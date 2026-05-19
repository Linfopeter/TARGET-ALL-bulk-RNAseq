### 11_ssGSEA_pathway_activity_PRIMARY_U18_TARGET_ALL_P2.R
### GOAL:
###   1) Calcular scores ssGSEA por muestra usando expresión VST de DESeq2
###   2) Evaluar si los clusters tienen actividad diferencial de vías por muestra
###   3) Usar MSigDB Hallmark como análisis principal, por ser compacto, interpretable y robusto
###   4) Generar tablas estadísticas, heatmaps y boxplots por cluster
###
### INPUT:
###   - subsets/primary_u18/final_analysis/07_differential_expression_clusters_deseq2/dds_cluster_de_ref2.rds
###   - subsets/primary_u18/final_analysis/07_differential_expression_clusters_deseq2/gene_annotation_table.csv
###
### OUTPUT:
###   - subsets/primary_u18/final_analysis/11_ssGSEA_pathway_activity/
###       * vst_expression_gene_symbol_collapsed.csv
###       * ssgsea_scores_hallmark.csv
###       * ssgsea_scores_hallmark_zscore.csv
###       * ssgsea_cluster_kruskal_results.csv
###       * ssgsea_cluster_pairwise_wilcox_results.csv
###       * ssgsea_cluster_mean_scores.csv
###       * ssgsea_cluster_mean_zscores.csv
###       * heatmap_ssgsea_hallmark_cluster_means.png
###       * heatmap_ssgsea_hallmark_samples_top_variable.png
###       * boxplot_ssgsea_top_pathways_by_cluster.png
###       * ssgsea_report.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(DESeq2)
  library(GSVA)
  library(msigdbr)
  library(pheatmap)
  library(RColorBrewer)
  library(ggplot2)
  library(tidyr)
})

cat("=== SCRIPT 11: ssGSEA PATHWAY ACTIVITY BY SAMPLE ===\n")

#---------------------------#
# 1) Paths y parámetros     #
#---------------------------#

indir <- "subsets/primary_u18/final_analysis/07_differential_expression_clusters_deseq2"
outdir <- "subsets/primary_u18/final_analysis/11_ssGSEA_pathway_activity"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

input_dds <- file.path(indir, "dds_cluster_de_ref2.rds")
input_annotation <- file.path(indir, "gene_annotation_table.csv")

# Parámetros ssGSEA
min_gene_set_size <- 10
max_gene_set_size <- 500
ssgsea_kcdf <- "Gaussian"

# Parámetros estadísticos
padj_cutoff <- 0.05
plot_top_n_pathways <- 25

# Paleta de clusters consistente con scripts previos
cluster_palette <- c(
  "1" = "#00BD00",
  "2" = "#FF6EB4",
  "3" = "#0000FF"
)

sample_type_palette_manual <- c(
  "Primary Blood Derived Cancer - Bone Marrow" = "#EEB422",
  "Primary Blood Derived Cancer - Peripheral Blood" = "#00B3EE"
)

#---------------------------#
# 2) Funciones auxiliares   #
#---------------------------#

clean_gene_symbol <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == "" | x == "NA" | grepl("^ENSG", x)] <- NA
  x
}

row_zscore <- function(mat) {
  z <- t(scale(t(mat)))
  z[is.na(z)] <- 0
  z
}

meta_with_sample_id <- function(meta_df) {
  meta_out <- as.data.frame(meta_df, check.names = FALSE)

  # Evita duplicar una columna sample_id preexistente en la metadata.
  if ("sample_id" %in% colnames(meta_out)) {
    meta_out <- meta_out[, colnames(meta_out) != "sample_id", drop = FALSE]
  }

  meta_out <- tibble::rownames_to_column(meta_out, "sample_id")
  return(meta_out)
}

collapse_expression_by_symbol <- function(expr_mat, annotation_df) {
  map_df <- annotation_df %>%
    dplyr::mutate(gene_symbol = clean_gene_symbol(gene_symbol)) %>%
    dplyr::filter(!is.na(gene_symbol)) %>%
    dplyr::select(ensembl_gene_id_original, gene_symbol) %>%
    dplyr::distinct(ensembl_gene_id_original, .keep_all = TRUE)

  common_genes <- intersect(rownames(expr_mat), map_df$ensembl_gene_id_original)

  if (length(common_genes) < 1000) {
    warning("Pocos genes comunes entre matriz VST y anotación: ", length(common_genes))
  }

  expr_sub <- expr_mat[common_genes, , drop = FALSE]
  map_sub <- map_df[match(common_genes, map_df$ensembl_gene_id_original), , drop = FALSE]

  expr_df <- as.data.frame(expr_sub, check.names = FALSE) %>%
    tibble::rownames_to_column("ensembl_gene_id_original") %>%
    dplyr::left_join(map_sub, by = "ensembl_gene_id_original") %>%
    dplyr::filter(!is.na(gene_symbol))

  sample_cols <- setdiff(colnames(expr_df), c("ensembl_gene_id_original", "gene_symbol"))

  # Colapsar símbolos duplicados usando el gen con mayor varianza entre muestras.
  variance_df <- expr_df %>%
    dplyr::mutate(row_variance = apply(dplyr::select(., dplyr::all_of(sample_cols)), 1, stats::var)) %>%
    dplyr::group_by(gene_symbol) %>%
    dplyr::arrange(dplyr::desc(row_variance), .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  collapsed_mat <- variance_df %>%
    dplyr::select(gene_symbol, dplyr::all_of(sample_cols)) %>%
    tibble::column_to_rownames("gene_symbol") %>%
    as.matrix()

  storage.mode(collapsed_mat) <- "numeric"
  return(collapsed_mat)
}

get_hallmark_gene_sets <- function() {
  cat("\n=== CARGANDO MSigDB Hallmark con msigdbr ===\n")

  msig_h <- tryCatch({
    msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  }, error = function(e1) {
    cat("[WARN] Falló msigdbr con category='H': ", conditionMessage(e1), "\n")
    tryCatch({
      msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
    }, error = function(e2) {
      stop("No se pudo cargar MSigDB Hallmark desde msigdbr: ", conditionMessage(e2))
    })
  })

  gene_col <- dplyr::case_when(
    "gene_symbol" %in% colnames(msig_h) ~ "gene_symbol",
    "human_gene_symbol" %in% colnames(msig_h) ~ "human_gene_symbol",
    TRUE ~ NA_character_
  )

  set_col <- dplyr::case_when(
    "gs_name" %in% colnames(msig_h) ~ "gs_name",
    TRUE ~ NA_character_
  )

  if (is.na(gene_col) || is.na(set_col)) {
    stop("No se encontraron columnas esperadas de msigdbr para genes o gene sets.")
  }

  gene_sets <- msig_h %>%
    dplyr::select(gs_name = dplyr::all_of(set_col), gene_symbol = dplyr::all_of(gene_col)) %>%
    dplyr::filter(!is.na(gene_symbol), gene_symbol != "") %>%
    dplyr::distinct(gs_name, gene_symbol) %>%
    split(x = .$gene_symbol, f = .$gs_name)

  gene_sets <- lapply(gene_sets, unique)

  cat("Gene sets Hallmark cargados:", length(gene_sets), "\n")
  return(gene_sets)
}

run_ssgsea_compatible <- function(expr_mat, gene_sets) {
  cat("\n=== EJECUTANDO ssGSEA ===\n")

  # Compatibilidad con versiones nuevas y antiguas de GSVA.
  res <- tryCatch({
    param <- GSVA::ssgseaParam(
      exprData = expr_mat,
      geneSets = gene_sets,
      minSize = min_gene_set_size,
      maxSize = max_gene_set_size,
      alpha = 0.25,
      normalize = TRUE
    )
    GSVA::gsva(param)
  }, error = function(e_new) {
    cat("[WARN] Falló API nueva ssgseaParam(); intentando API clásica gsva(..., method='ssgsea')\n")
    cat("Mensaje API nueva:", conditionMessage(e_new), "\n")

    tryCatch({
      GSVA::gsva(
        expr = expr_mat,
        gset.idx.list = gene_sets,
        method = "ssgsea",
        kcdf = ssgsea_kcdf,
        min.sz = min_gene_set_size,
        max.sz = max_gene_set_size,
        abs.ranking = FALSE,
        ssgsea.norm = TRUE,
        verbose = TRUE
      )
    }, error = function(e_old) {
      stop("Falló ssGSEA con API nueva y clásica de GSVA: ", conditionMessage(e_old))
    })
  })

  res <- as.matrix(res)
  return(res)
}

kruskal_by_cluster <- function(score_mat, meta_df) {
  score_df <- as.data.frame(t(score_mat), check.names = FALSE) %>%
    tibble::rownames_to_column("sample_id") %>%
    dplyr::left_join(
      meta_with_sample_id(meta_df),
      by = "sample_id"
    )

  pathway_cols <- rownames(score_mat)

  out <- lapply(pathway_cols, function(pw) {
    df <- score_df %>%
      dplyr::select(sample_id, cluster_k3, score = dplyr::all_of(pw)) %>%
      dplyr::filter(!is.na(cluster_k3), !is.na(score))

    if (length(unique(df$cluster_k3)) < 2) {
      return(tibble(pathway = pw, statistic = NA_real_, p.value = NA_real_))
    }

    kt <- kruskal.test(score ~ cluster_k3, data = df)

    tibble(
      pathway = pw,
      statistic = as.numeric(kt$statistic),
      p.value = kt$p.value
    )
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(p.adjust = p.adjust(p.value, method = "BH")) %>%
    dplyr::arrange(p.adjust)

  return(out)
}

pairwise_wilcox_by_cluster <- function(score_mat, meta_df) {
  score_df <- as.data.frame(t(score_mat), check.names = FALSE) %>%
    tibble::rownames_to_column("sample_id") %>%
    dplyr::left_join(
      meta_with_sample_id(meta_df),
      by = "sample_id"
    )

  pathway_cols <- rownames(score_mat)

  all_results <- lapply(pathway_cols, function(pw) {
    df <- score_df %>%
      dplyr::select(sample_id, cluster_k3, score = dplyr::all_of(pw)) %>%
      dplyr::filter(!is.na(cluster_k3), !is.na(score))

    clusters <- sort(unique(as.character(df$cluster_k3)))
    pairs <- combn(clusters, 2, simplify = FALSE)

    pair_results <- lapply(pairs, function(pr) {
      a <- pr[1]
      b <- pr[2]

      df_pair <- df %>% dplyr::filter(as.character(cluster_k3) %in% c(a, b))
      wt <- wilcox.test(score ~ cluster_k3, data = df_pair, exact = FALSE)

      mean_a <- mean(df_pair$score[df_pair$cluster_k3 == a], na.rm = TRUE)
      mean_b <- mean(df_pair$score[df_pair$cluster_k3 == b], na.rm = TRUE)

      tibble(
        pathway = pw,
        comparison = paste0("cluster_", a, "_vs_", b),
        cluster_a = a,
        cluster_b = b,
        mean_cluster_a = mean_a,
        mean_cluster_b = mean_b,
        delta_mean_a_minus_b = mean_a - mean_b,
        enriched_toward = ifelse(mean_a > mean_b, paste0("cluster_", a), paste0("cluster_", b)),
        p.value = wt$p.value
      )
    })

    dplyr::bind_rows(pair_results)
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(comparison) %>%
    dplyr::mutate(p.adjust = p.adjust(p.value, method = "BH")) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(comparison, p.adjust)

  return(all_results)
}

make_cluster_mean_heatmap <- function(mean_z_mat, outfile) {
  ann_col <- data.frame(
    cluster_k3 = factor(colnames(mean_z_mat), levels = c("1", "2", "3"))
  )
  rownames(ann_col) <- colnames(mean_z_mat)

  ann_colors <- list(cluster_k3 = cluster_palette[names(cluster_palette) %in% colnames(mean_z_mat)])

  pheatmap::pheatmap(
    mean_z_mat,
    scale = "none",
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
    breaks = seq(-2, 2, length.out = 101),
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    border_color = NA,
    fontsize_row = 8,
    fontsize_col = 11,
    main = "ssGSEA Hallmark pathway activity by cluster\n(mean sample-level Z-score)",
    filename = outfile,
    width = 8,
    height = 11
  )
}

make_sample_heatmap_top_variable <- function(score_z_mat, meta_df, outfile, top_n = 25) {
  pathway_var <- apply(score_z_mat, 1, stats::var)
  top_pathways <- names(sort(pathway_var, decreasing = TRUE))[seq_len(min(top_n, length(pathway_var)))]

  mat <- score_z_mat[top_pathways, , drop = FALSE]

  meta_df <- meta_df[colnames(mat), , drop = FALSE]

  annotation_col <- data.frame(
    cluster_k3 = factor(meta_df$cluster_k3)
  )
  rownames(annotation_col) <- colnames(mat)

  if ("sample_type.x" %in% colnames(meta_df)) {
    annotation_col$sample_type <- meta_df$sample_type.x
  } else if ("sample_type" %in% colnames(meta_df)) {
    annotation_col$sample_type <- meta_df$sample_type
  }

  ann_colors <- list()
  ann_colors$cluster_k3 <- cluster_palette[levels(annotation_col$cluster_k3)]

  if ("sample_type" %in% colnames(annotation_col)) {
    sample_type_levels <- unique(as.character(annotation_col$sample_type))
    sample_cols <- sample_type_palette_manual
    missing_levels <- setdiff(sample_type_levels, names(sample_cols))
    if (length(missing_levels) > 0) {
      extra_cols <- colorRampPalette(c("#7E6148", "#B09C85", "#6A3D9A"))(length(missing_levels))
      names(extra_cols) <- missing_levels
      sample_cols <- c(sample_cols, extra_cols)
    }
    ann_colors$sample_type <- sample_cols[sample_type_levels]
  }

  pheatmap::pheatmap(
    mat,
    scale = "none",
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
    breaks = seq(-2, 2, length.out = 101),
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    show_colnames = FALSE,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    border_color = NA,
    fontsize_row = 8,
    main = paste0("Top ", top_n, " most variable Hallmark ssGSEA pathways"),
    filename = outfile,
    width = 13,
    height = 9
  )
}

make_boxplots_top_pathways <- function(score_mat, meta_df, kruskal_df, outfile, top_n = 12) {
  top_pathways <- kruskal_df %>%
    dplyr::filter(!is.na(p.adjust)) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::pull(pathway)

  if (length(top_pathways) == 0) {
    cat("[INFO] No hay pathways para boxplot.\n")
    return(NULL)
  }

  plot_df <- as.data.frame(t(score_mat[top_pathways, , drop = FALSE]), check.names = FALSE) %>%
    tibble::rownames_to_column("sample_id") %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(top_pathways),
      names_to = "pathway",
      values_to = "score"
    ) %>%
    dplyr::left_join(
      meta_with_sample_id(meta_df),
      by = "sample_id"
    ) %>%
    dplyr::mutate(
      pathway = stringr::str_replace(pathway, "^HALLMARK_", ""),
      pathway = stringr::str_replace_all(pathway, "_", " "),
      pathway = stringr::str_to_title(pathway),
      cluster_k3 = factor(cluster_k3, levels = c("1", "2", "3"))
    )

  p <- ggplot(plot_df, aes(x = cluster_k3, y = score, fill = cluster_k3)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.15, alpha = 0.45, size = 0.9) +
    facet_wrap(~ pathway, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = cluster_palette) +
    labs(
      title = "Top differential Hallmark ssGSEA scores by cluster",
      x = "Cluster",
      y = "ssGSEA score",
      fill = "Cluster"
    ) +
    theme_bw(base_size = 12) +
    theme(
      strip.text = element_text(size = 9, face = "bold"),
      legend.position = "bottom"
    )

  ggsave(outfile, plot = p, width = 13, height = 10, dpi = 300)
}

#---------------------------#
# 3) Cargar datos           #
#---------------------------#

cat("\n=== CARGA DE DATOS ===\n")

if (!file.exists(input_dds)) stop("No existe dds: ", input_dds)
if (!file.exists(input_annotation)) stop("No existe anotación: ", input_annotation)

dds <- readRDS(input_dds)
annotation_df <- read.csv(input_annotation, stringsAsFactors = FALSE, check.names = FALSE)

if (!"cluster_k3" %in% colnames(colData(dds))) {
  stop("dds no contiene cluster_k3 en colData.")
}

meta_df <- as.data.frame(colData(dds))
meta_df$cluster_k3 <- factor(meta_df$cluster_k3, levels = sort(unique(as.character(meta_df$cluster_k3))))

cat("Muestras:", ncol(dds), "\n")
cat("Genes:", nrow(dds), "\n")
cat("Distribución de clusters:\n")
print(table(meta_df$cluster_k3))

#---------------------------#
# 4) VST y colapso a símbolos#
#---------------------------#

cat("\n=== VST EXPRESSION ===\n")

vsd <- DESeq2::vst(dds, blind = FALSE)
vst_mat <- assay(vsd)

expr_symbol_mat <- collapse_expression_by_symbol(vst_mat, annotation_df)

cat("Genes con símbolo HGNC tras colapso:", nrow(expr_symbol_mat), "\n")
cat("Muestras en matriz:", ncol(expr_symbol_mat), "\n")

write.csv(
  data.frame(gene_symbol = rownames(expr_symbol_mat), expr_symbol_mat, check.names = FALSE),
  file = file.path(outdir, "vst_expression_gene_symbol_collapsed.csv"),
  row.names = FALSE
)

#---------------------------#
# 5) Gene sets Hallmark     #
#---------------------------#

gene_sets_hallmark <- get_hallmark_gene_sets()

gene_sets_hallmark <- gene_sets_hallmark[lengths(gene_sets_hallmark) >= min_gene_set_size]

gene_set_overlap <- sapply(gene_sets_hallmark, function(gs) {
  sum(gs %in% rownames(expr_symbol_mat))
})

gene_sets_hallmark <- gene_sets_hallmark[
  gene_set_overlap >= min_gene_set_size & gene_set_overlap <= max_gene_set_size
]

cat("Gene sets Hallmark tras filtro por overlap:", length(gene_sets_hallmark), "\n")

if (length(gene_sets_hallmark) < 5) {
  stop("Muy pocos gene sets Hallmark disponibles tras filtro. Revisa símbolos génicos.")
}

#---------------------------#
# 6) ssGSEA                 #
#---------------------------#

ssgsea_scores <- run_ssgsea_compatible(expr_symbol_mat, gene_sets_hallmark)

# Alinear metadata
meta_df <- meta_df[colnames(ssgsea_scores), , drop = FALSE]

write.csv(
  data.frame(pathway = rownames(ssgsea_scores), ssgsea_scores, check.names = FALSE),
  file = file.path(outdir, "ssgsea_scores_hallmark.csv"),
  row.names = FALSE
)

ssgsea_scores_z <- row_zscore(ssgsea_scores)

write.csv(
  data.frame(pathway = rownames(ssgsea_scores_z), ssgsea_scores_z, check.names = FALSE),
  file = file.path(outdir, "ssgsea_scores_hallmark_zscore.csv"),
  row.names = FALSE
)

#---------------------------#
# 7) Estadística por cluster#
#---------------------------#

cat("\n=== ESTADÍSTICA ssGSEA POR CLUSTER ===\n")

kruskal_df <- kruskal_by_cluster(ssgsea_scores, meta_df)
write.csv(kruskal_df, file.path(outdir, "ssgsea_cluster_kruskal_results.csv"), row.names = FALSE)

pairwise_df <- pairwise_wilcox_by_cluster(ssgsea_scores, meta_df)
write.csv(pairwise_df, file.path(outdir, "ssgsea_cluster_pairwise_wilcox_results.csv"), row.names = FALSE)

cat("Pathways Hallmark significativos por Kruskal p.adjust <", padj_cutoff, ":", sum(kruskal_df$p.adjust < padj_cutoff, na.rm = TRUE), "\n")

#---------------------------#
# 8) Promedios por cluster  #
#---------------------------#

score_long <- as.data.frame(t(ssgsea_scores), check.names = FALSE) %>%
  tibble::rownames_to_column("sample_id") %>%
  tidyr::pivot_longer(
    cols = -sample_id,
    names_to = "pathway",
    values_to = "score"
  ) %>%
  dplyr::left_join(meta_with_sample_id(meta_df), by = "sample_id")

cluster_mean_df <- score_long %>%
  dplyr::group_by(pathway, cluster_k3) %>%
  dplyr::summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop")

write.csv(cluster_mean_df, file.path(outdir, "ssgsea_cluster_mean_scores.csv"), row.names = FALSE)

cluster_mean_mat <- cluster_mean_df %>%
  tidyr::pivot_wider(names_from = cluster_k3, values_from = mean_score) %>%
  tibble::column_to_rownames("pathway") %>%
  as.matrix()

cluster_mean_z_mat <- row_zscore(cluster_mean_mat)

write.csv(
  data.frame(pathway = rownames(cluster_mean_z_mat), cluster_mean_z_mat, check.names = FALSE),
  file = file.path(outdir, "ssgsea_cluster_mean_zscores.csv"),
  row.names = FALSE
)

#---------------------------#
# 9) Figuras                #
#---------------------------#

cat("\n=== FIGURAS ssGSEA ===\n")

make_cluster_mean_heatmap(
  mean_z_mat = cluster_mean_z_mat,
  outfile = file.path(outdir, "heatmap_ssgsea_hallmark_cluster_means.png")
)

make_sample_heatmap_top_variable(
  score_z_mat = ssgsea_scores_z,
  meta_df = meta_df,
  outfile = file.path(outdir, "heatmap_ssgsea_hallmark_samples_top_variable.png"),
  top_n = plot_top_n_pathways
)

make_boxplots_top_pathways(
  score_mat = ssgsea_scores,
  meta_df = meta_df,
  kruskal_df = kruskal_df,
  outfile = file.path(outdir, "boxplot_ssgsea_top_pathways_by_cluster.png"),
  top_n = 12
)

#---------------------------#
# 10) Reporte               #
#---------------------------#

cat("\n=== GENERANDO REPORTE TXT ===\n")

report_file <- file.path(outdir, "ssgsea_report.txt")
sink(report_file)

cat("REPORTE ssGSEA POR MUESTRA\n")
cat("Fecha:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("INPUTS:\n")
cat("- DDS:", input_dds, "\n")
cat("- Anotación:", input_annotation, "\n\n")

cat("DECISIÓN METODOLÓGICA:\n")
cat("- ssGSEA se usa como validación por muestra de actividad funcional de vías.\n")
cat("- No reemplaza GSEA pre-ranked del script 09; lo complementa.\n")
cat("- GSEA responde qué vías están enriquecidas en una comparación.\n")
cat("- ssGSEA responde si cada muestra individual tiene mayor o menor actividad de una vía.\n")
cat("- Se usa MSigDB Hallmark como colección principal por ser compacta, interpretable y menos redundante que GO BP completo.\n\n")

cat("PARÁMETROS:\n")
cat("- min_gene_set_size:", min_gene_set_size, "\n")
cat("- max_gene_set_size:", max_gene_set_size, "\n")
cat("- padj_cutoff:", padj_cutoff, "\n")
cat("- score input: DESeq2 VST, colapsado a símbolo HGNC por mayor varianza.\n\n")

cat("DIMENSIONES:\n")
cat("- Genes VST originales:", nrow(vst_mat), "\n")
cat("- Genes HGNC tras colapso:", nrow(expr_symbol_mat), "\n")
cat("- Muestras:", ncol(expr_symbol_mat), "\n")
cat("- Gene sets Hallmark usados:", length(gene_sets_hallmark), "\n\n")

cat("DISTRIBUCIÓN DE CLUSTERS:\n")
print(table(meta_df$cluster_k3))
cat("\n")

cat("RESULTADOS GLOBALES KRUSKAL-WALLIS:\n")
cat("- Pathways significativos p.adjust <", padj_cutoff, ":", sum(kruskal_df$p.adjust < padj_cutoff, na.rm = TRUE), "\n\n")

cat("TOP 20 PATHWAYS DIFERENCIALES POR CLUSTER:\n")
print(utils::head(kruskal_df, 20))
cat("\n")

cat("INTERPRETACIÓN RECOMENDADA:\n")
cat("- Usar heatmap de medias por cluster para resumir programas funcionales dominantes.\n")
cat("- Usar heatmap por muestra para verificar que la señal no dependa de pocas muestras.\n")
cat("- Usar boxplots para confirmar distribución por cluster.\n")
cat("- Priorizar vías con p.adjust < 0.05 y separación visual consistente entre clusters.\n\n")

cat("ARCHIVOS GENERADOS:\n")
cat("- vst_expression_gene_symbol_collapsed.csv\n")
cat("- ssgsea_scores_hallmark.csv\n")
cat("- ssgsea_scores_hallmark_zscore.csv\n")
cat("- ssgsea_cluster_kruskal_results.csv\n")
cat("- ssgsea_cluster_pairwise_wilcox_results.csv\n")
cat("- ssgsea_cluster_mean_scores.csv\n")
cat("- ssgsea_cluster_mean_zscores.csv\n")
cat("- heatmap_ssgsea_hallmark_cluster_means.png\n")
cat("- heatmap_ssgsea_hallmark_samples_top_variable.png\n")
cat("- boxplot_ssgsea_top_pathways_by_cluster.png\n")
cat("- ssgsea_report.txt\n")

sink()

cat("\n=== SCRIPT 11 COMPLETADO ===\n")
cat("Archivos generados en:\n")
cat(outdir, "\n")
