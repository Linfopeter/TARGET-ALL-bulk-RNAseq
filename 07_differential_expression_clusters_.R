### 07_differential_expression_clusters_deseq2_annotated_apeglm_PRIMARY_U18_TARGET_ALL_P2.R
### GOAL:
###   1) Realizar expresión diferencial entre los 3 clusters obtenidos en script 06
###   2) Usar apeglm en las 3 comparaciones para shrinkage de log2FC
###   3) Anotar genes ENSEMBL -> símbolo génico ANTES de exportar tablas o graficar
###   4) Generar tablas completas de resultados DESeq2
###   5) Visualizar resultados con volcano plots + heatmaps
###
### INPUT:
###   - subsets/primary_u18/final_analysis/TARGET_ALL_P2_PRIMARY_U18_FINAL_FOR_ANALYSIS.RData
###   - subsets/primary_u18/final_analysis/06_qc_normalization_deseq2/meta_with_clusters.rds
###
### OUTPUT:
###   - subsets/primary_u18/final_analysis/07_differential_expression_clusters_deseq2/
###       * gene_annotation_table.csv
###       * dds_cluster_de_ref2.rds
###       * dds_cluster_de_ref3.rds
###       * normalized_counts_cluster_de_annotated.csv
###       * results_cluster_1_vs_2.csv
###       * results_cluster_1_vs_3.csv
###       * results_cluster_2_vs_3.csv
###       * sig_cluster_1_vs_2.csv
###       * sig_cluster_1_vs_3.csv
###       * sig_cluster_2_vs_3.csv
###       * up_cluster_1_vs_2.csv
###       * up_cluster_1_vs_3.csv
###       * up_cluster_2_vs_3.csv
###       * down_cluster_1_vs_2.csv
###       * down_cluster_1_vs_3.csv
###       * down_cluster_2_vs_3.csv
###       * volcano_cluster_1_vs_2.png
###       * volcano_cluster_1_vs_3.png
###       * volcano_cluster_2_vs_3.png
###       * heatmap_topgenes_cluster_1_vs_2.png
###       * heatmap_topgenes_cluster_1_vs_3.png
###       * heatmap_topgenes_cluster_2_vs_3.png
###       * differential_expression_clusters_report.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(S4Vectors)
  library(biomaRt)
  library(apeglm)
})

#---------------------------#
# 1) Paths y directorios    #
#---------------------------#
input_counts_file <- "subsets/primary_u18/final_analysis/TARGET_ALL_P2_PRIMARY_U18_FINAL_FOR_ANALYSIS.RData"
input_meta_clusters <- "subsets/primary_u18/final_analysis/06_qc_normalization_deseq2/meta_with_clusters.rds"
outdir <- "subsets/primary_u18/final_analysis/07_differential_expression_clusters_deseq2"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("=== SCRIPT 07: DIFFERENTIAL EXPRESSION BETWEEN CLUSTERS (ANNOTATED + APEGLM) ===\n")

#---------------------------#
# 2) Cargar datos           #
#---------------------------#
cat("\n=== CARGA DE DATOS ===\n")

load(input_counts_file)
meta_clusters <- readRDS(input_meta_clusters)

cat("\nDimensiones iniciales:\n")
cat("- counts_final filas:", nrow(counts_final), "\n")
cat("- counts_final columnas:", ncol(counts_final), "\n")
cat("- meta_clusters filas:", nrow(meta_clusters), "\n")
cat("- meta_clusters columnas:", ncol(meta_clusters), "\n")

if (!"cluster_k3" %in% colnames(meta_clusters)) {
  stop("La metadata cargada no contiene la columna 'cluster_k3'.")
}

cat("\nChequeo de clusters:\n")
print(table(meta_clusters$cluster_k3, useNA = "ifany"))

#---------------------------#
# 3) Limpiar counts         #
#---------------------------#
cat("\n=== LIMPIEZA DE MATRIZ DE COUNTS ===\n")

star_technical_rows <- c(
  "N_unmapped",
  "N_multimapping",
  "N_noFeature",
  "N_ambiguous"
)

present_star_rows <- intersect(rownames(counts_final), star_technical_rows)

cat("\nFilas técnicas STAR encontradas:\n")
print(present_star_rows)

counts_gene_only <- counts_final[!rownames(counts_final) %in% star_technical_rows, , drop = FALSE]

cat("\nDimensiones después de remover filas técnicas:\n")
cat("- Filas antes:", nrow(counts_final), "\n")
cat("- Filas después:", nrow(counts_gene_only), "\n")
cat("- Columnas:", ncol(counts_gene_only), "\n")

if (!all(abs(counts_gene_only - round(counts_gene_only)) < .Machine$double.eps^0.5)) {
  stop("La matriz counts_gene_only contiene valores no enteros.")
}

counts_gene_only <- round(counts_gene_only)
storage.mode(counts_gene_only) <- "integer"

#---------------------------#
# 4) Alinear metadata       #
#---------------------------#
cat("\n=== ALINEACIÓN DE METADATA ===\n")

meta_deseq <- meta_clusters
meta_deseq <- meta_deseq[, !sapply(meta_deseq, is.list), drop = FALSE]

if (!"counts_colname" %in% colnames(meta_deseq)) {
  stop("No existe la columna counts_colname en meta_deseq.")
}

meta_deseq <- meta_deseq[match(colnames(counts_gene_only), meta_deseq$counts_colname), , drop = FALSE]

if (any(is.na(meta_deseq$counts_colname))) {
  stop("Hay muestras en counts_gene_only que no pudieron alinearse con meta_deseq.")
}

cat("\nChequeo de alineación:\n")
cat("¿colnames(counts_gene_only) == meta_deseq$counts_colname?: ",
    all(colnames(counts_gene_only) == meta_deseq$counts_colname), "\n")

meta_deseq$cluster_k3 <- as.factor(meta_deseq$cluster_k3)
meta_deseq$cluster_k3 <- droplevels(meta_deseq$cluster_k3)

cat("\nNiveles de cluster_k3:\n")
print(levels(meta_deseq$cluster_k3))

if (length(levels(meta_deseq$cluster_k3)) != 3) {
  stop("Se esperaban exactamente 3 clusters en cluster_k3.")
}

rownames(meta_deseq) <- meta_deseq$counts_colname

#---------------------------#
# 5) Crear DDS base         #
#---------------------------#
cat("\n=== CREACIÓN DE DESeqDataSet BASE ===\n")

coldata_s4 <- S4Vectors::DataFrame(meta_deseq)

dds_base <- DESeqDataSetFromMatrix(
  countData = counts_gene_only,
  colData = coldata_s4,
  design = ~ cluster_k3
)

cat("\nDimensiones iniciales del DDS:\n")
cat("- Genes:", nrow(dds_base), "\n")
cat("- Muestras:", ncol(dds_base), "\n")

#---------------------------#
# 6) Filtrado baja expresión#
#---------------------------#
cat("\n=== FILTRADO DE BAJA EXPRESIÓN ===\n")

keep_genes <- rowSums(counts(dds_base) >= 10) >= 10

cat("\nGenes antes del filtrado:", nrow(dds_base), "\n")
cat("Genes retenidos:", sum(keep_genes), "\n")
cat("Genes removidos:", sum(!keep_genes), "\n")

dds_base <- dds_base[keep_genes, ]

cat("\nDimensiones después del filtrado:\n")
cat("- Genes:", nrow(dds_base), "\n")
cat("- Muestras:", ncol(dds_base), "\n")

#---------------------------#
# 7) Anotación génica       #
#---------------------------#
cat("\n=== ANOTACIÓN ENSEMBL -> GENE SYMBOL ===\n")

ensembl_ids_original <- rownames(dds_base)
ensembl_ids_clean <- sub("\\..*$", "", ensembl_ids_original)

annotation_df <- data.frame(
  ensembl_gene_id_original = ensembl_ids_original,
  ensembl_gene_id = ensembl_ids_clean,
  stringsAsFactors = FALSE
)

connect_to_ensembl_biomart <- function() {
  mart <- NULL

  mirrors_to_try <- c("www", "useast", "asia")

  for (m in mirrors_to_try) {
    cat("\nIntentando conexión con mirror:", m, "\n")
    mart <- tryCatch({
      useEnsembl(
        biomart = "genes",
        dataset = "hsapiens_gene_ensembl",
        mirror = m
      )
    }, error = function(e) {
      cat("  Falló mirror", m, "->", conditionMessage(e), "\n")
      NULL
    })

    if (!is.null(mart)) {
      cat("  Conexión exitosa con mirror:", m, "\n")
      return(mart)
    }
  }

  archive_hosts <- c(
    "https://www.ensembl.org",
    "https://apr2024.archive.ensembl.org",
    "https://oct2024.archive.ensembl.org",
    "https://may2025.archive.ensembl.org"
  )

  for (h in archive_hosts) {
    cat("\nIntentando conexión con host:", h, "\n")
    mart <- tryCatch({
      useEnsembl(
        biomart = "genes",
        dataset = "hsapiens_gene_ensembl",
        host = h
      )
    }, error = function(e) {
      cat("  Falló host", h, "->", conditionMessage(e), "\n")
      NULL
    })

    if (!is.null(mart)) {
      cat("  Conexión exitosa con host:", h, "\n")
      return(mart)
    }
  }

  cat("\nIntentando ajuste SSL de biomaRt...\n")

  tryCatch({
    setEnsemblSSL(list(
      ssl_verifypeer = FALSE,
      ssl_cipher_list = "DEFAULT@SECLEVEL=1"
    ))
  }, error = function(e) {
    cat("  No se pudo aplicar setEnsemblSSL():", conditionMessage(e), "\n")
  })

  for (m in mirrors_to_try) {
    cat("\nReintentando con mirror y SSL ajustado:", m, "\n")
    mart <- tryCatch({
      useEnsembl(
        biomart = "genes",
        dataset = "hsapiens_gene_ensembl",
        mirror = m
      )
    }, error = function(e) {
      cat("  Falló mirror", m, "con SSL ajustado ->", conditionMessage(e), "\n")
      NULL
    })

    if (!is.null(mart)) {
      cat("  Conexión exitosa con mirror:", m, "tras ajuste SSL\n")
      return(mart)
    }
  }

  return(NULL)
}

cat("\nConectando a Ensembl biomart...\n")
mart <- connect_to_ensembl_biomart()

if (is.null(mart)) {
  stop(
    paste(
      "No fue posible conectar con Ensembl/BioMart tras probar mirrors, hosts archivados y ajuste SSL.",
      "Puedes volver a correr más tarde o usar una tabla local de anotación."
    )
  )
}

cat("\nRecuperando anotación con getBM...\n")

anno_bm <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "description"),
  filters = "ensembl_gene_id",
  values = unique(annotation_df$ensembl_gene_id),
  mart = mart
)

anno_bm <- anno_bm %>%
  dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)

annotation_df <- annotation_df %>%
  dplyr::left_join(anno_bm, by = "ensembl_gene_id") %>%
  dplyr::mutate(
    hgnc_symbol = ifelse(is.na(hgnc_symbol), "", hgnc_symbol),
    description = ifelse(is.na(description), "", description),
    gene_symbol = ifelse(hgnc_symbol == "", NA, hgnc_symbol),
    gene_label  = ifelse(is.na(gene_symbol), NA, gene_symbol)
  )

write.csv(
  annotation_df,
  file = file.path(outdir, "gene_annotation_table.csv"),
  row.names = FALSE
)

cat("\nResumen de anotación:\n")
cat("- Genes filtrados:", nrow(annotation_df), "\n")
cat("- Genes con símbolo HGNC:", sum(annotation_df$hgnc_symbol != ""), "\n")
cat("- Genes sin símbolo HGNC:", sum(annotation_df$hgnc_symbol == ""), "\n")

#---------------------------#
# 8) Funciones auxiliares   #
#---------------------------#

annotate_results_table <- function(res_df, annotation_df) {
  res_df$ensembl_gene_id_original <- rownames(res_df)

  res_df <- as.data.frame(res_df) %>%
    dplyr::left_join(
      annotation_df %>%
        dplyr::select(
          ensembl_gene_id_original,
          ensembl_gene_id,
          gene_symbol,
          description
        ),
      by = "ensembl_gene_id_original"
    ) %>%
    dplyr::relocate(ensembl_gene_id_original, ensembl_gene_id, gene_symbol, description) %>%
    dplyr::arrange(padj)

  return(res_df)
}

make_volcano_plot <- function(res_df, comparison_label, outfile) {
  plot_df <- res_df %>%
    dplyr::mutate(
      padj_plot = ifelse(is.na(padj) | padj <= 0, 1, padj),
      neglog10_padj = -log10(padj_plot),
      status = dplyr::case_when(
        !is.na(padj) & padj < 0.05 & log2FoldChange >= 1  ~ "Up",
        !is.na(padj) & padj < 0.05 & log2FoldChange <= -1 ~ "Down",
        TRUE ~ "NS"
      )
    )

  p <- ggplot(plot_df, aes(x = log2FoldChange, y = neglog10_padj, color = status)) +
    geom_point(alpha = 0.7, size = 1.8) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    scale_color_manual(values = c("Up" = "#D73027", "Down" = "#4575B4", "NS" = "grey70")) +
    labs(
      title = paste0("Volcano plot: ", comparison_label),
      x = "log2(Fold Change)",
      y = "-log10(adjusted p-value)",
      color = "Clasificación"
    ) +
    theme_bw(base_size = 12)

  ggsave(
    filename = outfile,
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}

make_heatmap <- function(dds_obj, sig_df, meta_df, annotation_df, comparison_label,
                         cluster_a, cluster_b, outfile, top_n = 30) {

  sig_df2 <- sig_df %>%
    dplyr::filter(
      !is.na(padj),
      !is.na(gene_symbol),
      gene_symbol != "",
      !grepl("^ENSG", gene_symbol)
    ) %>%
    dplyr::arrange(padj, dplyr::desc(abs(log2FoldChange)))

  top_df <- head(sig_df2, top_n)

  if (nrow(top_df) < 2) {
    cat("\n[INFO] No hay suficientes genes con símbolo HGNC válido para heatmap en", comparison_label, "\n")
    return(NULL)
  }

  vsd <- vst(dds_obj, blind = FALSE)
  vst_mat <- assay(vsd)

  selected_samples <- rownames(meta_df)[meta_df$cluster_k3 %in% c(cluster_a, cluster_b)]
  selected_samples <- intersect(selected_samples, colnames(vst_mat))

  selected_genes <- intersect(top_df$ensembl_gene_id_original, rownames(vst_mat))

  if (length(selected_genes) < 2) {
    cat("\n[INFO] No hay suficientes genes presentes en la matriz VST para", comparison_label, "\n")
    return(NULL)
  }

  mat <- vst_mat[selected_genes, selected_samples, drop = FALSE]

  row_labels <- top_df$gene_symbol[
    match(rownames(mat), top_df$ensembl_gene_id_original)
  ]

  keep_rows <- !is.na(row_labels) & row_labels != "" & !grepl("^ENSG", row_labels)

  mat <- mat[keep_rows, , drop = FALSE]
  row_labels <- row_labels[keep_rows]

  if (nrow(mat) < 2) {
    cat("\n[INFO] Tras remover genes sin símbolo HGNC válido quedan <2 genes en", comparison_label, "\n")
    return(NULL)
  }

  rownames(mat) <- make.unique(row_labels)

  annotation_col <- data.frame(
    cluster_k3 = meta_df[selected_samples, "cluster_k3", drop = TRUE]
  )
  rownames(annotation_col) <- selected_samples

  pheatmap(
    mat,
    scale = "row",
    annotation_col = annotation_col,
    show_rownames = TRUE,
    show_colnames = FALSE,
    clustering_method = "complete",
    color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
    fontsize_row = 8,
    main = paste0("Top genes DE con símbolo HGNC válido: ", comparison_label),
    filename = outfile,
    width = 8,
    height = 10
  )
}

run_deseq_with_reference <- function(dds_input, ref_level) {
  dds_tmp <- dds_input
  dds_tmp$cluster_k3 <- relevel(dds_tmp$cluster_k3, ref = ref_level)
  dds_tmp <- DESeq(dds_tmp)
  return(dds_tmp)
}

extract_results_apeglm <- function(dds_obj, coef_name, comparison_label, outdir, meta_df,
                                   annotation_df, cluster_a, cluster_b) {

  cat("\n=== PROCESANDO:", comparison_label, "===\n")
  cat("Coeficiente usado:", coef_name, "\n")

  res <- results(dds_obj, name = coef_name, alpha = 0.05)
  res_shrunk <- lfcShrink(dds_obj, coef = coef_name, res = res, type = "apeglm")

  res_df <- annotate_results_table(res_shrunk, annotation_df)

  comp_name <- gsub(" ", "_", tolower(comparison_label))
  comp_name <- gsub("-", "_", comp_name)

  write.csv(
    res_df,
    file = file.path(outdir, paste0("results_", comp_name, ".csv")),
    row.names = FALSE
  )

  sig_df <- res_df %>%
    dplyr::filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1)

  up_df <- sig_df %>%
    dplyr::filter(log2FoldChange >= 1) %>%
    dplyr::arrange(padj)

  down_df <- sig_df %>%
    dplyr::filter(log2FoldChange <= -1) %>%
    dplyr::arrange(padj)

  write.csv(
    sig_df,
    file = file.path(outdir, paste0("sig_", comp_name, ".csv")),
    row.names = FALSE
  )

  write.csv(
    up_df,
    file = file.path(outdir, paste0("up_", comp_name, ".csv")),
    row.names = FALSE
  )

  write.csv(
    down_df,
    file = file.path(outdir, paste0("down_", comp_name, ".csv")),
    row.names = FALSE
  )

  make_volcano_plot(
    res_df = res_df,
    comparison_label = comparison_label,
    outfile = file.path(outdir, paste0("volcano_", comp_name, ".png"))
  )

  make_heatmap(
    dds_obj = dds_obj,
    sig_df = sig_df,
    meta_df = meta_df,
    annotation_df = annotation_df,
    comparison_label = comparison_label,
    cluster_a = cluster_a,
    cluster_b = cluster_b,
    outfile = file.path(outdir, paste0("heatmap_topgenes_", comp_name, ".png")),
    top_n = 30
  )

  summary_list <- list(
    comparison = comparison_label,
    coefficient = coef_name,
    total_genes_tested = nrow(res_df),
    significant_genes = nrow(sig_df),
    up_genes = nrow(up_df),
    down_genes = nrow(down_df)
  )

  return(list(
    res_df = res_df,
    sig_df = sig_df,
    up_df = up_df,
    down_df = down_df,
    summary = summary_list
  ))
}

#---------------------------#
# 9) Correr DESeq2          #
#---------------------------#
cat("\n=== EJECUTANDO DESEQ2 CON REFERENCIAS NECESARIAS PARA APEGLM ===\n")

dds_ref2 <- run_deseq_with_reference(dds_base, ref_level = "2")
cat("\nCoeficientes disponibles con ref=2:\n")
print(resultsNames(dds_ref2))

dds_ref3 <- run_deseq_with_reference(dds_base, ref_level = "3")
cat("\nCoeficientes disponibles con ref=3:\n")
print(resultsNames(dds_ref3))

saveRDS(dds_ref2, file = file.path(outdir, "dds_cluster_de_ref2.rds"))
saveRDS(dds_ref3, file = file.path(outdir, "dds_cluster_de_ref3.rds"))

normalized_counts <- counts(dds_ref2, normalized = TRUE)

normalized_counts_df <- data.frame(
  ensembl_gene_id_original = rownames(normalized_counts),
  normalized_counts,
  check.names = FALSE
) %>%
  dplyr::left_join(
    annotation_df %>%
      dplyr::select(ensembl_gene_id_original, ensembl_gene_id, gene_symbol),
    by = "ensembl_gene_id_original"
  ) %>%
  dplyr::relocate(ensembl_gene_id_original, ensembl_gene_id, gene_symbol)

write.csv(
  normalized_counts_df,
  file = file.path(outdir, "normalized_counts_cluster_de_annotated.csv"),
  row.names = FALSE
)

#---------------------------#
# 10) Comparaciones         #
#---------------------------#
cat("\n=== COMPARACIONES ENTRE CLUSTERS ===\n")

meta_for_plots_ref2 <- as.data.frame(colData(dds_ref2))
meta_for_plots_ref2$cluster_k3 <- as.factor(meta_for_plots_ref2$cluster_k3)

meta_for_plots_ref3 <- as.data.frame(colData(dds_ref3))
meta_for_plots_ref3$cluster_k3 <- as.factor(meta_for_plots_ref3$cluster_k3)

res_1_vs_2 <- extract_results_apeglm(
  dds_obj = dds_ref2,
  coef_name = "cluster_k3_1_vs_2",
  comparison_label = "cluster 1 vs 2",
  outdir = outdir,
  meta_df = meta_for_plots_ref2,
  annotation_df = annotation_df,
  cluster_a = "1",
  cluster_b = "2"
)

res_1_vs_3 <- extract_results_apeglm(
  dds_obj = dds_ref3,
  coef_name = "cluster_k3_1_vs_3",
  comparison_label = "cluster 1 vs 3",
  outdir = outdir,
  meta_df = meta_for_plots_ref3,
  annotation_df = annotation_df,
  cluster_a = "1",
  cluster_b = "3"
)

res_2_vs_3 <- extract_results_apeglm(
  dds_obj = dds_ref3,
  coef_name = "cluster_k3_2_vs_3",
  comparison_label = "cluster 2 vs 3",
  outdir = outdir,
  meta_df = meta_for_plots_ref3,
  annotation_df = annotation_df,
  cluster_a = "2",
  cluster_b = "3"
)

#---------------------------#
# 11) Reporte TXT           #
#---------------------------#
cat("\n=== GENERANDO REPORTE TXT ===\n")

report_file <- file.path(outdir, "differential_expression_clusters_report.txt")
sink(report_file)

cat("REPORTE DE EXPRESIÓN DIFERENCIAL ENTRE CLUSTERS (DESEQ2 + ANOTACIÓN + APEGLM)\n")
cat("Fecha:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("INPUT COUNTS FILE:\n")
cat(input_counts_file, "\n\n")

cat("INPUT METADATA WITH CLUSTERS:\n")
cat(input_meta_clusters, "\n\n")

cat("DIMENSIONES INICIALES:\n")
cat("- counts_final filas:", nrow(counts_final), "\n")
cat("- counts_final columnas:", ncol(counts_final), "\n")
cat("- meta_clusters filas:", nrow(meta_clusters), "\n")
cat("- meta_clusters columnas:", ncol(meta_clusters), "\n\n")

cat("FILAS TÉCNICAS STAR REMOVIDAS:\n")
print(present_star_rows)
cat("\n")

cat("DIMENSIONES DESPUÉS DE REMOVER FILAS TÉCNICAS:\n")
cat("- counts_gene_only filas:", nrow(counts_gene_only), "\n")
cat("- counts_gene_only columnas:", ncol(counts_gene_only), "\n\n")

cat("FILTRADO DE BAJA EXPRESIÓN:\n")
cat("- genes antes:", length(keep_genes), "\n")
cat("- genes retenidos:", sum(keep_genes), "\n")
cat("- genes removidos:", sum(!keep_genes), "\n\n")

cat("DISTRIBUCIÓN DE CLUSTERS:\n")
print(table(meta_deseq$cluster_k3))
cat("\n")

cat("ANOTACIÓN GÉNICA:\n")
cat("- genes filtrados:", nrow(annotation_df), "\n")
cat("- con símbolo HGNC:", sum(annotation_df$hgnc_symbol != ""), "\n")
cat("- sin símbolo HGNC:", sum(annotation_df$hgnc_symbol == ""), "\n\n")

cat("COEFICIENTES USADOS PARA APEGLM:\n")
cat("- cluster 1 vs 2 -> cluster_k3_1_vs_2 (ref=2)\n")
cat("- cluster 1 vs 3 -> cluster_k3_1_vs_3 (ref=3)\n")
cat("- cluster 2 vs 3 -> cluster_k3_2_vs_3 (ref=3)\n\n")

cat("RESULTADOS POR COMPARACIÓN:\n\n")

for (x in list(res_1_vs_2, res_1_vs_3, res_2_vs_3)) {
  cat("Comparación:", x$summary$comparison, "\n")
  cat("- coeficiente:", x$summary$coefficient, "\n")
  cat("- genes evaluados:", x$summary$total_genes_tested, "\n")
  cat("- genes significativos (padj < 0.05 y |log2FC| >= 1):", x$summary$significant_genes, "\n")
  cat("- genes up:", x$summary$up_genes, "\n")
  cat("- genes down:", x$summary$down_genes, "\n\n")
}

cat("ARCHIVOS GENERADOS:\n")
cat("- gene_annotation_table.csv\n")
cat("- dds_cluster_de_ref2.rds\n")
cat("- dds_cluster_de_ref3.rds\n")
cat("- normalized_counts_cluster_de_annotated.csv\n")
cat("- results_cluster_1_vs_2.csv\n")
cat("- results_cluster_1_vs_3.csv\n")
cat("- results_cluster_2_vs_3.csv\n")
cat("- sig_cluster_1_vs_2.csv\n")
cat("- sig_cluster_1_vs_3.csv\n")
cat("- sig_cluster_2_vs_3.csv\n")
cat("- up_cluster_1_vs_2.csv\n")
cat("- up_cluster_1_vs_3.csv\n")
cat("- up_cluster_2_vs_3.csv\n")
cat("- down_cluster_1_vs_2.csv\n")
cat("- down_cluster_1_vs_3.csv\n")
cat("- down_cluster_2_vs_3.csv\n")
cat("- volcano_cluster_1_vs_2.png\n")
cat("- volcano_cluster_1_vs_3.png\n")
cat("- volcano_cluster_2_vs_3.png\n")
cat("- heatmap_topgenes_cluster_1_vs_2.png\n")
cat("- heatmap_topgenes_cluster_1_vs_3.png\n")
cat("- heatmap_topgenes_cluster_2_vs_3.png\n")
cat("- differential_expression_clusters_report.txt\n")

sink()

#---------------------------#
# 12) Mensaje final         #
#---------------------------#
cat("\n=== SCRIPT 07 COMPLETADO ===\n")
cat("\nArchivos generados en:\n")
cat(outdir, "\n")