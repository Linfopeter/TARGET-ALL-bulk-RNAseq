### 09_functional_enrichment_PRIMARY_U18_TARGET_ALL_P2.R
### GOAL:
###   1) Realizar enriquecimiento funcional robusto para comparaciones entre clusters
###   2) Priorizar GSEA pre-ranked sobre GO/KEGG ORA para evitar sesgo por umbrales arbitrarios
###   3) Ejecutar ORA GO/KEGG como análisis complementario para genes significativos Up/Down
###   4) Generar tablas, dotplots y reportes por comparación
###
### INPUT:
###   - subsets/primary_u18/final_analysis/07_differential_expression_clusters_deseq2/results_cluster_1_vs_2.csv
###   - subsets/primary_u18/final_analysis/07_differential_expression_clusters_deseq2/results_cluster_1_vs_3.csv
###   - subsets/primary_u18/final_analysis/07_differential_expression_clusters_deseq2/results_cluster_2_vs_3.csv
###   - subsets/primary_u18/final_analysis/07_differential_expression_clusters_deseq2/gene_annotation_table.csv
###
### OUTPUT:
###   - subsets/primary_u18/final_analysis/09_functional_enrichment/
###       * gsea_GO_BP_*.csv / .png
###       * gsea_GO_MF_*.csv / .png
###       * gsea_GO_CC_*.csv / .png
###       * gsea_KEGG_*.csv / .png
###       * ora_GO_BP_up/down_*.csv / .png
###       * ora_KEGG_up/down_*.csv / .png
###       * functional_enrichment_report.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(DOSE)
})

cat("=== SCRIPT 09: FUNCTIONAL ENRICHMENT BETWEEN CLUSTERS ===\n")

#---------------------------#
# 1) Paths y parámetros     #
#---------------------------#

indir <- "subsets/primary_u18/final_analysis/07_differential_expression_clusters_deseq2"
outdir <- "subsets/primary_u18/final_analysis/09_functional_enrichment"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

input_annotation <- file.path(indir, "gene_annotation_table.csv")

comparison_files <- list(
  cluster_1_vs_2 = file.path(indir, "results_cluster_1_vs_2.csv"),
  cluster_1_vs_3 = file.path(indir, "results_cluster_1_vs_3.csv"),
  cluster_2_vs_3 = file.path(indir, "results_cluster_2_vs_3.csv")
)

padj_cutoff <- 0.05
lfc_cutoff <- 1
min_gene_set_size <- 10
max_gene_set_size <- 500
gsea_pvalue_cutoff <- 1
ora_pvalue_cutoff <- 1
report_top_n <- 20
plot_top_n <- 20

#---------------------------#
# 2) Utilidades             #
#---------------------------#

safe_write_csv <- function(x, file) {
  if (is.null(x)) return(invisible(NULL))
  readr::write_csv(as.data.frame(x), file)
}

clean_gene_symbol <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == "" | x == "NA" | grepl("^ENSG", x)] <- NA
  x
}

make_safe_stub <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", x)
}

save_dotplot <- function(enrich_obj, outfile, title, top_n = 20, width = 10, height = 7) {
  df <- as.data.frame(enrich_obj)
  if (is.null(df) || nrow(df) == 0) {
    cat("[INFO] Sin términos para graficar:", title, "\n")
    return(NULL)
  }

  p <- dotplot(enrich_obj, showCategory = min(top_n, nrow(df))) +
    ggtitle(title) +
    theme_bw(base_size = 12)

  ggsave(outfile, plot = p, width = width, height = height, dpi = 300)
}

save_gseaplot_top_terms <- function(gsea_obj, outfile_prefix, top_n = 5) {
  df <- as.data.frame(gsea_obj)
  if (is.null(df) || nrow(df) == 0) return(NULL)

  df <- df %>%
    dplyr::arrange(p.adjust) %>%
    utils::head(top_n)

  for (i in seq_len(nrow(df))) {
    term_id <- df$ID[i]
    safe_id <- make_safe_stub(term_id)
    outfile <- paste0(outfile_prefix, "_", i, "_", safe_id, ".png")

    p <- tryCatch({
      gseaplot2(gsea_obj, geneSetID = term_id, title = df$Description[i])
    }, error = function(e) {
      cat("[WARN] No se pudo generar gseaplot para", term_id, ":", conditionMessage(e), "\n")
      NULL
    })

    if (!is.null(p)) {
      ggsave(outfile, plot = p, width = 10, height = 6, dpi = 300)
    }
  }
}

# Ranking recomendado para GSEA:
# - signo: dirección del cambio del cluster A vs cluster B
# - magnitud: combinación de tamaño de efecto y significancia
# - evita usar solo genes significativos, preserva señal distribuida en vías
make_ranked_gene_list <- function(res_df) {
  rank_df <- res_df %>%
    dplyr::mutate(
      gene_symbol = clean_gene_symbol(gene_symbol),
      padj_rank = ifelse(is.na(padj), 1, pmax(padj, .Machine$double.xmin)),
      rank_metric = sign(log2FoldChange) * (-log10(padj_rank)) * abs(log2FoldChange)
    ) %>%
    dplyr::filter(
      !is.na(gene_symbol),
      !is.na(rank_metric),
      is.finite(rank_metric)
    ) %>%
    dplyr::group_by(gene_symbol) %>%
    dplyr::arrange(desc(abs(rank_metric)), .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(rank_metric))

  gene_list <- rank_df$rank_metric
  names(gene_list) <- rank_df$gene_symbol
  gene_list <- sort(gene_list, decreasing = TRUE)
  return(gene_list)
}

symbol_to_entrez <- function(symbols) {
  symbols <- unique(na.omit(symbols))
  if (length(symbols) == 0) return(character(0))

  mapped <- suppressMessages(
    bitr(
      symbols,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
  )

  mapped <- mapped %>%
    dplyr::filter(!is.na(ENTREZID)) %>%
    dplyr::distinct(SYMBOL, .keep_all = TRUE)

  return(mapped)
}

make_entrez_ranked_list <- function(gene_list_symbol) {
  map_df <- symbol_to_entrez(names(gene_list_symbol))
  if (nrow(map_df) == 0) return(numeric(0))

  rank_df <- tibble(
    SYMBOL = names(gene_list_symbol),
    rank_metric = as.numeric(gene_list_symbol)
  ) %>%
    dplyr::inner_join(map_df, by = "SYMBOL") %>%
    dplyr::group_by(ENTREZID) %>%
    dplyr::arrange(desc(abs(rank_metric)), .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(rank_metric))

  gene_list <- rank_df$rank_metric
  names(gene_list) <- rank_df$ENTREZID
  gene_list <- sort(gene_list, decreasing = TRUE)
  return(gene_list)
}

run_gsea_go <- function(gene_list_symbol, ont, comp_stub) {
  cat("[GSEA GO", ont, "]", comp_stub, "\n")

  if (length(gene_list_symbol) < min_gene_set_size) return(NULL)

  gsea <- tryCatch({
    gseGO(
      geneList = gene_list_symbol,
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = ont,
      minGSSize = min_gene_set_size,
      maxGSSize = max_gene_set_size,
      pvalueCutoff = gsea_pvalue_cutoff,
      pAdjustMethod = "BH",
      verbose = FALSE,
      seed = TRUE
    )
  }, error = function(e) {
    cat("[WARN] Falló gseGO", ont, comp_stub, ":", conditionMessage(e), "\n")
    NULL
  })

  if (is.null(gsea)) return(NULL)

  out_csv <- file.path(outdir, paste0("gsea_GO_", ont, "_", comp_stub, ".csv"))
  out_png <- file.path(outdir, paste0("gsea_GO_", ont, "_dotplot_", comp_stub, ".png"))
  safe_write_csv(as.data.frame(gsea), out_csv)
  save_dotplot(gsea, out_png, paste0("GSEA GO ", ont, ": ", comp_stub), plot_top_n)
  save_gseaplot_top_terms(gsea, file.path(outdir, paste0("gseaplot_GO_", ont, "_", comp_stub)), top_n = 5)

  return(gsea)
}

run_gsea_kegg <- function(gene_list_entrez, comp_stub) {
  cat("[GSEA KEGG]", comp_stub, "\n")

  if (length(gene_list_entrez) < min_gene_set_size) return(NULL)

  gsea <- tryCatch({
    gseKEGG(
      geneList = gene_list_entrez,
      organism = "hsa",
      minGSSize = min_gene_set_size,
      maxGSSize = max_gene_set_size,
      pvalueCutoff = gsea_pvalue_cutoff,
      pAdjustMethod = "BH",
      verbose = FALSE,
      seed = TRUE
    )
  }, error = function(e) {
    cat("[WARN] Falló gseKEGG", comp_stub, ":", conditionMessage(e), "\n")
    NULL
  })

  if (is.null(gsea)) return(NULL)

  out_csv <- file.path(outdir, paste0("gsea_KEGG_", comp_stub, ".csv"))
  out_png <- file.path(outdir, paste0("gsea_KEGG_dotplot_", comp_stub, ".png"))
  safe_write_csv(as.data.frame(gsea), out_csv)
  save_dotplot(gsea, out_png, paste0("GSEA KEGG: ", comp_stub), plot_top_n)
  save_gseaplot_top_terms(gsea, file.path(outdir, paste0("gseaplot_KEGG_", comp_stub)), top_n = 5)

  return(gsea)
}

run_ora_go <- function(gene_symbols, universe_symbols, ont, direction, comp_stub) {
  cat("[ORA GO", ont, direction, "]", comp_stub, "\n")

  gene_symbols <- unique(na.omit(gene_symbols))
  universe_symbols <- unique(na.omit(universe_symbols))

  if (length(gene_symbols) < 5) {
    cat("[INFO] Muy pocos genes para ORA GO", ont, direction, comp_stub, "\n")
    return(NULL)
  }

  ora <- tryCatch({
    enrichGO(
      gene = gene_symbols,
      universe = universe_symbols,
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = ont,
      pAdjustMethod = "BH",
      pvalueCutoff = ora_pvalue_cutoff,
      qvalueCutoff = ora_pvalue_cutoff,
      readable = TRUE
    )
  }, error = function(e) {
    cat("[WARN] Falló enrichGO", ont, direction, comp_stub, ":", conditionMessage(e), "\n")
    NULL
  })

  if (is.null(ora)) return(NULL)

  out_csv <- file.path(outdir, paste0("ora_GO_", ont, "_", direction, "_", comp_stub, ".csv"))
  out_png <- file.path(outdir, paste0("ora_GO_", ont, "_", direction, "_dotplot_", comp_stub, ".png"))
  safe_write_csv(as.data.frame(ora), out_csv)
  save_dotplot(ora, out_png, paste0("ORA GO ", ont, " ", direction, ": ", comp_stub), plot_top_n)

  return(ora)
}

run_ora_kegg <- function(gene_symbols, universe_symbols, direction, comp_stub) {
  cat("[ORA KEGG", direction, "]", comp_stub, "\n")

  gene_map <- symbol_to_entrez(gene_symbols)
  universe_map <- symbol_to_entrez(universe_symbols)

  if (nrow(gene_map) < 5 || nrow(universe_map) < 10) {
    cat("[INFO] Muy pocos genes mapeados para ORA KEGG", direction, comp_stub, "\n")
    return(NULL)
  }

  ora <- tryCatch({
    enrichKEGG(
      gene = unique(gene_map$ENTREZID),
      universe = unique(universe_map$ENTREZID),
      organism = "hsa",
      pAdjustMethod = "BH",
      pvalueCutoff = ora_pvalue_cutoff,
      qvalueCutoff = ora_pvalue_cutoff
    )
  }, error = function(e) {
    cat("[WARN] Falló enrichKEGG", direction, comp_stub, ":", conditionMessage(e), "\n")
    NULL
  })

  if (is.null(ora)) return(NULL)

  ora <- tryCatch({
    setReadable(ora, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  }, error = function(e) ora)

  out_csv <- file.path(outdir, paste0("ora_KEGG_", direction, "_", comp_stub, ".csv"))
  out_png <- file.path(outdir, paste0("ora_KEGG_", direction, "_dotplot_", comp_stub, ".png"))
  safe_write_csv(as.data.frame(ora), out_csv)
  save_dotplot(ora, out_png, paste0("ORA KEGG ", direction, ": ", comp_stub), plot_top_n)

  return(ora)
}

#---------------------------#
# 3) Carga y validación     #
#---------------------------#

cat("\n=== CARGA DE DATOS ===\n")

if (!file.exists(input_annotation)) {
  stop("No existe gene_annotation_table.csv en: ", input_annotation)
}

annotation_df <- read.csv(input_annotation, stringsAsFactors = FALSE, check.names = FALSE)

for (nm in names(comparison_files)) {
  if (!file.exists(comparison_files[[nm]])) {
    stop("No existe archivo de resultados para ", nm, ": ", comparison_files[[nm]])
  }
}

required_cols <- c("ensembl_gene_id_original", "gene_symbol", "log2FoldChange", "padj")

#---------------------------#
# 4) Ejecutar análisis      #
#---------------------------#

all_summaries <- list()

for (comp_stub in names(comparison_files)) {
  cat("\n========================================\n")
  cat("PROCESANDO COMPARACIÓN:", comp_stub, "\n")
  cat("========================================\n")

  res_df <- read.csv(comparison_files[[comp_stub]], stringsAsFactors = FALSE, check.names = FALSE)

  missing_cols <- setdiff(required_cols, colnames(res_df))
  if (length(missing_cols) > 0) {
    stop("Faltan columnas en ", comp_stub, ": ", paste(missing_cols, collapse = ", "))
  }

  res_df <- res_df %>%
    dplyr::mutate(
      gene_symbol = clean_gene_symbol(gene_symbol),
      direction = case_when(
        !is.na(padj) & padj < padj_cutoff & log2FoldChange >= lfc_cutoff ~ "up",
        !is.na(padj) & padj < padj_cutoff & log2FoldChange <= -lfc_cutoff ~ "down",
        TRUE ~ "ns"
      )
    )

  universe_symbols <- res_df %>%
    dplyr::filter(!is.na(gene_symbol)) %>%
    dplyr::pull(gene_symbol) %>%
    unique()

  up_symbols <- res_df %>%
    dplyr::filter(direction == "up", !is.na(gene_symbol)) %>%
    dplyr::pull(gene_symbol) %>%
    unique()

  down_symbols <- res_df %>%
    dplyr::filter(direction == "down", !is.na(gene_symbol)) %>%
    dplyr::pull(gene_symbol) %>%
    unique()

  gene_list_symbol <- make_ranked_gene_list(res_df)
  gene_list_entrez <- make_entrez_ranked_list(gene_list_symbol)

  rank_table <- tibble(
    gene_symbol = names(gene_list_symbol),
    rank_metric = as.numeric(gene_list_symbol)
  )
  write_csv(rank_table, file.path(outdir, paste0("gsea_ranked_gene_list_SYMBOL_", comp_stub, ".csv")))

  rank_table_entrez <- tibble(
    entrez_id = names(gene_list_entrez),
    rank_metric = as.numeric(gene_list_entrez)
  )
  write_csv(rank_table_entrez, file.path(outdir, paste0("gsea_ranked_gene_list_ENTREZ_", comp_stub, ".csv")))

  gsea_go_bp <- run_gsea_go(gene_list_symbol, "BP", comp_stub)
  gsea_go_mf <- run_gsea_go(gene_list_symbol, "MF", comp_stub)
  gsea_go_cc <- run_gsea_go(gene_list_symbol, "CC", comp_stub)
  gsea_kegg <- run_gsea_kegg(gene_list_entrez, comp_stub)

  ora_go_bp_up <- run_ora_go(up_symbols, universe_symbols, "BP", "up", comp_stub)
  ora_go_bp_down <- run_ora_go(down_symbols, universe_symbols, "BP", "down", comp_stub)

  ora_kegg_up <- run_ora_kegg(up_symbols, universe_symbols, "up", comp_stub)
  ora_kegg_down <- run_ora_kegg(down_symbols, universe_symbols, "down", comp_stub)

  all_summaries[[comp_stub]] <- data.frame(
    comparison = comp_stub,
    genes_in_results = nrow(res_df),
    genes_with_valid_symbol = length(universe_symbols),
    significant_up_symbols = length(up_symbols),
    significant_down_symbols = length(down_symbols),
    ranked_symbols_for_gsea = length(gene_list_symbol),
    ranked_entrez_for_kegg = length(gene_list_entrez),
    gsea_go_bp_terms = ifelse(is.null(gsea_go_bp), 0, nrow(as.data.frame(gsea_go_bp))),
    gsea_go_mf_terms = ifelse(is.null(gsea_go_mf), 0, nrow(as.data.frame(gsea_go_mf))),
    gsea_go_cc_terms = ifelse(is.null(gsea_go_cc), 0, nrow(as.data.frame(gsea_go_cc))),
    gsea_kegg_terms = ifelse(is.null(gsea_kegg), 0, nrow(as.data.frame(gsea_kegg))),
    ora_go_bp_up_terms = ifelse(is.null(ora_go_bp_up), 0, nrow(as.data.frame(ora_go_bp_up))),
    ora_go_bp_down_terms = ifelse(is.null(ora_go_bp_down), 0, nrow(as.data.frame(ora_go_bp_down))),
    ora_kegg_up_terms = ifelse(is.null(ora_kegg_up), 0, nrow(as.data.frame(ora_kegg_up))),
    ora_kegg_down_terms = ifelse(is.null(ora_kegg_down), 0, nrow(as.data.frame(ora_kegg_down)))
  )
}

summary_df <- dplyr::bind_rows(all_summaries)
write_csv(summary_df, file.path(outdir, "functional_enrichment_summary.csv"))

#---------------------------#
# 5) Reporte                #
#---------------------------#

report_file <- file.path(outdir, "functional_enrichment_report.txt")
sink(report_file)

cat("REPORTE DE ENRIQUECIMIENTO FUNCIONAL\n")
cat("Fecha:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("INPUT DIRECTORY:\n")
cat(indir, "\n\n")

cat("OUTPUT DIRECTORY:\n")
cat(outdir, "\n\n")

cat("DECISIÓN METODOLÓGICA PRINCIPAL:\n")
cat("- Análisis primario recomendado: GSEA pre-ranked usando todos los genes evaluados con símbolo HGNC válido.\n")
cat("- Análisis complementario: ORA GO/KEGG separado para genes Up y Down significativos.\n")
cat("- Motivo: los clusters fueron definidos por transcriptoma; usar solo listas significativas puede amplificar sesgo por umbral.\n")
cat("- GSEA conserva la estructura continua del ranking de genes y detecta cambios coordinados moderados.\n")
cat("- KEGG se incluye por interpretabilidad de vías, pero no debe ser la única fuente por cobertura limitada.\n")
cat("- GO BP se considera el resultado biológico más amplio; GO MF/CC son complementarios.\n\n")

cat("PARÁMETROS:\n")
cat("- padj_cutoff:", padj_cutoff, "\n")
cat("- lfc_cutoff:", lfc_cutoff, "\n")
cat("- min_gene_set_size:", min_gene_set_size, "\n")
cat("- max_gene_set_size:", max_gene_set_size, "\n")
cat("- pAdjustMethod: BH\n")
cat("- GSEA rank metric: sign(log2FC) * -log10(padj) * abs(log2FC)\n\n")

cat("RESUMEN POR COMPARACIÓN:\n")
print(summary_df)
cat("\n")

cat("INTERPRETACIÓN DE DIRECCIÓN EN GSEA:\n")
cat("- NES positivo: vía enriquecida hacia el primer cluster de la comparación.\n")
cat("- NES negativo: vía enriquecida hacia el segundo cluster de la comparación.\n")
cat("Ejemplo: cluster_1_vs_2 con NES positivo implica enriquecimiento hacia cluster 1; NES negativo hacia cluster 2.\n\n")

cat("ARCHIVOS PRINCIPALES GENERADOS:\n")
cat("- functional_enrichment_summary.csv\n")
cat("- gsea_ranked_gene_list_SYMBOL_*.csv\n")
cat("- gsea_ranked_gene_list_ENTREZ_*.csv\n")
cat("- gsea_GO_BP_*.csv y dotplots\n")
cat("- gsea_GO_MF_*.csv y dotplots\n")
cat("- gsea_GO_CC_*.csv y dotplots\n")
cat("- gsea_KEGG_*.csv y dotplots\n")
cat("- gseaplot_GO_* para términos top\n")
cat("- gseaplot_KEGG_* para términos top\n")
cat("- ora_GO_BP_up/down_*.csv y dotplots\n")
cat("- ora_KEGG_up/down_*.csv y dotplots\n")
cat("- functional_enrichment_report.txt\n")

sink()

cat("\n=== SCRIPT 09 COMPLETADO ===\n")
cat("Archivos generados en:\n")
cat(outdir, "\n")
