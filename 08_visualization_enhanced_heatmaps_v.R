### 08_visualization_enhanced_heatmaps_volcano_PRIMARY_U18_TARGET_ALL_P2.R
### Versión corregida:
### - Volcano plots: no etiquetan ENSG sin símbolo HGNC válido
### - Heatmaps: no muestran ENSG como rownames
### - Genes sin símbolo HGNC se conservan en tablas, pero no se usan como etiquetas de figuras principales

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(DESeq2)
  library(pheatmap)
  library(RColorBrewer)
  library(EnhancedVolcano)
  library(grid)
})

#---------------------------#
# 1) Paths y directorios    #
#---------------------------#

indir <- "subsets/primary_u18/final_analysis/07_differential_expression_clusters_deseq2"
outdir <- "subsets/primary_u18/final_analysis/08_visualization_enhanced"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

input_dds_ref2 <- file.path(indir, "dds_cluster_de_ref2.rds")
input_dds_ref3 <- file.path(indir, "dds_cluster_de_ref3.rds")

input_res_1_vs_2 <- file.path(indir, "results_cluster_1_vs_2.csv")
input_res_1_vs_3 <- file.path(indir, "results_cluster_1_vs_3.csv")
input_res_2_vs_3 <- file.path(indir, "results_cluster_2_vs_3.csv")

cat("=== SCRIPT 08: ENHANCED VISUALIZATION OF DIFFERENTIAL EXPRESSION ===\n")

#---------------------------#
# 2) Parámetros             #
#---------------------------#

padj_cutoff <- 0.05
lfc_cutoff <- 1

n_label_up <- 20
n_label_down <- 20

heatmap_n_up <- 20
heatmap_n_down <- 20
heatmap_selection_mode <- "mixed"

heatmap_width <- 16
heatmap_height <- 9
volcano_width <- 12
volcano_height <- 9

show_sample_names <- FALSE
heatmap_scale <- "row"

volcano_xlim <- c(-5, 5)
volcano_point_size <- 2
volcano_lab_size <- 3.2
volcano_max_overlaps <- 200
volcano_connector_width <- 0.5

cluster_palette <- c(
  "1" = "#00BD00",
  "2" = "#FF6EB4",
  "3" = "#0000FF"
)

sample_type_palette_manual <- c(
  "Primary Blood Derived Cancer - Bone Marrow" = "#EEB422",
  "Primary Blood Derived Cancer - Peripheral Blood" = "#00B3EE"
)

heatmap_breaks <- seq(-2, 2, length.out = 201)
heatmap_legend_breaks <- c(-2, -1, 0, 1, 2)
heatmap_legend_labels <- c("-2", "-1", "0", "1", "2")

#---------------------------#
# 3) Cargar datos           #
#---------------------------#

cat("\n=== CARGA DE DATOS ===\n")

dds_ref2 <- readRDS(input_dds_ref2)
dds_ref3 <- readRDS(input_dds_ref3)

res_1_vs_2 <- read.csv(input_res_1_vs_2, stringsAsFactors = FALSE, check.names = FALSE)
res_1_vs_3 <- read.csv(input_res_1_vs_3, stringsAsFactors = FALSE, check.names = FALSE)
res_2_vs_3 <- read.csv(input_res_2_vs_3, stringsAsFactors = FALSE, check.names = FALSE)

cat("- dds_ref2 genes:", nrow(dds_ref2), " muestras:", ncol(dds_ref2), "\n")
cat("- dds_ref3 genes:", nrow(dds_ref3), " muestras:", ncol(dds_ref3), "\n")
cat("- res_1_vs_2 filas:", nrow(res_1_vs_2), "\n")
cat("- res_1_vs_3 filas:", nrow(res_1_vs_3), "\n")
cat("- res_2_vs_3 filas:", nrow(res_2_vs_3), "\n")

#---------------------------#
# 4) Funciones auxiliares   #
#---------------------------#

standardize_gene_symbol <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == "" | x == "NA"] <- NA
  x
}

has_valid_gene_symbol <- function(df) {
  !is.na(df$gene_symbol) &
    df$gene_symbol != "" &
    !grepl("^ENSG", df$gene_symbol)
}

make_unique_label_field <- function(res_df) {
  gene_symbol_clean <- standardize_gene_symbol(res_df$gene_symbol)

  duplicated_symbol <- duplicated(gene_symbol_clean) |
    duplicated(gene_symbol_clean, fromLast = TRUE)

  label_unique <- ifelse(
    duplicated_symbol | is.na(gene_symbol_clean),
    NA,
    gene_symbol_clean
  )

  return(label_unique)
}

prepare_results_for_plot <- function(res_df) {
  out <- res_df %>%
    mutate(
      gene_symbol = standardize_gene_symbol(gene_symbol),
      padj_plot = ifelse(is.na(padj), NA, pmax(padj, .Machine$double.xmin)),
      neglog10_padj = -log10(padj_plot),
      regulation = case_when(
        !is.na(padj) & padj < padj_cutoff & log2FoldChange >= lfc_cutoff  ~ "Up",
        !is.na(padj) & padj < padj_cutoff & log2FoldChange <= -lfc_cutoff ~ "Down",
        TRUE ~ "NS"
      )
    )

  out$label_unique <- make_unique_label_field(out)
  out
}

order_genes_by_mode <- function(df, mode = "mixed") {
  if (nrow(df) == 0) return(df)

  if (mode == "padj") {
    df %>% arrange(padj)
  } else if (mode == "lfc") {
    df %>% arrange(desc(abs(log2FoldChange)), padj)
  } else {
    df %>% arrange(padj, desc(abs(log2FoldChange)))
  }
}

select_heatmap_genes_balanced <- function(res_df, n_up = 20, n_down = 20, mode = "mixed") {

  res_df <- res_df %>%
    mutate(gene_symbol = standardize_gene_symbol(gene_symbol))

  sig_df <- res_df %>%
    filter(
      !is.na(padj),
      padj < padj_cutoff,
      abs(log2FoldChange) >= lfc_cutoff,
      has_valid_gene_symbol(.)
    )

  if (nrow(sig_df) == 0) return(sig_df)

  up_df <- sig_df %>%
    filter(log2FoldChange >= lfc_cutoff) %>%
    order_genes_by_mode(mode) %>%
    head(n_up)

  down_df <- sig_df %>%
    filter(log2FoldChange <= -lfc_cutoff) %>%
    order_genes_by_mode(mode) %>%
    head(n_down)

  bind_rows(up_df, down_df)
}

select_volcano_labels <- function(res_df, n_up = 20, n_down = 20) {

  res_df <- res_df %>%
    mutate(gene_symbol = standardize_gene_symbol(gene_symbol)) %>%
    filter(has_valid_gene_symbol(.))

  up_df <- res_df %>%
    filter(!is.na(padj), padj < padj_cutoff, log2FoldChange >= lfc_cutoff) %>%
    mutate(neglog10_padj = -log10(padj_plot)) %>%
    arrange(desc(abs(log2FoldChange)), desc(neglog10_padj)) %>%
    head(n_up)

  down_df <- res_df %>%
    filter(!is.na(padj), padj < padj_cutoff, log2FoldChange <= -lfc_cutoff) %>%
    mutate(neglog10_padj = -log10(padj_plot)) %>%
    arrange(desc(abs(log2FoldChange)), desc(neglog10_padj)) %>%
    head(n_down)

  bind_rows(up_df, down_df)
}

add_heatmap_colorbar_label <- function(label = "Row Z-score of VST expression",
                                       x = 0.82,
                                       y = 0.60,
                                       rot = 0,
                                       gp = grid::gpar(fontsize = 10.5, fontface = "plain")) {
  grid::grid.text(
    label = label,
    x = grid::unit(x, "npc"),
    y = grid::unit(y, "npc"),
    rot = rot,
    gp = gp
  )
}

#---------------------------#
# 5) Heatmap corregido      #
#---------------------------#

make_publication_heatmap <- function(dds_obj,
                                     res_df,
                                     comparison_label,
                                     cluster_a,
                                     cluster_b,
                                     outfile_png,
                                     outfile_gene_table,
                                     n_up = 20,
                                     n_down = 20,
                                     selection_mode = "mixed") {

  cat("\n=== HEATMAP:", comparison_label, "===\n")

  top_df <- select_heatmap_genes_balanced(
    res_df = res_df,
    n_up = n_up,
    n_down = n_down,
    mode = selection_mode
  )

  if (nrow(top_df) < 2) {
    cat("[INFO] No hay suficientes genes con símbolo HGNC válido para heatmap en", comparison_label, "\n")
    return(NULL)
  }

  write.csv(top_df, file = outfile_gene_table, row.names = FALSE)

  vsd <- vst(dds_obj, blind = FALSE)
  vst_mat <- assay(vsd)

  meta_df <- as.data.frame(colData(dds_obj))
  meta_df$cluster_k3 <- as.factor(meta_df$cluster_k3)

  selected_samples <- rownames(meta_df)[meta_df$cluster_k3 %in% c(cluster_a, cluster_b)]
  selected_samples <- intersect(selected_samples, colnames(vst_mat))

  selected_genes <- intersect(top_df$ensembl_gene_id_original, rownames(vst_mat))

  if (length(selected_genes) < 2) {
    cat("[INFO] No hay suficientes genes presentes en la matriz VST para", comparison_label, "\n")
    return(NULL)
  }

  mat <- vst_mat[selected_genes, selected_samples, drop = FALSE]

  gene_labels <- top_df$gene_symbol[
    match(rownames(mat), top_df$ensembl_gene_id_original)
  ]

  gene_labels <- standardize_gene_symbol(gene_labels)

  keep_label <- !is.na(gene_labels) & gene_labels != "" & !grepl("^ENSG", gene_labels)

  mat <- mat[keep_label, , drop = FALSE]
  gene_labels <- gene_labels[keep_label]

  if (nrow(mat) < 2) {
    cat("[INFO] Después de remover genes sin símbolo HGNC válido quedan <2 genes en", comparison_label, "\n")
    return(NULL)
  }

  rownames(mat) <- make.unique(gene_labels)

  annotation_col <- data.frame(row.names = selected_samples)

if ("cluster_k3" %in% colnames(meta_df)) {
  annotation_col$Clusters <- as.factor(meta_df[selected_samples, "cluster_k3"])
} else if ("cluster_k3" %in% names(colData(dds_obj))) {
  annotation_col$Clusters <- as.factor(colData(dds_obj)$cluster_k3[selected_samples])
} else {
  stop("No se encontró cluster_k3 en colData(dds_obj). Revisa dds_ref2/dds_ref3.")
}

if ("sample_type.x" %in% colnames(meta_df)) {
  annotation_col$`Sample type` <- meta_df[selected_samples, "sample_type.x"]
}
  ann_colors <- list()

  if ("Clusters" %in% colnames(annotation_col)) {
    cluster_levels <- sort(unique(as.character(annotation_col$Clusters)))
    ann_colors$Clusters <- cluster_palette[cluster_levels]
  }

  if ("Sample type" %in% colnames(annotation_col)) {
    sample_type_levels <- unique(as.character(annotation_col$`Sample type`))

    sample_type_palette <- sample_type_palette_manual
    missing_levels <- setdiff(sample_type_levels, names(sample_type_palette))

    if (length(missing_levels) > 0) {
      extra_cols <- colorRampPalette(c("#7E6148", "#B09C85", "#6A3D9A"))(length(missing_levels))
      names(extra_cols) <- missing_levels
      sample_type_palette <- c(sample_type_palette, extra_cols)
    }

    ann_colors$`Sample type` <- sample_type_palette[sample_type_levels]
  }

  png(filename = outfile_png, width = heatmap_width, height = heatmap_height, units = "in", res = 300)

  pheatmap(
    mat,
    scale = heatmap_scale,
    breaks = heatmap_breaks,
    color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(200),
    legend_breaks = heatmap_legend_breaks,
    legend_labels = heatmap_legend_labels,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    show_rownames = TRUE,
    show_colnames = show_sample_names,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_method = "complete",
    fontsize_row = 10,
    fontsize_col = 7,
    border_color = NA,
    main = paste0(
      "Heatmap de genes diferenciales: ", comparison_label,
      "\nTop genes con símbolo HGNC válido (", n_up, " + ", n_down, "; ", selection_mode, ")"
    )
  )

  add_heatmap_colorbar_label(
    label = "Row Z-score of VST expression",
    x = 0.82,
    y = 0.60,
    rot = 0,
    gp = grid::gpar(fontsize = 10.5)
  )

  dev.off()

  cat("[OK] Heatmap generado:", outfile_png, "\n")
}

#---------------------------#
# 6) Volcano corregido      #
#---------------------------#

make_enhanced_volcano <- function(res_df,
                                  comparison_label,
                                  outfile_png,
                                  outfile_label_table,
                                  n_up = 20,
                                  n_down = 20) {

  cat("\n=== ENHANCED VOLCANO:", comparison_label, "===\n")

  plot_df <- prepare_results_for_plot(res_df)

  label_df <- select_volcano_labels(
    res_df = plot_df,
    n_up = n_up,
    n_down = n_down
  )

  write.csv(label_df, file = outfile_label_table, row.names = FALSE)

  max_y_value <- max(-log10(plot_df$padj_plot), na.rm = TRUE) + 2

  png(filename = outfile_png, width = volcano_width, height = volcano_height, units = "in", res = 300)

  print(
    EnhancedVolcano(
      plot_df,
      lab = plot_df$label_unique,
      x = "log2FoldChange",
      y = "padj_plot",
      selectLab = unique(na.omit(label_df$label_unique)),
      xlim = volcano_xlim,
      ylim = c(0, max_y_value),
      title = paste0("Volcano Plot: ", comparison_label),
      subtitle = paste0(
        "Etiquetados: ", n_up, " sobreexpresados + ", n_down,
        " subexpresados con símbolo HGNC válido"
      ),
      pCutoff = padj_cutoff,
      FCcutoff = lfc_cutoff,
      pointSize = volcano_point_size,
      labSize = volcano_lab_size,
      boxedLabels = TRUE,
      drawConnectors = TRUE,
      widthConnectors = volcano_connector_width,
      max.overlaps = volcano_max_overlaps,
      col = c("grey30", "forestgreen", "royalblue", "red2"),
      legendLabels = c(
        "NS",
        paste0("|Log2 FC| >", lfc_cutoff),
        paste0("FDR < ", padj_cutoff),
        paste0("FDR < ", padj_cutoff, " & |log2FC| > ", lfc_cutoff)
      )
    ) +
      scale_x_continuous(
        labels = function(x) gsub("-", "\u2212", x)
      )
  )

  dev.off()

  cat("[OK] Enhanced volcano generado:", outfile_png, "\n")
}

#---------------------------#
# 7) Ejecutar pares         #
#---------------------------#

run_visualization_pair <- function(dds_obj,
                                   res_df,
                                   comparison_label,
                                   cluster_a,
                                   cluster_b,
                                   file_stub) {

  make_publication_heatmap(
    dds_obj = dds_obj,
    res_df = res_df,
    comparison_label = comparison_label,
    cluster_a = cluster_a,
    cluster_b = cluster_b,
    outfile_png = file.path(outdir, paste0("heatmap_publication_", file_stub, ".png")),
    outfile_gene_table = file.path(outdir, paste0("selected_heatmap_genes_", file_stub, ".csv")),
    n_up = heatmap_n_up,
    n_down = heatmap_n_down,
    selection_mode = heatmap_selection_mode
  )

  make_enhanced_volcano(
    res_df = res_df,
    comparison_label = comparison_label,
    outfile_png = file.path(outdir, paste0("enhanced_volcano_", file_stub, ".png")),
    outfile_label_table = file.path(outdir, paste0("selected_volcano_labels_", file_stub, ".csv")),
    n_up = n_label_up,
    n_down = n_label_down
  )
}

cat("\n=== GENERANDO VISUALIZACIONES ===\n")

run_visualization_pair(
  dds_obj = dds_ref2,
  res_df = res_1_vs_2,
  comparison_label = "cluster 1 vs 2",
  cluster_a = "1",
  cluster_b = "2",
  file_stub = "cluster_1_vs_2"
)

run_visualization_pair(
  dds_obj = dds_ref3,
  res_df = res_1_vs_3,
  comparison_label = "cluster 1 vs 3",
  cluster_a = "1",
  cluster_b = "3",
  file_stub = "cluster_1_vs_3"
)

run_visualization_pair(
  dds_obj = dds_ref3,
  res_df = res_2_vs_3,
  comparison_label = "cluster 2 vs 3",
  cluster_a = "2",
  cluster_b = "3",
  file_stub = "cluster_2_vs_3"
)

#---------------------------#
# 8) Reporte TXT            #
#---------------------------#

cat("\n=== GENERANDO REPORTE TXT ===\n")

report_file <- file.path(outdir, "visualization_report.txt")
sink(report_file)

cat("REPORTE DE VISUALIZACIÓN AVANZADA\n")
cat("Fecha:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("INPUT DIRECTORY:\n")
cat(indir, "\n\n")

cat("OUTPUT DIRECTORY:\n")
cat(outdir, "\n\n")

cat("PARÁMETROS GENERALES:\n")
cat("- padj_cutoff:", padj_cutoff, "\n")
cat("- lfc_cutoff:", lfc_cutoff, "\n")
cat("- n_label_up:", n_label_up, "\n")
cat("- n_label_down:", n_label_down, "\n")
cat("- heatmap_n_up:", heatmap_n_up, "\n")
cat("- heatmap_n_down:", heatmap_n_down, "\n")
cat("- heatmap_selection_mode:", heatmap_selection_mode, "\n")
cat("- heatmap_scale:", heatmap_scale, "\n")
cat("- volcano_xlim:", paste(volcano_xlim, collapse = " to "), "\n\n")

cat("CORRECCIÓN DE ANOTACIÓN PARA FIGURAS:\n")
cat("- Los resultados DESeq2 conservan todos los genes evaluados, incluyendo loci sin símbolo HGNC.\n")
cat("- Para las figuras principales, las etiquetas se restringieron a genes con símbolo HGNC válido.\n")
cat("- No se usó el identificador Ensembl original como fallback en volcano plots ni en heatmaps.\n")
cat("- Esta decisión evita mostrar identificadores Ensembl versionados, obsoletos o no mapeables en la anotación actual.\n")
cat("- Los Ensembl IDs sin símbolo no fueron eliminados del análisis estadístico; solo se excluyeron de las etiquetas visuales principales.\n\n")

cat("FUNDAMENTO METODOLÓGICO:\n")
cat("- Los archivos STAR counts de GDC pueden contener Ensembl gene IDs versionados, derivados de la anotación usada por el pipeline original.\n")
cat("- Algunos IDs Ensembl pueden no estar presentes en releases actuales porque los modelos génicos cambian entre versiones.\n")
cat("- Para evitar interpretaciones ambiguas en figuras de publicación, se priorizaron símbolos HGNC válidos.\n")
cat("- Los genes no mapeados permanecen disponibles en las tablas completas y pueden reportarse en material suplementario.\n")
cat("- Esta estrategia preserva la validez estadística del análisis DESeq2 y mejora la interpretabilidad biológica de las figuras.\n\n")

cat("VOLCANO LABEL CRITERIA:\n")
cat("- Primero: mayor |log2FoldChange|\n")
cat("- Segundo: mayor -log10(padj)\n")
cat("- Solo se etiquetan genes con gene_symbol válido.\n")
cat("- Genes duplicados por símbolo no se etiquetan para evitar asignaciones ambiguas.\n\n")

cat("HEATMAP GENE SELECTION:\n")
cat("- Se seleccionan genes significativos con padj < 0.05 y |log2FC| >= 1.\n")
cat("- La selección es balanceada: hasta 20 genes sobreexpresados y 20 subexpresados.\n")
cat("- Solo se incluyen genes con símbolo HGNC válido.\n")
cat("- No se reemplazan símbolos ausentes por Ensembl IDs.\n")
cat("- Por tanto, los heatmaps corregidos no deben mostrar ENSG como nombres de fila.\n\n")

cat("VOLCANO PADJ CORRECTION:\n")
cat("- Los padj == 0 no se convierten a 1.\n")
cat("- Se usa pmax(padj, .Machine$double.xmin) para evitar valores infinitos en -log10(padj).\n\n")

cat("HEATMAP COLOR BAR LABEL:\n")
cat("- Row Z-score of VST expression\n\n")

cat("ARCHIVOS GENERADOS:\n")
cat("- heatmap_publication_cluster_1_vs_2.png\n")
cat("- heatmap_publication_cluster_1_vs_3.png\n")
cat("- heatmap_publication_cluster_2_vs_3.png\n")
cat("- enhanced_volcano_cluster_1_vs_2.png\n")
cat("- enhanced_volcano_cluster_1_vs_3.png\n")
cat("- enhanced_volcano_cluster_2_vs_3.png\n")
cat("- selected_heatmap_genes_cluster_1_vs_2.csv\n")
cat("- selected_heatmap_genes_cluster_1_vs_3.csv\n")
cat("- selected_heatmap_genes_cluster_2_vs_3.csv\n")
cat("- selected_volcano_labels_cluster_1_vs_2.csv\n")
cat("- selected_volcano_labels_cluster_1_vs_3.csv\n")
cat("- selected_volcano_labels_cluster_2_vs_3.csv\n")
cat("- visualization_report.txt\n")

sink()

cat("\n=== SCRIPT 08 COMPLETADO ===\n")
cat("\nArchivos generados en:\n")
cat(outdir, "\n")