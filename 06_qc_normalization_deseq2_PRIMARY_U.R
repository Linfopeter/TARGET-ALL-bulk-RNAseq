### 06_qc_normalization_deseq2_PRIMARY_U18_TARGET_ALL_P2.R
### GOAL:
###   1) Limpiar la matriz de counts eliminando filas técnicas de STAR
###   2) Normalizar con DESeq2 y generar PCA + heatmap de genes variables
###
### INPUT:
###   - subsets/primary_u18/final_analysis/TARGET_ALL_P2_PRIMARY_U18_FINAL_FOR_ANALYSIS.RData
###
### OUTPUT:
###   - subsets/primary_u18/final_analysis/06_qc_normalization_deseq2/
###       * counts_gene_only_raw.RData
###       * dds_filtered.rds
###       * vst_object.rds
###       * vst_matrix.csv
###       * sample_qc_metadata.csv
###       * pca_samples_deseq2.png
###       * heatmap_top_variable_genes_deseq2.png
###       * top_variable_genes_deseq2.csv
###       * qc_normalization_deseq2_report.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(S4Vectors)
  library(uwot)
})

#---------------------------#
# 1) Paths y directorios    #
#---------------------------#
input_file <- "subsets/primary_u18/final_analysis/TARGET_ALL_P2_PRIMARY_U18_FINAL_FOR_ANALYSIS.RData"
outdir <- "subsets/primary_u18/final_analysis/06_qc_normalization_deseq2"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("=== SCRIPT 06: QC + NORMALIZATION WITH DESEQ2 ===\n")

#---------------------------#
# 2) Cargar dataset final   #
#---------------------------#
load(input_file)
# Esperado:
#   - counts_final
#   - meta_final

cat("\nDimensiones iniciales:\n")
cat("- Genes/filas en counts_final:", nrow(counts_final), "\n")
cat("- Muestras/columnas en counts_final:", ncol(counts_final), "\n")
cat("- Filas en meta_final:", nrow(meta_final), "\n")
cat("- Columnas en meta_final:", ncol(meta_final), "\n")

cat("\nChequeo de alineación inicial:\n")
cat("¿colnames(counts_final) == meta_final$counts_colname?: ",
    all(colnames(counts_final) == meta_final$counts_colname), "\n")

#---------------------------#
# 3) Preparar metadata plana#
#---------------------------#
cat("\n=== PREPARACIÓN DE METADATA PLANA ===\n")

list_cols <- names(meta_final)[sapply(meta_final, is.list)]

cat("\nColumnas tipo list detectadas en meta_final:\n")
print(list_cols)

meta_flat <- meta_final[, !sapply(meta_final, is.list), drop = FALSE]

cat("\nDimensiones de meta_flat:\n")
cat("- Filas:", nrow(meta_flat), "\n")
cat("- Columnas:", ncol(meta_flat), "\n")

#---------------------------#
# 4) Remover filas STAR     #
#---------------------------#
cat("\n=== LIMPIEZA DE FILAS TÉCNICAS DE STAR ===\n")

star_technical_rows <- c(
  "N_unmapped",
  "N_multimapping",
  "N_noFeature",
  "N_ambiguous"
)

present_star_rows <- intersect(rownames(counts_final), star_technical_rows)
missing_star_rows <- setdiff(star_technical_rows, rownames(counts_final))

cat("\nFilas técnicas encontradas:\n")
print(present_star_rows)

cat("\nFilas técnicas no encontradas:\n")
print(missing_star_rows)

technical_counts_summary <- NULL
if (length(present_star_rows) > 0) {
  technical_counts_summary <- counts_final[present_star_rows, , drop = FALSE]
}

counts_gene_only <- counts_final[!rownames(counts_final) %in% star_technical_rows, , drop = FALSE]

cat("\nDimensiones después de remover filas técnicas:\n")
cat("- Filas antes:", nrow(counts_final), "\n")
cat("- Filas después:", nrow(counts_gene_only), "\n")
cat("- Columnas después:", ncol(counts_gene_only), "\n")

#---------------------------#
# 5) QC básico por muestra  #
#---------------------------#
cat("\n=== QC BÁSICO POR MUESTRA ===\n")

library_size_raw <- colSums(counts_final, na.rm = TRUE)
library_size_gene_only <- colSums(counts_gene_only, na.rm = TRUE)

if (!is.null(technical_counts_summary)) {
  technical_total <- colSums(technical_counts_summary, na.rm = TRUE)
} else {
  technical_total <- rep(0, ncol(counts_final))
  names(technical_total) <- colnames(counts_final)
}

pct_technical <- 100 * technical_total / library_size_raw
genes_detected_gt0 <- colSums(counts_gene_only > 0, na.rm = TRUE)

sample_qc <- meta_flat %>%
  mutate(
    library_size_raw = library_size_raw[counts_colname],
    library_size_gene_only = library_size_gene_only[counts_colname],
    technical_counts = technical_total[counts_colname],
    pct_technical = pct_technical[counts_colname],
    genes_detected_gt0 = genes_detected_gt0[counts_colname]
  )

cat("\nResumen de library_size_raw:\n")
print(summary(sample_qc$library_size_raw))

cat("\nResumen de pct_technical:\n")
print(summary(sample_qc$pct_technical))

cat("\nResumen de genes_detected_gt0:\n")
print(summary(sample_qc$genes_detected_gt0))

write.csv(
  sample_qc,
  file = file.path(outdir, "sample_qc_metadata.csv"),
  row.names = FALSE
)

#---------------------------#
# 6) Preparar matriz counts #
#---------------------------#
cat("\n=== PREPARACIÓN DE MATRIZ DE COUNTS ===\n")

if (!all(abs(counts_gene_only - round(counts_gene_only)) < .Machine$double.eps^0.5)) {
  stop("La matriz counts_gene_only contiene valores no enteros. DESeq2 requiere counts enteros.")
}

counts_gene_only <- round(counts_gene_only)
storage.mode(counts_gene_only) <- "integer"

meta_deseq <- meta_flat
meta_deseq <- meta_deseq[match(colnames(counts_gene_only), meta_deseq$counts_colname), , drop = FALSE]

if (any(is.na(meta_deseq$counts_colname))) {
  stop("Hay muestras en counts_gene_only que no pudieron alinearse con meta_deseq.")
}

cat("\nChequeo de alineación antes de DESeq2:\n")
cat("¿colnames(counts_gene_only) == meta_deseq$counts_colname?: ",
    all(colnames(counts_gene_only) == meta_deseq$counts_colname), "\n")

rownames(meta_deseq) <- meta_deseq$counts_colname

#---------------------------#
# 7) DESeqDataSet           #
#---------------------------#
cat("\n=== CREACIÓN DE DESeqDataSet ===\n")

coldata_s4 <- S4Vectors::DataFrame(meta_deseq)

dds <- DESeqDataSetFromMatrix(
  countData = counts_gene_only,
  colData = coldata_s4,
  design = ~ 1
)

cat("\nDimensiones iniciales del DDS:\n")
cat("- Genes:", nrow(dds), "\n")
cat("- Muestras:", ncol(dds), "\n")

#---------------------------#
# 8) Filtrado baja expresión#
#---------------------------#
cat("\n=== FILTRADO DE BAJA EXPRESIÓN ===\n")

keep_genes <- rowSums(counts(dds) >= 10) >= 10

cat("\nGenes antes del filtrado:", nrow(dds), "\n")
cat("Genes retenidos:", sum(keep_genes), "\n")
cat("Genes removidos:", sum(!keep_genes), "\n")

dds <- dds[keep_genes, ]

cat("\nDimensiones después del filtrado:\n")
cat("- Genes:", nrow(dds), "\n")
cat("- Muestras:", ncol(dds), "\n")

#---------------------------#
# 9) Normalización DESeq2   #
#---------------------------#
cat("\n=== NORMALIZACIÓN CON DESEQ2 ===\n")

dds <- estimateSizeFactors(dds)

cat("\nSize factors:\n")
print(sizeFactors(dds))

normalized_counts <- counts(dds, normalized = TRUE)

cat("\nResumen de size factors:\n")
print(summary(sizeFactors(dds)))

#---------------------------#
# 10) Transformación VST    #
#---------------------------#
cat("\n=== TRANSFORMACIÓN VST ===\n")

vsd <- vst(dds, blind = TRUE)
vst_mat <- assay(vsd)

cat("\nDimensiones de la matriz VST:\n")
cat("- Genes:", nrow(vst_mat), "\n")
cat("- Muestras:", ncol(vst_mat), "\n")

#---------------------------#
# 11) PCA corregido         #
#---------------------------#
cat("\n=== PCA CORREGIDO ===\n")

pca <- prcomp(t(vst_mat), center = TRUE, scale. = FALSE)

percent_var <- (pca$sdev^2) / sum(pca$sdev^2)
percent_var <- round(100 * percent_var, 2)

pca_df <- data.frame(
  sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  stringsAsFactors = FALSE
)

plot_cols <- c("counts_colname", "sample_type.x", "age_at_diagnosis_years", "year_of_diagnosis")
plot_cols <- plot_cols[plot_cols %in% colnames(sample_qc)]

meta_for_plot <- sample_qc[, plot_cols, drop = FALSE]

pca_df <- pca_df %>%
  left_join(meta_for_plot, by = c("sample" = "counts_colname"))

if ("sample_type.x" %in% colnames(pca_df)) {
  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = sample_type.x)) +
    geom_point(size = 3, alpha = 0.9) +
    labs(
      title = "PCA de muestras (DESeq2 VST, todos los genes)",
      x = paste0("PC1 (", percent_var[1], "%)"),
      y = paste0("PC2 (", percent_var[2], "%)"),
      color = "Sample type"
    ) +
    theme_bw(base_size = 12)
} else {
  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
    geom_point(size = 3, alpha = 0.9) +
    labs(
      title = "PCA de muestras (DESeq2 VST, todos los genes)",
      x = paste0("PC1 (", percent_var[1], "%)"),
      y = paste0("PC2 (", percent_var[2], "%)")
    ) +
    theme_bw(base_size = 12)
}

ggsave(
  filename = file.path(outdir, "pca_samples_deseq2.png"),
  plot = p_pca,
  width = 8,
  height = 6,
  dpi = 300
)

#---------------------------#
# 12) UMAP                  #
#---------------------------#
cat("\n=== UMAP ===\n")

set.seed(123)

umap_res <- uwot::umap(
  t(vst_mat),
  n_neighbors = 30,
  min_dist = 0.30,
  metric = "euclidean",
  scale = FALSE,
  verbose = TRUE
)

umap_df <- data.frame(
  sample = colnames(vst_mat),
  UMAP1 = umap_res[, 1],
  UMAP2 = umap_res[, 2],
  stringsAsFactors = FALSE
)

umap_df <- umap_df %>%
  left_join(meta_for_plot, by = c("sample" = "counts_colname"))

if ("sample_type.x" %in% colnames(umap_df)) {
  p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = sample_type.x)) +
    geom_point(size = 3, alpha = 0.9) +
    labs(
      title = "UMAP de muestras (DESeq2 VST, todos los genes)",
      x = "UMAP1",
      y = "UMAP2",
      color = "Sample type"
    ) +
    theme_bw(base_size = 12)
} else {
  p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
    geom_point(size = 3, alpha = 0.9) +
    labs(
      title = "UMAP de muestras (DESeq2 VST, todos los genes)",
      x = "UMAP1",
      y = "UMAP2"
    ) +
    theme_bw(base_size = 12)
}

ggsave(
  filename = file.path(outdir, "umap_samples_deseq2.png"),
  plot = p_umap,
  width = 8,
  height = 6,
  dpi = 300
)

#---------------------------#
# 13) Guardado de objetos   #
#---------------------------#
cat("\n=== GUARDADO DE OBJETOS ===\n")

save(
  counts_gene_only,
  technical_counts_summary,
  file = file.path(outdir, "counts_gene_only_raw.RData")
)

saveRDS(
  dds,
  file = file.path(outdir, "dds_filtered.rds")
)

saveRDS(
  vsd,
  file = file.path(outdir, "vst_object.rds")
)

write.csv(
  data.frame(gene = rownames(vst_mat), vst_mat, check.names = FALSE),
  file = file.path(outdir, "vst_matrix.csv"),
  row.names = FALSE
)

#---------------------------#
# 14) Reporte TXT           #
#---------------------------#
cat("\n=== GENERANDO REPORTE TXT ===\n")

report_file <- file.path(outdir, "qc_normalization_deseq2_report.txt")
sink(report_file)

cat("REPORTE QC + NORMALIZATION (DESEQ2)\n")
cat("Fecha:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("INPUT FILE:\n")
cat(input_file, "\n\n")

cat("DIMENSIONES INICIALES:\n")
cat("- counts_final filas:", nrow(counts_final), "\n")
cat("- counts_final columnas:", ncol(counts_final), "\n")
cat("- meta_final filas:", nrow(meta_final), "\n")
cat("- meta_final columnas:", ncol(meta_final), "\n\n")

cat("COLUMNAS LIST EXCLUIDAS DE META_FLAT:\n")
print(list_cols)
cat("\n")

cat("FILAS TÉCNICAS STAR ENCONTRADAS:\n")
print(present_star_rows)
cat("\n")

cat("DIMENSIONES DESPUÉS DE LIMPIEZA:\n")
cat("- counts_gene_only filas:", nrow(counts_gene_only), "\n")
cat("- counts_gene_only columnas:", ncol(counts_gene_only), "\n\n")

cat("RESUMEN QC DE MUESTRAS:\n")
cat("\nlibrary_size_raw:\n")
print(summary(sample_qc$library_size_raw))

cat("\npct_technical:\n")
print(summary(sample_qc$pct_technical))

cat("\ngenes_detected_gt0:\n")
print(summary(sample_qc$genes_detected_gt0))

cat("\nFILTRADO DE BAJA EXPRESIÓN:\n")
cat("- genes antes:", length(keep_genes), "\n")
cat("- genes retenidos:", sum(keep_genes), "\n")
cat("- genes removidos:", sum(!keep_genes), "\n\n")

cat("SIZE FACTORS:\n")
print(sizeFactors(dds))
cat("\n")

cat("RESUMEN SIZE FACTORS:\n")
print(summary(sizeFactors(dds)))
cat("\n")

cat("DIMENSIONES VST:\n")
cat("- filas:", nrow(vst_mat), "\n")
cat("- columnas:", ncol(vst_mat), "\n\n")

cat("VARIANZA EXPLICADA PCA:\n")
print(percent_var[1:10])
cat("\n")

cat("PARÁMETROS PCA:\n")
cat("- input: matriz VST\n")
cat("- genes usados: todos los genes filtrados\n")
cat("- center = TRUE\n")
cat("- scale = FALSE\n\n")

cat("PARÁMETROS UMAP:\n")
cat("- input: matriz VST\n")
cat("- genes usados: todos los genes filtrados\n")
cat("- n_neighbors = 15\n")
cat("- min_dist = 0.30\n")
cat("- metric = euclidean\n")
cat("- scale = FALSE\n")
cat("- seed = 123\n\n")

cat("ARCHIVOS GENERADOS:\n")
cat("- counts_gene_only_raw.RData\n")
cat("- dds_filtered.rds\n")
cat("- vst_object.rds\n")
cat("- vst_matrix.csv\n")
cat("- sample_qc_metadata.csv\n")
cat("- pca_samples_deseq2.png\n")
cat("- umap_samples_deseq2.png\n")
cat("- qc_normalization_deseq2_report.txt\n")

sink()

#---------------------------#
# 15) Mensaje final         #
#---------------------------#
cat("\n=== SCRIPT 06 COMPLETADO ===\n")
cat("\nArchivos generados en:\n")
cat(outdir, "\n")



#---------------------------#
# 16) CLUSTERING + VALIDACIÓN
#---------------------------#
cat("\n=== CLUSTERING (PCA + KMEANS) ===\n")

library(cluster)

# 1. PCA (ya lo tienes, pero lo rehacemos por claridad)
pca <- prcomp(t(vst_mat), center = TRUE, scale. = FALSE)

# usar primeras 10 PCs
pc_mat <- pca$x[, 1:10]

cat("\nDimensiones PCA usadas para clustering:\n")
print(dim(pc_mat))

#---------------------------#
# 2. K-means (k = 3)
#---------------------------#
set.seed(123)

k <- 3
kmeans_res <- kmeans(pc_mat, centers = k, nstart = 50)

clusters <- kmeans_res$cluster

cat("\nTamaño de clusters:\n")
print(table(clusters))

#---------------------------#
# 3. Silhouette score
#---------------------------#
cat("\n=== SILHOUETTE ===\n")

dist_mat <- dist(pc_mat)

sil <- silhouette(clusters, dist_mat)

sil_mean <- mean(sil[, 3])

cat("\nSilhouette promedio:\n")
print(sil_mean)

# guardar silhouette por muestra
sil_df <- data.frame(
  sample = rownames(pc_mat),
  cluster = clusters,
  silhouette = sil[, 3]
)

write.csv(
  sil_df,
  file = file.path(outdir, "silhouette_scores.csv"),
  row.names = FALSE
)

#---------------------------#
# 4. Añadir clusters a metadata
#---------------------------#
cat("\n=== AÑADIENDO CLUSTERS A METADATA ===\n")

meta_final$cluster_k3 <- clusters[meta_final$counts_colname]

# versión factor (mejor para plots)
meta_final$cluster_k3 <- as.factor(meta_final$cluster_k3)

# guardar metadata actualizada
write.csv(
  meta_final[, !sapply(meta_final, is.list)],
  file = file.path(outdir, "metadata_with_clusters.csv"),
  row.names = FALSE
)

#---------------------------#
# PCA con clusters
#---------------------------#
cat("\n=== PCA CON CLUSTERS ===\n")

pca_df$cluster_k3 <- meta_final$cluster_k3[
  match(pca_df$sample, meta_final$counts_colname)
]

p_pca_cluster <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster_k3)) +
  geom_point(size = 3, alpha = 0.9) +
  labs(
    title = "PCA con clusters k-means (k=3)",
    x = paste0("PC1 (", percent_var[1], "%)"),
    y = paste0("PC2 (", percent_var[2], "%)"),
    color = "Cluster"
  ) +
  theme_bw(base_size = 12)

ggsave(
  filename = file.path(outdir, "pca_clusters_k3.png"),
  plot = p_pca_cluster,
  width = 8,
  height = 6,
  dpi = 300
)

#ALTERNATIVO
#---------------------------#
# PCA con clusters
#---------------------------#
cat("\n=== PCA CON CLUSTERS ===\n")

pca_df$cluster_k3 <- meta_final$cluster_k3[
  match(pca_df$sample, meta_final$counts_colname)
]

p_pca_cluster <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster_k3)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = c(
    "1" = "#00BA38",   # verde
    "2" = "#F8766D",   # salmon/rojizo
    "3" = "#619CFF"    # azul
  )) +
  labs(
    title = "PCA con clusters k-means (k=3)",
    x = paste0("PC1 (", percent_var[1], "%)"),
    y = paste0("PC2 (", percent_var[2], "%)"),
    color = "Cluster"
  ) +
  theme_bw(base_size = 12)

ggsave(
  filename = file.path(outdir, "pca_clusters_k3_ALT.png"),
  plot = p_pca_cluster,
  width = 8,
  height = 6,
  dpi = 300
)

#---------------------------#
# 5. Añadir clusters a UMAP
#---------------------------#
cat("\n=== UMAP CON CLUSTERS ===\n")

umap_df$cluster_k3 <- meta_final$cluster_k3[match(umap_df$sample, meta_final$counts_colname)]

p_umap_cluster <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = cluster_k3)) +
  geom_point(size = 3, alpha = 0.9) +
  labs(
    title = "UMAP con clusters k-means (k=3)",
    x = "UMAP1",
    y = "UMAP2",
    color = "Cluster"
  ) +
  theme_bw(base_size = 12)

ggsave(
  filename = file.path(outdir, "umap_clusters_k3.png"),
  plot = p_umap_cluster,
  width = 8,
  height = 6,
  dpi = 300
)


#ALTERNATIVO
#---------------------------#
# 5. Añadir clusters a UMAP
#---------------------------#
cat("\n=== UMAP CON CLUSTERS ===\n")

umap_df$cluster_k3 <- meta_final$cluster_k3[
  match(umap_df$sample, meta_final$counts_colname)
]

p_umap_cluster <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = cluster_k3)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = c(
    "1" = "#00BA38",   # verde
    "2" = "#F8766D",   # salmon/rojizo
    "3" = "#619CFF"    # azul
  )) +
  labs(
    title = "UMAP con clusters k-means (k=3)",
    x = "UMAP1",
    y = "UMAP2",
    color = "Cluster"
  ) +
  theme_bw(base_size = 12)

ggsave(
  filename = file.path(outdir, "umap_clusters_k3_ALT.png"),
  plot = p_umap_cluster,
  width = 8,
  height = 6,
  dpi = 300
)




#---------------------------#
# 6. Comparación con sample_type
#---------------------------#
cat("\n=== CRUCE CLUSTER vs SAMPLE TYPE ===\n")

if ("sample_type.x" %in% colnames(meta_final)) {
  print(table(meta_final$cluster_k3, meta_final$sample_type.x))
}

#---------------------------#
# 7. Guardado final
#---------------------------#
saveRDS(
  meta_final,
  file = file.path(outdir, "meta_with_clusters.rds")
)

cat("\n=== CLUSTERING COMPLETADO ===\n")
cat("\nSilhouette promedio:", sil_mean, "\n")