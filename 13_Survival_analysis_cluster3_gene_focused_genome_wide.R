############################################################
### 13_SURVIVAL_ANALYSIS_CLUSTER3_GENE_FOCUSED.R
### GOAL:
###   Survival analysis restricted to Cluster 3 only
###   Gene-focused Cox and KM analysis
###
### INPUT:
###   - TARGET_ALL_P2_PRIMARY_U18_FINAL_FOR_ANALYSIS.RData
###   - meta_with_clusters.rds
###
### OUTPUT:
###   - Cox gene-wise survival tables
###   - KM plots per gene
###   - Forest plot
###   - PH assumption tables
###   - C-index comparison
###   - TXT report
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(DESeq2)
  library(survival)
  library(survminer)
  library(broom)
  library(ggplot2)
  library(patchwork)
  library(biomaRt)
})

cat("=== SCRIPT 13: CLUSTER 3 GENE-FOCUSED SURVIVAL ANALYSIS ===\n")

#---------------------------#
# 1) Editable gene list     #
#---------------------------#

target_genes <- c(
  "BTLA",
  "ENTPD1",   # CD39
  "NT5E",     # CD73
  
  "PDCD1",
  "CD274",
  "CTLA4",
  "LAG3",
  "HAVCR2",
  "TIGIT",
  "VSIR",
  
  "CD47",
  "SIRPA",
  "LGALS9",
  
  "IDO1",
  "ARG1",
  "IL10",
  "TGFB1",
  "FOXP3",
  
  "CXCR4",
  "CCR7",
  "SELL",
  "MKI67"
)

target_genes <- unique(target_genes)

#---------------------------#
# 2) Paths                  #
#---------------------------#

input_counts_file <- "subsets/primary_u18/final_analysis/TARGET_ALL_P2_PRIMARY_U18_FINAL_FOR_ANALYSIS.RData"
input_meta_clusters <- "subsets/primary_u18/final_analysis/06_qc_normalization_deseq2/meta_with_clusters.rds"

outdir <- "subsets/primary_u18/final_analysis/13_survival_cluster3_gene_focused_genome_wide"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

#---------------------------#
# 3) Load data              #
#---------------------------#

cat("\n=== LOADING DATA ===\n")

load(input_counts_file)
meta_clusters <- readRDS(input_meta_clusters)

cat("- Raw rows in counts_final:", nrow(counts_final), "\n")
cat("- Raw samples in counts_final:", ncol(counts_final), "\n")

#---------------------------#
# 4) Remove STAR technical rows
#    and clean ENSEMBL IDs
#---------------------------#

cat("\n=== CLEANING COUNT MATRIX ===\n")

technical_rows <- c(
  "N_unmapped",
  "N_multimapping",
  "N_noFeature",
  "N_ambiguous"
)

counts_final <- counts_final[
  !rownames(counts_final) %in% technical_rows,
  ,
  drop = FALSE
]

# Remove ENSEMBL version suffix: ENSG000001234.5 -> ENSG000001234
rownames(counts_final) <- sub("\\..*$", "", rownames(counts_final))

# Collapse duplicated ENSEMBL IDs if any appeared after removing version suffix
if (anyDuplicated(rownames(counts_final)) > 0) {
  cat("Duplicated ENSEMBL IDs detected after version removal. Collapsing by sum.\n")
  
  counts_final <- rowsum(
    as.matrix(counts_final),
    group = rownames(counts_final)
  )
}

cat("- Rows after removing technical rows:", nrow(counts_final), "\n")
cat("- Example cleaned rownames:\n")
print(head(rownames(counts_final)))

#---------------------------#
# 5) Align metadata         #
#---------------------------#

cat("\n=== ALIGNING METADATA ===\n")

meta_clusters <- meta_clusters[
  match(colnames(counts_final), meta_clusters$counts_colname),
  ,
  drop = FALSE
]

stopifnot(all(colnames(counts_final) == meta_clusters$counts_colname))

cat("- Samples after metadata alignment:", nrow(meta_clusters), "\n")

#---------------------------#
# 6) Build survival metadata
#    restricted to Cluster 3
#---------------------------#

cat("\n=== BUILDING CLUSTER 3 SURVIVAL DATASET ===\n")

meta_surv <- meta_clusters %>%
  mutate(
    OS_days = case_when(
      vital_status == "Dead"  ~ suppressWarnings(as.numeric(days_to_death)),
      vital_status == "Alive" ~ suppressWarnings(as.numeric(days_to_last_follow_up)),
      TRUE ~ NA_real_
    ),
    OS_event = case_when(
      vital_status == "Dead" ~ 1,
      vital_status == "Alive" ~ 0,
      TRUE ~ NA_real_
    ),
    cluster_k3 = factor(cluster_k3),
    age = suppressWarnings(as.numeric(age_at_diagnosis_years)),
    year = suppressWarnings(as.numeric(year_of_diagnosis))
  ) %>%
  filter(
    cluster_k3 == "3",
    !is.na(OS_days),
    OS_days > 0,
    !is.na(OS_event),
    !is.na(age),
    !is.na(year)
  )

cat("\nSamples in Cluster 3 survival analysis:", nrow(meta_surv), "\n")

cat("\nEvents:\n")
print(table(meta_surv$OS_event, useNA = "ifany"))

cat("\nOS_days summary:\n")
print(summary(meta_surv$OS_days))

cat("\nAge summary:\n")
print(summary(meta_surv$age))

cat("\nYear summary:\n")
print(summary(meta_surv$year))

write.csv(
  meta_surv[, !sapply(meta_surv, is.list), drop = FALSE],
  file = file.path(outdir, "cluster3_survival_metadata_used.csv"),
  row.names = FALSE
)

saveRDS(
  meta_surv,
  file = file.path(outdir, "cluster3_survival_metadata_used.rds")
)

#---------------------------#
# 7) Subset counts to Cluster 3
#---------------------------#

cat("\n=== SUBSETTING COUNTS TO CLUSTER 3 ===\n")

counts_cluster3 <- counts_final[
  ,
  meta_surv$counts_colname,
  drop = FALSE
]

stopifnot(all(colnames(counts_cluster3) == meta_surv$counts_colname))

cat("- Genes:", nrow(counts_cluster3), "\n")
cat("- Cluster 3 samples:", ncol(counts_cluster3), "\n")

#---------------------------#
# 8) VST normalization      #
#---------------------------#

cat("\n=== VST NORMALIZATION ===\n")

dds <- DESeqDataSetFromMatrix(
  countData = round(counts_cluster3),
  colData = meta_surv,
  design = ~ 1
)

dds <- estimateSizeFactors(dds)

vst_mat <- assay(
  vst(dds, blind = TRUE)
)

cat("- VST matrix generated.\n")
cat("- VST genes:", nrow(vst_mat), "\n")
cat("- VST samples:", ncol(vst_mat), "\n")

#---------------------------#
# 9) ENSEMBL to HGNC annotation
#---------------------------#

cat("\n=== ANNOTATING ENSEMBL IDS TO HGNC SYMBOLS ===\n")

ensembl_ids <- rownames(vst_mat)

annotation_df <- data.frame(
  ensembl_gene_id = ensembl_ids,
  original_rowname = ensembl_ids,
  stringsAsFactors = FALSE
)

cat("\nRecuperando anotación con getBM...\n")

mart <- NULL

mart_attempts <- list(
  list(host = "https://www.ensembl.org"),
  list(host = "https://useast.ensembl.org"),
  list(host = "https://uswest.ensembl.org"),
  list(host = "https://asia.ensembl.org")
)

for (attempt in mart_attempts) {
  cat("Trying Ensembl host:", attempt$host, "\n")
  
  mart <- tryCatch(
    {
      useMart(
        biomart = "ENSEMBL_MART_ENSEMBL",
        dataset = "hsapiens_gene_ensembl",
        host = attempt$host
      )
    },
    error = function(e) {
      cat("Failed:", conditionMessage(e), "\n")
      NULL
    }
  )
  
  if (!is.null(mart)) {
    cat("Connected successfully to:", attempt$host, "\n")
    break
  }
}

if (is.null(mart)) {
  stop("Could not connect to any Ensembl mirror. Check internet connection or try again later.")
}

anno_bm <- tryCatch(
  {
    getBM(
      attributes = c(
        "ensembl_gene_id",
        "hgnc_symbol",
        "description"
      ),
      filters = "ensembl_gene_id",
      values = unique(annotation_df$ensembl_gene_id),
      mart = mart
    )
  },
  error = function(e) {
    stop("getBM failed even after connecting to Ensembl: ", conditionMessage(e))
  }
)

anno_bm <- anno_bm %>%
  dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)

annotation_df <- annotation_df %>%
  dplyr::left_join(anno_bm, by = "ensembl_gene_id") %>%
  dplyr::mutate(
    hgnc_symbol = ifelse(is.na(hgnc_symbol), "", hgnc_symbol),
    description = ifelse(is.na(description), "", description),
    gene_symbol = ifelse(hgnc_symbol == "", NA, hgnc_symbol),
    gene_label  = ifelse(is.na(gene_symbol), ensembl_gene_id, gene_symbol)
  )

write.csv(
  annotation_df,
  file = file.path(outdir, "gene_annotation_table.csv"),
  row.names = FALSE
)

gene_names_final <- annotation_df$gene_label
names(gene_names_final) <- annotation_df$ensembl_gene_id

rownames(vst_mat) <- gene_names_final[rownames(vst_mat)]

# Safety check
rownames(vst_mat) <- ifelse(
  is.na(rownames(vst_mat)) | rownames(vst_mat) == "",
  ensembl_ids,
  rownames(vst_mat)
)

cat("- Annotated genes with HGNC symbols:", sum(!is.na(annotation_df$gene_symbol)), "\n")

# Collapse duplicated HGNC symbols by mean VST expression
if (anyDuplicated(rownames(vst_mat)) > 0) {
  cat("Duplicated HGNC symbols detected. Collapsing by mean VST expression.\n")
  
  vst_mat <- rowsum(
    vst_mat,
    group = rownames(vst_mat)
  ) / as.vector(table(rownames(vst_mat)))
}

cat("- Final VST matrix rows after annotation/collapse:", nrow(vst_mat), "\n")

#---------------------------#
# 10) Genome-wide gene list #
#---------------------------#

cat("\n=== PREPARING GENOME-WIDE GENE SURVIVAL SCREEN ===\n")

all_genes <- rownames(vst_mat)

# Remove ENSEMBL-like unnamed genes if desired
# Keep this as FALSE if you want absolutely everything
keep_only_hgnc_symbols <- TRUE

if (keep_only_hgnc_symbols) {
  all_genes <- all_genes[!grepl("^ENSG", all_genes)]
}

# Remove genes with zero or near-zero variance
gene_sd <- apply(vst_mat[all_genes, , drop = FALSE], 1, sd, na.rm = TRUE)

all_genes <- all_genes[
  !is.na(gene_sd) &
    is.finite(gene_sd) &
    gene_sd > 0.1
]

cat("- Genes tested after filtering:", length(all_genes), "\n")

if (length(all_genes) == 0) {
  stop("No genes available for genome-wide survival analysis.")
}

#---------------------------#
# 11) Cox analysis all genes#
#---------------------------#

cat("\n=== RUNNING GENOME-WIDE COX SURVIVAL SCREEN ===\n")

cox_results <- list()
zph_results <- list()
cindex_results <- list()

cox_base <- coxph(
  Surv(OS_days, OS_event) ~ age + year,
  data = meta_surv
)

base_cindex <- summary(cox_base)$concordance[1]

for (gene in all_genes) {
  
  gene_expr <- as.numeric(vst_mat[gene, meta_surv$counts_colname])
  
  df_gene <- meta_surv %>%
    dplyr::mutate(
      gene = gene,
      expr = gene_expr,
      expr_z = as.numeric(scale(expr))
    ) %>%
    dplyr::filter(
      !is.na(expr_z),
      is.finite(expr_z)
    )
  
  if (nrow(df_gene) < 20) next
  if (sd(df_gene$expr_z, na.rm = TRUE) == 0) next
  
  cox_adj <- tryCatch(
    {
      coxph(
        Surv(OS_days, OS_event) ~ expr_z + age + year,
        data = df_gene
      )
    },
    error = function(e) NULL
  )
  
  if (is.null(cox_adj)) next
  
  tidy_adj <- tryCatch(
    {
      broom::tidy(
        cox_adj,
        exponentiate = TRUE,
        conf.int = TRUE
      ) %>%
        dplyr::filter(term == "expr_z") %>%
        dplyr::mutate(
          gene = gene,
          model = "Adjusted age + year",
          interpretation = ifelse(
            estimate > 1,
            "Higher expression = higher risk",
            "Higher expression = lower risk"
          )
        )
    },
    error = function(e) NULL
  )
  
  if (is.null(tidy_adj)) next
  if (nrow(tidy_adj) == 0) next
  
  cox_results[[gene]] <- tidy_adj
  
  gene_cindex <- tryCatch(
    {
      summary(cox_adj)$concordance[1]
    },
    error = function(e) NA_real_
  )
  
  cindex_results[[gene]] <- data.frame(
    gene = gene,
    base_model_cindex = base_cindex,
    gene_model_cindex = gene_cindex,
    delta_cindex = gene_cindex - base_cindex
  )
}

if (length(cox_results) == 0) {
  stop("No Cox models were successfully fitted.")
}

#---------------------------#
# 12) Combine genome-wide results
#---------------------------#

cat("\n=== COMBINING GENOME-WIDE RESULTS ===\n")

top_n_genes <- 50

cox_table <- dplyr::bind_rows(cox_results) %>%
  dplyr::mutate(
    FDR = p.adjust(p.value, method = "BH"),
    HR_CI = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high),
    p_value_clean = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)),
    FDR_clean = ifelse(FDR < 0.001, "<0.001", sprintf("%.3f", FDR)),
    significance = dplyr::case_when(
      FDR < 0.05 & estimate > 1 ~ "Significant higher risk",
      FDR < 0.05 & estimate < 1 ~ "Significant protective",
      p.value < 0.05 ~ "Nominal only",
      TRUE ~ "Not significant"
    )
  ) %>%
  dplyr::arrange(FDR, p.value)

cindex_table <- dplyr::bind_rows(cindex_results) %>%
  dplyr::arrange(dplyr::desc(delta_cindex))

significant_genes_all <- cox_table %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::arrange(FDR, p.value)

significant_genes <- significant_genes_all %>%
  dplyr::slice_head(n = top_n_genes)

nominal_genes <- cox_table %>%
  dplyr::filter(p.value < 0.05, FDR >= 0.05) %>%
  dplyr::arrange(p.value)

cat("\nFDR-significant genes total:", nrow(significant_genes_all), "\n")
cat("Top genes selected for plots:", nrow(significant_genes), "\n")
cat("Nominal genes p < 0.05 but FDR >= 0.05:", nrow(nominal_genes), "\n")

write.csv(
  cox_table,
  file = file.path(outdir, "cluster3_genomewide_survival_cox_all_genes.csv"),
  row.names = FALSE
)

write.csv(
  significant_genes_all,
  file = file.path(outdir, "cluster3_genomewide_survival_ALL_FDR_significant_genes.csv"),
  row.names = FALSE
)

write.csv(
  significant_genes,
  file = file.path(outdir, "cluster3_genomewide_survival_TOP50_FDR_significant_genes.csv"),
  row.names = FALSE
)

write.csv(
  nominal_genes,
  file = file.path(outdir, "cluster3_genomewide_survival_nominal_genes.csv"),
  row.names = FALSE
)

write.csv(
  cindex_table,
  file = file.path(outdir, "cluster3_genomewide_survival_cindex_comparison.csv"),
  row.names = FALSE
)
#---------------------------#
# 13) Forest plot significant genes
#---------------------------#

cat("\n=== GENERATING FOREST PLOT FOR SIGNIFICANT GENES ===\n")

if (nrow(significant_genes) > 0) {
  
  forest_df <- significant_genes %>%
    dplyr::mutate(
      gene = factor(gene, levels = rev(gene)),
      status = dplyr::case_when(
        estimate > 1 ~ "Higher risk",
        estimate < 1 ~ "Protective",
        TRUE ~ "Neutral"
      )
    )
  
  p_forest <- ggplot(forest_df, aes(x = estimate, y = gene)) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      color = "gray40"
    ) +
    geom_errorbar(
      aes(xmin = conf.low, xmax = conf.high),
      orientation = "y",
      width = 0.2,
      linewidth = 0.7
    ) +
    geom_point(size = 3) +
    scale_x_log10() +
    labs(
      title = "Cluster 3 genome-wide survival screen",
      subtitle = "FDR-significant genes; adjusted Cox: expression + age + year",
      x = "Hazard ratio per 1 SD increase in VST expression",
      y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.y = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  ggsave(
    filename = file.path(outdir, "cluster3_genomewide_survival_forest_FDR_significant.png"),
    plot = p_forest,
    width = 9,
    height = max(5, nrow(significant_genes) * 0.35),
    dpi = 300
  )
  
  ggsave(
    filename = file.path(outdir, "cluster3_genomewide_survival_forest_FDR_significant.pdf"),
    plot = p_forest,
    width = 9,
    height = max(5, nrow(significant_genes) * 0.35)
  )
}

#---------------------------#
# 14) KM plots only significant genes
#---------------------------#

cat("\n=== GENERATING KM PLOTS FOR SIGNIFICANT GENES ONLY ===\n")

if (nrow(significant_genes) > 0) {
  
  km_dir <- file.path(outdir, "KM_FDR_significant_genes")
  dir.create(km_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (gene in significant_genes$gene) {
    
    gene_expr <- as.numeric(vst_mat[gene, meta_surv$counts_colname])
    
    df_gene <- meta_surv %>%
      dplyr::mutate(
        expr = gene_expr,
        expr_group = ifelse(
          expr >= median(expr, na.rm = TRUE),
          "High",
          "Low"
        ),
        expr_group = factor(expr_group, levels = c("Low", "High"))
      )
    
    fit_km <- survfit(
      Surv(OS_days, OS_event) ~ expr_group,
      data = df_gene
    )
    
    km_plot <- ggsurvplot(
      fit_km,
      data = df_gene,
      pval = TRUE,
      risk.table = TRUE,
      conf.int = TRUE,
      title = paste0("Cluster 3 OS by ", gene, " expression"),
      xlab = "Days",
      ylab = "Overall survival probability",
      legend.title = gene,
      legend.labs = c("Low expression", "High expression"),
      risk.table.height = 0.25,
      ggtheme = theme_bw(base_size = 12)
    )
    
    png(
      filename = file.path(km_dir, paste0("KM_cluster3_", gene, "_high_vs_low.png")),
      width = 2400,
      height = 2200,
      res = 300
    )
    print(km_plot)
    dev.off()
    
    pdf(
      file = file.path(km_dir, paste0("KM_cluster3_", gene, "_high_vs_low.pdf")),
      width = 8,
      height = 7
    )
    print(km_plot)
    dev.off()
  }
}

#---------------------------#
# 15) TXT report
#---------------------------#

cat("\n=== GENERATING TXT REPORT ===\n")

report_file <- file.path(outdir, "cluster3_genomewide_gene_survival_report.txt")

sink(report_file)

cat("CLUSTER 3 GENOME-WIDE GENE-FOCUSED SURVIVAL ANALYSIS REPORT\n\n")

cat("Dataset: TARGET-ALL-P2 PRIMARY U18\n")
cat("Analysis restricted to: cluster_k3 == 3\n\n")

cat("Survival definition:\n")
cat("- OS_days = days_to_death if Dead\n")
cat("- OS_days = days_to_last_follow_up if Alive\n")
cat("- OS_event: Dead = 1, Alive = 0\n\n")

cat("Samples used:", nrow(meta_surv), "\n")
cat("Events:", sum(meta_surv$OS_event == 1), "\n")
cat("Censored:", sum(meta_surv$OS_event == 0), "\n\n")

cat("Genes tested after filtering:", length(all_genes), "\n")
cat("Cox models successfully fitted:", nrow(cox_table), "\n")
cat("FDR-significant genes:", nrow(significant_genes), "\n")
cat("Nominal genes:", nrow(nominal_genes), "\n\n")

cat("Top FDR-significant genes:\n")
print(significant_genes)
cat("\n\n")

cat("Top 50 genes by FDR:\n")
print(head(cox_table, 50))
cat("\n\n")

cat("Interpretation guide:\n")
cat("- HR > 1: higher gene expression is associated with increased risk of death.\n")
cat("- HR < 1: higher gene expression is associated with reduced risk of death.\n")
cat("- HR is interpreted per 1 SD increase in VST-normalized expression.\n")
cat("- FDR < 0.05: significant after genome-wide multiple-testing correction.\n")
cat("- p < 0.05 but FDR >= 0.05: nominal association only.\n")
cat("- KM plots are generated only for FDR-significant genes.\n")

sink()

cat("\n=== GENOME-WIDE CLUSTER 3 SURVIVAL ANALYSIS COMPLETED ===\n")
