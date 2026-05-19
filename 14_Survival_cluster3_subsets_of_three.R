############################################################
### 14_SURVIVAL_CLUSTER3_THREE_GENE_MODEL.R
### GOAL:
###   Survival analysis restricted to Cluster 3 only
###   Gene-only multigene Cox model using up to 3 selected genes
###
### MODEL:
###   Surv(OS_days, OS_event) ~ gene1_z + gene2_z + gene3_z
###
### OUTPUT:
###   - Multigene Cox table
###   - Forest plot
###   - PH assumption table
###   - C-index
###   - Risk score
###   - KM plot by risk score median split
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

cat("=== SCRIPT 14: CLUSTER 3 THREE-GENE SURVIVAL MODEL ===\n")

#---------------------------#
# 1) Editable gene list
#---------------------------#

# Use maximum 3 genes
target_genes <- c(
  "GPR176", "ODC1", 
  "NGRN" 
  #"TUNAR" 
)

target_genes <- unique(target_genes)

if (length(target_genes) == 0) {
  stop("You must define 1 to 3 target genes before running the script.")
}

if (length(target_genes) > 3) {
  stop("Use a maximum of 3 genes for this multigene Cox model.")
}

#---------------------------#
# 2) Settings
#---------------------------#

min_count <- 10
min_prop_samples <- 0.20

#---------------------------#
# 3) Paths
#---------------------------#

input_counts_file <- "subsets/primary_u18/final_analysis/TARGET_ALL_P2_PRIMARY_U18_FINAL_FOR_ANALYSIS.RData"
input_meta_clusters <- "subsets/primary_u18/final_analysis/06_qc_normalization_deseq2/meta_with_clusters.rds"

outdir <- "subsets/primary_u18/final_analysis/14_survival_cluster3_three_gene_model"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

#---------------------------#
# 4) Load data
#---------------------------#

cat("\n=== LOADING DATA ===\n")

load(input_counts_file)
meta_clusters <- readRDS(input_meta_clusters)

cat("- Raw rows in counts_final:", nrow(counts_final), "\n")
cat("- Raw samples in counts_final:", ncol(counts_final), "\n")

#---------------------------#
# 5) Clean count matrix
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

rownames(counts_final) <- sub("\\..*$", "", rownames(counts_final))

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
# 6) Align metadata
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
# 7) Build Cluster 3 survival metadata
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
    cluster_k3 = factor(cluster_k3)
  ) %>%
  filter(
    cluster_k3 == "3",
    !is.na(OS_days),
    OS_days > 0,
    !is.na(OS_event)
  )

cat("\nSamples in Cluster 3 survival analysis:", nrow(meta_surv), "\n")

cat("\nEvents:\n")
print(table(meta_surv$OS_event, useNA = "ifany"))

cat("\nOS_days summary:\n")
print(summary(meta_surv$OS_days))

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
# 8) Subset counts to Cluster 3
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
# 9) Expression filtering
#---------------------------#

cat("\n=== FILTERING GENES BY MINIMUM EXPRESSION ===\n")

min_samples <- ceiling(min_prop_samples * ncol(counts_cluster3))

keep_expr <- rowSums(counts_cluster3 >= min_count, na.rm = TRUE) >= min_samples

counts_cluster3_filtered <- counts_cluster3[
  keep_expr,
  ,
  drop = FALSE
]

cat("- Genes before filtering:", nrow(counts_cluster3), "\n")
cat("- Genes after filtering:", nrow(counts_cluster3_filtered), "\n")
cat("- Minimum samples required with count >=", min_count, ":", min_samples, "\n")

#---------------------------#
# 10) VST normalization
#---------------------------#

cat("\n=== VST NORMALIZATION ===\n")

dds_filtered <- DESeqDataSetFromMatrix(
  countData = round(counts_cluster3_filtered),
  colData = meta_surv,
  design = ~ 1
)

dds_filtered <- estimateSizeFactors(dds_filtered)

vst_mat <- assay(
  vst(dds_filtered, blind = TRUE)
)

cat("- Filtered VST genes:", nrow(vst_mat), "\n")
cat("- Filtered VST samples:", ncol(vst_mat), "\n")

#---------------------------#
# 11) ENSEMBL to HGNC annotation
#---------------------------#

cat("\n=== ANNOTATING ENSEMBL IDS TO HGNC SYMBOLS ===\n")

ensembl_ids <- rownames(vst_mat)

annotation_df <- data.frame(
  ensembl_gene_id = ensembl_ids,
  original_rowname = ensembl_ids,
  stringsAsFactors = FALSE
)

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
  stop("Could not connect to any Ensembl mirror.")
}

anno_bm <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "hgnc_symbol",
    "description"
  ),
  filters = "ensembl_gene_id",
  values = unique(annotation_df$ensembl_gene_id),
  mart = mart
)

anno_bm <- anno_bm %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

annotation_df <- annotation_df %>%
  left_join(anno_bm, by = "ensembl_gene_id") %>%
  mutate(
    hgnc_symbol = ifelse(is.na(hgnc_symbol), "", hgnc_symbol),
    description = ifelse(is.na(description), "", description),
    gene_symbol = ifelse(hgnc_symbol == "", NA, hgnc_symbol),
    gene_label = ifelse(is.na(gene_symbol), ensembl_gene_id, gene_symbol)
  )

write.csv(
  annotation_df,
  file = file.path(outdir, "gene_annotation_table.csv"),
  row.names = FALSE
)

gene_names_final <- annotation_df$gene_label
names(gene_names_final) <- annotation_df$ensembl_gene_id

new_gene_names <- gene_names_final[rownames(vst_mat)]

new_gene_names <- ifelse(
  is.na(new_gene_names) | new_gene_names == "",
  rownames(vst_mat),
  new_gene_names
)

rownames(vst_mat) <- new_gene_names

if (anyDuplicated(rownames(vst_mat)) > 0) {
  cat("Duplicated HGNC symbols detected. Collapsing by mean VST expression.\n")
  
  vst_mat <- rowsum(
    vst_mat,
    group = rownames(vst_mat)
  ) / as.vector(table(rownames(vst_mat)))
}

cat("- Final VST genes after annotation/collapse:", nrow(vst_mat), "\n")

#---------------------------#
# 12) Check target genes
#---------------------------#

cat("\n=== CHECKING TARGET GENES ===\n")

available_genes <- rownames(vst_mat)

missing_genes <- setdiff(target_genes, available_genes)
genes_to_model <- intersect(target_genes, available_genes)

if (length(missing_genes) > 0) {
  cat("\nWARNING: These target genes were not found after filtering/annotation:\n")
  print(missing_genes)
}

if (length(genes_to_model) < 1) {
  stop("None of the selected target genes are available after expression filtering and annotation.")
}

if (length(genes_to_model) > 3) {
  stop("More than 3 genes available. This script is intended for a maximum of 3 genes.")
}

cat("\nGenes used in final model:\n")
print(genes_to_model)

#---------------------------#
# 13) Build modeling dataframe
#---------------------------#

cat("\n=== BUILDING MODELING DATAFRAME ===\n")

expr_df <- t(vst_mat[genes_to_model, meta_surv$counts_colname, drop = FALSE]) %>%
  as.data.frame() %>%
  rownames_to_column("counts_colname")

colnames(expr_df) <- make.names(colnames(expr_df))

genes_model_safe <- make.names(genes_to_model)

model_df <- meta_surv %>%
  left_join(expr_df, by = "counts_colname")

for (g in genes_model_safe) {
  model_df[[paste0(g, "_z")]] <- as.numeric(scale(model_df[[g]]))
}

gene_z_terms <- paste0(genes_model_safe, "_z")

model_df <- model_df %>%
  filter(
    if_all(all_of(gene_z_terms), ~ !is.na(.) & is.finite(.))
  )

cat("- Samples in final model dataframe:", nrow(model_df), "\n")

#---------------------------#
# 14) Statistical power summary
#---------------------------#

cat("\n=== STATISTICAL POWER SUMMARY ===\n")

n_samples <- nrow(model_df)
n_events <- sum(model_df$OS_event == 1)
n_predictors <- length(gene_z_terms)
events_per_variable <- n_events / n_predictors

cat("- Samples:", n_samples, "\n")
cat("- Events:", n_events, "\n")
cat("- Censored:", sum(model_df$OS_event == 0), "\n")
cat("- Predictors:", n_predictors, "\n")
cat("- Events per variable:", round(events_per_variable, 2), "\n")

if (events_per_variable < 10) {
  cat("\nWARNING: EPV < 10. Interpret this multigene model cautiously.\n")
}

#---------------------------#
# 15) Multigene Cox model
#---------------------------#

cat("\n=== RUNNING THREE-GENE COX MODEL ===\n")

formula_multigene <- as.formula(
  paste(
    "Surv(OS_days, OS_event) ~",
    paste(gene_z_terms, collapse = " + ")
  )
)

cox_multigene <- coxph(
  formula_multigene,
  data = model_df
)

cox_summary <- summary(cox_multigene)

cox_table <- broom::tidy(
  cox_multigene,
  exponentiate = TRUE,
  conf.int = TRUE
) %>%
  mutate(
    gene = gsub("_z$", "", term),
    gene = genes_to_model[match(gene, genes_model_safe)],
    HR_CI = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high),
    p_value_clean = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)),
    FDR = p.adjust(p.value, method = "BH"),
    FDR_clean = ifelse(FDR < 0.001, "<0.001", sprintf("%.3f", FDR)),
    significance = case_when(
      FDR < 0.05 & estimate > 1 ~ "FDR-significant higher risk",
      FDR < 0.05 & estimate < 1 ~ "FDR-significant protective",
      p.value < 0.05 ~ "Nominal only",
      TRUE ~ "Not significant"
    ),
    interpretation = case_when(
      estimate > 1 ~ "Higher VST-normalized expression associated with higher risk",
      estimate < 1 ~ "Higher VST-normalized expression associated with lower risk",
      TRUE ~ "Neutral"
    )
  )

write.csv(
  cox_table,
  file = file.path(outdir, "cluster3_THREE_GENE_cox_results.csv"),
  row.names = FALSE
)

cat("\nMultigene Cox results:\n")
print(cox_table)

#---------------------------#
# 16) PH assumption
#---------------------------#

cat("\n=== PROPORTIONAL HAZARDS TEST ===\n")

zph <- cox.zph(cox_multigene)

zph_table <- as.data.frame(zph$table) %>%
  rownames_to_column("term") %>%
  mutate(
    interpretation = case_when(
      p < 0.05 ~ "Possible PH violation",
      p >= 0.05 ~ "No evidence of PH violation"
    )
  )

write.csv(
  zph_table,
  file = file.path(outdir, "cluster3_THREE_GENE_PH_assumption.csv"),
  row.names = FALSE
)

cat("\nPH assumption results:\n")
print(zph_table)

# Global PH plot
png(
  filename = file.path(outdir, "cluster3_THREE_GENE_cox_zph_GLOBAL.png"),
  width = 2400,
  height = 2000,
  res = 300
)
plot(zph)
dev.off()

pdf(
  file = file.path(outdir, "cluster3_THREE_GENE_cox_zph_GLOBAL.pdf"),
  width = 9,
  height = 7
)
plot(zph)
dev.off()

# Individual PH plots per gene
zph_dir <- file.path(outdir, "PH_plots_individual_genes")
dir.create(zph_dir, recursive = TRUE, showWarnings = FALSE)

for (g in gene_z_terms) {
  
  gene_label <- genes_to_model[match(gsub("_z$", "", g), genes_model_safe)]
  
  png(
    filename = file.path(
      zph_dir,
      paste0("PH_plot_", gene_label, ".png")
    ),
    width = 1800,
    height = 1600,
    res = 300
  )
  plot(
    zph[g],
    main = paste0("Proportional hazards test: ", gene_label)
  )
  abline(h = 0, lty = 2, col = "gray40")
  dev.off()
  
  pdf(
    file = file.path(
      zph_dir,
      paste0("PH_plot_", gene_label, ".pdf")
    ),
    width = 7,
    height = 6
  )
  plot(
    zph[g],
    main = paste0("Proportional hazards test: ", gene_label)
  )
  abline(h = 0, lty = 2, col = "gray40")
  dev.off()
}

#---------------------------#
# 17) C-index
#---------------------------#

cat("\n=== C-INDEX ===\n")

cindex <- cox_summary$concordance[1]
cindex_se <- cox_summary$concordance[2]

cindex_table <- data.frame(
  model = "Three-gene Cox model",
  genes = paste(genes_to_model, collapse = " + "),
  c_index = cindex,
  c_index_se = cindex_se,
  n_samples = n_samples,
  n_events = n_events,
  n_predictors = n_predictors,
  events_per_variable = events_per_variable
)

write.csv(
  cindex_table,
  file = file.path(outdir, "cluster3_THREE_GENE_cindex.csv"),
  row.names = FALSE
)

print(cindex_table)

#---------------------------#
# 18) Risk score
#---------------------------#

cat("\n=== COMPUTING RISK SCORE ===\n")

model_df$risk_score <- predict(
  cox_multigene,
  type = "lp"
)

model_df$risk_group <- ifelse(
  model_df$risk_score >= median(model_df$risk_score, na.rm = TRUE),
  "High risk",
  "Low risk"
)

model_df$risk_group <- factor(
  model_df$risk_group,
  levels = c("Low risk", "High risk")
)

write.csv(
  model_df[, !sapply(model_df, is.list), drop = FALSE],
  file = file.path(outdir, "cluster3_THREE_GENE_model_dataframe_with_risk_score.csv"),
  row.names = FALSE
)

#---------------------------#
# 19) Forest plot
#---------------------------#

cat("\n=== GENERATING FOREST PLOT ===\n")

forest_df <- cox_table %>%
  mutate(
    term_clean = gene,
    status = case_when(
      p.value < 0.05 & estimate > 1 ~ "Increased Risk",
      p.value < 0.05 & estimate < 1 ~ "Protective",
      TRUE ~ "Not Significant"
    ),
    HR_CI = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high),
    p_value_clean = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
  )

forest_df$term_clean <- factor(
  forest_df$term_clean,
  levels = rev(forest_df$term_clean)
)

p_forest <- ggplot(
  forest_df,
  aes(x = estimate, y = term_clean)
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "darkgrey",
    linewidth = 0.8
  ) +
  geom_errorbarh(
    aes(
      xmin = conf.low,
      xmax = conf.high,
      color = status
    ),
    height = 0.2,
    linewidth = 0.8
  ) +
  geom_point(
    aes(color = status),
    size = 3.5
  ) +
  scale_x_log10(
    breaks = c(0.1, 0.5, 1, 2, 5, 10)
  ) +
  scale_color_manual(
    values = c(
      "Increased Risk" = "#D55E00",
      "Protective" = "#0072B2",
      "Not Significant" = "black"
    )
  ) +
  labs(
    title = "Cluster 3 three-gene Cox survival model",
    subtitle = paste0(
      "Gene-only model: ",
      paste(genes_to_model, collapse = " + ")
    ),
    x = "Hazard Ratio per 1 SD increase in VST-normalized expression",
    y = NULL
  ) +
  theme_bw(base_size = 16) + #here was 12
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "bold", size = 18), #HERE
    plot.margin = margin(5, 10, 5, 5)
  )

p_table <- ggplot(
  forest_df,
  aes(y = term_clean)
) +
  geom_text(
    aes(x = 0, label = HR_CI),
    size = 6.5, #here
    hjust = 0.5
  ) +
  geom_text(
    aes(x = 1.2, label = p_value_clean),
    size = 6.5, #here
    hjust = 0.5
  ) +
  annotate(
    "text",
    x = 0,
    y = Inf,
    label = "HR (95% CI)",
    fontface = "bold",
    vjust = 2,
    size = 6.5  #HERE
  ) +
  annotate(
    "text",
    x = 1.2,
    y = Inf,
    label = "p-value",
    fontface = "bold",
    vjust = 2,
    size = 6.5 #HERE
  ) +
  scale_x_continuous(
    limits = c(-0.8, 2)
  ) +
  theme_void() +
  theme(
    plot.margin = margin(l = 20, r = 20)
  )

combined_forest <- p_forest +
  p_table +
  patchwork::plot_layout(widths = c(2, 1))

ggsave(
  filename = file.path(outdir, "cluster3_THREE_GENE_forest_final.png"),
  plot = combined_forest,
  width = 14, #here was 12
  height = 5,
  dpi = 300
)

ggsave(
  filename = file.path(outdir, "cluster3_THREE_GENE_forest_final.pdf"),
  plot = combined_forest,
  width = 12,
  height = 5
)

write.csv(
  forest_df,
  file = file.path(outdir, "cluster3_THREE_GENE_forest_table_clean.csv"),
  row.names = FALSE
)

print(combined_forest)

#---------------------------#
# 20) KM plot by risk score
#---------------------------#

cat("\n=== GENERATING KM PLOT BY RISK SCORE ===\n")

fit_km_risk <- survfit(
  Surv(OS_days / 365.25, OS_event) ~ risk_group,
  data = model_df
)

km_risk_plot <- ggsurvplot(
  fit_km_risk,
  data = model_df,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = TRUE,
  xlim = c(0, 11),
  break.time.by = 2,
  title = "Cluster 3 overall survival by three-gene risk score",
  subtitle = paste0(
    "Genes: ",
    paste(genes_to_model, collapse = " + "),
    "; median split for visualization only"
  ),
  xlab = "Years",
  ylab = "Overall survival probability",
  legend.title = "Risk score group",
  legend.labs = c("Low risk", "High risk"),
  risk.table.height = 0.25,
  ggtheme = theme_bw(base_size = 16) #HERE
)

png(
  filename = file.path(outdir, "KM_cluster3_THREE_GENE_risk_score_median_split.png"),
  width = 2400,
  height = 2200,
  res = 300
)
print(km_risk_plot)
dev.off()

pdf(
  file = file.path(outdir, "KM_cluster3_THREE_GENE_risk_score_median_split.pdf"),
  width = 8,
  height = 7
)
print(km_risk_plot)
dev.off()


#---------------------------#
# 20B) Individual KM plots per gene
#---------------------------#

cat("\n=== GENERATING INDIVIDUAL KM PLOTS PER GENE ===\n")

km_gene_dir <- file.path(outdir, "KM_individual_genes_median_split")
dir.create(km_gene_dir, recursive = TRUE, showWarnings = FALSE)

for (gene in genes_to_model) {
  
  gene_safe <- make.names(gene)
  
  df_gene <- model_df %>%
    mutate(
      expr_value = .data[[gene_safe]],
      expr_group = ifelse(
        expr_value >= median(expr_value, na.rm = TRUE),
        "High expression",
        "Low expression"
      ),
      expr_group = factor(
        expr_group,
        levels = c("Low expression", "High expression")
      )
    )
  
  fit_km_gene <- survfit(
    Surv(OS_days, OS_event) ~ expr_group,
    data = df_gene
  )
  
  km_gene_plot <- ggsurvplot(
    fit_km_gene,
    data = df_gene,
    pval = TRUE,
    risk.table = TRUE,
    conf.int = TRUE,
    title = paste0("Cluster 3 overall survival by ", gene, " expression"),
    subtitle = "Median split; visualization only",
    xlab = "Days",
    ylab = "Overall survival probability",
    legend.title = gene,
    legend.labs = c("Low expression", "High expression"),
    risk.table.height = 0.25,
    ggtheme = theme_bw(base_size = 12)
  )
  
  png(
    filename = file.path(
      km_gene_dir,
      paste0("KM_cluster3_", gene, "_median_split.png")
    ),
    width = 2400,
    height = 2200,
    res = 300
  )
  print(km_gene_plot)
  dev.off()
  
  pdf(
    file = file.path(
      km_gene_dir,
      paste0("KM_cluster3_", gene, "_median_split.pdf")
    ),
    width = 8,
    height = 7
  )
  print(km_gene_plot)
  dev.off()
}


#---------------------------#
# 21) Save R objects
#---------------------------#

saveRDS(
  cox_multigene,
  file = file.path(outdir, "cluster3_THREE_GENE_cox_model.rds")
)

saveRDS(
  model_df,
  file = file.path(outdir, "cluster3_THREE_GENE_model_dataframe.rds")
)

#---------------------------#
# 22) TXT report
#---------------------------#

cat("\n=== GENERATING TXT REPORT ===\n")

report_file <- file.path(outdir, "cluster3_THREE_GENE_survival_report.txt")

sink(report_file)

cat("CLUSTER 3 THREE-GENE SURVIVAL MODEL REPORT\n\n")

cat("Dataset: TARGET-ALL-P2 PRIMARY U18\n")
cat("Analysis restricted to: cluster_k3 == 3\n")
cat("Model type: gene-only multigene Cox proportional hazards model\n\n")

cat("Genes requested:\n")
print(target_genes)
cat("\n")

cat("Genes used in final model:\n")
print(genes_to_model)
cat("\n\n")

cat("Survival definition:\n")
cat("- OS_days = days_to_death if Dead\n")
cat("- OS_days = days_to_last_follow_up if Alive\n")
cat("- OS_event: Dead = 1, Alive = 0\n\n")

cat("Samples used:", n_samples, "\n")
cat("Events:", n_events, "\n")
cat("Censored:", sum(model_df$OS_event == 0), "\n\n")

cat("Expression filtering:\n")
cat("- Minimum raw count:", min_count, "\n")
cat("- Minimum proportion of samples:", min_prop_samples, "\n")
cat("- Minimum number of samples:", min_samples, "\n\n")

cat("Statistical model:\n")
cat(deparse(formula_multigene), "\n\n")

cat("Important modeling note:\n")
cat("- Age and year of diagnosis were intentionally not included.\n")
cat("- The analysis is restricted to Cluster 3 and pediatric/U18 samples.\n")
cat("- The goal is to evaluate a focused transcriptomic risk model inside the poor-survival cluster.\n\n")

cat("Statistical power:\n")
cat("- Predictors:", n_predictors, "\n")
cat("- Events per variable:", round(events_per_variable, 2), "\n")
if (events_per_variable < 10) {
  cat("- WARNING: EPV < 10. Interpret cautiously.\n")
}
cat("\n")

cat("Multigene Cox results:\n")
print(cox_table)
cat("\n\n")

cat("Proportional hazards assumption:\n")
print(zph_table)
cat("\n\n")

cat("C-index:\n")
print(cindex_table)
cat("\n\n")

cat("Interpretation guide:\n")
cat("- HR > 1: higher VST-normalized expression is associated with increased risk of death, adjusted for the other genes in the model.\n")
cat("- HR < 1: higher VST-normalized expression is associated with reduced risk of death, adjusted for the other genes in the model.\n")
cat("- HR is interpreted per 1 SD increase in VST-normalized expression.\n")
cat("- FDR is corrected only across the selected genes, not genome-wide.\n")
cat("- Risk score is the Cox linear predictor from the three-gene model.\n")
cat("- KM risk groups use median split only for visualization.\n")
cat("- Cox model with continuous standardized expression is the primary statistical result.\n")

sink()

#---------------------------#
# 23) Final message
#---------------------------#

cat("\n=== CLUSTER 3 THREE-GENE SURVIVAL ANALYSIS COMPLETED ===\n")
cat("Results saved in:\n")
cat(outdir, "\n")

cat("\nGenes used:\n")
print(genes_to_model)

cat("\nCox results:\n")
print(cox_table)

cat("\nC-index:\n")
print(cindex_table)

cat("\nPH assumption:\n")
print(zph_table)