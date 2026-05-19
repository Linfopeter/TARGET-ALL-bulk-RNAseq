############################################################
### 15_SURVIVAL_ANALYSIS_CLUSTER3_GENOMEWIDE_BY_SAMPLE_TYPE.R
### GOAL:
###   Genome-wide survival screening restricted to Cluster 3
###   separated by sample type:
###     - bone_marrow
###     - PBMCs
###
### MODEL PER GENE:
###   Surv(OS_days, OS_event) ~ expr_z
###
### PURPOSE:
###   Exploratory genome-wide gene-only Cox screening to identify
###   candidate survival-associated genes for final focused models.
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

cat("=== SCRIPT 15: CLUSTER 3 GENOME-WIDE GENE-ONLY SURVIVAL SCREENING BY SAMPLE TYPE ===\n")

#---------------------------#
# 1) Settings
#---------------------------#

top_n_forest <- 25
top_n_km <- 10
top_n_ph <- 25

min_count <- 10
min_prop_samples <- 0.20

#---------------------------#
# 2) Paths
#---------------------------#

input_counts_file <- "subsets/primary_u18/final_analysis/TARGET_ALL_P2_PRIMARY_U18_FINAL_FOR_ANALYSIS.RData"
input_meta_clusters <- "subsets/primary_u18/final_analysis/06_qc_normalization_deseq2/meta_with_clusters.rds"

outdir <- "subsets/primary_u18/final_analysis/15_survival_cluster3_genomewide_by_sample_type_gene_only"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

#---------------------------#
# 3) Load data
#---------------------------#

cat("\n=== LOADING DATA ===\n")

load(input_counts_file)
meta_clusters <- readRDS(input_meta_clusters)

cat("- Raw rows in counts_final:", nrow(counts_final), "\n")
cat("- Raw samples in counts_final:", ncol(counts_final), "\n")

#---------------------------#
# 4) Clean count matrix
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
# 5) Align metadata
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
# 6) Define sample groups
#---------------------------#

meta_clusters <- meta_clusters %>%
  mutate(
    sample_group = case_when(
      sample_type.x == "Primary Blood Derived Cancer - Bone Marrow" ~ "bone_marrow",
      sample_type.x == "Primary Blood Derived Cancer - Peripheral Blood" ~ "PBMCs",
      TRUE ~ NA_character_
    )
  )

cat("\nSample groups in metadata:\n")
print(table(meta_clusters$sample_group, useNA = "ifany"))

#---------------------------#
# 7) Global annotation
#---------------------------#

cat("\n=== ANNOTATING ENSEMBL IDS TO HGNC SYMBOLS ===\n")

ensembl_ids <- rownames(counts_final)

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
  file = file.path(outdir, "gene_annotation_table_global.csv"),
  row.names = FALSE
)

cat("- Annotated genes with HGNC symbols:", sum(!is.na(annotation_df$gene_symbol)), "\n")

#---------------------------#
# 8) Function
#---------------------------#

run_cluster3_genomewide_survival <- function(sample_group_name) {
  
  cat("\n============================================================\n")
  cat("RUNNING CLUSTER 3 GENOME-WIDE GENE-ONLY SURVIVAL SCREENING FOR:", sample_group_name, "\n")
  cat("============================================================\n")
  
  outdir_group <- file.path(outdir, sample_group_name)
  dir.create(outdir_group, recursive = TRUE, showWarnings = FALSE)
  
  #---------------------------#
  # Survival metadata
  #---------------------------#
  
  meta_surv <- meta_clusters %>%
    mutate(
      OS_days = case_when(
        vital_status == "Dead" ~ suppressWarnings(as.numeric(days_to_death)),
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
      sample_group == sample_group_name,
      !is.na(OS_days),
      OS_days > 0,
      !is.na(OS_event)
    )
  
  cat("\nSamples in Cluster 3 -", sample_group_name, ":", nrow(meta_surv), "\n")
  cat("\nEvents:\n")
  print(table(meta_surv$OS_event, useNA = "ifany"))
  
  if (nrow(meta_surv) < 20) {
    cat("\nWARNING: fewer than 20 samples. Skipping:", sample_group_name, "\n")
    return(NULL)
  }
  
  if (sum(meta_surv$OS_event == 1) < 5) {
    cat("\nWARNING: fewer than 5 events. Skipping:", sample_group_name, "\n")
    return(NULL)
  }
  
  write.csv(
    meta_surv[, !sapply(meta_surv, is.list), drop = FALSE],
    file = file.path(outdir_group, paste0("cluster3_", sample_group_name, "_survival_metadata_used.csv")),
    row.names = FALSE
  )
  
  saveRDS(
    meta_surv,
    file = file.path(outdir_group, paste0("cluster3_", sample_group_name, "_survival_metadata_used.rds"))
  )
  
  #---------------------------#
  # Subset counts
  #---------------------------#
  
  counts_cluster3 <- counts_final[
    ,
    meta_surv$counts_colname,
    drop = FALSE
  ]
  
  stopifnot(all(colnames(counts_cluster3) == meta_surv$counts_colname))
  
  #---------------------------#
  # Expression filtering
  #---------------------------#
  
  min_samples <- ceiling(min_prop_samples * ncol(counts_cluster3))
  
  keep_expr <- rowSums(counts_cluster3 >= min_count, na.rm = TRUE) >= min_samples
  
  counts_cluster3_filtered <- counts_cluster3[
    keep_expr,
    ,
    drop = FALSE
  ]
  
  cat("\nGenes before filtering:", nrow(counts_cluster3), "\n")
  cat("Genes after filtering:", nrow(counts_cluster3_filtered), "\n")
  cat("Minimum samples required with count >=", min_count, ":", min_samples, "\n")
  
  #---------------------------#
  # VST
  #---------------------------#
  
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
  # Annotate filtered matrix
  #---------------------------#
  
  annotation_filtered <- annotation_df %>%
    filter(ensembl_gene_id %in% rownames(vst_mat)) %>%
    distinct(ensembl_gene_id, .keep_all = TRUE)
  
  gene_names_final <- annotation_filtered$gene_label
  names(gene_names_final) <- annotation_filtered$ensembl_gene_id
  
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
  # Power summary
  #---------------------------#
  
  n_samples <- nrow(meta_surv)
  n_events <- sum(meta_surv$OS_event == 1)
  n_predictors_per_model <- 1
  events_per_variable <- n_events / n_predictors_per_model
  
  cat("\n=== STATISTICAL POWER SUMMARY ===\n")
  cat("- Samples:", n_samples, "\n")
  cat("- Events:", n_events, "\n")
  cat("- Censored:", sum(meta_surv$OS_event == 0), "\n")
  cat("- Predictors per Cox model: expr_z =", n_predictors_per_model, "\n")
  cat("- Events per variable:", round(events_per_variable, 2), "\n")
  
  #---------------------------#
  # Genome-wide Cox analysis
  #---------------------------#
  
  cat("\n=== RUNNING GENOME-WIDE GENE-ONLY COX MODELS ===\n")
  
  genes_to_test <- rownames(vst_mat)
  
  cat("- Genes to test:", length(genes_to_test), "\n")
  
  cox_results <- list()
  cindex_results <- list()
  zph_results_top <- list()
  
  for (i in seq_along(genes_to_test)) {
    
    gene <- genes_to_test[i]
    
    if (i %% 1000 == 0) {
      cat("Processed", i, "of", length(genes_to_test), "genes\n")
    }
    
    gene_expr <- as.numeric(vst_mat[gene, meta_surv$counts_colname])
    
    df_gene <- meta_surv %>%
      mutate(
        gene = gene,
        expr = gene_expr,
        expr_z = as.numeric(scale(expr))
      ) %>%
      filter(
        !is.na(expr_z),
        is.finite(expr_z)
      )
    
    if (nrow(df_gene) < 20) next
    if (sd(df_gene$expr_z, na.rm = TRUE) == 0) next
    
    cox_gene <- tryCatch(
      suppressWarnings(
        coxph(
          Surv(OS_days, OS_event) ~ expr_z,
          data = df_gene
        )
      ),
      error = function(e) NULL
    )
    
    if (is.null(cox_gene)) next
    
    tidy_gene <- tryCatch(
      broom::tidy(
        cox_gene,
        exponentiate = TRUE,
        conf.int = TRUE
      ) %>%
        filter(term == "expr_z") %>%
        mutate(
          gene = gene,
          sample_group = sample_group_name,
          model = "Gene-only Cox model",
          interpretation = ifelse(
            estimate > 1,
            "Higher VST-normalized expression associated with higher risk",
            "Higher VST-normalized expression associated with lower risk"
          )
        ),
      error = function(e) NULL
    )
    
    if (is.null(tidy_gene)) next
    if (nrow(tidy_gene) == 0) next
    
    cox_results[[gene]] <- tidy_gene
    
    gene_cindex <- tryCatch(
      summary(cox_gene)$concordance[1],
      error = function(e) NA_real_
    )
    
    cindex_results[[gene]] <- data.frame(
      gene = gene,
      sample_group = sample_group_name,
      gene_model_cindex = gene_cindex
    )
  }
  
  if (length(cox_results) == 0) {
    cat("\nNo Cox models were successfully fitted for:", sample_group_name, "\n")
    return(NULL)
  }
  
  #---------------------------#
  # Combine results
  #---------------------------#
  
  cat("\n=== COMBINING GENOME-WIDE RESULTS ===\n")
  
  cox_table <- bind_rows(cox_results) %>%
    mutate(
      FDR = p.adjust(p.value, method = "BH"),
      HR_CI = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high),
      p_value_clean = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)),
      FDR_clean = ifelse(FDR < 0.001, "<0.001", sprintf("%.3f", FDR)),
      significance = case_when(
        FDR < 0.05 & estimate > 1 ~ "FDR-significant higher risk",
        FDR < 0.05 & estimate < 1 ~ "FDR-significant protective",
        p.value < 0.05 ~ "Nominal only",
        TRUE ~ "Not significant"
      )
    ) %>%
    arrange(FDR, p.value)
  
  cindex_table <- bind_rows(cindex_results) %>%
    arrange(desc(gene_model_cindex))
  
  significant_genes <- cox_table %>%
    filter(FDR < 0.05)
  
  nominal_genes <- cox_table %>%
    filter(p.value < 0.05)
  
  top100_genes <- cox_table %>%
    slice_head(n = 100)
  
  write.csv(
    cox_table,
    file = file.path(outdir_group, paste0("cluster3_", sample_group_name, "_GENOMEWIDE_gene_only_cox_all_genes.csv")),
    row.names = FALSE
  )
  
  write.csv(
    significant_genes,
    file = file.path(outdir_group, paste0("cluster3_", sample_group_name, "_GENOMEWIDE_gene_only_FDR_significant_genes.csv")),
    row.names = FALSE
  )
  
  write.csv(
    nominal_genes,
    file = file.path(outdir_group, paste0("cluster3_", sample_group_name, "_GENOMEWIDE_gene_only_nominal_p005_genes.csv")),
    row.names = FALSE
  )
  
  write.csv(
    top100_genes,
    file = file.path(outdir_group, paste0("cluster3_", sample_group_name, "_GENOMEWIDE_gene_only_top100_genes.csv")),
    row.names = FALSE
  )
  
  write.csv(
    cindex_table,
    file = file.path(outdir_group, paste0("cluster3_", sample_group_name, "_GENOMEWIDE_gene_only_cindex.csv")),
    row.names = FALSE
  )
  
  cat("\nTop 20 survival-associated genes:\n")
  print(head(cox_table, 20))
  
  cat("\nFDR-significant genes:", nrow(significant_genes), "\n")
  cat("Nominal p < 0.05 genes:", nrow(nominal_genes), "\n")
  
  #---------------------------#
  # PH assumption for top genes
  #---------------------------#
  
  cat("\n=== TESTING PH ASSUMPTION FOR TOP GENES ===\n")
  
  ph_dir <- file.path(outdir_group, "PH_plots_top_genes")
  dir.create(ph_dir, recursive = TRUE, showWarnings = FALSE)
  
  top_ph_genes <- cox_table %>%
    slice_head(n = min(top_n_ph, nrow(cox_table))) %>%
    pull(gene)
  
  for (gene in top_ph_genes) {
    
    gene_expr <- as.numeric(vst_mat[gene, meta_surv$counts_colname])
    
    df_gene <- meta_surv %>%
      mutate(
        expr = gene_expr,
        expr_z = as.numeric(scale(expr))
      ) %>%
      filter(
        !is.na(expr_z),
        is.finite(expr_z)
      )
    
    cox_gene <- tryCatch(
      suppressWarnings(
        coxph(
          Surv(OS_days, OS_event) ~ expr_z,
          data = df_gene
        )
      ),
      error = function(e) NULL
    )
    
    if (is.null(cox_gene)) next
    
    zph <- tryCatch(
      cox.zph(cox_gene),
      error = function(e) NULL
    )
    
    if (!is.null(zph)) {
      
      zph_table <- as.data.frame(zph$table)
      zph_table$term <- rownames(zph_table)
      zph_table$gene <- gene
      zph_table$sample_group <- sample_group_name
      rownames(zph_table) <- NULL
      zph_results_top[[gene]] <- zph_table
      
      png(
        filename = file.path(ph_dir, paste0("PH_plot_", gene, ".png")),
        width = 1800,
        height = 1600,
        res = 300
      )
      plot(
        zph["expr_z"],
        main = paste0("Proportional hazards test: ", gene)
      )
      abline(h = 0, lty = 2, col = "gray40")
      dev.off()
      
      pdf(
        file = file.path(ph_dir, paste0("PH_plot_", gene, ".pdf")),
        width = 7,
        height = 6
      )
      plot(
        zph["expr_z"],
        main = paste0("Proportional hazards test: ", gene)
      )
      abline(h = 0, lty = 2, col = "gray40")
      dev.off()
    }
  }
  
  zph_table_top <- bind_rows(zph_results_top)
  
  write.csv(
    zph_table_top,
    file = file.path(outdir_group, paste0("cluster3_", sample_group_name, "_GENOMEWIDE_gene_only_top", top_n_ph, "_PH_assumption.csv")),
    row.names = FALSE
  )
  
  #---------------------------#
  # Forest plot top genes
  #---------------------------#
  
  cat("\n=== GENERATING FOREST PLOT ===\n")
  
  forest_df <- cox_table %>%
    slice_head(n = min(top_n_forest, nrow(cox_table))) %>%
    mutate(
      term_clean = gene,
      status = case_when(
        p.value < 0.05 & estimate > 1 ~ "Increased Risk",
        p.value < 0.05 & estimate < 1 ~ "Protective",
        TRUE ~ "Not Significant"
      )
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
    geom_errorbar(
      aes(
        xmin = conf.low,
        xmax = conf.high,
        color = status
      ),
      orientation = "y",
      width = 0.2,
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
      title = paste0("Top survival-associated genes in Cluster 3: ", sample_group_name),
      subtitle = "Genome-wide gene-only Cox screening",
      x = "Hazard Ratio per 1 SD increase in VST-normalized expression",
      y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(face = "bold", size = 9),
      plot.margin = margin(5, 10, 5, 5)
    )
  
  p_table <- ggplot(
    forest_df,
    aes(y = term_clean)
  ) +
    geom_text(
      aes(x = 0, label = HR_CI),
      size = 3.3,
      hjust = 0.5
    ) +
    geom_text(
      aes(x = 1.2, label = p_value_clean),
      size = 3.3,
      hjust = 0.5
    ) +
    geom_text(
      aes(x = 2.2, label = FDR_clean),
      size = 3.3,
      hjust = 0.5
    ) +
    annotate(
      "text",
      x = 0,
      y = Inf,
      label = "HR (95% CI)",
      fontface = "bold",
      vjust = 2,
      size = 3.8
    ) +
    annotate(
      "text",
      x = 1.2,
      y = Inf,
      label = "p-value",
      fontface = "bold",
      vjust = 2,
      size = 3.8
    ) +
    annotate(
      "text",
      x = 2.2,
      y = Inf,
      label = "FDR",
      fontface = "bold",
      vjust = 2,
      size = 3.8
    ) +
    scale_x_continuous(
      limits = c(-0.8, 3)
    ) +
    theme_void() +
    theme(
      plot.margin = margin(l = 20, r = 20)
    )
  
  combined_forest <- p_forest +
    p_table +
    patchwork::plot_layout(widths = c(2, 1.2))
  
  ggsave(
    filename = file.path(outdir_group, paste0("cluster3_", sample_group_name, "_GENOMEWIDE_gene_only_top_genes_forest_final.png")),
    plot = combined_forest,
    width = 14,
    height = max(6, nrow(forest_df) * 0.35),
    dpi = 300
  )
  
  ggsave(
    filename = file.path(outdir_group, paste0("cluster3_", sample_group_name, "_GENOMEWIDE_gene_only_top_genes_forest_final.pdf")),
    plot = combined_forest,
    width = 14,
    height = max(6, nrow(forest_df) * 0.35)
  )
  
  write.csv(
    forest_df,
    file = file.path(outdir_group, paste0("cluster3_", sample_group_name, "_GENOMEWIDE_gene_only_forest_table_clean.csv")),
    row.names = FALSE
  )
  
  #---------------------------#
  # KM plots for top genes
  #---------------------------#
  
  cat("\n=== GENERATING KM PLOTS FOR TOP GENES ===\n")
  
  km_dir <- file.path(outdir_group, "KM_GENOMEWIDE_top_genes_median_split")
  dir.create(km_dir, recursive = TRUE, showWarnings = FALSE)
  
  km_genes <- cox_table %>%
    slice_head(n = min(top_n_km, nrow(cox_table))) %>%
    pull(gene)
  
  for (gene in km_genes) {
    
    gene_expr <- as.numeric(vst_mat[gene, meta_surv$counts_colname])
    
    df_gene <- meta_surv %>%
      mutate(
        expr = gene_expr,
        expr_group = ifelse(
          expr >= median(expr, na.rm = TRUE),
          "High expression",
          "Low expression"
        ),
        expr_group = factor(
          expr_group,
          levels = c("Low expression", "High expression")
        )
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
      title = paste0("Cluster 3 ", sample_group_name, " OS by ", gene, " expression"),
      subtitle = "Median split; visualization only. Cox continuous model is primary.",
      xlab = "Days",
      ylab = "Overall survival probability",
      legend.title = gene,
      legend.labs = c("Low expression", "High expression"),
      risk.table.height = 0.25,
      ggtheme = theme_bw(base_size = 12)
    )
    
    png(
      filename = file.path(km_dir, paste0("KM_cluster3_", sample_group_name, "_", gene, "_median_split.png")),
      width = 2400,
      height = 2200,
      res = 300
    )
    print(km_plot)
    dev.off()
    
    pdf(
      file = file.path(km_dir, paste0("KM_cluster3_", sample_group_name, "_", gene, "_median_split.pdf")),
      width = 8,
      height = 7
    )
    print(km_plot)
    dev.off()
  }
  
  #---------------------------#
  # TXT report
  #---------------------------#
  
  cat("\n=== GENERATING TXT REPORT ===\n")
  
  report_file <- file.path(outdir_group, paste0("cluster3_", sample_group_name, "_GENOMEWIDE_gene_only_survival_report.txt"))
  
  sink(report_file)
  
  cat("CLUSTER 3 GENOME-WIDE GENE-ONLY SURVIVAL SCREENING REPORT\n\n")
  
  cat("Dataset: TARGET-ALL-P2 PRIMARY U18\n")
  cat("Analysis restricted to: cluster_k3 == 3\n")
  cat("Sample group:", sample_group_name, "\n\n")
  
  cat("Original sample_type.x values:\n")
  print(unique(meta_surv$sample_type.x))
  cat("\n")
  
  cat("Survival definition:\n")
  cat("- OS_days = days_to_death if Dead\n")
  cat("- OS_days = days_to_last_follow_up if Alive\n")
  cat("- OS_event: Dead = 1, Alive = 0\n\n")
  
  cat("Samples used:", nrow(meta_surv), "\n")
  cat("Events:", sum(meta_surv$OS_event == 1), "\n")
  cat("Censored:", sum(meta_surv$OS_event == 0), "\n\n")
  
  cat("Expression filtering:\n")
  cat("- Minimum raw count:", min_count, "\n")
  cat("- Minimum proportion of samples:", min_prop_samples, "\n")
  cat("- Minimum number of samples:", min_samples, "\n")
  cat("- Genes tested after filtering:", length(genes_to_test), "\n\n")
  
  cat("Statistical model:\n")
  cat("- Cox model per gene: Surv(OS_days, OS_event) ~ expr_z\n")
  cat("- expr_z = scaled VST-normalized expression\n")
  cat("- HR is interpreted per 1 SD increase in VST-normalized expression\n")
  cat("- Age and year were intentionally not included.\n\n")
  
  cat("Statistical power:\n")
  cat("- Cox model predictors per gene: expr_z only\n")
  cat("- Number of predictors per model:", n_predictors_per_model, "\n")
  cat("- Events per variable:", round(events_per_variable, 2), "\n\n")
  
  cat("Genome-wide results summary:\n")
  cat("- Total genes tested:", nrow(cox_table), "\n")
  cat("- FDR-significant genes:", nrow(significant_genes), "\n")
  cat("- Nominal p < 0.05 genes:", nrow(nominal_genes), "\n\n")
  
  cat("Top 50 survival-associated genes:\n")
  print(head(cox_table, 50))
  cat("\n\n")
  
  cat("Top C-index genes:\n")
  print(head(cindex_table, 50))
  cat("\n\n")
  
  cat("PH assumption for top genes:\n")
  print(zph_table_top)
  cat("\n\n")
  
  cat("Interpretation guide:\n")
  cat("- HR > 1: higher VST-normalized gene expression is associated with increased risk of death.\n")
  cat("- HR < 1: higher VST-normalized gene expression is associated with reduced risk of death.\n")
  cat("- FDR is Benjamini-Hochberg corrected across all tested genes within this sample group.\n")
  cat("- KM plots use median split only for visualization.\n")
  cat("- Cox models use continuous standardized expression and are the primary statistical result.\n")
  cat("- This script is an exploratory genome-wide screening step; final focused models should be built separately.\n")
  cat("- If no FDR-significant genes are found, nominal p < 0.05 genes should be reported only as exploratory.\n")
  
  sink()
  
  cat("\nCompleted:", sample_group_name, "\n")
  
  return(list(
    sample_group = sample_group_name,
    meta_surv = meta_surv,
    cox_table = cox_table,
    significant_genes = significant_genes,
    nominal_genes = nominal_genes,
    cindex_table = cindex_table,
    zph_table_top = zph_table_top
  ))
}

#---------------------------#
# 9) Run both groups
#---------------------------#

results_by_sample_group <- list()

for (sg in c("bone_marrow", "PBMCs")) {
  results_by_sample_group[[sg]] <- run_cluster3_genomewide_survival(sg)
}

saveRDS(
  results_by_sample_group,
  file = file.path(outdir, "results_by_sample_group_GENOMEWIDE_gene_only.rds")
)

cat("\n=== CLUSTER 3 GENOME-WIDE SAMPLE-TYPE STRATIFIED GENE-ONLY SURVIVAL SCREENING COMPLETED ===\n")
cat("Results saved in:\n")
cat(outdir, "\n")