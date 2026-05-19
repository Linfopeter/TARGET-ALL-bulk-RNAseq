### 12_survival_analysis_PRIMARY_U18_TARGET_ALL_P2.R
### GOAL:
###   Kaplan-Meier + Cox survival analysis by cluster_k3
### OUTPUT:
###   - KM curve PNG/PDF
###   - Cox univariate table
###   - Cox multivariate table
###   - cox.zph proportional hazards tests
###   - TXT report

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(survminer)
  library(broom)
})

cat("=== SCRIPT 12: SURVIVAL ANALYSIS FINAL ===\n")

#---------------------------#
# 1) Paths                  #
#---------------------------#
input_counts_file <- "subsets/primary_u18/final_analysis/TARGET_ALL_P2_PRIMARY_U18_FINAL_FOR_ANALYSIS.RData"
input_meta_clusters <- "subsets/primary_u18/final_analysis/06_qc_normalization_deseq2/meta_with_clusters.rds"

outdir <- "subsets/primary_u18/final_analysis/12_survival_analysis"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

#---------------------------#
# 2) Load data              #
#---------------------------#
cat("\n=== CARGANDO DATOS ===\n")

load(input_counts_file)
meta_clusters <- readRDS(input_meta_clusters)

meta_clusters <- meta_clusters[
  match(colnames(counts_final), meta_clusters$counts_colname),
  ,
  drop = FALSE
]

stopifnot(all(colnames(counts_final) == meta_clusters$counts_colname))

cat("- Muestras:", nrow(meta_clusters), "\n")
cat("- Variables:", ncol(meta_clusters), "\n")

#---------------------------#
# 3) Build survival dataset #
#---------------------------#
cat("\n=== CONSTRUYENDO DATASET DE SUPERVIVENCIA ===\n")

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
    !is.na(OS_days),
    OS_days > 0,
    !is.na(OS_event),
    !is.na(cluster_k3)
  )

meta_surv$cluster_k3 <- relevel(meta_surv$cluster_k3, ref = "1")

cat("\nMuestras usadas en supervivencia:", nrow(meta_surv), "\n")

cat("\nEventos OS:\n")
print(table(meta_surv$OS_event, useNA = "ifany"))

cat("\nClusters:\n")
print(table(meta_surv$cluster_k3, useNA = "ifany"))

cat("\nResumen OS_days:\n")
print(summary(meta_surv$OS_days))

# Guardar metadata completa en RDS
saveRDS(
  meta_surv,
  file = file.path(outdir, "survival_metadata_used_full.rds")
)

# Guardar metadata plana en CSV
meta_surv_csv <- meta_surv[, !sapply(meta_surv, is.list), drop = FALSE]

write.csv(
  meta_surv_csv,
  file = file.path(outdir, "survival_metadata_used.csv"),
  row.names = FALSE
)
#---------------------------#
# 4) Kaplan-Meier           #
#---------------------------#
cat("\n=== KAPLAN-MEIER ===\n")

# Convertimos los días a años antes de crear el objeto Surv
meta_surv$OS_years <- meta_surv$OS_days / 365.25

# Usamos OS_years en lugar de OS_days
surv_obj <- Surv(time = meta_surv$OS_years, event = meta_surv$OS_event)

fit_km <- survfit(surv_obj ~ cluster_k3, data = meta_surv)

logrank <- survdiff(surv_obj ~ cluster_k3, data = meta_surv)
logrank_p <- 1 - pchisq(logrank$chisq, df = length(logrank$n) - 1)

km_plot <- ggsurvplot(
  fit_km,
  data = meta_surv,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = TRUE,
  title = "Overall Survival by cluster_k3",
  xlab = "Years", # Cambiado de Days a Years
  ylab = "Overall survival probability",
  legend.title = "Cluster",
  legend.labs = paste("Cluster", levels(meta_surv$cluster_k3)),
  risk.table.height = 0.25,
  ggtheme = theme_bw(base_size = 12)
)

# Guardar PNG
png(
  filename = file.path(outdir, "kaplan_meier_OS_by_cluster_k3.png"),
  width = 2400,
  height = 2200,
  res = 300
)
print(km_plot)
dev.off()

# Guardar PDF
pdf(
  file = file.path(outdir, "kaplan_meier_OS_by_cluster_k3.pdf"),
  width = 8,
  height = 7
)
print(km_plot)
dev.off()

#---------------------------#
# 5) Cox univariate         #
#---------------------------#
cat("\n=== COX UNIVARIADO ===\n")

cox_uni <- coxph(surv_obj ~ cluster_k3, data = meta_surv)
cox_uni_summary <- summary(cox_uni)

cox_uni_table <- broom::tidy(
  cox_uni,
  exponentiate = TRUE,
  conf.int = TRUE
)

write.csv(
  cox_uni_table,
  file = file.path(outdir, "cox_univariate_cluster_k3.csv"),
  row.names = FALSE
)

#---------------------------#
# 6) Cox multivariate       #
#---------------------------#
cat("\n=== COX MULTIVARIADO ===\n")

cox_multi <- coxph(
  surv_obj ~ cluster_k3 + age + year,
  data = meta_surv
)

cox_multi_summary <- summary(cox_multi)

cox_multi_table <- broom::tidy(
  cox_multi,
  exponentiate = TRUE,
  conf.int = TRUE
)

write.csv(
  cox_multi_table,
  file = file.path(outdir, "cox_multivariate_cluster_k3_age_year.csv"),
  row.names = FALSE
)

#---------------------------#
# 7) Proportional hazards   #
#---------------------------#
cat("\n=== TEST DE PROPORCIONALIDAD DE RIESGOS ===\n")

zph_uni <- cox.zph(cox_uni)
zph_multi <- cox.zph(cox_multi)

zph_uni_table <- as.data.frame(zph_uni$table)
zph_multi_table <- as.data.frame(zph_multi$table)

write.csv(
  zph_uni_table,
  file = file.path(outdir, "cox_zph_univariate.csv"),
  row.names = TRUE
)

write.csv(
  zph_multi_table,
  file = file.path(outdir, "cox_zph_multivariate.csv"),
  row.names = TRUE
)

png(
  filename = file.path(outdir, "cox_zph_univariate.png"),
  width = 2000,
  height = 1600,
  res = 300
)
plot(zph_uni)
dev.off()

png(
  filename = file.path(outdir, "cox_zph_multivariate.png"),
  width = 2400,
  height = 2000,
  res = 300
)
plot(zph_multi)
dev.off()

#---------------------------#
# 8) Save R objects         #
#---------------------------#
saveRDS(
  fit_km,
  file = file.path(outdir, "survfit_km_cluster_k3.rds")
)

saveRDS(
  cox_uni,
  file = file.path(outdir, "cox_univariate_cluster_k3.rds")
)

saveRDS(
  cox_multi,
  file = file.path(outdir, "cox_multivariate_cluster_k3_age_year.rds")
)

#---------------------------#
# 9) TXT report             #
#---------------------------#
cat("\n=== GENERANDO REPORTE TXT ===\n")

report_file <- file.path(outdir, "survival_analysis_report.txt")

sink(report_file)

cat("SURVIVAL ANALYSIS REPORT\n")
cat("Dataset: TARGET-ALL-P2 PRIMARY U18\n")
cat("Analysis: Overall Survival by cluster_k3\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("INPUT FILES\n")
cat("- ", input_counts_file, "\n", sep = "")
cat("- ", input_meta_clusters, "\n\n", sep = "")

cat("SURVIVAL DEFINITION\n")
cat("- Time variable: days_to_last_follow_up\n")
cat("- Event variable: vital_status\n")
cat("- Dead = 1\n")
cat("- Alive = 0\n")
cat("- Unknown / missing removed\n\n")

cat("SAMPLE SUMMARY\n")
cat("- Samples used:", nrow(meta_surv), "\n")
cat("- Events:", sum(meta_surv$OS_event == 1), "\n")
cat("- Censored:", sum(meta_surv$OS_event == 0), "\n\n")

cat("CLUSTER DISTRIBUTION\n")
print(table(meta_surv$cluster_k3, useNA = "ifany"))
cat("\n")

cat("EVENTS BY CLUSTER\n")
print(table(meta_surv$cluster_k3, meta_surv$OS_event))
cat("\n")

cat("OS DAYS SUMMARY\n")
print(summary(meta_surv$OS_days))
cat("\n")

cat("LOG-RANK TEST\n")
print(logrank)
cat("\nApproximate log-rank p-value:", signif(logrank_p, 4), "\n\n")

cat("COX UNIVARIATE MODEL\n")
cat("Formula: Surv(OS_days, OS_event) ~ cluster_k3\n\n")
print(cox_uni_summary)
cat("\n")

cat("COX UNIVARIATE TABLE\n")
print(cox_uni_table)
cat("\n\n")

cat("COX MULTIVARIATE MODEL\n")
cat("Formula: Surv(OS_days, OS_event) ~ cluster_k3 + age + year\n\n")
print(cox_multi_summary)
cat("\n")

cat("COX MULTIVARIATE TABLE\n")
print(cox_multi_table)
cat("\n\n")

cat("PROPORTIONAL HAZARDS TEST - UNIVARIATE\n")
print(zph_uni)
cat("\n")

cat("PROPORTIONAL HAZARDS TEST - MULTIVARIATE\n")
print(zph_multi)
cat("\n")

cat("BASIC INTERPRETATION\n")
cat("Kaplan-Meier and log-rank analyses evaluate whether overall survival differs across cluster_k3 groups.\n")
cat("Cox univariate analysis estimates the hazard ratio of each cluster compared with cluster 1.\n")
cat("Cox multivariate analysis adjusts this association by age at diagnosis and year of diagnosis.\n")
cat("A hazard ratio greater than 1 indicates increased risk of death, while a hazard ratio below 1 indicates reduced risk.\n")
cat("The cox.zph test evaluates whether the proportional hazards assumption is violated.\n")
cat("For cox.zph, p < 0.05 suggests possible violation of the proportional hazards assumption.\n\n")

cat("FILES GENERATED\n")
cat("- kaplan_meier_OS_by_cluster_k3.png\n")
cat("- kaplan_meier_OS_by_cluster_k3.pdf\n")
cat("- survival_metadata_used.csv\n")
cat("- cox_univariate_cluster_k3.csv\n")
cat("- cox_multivariate_cluster_k3_age_year.csv\n")
cat("- cox_zph_univariate.csv\n")
cat("- cox_zph_multivariate.csv\n")
cat("- cox_zph_univariate.png\n")
cat("- cox_zph_multivariate.png\n")
cat("- survfit_km_cluster_k3.rds\n")
cat("- cox_univariate_cluster_k3.rds\n")
cat("- cox_multivariate_cluster_k3_age_year.rds\n")
cat("- survival_analysis_report.txt\n")

sink()

#---------------------------#
# 10) Final message         #
#---------------------------#
cat("\n=== SCRIPT 12 COMPLETADO ===\n")
cat("\nArchivos guardados en:\n")
cat(outdir, "\n")

cat("\nPrincipales resultados:\n")
cat("- Log-rank p-value:", signif(logrank_p, 4), "\n")

cat("\nCox univariado:\n")
print(cox_uni_table)

cat("\nCox multivariado:\n")
print(cox_multi_table)

cat("\nProportional hazards - multivariado:\n")
print(zph_multi)





#-----------------------------------------------------------#
# Forest plot Cox (Enhanced) + PH Table
#-----------------------------------------------------------#

library(ggplot2)
library(dplyr)
library(broom)
library(patchwork)

# 1. Limpieza de datos y preparación de etiquetas
cox_forest_df <- broom::tidy(
  cox_multi,
  exponentiate = TRUE,
  conf.int = TRUE
) %>%
  mutate(
    term_clean = case_when(
      term == "cluster_k32" ~ "Cluster 2 vs Cluster 1",
      term == "cluster_k33" ~ "Cluster 3 vs Cluster 1",
      term == "age" ~ "Age at diagnosis",
      term == "year" ~ "Year of diagnosis",
      TRUE ~ term
    ),
    # Clasificación para colores según impacto y significancia
    status = case_when(
      p.value < 0.05 & estimate > 1 ~ "Increased Risk",
      p.value < 0.05 & estimate < 1 ~ "Protective",
      TRUE ~ "Not Significant"
    ),
    HR_CI = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high),
    p_value_clean = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
  )

# 2. Guardar tabla Cox bonita (CSV)
write.csv(
  cox_forest_df[, c("term_clean", "estimate", "conf.low", "conf.high", "HR_CI", "p_value_clean")],
  file = file.path(outdir, "cox_multivariate_forest_table_clean.csv"),
  row.names = FALSE
)

# 3. Creación del Forest Plot (Panel Izquierdo)
p_forest <- ggplot(cox_forest_df, aes(x = estimate, y = reorder(term_clean, estimate))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "darkgrey", size = 0.8) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = status), height = 0.2, size = 0.8) +
  geom_point(aes(color = status), size = 3.5) +
  scale_x_log10(breaks = c(0.1, 0.5, 1, 2, 5, 10)) +
  scale_color_manual(values = c(
    "Increased Risk" = "#D55E00", # Naranja
    "Protective" = "#0072B2",     # Azul
    "Not Significant" = "black"
  )) +
  labs(
    title = "Multivariate Cox Regression: Overall Survival",
    x = "Hazard Ratio (95% CI, log scale)",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "bold", size = 11),
    plot.margin = margin(5, 10, 5, 5) 
  )

# 4. Creación de la Tabla de Datos (Panel Derecho) con espaciado corregido
p_table <- ggplot(cox_forest_df, aes(y = reorder(term_clean, estimate))) +
  geom_text(aes(x = 0, label = HR_CI), size = 4, hjust = 0.5) +
  geom_text(aes(x = 1.2, label = p_value_clean), size = 4, hjust = 0.5) +
  annotate("text", x = 0, y = Inf, label = "HR (95% CI)", fontface = "bold", vjust = 2, size = 4.5) +
  annotate("text", x = 1.2, y = Inf, label = "p-value", fontface = "bold", vjust = 2, size = 4.5) +
  scale_x_continuous(limits = c(-0.8, 2)) + 
  theme_void() +
  theme(
    plot.margin = margin(l = 20, r = 20)
  )

# 5. Ensamble Final con patchwork (Proporción 2:1 para más espacio a la derecha)
combined_forest <- p_forest + p_table + plot_layout(widths = c(2, 1))

# 6. Guardado del gráfico en alta resolución
ggsave(
  filename = file.path(outdir, "cox_multivariate_forest_final.png"),
  plot = combined_forest,
  width = 12, 
  height = 5,
  dpi = 300
)

ggsave(
  filename = file.path(outdir, "cox_multivariate_forest_final.pdf"),
  plot = combined_forest,
  width = 12,
  height = 5
)

# 7. Tabla de proporcionalidad de riesgos (zph)
ph_table <- as.data.frame(zph_multi$table) %>%
  tibble::rownames_to_column("variable") %>%
  mutate(
    p_interpretation = case_when(
      p < 0.05 ~ "Possible PH violation",
      p >= 0.05 ~ "No evidence of PH violation"
    )
  )

write.csv(
  ph_table,
  file = file.path(outdir, "cox_zph_multivariate_clean.csv"),
  row.names = FALSE
)

# Outputs finales en consola
print(combined_forest)
print(ph_table)





#BiocManager::install("survcomp")
library(survcomp)

cindex <- concordance.index(
  x = predict(cox_multi),
  surv.time = meta_surv$OS_days,
  surv.event = meta_surv$OS_event
)

cindex$c.index



#---------------------------#
# COMPARACIÓN DE MODELOS (C-INDEX)
#---------------------------#

cat("\n=== COMPARACIÓN DE MODELOS (C-INDEX) ===\n")

# Modelo BASE (sin clusters)
cox_base <- coxph(
  Surv(OS_days, OS_event) ~ age + year,
  data = meta_surv
)

# C-index modelo base
cindex_base <- summary(cox_base)$concordance[1]

# C-index modelo completo (ya lo tienes)
cindex_full <- summary(cox_multi)$concordance[1]

cat("\nC-index modelo BASE (age + year):", round(cindex_base, 3), "\n")
cat("C-index modelo COMPLETO (+ cluster):", round(cindex_full, 3), "\n")

cat("\nDiferencia:", round(cindex_full - cindex_base, 3), "\n")

# Guardar en archivo
cindex_df <- data.frame(
  Model = c("Base (age + year)", "Full (+ cluster)"),
  C_index = c(cindex_base, cindex_full)
)

write.csv(
  cindex_df,
  file = file.path(outdir, "cindex_model_comparison.csv"),
  row.names = FALSE
)