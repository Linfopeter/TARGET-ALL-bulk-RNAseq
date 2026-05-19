### 04_filter_PRIMARY_U18_TARGET_ALL_P2.R
### GOAL: From the U18 subset, remove recurrent samples and keep ONLY Primary samples.
### INPUT:
###   - subsets/U18/TARGET_ALL_P2_U18_counts_and_metadata.RData
### OUTPUT (saved inside folder subsets/primary_u18/):
###   - TARGET_ALL_P2_PRIMARY_U18_counts_and_metadata.RData
###   - TARGET_ALL_P2_PRIMARY_U18_metadata_per_file.rds

suppressPackageStartupMessages({
  library(dplyr)
})

#---------------------------#
# 0) Output folder          #
#---------------------------#
outdir <- file.path("subsets", "primary_u18")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

#---------------------------#
# 1) Load U18 subset        #
#---------------------------#
load("subsets/U18/TARGET_ALL_P2_U18_counts_and_metadata.RData")  # loads counts_u18, final_meta_u18

#---------------------------#
# 2) Sanity: alignment      #
#---------------------------#
stopifnot(ncol(counts_u18) == nrow(final_meta_u18))
stopifnot(all(colnames(counts_u18) == final_meta_u18$counts_colname))

#---------------------------#
# 3) Filter: keep only PRIMARY
#---------------------------#
# We will use sample_type.x as the main label (it is complete and file-level).
stopifnot("sample_type.x" %in% names(final_meta_u18))

keep_primary <- grepl("^Primary", final_meta_u18$sample_type.x)

final_meta_primary_u18 <- final_meta_u18[keep_primary, , drop = FALSE]
counts_primary_u18 <- counts_u18[, final_meta_primary_u18$counts_colname, drop = FALSE]

#---------------------------#
# 4) Post-filter checks     #
#---------------------------#
stopifnot(ncol(counts_primary_u18) == nrow(final_meta_primary_u18))
stopifnot(all(colnames(counts_primary_u18) == final_meta_primary_u18$counts_colname))

cat("=== FILTRO PRIMARY_U18 ===\n")
cat("Muestras U18 originales:", ncol(counts_u18), "\n")
cat("Muestras PRIMARY_U18:", ncol(counts_primary_u18), "\n")

cat("\nSample types (PRIMARY_U18):\n")
print(table(final_meta_primary_u18$sample_type.x, useNA = "ifany"))

cat("\nEdad (años) PRIMARY_U18 - resumen:\n")
print(summary(final_meta_primary_u18$age_at_diagnosis_years))

#---------------------------#
# 5) Save outputs           #
#---------------------------#
save(
  counts_primary_u18,
  final_meta_primary_u18,
  file = file.path(outdir, "TARGET_ALL_P2_PRIMARY_U18_counts_and_metadata.RData")
)

saveRDS(
  final_meta_primary_u18,
  file = file.path(outdir, "TARGET_ALL_P2_PRIMARY_U18_metadata_per_file.rds")
)

cat("\nGuardado en carpeta:", outdir, "\n")
cat("- TARGET_ALL_P2_PRIMARY_U18_counts_and_metadata.RData\n")
cat("- TARGET_ALL_P2_PRIMARY_U18_metadata_per_file.rds\n")

