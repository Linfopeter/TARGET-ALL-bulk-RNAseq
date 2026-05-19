### 03_filter_U18_TARGET_ALL_P2.R
### GOAL: Subset TARGET-ALL-P2 RNA-seq to patients < 18 years at diagnosis
### INPUT:
###   - TARGET_ALL_P2_STAR_counts.RData  (object: counts)
###   - TARGET_ALL_P2_metadata_FULL_per_file.rds (object: final_meta)
### OUTPUT (saved inside folder subsets/U18/):
###   - TARGET_ALL_P2_U18_counts_and_metadata.RData
###   - TARGET_ALL_P2_U18_metadata_per_file.rds

suppressPackageStartupMessages({
  library(dplyr)
})

#---------------------------#
# 0) Output folder          #
#---------------------------#
outdir <- file.path("subsets", "U18")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

#---------------------------#
# 1) Load inputs            #
#---------------------------#
load("TARGET_ALL_P2_STAR_counts.RData")  # loads 'counts'
final_meta <- readRDS("TARGET_ALL_P2_metadata_FULL_per_file.rds")

#---------------------------#
# 2) Sanity: alignment      #
#---------------------------#
stopifnot(all(colnames(counts) == final_meta$counts_colname))

#---------------------------#
# 3) Age handling           #
#---------------------------#
# GDC 'age_at_diagnosis' is typically in DAYS.
# Convert to years; keep original column.
final_meta <- final_meta %>%
  mutate(
    age_at_diagnosis_num = suppressWarnings(as.numeric(age_at_diagnosis)),
    age_at_diagnosis_years = age_at_diagnosis_num / 365.25
  )

# Report missing ages
cat("N muestras:", nrow(final_meta), "\n")
cat("N con age_at_diagnosis NA:", sum(is.na(final_meta$age_at_diagnosis_years)), "\n")

#---------------------------#
# 4) Filter <18 years       #
#---------------------------#
keep <- !is.na(final_meta$age_at_diagnosis_years) & final_meta$age_at_diagnosis_years < 18

final_meta_u18 <- final_meta[keep, , drop = FALSE]
counts_u18 <- counts[, final_meta_u18$counts_colname, drop = FALSE]

#---------------------------#
# 5) Post-filter checks     #
#---------------------------#
stopifnot(ncol(counts_u18) == nrow(final_meta_u18))
stopifnot(all(colnames(counts_u18) == final_meta_u18$counts_colname))

cat("=== FILTRO U18 ===\n")
cat("Muestras originales:", ncol(counts), "\n")
cat("Muestras U18:", ncol(counts_u18), "\n")
cat("Edad (años) U18 - resumen:\n")
print(summary(final_meta_u18$age_at_diagnosis_years))

cat("\nTipos de muestra (sample_type.x) en U18:\n")
if ("sample_type.x" %in% names(final_meta_u18)) {
  print(table(final_meta_u18$sample_type.x, useNA = "ifany"))
} else {
  cat("No existe sample_type.x en final_meta_u18\n")
}

#---------------------------#
# 6) Save outputs           #
#---------------------------#
save(
  counts_u18,
  final_meta_u18,
  file = file.path(outdir, "TARGET_ALL_P2_U18_counts_and_metadata.RData")
)

saveRDS(
  final_meta_u18,
  file = file.path(outdir, "TARGET_ALL_P2_U18_metadata_per_file.rds")
)

cat("\nGuardado en carpeta:", outdir, "\n")
cat("- TARGET_ALL_P2_U18_counts_and_metadata.RData\n")
cat("- TARGET_ALL_P2_U18_metadata_per_file.rds\n")

#VERIFICAICON DEL FILTRADO 
load("subsets/U18/TARGET_ALL_P2_U18_counts_and_metadata.RData")

cat("=== QC POST-GUARDADO U18 ===\n")
cat("Dim counts_u18:", paste(dim(counts_u18), collapse=" x "), "\n")
cat("Dim final_meta_u18:", paste(dim(final_meta_u18), collapse=" x "), "\n")

cat("Alineación colnames:", all(colnames(counts_u18) == final_meta_u18$counts_colname), "\n")
cat("Edad NA en U18:", sum(is.na(final_meta_u18$age_at_diagnosis_years)), "\n")
cat("Max edad U18:", max(final_meta_u18$age_at_diagnosis_years, na.rm=TRUE), "\n")

cat("\nSample types:\n")
print(table(final_meta_u18$sample_type.x, useNA="ifany"))
