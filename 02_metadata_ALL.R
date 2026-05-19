### 02_metadata_ALL.R
### GOAL: Build the most complete metadata table possible (per RNA-seq file) for TARGET-ALL-P2,
###       aligned 1:1 with the columns of `counts` (STAR - Counts)
### OUTPUT: RDS (recommended). Optional CSV export (commented at end).

library(TCGAbiolinks)
library(dplyr)

#---------------------------#
# 0) Load counts            #
#---------------------------#
load("TARGET_ALL_P2_STAR_counts.RData")  # loads 'counts'

#---------------------------#
# 1) File-level metadata    #
#---------------------------#
query_meta <- GDCquery(
  project = "TARGET-ALL-P2",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

meta_df <- as.data.frame(getResults(query_meta))

# Keep only files present in counts (should be 532)
meta_sub <- meta_df %>%
  filter(file_name %in% colnames(counts))

cat("Archivos en counts:", ncol(counts), "\n")
cat("Archivos matcheados en meta_sub:", nrow(meta_sub), "\n")

#---------------------------#
# 2) Clinical + biospecimen  #
#---------------------------#
clinical <- GDCquery_clinic(project = "TARGET-ALL-P2", type = "clinical")
biospecimen <- GDCquery_clinic(project = "TARGET-ALL-P2", type = "biospecimen")

# Ensure unique column names (clinical/biospecimen can contain duplicates)
names(clinical) <- make.unique(names(clinical))
names(biospecimen) <- make.unique(names(biospecimen))

# Identify submitter_id columns after make.unique()
clinical_submit_col <- if ("submitter_id" %in% names(clinical)) {
  "submitter_id"
} else {
  grep("^submitter_id(\\.|$)", names(clinical), value = TRUE)[1]
}

biospec_submit_col <- if ("submitter_id" %in% names(biospecimen)) {
  "submitter_id"
} else {
  grep("^submitter_id(\\.|$)", names(biospecimen), value = TRUE)[1]
}

stopifnot(!is.na(clinical_submit_col), !is.na(biospec_submit_col))

#---------------------------#
# 3) Merge all metadata      #
#---------------------------#
final_meta <- meta_sub %>%
  mutate(counts_colname = file_name) %>%
  left_join(
    clinical,
    by = setNames(clinical_submit_col, "cases.submitter_id")
  ) %>%
  left_join(
    biospecimen,
    by = setNames(biospec_submit_col, "sample.submitter_id")
  )

# Make final names unique (safe for R objects and downstream handling)
names(final_meta) <- make.unique(names(final_meta))

#---------------------------#
# 4) Reorder to match counts #
#---------------------------#
final_meta <- final_meta %>%
  slice(match(colnames(counts), counts_colname))

cat("Dimensiones final_meta:", paste(dim(final_meta), collapse = " x "), "\n")
cat("NAs en counts_colname:", sum(is.na(final_meta$counts_colname)), "\n")

#---------------------------#
# 5) Save outputs (RDS)      #
#---------------------------#
saveRDS(final_meta, file = "TARGET_ALL_P2_metadata_FULL_per_file.rds")
save(final_meta, clinical, biospecimen, meta_sub,
     file = "TARGET_ALL_P2_metadata_FULL_objects.RData")

cat("Listo: TARGET_ALL_P2_metadata_FULL_per_file.rds\n")

#-------------------------------------------------------------------------------
# OPTIONAL: Export a flattened CSV (ONLY if you need Excel/sharing).
# Note: CSV cannot store list-columns; this converts any list column to text.
#
# is_list_col <- vapply(final_meta, is.list, logical(1))
# list_cols <- names(final_meta)[is_list_col]
# final_meta_csv <- final_meta
#
# if (length(list_cols) > 0) {
#   for (cc in list_cols) {
#     final_meta_csv[[cc]] <- vapply(final_meta_csv[[cc]], function(x) {
#       if (is.null(x) || length(x) == 0) return(NA_character_)
#       paste(unlist(x), collapse = "|")
#     }, character(1))
#   }
# }
#
# write.csv(final_meta_csv, "TARGET_ALL_P2_metadata_FULL_per_file.csv", row.names = FALSE)
# cat("Listo (opcional): TARGET_ALL_P2_metadata_FULL_per_file.csv\n")
#-------------------------------------------------------------------------------



#VERIFICACION DE DATASET (congruencia muestra -- metadata)
## =========================
## QC congruencia counts ↔ metadata (TARGET-ALL-P2)
## Pega y corre todo este bloque
## =========================

suppressPackageStartupMessages(library(dplyr))

cat("\n=== (1) Alineación counts ↔ final_meta ===\n")
cat("colnames(counts) == final_meta$counts_colname ? -> ",
    all(colnames(counts) == final_meta$counts_colname), "\n")

cat("\n=== (2) Completitud de llaves críticas (NA counts) ===\n")
key_na <- c(
  na_case        = sum(is.na(final_meta$cases.submitter_id)),
  na_sample      = sum(is.na(final_meta$sample.submitter_id)),
  na_file_name   = sum(is.na(final_meta$file_name)),
  na_id_uuid     = sum(is.na(final_meta$id)),
  na_counts_name = sum(is.na(final_meta$counts_colname))
)
print(key_na)

cat("\n=== (3) Unicidad / duplicados (deberían ser 0) ===\n")
dup_check <- c(
  duplicated_file_name = sum(duplicated(final_meta$file_name)),
  duplicated_id_uuid   = sum(duplicated(final_meta$id)),
  duplicated_countscol = sum(duplicated(final_meta$counts_colname))
)
print(dup_check)

cat("\n=== (4) sample_type: existencia + comparación x vs y ===\n")
st_cols <- grep("^sample_type(\\.|$)|sample.*type|type.*sample",
                names(final_meta), value = TRUE, ignore.case = TRUE)
cat("Columnas sample_type detectadas:\n")
print(st_cols)

if ("sample_type.x" %in% names(final_meta) || "sample_type.y" %in% names(final_meta)) {
  if ("sample_type.x" %in% names(final_meta)) {
    cat("\nTabla sample_type.x:\n")
    print(table(final_meta$sample_type.x, useNA = "ifany"))
  } else {
    cat("\nNo existe sample_type.x\n")
  }
  
  if ("sample_type.y" %in% names(final_meta)) {
    cat("\nTabla sample_type.y:\n")
    print(table(final_meta$sample_type.y, useNA = "ifany"))
  } else {
    cat("\nNo existe sample_type.y\n")
  }
  
  if ("sample_type.x" %in% names(final_meta) && "sample_type.y" %in% names(final_meta)) {
    cat("\nConcordancia sample_type.x vs sample_type.y (TRUE = iguales):\n")
    print(table(final_meta$sample_type.x == final_meta$sample_type.y, useNA = "ifany"))
  }
} else {
  cat("No se encontró sample_type.x / sample_type.y en final_meta.\n")
}

cat("\n=== (5) Distribución de #muestras por paciente (cases.submitter_id) ===\n")
if ("cases.submitter_id" %in% names(final_meta)) {
  per_case <- table(final_meta$cases.submitter_id)
  cat("Resumen (#muestras por paciente):\n")
  print(summary(as.numeric(per_case)))
  cat("\nTop 20 pacientes con más muestras:\n")
  print(head(sort(per_case, decreasing = TRUE), 20))
  
  cat("\nFrecuencia: cuántos pacientes tienen 1,2,3,... muestras:\n")
  print(head(table(as.numeric(per_case)), 20))
} else {
  cat("No existe cases.submitter_id en final_meta.\n")
}

cat("\n=== (6) Chequeo rápido: file_name -> case y sample (debe ser 1:1) ===\n")
if (all(c("file_name","cases.submitter_id","sample.submitter_id") %in% names(final_meta))) {
  tmp <- final_meta %>%
    group_by(file_name) %>%
    summarise(
      n_cases  = n_distinct(cases.submitter_id),
      n_samples = n_distinct(sample.submitter_id),
      .groups = "drop"
    )
  cat("Archivos con >1 case asociado:\n")
  print(sum(tmp$n_cases > 1))
  cat("Archivos con >1 sample asociado:\n")
  print(sum(tmp$n_samples > 1))
  if (sum(tmp$n_cases > 1) > 0 || sum(tmp$n_samples > 1) > 0) {
    cat("\nEjemplos problemáticos (primeros 10):\n")
    print(head(tmp %>% filter(n_cases > 1 | n_samples > 1), 10))
  }
} else {
  cat("Faltan columnas para este chequeo (file_name/cases.submitter_id/sample.submitter_id).\n")
}

cat("\n=== FIN QC ===\n")


#grep("age|days_to_birth|days_to_diagnosis|year", names(final_meta), value = TRUE, ignore.case = TRUE)
