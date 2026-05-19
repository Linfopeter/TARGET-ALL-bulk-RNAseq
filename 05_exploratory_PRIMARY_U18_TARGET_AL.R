### 05_exploratory_PRIMARY_U18_TARGET_ALL_P2.R
### GOAL: Exploratory analysis of the final filtered dataset (PRIMARY + U18)
### INPUT:
###   - subsets/primary_u18/TARGET_ALL_P2_PRIMARY_U18_counts_and_metadata.RData
### OUTPUT:
###   - Console summaries (dimensions, variables, structure)

suppressPackageStartupMessages({
  library(dplyr)
})

#---------------------------#
# 1) Load final dataset     #
#---------------------------#
load("subsets/primary_u18/TARGET_ALL_P2_PRIMARY_U18_counts_and_metadata.RData")
# loads:
#   - counts_primary_u18
#   - final_meta_primary_u18

#---------------------------#
# 2) Basic dimensions       #
#---------------------------#
cat("=== DIMENSIONES DEL DATASET ===\n")

cat("\nMatriz de expresión (counts):\n")
cat("- Genes (filas):", nrow(counts_primary_u18), "\n")
cat("- Muestras (columnas):", ncol(counts_primary_u18), "\n")

cat("\nMetadata:\n")
cat("- Filas (muestras):", nrow(final_meta_primary_u18), "\n")
cat("- Columnas (variables):", ncol(final_meta_primary_u18), "\n")

#---------------------------#
# 3) Alignment check        #
#---------------------------#
cat("\n=== CHEQUEO DE ALINEACIÓN ===\n")
cat("¿colnames(counts) == metadata$counts_colname?: ",
    all(colnames(counts_primary_u18) == final_meta_primary_u18$counts_colname), "\n")

#---------------------------#
# 4) Variables in metadata  #
#---------------------------#
cat("\n=== VARIABLES EN METADATA ===\n")
vars <- colnames(final_meta_primary_u18)

cat("Número total de variables:", length(vars), "\n\n")
print(vars)

#---------------------------#
# 5) Structure overview     #
#---------------------------#
cat("\n=== ESTRUCTURA GENERAL ===\n")

cat("\nEstructura de counts:\n")
str(counts_primary_u18)

cat("\nEstructura de metadata:\n")
str(final_meta_primary_u18)

#---------------------------#
# 6) Quick summary stats    #
#---------------------------#
cat("\n=== RESUMEN RÁPIDO ===\n")

cat("\nDistribución de edad (años):\n")
if ("age_at_diagnosis_years" %in% names(final_meta_primary_u18)) {
  print(summary(final_meta_primary_u18$age_at_diagnosis_years))
} else {
  cat("No existe age_at_diagnosis_years\n")
}

cat("\nTipos de muestra:\n")
if ("sample_type.x" %in% names(final_meta_primary_u18)) {
  print(table(final_meta_primary_u18$sample_type.x, useNA = "ifany"))
} else {
  cat("No existe sample_type.x\n")
}

#---------------------------#
# 7) Final message          #
#---------------------------#
cat("\n=== EXPLORACIÓN COMPLETADA ===\n")

#---------------------------#
# 8) Variantes diagnósticas #
#---------------------------#
cat("\n=== VARIANTES DIAGNÓSTICAS ===\n")

clean_values <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "N/A", "Not Reported", "not reported", "unknown", "Unknown")] <- NA
  x
}

diag_cols <- c("primary_diagnosis", "icd_10_code", "morphology")

missing_cols <- setdiff(diag_cols, colnames(final_meta_primary_u18))
if (length(missing_cols) > 0) {
  cat("Faltan estas columnas en metadata:\n")
  print(missing_cols)
} else {
  
  final_meta_primary_u18 <- final_meta_primary_u18 %>%
    mutate(
      primary_diagnosis_clean = clean_values(primary_diagnosis),
      icd_10_code_clean       = clean_values(icd_10_code),
      morphology_clean        = clean_values(morphology)
    )
  
  #---------------------------#
  # A) Número de variantes    #
  #---------------------------#
  cat("\n--- Número de variantes por columna ---\n")
  cat("primary_diagnosis:", dplyr::n_distinct(final_meta_primary_u18$primary_diagnosis_clean, na.rm = TRUE), "\n")
  cat("icd_10_code:",       dplyr::n_distinct(final_meta_primary_u18$icd_10_code_clean, na.rm = TRUE), "\n")
  cat("morphology:",        dplyr::n_distinct(final_meta_primary_u18$morphology_clean, na.rm = TRUE), "\n")
  
  cat("\nNAs por columna:\n")
  cat("primary_diagnosis:", sum(is.na(final_meta_primary_u18$primary_diagnosis_clean)), "\n")
  cat("icd_10_code:",       sum(is.na(final_meta_primary_u18$icd_10_code_clean)), "\n")
  cat("morphology:",        sum(is.na(final_meta_primary_u18$morphology_clean)), "\n")
  
  #---------------------------#
  # B) Frecuencias por campo  #
  #---------------------------#
  cat("\n--- Variantes en primary_diagnosis ---\n")
  print(sort(table(final_meta_primary_u18$primary_diagnosis_clean, useNA = "ifany"), decreasing = TRUE))
  
  cat("\n--- Variantes en icd_10_code ---\n")
  print(sort(table(final_meta_primary_u18$icd_10_code_clean, useNA = "ifany"), decreasing = TRUE))
  
  cat("\n--- Variantes en morphology ---\n")
  print(sort(table(final_meta_primary_u18$morphology_clean, useNA = "ifany"), decreasing = TRUE))
  
  #---------------------------#
  # C) Combinaciones únicas   #
  #---------------------------#
  cat("\n--- Combinaciones únicas entre primary_diagnosis + icd_10_code + morphology ---\n")
  combo_tab <- final_meta_primary_u18 %>%
    count(
      primary_diagnosis_clean,
      icd_10_code_clean,
      morphology_clean,
      sort = TRUE
    )
  
  print(combo_tab)
}

#---------------------------#
# 9) Año de diagnóstico     #
#---------------------------#
cat("\n=== AÑO DE DIAGNÓSTICO ===\n")

if ("year_of_diagnosis" %in% names(final_meta_primary_u18)) {
  
  year_diag <- final_meta_primary_u18$year_of_diagnosis
  
  # limpiar
  year_diag_clean <- as.numeric(trimws(as.character(year_diag)))
  
  cat("\nResumen numérico:\n")
  print(summary(year_diag_clean))
  
  cat("\nNAs:", sum(is.na(year_diag_clean)), "\n")
  
  cat("\nFrecuencia por año:\n")
  print(sort(table(year_diag_clean, useNA = "ifany")))
  
} else {
  cat("No existe year_of_diagnosis\n")
}


#---------------------------#
# 10) Follow-up             #
#---------------------------#
cat("\n=== DAYS TO LAST FOLLOW-UP ===\n")

if ("days_to_last_follow_up" %in% names(final_meta_primary_u18)) {
  
  follow_up <- final_meta_primary_u18$days_to_last_follow_up
  
  # limpiar
  follow_up_clean <- as.numeric(trimws(as.character(follow_up)))
  
  cat("\nResumen numérico (días):\n")
  print(summary(follow_up_clean))
  
  cat("\nNAs:", sum(is.na(follow_up_clean)), "\n")
  
  cat("\nFollow-up en años (aprox):\n")
  follow_up_years <- follow_up_clean / 365
  print(summary(follow_up_years))
  
  # categorías útiles
  cat("\nDistribución por categorías de follow-up:\n")
  follow_up_cat <- cut(
    follow_up_years,
    breaks = c(0, 1, 3, 5, 10, Inf),
    labels = c("<1 año", "1-3 años", "3-5 años", "5-10 años", ">10 años"),
    include.lowest = TRUE
  )
  
  print(table(follow_up_cat, useNA = "ifany"))
  
} else {
  cat("No existe days_to_last_follow_up\n")
}

#---------------------------#
# 11) Investigación de NAs  #
#---------------------------#
cat("\n=== INVESTIGACIÓN DE NAs EN FOLLOW-UP ===\n")

if ("days_to_last_follow_up" %in% names(final_meta_primary_u18)) {
  
  raw_follow_up <- final_meta_primary_u18$days_to_last_follow_up
  
  # ver valores originales problemáticos
  problematic_idx <- which(is.na(as.numeric(trimws(as.character(raw_follow_up)))))
  
  cat("\nÍndices con problema:\n")
  print(problematic_idx)
  
  cat("\nValores originales (crudos):\n")
  print(raw_follow_up[problematic_idx])
  
  cat("\nConteo de valores únicos problemáticos:\n")
  print(table(raw_follow_up[problematic_idx], useNA = "ifany"))
  
  # ver filas completas para contexto
  cat("\nFilas completas en metadata (para inspección):\n")
  print(final_meta_primary_u18[problematic_idx, 
                              c("cases", "days_to_last_follow_up", "vital_status", "year_of_diagnosis")])
  
} else {
  cat("No existe days_to_last_follow_up\n")
}

#---------------------------#
# 12) Filtrado final        #
#---------------------------#
cat("\n=== FILTRADO FINAL DE MUESTRAS PROBLEMÁTICAS ===\n")

# Copias de trabajo
counts_final <- counts_primary_u18
meta_final   <- final_meta_primary_u18

# Limpiar days_to_last_follow_up usando la misma lógica que en diagnósticos
follow_up_raw <- clean_values(meta_final$days_to_last_follow_up)
follow_up_num <- suppressWarnings(as.numeric(follow_up_raw))

# Detectar muestras problemáticas:
# aquí se eliminan las que no tienen follow-up numérico válido
problematic_idx <- which(is.na(follow_up_num))

cat("\nNúmero de muestras problemáticas detectadas:", length(problematic_idx), "\n")

if (length(problematic_idx) > 0) {
  cat("\nMuestras problemáticas:\n")
  print(meta_final[problematic_idx,
                   c("cases", "counts_colname", "days_to_last_follow_up", "vital_status", "year_of_diagnosis")])
  
  # Filtrar metadata
  meta_final <- meta_final[-problematic_idx, , drop = FALSE]
  
  # Filtrar counts respetando alineación por columnas
  counts_final <- counts_final[, meta_final$counts_colname, drop = FALSE]
}

cat("\nDimensiones después del filtrado:\n")
cat("- Genes (filas):", nrow(counts_final), "\n")
cat("- Muestras (columnas):", ncol(counts_final), "\n")
cat("- Filas en metadata:", nrow(meta_final), "\n")

cat("\nChequeo de alineación post-filtrado:\n")
cat("¿colnames(counts_final) == meta_final$counts_colname?: ",
    all(colnames(counts_final) == meta_final$counts_colname), "\n")


#---------------------------#
# 13) Re-exploración final  #
#---------------------------#
cat("\n=== RE-EXPLORACIÓN FINAL TRAS FILTRADO ===\n")

# A) Variables en metadata
cat("\n=== VARIABLES EN METADATA FINAL ===\n")
vars_final <- colnames(meta_final)
cat("Número total de variables:", length(vars_final), "\n\n")
print(vars_final)

# B) Estructura general
cat("\n=== ESTRUCTURA GENERAL FINAL ===\n")

cat("\nEstructura de counts_final:\n")
str(counts_final)

cat("\nEstructura de meta_final:\n")
str(meta_final)

# C) Edad
# C) Edad
cat("\n=== RESUMEN DE EDAD FINAL ===\n")
if ("age_at_diagnosis_years" %in% names(meta_final)) {
  print(summary(meta_final$age_at_diagnosis_years))
  cat("NAs en age_at_diagnosis_years:", sum(is.na(meta_final$age_at_diagnosis_years)), "\n")
  
  # 🔥 NUEVO BLOQUE
  cat("\nFrecuencia exacta por edad (años):\n")
  
  age_clean <- meta_final$age_at_diagnosis_years
  
  # si quieres redondear (recomendado porque vienen con decimales)
  age_rounded <- floor(age_clean)
  
  print(sort(table(age_rounded, useNA = "ifany")))
  
  # opcional: tabla ordenada tipo data.frame
  cat("\nTabla ordenada de edades:\n")
  age_table <- as.data.frame(table(age_rounded))
  colnames(age_table) <- c("Edad", "Frecuencia")
  age_table <- age_table[order(as.numeric(as.character(age_table$Edad))), ]
  
  print(age_table)
  
} else {
  cat("No existe age_at_diagnosis_years\n")
}

# D) Tipos de muestra
cat("\n=== TIPOS DE MUESTRA FINAL ===\n")
if ("sample_type.x" %in% names(meta_final)) {
  print(table(meta_final$sample_type.x, useNA = "ifany"))
} else {
  cat("No existe sample_type.x\n")
}

# E) Variantes diagnósticas
cat("\n=== VARIANTES DIAGNÓSTICAS FINAL ===\n")

diag_cols <- c("primary_diagnosis", "icd_10_code", "morphology")
missing_cols_final <- setdiff(diag_cols, colnames(meta_final))

if (length(missing_cols_final) > 0) {
  cat("Faltan estas columnas en metadata final:\n")
  print(missing_cols_final)
} else {
  
  meta_final <- meta_final %>%
    mutate(
      primary_diagnosis_clean = clean_values(primary_diagnosis),
      icd_10_code_clean       = clean_values(icd_10_code),
      morphology_clean        = clean_values(morphology)
    )
  
  cat("\n--- Número de variantes por columna ---\n")
  cat("primary_diagnosis:", dplyr::n_distinct(meta_final$primary_diagnosis_clean, na.rm = TRUE), "\n")
  cat("icd_10_code:",       dplyr::n_distinct(meta_final$icd_10_code_clean, na.rm = TRUE), "\n")
  cat("morphology:",        dplyr::n_distinct(meta_final$morphology_clean, na.rm = TRUE), "\n")
  
  cat("\nNAs por columna:\n")
  cat("primary_diagnosis:", sum(is.na(meta_final$primary_diagnosis_clean)), "\n")
  cat("icd_10_code:",       sum(is.na(meta_final$icd_10_code_clean)), "\n")
  cat("morphology:",        sum(is.na(meta_final$morphology_clean)), "\n")
  
  cat("\n--- Variantes en primary_diagnosis ---\n")
  print(sort(table(meta_final$primary_diagnosis_clean, useNA = "ifany"), decreasing = TRUE))
  
  cat("\n--- Variantes en icd_10_code ---\n")
  print(sort(table(meta_final$icd_10_code_clean, useNA = "ifany"), decreasing = TRUE))
  
  cat("\n--- Variantes en morphology ---\n")
  print(sort(table(meta_final$morphology_clean, useNA = "ifany"), decreasing = TRUE))
  
  cat("\n--- Combinaciones únicas entre primary_diagnosis + icd_10_code + morphology ---\n")
  combo_tab_final <- meta_final %>%
    count(
      primary_diagnosis_clean,
      icd_10_code_clean,
      morphology_clean,
      sort = TRUE
    )
  
  print(combo_tab_final)
}

# F) Año de diagnóstico
cat("\n=== AÑO DE DIAGNÓSTICO FINAL ===\n")

if ("year_of_diagnosis" %in% names(meta_final)) {
  year_diag_final <- as.numeric(clean_values(meta_final$year_of_diagnosis))
  
  cat("\nResumen numérico:\n")
  print(summary(year_diag_final))
  
  cat("\nNAs:", sum(is.na(year_diag_final)), "\n")
  
  cat("\nFrecuencia por año:\n")
  print(sort(table(year_diag_final, useNA = "ifany")))
  
} else {
  cat("No existe year_of_diagnosis\n")
}

# G) Follow-up
cat("\n=== DAYS TO LAST FOLLOW-UP FINAL ===\n")

if ("days_to_last_follow_up" %in% names(meta_final)) {
  
  follow_up_final <- as.numeric(clean_values(meta_final$days_to_last_follow_up))
  
  cat("\nResumen numérico (días):\n")
  print(summary(follow_up_final))
  
  cat("\nNAs:", sum(is.na(follow_up_final)), "\n")
  
  cat("\nFollow-up en años (aprox):\n")
  follow_up_final_years <- follow_up_final / 365
  print(summary(follow_up_final_years))
  
  cat("\nDistribución por categorías de follow-up:\n")
  follow_up_final_cat <- cut(
    follow_up_final_years,
    breaks = c(0, 1, 3, 5, 10, Inf),
    labels = c("<1 año", "1-3 años", "3-5 años", "5-10 años", ">10 años"),
    include.lowest = TRUE
  )
  
  print(table(follow_up_final_cat, useNA = "ifany"))
  
} else {
  cat("No existe days_to_last_follow_up\n")
}


#---------------------------#
# 14) Guardado final        #
#---------------------------#
cat("\n=== GUARDANDO DATASET FINAL ===\n")

dir.create("subsets/primary_u18/final_analysis", recursive = TRUE, showWarnings = FALSE)

# Guardado completo en RData (SIN problemas)
save(
  counts_final,
  meta_final,
  file = "subsets/primary_u18/final_analysis/TARGET_ALL_P2_PRIMARY_U18_FINAL_FOR_ANALYSIS.RData"
)

# 🔥 ELIMINAR columnas tipo list para CSV
meta_final_csv <- meta_final[, !sapply(meta_final, is.list)]

cat("\nColumnas eliminadas para CSV (tipo list):\n")
print(names(meta_final)[sapply(meta_final, is.list)])

# Guardar CSV limpio
write.csv(
  meta_final_csv,
  file = "subsets/primary_u18/final_analysis/TARGET_ALL_P2_PRIMARY_U18_FINAL_METADATA.csv",
  row.names = FALSE
)

cat("\nArchivos guardados:\n")
cat("- RData (completo, incluye listas)\n")
cat("- CSV (sin columnas tipo list)\n")

cat("\n=== PROCESO FINAL COMPLETADO ===\n")


#---------------------------#
# 15) Reporte TXT final     #
#---------------------------#
cat("\n=== GENERANDO REPORTE TXT FINAL ===\n")

report_file <- "subsets/primary_u18/final_analysis/TARGET_ALL_P2_PRIMARY_U18_FULL_REPORT.txt"

# Objetos "antes del filtrado"
counts_before <- counts_primary_u18
meta_before   <- final_meta_primary_u18

# Objetos "después del filtrado"
# counts_final y meta_final ya existen en tu script

# Recalcular índices problemáticos en el dataset original
follow_up_before_raw <- clean_values(meta_before$days_to_last_follow_up)
follow_up_before_num <- suppressWarnings(as.numeric(follow_up_before_raw))
problematic_idx_before <- which(is.na(follow_up_before_num))

# Utilidad para imprimir secciones
write_section <- function(title) {
  cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
  cat(title, "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n", sep = "")
}

sink(report_file)

cat("REPORTE COMPLETO DEL DATASET TARGET ALL P2 PRIMARY U18\n")
cat("Generado el:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

#------------------------------------------------#
# 1) RESUMEN GENERAL ANTES DEL FILTRADO          #
#------------------------------------------------#
write_section("1) RESUMEN GENERAL ANTES DEL FILTRADO")

cat("\nDimensiones:\n")
cat("- Genes (counts):", nrow(counts_before), "\n")
cat("- Muestras (counts):", ncol(counts_before), "\n")
cat("- Filas metadata:", nrow(meta_before), "\n")
cat("- Columnas metadata:", ncol(meta_before), "\n")

cat("\nChequeo de alineación:\n")
cat("- colnames(counts_before) == meta_before$counts_colname: ",
    all(colnames(counts_before) == meta_before$counts_colname), "\n", sep = "")

cat("\nVariables en metadata:\n")
print(colnames(meta_before))

cat("\nEstructura de counts_before:\n")
str(counts_before)

cat("\nEstructura de meta_before:\n")
str(meta_before)

#------------------------------------------------#
# 2) EDAD ANTES DEL FILTRADO                     #
#------------------------------------------------#
write_section("2) EDAD ANTES DEL FILTRADO")

if ("age_at_diagnosis_years" %in% names(meta_before)) {
  age_before <- meta_before$age_at_diagnosis_years
  
  cat("\nResumen numérico de edad:\n")
  print(summary(age_before))
  
  cat("\nNAs en edad:", sum(is.na(age_before)), "\n")
  
  cat("\nFrecuencia exacta por edad (floor):\n")
  age_before_floor <- floor(age_before)
  print(sort(table(age_before_floor, useNA = "ifany")))
  
} else {
  cat("No existe age_at_diagnosis_years\n")
}

#------------------------------------------------#
# 3) TIPOS DE MUESTRA ANTES DEL FILTRADO         #
#------------------------------------------------#
write_section("3) TIPOS DE MUESTRA ANTES DEL FILTRADO")

if ("sample_type.x" %in% names(meta_before)) {
  print(table(meta_before$sample_type.x, useNA = "ifany"))
} else {
  cat("No existe sample_type.x\n")
}

#------------------------------------------------#
# 4) DIAGNÓSTICO ANTES DEL FILTRADO              #
#------------------------------------------------#
write_section("4) DIAGNÓSTICO ANTES DEL FILTRADO")

diag_cols <- c("primary_diagnosis", "icd_10_code", "morphology")
missing_before <- setdiff(diag_cols, colnames(meta_before))

if (length(missing_before) > 0) {
  cat("Faltan estas columnas:\n")
  print(missing_before)
} else {
  meta_before_report <- meta_before %>%
    mutate(
      primary_diagnosis_clean_report = clean_values(primary_diagnosis),
      icd_10_code_clean_report       = clean_values(icd_10_code),
      morphology_clean_report        = clean_values(morphology)
    )
  
  cat("\nNúmero de variantes por columna:\n")
  cat("- primary_diagnosis:",
      dplyr::n_distinct(meta_before_report$primary_diagnosis_clean_report, na.rm = TRUE), "\n")
  cat("- icd_10_code:",
      dplyr::n_distinct(meta_before_report$icd_10_code_clean_report, na.rm = TRUE), "\n")
  cat("- morphology:",
      dplyr::n_distinct(meta_before_report$morphology_clean_report, na.rm = TRUE), "\n")
  
  cat("\nNAs por columna:\n")
  cat("- primary_diagnosis:", sum(is.na(meta_before_report$primary_diagnosis_clean_report)), "\n")
  cat("- icd_10_code:", sum(is.na(meta_before_report$icd_10_code_clean_report)), "\n")
  cat("- morphology:", sum(is.na(meta_before_report$morphology_clean_report)), "\n")
  
  cat("\nFrecuencia de primary_diagnosis:\n")
  print(sort(table(meta_before_report$primary_diagnosis_clean_report, useNA = "ifany"), decreasing = TRUE))
  
  cat("\nFrecuencia de icd_10_code:\n")
  print(sort(table(meta_before_report$icd_10_code_clean_report, useNA = "ifany"), decreasing = TRUE))
  
  cat("\nFrecuencia de morphology:\n")
  print(sort(table(meta_before_report$morphology_clean_report, useNA = "ifany"), decreasing = TRUE))
  
  cat("\nCombinaciones únicas:\n")
  combo_before <- meta_before_report %>%
    count(
      primary_diagnosis_clean_report,
      icd_10_code_clean_report,
      morphology_clean_report,
      sort = TRUE
    )
  print(combo_before)
}

#------------------------------------------------#
# 5) AÑO DE DIAGNÓSTICO ANTES DEL FILTRADO       #
#------------------------------------------------#
write_section("5) AÑO DE DIAGNÓSTICO ANTES DEL FILTRADO")

if ("year_of_diagnosis" %in% names(meta_before)) {
  year_before <- suppressWarnings(as.numeric(clean_values(meta_before$year_of_diagnosis)))
  
  cat("\nResumen numérico:\n")
  print(summary(year_before))
  
  cat("\nNAs:", sum(is.na(year_before)), "\n")
  
  cat("\nFrecuencia por año:\n")
  print(sort(table(year_before, useNA = "ifany")))
  
} else {
  cat("No existe year_of_diagnosis\n")
}

#------------------------------------------------#
# 6) FOLLOW-UP ANTES DEL FILTRADO                #
#------------------------------------------------#
write_section("6) FOLLOW-UP ANTES DEL FILTRADO")

if ("days_to_last_follow_up" %in% names(meta_before)) {
  cat("\nResumen numérico (días):\n")
  print(summary(follow_up_before_num))
  
  cat("\nNAs:", sum(is.na(follow_up_before_num)), "\n")
  
  cat("\nFollow-up en años:\n")
  follow_up_before_years <- follow_up_before_num / 365
  print(summary(follow_up_before_years))
  
  cat("\nDistribución por categorías:\n")
  follow_up_before_cat <- cut(
    follow_up_before_years,
    breaks = c(0, 1, 3, 5, 10, Inf),
    labels = c("<1 año", "1-3 años", "3-5 años", "5-10 años", ">10 años"),
    include.lowest = TRUE
  )
  print(table(follow_up_before_cat, useNA = "ifany"))
  
} else {
  cat("No existe days_to_last_follow_up\n")
}

#------------------------------------------------#
# 7) MUESTRAS PROBLEMÁTICAS FILTRADAS            #
#------------------------------------------------#
write_section("7) MUESTRAS PROBLEMÁTICAS FILTRADAS")

cat("\nNúmero de muestras problemáticas detectadas:", length(problematic_idx_before), "\n")

if (length(problematic_idx_before) > 0) {
  print(meta_before[problematic_idx_before,
                    c("cases", "counts_colname", "days_to_last_follow_up", "vital_status", "year_of_diagnosis")])
} else {
  cat("No se detectaron muestras problemáticas.\n")
}

#------------------------------------------------#
# 8) RESUMEN GENERAL DESPUÉS DEL FILTRADO        #
#------------------------------------------------#
write_section("8) RESUMEN GENERAL DESPUÉS DEL FILTRADO")

cat("\nDimensiones:\n")
cat("- Genes (counts_final):", nrow(counts_final), "\n")
cat("- Muestras (counts_final):", ncol(counts_final), "\n")
cat("- Filas meta_final:", nrow(meta_final), "\n")
cat("- Columnas meta_final:", ncol(meta_final), "\n")

cat("\nChequeo de alineación:\n")
cat("- colnames(counts_final) == meta_final$counts_colname: ",
    all(colnames(counts_final) == meta_final$counts_colname), "\n", sep = "")

cat("\nVariables en metadata final:\n")
print(colnames(meta_final))

cat("\nEstructura de counts_final:\n")
str(counts_final)

cat("\nEstructura de meta_final:\n")
str(meta_final)

#------------------------------------------------#
# 9) EDAD DESPUÉS DEL FILTRADO                   #
#------------------------------------------------#
write_section("9) EDAD DESPUÉS DEL FILTRADO")

if ("age_at_diagnosis_years" %in% names(meta_final)) {
  age_final <- meta_final$age_at_diagnosis_years
  
  cat("\nResumen numérico de edad:\n")
  print(summary(age_final))
  
  cat("\nNAs en edad:", sum(is.na(age_final)), "\n")
  
  cat("\nFrecuencia exacta por edad (floor):\n")
  age_final_floor <- floor(age_final)
  print(sort(table(age_final_floor, useNA = "ifany")))
  
} else {
  cat("No existe age_at_diagnosis_years\n")
}

#------------------------------------------------#
# 10) TIPOS DE MUESTRA DESPUÉS DEL FILTRADO      #
#------------------------------------------------#
write_section("10) TIPOS DE MUESTRA DESPUÉS DEL FILTRADO")

if ("sample_type.x" %in% names(meta_final)) {
  print(table(meta_final$sample_type.x, useNA = "ifany"))
} else {
  cat("No existe sample_type.x\n")
}

#------------------------------------------------#
# 11) DIAGNÓSTICO DESPUÉS DEL FILTRADO           #
#------------------------------------------------#
write_section("11) DIAGNÓSTICO DESPUÉS DEL FILTRADO")

missing_after <- setdiff(diag_cols, colnames(meta_final))

if (length(missing_after) > 0) {
  cat("Faltan estas columnas:\n")
  print(missing_after)
} else {
  meta_final_report <- meta_final %>%
    mutate(
      primary_diagnosis_clean_report = clean_values(primary_diagnosis),
      icd_10_code_clean_report       = clean_values(icd_10_code),
      morphology_clean_report        = clean_values(morphology)
    )
  
  cat("\nNúmero de variantes por columna:\n")
  cat("- primary_diagnosis:",
      dplyr::n_distinct(meta_final_report$primary_diagnosis_clean_report, na.rm = TRUE), "\n")
  cat("- icd_10_code:",
      dplyr::n_distinct(meta_final_report$icd_10_code_clean_report, na.rm = TRUE), "\n")
  cat("- morphology:",
      dplyr::n_distinct(meta_final_report$morphology_clean_report, na.rm = TRUE), "\n")
  
  cat("\nNAs por columna:\n")
  cat("- primary_diagnosis:", sum(is.na(meta_final_report$primary_diagnosis_clean_report)), "\n")
  cat("- icd_10_code:", sum(is.na(meta_final_report$icd_10_code_clean_report)), "\n")
  cat("- morphology:", sum(is.na(meta_final_report$morphology_clean_report)), "\n")
  
  cat("\nFrecuencia de primary_diagnosis:\n")
  print(sort(table(meta_final_report$primary_diagnosis_clean_report, useNA = "ifany"), decreasing = TRUE))
  
  cat("\nFrecuencia de icd_10_code:\n")
  print(sort(table(meta_final_report$icd_10_code_clean_report, useNA = "ifany"), decreasing = TRUE))
  
  cat("\nFrecuencia de morphology:\n")
  print(sort(table(meta_final_report$morphology_clean_report, useNA = "ifany"), decreasing = TRUE))
  
  cat("\nCombinaciones únicas:\n")
  combo_after <- meta_final_report %>%
    count(
      primary_diagnosis_clean_report,
      icd_10_code_clean_report,
      morphology_clean_report,
      sort = TRUE
    )
  print(combo_after)
}

#------------------------------------------------#
# 12) AÑO DE DIAGNÓSTICO DESPUÉS DEL FILTRADO    #
#------------------------------------------------#
write_section("12) AÑO DE DIAGNÓSTICO DESPUÉS DEL FILTRADO")

if ("year_of_diagnosis" %in% names(meta_final)) {
  year_after <- suppressWarnings(as.numeric(clean_values(meta_final$year_of_diagnosis)))
  
  cat("\nResumen numérico:\n")
  print(summary(year_after))
  
  cat("\nNAs:", sum(is.na(year_after)), "\n")
  
  cat("\nFrecuencia por año:\n")
  print(sort(table(year_after, useNA = "ifany")))
  
} else {
  cat("No existe year_of_diagnosis\n")
}

#------------------------------------------------#
# 13) FOLLOW-UP DESPUÉS DEL FILTRADO             #
#------------------------------------------------#
write_section("13) FOLLOW-UP DESPUÉS DEL FILTRADO")

if ("days_to_last_follow_up" %in% names(meta_final)) {
  follow_up_after <- suppressWarnings(as.numeric(clean_values(meta_final$days_to_last_follow_up)))
  
  cat("\nResumen numérico (días):\n")
  print(summary(follow_up_after))
  
  cat("\nNAs:", sum(is.na(follow_up_after)), "\n")
  
  cat("\nFollow-up en años:\n")
  follow_up_after_years <- follow_up_after / 365
  print(summary(follow_up_after_years))
  
  cat("\nDistribución por categorías:\n")
  follow_up_after_cat <- cut(
    follow_up_after_years,
    breaks = c(0, 1, 3, 5, 10, Inf),
    labels = c("<1 año", "1-3 años", "3-5 años", "5-10 años", ">10 años"),
    include.lowest = TRUE
  )
  print(table(follow_up_after_cat, useNA = "ifany"))
  
} else {
  cat("No existe days_to_last_follow_up\n")
}

#------------------------------------------------#
# 14) ARCHIVOS GENERADOS                         #
#------------------------------------------------#
write_section("14) ARCHIVOS GENERADOS")

cat("\nArchivo RData final:\n")
cat("subsets/primary_u18/final_analysis/TARGET_ALL_P2_PRIMARY_U18_FINAL_FOR_ANALYSIS.RData\n")

cat("\nArchivo CSV final:\n")
cat("subsets/primary_u18/final_analysis/TARGET_ALL_P2_PRIMARY_U18_FINAL_METADATA.csv\n")

cat("\nArchivo TXT de reporte:\n")
cat(report_file, "\n")

cat("\nColumnas eliminadas para CSV por ser tipo list:\n")
list_cols <- names(meta_final)[sapply(meta_final, is.list)]
if (length(list_cols) > 0) {
  print(list_cols)
} else {
  cat("Ninguna\n")
}

sink()

cat("Reporte TXT guardado en:\n")
cat(report_file, "\n")
cat("\n=== REPORTE TXT COMPLETADO ===\n")