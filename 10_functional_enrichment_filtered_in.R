### 10_functional_enrichment_filtered_interpretation_PRIMARY_U18_TARGET_ALL_P2.R
### GOAL:
###   1) Filtrar resultados del script 09 por significancia estadística objetiva
###   2) Generar tablas finales estrictas y exploratorias
###   3) Generar dotplots finales SOLO con términos significativos
###   4) Separar GSEA positivo/negativo para interpretación por cluster
###   5) Crear tablas resumen listas para interpretación biológica / paper
###
### INPUT:
###   - subsets/primary_u18/final_analysis/09_functional_enrichment/*.csv
###
### OUTPUT:
###   - subsets/primary_u18/final_analysis/10_functional_enrichment_filtered/
###       * strict_padj_0.05/*.csv
###       * exploratory_padj_0.10/*.csv
###       * final_dotplots/*.png
###       * top_terms_summary_*.csv
###       * functional_enrichment_filtered_report.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(forcats)
})

cat("=== SCRIPT 10: FILTERED FUNCTIONAL ENRICHMENT INTERPRETATION ===\n")

#---------------------------#
# 1) Paths y parámetros     #
#---------------------------#

indir <- "subsets/primary_u18/final_analysis/09_functional_enrichment"
outdir <- "subsets/primary_u18/final_analysis/10_functional_enrichment_filtered"

strict_dir <- file.path(outdir, "strict_padj_0.05")
exploratory_dir <- file.path(outdir, "exploratory_padj_0.10")
plot_dir <- file.path(outdir, "final_dotplots")

for (d in c(outdir, strict_dir, exploratory_dir, plot_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

strict_padj_cutoff <- 0.05
exploratory_padj_cutoff <- 0.10
plot_top_n <- 20
min_terms_to_plot <- 1

# Dotplot parameters
plot_width <- 11
plot_height <- 8
plot_dpi <- 300

#---------------------------#
# 2) Funciones auxiliares   #
#---------------------------#

safe_read_csv <- function(file) {
  if (!file.exists(file)) {
    warning("Archivo no encontrado: ", file)
    return(NULL)
  }

  df <- tryCatch({
    readr::read_csv(file, show_col_types = FALSE)
  }, error = function(e) {
    warning("No se pudo leer: ", file, " | ", conditionMessage(e))
    NULL
  })

  df
}

safe_write_csv <- function(df, file) {
  if (is.null(df)) return(invisible(NULL))
  readr::write_csv(as.data.frame(df), file)
}

make_safe_stub <- function(x) {
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

parse_enrichment_filename <- function(file) {
  fname <- basename(file)
  fname <- sub("\\.csv$", "", fname)

  analysis_type <- dplyr::case_when(
    startsWith(fname, "gsea_GO_BP_") ~ "GSEA_GO_BP",
    startsWith(fname, "gsea_GO_MF_") ~ "GSEA_GO_MF",
    startsWith(fname, "gsea_GO_CC_") ~ "GSEA_GO_CC",
    startsWith(fname, "gsea_KEGG_") ~ "GSEA_KEGG",
    startsWith(fname, "ora_GO_BP_up_") ~ "ORA_GO_BP_UP",
    startsWith(fname, "ora_GO_BP_down_") ~ "ORA_GO_BP_DOWN",
    startsWith(fname, "ora_KEGG_up_") ~ "ORA_KEGG_UP",
    startsWith(fname, "ora_KEGG_down_") ~ "ORA_KEGG_DOWN",
    TRUE ~ NA_character_
  )

  comparison <- fname %>%
    str_replace("^gsea_GO_BP_", "") %>%
    str_replace("^gsea_GO_MF_", "") %>%
    str_replace("^gsea_GO_CC_", "") %>%
    str_replace("^gsea_KEGG_", "") %>%
    str_replace("^ora_GO_BP_up_", "") %>%
    str_replace("^ora_GO_BP_down_", "") %>%
    str_replace("^ora_KEGG_up_", "") %>%
    str_replace("^ora_KEGG_down_", "")

  tibble(
    file = file,
    filename = basename(file),
    analysis_type = analysis_type,
    comparison = comparison
  )
}

add_direction_interpretation <- function(df, analysis_type, comparison) {
  if (is.null(df) || nrow(df) == 0) return(df)

  # Para GSEA:
  # NES positivo = enriquecido hacia el primer cluster de la comparación.
  # NES negativo = enriquecido hacia el segundo cluster de la comparación.
  cluster_a <- str_match(comparison, "cluster_([0-9]+)_vs_([0-9]+)")[, 2]
  cluster_b <- str_match(comparison, "cluster_([0-9]+)_vs_([0-9]+)")[, 3]

  if (grepl("^GSEA", analysis_type) && "NES" %in% colnames(df)) {
    df <- df %>%
      dplyr::mutate(
        enriched_toward = dplyr::case_when(
          NES > 0 ~ paste0("cluster_", cluster_a),
          NES < 0 ~ paste0("cluster_", cluster_b),
          TRUE ~ "neutral"
        ),
        direction = dplyr::case_when(
          NES > 0 ~ "positive_NES",
          NES < 0 ~ "negative_NES",
          TRUE ~ "neutral"
        )
      )
  } else if (grepl("_UP$", analysis_type)) {
    df <- df %>%
      dplyr::mutate(
        enriched_toward = paste0("cluster_", cluster_a),
        direction = "up_gene_ORA"
      )
  } else if (grepl("_DOWN$", analysis_type)) {
    df <- df %>%
      dplyr::mutate(
        enriched_toward = paste0("cluster_", cluster_b),
        direction = "down_gene_ORA"
      )
  }

  df
}

filter_enrichment <- function(df, padj_cutoff) {
  if (is.null(df) || nrow(df) == 0) return(df)

  if (!"p.adjust" %in% colnames(df)) {
    warning("La tabla no tiene columna p.adjust")
    return(df[0, , drop = FALSE])
  }

  df %>%
    dplyr::filter(!is.na(p.adjust), p.adjust < padj_cutoff) %>%
    dplyr::arrange(p.adjust)
}

make_final_dotplot <- function(df, outfile, title, top_n = 20) {
  if (is.null(df) || nrow(df) < min_terms_to_plot) {
    cat("[INFO] No hay términos suficientes para graficar:", title, "\n")
    return(NULL)
  }

  plot_df <- df %>%
    dplyr::arrange(p.adjust) %>%
    utils::head(top_n) %>%
    dplyr::mutate(
      Description = stringr::str_wrap(Description, width = 55),
      Description = forcats::fct_reorder(Description, -log10(p.adjust))
    )

  # Si existe NES, usarlo para mostrar dirección en GSEA.
  if ("NES" %in% colnames(plot_df)) {
    p <- ggplot(plot_df, aes(x = NES, y = Description)) +
      geom_point(aes(size = setSize, color = p.adjust), alpha = 0.9) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      scale_color_gradient(low = "red", high = "blue") +
      labs(
        title = title,
        x = "Normalized Enrichment Score (NES)",
        y = NULL,
        size = "Gene set size",
        color = "p.adjust"
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 10)
      )
  } else {
    # ORA: usar GeneRatio si existe; si viene como "12/100", convertir a numérico.
    if ("GeneRatio" %in% colnames(plot_df)) {
      plot_df <- plot_df %>%
        dplyr::mutate(
          GeneRatio_num = sapply(GeneRatio, function(x) {
            if (is.na(x)) return(NA_real_)
            if (grepl("/", x)) {
              parts <- strsplit(x, "/")[[1]]
              return(as.numeric(parts[1]) / as.numeric(parts[2]))
            }
            as.numeric(x)
          })
        )
      xvar <- "GeneRatio_num"
      xlabel <- "GeneRatio"
    } else {
      plot_df <- plot_df %>% dplyr::mutate(GeneRatio_num = -log10(p.adjust))
      xvar <- "GeneRatio_num"
      xlabel <- "-log10(p.adjust)"
    }

    size_var <- if ("Count" %in% colnames(plot_df)) "Count" else NULL

    p <- ggplot(plot_df, aes(x = .data[[xvar]], y = Description)) +
      geom_point(aes(size = if (!is.null(size_var)) .data[[size_var]] else 3,
                     color = p.adjust), alpha = 0.9) +
      scale_color_gradient(low = "red", high = "blue") +
      labs(
        title = title,
        x = xlabel,
        y = NULL,
        size = ifelse(is.null(size_var), "", "Count"),
        color = "p.adjust"
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 10)
      )
  }

  ggsave(outfile, plot = p, width = plot_width, height = plot_height, dpi = plot_dpi)
  cat("[OK] Dotplot final:", outfile, "\n")
}

#---------------------------#
# 3) Buscar archivos        #
#---------------------------#

cat("\n=== BUSCANDO RESULTADOS DEL SCRIPT 09 ===\n")

if (!dir.exists(indir)) {
  stop("No existe el directorio del script 09: ", indir)
}

all_csv <- list.files(indir, pattern = "\\.csv$", full.names = TRUE)

# Excluir listas rankeadas y summaries generales; aquí solo queremos tablas de enriquecimiento.
enrichment_csv <- all_csv[
  grepl("^(gsea_GO_|gsea_KEGG_|ora_GO_|ora_KEGG_)", basename(all_csv)) &
    !grepl("ranked_gene_list|summary", basename(all_csv))
]

if (length(enrichment_csv) == 0) {
  stop("No se encontraron tablas de enriquecimiento en: ", indir)
}

file_index <- dplyr::bind_rows(lapply(enrichment_csv, parse_enrichment_filename)) %>%
  dplyr::filter(!is.na(analysis_type))

write_csv(file_index, file.path(outdir, "input_enrichment_files_index.csv"))

cat("Archivos de enriquecimiento encontrados:", nrow(file_index), "\n")
print(file_index %>% dplyr::count(analysis_type))

#---------------------------#
# 4) Filtrar resultados     #
#---------------------------#

cat("\n=== FILTRANDO RESULTADOS ===\n")

summary_list <- list()
all_strict <- list()
all_exploratory <- list()

for (i in seq_len(nrow(file_index))) {
  meta <- file_index[i, ]
  df <- safe_read_csv(meta$file)
  if (is.null(df)) next

  df <- add_direction_interpretation(
    df = df,
    analysis_type = meta$analysis_type,
    comparison = meta$comparison
  )

  strict_df <- filter_enrichment(df, strict_padj_cutoff)
  exploratory_df <- filter_enrichment(df, exploratory_padj_cutoff)

  strict_outfile <- file.path(
    strict_dir,
    paste0(sub("\\.csv$", "", meta$filename), "_padj_lt_0.05.csv")
  )

  exploratory_outfile <- file.path(
    exploratory_dir,
    paste0(sub("\\.csv$", "", meta$filename), "_padj_lt_0.10.csv")
  )

  safe_write_csv(strict_df, strict_outfile)
  safe_write_csv(exploratory_df, exploratory_outfile)

  if (!is.null(strict_df) && nrow(strict_df) > 0) {
    strict_df$source_file <- meta$filename
    strict_df$analysis_type <- meta$analysis_type
    strict_df$comparison <- meta$comparison
    all_strict[[length(all_strict) + 1]] <- strict_df
  }

  if (!is.null(exploratory_df) && nrow(exploratory_df) > 0) {
    exploratory_df$source_file <- meta$filename
    exploratory_df$analysis_type <- meta$analysis_type
    exploratory_df$comparison <- meta$comparison
    all_exploratory[[length(all_exploratory) + 1]] <- exploratory_df
  }

  summary_list[[length(summary_list) + 1]] <- tibble(
    filename = meta$filename,
    analysis_type = meta$analysis_type,
    comparison = meta$comparison,
    total_terms = nrow(df),
    strict_terms_padj_lt_0.05 = nrow(strict_df),
    exploratory_terms_padj_lt_0.10 = nrow(exploratory_df)
  )

  # Figuras finales estrictas.
  if (nrow(strict_df) > 0) {
    plot_stub <- make_safe_stub(paste0(meta$analysis_type, "_", meta$comparison, "_strict_padj_0_05"))
    plot_file <- file.path(plot_dir, paste0("dotplot_", plot_stub, ".png"))
    plot_title <- paste0(meta$analysis_type, " | ", meta$comparison, " | p.adjust < 0.05")
    make_final_dotplot(strict_df, plot_file, plot_title, plot_top_n)
  }
}

summary_df <- dplyr::bind_rows(summary_list)
write_csv(summary_df, file.path(outdir, "filtering_summary_by_file.csv"))

strict_master <- if (length(all_strict) > 0) dplyr::bind_rows(all_strict) else tibble()
exploratory_master <- if (length(all_exploratory) > 0) dplyr::bind_rows(all_exploratory) else tibble()

write_csv(strict_master, file.path(outdir, "all_significant_terms_strict_padj_lt_0.05.csv"))
write_csv(exploratory_master, file.path(outdir, "all_significant_terms_exploratory_padj_lt_0.10.csv"))

#---------------------------#
# 5) Tablas top finales     #
#---------------------------#

cat("\n=== GENERANDO TABLAS TOP ===\n")

if (nrow(strict_master) > 0) {
  top_by_analysis <- strict_master %>%
    dplyr::group_by(comparison, analysis_type) %>%
    dplyr::arrange(p.adjust, .by_group = TRUE) %>%
    dplyr::slice_head(n = 20) %>%
    dplyr::ungroup()

  write_csv(top_by_analysis, file.path(outdir, "top20_strict_terms_by_comparison_and_analysis.csv"))

  # Tabla especialmente útil para GSEA: GO BP + KEGG, separada por dirección.
  top_gsea_interpretable <- strict_master %>%
    dplyr::filter(analysis_type %in% c("GSEA_GO_BP", "GSEA_KEGG")) %>%
    dplyr::arrange(comparison, enriched_toward, p.adjust) %>%
    dplyr::group_by(comparison, analysis_type, enriched_toward) %>%
    dplyr::slice_head(n = 15) %>%
    dplyr::ungroup()

  write_csv(top_gsea_interpretable, file.path(outdir, "top_GSEA_GO_BP_KEGG_by_cluster_direction_strict.csv"))

  # Conteo de términos por cluster/dirección.
  direction_summary <- strict_master %>%
    dplyr::count(comparison, analysis_type, enriched_toward, name = "n_significant_terms") %>%
    dplyr::arrange(comparison, analysis_type, enriched_toward)

  write_csv(direction_summary, file.path(outdir, "significant_terms_count_by_direction.csv"))
} else {
  cat("[INFO] No hubo términos estrictos con p.adjust < 0.05.\n")
}

#---------------------------#
# 6) Reporte                #
#---------------------------#

cat("\n=== GENERANDO REPORTE TXT ===\n")

report_file <- file.path(outdir, "functional_enrichment_filtered_report.txt")
sink(report_file)

cat("REPORTE DE FILTRADO E INTERPRETACIÓN DE ENRIQUECIMIENTO FUNCIONAL\n")
cat("Fecha:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("INPUT DIRECTORY:\n")
cat(indir, "\n\n")

cat("OUTPUT DIRECTORY:\n")
cat(outdir, "\n\n")

cat("CRITERIOS DE FILTRADO:\n")
cat("- Estricto: p.adjust <", strict_padj_cutoff, "\n")
cat("- Exploratorio: p.adjust <", exploratory_padj_cutoff, "\n\n")

cat("INTERPRETACIÓN RECOMENDADA:\n")
cat("- Para resultados principales/paper: usar tablas strict_padj_0.05 y figuras en final_dotplots.\n")
cat("- Para discusión exploratoria: revisar exploratory_padj_0.10, pero no presentarlo como evidencia definitiva.\n")
cat("- No seleccionar términos solo por color rojo/azul del gráfico; seleccionar por p.adjust numérico.\n")
cat("- En GSEA, NES positivo indica enriquecimiento hacia el primer cluster de la comparación.\n")
cat("- En GSEA, NES negativo indica enriquecimiento hacia el segundo cluster de la comparación.\n")
cat("- En ORA Up, los términos corresponden a genes aumentados hacia el primer cluster de la comparación.\n")
cat("- En ORA Down, los términos corresponden a genes aumentados hacia el segundo cluster de la comparación.\n\n")

cat("RESUMEN POR ARCHIVO:\n")
print(summary_df)
cat("\n")

cat("RESUMEN GLOBAL:\n")
cat("- Archivos analizados:", nrow(file_index), "\n")
cat("- Términos estrictos totales:", nrow(strict_master), "\n")
cat("- Términos exploratorios totales:", nrow(exploratory_master), "\n\n")

if (nrow(strict_master) > 0) {
  cat("TOP 30 TÉRMINOS ESTRICTOS GLOBALES:\n")
  print(
    strict_master %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::select(any_of(c("comparison", "analysis_type", "enriched_toward", "ID", "Description", "NES", "GeneRatio", "Count", "p.adjust"))) %>%
      utils::head(30)
  )
  cat("\n")
}

cat("ARCHIVOS PRINCIPALES GENERADOS:\n")
cat("- filtering_summary_by_file.csv\n")
cat("- all_significant_terms_strict_padj_lt_0.05.csv\n")
cat("- all_significant_terms_exploratory_padj_lt_0.10.csv\n")
cat("- top20_strict_terms_by_comparison_and_analysis.csv\n")
cat("- top_GSEA_GO_BP_KEGG_by_cluster_direction_strict.csv\n")
cat("- significant_terms_count_by_direction.csv\n")
cat("- final_dotplots/*.png\n")
cat("- functional_enrichment_filtered_report.txt\n")

sink()

cat("\n=== SCRIPT 10 COMPLETADO ===\n")
cat("Archivos generados en:\n")
cat(outdir, "\n")
