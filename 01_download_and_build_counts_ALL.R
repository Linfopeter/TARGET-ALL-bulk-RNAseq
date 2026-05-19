### 01_download_and_build_counts_ALL.R
### GOAL: Download TARGET-ALL-P2 STAR counts and build counts matrix

# BiocManager::install("TCGAbiolinks")
# BiocManager::install("SummarizedExperiment")

library(TCGAbiolinks)
library(dplyr)
library(readr)

#---------------------------#
# 1) GDC query + download   #
#---------------------------#

query <- GDCquery(
  project = "TARGET-ALL-P2",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"   # o "HTSeq - Counts" si cambias
)

GDCdownload(query)

#---------------------------#
# 2) Locate downloaded TSVs #
#---------------------------#

files <- list.files(
  "GDCdata/TARGET-ALL-P2/Transcriptome_Profiling/Gene_Expression_Quantification",
  pattern = "tsv",
  full.names = TRUE,
  recursive = TRUE
)

cat("N archivos encontrados:", length(files), "\n")
head(files)

#---------------------------#
# 3) Function to read one   #
#---------------------------#

read_counts <- function(file) {
  df <- read_tsv(
    file,
    comment = "#",
    col_names = TRUE,
    show_col_types = FALSE
  )
  
  df %>% 
    select(gene_id, unstranded)
}

# Test on first file
one_file <- files[1]
test_read <- read_counts(one_file)
head(test_read)

genes <- test_read$gene_id
length(genes)
head(genes)

#---------------------------#
# 4) Build counts matrix    #
#---------------------------#

mat <- matrix(
  NA_integer_,
  nrow = length(genes),
  ncol = length(files),
  dimnames = list(genes, basename(files))
)

for (i in seq_along(files)) {
  message("Leyendo archivo ", i, " de ", length(files), " ...")
  
  df <- read_counts(files[i])
  
  if (!all(df$gene_id == genes)) {
    stop("Error: los genes no están en el mismo orden en el archivo ", files[i])
  }
  
  mat[, i] <- df$unstranded
}

counts <- mat

cat("Dimensiones de counts:\n")
print(dim(counts))

str(counts)


#---------------------------#
# 5) Save objects           #
#---------------------------#

save(counts, file = "TARGET_ALL_P2_STAR_counts.RData")
cat("Archivo guardado: TARGET_ALL_P2_STAR_counts.RData\n")
