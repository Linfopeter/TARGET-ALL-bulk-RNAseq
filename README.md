# TARGET-BALL-bulk-RNAseq

Reproducible bulk RNA-seq pipeline for pediatric B-cell acute lymphoblastic leukemia using TARGET-ALL-P2 STAR-count data.

## Overview

This repository contains a complete R/Bioconductor workflow for downloading, curating, normalizing, and analyzing TARGET-ALL-P2 bulk RNA-seq data. The pipeline includes metadata harmonization, pediatric and primary-sample filtering, DESeq2 normalization, unsupervised transcriptomic stratification, differential expression analysis, pathway enrichment, and survival modeling.

## Main analyses

- GDC/TARGET data acquisition using TCGAbiolinks
- Metadata harmonization and quality control
- Pediatric cohort filtering
- Primary sample selection
- DESeq2 normalization and variance-stabilizing transformation
- PCA, UMAP, and k-means clustering
- Differential expression analysis
- Functional enrichment analysis
- Kaplan–Meier and Cox survival modeling

## Data availability

Raw TARGET-ALL-P2 RNA-seq data are publicly available through the Genomic Data Commons. This repository does not redistribute raw sequencing files or processed count matrices.

## Reproducibility

Run scripts sequentially from `01_download_and_build_counts_ALL.R` to downstream analysis scripts.
