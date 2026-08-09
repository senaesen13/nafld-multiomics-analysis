#!/usr/bin/env Rscript
# Download GSE162694 - NAFLD liver bulk RNA-seq, discovery cohort.
# Run from the project root: Rscript code/00_download/download_GSE162694.R

suppressPackageStartupMessages(library(GEOquery))
source(here::here("code", "utils", "paths.R"))

out_dir <- dataset_path("RNAseq", "GSE162694")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Raw counts matrix (already gene-level, no alignment/quantification needed)
counts_url  <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162694/suppl/GSE162694_raw_counts.csv.gz"
counts_dest <- file.path(out_dir, "GSE162694_raw_counts.csv.gz")
if (!file.exists(counts_dest)) download.file(counts_url, counts_dest, mode = "wb")

# Sample-level metadata (condition, age, sex, etc.) from the GEO series matrix
gse  <- getGEO("GSE162694", GSEMatrix = TRUE, destdir = out_dir, AnnotGPL = FALSE)
meta <- pData(gse[[1]])
saveRDS(meta, file.path(out_dir, "meta_gse162694.rds"))
write.csv(meta, file.path(out_dir, "meta_gse162694.csv"), row.names = TRUE)

cat(sprintf("GSE162694: %d samples, counts at %s\n", nrow(meta), counts_dest))
