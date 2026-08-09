#!/usr/bin/env Rscript
# Download GSE126848 - Suppli et al. 2019, NAFLD liver bulk RNA-seq.
# 26 healthy controls (14 normal-weight + 12 obese) vs 15 NAFL + 16 NASH -
# added as a properly-powered replacement for the under-controlled third
# cohort used earlier in this project.
# Run from the project root: Rscript code/00_download/download_GSE126848.R

suppressPackageStartupMessages(library(GEOquery))
source(here::here("code", "utils", "paths.R"))

out_dir <- dataset_path("RNAseq", "GSE126848")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts_url  <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126848/suppl/GSE126848_Gene_counts_raw.txt.gz"
counts_dest <- file.path(out_dir, "GSE126848_Gene_counts_raw.txt.gz")
if (!file.exists(counts_dest)) download.file(counts_url, counts_dest, mode = "wb")

gse  <- getGEO("GSE126848", GSEMatrix = TRUE, destdir = out_dir, AnnotGPL = FALSE)
meta <- pData(gse[[1]])
saveRDS(meta, file.path(out_dir, "meta_gse126848.rds"))
write.csv(meta, file.path(out_dir, "meta_gse126848.csv"), row.names = TRUE)

cat(sprintf("GSE126848: %d samples, counts at %s\n", nrow(meta), counts_dest))
