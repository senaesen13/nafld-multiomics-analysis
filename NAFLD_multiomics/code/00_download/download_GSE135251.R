#!/usr/bin/env Rscript
# Download GSE135251 - NAFLD liver bulk RNA-seq, validation cohort.
# Unlike GSE162694/GSE126848, this series deposits per-sample HTSeq count
# files inside one RAW.tar rather than a single assembled matrix, so this
# script also builds the gene x sample count matrix.
# Run from the project root: Rscript code/00_download/download_GSE135251.R

suppressPackageStartupMessages(library(GEOquery))
source(here::here("code", "utils", "paths.R"))

out_dir <- dataset_path("RNAseq", "GSE135251")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- sample metadata ---------------------------------------------------
gse  <- getGEO("GSE135251", GSEMatrix = TRUE, destdir = out_dir, AnnotGPL = FALSE)
meta <- pData(gse[[1]])
saveRDS(meta, file.path(out_dir, "meta_gse135251.rds"))
write.csv(meta, file.path(out_dir, "meta_gse135251.csv"), row.names = TRUE)

# --- per-sample HTSeq count files (bundled as one tar) ------------------
tar_path <- file.path(out_dir, "GSE135251_RAW.tar")
if (!file.exists(tar_path)) {
  download.file(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135251/suppl/GSE135251_RAW.tar",
    destfile = tar_path, mode = "wb"
  )
}
untar(tar_path, exdir = out_dir)
count_files <- list.files(out_dir, pattern = "\\.counts\\.txt\\.gz$", full.names = TRUE)
cat("Count files found:", length(count_files), "\n")

# --- assemble into one gene x sample matrix ------------------------------
# HTSeq output has a few summary rows (__no_feature, __ambiguous, ...)
# that are not genes and must be dropped before building the matrix.
gsm_ids  <- sub("_(.*)\\.counts\\.txt\\.gz$", "", basename(count_files))
first    <- read.table(gzfile(count_files[1]), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
gene_ids <- first[!grepl("^__", first[, 1]), 1]

count_mat <- matrix(0L, nrow = length(gene_ids), ncol = length(count_files),
                     dimnames = list(gene_ids, gsm_ids))
for (i in seq_along(count_files)) {
  d <- read.table(gzfile(count_files[i]), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  d <- d[!grepl("^__", d[, 1]), ]
  count_mat[d[, 1], i] <- as.integer(d[, 2])
}

write.csv(count_mat, file.path(out_dir, "GSE135251_count_matrix.csv"))
cat(sprintf("GSE135251: %d genes x %d samples, matrix at %s\n",
            nrow(count_mat), ncol(count_mat),
            file.path(out_dir, "GSE135251_count_matrix.csv")))
