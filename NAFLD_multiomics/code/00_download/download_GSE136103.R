#!/usr/bin/env Rscript
# Download GSE136103 - Ramachandran et al. 2019, human liver scRNA-seq
# (healthy + cirrhotic/NAFLD-spectrum liver, ~35,000 cells, 10x Genomics).
# Each library is deposited as a standard 10x triplet (barcodes/genes/matrix),
# bundled together in one RAW.tar.
# Run from the project root: Rscript code/00_download/download_GSE136103.R

suppressPackageStartupMessages(library(GEOquery))
source(here::here("code", "utils", "paths.R"))

out_dir <- dataset_path("single_cell", "GSE136103")
raw_dir <- file.path(out_dir, "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

gse  <- getGEO("GSE136103", GSEMatrix = TRUE, destdir = out_dir, AnnotGPL = FALSE)
meta <- pData(gse[[1]])
saveRDS(meta, file.path(out_dir, "meta_gse136103.rds"))
write.csv(meta, file.path(out_dir, "meta_gse136103.csv"), row.names = TRUE)

tar_path <- file.path(out_dir, "GSE136103_RAW.tar")
if (!file.exists(tar_path)) {
  cat("Downloading GSE136103_RAW.tar (~436 MB, this will take a few minutes)...\n")
  download.file(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136103/suppl/GSE136103_RAW.tar",
    destfile = tar_path, mode = "wb"
  )
}
untar(tar_path, exdir = raw_dir)

libraries <- unique(sub("_(barcodes|genes|matrix).*$", "", list.files(raw_dir)))
cat(sprintf("GSE136103: %d libraries extracted to %s\n", length(libraries), raw_dir))
