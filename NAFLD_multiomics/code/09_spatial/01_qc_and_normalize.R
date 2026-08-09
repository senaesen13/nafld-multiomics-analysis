#!/usr/bin/env Rscript
# Spatial stage 1: per-array QC and normalization - Vu et al. 2025 MASLD
# Visium CytAssist dataset.
#
# NOT YET RUN in this project - the input data itself is not downloaded yet
# (data/spatial/Vu2025_Visium/ requires a manual step, see
# code/00_download/download_spatial_Vu2025.R). This script is written and
# ready to run once that data is in place.
#
# Usage: Rscript code/09_spatial/01_qc_and_normalize.R
#
# QC thresholds and the per-array-then-merge normalization strategy match
# what was already established and validated on this same dataset in this
# project's predecessor (bioinformatics-learning/week8-day1-masld-spatial-setup) -
# kept identical here rather than re-derived, since that choice was already
# reasoned through there (see that project's NOTES.md for the full
# rationale on why per-array SCTransform, not post-merge).

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})
source(here::here("code", "utils", "paths.R"))

raw_dir <- dataset_path("spatial", "Vu2025_Visium", "raw")
out_dir <- results_path("spatial", "Vu2025_Visium")

MIN_GENES_PER_SPOT <- 200   # CytAssist FFPE has lower UMI depth than fresh-frozen
MAX_PCT_MITO       <- 10    # compromised/dying cells or FFPE extraction artefact above this
MIN_SPOTS_PER_GENE <- 3     # removes probe artefacts and very lowly detected transcripts

array_dirs <- list.dirs(raw_dir, recursive = FALSE)
if (length(array_dirs) == 0) {
  stop("No array folders found in ", raw_dir, " - run code/00_download/download_spatial_Vu2025.R ",
       "first (it requires a manual download step, see the comment at the top of that script).")
}

qc_summary <- list()
seurat_list <- lapply(array_dirs, function(array_dir) {
  array_id <- basename(array_dir)
  obj <- Load10X_Spatial(array_dir)
  obj$array <- array_id
  n_before <- ncol(obj)

  obj[["percent_mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj <- subset(obj, subset = nFeature_Spatial >= MIN_GENES_PER_SPOT & percent_mt < MAX_PCT_MITO)
  obj <- obj[rowSums(GetAssayData(obj, layer = "counts") > 0) >= MIN_SPOTS_PER_GENE, ]

  qc_summary[[array_id]] <<- data.frame(
    array = array_id, spots_before = n_before, spots_after = ncol(obj),
    pct_kept = round(100 * ncol(obj) / n_before, 1),
    median_genes = median(obj$nFeature_Spatial), median_umis = median(obj$nCount_Spatial)
  )

  # SCTransform per array (not after merging) so that per-array sequencing-
  # depth normalisation does not get confounded with genuine cross-array
  # biological variation - see the predecessor project's NOTES.md for the
  # full reasoning.
  SCTransform(obj, assay = "Spatial", verbose = FALSE)
})

write.csv(bind_rows(qc_summary), file.path(out_dir, "01_qc_summary.csv"), row.names = FALSE)

features <- SelectIntegrationFeatures(seurat_list, nfeatures = 3000)
merged <- merge(seurat_list[[1]], y = seurat_list[-1])
VariableFeatures(merged) <- features

saveRDS(merged, file.path(out_dir, "01_normalized_merged.rds"))
cat(sprintf("Merged %d arrays: %d spots x %d genes, %d variable features\n",
            length(array_dirs), ncol(merged), nrow(merged), length(features)))
