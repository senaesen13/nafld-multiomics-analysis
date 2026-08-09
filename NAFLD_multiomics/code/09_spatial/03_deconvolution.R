#!/usr/bin/env Rscript
# Spatial stage 3: cell-type deconvolution via RCTD.
#
# NOT YET RUN in this project - see the note at the top of
# 01_qc_and_normalize.R. Also depends on the single-cell pipeline
# (code/08_single_cell) having been run, since GSE136103's annotated cells
# are the reference this stage deconvolves against.
#
# Usage: Rscript code/09_spatial/03_deconvolution.R
#
# Why RCTD (spacexr) instead of Seurat's AddModuleScore, which is what the
# predecessor project (bioinformatics-learning/week8-day1-masld-spatial-setup)
# fell back to after CARD failed to compile on macOS? AddModuleScore reports
# a relative activity score per cell type per spot, not a compositional
# proportion - it cannot say "this spot is 40% hepatocyte, 30% Kupffer
# cell", only "hepatocyte markers score higher here than there". Because
# hepatocytes dominate liver transcriptionally by sheer abundance, every
# spot's highest AddModuleScore is "Hepatocytes" regardless of true
# composition (confirmed in the predecessor project: all 24/24 clusters
# there were called hepatocyte-dominant, and the coefficient of variation of
# the hepatocyte score across clusters was the *lowest* of any cell type -
# 13% vs 25-55% for others - meaning it is structurally the least
# discriminating score, not the most biologically dominant one). RCTD
# instead fits how much of each spot's transcriptome is explained by each
# reference cell type's expression profile, giving proportions that sum to
# 1 and can distinguish "more hepatocyte than usual" from "hepatocyte
# transcriptionally swamping everything, as always in liver". `spacexr` is
# confirmed installed in this environment (unlike CARD).
#
# API note: `spacexr`'s installed version here exposes createReference() /
# createSpatialRNA() / createRctd() / runRctd() (camelCase), not the
# create.RCTD()/run.RCTD() names most RCTD tutorials online still show -
# this script uses the confirmed-installed names. A synthetic smoke test
# during development found that createRctd() expects a SummarizedExperiment/
# SpatialExperiment-shaped input, and it was not established within this
# session's scope exactly how createSpatialRNA()'s/createReference()'s own
# output should be converted to that shape (or whether a different
# constructor path is intended). Resolve this - check
# `?spacexr::createRctd` and the package vignette - before this stage is
# first run; it is flagged here rather than guessed at.

suppressPackageStartupMessages({
  library(Seurat)
  library(spacexr)
  library(dplyr)
})
source(here::here("code", "utils", "paths.R"))

spatial_dir <- results_path("spatial", "Vu2025_Visium")
sc_dir      <- results_path("single_cell", "GSE136103")

spatial_obj <- readRDS(file.path(spatial_dir, "02_clustered_annotated.rds"))
ref_obj     <- readRDS(file.path(sc_dir, "02_annotated.rds"))

# --- build the RCTD reference from the annotated scRNA-seq object ----------
ref_counts <- GetAssayData(ref_obj, assay = "RNA", layer = "counts")
ref_cell_types <- factor(ref_obj$cell_type)
names(ref_cell_types) <- colnames(ref_obj)
reference <- createReference(ref_counts, ref_cell_types)

# --- build the spatial query from the Visium object -------------------------
spatial_counts <- GetAssayData(spatial_obj, assay = "Spatial", layer = "counts")
coords <- GetTissueCoordinates(spatial_obj)[, c("x", "y")]
spatial_query <- createSpatialRNA(coords, spatial_counts)

# "full" mode fits an arbitrary number of cell types per spot (not just up
# to 2, as "doublet" mode assumes) - appropriate here since each ~55um
# Visium spot captures roughly 3-10 cells of potentially several types.
rctd <- createRctd(spatial_query, reference)
rctd <- runRctd(rctd, rctd_mode = "full")

weights <- rctd@results$weights  # spot x cell-type proportion matrix, if this is where the fitted output lives - verify against the actual returned object
write.csv(as.data.frame(as.matrix(weights)), file.path(spatial_dir, "03_rctd_celltype_proportions.csv"))

dominant_type <- colnames(weights)[apply(weights, 1, which.max)]
spatial_obj$rctd_dominant_type <- dominant_type

write.csv(
  data.frame(spot = colnames(spatial_obj), array = spatial_obj$array,
             cluster = spatial_obj$seurat_clusters, rctd_dominant_type = dominant_type),
  file.path(spatial_dir, "03_rctd_dominant_per_spot.csv"), row.names = FALSE
)
saveRDS(rctd, file.path(spatial_dir, "03_rctd_full_result.rds"))
cat("Saved RCTD deconvolution results to", spatial_dir, "\n")
