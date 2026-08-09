#!/usr/bin/env Rscript
# Single-cell stage 2: normalization, clustering, and cell-type annotation.
#
# NOT YET RUN in this project - see the note at the top of
# 01_qc_and_doublets.R. Depends on that stage's output.
#
# Usage: Rscript code/08_single_cell/02_normalize_cluster_annotate.R

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})
source(here::here("code", "utils", "paths.R"))

out_dir <- results_path("single_cell", "GSE136103")
obj <- readRDS(file.path(out_dir, "01_qc_singlets.rds"))

# --- normalization + dimensionality reduction --------------------------------
# SCTransform (regularised negative binomial regression) rather than
# log-normalisation + scaling: it models the mean-variance relationship of
# UMI counts directly and is the current Seurat-recommended default for
# datasets, like this one, assembled from several libraries prepared at
# different times (it is more robust to between-library differences in
# sequencing depth than log-normalisation).
obj <- SCTransform(obj, vars.to.regress = "percent_mt", verbose = FALSE)
obj <- RunPCA(obj, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.8, verbose = FALSE)

# --- marker genes per cluster --------------------------------------------------
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.4, verbose = FALSE)
write.csv(markers, file.path(out_dir, "02_cluster_markers.csv"), row.names = FALSE)

# --- cell-type annotation by marker-set scoring -------------------------------
# Shared with the spatial pipeline (code/09_spatial) - see
# code/utils/liver_cell_markers.R for why this is factored out rather than
# defined separately in each place.
source(here::here("code", "utils", "liver_cell_markers.R"))
cluster_labels <- annotate_clusters_by_markers(markers, cluster_col = "cluster", gene_col = "gene")
obj$cell_type <- cluster_labels[as.character(obj$seurat_clusters)]

composition <- obj@meta.data |>
  count(condition, cell_type) |>
  group_by(condition) |>
  mutate(pct = round(100 * n / sum(n), 2)) |>
  ungroup()
write.csv(composition, file.path(out_dir, "02_celltype_composition.csv"), row.names = FALSE)

saveRDS(obj, file.path(out_dir, "02_annotated.rds"))
cat(sprintf("Annotated %d clusters into %d cell types across %d cells\n",
            length(unique(obj$seurat_clusters)), length(unique(obj$cell_type)), ncol(obj)))
