#!/usr/bin/env Rscript
# Spatial stage 2: dimensionality reduction, clustering, and biological
# annotation of spatial domains.
#
# NOT YET RUN in this project - see the note at the top of
# 01_qc_and_normalize.R. Depends on that stage's output.
#
# Usage: Rscript code/09_spatial/02_clustering_and_annotation.R

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})
source(here::here("code", "utils", "paths.R"))
source(here::here("code", "utils", "liver_cell_markers.R"))

out_dir <- results_path("spatial", "Vu2025_Visium")
obj <- readRDS(file.path(out_dir, "01_normalized_merged.rds"))

# PCA/UMAP/clustering parameters match what was already validated on this
# same dataset in the predecessor project (30 PCs, 20 dims for
# neighbours/UMAP, resolution 0.5 - chosen there from an elbow plot showing
# eigenvalues levelling around PC15-18).
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)

markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE)
write.csv(markers, file.path(out_dir, "02_cluster_markers.csv"), row.names = FALSE)

# Same shared marker-set scoring as the single-cell pipeline (code/utils/
# liver_cell_markers.R) - a spatial cluster and a scRNA-seq cluster get the
# same label from the same evidence, so results from the two modalities can
# be compared directly (e.g. "is the TREM2+/SPP1+ macrophage niche seen in
# scRNA-seq spatially co-localised with fibrotic regions?").
# NB: a caveat that does not apply to single-cell - each Visium spot
# captures several cells (~3-10, mixed cell types), so a cluster labelled
# e.g. "Hepatocytes" reflects the dominant transcriptional signal in that
# spot, not a pure hepatocyte population the way a single-cell cluster
# would. Read cluster labels here as "hepatocyte-dominated" rather than
# "hepatocyte-only".
cluster_labels <- annotate_clusters_by_markers(markers, cluster_col = "cluster", gene_col = "gene")
obj$cell_type_dominant <- cluster_labels[as.character(obj$seurat_clusters)]

write.csv(
  obj@meta.data |> count(array, seurat_clusters, cell_type_dominant),
  file.path(out_dir, "02_cluster_by_array.csv"), row.names = FALSE
)
saveRDS(obj, file.path(out_dir, "02_clustered_annotated.rds"))
cat(sprintf("%d spatial domains identified across %d spots\n",
            length(unique(obj$seurat_clusters)), ncol(obj)))
