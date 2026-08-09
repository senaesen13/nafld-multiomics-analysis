# Curated liver cell-type marker genes, used for both single-cell
# (code/08_single_cell) and spatial (code/09_spatial) cluster annotation.
#
# Kept in one shared file rather than copied into each stage so that a
# cluster labelled "Kupffer cells" means the same thing - the same marker
# genes - whether it came from the scRNA-seq or the spatial pipeline. This
# is what makes it possible to compare the two later (e.g. "does the
# TREM2+ macrophage niche found spatially match the Kupffer cell cluster
# found in scRNA-seq?").
#
# Source: bioinformatics-learning/week6-day1-nafld-scrna-gse136103/scripts/
# 06_doublet_detection.R (established and validated on GSE136103 in this
# project's predecessor).

LIVER_CELL_MARKERS <- list(
  "CD4+ T cells"              = c("CD3D", "CD3G", "CD4", "IL7R", "LTB", "CD40LG"),
  "CD8+ T cells"               = c("CD8A", "CD8B", "CD3G"),
  "CD8+ T cells (exhausted)"    = c("CD8A", "LAG3", "CRTAM", "TIGIT", "HAVCR2"),
  "Naive T cells"                = c("MAL", "LEF1", "CCR7", "TCF7", "SELL"),
  "NK/NKT cells"                  = c("XCL1", "TOX2", "KLRC1", "KLRC2"),
  "NK cells (cytotoxic)"            = c("GNLY", "GZMB", "FGFBP2", "NKG7", "PRF1"),
  "NK cells (liver-resident)"        = c("IL2RB", "CD160", "NCR1", "CXCR6"),
  "gdT/NK cells"                       = c("GZMH", "TRGC2", "CX3CR1", "FCGR3A"),
  "Monocytes"                            = c("S100A8", "S100A9", "S100A12", "FCN1", "LYZ"),
  "Dendritic cells"                        = c("CD1C", "FCER1A", "CLEC10A", "HLA-DQA1"),
  "Kupffer cells"                            = c("C1QB", "C1QC", "CD5L", "CD163", "GPNMB", "TIMD4"),
  "B cells"                                    = c("CD79A", "CD19", "MS4A1", "BANK1"),
  "Plasma cells"                                 = c("IGHGP", "IGLL5", "IGHA2", "MZB1", "JCHAIN"),
  "Endothelial cells"                              = c("GPIHBP1", "PODXL", "AQP7", "PECAM1"),
  "LSEC"                                             = c("CLEC4G", "FCN2", "FCN3", "OIT3", "STAB2"),
  "Hepatic stellate cells"                             = c("DCN", "TCF21", "RGS5", "COL1A1", "ACTA2"),
  "Hepatocytes"                                          = c("UGT2B15", "UGT2A3", "CYP3A4", "ALB", "APOB"),
  "Proliferating cells"                                    = c("CENPA", "RRM2", "TYMS", "MKI67", "TOP2A"),
  "Cholangiocytes"                                           = c("SCT", "PTCRA", "LRRC26", "KRT7", "KRT19"),
  "Mast cells"                                                 = c("TPSAB1", "TPSB2", "CPA3", "GATA2")
)

#' Assign each cluster the cell type whose marker set best overlaps that
#' cluster's own top marker genes (from Seurat::FindAllMarkers output).
annotate_clusters_by_markers <- function(markers_df, cluster_col = "cluster", gene_col = "gene") {
  clusters <- unique(markers_df[[cluster_col]])
  labels <- sapply(clusters, function(cl) {
    cl_genes <- markers_df[[gene_col]][markers_df[[cluster_col]] == cl]
    scores <- sapply(LIVER_CELL_MARKERS, function(m) sum(m %in% cl_genes))
    if (max(scores) == 0) paste0("Unknown_", cl) else names(which.max(scores))
  })
  setNames(labels, clusters)
}
