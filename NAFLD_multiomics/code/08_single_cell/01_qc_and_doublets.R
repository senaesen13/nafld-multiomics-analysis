#!/usr/bin/env Rscript
# Single-cell stage 1: load, QC, and doublet removal - GSE136103 (Ramachandran
# et al. 2019, human liver scRNA-seq, 10x Genomics).
#
# NOT YET RUN in this project - written and verified against the real
# downloaded metadata (data/single_cell/GSE136103/meta_gse136103.csv), but
# not executed, per instruction. Requires `scDblFinder`
# (BiocManager::install("scDblFinder")), which is not yet installed.
#
# Usage: Rscript code/08_single_cell/01_qc_and_doublets.R

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(dplyr)
})
source(here::here("code", "utils", "paths.R"))

raw_dir <- dataset_path("single_cell", "GSE136103", "raw")
out_dir <- results_path("single_cell", "GSE136103")

# --- restrict to the NAFLD-vs-healthy comparison ----------------------------
# GSE136103 is not an NAFLD-specific dataset - it profiles cirrhotic liver of
# five different etiologies plus healthy controls plus peripheral blood
# (PBMC) samples. Per this project's control-vs-disease convention (NAFLD is
# the disease of interest), only two groups are usable here:
#   - "healthy": normal liver, no disease
#   - "cause of liver disease: NAFLD": the only cirrhosis etiology that is
#     actually NAFLD (the metadata also separately lists Alcohol, Hereditary
#     Haemochromatosis, and PBC as causes in other patients - these must be
#     excluded, they are liver disease but not NAFLD)
# PBMC samples are excluded regardless of their disease label: they are
# peripheral blood, not liver tissue, and this project's blood-donor IDs do
# not reliably correspond to the same patients as the liver biopsies (3 of 4
# PBMC samples are labelled "NAFLD" cause despite only 2 of 5 cirrhotic
# *liver* donors having that etiology) - a data-provenance mismatch worth
# not building on.
meta <- read.csv(file.path(dataset_path("single_cell", "GSE136103"), "meta_gse136103.csv"),
                  row.names = 1, check.names = FALSE)
meta$sample_id <- meta$title
meta$condition <- ifelse(meta$"disease status:ch1" == "healthy", "Healthy", "NAFLD")

keep <- meta$"population:ch1" != "PBMC" &
  (meta$"disease status:ch1" == "healthy" | meta$"cause of liver disease:ch1" == "NAFLD")
meta_keep <- meta[keep, ]
cat(sprintf("Keeping %d / %d libraries (%d Healthy, %d NAFLD), dropping PBMC and non-NAFLD cirrhosis\n",
            nrow(meta_keep), nrow(meta), sum(meta_keep$condition == "Healthy"),
            sum(meta_keep$condition == "NAFLD")))

# --- load each kept library as a 10x triplet and tag with condition/donor ---
donor_of <- function(title) sub("_(Cd45|CD45).*", "", title, ignore.case = TRUE)

seurat_list <- lapply(seq_len(nrow(meta_keep)), function(i) {
  lib_title <- meta_keep$title[i]
  gsm       <- meta_keep$geo_accession[i]
  # All 26 libraries' 10x triplets sit flat in one directory, named
  # <GSM>_<samplename>_{barcodes,genes,matrix}.*.gz (confirmed from the
  # actual extracted files - this is GEO's deposit layout, not the standard
  # per-library-subfolder layout Read10X() expects), so ReadMtx() with
  # explicit per-file paths is required instead of Read10X(data.dir=...).
  lib_files <- list.files(raw_dir, pattern = paste0("^", gsm, "_"), full.names = TRUE)
  counts <- ReadMtx(
    mtx      = grep("matrix",   lib_files, value = TRUE),
    cells    = grep("barcodes", lib_files, value = TRUE),
    features = grep("genes",    lib_files, value = TRUE)
  )
  obj <- CreateSeuratObject(counts = counts, project = lib_title, min.cells = 3, min.features = 200)
  obj$library   <- lib_title
  obj$donor     <- donor_of(lib_title)
  obj$condition <- meta_keep$condition[i]
  obj
})
merged <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = meta_keep$title)

# --- QC filtering -------------------------------------------------------------
merged[["percent_mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")
merged <- subset(merged, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent_mt < 20)
cat(sprintf("After QC: %d cells across %d libraries\n", ncol(merged), length(unique(merged$library))))

# --- doublet detection (per library, as doublet rate depends on loading
#     density which differs by library) --------------------------------------
sce <- as.SingleCellExperiment(merged)
sce <- scDblFinder(sce, samples = "library")
merged$doublet_class <- sce$scDblFinder.class
merged_clean <- subset(merged, subset = doublet_class == "singlet")
cat(sprintf("After doublet removal: %d cells (%d doublets removed, %.1f%%)\n",
            ncol(merged_clean), sum(merged$doublet_class == "doublet"),
            100 * mean(merged$doublet_class == "doublet")))

saveRDS(merged_clean, file.path(out_dir, "01_qc_singlets.rds"))
write.csv(meta_keep, file.path(out_dir, "01_libraries_used.csv"), row.names = FALSE)
cat("Saved QC'd, doublet-free Seurat object to", out_dir, "\n")
