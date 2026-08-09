#!/usr/bin/env Rscript
# Download the Vu et al. 2025 MASLD Visium CytAssist spatial transcriptomics
# dataset (JHEP Reports, DOI 10.1016/j.jhepr.2025.101657).
#
# MANUAL STEP REQUIRED: unlike the RNA-seq/single-cell datasets in this
# project, the UQ eSpace repository (DOI 10.48610/e95155f) returns 403
# Forbidden on a direct URL - the item page requires a click-through before
# the file becomes downloadable. This script cannot complete that step
# automatically. To populate this folder:
#   1. Open https://doi.org/10.48610/e95155f in a browser
#   2. Download Visium.zip (~336 MB) - CODEX.zip is not used by this project
#   3. Unzip it into data/spatial/Vu2025_Visium/raw/
# This script then organises whatever per-array folders it finds.

source(here::here("code", "utils", "paths.R"))

out_dir <- dataset_path("spatial", "Vu2025_Visium")
raw_dir <- file.path(out_dir, "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

zip_path <- file.path(out_dir, "Visium.zip")
if (!file.exists(zip_path) && length(list.files(raw_dir)) == 0) {
  stop(
    "Visium.zip not found. This dataset requires a manual download step ",
    "(see comment at the top of this script) - open ",
    "https://doi.org/10.48610/e95155f in a browser, download Visium.zip, ",
    "and place it at ", zip_path, " before re-running this script."
  )
}

if (file.exists(zip_path) && length(list.files(raw_dir)) == 0) {
  unzip(zip_path, exdir = raw_dir)
}

arrays <- list.dirs(raw_dir, recursive = FALSE)
cat(sprintf("Vu2025 spatial: %d array folders found in %s\n", length(arrays), raw_dir))
