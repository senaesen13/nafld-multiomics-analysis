# Project-wide path helpers.
#
# Every analysis script starts with source("code/utils/paths.R") (or a
# relative path to it) and then refers to data with dataset_path(), and
# writes output with results_path(). This means a script never hardcodes
# "/Users/..." anywhere, and the whole project can be moved, cloned, or
# run on a different machine without editing a single file path.
#
# how it works: the `here` package finds the project root by walking up
# from the current file until it finds the `.here` marker file that sits
# in NAFLD_multiomics/. Every path below is built from that root, so
# scripts can be run from any working directory.

suppressPackageStartupMessages(library(here))

#' Path to a dataset's raw data folder, e.g. dataset_path("RNAseq", "GSE162694")
dataset_path <- function(...) here::here("data", ...)

#' Path to a results folder, created if it doesn't exist yet
results_path <- function(..., create = TRUE) {
  p <- here::here("results", ...)
  if (create) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  p
}

# The bulk RNA-seq cohorts used throughout the project, with a short label
# used in plots/tables and the one design choice (which sample group counts
# as the control arm) that differs between them.
#
# Originally 4 cohorts (GSE162694, GSE135251, GSE126848, EMTAB12807), then
# narrowed to GSE135251+GSE126848 (largest raw DEG overlap, 178 genes), then
# swapped to this pairing instead: GSE162694+GSE135251 has slightly less raw
# overlap (138 genes) but is the strongest pair on every cross-cohort
# consistency metric - genome-wide LFC correlation (0.382, vs -0.065 for the
# GSE135251+GSE126848 pairing), direction concordance (73.2% vs 65.2%), and
# Fisher's exact test significance (9.2e-54 vs 4.4e-30). It is also Sena's
# own original discovery+validation pairing. See
# backup_data/results/RNAseq/cross_cohort_4dataset_archive/ for the full
# original 6-pair comparison this was decided from. GSE126848 and
# EMTAB12807 are archived, not deleted, at backup_data/. To restore a
# cohort, move its data/RNAseq, results/RNAseq, and
# code/00_download/download_<id>.R folders back out of backup_data/ and add
# its entry back here.
RNASEQ_COHORTS <- list(
  GSE162694 = list(label = "Discovery (GSE162694)",  control_group = "Normal"),
  GSE135251 = list(label = "Validation (GSE135251)", control_group = "Normal")
)
