#!/usr/bin/env Rscript
# DEG overlap across all 4 original cohorts using degs_list_overlap()
# (code/utils/degs_list_overlap.R), extracted from
# /Users/k2254978/Desktop/Work/01_Projects/Microarray_Urbani_Clean/source_code/Codes/DEGs_list_overlap.R -
# the original glioblastoma TCGA/GEO/CGGA survival-analysis project this
# was written for, not an approximation of it.
#
# Usage: Rscript code/05_cross_cohort/degs_list_overlap_4way.R
#
# Reads all 4 cohorts directly from wherever they currently sit (active
# results/RNAseq/ or archived backup_data/results/RNAseq/) without moving
# any files.

suppressPackageStartupMessages(library(dplyr))
source(here::here("code", "utils", "paths.R"))
source(here::here("code", "utils", "degs_list_overlap.R"))

ALL_FOUR <- c("GSE162694", "GSE135251", "GSE126848", "EMTAB12807")

locate_results <- function(ds) {
  active_path <- here::here("results", "RNAseq", ds, "deseq2_object.rds")
  backup_path <- here::here("backup_data", "results", "RNAseq", ds, "deseq2_object.rds")
  if (file.exists(active_path)) active_path
  else if (file.exists(backup_path)) backup_path
  else stop("Cannot find deseq2_object.rds for ", ds)
}

res_list <- lapply(ALL_FOUR, function(ds) readRDS(locate_results(ds))$res_df)
names(res_list) <- ALL_FOUR

up_of   <- function(ds) res_list[[ds]]$gene[res_list[[ds]]$significant & res_list[[ds]]$log2FoldChange > 0]
down_of <- function(ds) res_list[[ds]]$gene[res_list[[ds]]$significant & res_list[[ds]]$log2FoldChange < 0]

pairs <- combn(ALL_FOUR, 2, simplify = FALSE)
out_dir <- results_path("RNAseq", "cross_cohort")

comparison <- bind_rows(lapply(pairs, function(p) {
  bg <- intersect(res_list[[p[1]]]$gene, res_list[[p[2]]]$gene)
  res <- degs_list_overlap(up_of(p[1]), down_of(p[1]), up_of(p[2]), down_of(p[2]), bg)
  cbind(data.frame(cohort_a = p[1], cohort_b = p[2]), res$over_result)
}))

write.csv(comparison, file.path(out_dir, "degs_list_overlap_4way_comparison.csv"), row.names = FALSE)
cat("=== DEG overlap, all 6 pairs, degs_list_overlap() (direction-aware hypergeometric + binomial) ===\n")
print(comparison)
cat("\nSaved to", file.path(out_dir, "degs_list_overlap_4way_comparison.csv"), "\n")
