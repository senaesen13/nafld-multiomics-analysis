#!/usr/bin/env Rscript
# DEG overlap across all 4 original cohorts (2 active + 2 archived),
# using the hypergeometric + binomial-concordance methodology from
# bioinformatics-learning/improvements/01_deg_overlap_jaccard_ak.R
# (hypergeometric overlap test) and
# bioinformatics-learning/week4-day4-nafld-cross-cohort-comparison/scripts/
# pairwise_comparison.R (binomial direction-concordance test) - not the
# Fisher's-exact-only version code/05_cross_cohort/run_cross_cohort.R uses.
#
# Usage: Rscript code/05_cross_cohort/hypergeometric_binomial_4way.R
#
# Reads all 4 cohorts directly from wherever they currently sit (active
# results/RNAseq/ or archived backup_data/results/RNAseq/) without moving
# any files - this is a read-only comparison, independent of which cohorts
# are "active" in code/utils/paths.R right now.

suppressPackageStartupMessages(library(dplyr))
source(here::here("code", "utils", "paths.R"))
source(here::here("code", "utils", "analysis_helpers.R"))

ALL_FOUR <- c("GSE162694", "GSE135251", "GSE126848", "EMTAB12807")

locate_results <- function(ds) {
  active_path <- here::here("results", "RNAseq", ds, "deseq2_object.rds")
  backup_path <- here::here("backup_data", "results", "RNAseq", ds, "deseq2_object.rds")
  if (file.exists(active_path)) active_path
  else if (file.exists(backup_path)) backup_path
  else stop("Cannot find deseq2_object.rds for ", ds, " in results/ or backup_data/results/")
}

res_list <- lapply(ALL_FOUR, function(ds) readRDS(locate_results(ds))$res_df)
names(res_list) <- ALL_FOUR
sig_sets <- lapply(res_list, function(df) df$gene[df$significant])

pairs <- combn(ALL_FOUR, 2, simplify = FALSE)

comparison <- bind_rows(lapply(pairs, function(p) {
  hyper <- hypergeometric_overlap_test(sig_sets[[p[1]]], sig_sets[[p[2]]])

  shared <- intersect(sig_sets[[p[1]]], sig_sets[[p[2]]])
  lfc_a <- res_list[[p[1]]]$log2FoldChange[match(shared, res_list[[p[1]]]$gene)]
  lfc_b <- res_list[[p[2]]]$log2FoldChange[match(shared, res_list[[p[2]]]$gene)]
  binom <- binomial_concordance_test(lfc_a, lfc_b)

  data.frame(
    cohort_a = p[1], cohort_b = p[2],
    n_sig_a = hyper$n_a, n_sig_b = hyper$n_b,
    overlap = hyper$overlap,
    jaccard = round(jaccard_index(sig_sets[[p[1]]], sig_sets[[p[2]]]), 4),
    hypergeometric_p = signif(hyper$p_value, 4),
    n_concordant = binom$n_concordant,
    pct_concordant = binom$pct_concordant,
    binomial_p = signif(binom$p_value, 4)
  )
}))

out_dir <- results_path("RNAseq", "cross_cohort")
write.csv(comparison, file.path(out_dir, "hypergeometric_binomial_4way_comparison.csv"), row.names = FALSE)

cat("=== DEG overlap, all 6 pairs, hypergeometric + binomial concordance ===\n")
print(comparison)
cat("\nSaved to", file.path(out_dir, "hypergeometric_binomial_4way_comparison.csv"), "\n")
