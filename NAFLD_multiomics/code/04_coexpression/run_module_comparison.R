#!/usr/bin/env Rscript
# Compare WGCNA modules between two cohorts via pairwise Jaccard index -
# "module preservation": does a co-expression module found in one cohort
# correspond to a similar gene set in another, independently-clustered
# cohort?
#
# Usage: Rscript code/04_coexpression/run_module_comparison.R --a GSE162694 --b GSE135251
# (defaults to the two active cohorts in code/utils/paths.R's
# RNASEQ_COHORTS - see that file for the current pairing and why)
#
# Note this is a different question from stage 05's DEG-set overlap: two
# modules can be well-preserved (similar gene membership) even if most of
# their genes are not individually DE-significant, and a module can fail to
# be preserved even if its most-significant DEGs are shared, if the rest of
# the module's membership differs between cohorts. Module preservation asks
# about co-regulation structure, not individual-gene significance.

suppressPackageStartupMessages(library(dplyr))
source(here::here("code", "utils", "paths.R"))
source(here::here("code", "utils", "analysis_helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
# Defaults read from RNASEQ_COHORTS rather than hardcoded here a second
# time, so changing the active pair in code/utils/paths.R is the only edit
# needed - this script cannot silently fall back to an archived cohort.
active_cohorts <- names(RNASEQ_COHORTS)
ds_a <- if ("--a" %in% args) args[which(args == "--a") + 1] else active_cohorts[1]
ds_b <- if ("--b" %in% args) args[which(args == "--b") + 1] else active_cohorts[2]

cat(sprintf("=== WGCNA module comparison: %s vs %s ===\n", ds_a, ds_b))
mod_a <- read.csv(file.path(results_path("RNAseq", ds_a), "wgcna_gene_modules.csv")) |> filter(module != "grey")
mod_b <- read.csv(file.path(results_path("RNAseq", ds_b), "wgcna_gene_modules.csv")) |> filter(module != "grey")

modules_a <- split(mod_a$gene, mod_a$module)
modules_b <- split(mod_b$gene, mod_b$module)

jac_table <- expand.grid(module_a = names(modules_a), module_b = names(modules_b), stringsAsFactors = FALSE) |>
  rowwise() |>
  mutate(
    size_a  = length(modules_a[[module_a]]),
    size_b  = length(modules_b[[module_b]]),
    overlap = length(intersect(modules_a[[module_a]], modules_b[[module_b]])),
    jaccard = jaccard_index(modules_a[[module_a]], modules_b[[module_b]])
  ) |>
  ungroup() |>
  arrange(desc(jaccard))

out_dir <- results_path("RNAseq", "cross_cohort")
write.csv(jac_table, file.path(out_dir, sprintf("module_jaccard_%s_vs_%s.csv", ds_a, ds_b)), row.names = FALSE)

# Best-matching partner for each module_a - the module preservation summary
best_match <- jac_table |> group_by(module_a) |> slice_max(jaccard, n = 1, with_ties = FALSE) |> ungroup()
write.csv(best_match, file.path(out_dir, sprintf("module_best_match_%s_vs_%s.csv", ds_a, ds_b)), row.names = FALSE)

cat(sprintf("%d modules (%s) x %d modules (%s) compared\n", length(modules_a), ds_a, length(modules_b), ds_b))
cat("\nBest-preserved module pairs (by Jaccard):\n")
print(best_match |> arrange(desc(jaccard)) |> select(module_a, module_b, size_a, size_b, overlap, jaccard))
