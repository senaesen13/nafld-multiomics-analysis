# Shared analysis functions used by more than one stage of the pipeline.
#
# Keeping these in one place means every script that needs "is this gene
# significant?" or "how do I rank genes for GSEA?" uses the exact same
# definition. That matters: if 02_differential_expression and
# 05_cross_cohort used slightly different significance thresholds, the
# gene counts reported in each stage would silently disagree.

suppressPackageStartupMessages({
  library(dplyr)
})

# ---------------------------------------------------------------------------
# Significance calling
# ---------------------------------------------------------------------------

#' Standard significance threshold used everywhere in this project:
#' adjusted p < 0.05 AND at least a 2-fold change (|log2FC| > 1).
#' Using the MLE (unshrunken) log2FC for the *filter* and the apeglm-shrunk
#' value for *reporting effect size* is deliberate: shrinkage pulls
#' low-count genes toward zero, which is desirable when quoting an effect
#' size but would make the significance filter itself count-dependent.
SIG_PADJ <- 0.05
SIG_LFC  <- 1

flag_significant <- function(df, padj_col = "padj", lfc_col = "log2FoldChange") {
  df %>%
    mutate(significant = !is.na(.data[[padj_col]]) &
             .data[[padj_col]] < SIG_PADJ &
             abs(.data[[lfc_col]]) > SIG_LFC)
}

# ---------------------------------------------------------------------------
# GSEA ranking
# ---------------------------------------------------------------------------

#' Standardized GSEA ranking metric: sign(log2FC) * -log10(p-value).
#'
#' Why not rank by log2FC alone? A gene can have a huge fold change purely
#' because it is barely expressed in one group (noisy estimate). Why not
#' rank by p-value alone? Then direction of change is lost. Combining both
#' - sign from the fold change, magnitude from statistical confidence -
#' gives a ranking that reflects "how confidently, and in which direction,
#' does this gene move", which is what GSEA's running-sum statistic needs.
#'
#' Uses the *unshrunken* (MLE) log2FC for the sign, and the nominal
#' (not BH-adjusted) p-value for the magnitude, following the standard
#' GSEA Preranked convention.
gsea_rank_metric <- function(log2fc, pvalue) {
  sign(log2fc) * -log10(pmax(pvalue, .Machine$double.xmin))
}

#' Build a ranked gene vector ready for clusterProfiler::GSEA() or export
#' to a GSEA Desktop .rnk file.
build_gsea_ranking <- function(df, gene_col, lfc_col = "log2FoldChange", pval_col = "pvalue") {
  df <- df[!is.na(df[[lfc_col]]) & !is.na(df[[pval_col]]) & !is.na(df[[gene_col]]), ]
  ranks <- gsea_rank_metric(df[[lfc_col]], df[[pval_col]])
  names(ranks) <- df[[gene_col]]
  sort(ranks, decreasing = TRUE)
}

#' Write a ranked vector to a GSEA Desktop-compatible .rnk file
write_rnk <- function(ranked_vec, path) {
  out <- data.frame(gene = names(ranked_vec), score = as.numeric(ranked_vec))
  write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# ---------------------------------------------------------------------------
# Set overlap (used by cross-cohort comparison)
# ---------------------------------------------------------------------------

jaccard_index <- function(set_a, set_b) {
  length(intersect(set_a, set_b)) / length(union(set_a, set_b))
}

#' Fisher's exact test for whether two gene sets overlap more than expected
#' by chance, given a shared background universe (e.g. all genes tested in
#' both cohorts after expression filtering).
fisher_overlap_test <- function(set_a, set_b, universe) {
  set_a <- intersect(set_a, universe)
  set_b <- intersect(set_b, universe)
  both     <- length(intersect(set_a, set_b))
  only_a   <- length(setdiff(set_a, set_b))
  only_b   <- length(setdiff(set_b, set_a))
  neither  <- length(universe) - both - only_a - only_b
  tab <- matrix(c(both, only_a, only_b, neither), nrow = 2)
  ft <- fisher.test(tab)
  list(p_value = ft$p.value, odds_ratio = unname(ft$estimate), overlap = both)
}

#' Hypergeometric test for whether two gene sets overlap more than expected
#' by chance - the one-tailed test underlying the classic "gene set
#' enrichment against a background" calculation, computed directly with
#' phyper() rather than via fisher.test(). For a single pairwise 2x2 table
#' the two agree to many decimal places (a one-tailed Fisher's exact test
#' on a 2x2 table *is* a hypergeometric test), but phyper() is the direct,
#' standard formulation: "given a universe of universe_n genes containing
#' n_b genes from set B, what is the probability of drawing at least k of
#' them by chance in a random draw of n_a genes (the size of set A)?"
#'
#' Matches bioinformatics-learning/improvements/01_deg_overlap_jaccard_ak.R's
#' calc_overlap_stats() exactly: set_a/set_b sizes are used as-is (NOT
#' restricted to genes shared with the other cohort - each cohort's own
#' full significant-gene count is what n_a/n_b report, so they stay
#' consistent across every pairing that cohort appears in), and universe_n
#' is a fixed approximation of the protein-coding genome (default 20,000)
#' rather than a per-pair "genes tested in both" count.
hypergeometric_overlap_test <- function(set_a, set_b, universe_n = 20000) {
  n_a <- length(set_a)
  n_b <- length(set_b)
  k   <- length(intersect(set_a, set_b))
  p_value <- phyper(k - 1, n_b, universe_n - n_b, n_a, lower.tail = FALSE)
  list(p_value = p_value, overlap = k, n_a = n_a, n_b = n_b)
}

#' Binomial test for direction concordance: of the genes significant in
#' both cohorts, is the fraction moving in the *same* direction (both up or
#' both down) significantly greater than the 50% expected by chance alone?
#' This is a different question from the overlap test above - two cohorts
#' can share many significant genes (a strong hypergeometric hit) while
#' disagreeing about which way those genes move, which this test would
#' catch and the overlap test alone would not.
binomial_concordance_test <- function(lfc_a, lfc_b) {
  n_total <- length(lfc_a)
  n_concordant <- sum(sign(lfc_a) == sign(lfc_b))
  bt <- binom.test(n_concordant, n_total, p = 0.5, alternative = "greater")
  list(p_value = bt$p.value, n_concordant = n_concordant, n_total = n_total,
       pct_concordant = round(100 * n_concordant / n_total, 1))
}

# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

#' Consistent ggplot theme for all figures produced in this project.
theme_project <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right"
    )
}
