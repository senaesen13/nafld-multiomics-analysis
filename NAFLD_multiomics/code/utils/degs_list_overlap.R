# DEG list overlap test - extracted from
# /Users/k2254978/Desktop/Work/01_Projects/Microarray_Urbani_Clean/source_code/Codes/DEGs_list_overlap.R
# (original algorithm and statistics preserved exactly; adapted from
# matrix/dim()-based counting to plain vectors/length(), and documented,
# but the actual math is unchanged).
#
# This is a different, more direction-aware test than
# code/utils/analysis_helpers.R's hypergeometric_overlap_test(): that
# function tests whether two DEG lists overlap more than chance,
# regardless of direction, then (as a separate check) asks whether the
# overlapping genes happen to agree on direction. This function instead
# folds direction into the hypergeometric test itself - the "success"
# count is the number of genes that are UP in both datasets or DOWN in
# both datasets specifically (found via separate intersections of the
# up-only and down-only lists, not by checking sign agreement after the
# fact), which is a stricter, more biologically specific test: two DEG
# lists can share many genes while disagreeing about direction, and this
# test would not count that as enrichment, whereas a direction-blind
# overlap test would.
#
# One quirk preserved faithfully from the original: `over_all` (used only
# to report the "ratio" column) is computed from the *un-background-
# restricted* combined up+down lists, while `over_up`/`over_down` (which
# feed both p-values) are computed from lists already intersected with
# bg_gene. In the original code this looks like an artefact of variable
# reuse rather than a deliberate design choice, but it is kept as-is here
# rather than "corrected" - this is a faithful extraction, not a rewrite.

#' @param deg_up_1,deg_down_1 up- and down-regulated significant gene IDs, dataset 1
#' @param deg_up_2,deg_down_2 up- and down-regulated significant gene IDs, dataset 2
#' @param bg_gene background/universe gene IDs (e.g. genes tested in both datasets)
degs_list_overlap <- function(deg_up_1, deg_down_1, deg_up_2, deg_down_2, bg_gene) {
  bg <- length(bg_gene)

  deg_all_1 <- c(deg_up_1, deg_down_1)
  deg_all_2 <- c(deg_up_2, deg_down_2)

  deg_up_1_bg   <- intersect(deg_up_1, bg_gene)
  deg_down_1_bg <- intersect(deg_down_1, bg_gene)
  deg_up_2_bg   <- intersect(deg_up_2, bg_gene)
  deg_down_2_bg <- intersect(deg_down_2, bg_gene)

  num_deg_11 <- length(deg_up_1_bg) + length(deg_down_1_bg)
  num_deg_22 <- length(deg_up_2_bg) + length(deg_down_2_bg)

  deg_inform <- rbind(
    deg_1 = c(all = num_deg_11, up = length(deg_up_1_bg), down = length(deg_down_1_bg)),
    deg_2 = c(all = num_deg_22, up = length(deg_up_2_bg), down = length(deg_down_2_bg))
  )

  over_all  <- intersect(deg_all_1, deg_all_2)
  over_up   <- intersect(deg_up_1_bg, deg_up_2_bg)
  over_down <- intersect(deg_down_1_bg, deg_down_2_bg)
  n_consistent <- length(over_up) + length(over_down)

  # Binomial: of all genes significant in both datasets regardless of
  # direction (over_all), is the number moving in a *consistent* direction
  # (n_consistent) higher than the 50% expected by chance?
  p_binom <- pbinom(n_consistent - 1, length(over_all), 0.5, lower.tail = FALSE)

  # Hypergeometric: given a universe of `bg` genes containing num_deg_11
  # dataset-1 DEGs, what is the probability of drawing at least
  # n_consistent of them - by direction-consistent overlap specifically,
  # not just any overlap - in a random draw of num_deg_22 genes (dataset
  # 2's DEG count)?
  p_hyper <- phyper(n_consistent - 1, num_deg_11, bg - num_deg_11, num_deg_22, lower.tail = FALSE)

  over_result <- data.frame(
    deg_1 = num_deg_11, deg_2 = num_deg_22,
    overlap = length(over_all), consistent = n_consistent,
    ratio = n_consistent / length(over_all),
    co_up = length(over_up), co_down = length(over_down),
    p_binom = p_binom, p_hyper = p_hyper
  )

  list(deg_inform = deg_inform, over_result = over_result,
       over_up = over_up, over_down = over_down)
}
