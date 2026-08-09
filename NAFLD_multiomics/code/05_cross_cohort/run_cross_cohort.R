#!/usr/bin/env Rscript
# Stage 5: cross-cohort comparison across every active bulk RNA-seq dataset
# (code/utils/paths.R's RNASEQ_COHORTS - currently 2, previously 4; this
# script reads that list dynamically rather than assuming a fixed count).
#
# Usage: Rscript code/05_cross_cohort/run_cross_cohort.R
# (no --dataset flag - this stage is the one place all active cohorts meet)
#
# A gene found significant in one cohort could be a true NAFLD signal, or it
# could be a cohort-specific artefact (batch effects, a confounded control
# group, a sequencing-depth effect). Multiple independent cohorts, collected
# by different groups with different patients, make that much less likely
# to happen by chance in all of them at once. This stage asks three separate
# questions of the same data:
#   1. Set overlap  - do the same *genes* pass significance in each cohort?
#   2. Direction concordance - among the shared significant genes, do they
#      move the same direction (up in NAFLD in one cohort, also up in the
#      other)?
#   3. Genome-wide correlation - even for genes that are not individually
#      significant, is the overall fold-change pattern similar between
#      cohorts?
# All three can disagree with each other, and that disagreement is
# informative: e.g. low set overlap but high genome-wide LFC correlation
# usually means "same biology, different statistical power", not
# "different biology".

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})
source(here::here("code", "utils", "paths.R"))
source(here::here("code", "utils", "analysis_helpers.R"))

DATASETS <- names(RNASEQ_COHORTS)
CONSENSUS_MARKERS <- c("TREM2", "SPP1", "GPNMB", "CD68", "FABP4", "FASN", "PNPLA3")

cat("=== Cross-cohort comparison:", paste(DATASETS, collapse = ", "), "===\n")
res_list <- lapply(DATASETS, function(ds) {
  readRDS(file.path(results_path("RNAseq", ds), "deseq2_object.rds"))$res_df
})
names(res_list) <- DATASETS

out_dir <- results_path("RNAseq", "cross_cohort")

# --- 1. pairwise set overlap (Jaccard + Fisher's exact) ---------------------
sig_sets <- lapply(res_list, function(df) df$gene[df$significant])
pairs <- combn(DATASETS, 2, simplify = FALSE)

overlap_summary <- bind_rows(lapply(pairs, function(p) {
  universe <- intersect(res_list[[p[1]]]$gene, res_list[[p[2]]]$gene)
  ft <- fisher_overlap_test(sig_sets[[p[1]]], sig_sets[[p[2]]], universe)
  data.frame(
    cohort_a = p[1], cohort_b = p[2],
    n_sig_a = length(sig_sets[[p[1]]]), n_sig_b = length(sig_sets[[p[2]]]),
    overlap = ft$overlap,
    jaccard = jaccard_index(sig_sets[[p[1]]], sig_sets[[p[2]]]),
    fisher_odds_ratio = ft$odds_ratio, fisher_p = ft$p_value
  )
}))
write.csv(overlap_summary, file.path(out_dir, "pairwise_overlap_summary.csv"), row.names = FALSE)
cat("\nPairwise overlap:\n")
print(overlap_summary %>% select(cohort_a, cohort_b, n_sig_a, n_sig_b, overlap, jaccard, fisher_p))

# --- 2. direction concordance among shared significant genes ----------------
concordance_summary <- bind_rows(lapply(pairs, function(p) {
  shared <- intersect(sig_sets[[p[1]]], sig_sets[[p[2]]])
  if (length(shared) == 0) return(NULL)
  lfc_a <- res_list[[p[1]]]$log2FoldChange[match(shared, res_list[[p[1]]]$gene)]
  lfc_b <- res_list[[p[2]]]$log2FoldChange[match(shared, res_list[[p[2]]]$gene)]
  data.frame(
    cohort_a = p[1], cohort_b = p[2], n_shared_sig = length(shared),
    pct_concordant_direction = round(100 * mean(sign(lfc_a) == sign(lfc_b)), 1)
  )
}))
write.csv(concordance_summary, file.path(out_dir, "direction_concordance.csv"), row.names = FALSE)
cat("\nDirection concordance among shared significant genes:\n")
print(concordance_summary)

# --- 3. consensus gene set (significant in a strong majority of cohorts) ----
# Threshold scales with how many cohorts are active (originally ">=3 of 4",
# i.e. 75%; kept at the same 75% rate rather than hardcoded to "3" so this
# still means the same thing if a cohort is later restored from
# backup_data/ or a new one added). For the current 2 cohorts this reduces
# to "significant in both" - with only 2 cohorts there is no intermediate
# "majority but not all" case, so this is equivalent to the overlap count
# already in pairwise_overlap_summary.csv, not a distinct additional
# result - it is kept here mainly so the file/column names stay stable as
# the cohort count changes.
n_cohorts <- length(DATASETS)
consensus_threshold <- ceiling(n_cohorts * 0.75)
gene_hit_count <- table(unlist(sig_sets))
consensus_genes <- names(gene_hit_count[gene_hit_count >= consensus_threshold])
write.csv(data.frame(gene = consensus_genes, n_cohorts = as.integer(gene_hit_count[consensus_genes])),
          file.path(out_dir, sprintf("consensus_genes_%dof%d.csv", consensus_threshold, n_cohorts)),
          row.names = FALSE)
cat(sprintf("\nConsensus genes (significant in >=%d/%d cohorts): %d\n",
            consensus_threshold, n_cohorts, length(consensus_genes)))

# --- 4. genome-wide LFC correlation across all shared genes ------------------
shared_genes <- Reduce(intersect, lapply(res_list, function(df) df$gene))
lfc_matrix <- sapply(res_list, function(df) df$log2FoldChange[match(shared_genes, df$gene)])
rownames(lfc_matrix) <- shared_genes
cor_matrix <- cor(lfc_matrix, method = "spearman", use = "complete.obs")
write.csv(cor_matrix, file.path(out_dir, "genome_wide_lfc_correlation.csv"))
cat(sprintf("\nGenome-wide log2FC Spearman correlation (%d shared genes):\n", length(shared_genes)))
print(round(cor_matrix, 3))

# --- 5. core consensus marker tracking ---------------------------------------
marker_table <- bind_rows(lapply(DATASETS, function(ds) {
  df <- res_list[[ds]] %>% filter(gene %in% CONSENSUS_MARKERS)
  data.frame(cohort = ds, gene = df$gene, log2FC = round(df$log2FoldChange, 2),
             padj = signif(df$padj, 3), significant = df$significant)
}))
write.csv(marker_table, file.path(out_dir, "consensus_marker_tracking.csv"), row.names = FALSE)

# --- figure: Jaccard heatmap --------------------------------------------------
jac_matrix <- matrix(1, length(DATASETS), length(DATASETS), dimnames = list(DATASETS, DATASETS))
for (p in pairs) {
  j <- jaccard_index(sig_sets[[p[1]]], sig_sets[[p[2]]])
  jac_matrix[p[1], p[2]] <- j
  jac_matrix[p[2], p[1]] <- j
}
jac_long <- as.data.frame(as.table(jac_matrix))
p_heat <- ggplot(jac_long, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = round(Freq, 2)), colour = "white") +
  scale_fill_viridis_c(name = "Jaccard") +
  labs(title = "Pairwise DEG set overlap (Jaccard index)", x = NULL, y = NULL) +
  theme_project()
ggsave(file.path(out_dir, "jaccard_heatmap.png"), p_heat, width = 6, height = 5, dpi = 150)

# --- 6. meta-analysis: combine all active cohorts into one consensus signature
# "Significant in a majority of cohorts" (section 3) is a useful descriptive
# summary, but it throws away effect size and treats a small cohort the same
# as a large one. A proper fixed-effect meta-analysis (Stouffer's method,
# weighted by sqrt(sample size)) combines every active cohort's evidence for
# each gene into one signature, gives more powerful cohorts more say, and
# produces a single combined p-value per gene rather than a vote count. This
# combined signature is what downstream stages (enrichment, reporter
# metabolites, drug repositioning) are run on when treating the active
# cohorts as one consensus finding, via `--dataset cross_cohort`.
cat(sprintf("\n=== Meta-analysis: Stouffer-combined signature across %d active cohort(s) ===\n", n_cohorts))

cohort_n <- sapply(DATASETS, function(ds) {
  length(readRDS(file.path(results_path("RNAseq", ds), "preprocessed.rds"))$condition)
})
weights <- sqrt(cohort_n)

signed_z_matrix <- sapply(DATASETS, function(ds) {
  df <- res_list[[ds]][match(shared_genes, res_list[[ds]]$gene), ]
  sign(df$log2FoldChange) * qnorm(df$pvalue / 2, lower.tail = FALSE)
})
rownames(signed_z_matrix) <- shared_genes

combined_z <- as.numeric(signed_z_matrix %*% weights) / sqrt(sum(weights^2))
combined_lfc <- as.numeric(lfc_matrix %*% weights) / sum(weights)

meta_res_df <- data.frame(
  gene = shared_genes,
  log2FoldChange = combined_lfc,
  pvalue = 2 * pnorm(-abs(combined_z)),
  stringsAsFactors = FALSE
) %>%
  mutate(padj = p.adjust(pvalue, method = "BH")) %>%
  flag_significant() %>%
  arrange(padj)

meta_out_dir <- results_path("RNAseq", "cross_cohort")
write.csv(meta_res_df, file.path(meta_out_dir, "meta_analysis_results.csv"), row.names = FALSE)
# Saved in the same shape as every other cohort's deseq2_object.rds (a
# res_df with gene/log2FoldChange/pvalue/padj/significant) so that
# 03_enrichment, 06_drug_repositioning and 07_metabolic_modelling can all be
# run unmodified with --dataset cross_cohort. There is no vst_matrix here
# (no single expression matrix spans samples from 4 different cohorts
# without batch correction), which is why 04_coexpression is the one stage
# that cannot be run on this pseudo-dataset - WGCNA needs real per-sample
# expression, not a meta-analysis summary.
saveRDS(list(res_df = meta_res_df, dataset_id = "cross_cohort"),
        file.path(meta_out_dir, "deseq2_object.rds"))

cat(sprintf("Meta-analysis: %d / %d genes significant (padj<%.2f, |LFC|>%d): %d up, %d down\n",
            sum(meta_res_df$significant), nrow(meta_res_df), SIG_PADJ, SIG_LFC,
            sum(meta_res_df$significant & meta_res_df$log2FoldChange > 0),
            sum(meta_res_df$significant & meta_res_df$log2FoldChange < 0)))

cat("\nSaved cross-cohort comparison outputs to", out_dir, "\n")
