#!/usr/bin/env Rscript
# Metabolomics: Ji et al. 2022 plasma metabolomics, re-grouped to the same
# Control-vs-NAFLD convention used throughout this project.
#
# Usage: Rscript code/10_metabolomics/run_ji2022_metabolomics.R
#
# This dataset is deposited as *published group-level summary statistics*
# (mean and SD per metabolite per group: Control n=25, NAFL n=42, NASH
# n=19), not raw per-patient measurements - the paper's own p-values are
# for the pairwise Control-vs-NAFL and Control-vs-NASH contrasts only.
# Since this project treats NAFLD as the disease umbrella (NAFL + NASH
# merged, matching every RNA-seq cohort), the NAFL and NASH groups first
# need to be combined into one "NAFLD" group before testing - and because
# only summary statistics are available, that combination has to be done
# correctly:
#   - the combined mean is the sample-size-weighted mean of the two groups
#   - the combined SD is NOT the pooled within-group SD alone; it must also
#     account for the fact that NAFL and NASH have different means, or the
#     spread of the merged group will be systematically underestimated.
#     This is the standard formula for combining two summary-statistic
#     subgroups (see e.g. Cochrane Handbook 6.5.2.10).
# A new Control-vs-NAFLD significance test is then computed from the
# combined summary statistics via Welch's t-test (the standard test for two
# groups with unequal variance, computed here directly from means/SDs/Ns
# since the raw values are not available), followed by BH-FDR correction
# across all 79 metabolites for this specific (re-derived) contrast - the
# paper's own Qval_FDR column corrected for the original 3-way comparisons
# and does not apply to this regrouped one.

suppressPackageStartupMessages(library(dplyr))
source(here::here("code", "utils", "paths.R"))

N_CONTROL <- 25
N_NAFL    <- 42
N_NASH    <- 19

df <- read.csv(dataset_path("metabolomics", "Ji2022", "Table_S1_Metabolite_Levels.csv"))
out_dir <- results_path("metabolomics", "Ji2022")

# --- combine NAFL + NASH into one NAFLD group --------------------------------
pooled_mean <- (N_NAFL * df$NAFL_Mean + N_NASH * df$NASH_Mean) / (N_NAFL + N_NASH)
pooled_var <- (
  (N_NAFL - 1) * df$NAFL_SD^2 + (N_NASH - 1) * df$NASH_SD^2 +
    (N_NAFL * N_NASH / (N_NAFL + N_NASH)) * (df$NAFL_Mean - df$NASH_Mean)^2
) / (N_NAFL + N_NASH - 1)
df$NAFLD_Mean <- pooled_mean
df$NAFLD_SD   <- sqrt(pooled_var)
N_NAFLD <- N_NAFL + N_NASH

# --- Welch's t-test, Control vs combined NAFLD, from summary statistics -----
se_diff <- sqrt(df$Control_SD^2 / N_CONTROL + df$NAFLD_SD^2 / N_NAFLD)
t_stat  <- (df$NAFLD_Mean - df$Control_Mean) / se_diff
welch_df <- se_diff^4 / (
  (df$Control_SD^2 / N_CONTROL)^2 / (N_CONTROL - 1) +
    (df$NAFLD_SD^2 / N_NAFLD)^2 / (N_NAFLD - 1)
)
df$pvalue <- 2 * pt(-abs(t_stat), df = welch_df)
df$padj   <- p.adjust(df$pvalue, method = "BH")
df$log2FC_NAFLD_vs_Control <- log2(df$NAFLD_Mean / df$Control_Mean)
# No |log2FC| cutoff is applied (unlike the gene-level padj<0.05 & |LFC|>1
# convention used elsewhere in this project) - metabolite fold-changes are
# not directly comparable to gene expression fold-changes in magnitude, and
# an effect-size threshold calibrated for one is not automatically
# meaningful for the other. Significance here is padj<0.05 alone, consistent
# with how the original Ji et al. paper itself reported significance.
df$significant <- !is.na(df$padj) & df$padj < 0.05

result <- df |>
  select(Category, Metabolite, Control_Mean, Control_SD, NAFLD_Mean, NAFLD_SD,
         log2FC_NAFLD_vs_Control, pvalue, padj, significant) |>
  arrange(padj)

write.csv(result, file.path(out_dir, "ji2022_control_vs_nafld.csv"), row.names = FALSE)
cat(sprintf("Ji2022 (re-grouped Control n=%d vs NAFLD n=%d): %d / %d metabolites significant (padj<0.05)\n",
            N_CONTROL, N_NAFLD, sum(result$significant), nrow(result)))
print(head(result[result$significant, c("Metabolite", "log2FC_NAFLD_vs_Control", "padj")], 10))
