#!/usr/bin/env Rscript
# Metabolomics: MTBLS174 NMR serum metabolomics, High vs Low steatosis.
#
# Usage: Rscript code/10_metabolomics/run_mtbls174_metabolomics.R
#
# Unlike Ji2022 (stage 10's other script), every sample in this cohort is
# already NAFLD - there is no healthy control arm at all (MTBLS174's own
# sample sheet labels every subject's "Metabolic syndrome" factor as
# "Human non-alcoholic fatty liver disease"). So this cannot be forced into
# the project's Control-vs-NAFLD convention the way Ji2022 was: this stage
# instead reports its own thing - a within-disease severity gradient,
# comparing subjects with high liver fat content (steatosis) against
# subjects with low liver fat content, both still NAFLD. Do not read
# "Low steatosis" as a stand-in for "healthy".
#
# Method matches bioinformatics-learning/week4-day5-nafld-metabolomics/
# scripts/nafld_nmr_metabolomics.R (limma moderated t-test on log2-
# transformed NMR intensities - not DESeq2 or a t-test from summary
# statistics, because this cohort's raw per-sample metabolite values are
# available, unlike Ji2022's published-summary-only data), with one fix:
# that script's steatosis parser intended "<5" (below quantification
# limit) to be treated as 0 (see its own comment), but
# gsub("<", "0", "<5") actually produces the string "05" -> 5, not 0. Kept
# as an actual 0 here. In this specific dataset only one sample is
# affected and it does not change that sample's Low/Medium/High bucket
# (0 and 5 are both "Low"), but the exact steatosis_pct value reported for
# it does change.

suppressPackageStartupMessages({
  library(dplyr)
  library(limma)
  library(ggplot2)
})
source(here::here("code", "utils", "paths.R"))

data_dir <- dataset_path("metabolomics", "MTBLS174")
out_dir  <- results_path("metabolomics", "MTBLS174")

# --- metabolite abundance matrix --------------------------------------------
maf <- read.delim(file.path(data_dir, "m_MTBLS174_hna_fld_metabolite_profiling_NMR_spectroscopy_v2_maf.tsv"),
                   check.names = FALSE, stringsAsFactors = FALSE)
meta_cols  <- 18  # columns 1-18 are metabolite metadata (CHEBI ID, formula, ...); 19+ are per-sample abundances
sample_ids <- colnames(maf)[(meta_cols + 1):ncol(maf)]

abund <- as.matrix(maf[, sample_ids, drop = FALSE])
mode(abund) <- "numeric"
rownames(abund) <- maf$metabolite_identification
abund_log <- log2(abund + 1)

# --- sample metadata + steatosis grading ------------------------------------
meta_raw <- read.delim(file.path(data_dir, "s_MTBLS174.txt"), check.names = FALSE, stringsAsFactors = FALSE)
# The ISA-Tab sample sheet repeats "Term Source REF"/"Term Accession
# Number"/"Unit" as column headers once per Factor Value block, so the raw
# data.frame has several literally-duplicate column names - harmless for
# base R indexing, but dplyr's NSE can't disambiguate them even when they
# aren't among the columns being selected. make.unique() resolves this
# without touching the (already-unique) "Factor Value[...]" names actually
# used below.
colnames(meta_raw) <- make.unique(colnames(meta_raw))
meta <- meta_raw %>%
  transmute(
    sample_name = `Sample Name`,
    gender      = `Factor Value[Gender]`,
    bmi         = `Factor Value[BMI]`,
    age         = `Factor Value[Age]`,
    steatosis   = `Factor Value[Steatosis]`
  ) %>%
  filter(sample_name %in% sample_ids)

# Values are recorded as "<5" (below quantification limit -> 0), ">X", "X%",
# or a range "X-Y" (take the midpoint).
parse_steatosis <- function(x) {
  below_limit <- grepl("^<", x)
  x <- gsub("^<", "", x)
  x <- gsub(">", "", x)
  x <- gsub("%", "", x)
  val <- sapply(x, function(v) mean(as.numeric(strsplit(v, "-")[[1]]), na.rm = TRUE))
  ifelse(below_limit, 0, val)
}
meta$steatosis_pct <- parse_steatosis(meta$steatosis)

# Low <20%, Medium 20-40%, High >40% (same breakpoints as the original
# analysis) - Medium is excluded from the differential test below, kept
# only as a reference group shown on the PCA plot.
meta$steatosis_grade <- cut(meta$steatosis_pct, breaks = c(-Inf, 19.9, 40, Inf),
                             labels = c("Low", "Medium", "High"))
meta <- meta[match(sample_ids, meta$sample_name), ]
rownames(meta) <- meta$sample_name

cat("Samples:", nrow(meta), " Metabolites:", nrow(abund_log), "\n")
cat("Steatosis grades:\n")
print(table(meta$steatosis_grade, useNA = "ifany"))

# --- PCA (all grades, exploratory) ------------------------------------------
pca <- prcomp(t(abund_log), scale. = TRUE)
pca_df <- as.data.frame(pca$x[, 1:2]) %>%
  mutate(sample = rownames(.)) %>%
  left_join(meta, by = c("sample" = "sample_name"))
pct_var <- round(100 * summary(pca)$importance[2, 1:2], 1)

p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = steatosis_grade)) +
  geom_point(size = 3) +
  scale_colour_manual(values = c(Low = "#2196F3", Medium = "#FF9800", High = "#F44336"), na.value = "grey60") +
  labs(title = "PCA - MTBLS174 NMR metabolomics (NAFLD serum)",
       x = paste0("PC1 (", pct_var[1], "%)"), y = paste0("PC2 (", pct_var[2], "%)"),
       colour = "Steatosis grade") +
  theme_bw()
ggsave(file.path(out_dir, "pca_steatosis.png"), p_pca, width = 7, height = 5, dpi = 150)

# --- differential metabolite analysis: High vs Low steatosis ----------------
# limma's moderated t-test (empirical Bayes variance shrinkage across all
# metabolites) rather than a plain per-metabolite t-test: with only a
# double-digit number of metabolites and samples, per-metabolite variance
# estimates are noisy on their own, and moderation borrows strength across
# the whole panel the same way DESeq2's dispersion shrinkage does for genes.
keep <- meta$steatosis_grade %in% c("Low", "High")
design <- model.matrix(~ factor(meta$steatosis_grade[keep], levels = c("Low", "High")))
fit <- eBayes(lmFit(abund_log[, keep], design))

result <- topTable(fit, coef = 2, number = Inf, sort.by = "P") %>%
  mutate(metabolite = rownames(.), significant = adj.P.Val < 0.05 & abs(logFC) > 0.5) %>%
  select(metabolite, logFC, AveExpr, P.Value, adj.P.Val, significant)

write.csv(result, file.path(out_dir, "mtbls174_high_vs_low_steatosis.csv"), row.names = FALSE)
cat(sprintf("High vs Low steatosis (n=%d Low, %d High): %d / %d metabolites significant (padj<0.05, |logFC|>0.5)\n",
            sum(meta$steatosis_grade[keep] == "Low"), sum(meta$steatosis_grade[keep] == "High"),
            sum(result$significant), nrow(result)))
print(head(result[result$significant, ], 10))

cat("Saved MTBLS174 results to", out_dir, "\n")
