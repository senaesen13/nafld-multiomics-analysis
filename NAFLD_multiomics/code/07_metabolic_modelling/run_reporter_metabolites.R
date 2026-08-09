#!/usr/bin/env Rscript
# Stage 7: Reporter Metabolite analysis (Patil & Nielsen, 2005, PNAS).
#
# Usage: Rscript code/07_metabolic_modelling/run_reporter_metabolites.R --dataset GSE162694
#
# Differential expression tells you which genes changed. It does not, by
# itself, say which parts of metabolism are affected, because one metabolite
# is typically produced or consumed by several enzymes (genes), and a
# metabolite can be under coordinated regulatory pressure even if no single
# one of its neighbouring genes individually reaches significance. Reporter
# Metabolite analysis re-expresses gene-level statistical significance as a
# property of the metabolic network's nodes (metabolites) rather than genes:
#
#   1. Convert each gene's DE p-value to a Z-score (bigger = more
#      significant, regardless of direction).
#   2. For each metabolite, average the Z-scores of every gene encoding an
#      enzyme that produces or consumes it (its "neighbourhood" in the
#      genome-scale metabolic network).
#   3. That raw average is biased by how many neighbours a metabolite has
#      (hubs with many neighbours regress toward the mean p-value just by
#      the law of large numbers) - correct for this by comparing each
#      metabolite's score against 1,000 random gene sets of the same size,
#      drawn from the same network, and standardising against that null.
#
# The genome-scale network itself (which genes' products act on which
# metabolites) comes from Human-GEM (Robinson et al. 2020, Sci Signal),
# not from this project's own expression data - it encodes literature-
# curated human metabolic stoichiometry, independent of any cohort.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
source(here::here("code", "utils", "paths.R"))
source(here::here("code", "utils", "analysis_helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
dataset_id <- args[which(args == "--dataset") + 1]
if (is.na(dataset_id) || dataset_id == "") stop("Usage: --dataset <dataset_id>")

N_PERMUTATIONS <- 1000
MIN_NEIGHBOURS <- 3    # metabolites with too few neighbours give an unstable Z-score
MAX_NEIGHBOURS <- 100  # metabolites with too many (currency-like hubs) dilute any real signal

# Metabolites that participate in a very large fraction of all reactions
# (energy/redox/proton carriers) act as network hubs rather than pathway-
# specific nodes, and would dominate the "top hits" list without this
# exclusion despite not being biologically informative on their own.
CURRENCY_METABOLITES <- c("H2O", "ATP", "ADP", "AMP", "NAD+", "NADH", "NADP+", "NADPH",
                           "H+", "Pi", "PPi", "CoA", "CO2", "O2", "HCO3-", "Na+", "K+",
                           "Cl-", "water", "oxygen", "phosphate")

cat(sprintf("=== Reporter metabolites: %s ===\n", dataset_id))
out_dir <- results_path("RNAseq", dataset_id)
res_df  <- readRDS(file.path(out_dir, "deseq2_object.rds"))$res_df

gem_map <- read.csv(dataset_path("GEM", "human_gem_topology_network.csv")) |>
  filter(!is.na(metabolite_name), metabolite_name != "",
         !is.na(gene_symbol), gene_symbol != "",
         !tolower(metabolite_name) %in% tolower(CURRENCY_METABOLITES))
cat(sprintf("GEM network: %d metabolite-gene edges, %d unique metabolites\n",
            nrow(gem_map), length(unique(gem_map$metabolite_name))))

# --- gene-level Z-scores -----------------------------------------------------
# qnorm(1 - p/2) looks equivalent to qnorm(p/2, lower.tail=FALSE) but is not,
# numerically: doubles only carry ~16 significant digits, so once p/2 drops
# below ~1e-16, "1 - p/2" rounds to exactly 1.0 and qnorm() returns Inf. A
# well-powered cohort easily produces p-values that small, and one Inf gene
# Z-score is enough to poison every permutation that happens to draw it
# (sd() of a set containing Inf is NA). lower.tail=FALSE computes the upper
# tail directly, without subtracting from 1, so it stays finite.
gene_stats <- res_df |>
  filter(!is.na(pvalue)) |>
  mutate(pvalue = pmax(pvalue, .Machine$double.xmin),
         gene_z = qnorm(pvalue / 2, lower.tail = FALSE))
gene_z_vec    <- setNames(gene_stats$gene_z, gene_stats$gene)
gene_lfc_vec  <- setNames(gene_stats$log2FoldChange, gene_stats$gene)
gene_universe <- names(gene_z_vec)

# --- metabolite neighbourhoods within the size window ------------------------
metabolite_genes <- gem_map |>
  filter(gene_symbol %in% gene_universe) |>
  group_by(metabolite_name) |>
  summarise(genes = list(unique(gene_symbol)), .groups = "drop") |>
  mutate(k = lengths(genes)) |>
  filter(k >= MIN_NEIGHBOURS, k <= MAX_NEIGHBOURS)
cat(sprintf("Testing %d metabolites with %d-%d gene neighbours\n",
            nrow(metabolite_genes), MIN_NEIGHBOURS, MAX_NEIGHBOURS))

# --- Patil & Nielsen Z-score with permutation correction ---------------------
# The null distribution of "average Z-score of k randomly chosen genes" only
# depends on k, not on which metabolite is being tested - many metabolites
# share the same neighbourhood size, so the background only needs to be
# simulated once per *unique* k (a few dozen), not once per metabolite (a
# few thousand). This is the standard way to implement Patil & Nielsen
# efficiently, not just a speed hack: it gives the identical null
# distribution a metabolite-by-metabolite simulation would, for a fraction
# of the computation.
set.seed(42)
unique_ks <- sort(unique(metabolite_genes$k))
background_by_k <- setNames(
  lapply(unique_ks, function(k) {
    replicate(N_PERMUTATIONS, sum(gene_z_vec[sample(gene_universe, k)]) / sqrt(k))
  }),
  as.character(unique_ks)
)

reporter_results <- lapply(seq_len(nrow(metabolite_genes)), function(i) {
  g_set <- metabolite_genes$genes[[i]]
  k     <- metabolite_genes$k[i]

  z_raw  <- sum(gene_z_vec[g_set], na.rm = TRUE) / sqrt(k)
  perm_z <- background_by_k[[as.character(k)]]
  z_corr <- if (sd(perm_z) > 0) (z_raw - mean(perm_z)) / sd(perm_z) else z_raw

  data.frame(
    metabolite = metabolite_genes$metabolite_name[i],
    n_neighbours = k,
    mean_log2FC = round(mean(gene_lfc_vec[g_set], na.rm = TRUE), 3),
    reporter_z = round(z_corr, 3),
    pvalue = signif(pnorm(z_corr, lower.tail = FALSE), 4)
  )
}) |>
  bind_rows() |>
  mutate(padj = signif(p.adjust(pvalue, method = "BH"), 4)) |>
  arrange(desc(reporter_z))

write.csv(reporter_results, file.path(out_dir, "reporter_metabolites.csv"), row.names = FALSE)
# With thousands of metabolites tested, BH correction is often conservative
# enough that few or none pass padj<0.05 - the field convention (following
# Patil & Nielsen's own papers) is to interpret this method by rank, i.e.
# the top 10-20 metabolites in reporter_metabolites.csv, not a strict cutoff.
cat(sprintf("%d / %d metabolites pass padj<0.05; top hit: %s (Z=%.2f, p=%.1e)\n",
            sum(reporter_results$padj < 0.05), nrow(reporter_results),
            reporter_results$metabolite[1], reporter_results$reporter_z[1],
            reporter_results$pvalue[1]))

# --- figure: top 20 reporter metabolites --------------------------------------
top20 <- head(reporter_results, 20)
p <- ggplot(top20, aes(x = reorder(metabolite, reporter_z), y = reporter_z, fill = mean_log2FC)) +
  geom_col(colour = "black", width = 0.7) +
  coord_flip() +
  scale_fill_gradient2(low = "#1E88E5", mid = "grey90", high = "#D81B60", midpoint = 0,
                        name = "Mean log2FC") +
  labs(title = paste("Top 20 reporter metabolites -", dataset_id),
       x = NULL, y = "Corrected reporter Z-score") +
  theme_project()
ggsave(file.path(out_dir, "reporter_metabolites_top20.png"), p, width = 8, height = 6.5, dpi = 150)

cat("Saved reporter metabolite results to", out_dir, "\n")
