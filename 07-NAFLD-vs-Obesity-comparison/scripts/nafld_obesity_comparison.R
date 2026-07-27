## ============================================================
## NAFLD vs Obesity gene signature comparison
## NAFLD primary signature: 139-gene overlap (GSE162694 ∩ GSE135251)
## NAFLD secondary: D1 alone (485 sig genes) and D2 alone (1058 sig genes)
## Obesity signature: AZGP1, CLEC4C (GSE166047 adipose, padj < 0.05 & |LFC| > 1)
## Metrics: overlap count, Jaccard Index, Fisher's exact, direction concordance
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

if (interactive()) {
  tryCatch({
    proj_dir <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
    setwd(proj_dir)
  }, error = function(e) invisible(NULL))
}
cat("Working directory:", getwd(), "\n")
dir.create("results", showWarnings = FALSE)
dir.create("plots",   showWarnings = FALSE)

# ============================================================
# Paths (relative to 07-NAFLD-vs-Obesity-comparison/)
# ============================================================
d1_dir  <- "../01-NAFLD-discovery-cohort-GSE162694"
d2_dir  <- "../02-NAFLD-validation-cohort-GSE135251"
ob_dir  <- "data"   # obesity files copied locally (GSE166047)
nafld_overlap_file <- "../06-NAFLD-coexpression-analysis/results/module_assignments.csv"

# ============================================================
# 1. Load gene lists
# ============================================================
cat("\n=== 1. Loading gene lists ===\n")

ob_sig <- read.csv(file.path(ob_dir, "deseq2_significant_genes.csv"),
                   stringsAsFactors = FALSE)
ob_all <- read.csv(file.path(ob_dir, "deseq2_results_all.csv"),
                   stringsAsFactors = FALSE)

d1_all <- read.csv(file.path(d1_dir, "results/deseq2_results.csv"),
                   stringsAsFactors = FALSE)
d2_all <- read.csv(file.path(d2_dir, "results/gse135251_results.csv"),
                   stringsAsFactors = FALSE)

d1_sig <- read.csv(file.path(d1_dir, "results/significant_genes.csv"),
                   stringsAsFactors = FALSE)
d2_sig <- read.csv(file.path(d2_dir, "results/significant_genes.csv"),
                   stringsAsFactors = FALSE)

nafld_overlap <- read.csv(nafld_overlap_file, stringsAsFactors = FALSE)

cat(sprintf("Obesity sig genes:        %d (%s)\n",
            nrow(ob_sig), paste(ob_sig$gene_name, collapse = ", ")))
cat(sprintf("NAFLD D1 sig genes:       %d\n", nrow(d1_sig)))
cat(sprintf("NAFLD D2 sig genes:       %d\n", nrow(d2_sig)))
cat(sprintf("NAFLD 139-gene overlap:   %d\n", nrow(nafld_overlap)))

ob_syms    <- na.omit(unique(ob_sig$gene_name[ob_sig$gene_name != ""]))
d1_syms    <- na.omit(unique(d1_sig$gene_symbol[d1_sig$gene_symbol != ""]))
d2_syms    <- na.omit(unique(d2_sig$gene_symbol[d2_sig$gene_symbol != ""]))
nafld_syms <- na.omit(unique(nafld_overlap$gene_symbol[nafld_overlap$gene_symbol != ""]))

ob_all_syms <- na.omit(unique(ob_all$gene_name[ob_all$gene_name != ""]))
d1_all_syms <- na.omit(unique(d1_all$gene_symbol[d1_all$gene_symbol != ""]))
d2_all_syms <- na.omit(unique(d2_all$gene_symbol[d2_all$gene_symbol != ""]))

# ============================================================
# 2. Overlap counts and Jaccard Index
# ============================================================
cat("\n=== 2. Overlap and Jaccard Index ===\n")

jaccard <- function(a, b) {
  ov <- length(intersect(a, b))
  un <- length(union(a, b))
  list(overlap = ov, union = un, jaccard = if (un > 0) round(ov / un, 6) else 0)
}

j_139_ob <- jaccard(nafld_syms, ob_syms)
j_d1_ob  <- jaccard(d1_syms, ob_syms)
j_d2_ob  <- jaccard(d2_syms, ob_syms)

cat(sprintf("139-gene NAFLD ∩ Obesity (2): overlap = %d  Jaccard = %.6f\n",
            j_139_ob$overlap, j_139_ob$jaccard))
cat(sprintf("D1 (485) ∩ Obesity (2):       overlap = %d  Jaccard = %.6f\n",
            j_d1_ob$overlap, j_d1_ob$jaccard))
cat(sprintf("D2 (1058) ∩ Obesity (2):      overlap = %d  Jaccard = %.6f\n",
            j_d2_ob$overlap, j_d2_ob$jaccard))

# ============================================================
# 3. Fisher's exact test
# ============================================================
cat("\n=== 3. Fisher's exact test ===\n")

fisher_test <- function(sig_nafld, sig_ob, all_nafld, all_ob, label) {
  universe <- intersect(all_nafld, all_ob)
  n_u  <- length(universe)
  a    <- length(intersect(sig_nafld, sig_ob))
  b    <- length(setdiff(sig_nafld, sig_ob))
  cv   <- length(setdiff(sig_ob, sig_nafld))
  d    <- max(0, n_u - a - b - cv)
  ft   <- fisher.test(matrix(c(a, b, cv, d), nrow = 2), alternative = "greater")
  cat(sprintf("  %s: universe=%d | overlap=%d | OR=%.3f | p=%.3e\n",
              label, n_u, a, ft$estimate, ft$p.value))
  data.frame(comparison = label, universe = n_u, overlap = a,
             fisher_or = round(ft$estimate, 4),
             fisher_p  = signif(ft$p.value, 4),
             stringsAsFactors = FALSE)
}

ft_139 <- fisher_test(nafld_syms, ob_syms, c(d1_all_syms, d2_all_syms), ob_all_syms,
                      "139-gene NAFLD vs Obesity")
ft_d1  <- fisher_test(d1_syms, ob_syms, d1_all_syms, ob_all_syms,
                      "D1 (485 sig) vs Obesity")
ft_d2  <- fisher_test(d2_syms, ob_syms, d2_all_syms, ob_all_syms,
                      "D2 (1058 sig) vs Obesity")

# ============================================================
# 4. Direction concordance
# ============================================================
cat("\n=== 4. Direction of obesity genes in NAFLD datasets ===\n")

ob_gene_info <- ob_sig %>%
  dplyr::select(gene = gene_name, ob_lfc = MLE_log2FC, ob_padj = padj)

d1_lookup <- d1_all %>%
  filter(gene_symbol %in% ob_syms) %>%
  group_by(gene_symbol) %>% slice(1) %>% ungroup() %>%
  dplyr::select(gene = gene_symbol, d1_lfc = lfc_mle, d1_padj = padj_mle)

d2_lookup <- d2_all %>%
  filter(gene_symbol %in% ob_syms) %>%
  group_by(gene_symbol) %>% slice(1) %>% ungroup() %>%
  dplyr::select(gene = gene_symbol, d2_lfc = lfc_mle, d2_padj = padj_mle)

direction_df <- ob_gene_info %>%
  left_join(d1_lookup, by = "gene") %>%
  left_join(d2_lookup, by = "gene") %>%
  mutate(
    in_nafld_139 = gene %in% nafld_syms,
    in_d1_sig    = gene %in% d1_syms,
    in_d2_sig    = gene %in% d2_syms,
    d1_direction_concordant = !is.na(d1_lfc) & sign(ob_lfc) == sign(d1_lfc),
    d2_direction_concordant = !is.na(d2_lfc) & sign(ob_lfc) == sign(d2_lfc)
  )

print(direction_df %>%
  dplyr::select(gene, ob_lfc, ob_padj, d1_lfc, d1_padj, d2_lfc, d2_padj,
                in_nafld_139, in_d1_sig, in_d2_sig,
                d1_direction_concordant, d2_direction_concordant))

# ============================================================
# 5. Save results
# ============================================================
summary_df <- data.frame(
  comparison    = c("139-gene NAFLD vs Obesity", "D1 (485) vs Obesity", "D2 (1058) vs Obesity"),
  nafld_sig_n   = c(length(nafld_syms), length(d1_syms), length(d2_syms)),
  obesity_sig_n = c(length(ob_syms), length(ob_syms), length(ob_syms)),
  overlap_n     = c(j_139_ob$overlap, j_d1_ob$overlap, j_d2_ob$overlap),
  union_n       = c(j_139_ob$union, j_d1_ob$union, j_d2_ob$union),
  jaccard       = c(j_139_ob$jaccard, j_d1_ob$jaccard, j_d2_ob$jaccard),
  fisher_or     = c(ft_139$fisher_or, ft_d1$fisher_or, ft_d2$fisher_or),
  fisher_p      = c(ft_139$fisher_p, ft_d1$fisher_p, ft_d2$fisher_p),
  stringsAsFactors = FALSE
)

write.csv(summary_df,   "results/overlap_summary.csv",     row.names = FALSE)
write.csv(direction_df, "results/gene_direction_table.csv", row.names = FALSE)
cat("\nSaved: results/overlap_summary.csv\n")
cat("Saved: results/gene_direction_table.csv\n")
print(summary_df)

# ============================================================
# 6. Plots
# ============================================================

# Plot A: LFC comparison for AZGP1 and CLEC4C
lfc_plot_df <- direction_df %>%
  dplyr::select(gene, ob_lfc, d1_lfc, d2_lfc) %>%
  pivot_longer(cols = c(ob_lfc, d1_lfc, d2_lfc),
               names_to = "dataset", values_to = "lfc") %>%
  mutate(
    dataset = recode(dataset,
      ob_lfc = "Obesity\n(GSE166047 adipose)",
      d1_lfc = "NAFLD D1\n(GSE162694 liver)",
      d2_lfc = "NAFLD D2\n(GSE135251 liver)"
    ),
    dataset = factor(dataset, levels = c(
      "NAFLD D1\n(GSE162694 liver)",
      "NAFLD D2\n(GSE135251 liver)",
      "Obesity\n(GSE166047 adipose)"
    )),
    fill_color = case_when(
      dataset == "Obesity\n(GSE166047 adipose)" ~ "#E65100",
      TRUE ~ "#B0BEC5"
    )
  ) %>%
  filter(!is.na(lfc))

p_lfc <- ggplot(lfc_plot_df, aes(x = dataset, y = lfc, fill = fill_color)) +
  geom_col(width = 0.55, colour = "grey30", linewidth = 0.3) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_fill_identity() +
  facet_wrap(~ gene, scales = "free_y", ncol = 2) +
  labs(
    title    = "Obesity Genes (AZGP1, CLEC4C) — Log2FC in NAFLD vs Obesity",
    subtitle = "Orange = significant in obesity  |  Grey = tested but not significant in NAFLD\nNeither gene meets NAFLD significance threshold (padj < 0.05 AND |LFC| > 1)",
    x = NULL,
    y = "MLE log2FC"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9, colour = "grey40"),
    strip.text    = element_text(face = "bold", size = 11),
    axis.text.x   = element_text(size = 9)
  )

ggsave("plots/lfc_comparison.png", p_lfc, width = 9, height = 5, dpi = 150)

# Plot B: Set overlap summary
overlap_plot_df <- data.frame(
  comparison = c("139-gene\nNAFLD overlap\nvs Obesity",
                 "NAFLD D1\n(485 genes)\nvs Obesity",
                 "NAFLD D2\n(1058 genes)\nvs Obesity"),
  nafld_only = c(length(nafld_syms) - j_139_ob$overlap,
                 length(d1_syms)    - j_d1_ob$overlap,
                 length(d2_syms)    - j_d2_ob$overlap),
  shared     = c(j_139_ob$overlap, j_d1_ob$overlap, j_d2_ob$overlap),
  ob_only    = c(length(ob_syms) - j_139_ob$overlap,
                 length(ob_syms) - j_d1_ob$overlap,
                 length(ob_syms) - j_d2_ob$overlap)
) %>%
  pivot_longer(cols = c(nafld_only, shared, ob_only),
               names_to = "category", values_to = "n_genes") %>%
  mutate(category = factor(category,
    levels = c("nafld_only", "shared", "ob_only"),
    labels = c("NAFLD only", "Shared", "Obesity only")))

p_overlap <- ggplot(overlap_plot_df,
                    aes(x = comparison, y = n_genes, fill = category)) +
  geom_col(width = 0.55, colour = "grey30", linewidth = 0.3) +
  scale_fill_manual(values = c(
    "NAFLD only"   = "#1565C0",
    "Shared"       = "#6A1B9A",
    "Obesity only" = "#E65100"
  )) +
  labs(
    title    = "NAFLD vs Obesity Significant Gene Set Comparison",
    subtitle = "Jaccard Index = 0.000 for all three comparisons (zero overlap)",
    x = NULL, y = "Number of genes", fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, colour = "grey40"),
        legend.position = "bottom")

ggsave("plots/set_overlap_summary.png", p_overlap, width = 8, height = 5, dpi = 150)

cat("\nSaved: plots/lfc_comparison.png\n")
cat("Saved: plots/set_overlap_summary.png\n")
cat("\n============================================================\n")
cat("NAFLD vs Obesity comparison complete.\n")
cat("============================================================\n")
