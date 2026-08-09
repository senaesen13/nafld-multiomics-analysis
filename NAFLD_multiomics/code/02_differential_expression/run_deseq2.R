#!/usr/bin/env Rscript
# Stage 2: differential expression (NAFLD vs Normal) with DESeq2.
#
# Usage: Rscript code/02_differential_expression/run_deseq2.R --dataset GSE162694
#
# Two log2 fold-change estimates are kept side by side throughout this
# project, on purpose:
#   - MLE (maximum likelihood): the "raw" estimate. Used for the
#     significance filter and for GSEA ranking, because it does not depend
#     on how many other genes happen to be in the dataset.
#   - apeglm-shrunk: pulls noisy, low-count-driven estimates toward zero
#     (Zhu, Ibrahim & Love 2019). Used when *reporting* an effect size,
#     because the raw MLE for a low-count gene can be enormous and
#     essentially meaningless.
# Filtering on the shrunk value and reporting the raw value (or vice versa)
# both happen in the literature and give different gene lists - this
# project is consistent about which is used for which purpose everywhere.

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(ggplot2)
})
source(here::here("code", "utils", "paths.R"))
source(here::here("code", "utils", "analysis_helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
dataset_id <- args[which(args == "--dataset") + 1]
if (is.na(dataset_id) || dataset_id == "") stop("Usage: --dataset <dataset_id>")

cat(sprintf("=== Differential expression: %s ===\n", dataset_id))
prep <- readRDS(file.path(results_path("RNAseq", dataset_id), "preprocessed.rds"))

dds <- DESeqDataSetFromMatrix(
  countData = prep$counts,
  colData   = data.frame(condition = prep$condition, row.names = colnames(prep$counts)),
  design    = ~condition
)
dds <- DESeq(dds)

res_mle <- results(dds, contrast = c("condition", "NAFLD", "Normal"), alpha = 0.05)
res_shrunk <- lfcShrink(dds, coef = resultsNames(dds)[grepl("NAFLD", resultsNames(dds))], type = "apeglm")

res_df <- data.frame(
  gene           = rownames(res_mle),
  baseMean       = res_mle$baseMean,
  log2FoldChange = res_mle$log2FoldChange,   # MLE - used for significance calls + GSEA ranking
  lfcSE          = res_mle$lfcSE,
  pvalue         = res_mle$pvalue,
  padj           = res_mle$padj,
  lfc_shrunk     = res_shrunk$log2FoldChange  # apeglm - used only for reporting effect size
) %>% flag_significant()

out_dir <- results_path("RNAseq", dataset_id)
write.csv(res_df, file.path(out_dir, "deseq2_results_all.csv"), row.names = FALSE)
write.csv(subset(res_df, significant), file.path(out_dir, "deseq2_significant_genes.csv"), row.names = FALSE)

vsd <- vst(dds, blind = FALSE)
saveRDS(list(dds = dds, vst_matrix = assay(vsd), res_df = res_df), file.path(out_dir, "deseq2_object.rds"))

n_sig <- sum(res_df$significant, na.rm = TRUE)
cat(sprintf("%d / %d genes significant (padj<%.2f, |LFC|>%d): %d up, %d down\n",
            n_sig, nrow(res_df), SIG_PADJ, SIG_LFC,
            sum(res_df$significant & res_df$log2FoldChange > 0, na.rm = TRUE),
            sum(res_df$significant & res_df$log2FoldChange < 0, na.rm = TRUE)))

# --- standard figures --------------------------------------------------
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
p_pca <- ggplot(pca_data, aes(PC1, PC2, colour = condition)) +
  geom_point(size = 2.5, alpha = 0.8) +
  labs(title = paste("PCA -", dataset_id), x = "PC1", y = "PC2") +
  theme_project()
ggsave(file.path(out_dir, "pca.png"), p_pca, width = 6, height = 5, dpi = 150)

volcano_df <- res_df %>% mutate(neg_log10_padj = -log10(pmax(padj, 1e-300)))
p_volcano <- ggplot(volcano_df, aes(log2FoldChange, neg_log10_padj, colour = significant)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_colour_manual(values = c(`TRUE` = "#D64550", `FALSE` = "grey70")) +
  labs(title = paste("Volcano -", dataset_id), x = "log2 fold change (NAFLD / Normal)",
       y = expression(-log[10](adjusted~p))) +
  theme_project()
ggsave(file.path(out_dir, "volcano.png"), p_volcano, width = 6, height = 5, dpi = 150)

cat("Saved results, VST matrix, and figures to", out_dir, "\n")
