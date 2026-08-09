#!/usr/bin/env Rscript
# Stage 4: weighted gene co-expression network analysis (WGCNA).
#
# Usage: Rscript code/04_coexpression/run_wgcna.R --dataset GSE162694
#
# WGCNA looks for genes that rise and fall together across samples,
# independent of whether they individually reach DE significance - two
# genes can be tightly co-regulated (same module) even if only one of them
# has a fold-change large enough to pass the stage-02 significance filter.
# That is why this stage runs on the top most-variable genes genome-wide,
# not just the significant-gene list: restricting to an already-significant
# subset before building the network would only ever rediscover genes
# already known to be significant, rather than finding their regulatory
# neighbours.

suppressPackageStartupMessages({
  library(WGCNA)
  library(dplyr)
})
source(here::here("code", "utils", "paths.R"))
source(here::here("code", "utils", "analysis_helpers.R"))
allowWGCNAThreads()

args <- commandArgs(trailingOnly = TRUE)
dataset_id <- args[which(args == "--dataset") + 1]
if (is.na(dataset_id) || dataset_id == "") stop("Usage: --dataset <dataset_id>")

N_TOP_VAR_GENES <- 3000  # standard WGCNA convention: restrict to the most
                          # variable genes so the network reflects biology
                          # rather than noise from barely-expressed genes

cat(sprintf("=== WGCNA: %s ===\n", dataset_id))
out_dir <- results_path("RNAseq", dataset_id)
de <- readRDS(file.path(out_dir, "deseq2_object.rds"))

vst_mat <- de$vst_matrix
top_var_genes <- names(sort(apply(vst_mat, 1, var), decreasing = TRUE))[1:N_TOP_VAR_GENES]
expr <- t(vst_mat[top_var_genes, ])  # WGCNA expects samples x genes

# --- soft-thresholding power selection --------------------------------------
# The network is built as a weighted adjacency (correlation^power) rather
# than a hard cutoff; `power` is chosen as the smallest value at which the
# network approximates scale-free topology (R^2 >= 0.8), following the
# standard WGCNA tutorial recommendation.
sft <- pickSoftThreshold(expr, powerVector = c(1:20), verbose = 0)
power <- sft$fitIndices$Power[which(sft$fitIndices$SFT.R.sq >= 0.8)[1]]
if (is.na(power)) power <- 6  # fallback used by the WGCNA tutorial when no power reaches R^2 >= 0.8
cat("Soft-thresholding power selected:", power, "\n")
write.csv(sft$fitIndices, file.path(out_dir, "wgcna_soft_threshold.csv"), row.names = FALSE)

# --- network construction + module detection --------------------------------
net <- blockwiseModules(
  expr, power = power, TOMType = "signed", minModuleSize = 30,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE, verbose = 0
)
module_colors <- labels2colors(net$colors)

module_df <- data.frame(gene = colnames(expr), module = module_colors)
write.csv(module_df, file.path(out_dir, "wgcna_gene_modules.csv"), row.names = FALSE)
cat("Modules found:", length(unique(module_colors)), "\n")

# --- hub genes per module (highest intramodular connectivity) ---------------
adjacency <- adjacency(expr, power = power, type = "signed")
khub <- intramodularConnectivity(adjacency, module_colors)
khub$gene <- rownames(khub)
hub_genes <- khub %>%
  filter(module_colors != "grey") %>%
  mutate(module = module_colors[match(gene, colnames(expr))]) %>%
  group_by(module) %>%
  slice_max(kWithin, n = 5) %>%
  ungroup() %>%
  select(module, gene, kWithin)
write.csv(hub_genes, file.path(out_dir, "wgcna_hub_genes.csv"), row.names = FALSE)

# --- does any module overlap disproportionately with the DEG list? ----------
sig_genes <- de$res_df$gene[de$res_df$significant]
module_overlap <- module_df %>%
  group_by(module) %>%
  summarise(
    module_size = n(),
    n_sig = sum(gene %in% sig_genes),
    jaccard = jaccard_index(gene[gene %in% sig_genes], sig_genes),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sig))
write.csv(module_overlap, file.path(out_dir, "wgcna_module_deg_overlap.csv"), row.names = FALSE)

# --- dendrogram figure -------------------------------------------------------
png(file.path(out_dir, "wgcna_dendrogram.png"), width = 1000, height = 600, res = 120)
plotDendroAndColors(net$dendrograms[[1]], module_colors[net$blockGenes[[1]]],
                     "Module", dendroLabels = FALSE, main = paste("WGCNA -", dataset_id))
dev.off()

cat("Saved WGCNA modules, hub genes, and DEG overlap to", out_dir, "\n")
