#!/usr/bin/env Rscript
# Stage 3: pathway enrichment - GSEA (Hallmark, KEGG, GO:BP) and ORA
# (over-representation analysis on the significant gene lists).
#
# Usage: Rscript code/03_enrichment/run_enrichment.R --dataset GSE162694
#
# GSEA and ORA answer different questions and are both reported:
#   - GSEA uses every tested gene, ranked, and asks "are genes in this
#     pathway systematically shifted toward one end of the ranking?" - it
#     does not require a significance cutoff and can detect coordinated
#     small shifts that no single gene would reach significance for.
#   - ORA only uses the genes that already passed the significance filter
#     and asks "is this pathway over-represented in that fixed list,
#     relative to chance?" - simpler, but throws away the genes that were
#     close to but did not cross the threshold.

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(msigdbr)
  library(dplyr)
  library(ggplot2)
})
source(here::here("code", "utils", "paths.R"))
source(here::here("code", "utils", "analysis_helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
dataset_id <- args[which(args == "--dataset") + 1]
if (is.na(dataset_id) || dataset_id == "") stop("Usage: --dataset <dataset_id>")

cat(sprintf("=== Enrichment: %s ===\n", dataset_id))
out_dir <- results_path("RNAseq", dataset_id)
de <- readRDS(file.path(out_dir, "deseq2_object.rds"))
res_df <- de$res_df

# --- gene set collections (see code/utils/analysis_helpers.R comment on
#     msigdbr's collection/subcollection naming) -----------------------
gene_sets <- list(
  Hallmark = msigdbr(species = "human", collection = "H"),
  KEGG     = msigdbr(species = "human", collection = "C2", subcollection = "CP:KEGG_LEGACY"),
  GO_BP    = msigdbr(species = "human", collection = "C5", subcollection = "GO:BP")
) %>% lapply(\(x) x %>% select(gs_name, gene_symbol) %>% distinct())

# --- GSEA ----------------------------------------------------------------
ranking <- build_gsea_ranking(res_df, gene_col = "gene")
write_rnk(ranking, file.path(out_dir, "gsea_ranked_genes.rnk"))

gsea_results <- lapply(names(gene_sets), function(set_name) {
  fit <- GSEA(ranking, TERM2GENE = gene_sets[[set_name]], pvalueCutoff = 0.05,
              eps = 0, seed = TRUE)
  as.data.frame(fit) %>% mutate(collection = set_name)
})
gsea_all <- bind_rows(gsea_results)
write.csv(gsea_all, file.path(out_dir, "gsea_results.csv"), row.names = FALSE)
cat(sprintf("GSEA: %d significant gene sets across Hallmark/KEGG/GO:BP\n", nrow(gsea_all)))

# --- ORA (up- and down-regulated genes tested separately) -----------------
sig <- subset(res_df, significant)
up_genes   <- sig$gene[sig$log2FoldChange > 0]
down_genes <- sig$gene[sig$log2FoldChange < 0]
universe   <- res_df$gene

ora_results <- lapply(names(gene_sets), function(set_name) {
  t2g <- gene_sets[[set_name]]
  up   <- enricher(up_genes,   TERM2GENE = t2g, universe = universe, pvalueCutoff = 0.05)
  down <- enricher(down_genes, TERM2GENE = t2g, universe = universe, pvalueCutoff = 0.05)
  bind_rows(
    if (!is.null(up))   as.data.frame(up)   %>% mutate(direction = "Up in NAFLD")   else NULL,
    if (!is.null(down)) as.data.frame(down) %>% mutate(direction = "Down in NAFLD") else NULL
  ) %>% mutate(collection = set_name)
})
ora_all <- bind_rows(ora_results)
write.csv(ora_all, file.path(out_dir, "ora_results.csv"), row.names = FALSE)
cat(sprintf("ORA: %d significant terms across Hallmark/KEGG/GO:BP\n", nrow(ora_all)))

# --- figure: top GSEA hits by collection ----------------------------------
if (nrow(gsea_all) > 0) {
  top_hits <- gsea_all %>% group_by(collection) %>% slice_min(p.adjust, n = 8) %>% ungroup()
  p <- ggplot(top_hits, aes(x = NES, y = reorder(Description, NES), fill = collection)) +
    geom_col() +
    labs(title = paste("Top GSEA hits -", dataset_id), x = "Normalized enrichment score", y = NULL) +
    theme_project(base_size = 9)
  ggsave(file.path(out_dir, "gsea_top_hits.png"), p, width = 9, height = 7, dpi = 150)
}

cat("Saved GSEA + ORA results to", out_dir, "\n")
