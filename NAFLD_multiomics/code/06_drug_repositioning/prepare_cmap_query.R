#!/usr/bin/env Rscript
# Stage 6: prepare a Connectivity Map (CMap) drug-repositioning query.
#
# Usage: Rscript code/06_drug_repositioning/prepare_cmap_query.R --dataset GSE162694
#
# The idea behind CMap-based drug repositioning (Lamb et al. 2006, Science):
# submit the top up- and down-regulated genes of a disease signature to the
# Broad Institute's LINCS L1000 database, and look for compounds whose own
# transcriptional signature is the *mirror image* - genes the disease turns
# up, the compound turns down, and vice versa. A strong "reversal" is a
# candidate for repositioning.
#
# This script produces the two query files, correctly formatted for
# submission at https://clue.io/query. It does not produce connectivity
# scores itself: CMap querying is an interactive step against the Broad's
# own compound reference database (not an offline computation this pipeline
# can reproduce), so submitting the files and reading back the results is a
# manual step - see cmap_query_readme.txt for exactly what to do with the
# two .grp files.
#
# Cell line: L1000's core reference panel (Touchstone) has 9 cell lines with
# dense, reliable signature coverage, of which HEPG2 (human hepatocellular
# carcinoma) is the only liver-derived one - HUH7 and other liver lines are
# not in the core panel. clue.io's query tool reports both a HEPG2-specific
# score and a summary score pooled across all 9 core lines: a compound whose
# reversal is strong in HEPG2 specifically but weak in the pooled summary is
# a more liver-mechanism-specific candidate than one that reverses broadly
# across unrelated tissue types (skin, prostate, lung, etc. are among the
# other 8 core lines) - read both, not just the summary score.

suppressPackageStartupMessages(library(dplyr))
source(here::here("code", "utils", "paths.R"))
source(here::here("code", "utils", "analysis_helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
dataset_id <- args[which(args == "--dataset") + 1]
if (is.na(dataset_id) || dataset_id == "") stop("Usage: --dataset <dataset_id>")

N_GENES <- 150  # CMap's own submission guidance: 150 up + 150 down is the
                 # standard signature size used in the original methodology

cat(sprintf("=== CMap query preparation: %s ===\n", dataset_id))
out_dir <- results_path("RNAseq", dataset_id)
res_df  <- readRDS(file.path(out_dir, "deseq2_object.rds"))$res_df

sig <- res_df |> filter(significant)
# Must filter to the correct direction *before* ranking: if a cohort has
# fewer than N_GENES truly down-regulated significant genes, sorting the
# whole significant set ascending and taking head(N_GENES) would silently
# pad the "down" file with the least-upregulated genes instead - which are
# not down-regulated at all.
up_genes   <- sig |> filter(log2FoldChange > 0) |> arrange(desc(log2FoldChange)) |> pull(gene) |> head(N_GENES)
down_genes <- sig |> filter(log2FoldChange < 0) |> arrange(log2FoldChange)       |> pull(gene) |> head(N_GENES)

writeLines(up_genes,   file.path(out_dir, "cmap_up_genes.grp"))
writeLines(down_genes, file.path(out_dir, "cmap_down_genes.grp"))

writeLines(c(
  paste("CMap query files for", dataset_id),
  "",
  "To get real connectivity scores:",
  "  1. Go to https://clue.io/query (free registration required)",
  "  2. Upload cmap_up_genes.grp as the 'up' gene set",
  "  3. Upload cmap_down_genes.grp as the 'down' gene set",
  "  4. Submit against the L1000 dataset",
  "  5. In the results, check BOTH the HEPG2-specific score and the",
  "     pooled summary score across all 9 core cell lines (HEPG2 is the",
  "     only liver line in the core panel). Compounds with the most",
  "     negative connectivity score in HEPG2 specifically are the",
  "     strongest liver-relevant candidates: their transcriptional effect",
  "     is most opposite to the disease signature submitted here, in the",
  "     tissue context that matters for this disease.",
  "",
  sprintf("Signature size: %d up-regulated, %d down-regulated genes (padj<%.2f, |LFC|>%d, top %d by effect size).",
          length(up_genes), length(down_genes), SIG_PADJ, SIG_LFC, N_GENES)
), file.path(out_dir, "cmap_query_readme.txt"))

cat(sprintf("Wrote %d up / %d down gene CMap query files to %s\n",
            length(up_genes), length(down_genes), out_dir))
