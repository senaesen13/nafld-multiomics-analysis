#!/usr/bin/env Rscript
# Stage 1: load raw counts, restrict to protein-coding genes, collapse to
# gene symbol, and apply a minimum-expression filter.
#
# Usage: Rscript code/01_preprocessing/run_preprocessing.R --dataset GSE162694
#
# Why filter to protein-coding genes before differential expression? Non-
# coding RNAs and pseudogenes are real biology, but including them roughly
# doubles the number of hypothesis tests without adding statistical power,
# which costs significant genes at the multiple-testing correction step for
# no benefit if the downstream analysis (GSEA, WGCNA, GEM mapping) all
# operate on protein-coding gene sets anyway.
#
# Why collapse to gene symbol rather than keep Ensembl IDs? The four cohorts
# were quantified against different references (three against Ensembl gene
# IDs, one via transcript-level Salmon output aggregated with a different
# annotation). Gene symbol is the only identifier space common to all four,
# and cross-cohort comparison (stage 05) depends on that.

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(dplyr)
})
source(here::here("code", "utils", "paths.R"))
source(here::here("code", "01_preprocessing", "load_dataset.R"))

args <- commandArgs(trailingOnly = TRUE)
dataset_id <- args[which(args == "--dataset") + 1]
if (is.na(dataset_id) || dataset_id == "") stop("Usage: --dataset <GSE162694|GSE135251|GSE126848|EMTAB12807>")

MIN_MEAN_COUNT <- 10  # matches the expression filter used throughout this project

cat(sprintf("=== Preprocessing %s ===\n", dataset_id))
d <- load_aligned(dataset_id)
cat(sprintf("Loaded: %d genes x %d samples (%s)\n",
            nrow(d$counts), ncol(d$counts), paste(table(d$condition), collapse = " / ")))

# --- map gene IDs to symbol, restricted to protein-coding -------------------
is_ensembl <- grepl("^ENSG", rownames(d$counts)[1])
key_type   <- if (is_ensembl) "ENSEMBL" else "SYMBOL"

gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(rownames(d$counts)),
  keytype = key_type,
  columns = c("SYMBOL", "GENETYPE")
) %>%
  filter(GENETYPE == "protein-coding", !is.na(SYMBOL)) %>%
  distinct(.data[[key_type]], .keep_all = TRUE)

counts_pc <- d$counts[rownames(d$counts) %in% gene_map[[key_type]], , drop = FALSE]
symbol_for_row <- gene_map$SYMBOL[match(rownames(counts_pc), gene_map[[key_type]])]

# Multiple input IDs can map to the same symbol (e.g. retired/merged Ensembl
# IDs) - sum their counts rather than silently keeping only one.
counts_by_symbol <- rowsum(counts_pc, group = symbol_for_row)
cat(sprintf("After protein-coding filter + symbol collapse: %d genes\n", nrow(counts_by_symbol)))

# --- minimum expression filter ----------------------------------------------
keep <- rowMeans(counts_by_symbol) >= MIN_MEAN_COUNT
counts_filt <- counts_by_symbol[keep, , drop = FALSE]
cat(sprintf("After mean-count >= %d filter: %d genes\n", MIN_MEAN_COUNT, nrow(counts_filt)))

# --- save --------------------------------------------------------------------
out_dir <- results_path("RNAseq", dataset_id)
saveRDS(
  list(counts = round(counts_filt), condition = d$condition, dataset_id = dataset_id),
  file.path(out_dir, "preprocessed.rds")
)
cat("Saved:", file.path(out_dir, "preprocessed.rds"), "\n")
