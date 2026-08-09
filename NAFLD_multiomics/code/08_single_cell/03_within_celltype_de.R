#!/usr/bin/env Rscript
# Single-cell stage 3: within-cell-type differential expression, Healthy vs
# NAFLD, using pseudobulk aggregation.
#
# NOT YET RUN in this project - see the note at the top of
# 01_qc_and_doublets.R. Depends on stage 2's output.
#
# Usage: Rscript code/08_single_cell/03_within_celltype_de.R
#
# Why pseudobulk (sum counts per donor per cell type, then DESeq2) rather
# than testing cell-by-cell (e.g. Seurat's default Wilcoxon FindMarkers)?
# Cells from the same donor are not independent observations - they share
# that donor's genetics, disease severity, and technical batch. Treating
# each of a donor's thousands of cells as an independent replicate is
# pseudo-replication: it collapses "5 healthy donors vs 2 NAFLD donors" into
# an artificial "thousands of healthy cells vs thousands of NAFLD cells"
# comparison, which produces a flood of significant genes driven by
# between-donor variation, not the disease effect (Squair et al. 2021,
# Nature Communications, "Confronting false discoveries in single-cell
# differential expression"). Pseudobulk restores the correct unit of
# replication (donor), at the cost of only having as much statistical power
# as the number of donors - which in this cohort (5 Healthy, 2 NAFLD, see
# stage 1's comment on why the NAFLD arm is this small) is honestly limited,
# and any results from this stage should be read as hypothesis-generating,
# not confirmatory, until a larger single-cell NAFLD cohort is available.

suppressPackageStartupMessages({
  library(Seurat)
  library(DESeq2)
  library(dplyr)
})
source(here::here("code", "utils", "paths.R"))
source(here::here("code", "utils", "analysis_helpers.R"))

out_dir <- results_path("single_cell", "GSE136103")
obj <- readRDS(file.path(out_dir, "02_annotated.rds"))

MIN_DONORS_PER_CONDITION <- 2  # below this, DESeq2's dispersion estimate is unreliable

# --- build one pseudobulk count matrix per cell type --------------------------
cell_types <- unique(obj$cell_type)
donor_condition <- obj@meta.data |> distinct(donor, condition)

results_by_celltype <- list()
for (ct in cell_types) {
  sub <- subset(obj, subset = cell_type == ct)
  donors_present <- unique(sub$donor)
  cond_counts <- table(donor_condition$condition[donor_condition$donor %in% donors_present])
  if (length(cond_counts) < 2 || any(cond_counts < MIN_DONORS_PER_CONDITION)) {
    cat(sprintf("Skipping %s: insufficient donors per condition (%s)\n",
                ct, paste(names(cond_counts), cond_counts, sep = "=", collapse = ", ")))
    next
  }

  raw_counts <- GetAssayData(sub, assay = "RNA", layer = "counts")
  pseudobulk <- sapply(donors_present, function(d) {
    Matrix::rowSums(raw_counts[, sub$donor == d, drop = FALSE])
  })
  colnames(pseudobulk) <- donors_present

  condition <- donor_condition$condition[match(donors_present, donor_condition$donor)]
  dds <- DESeqDataSetFromMatrix(
    countData = round(pseudobulk),
    colData   = data.frame(condition = factor(condition, levels = c("Healthy", "NAFLD"))),
    design    = ~condition
  )
  dds <- dds[rowSums(counts(dds)) >= 10, ]  # same expression-level pre-filter used throughout this project
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", "NAFLD", "Healthy"))

  res_df <- data.frame(
    gene = rownames(res), log2FoldChange = res$log2FoldChange,
    pvalue = res$pvalue, padj = res$padj, cell_type = ct
  ) |> flag_significant()

  write.csv(res_df, file.path(out_dir, paste0("03_de_", gsub("[^A-Za-z0-9]+", "_", ct), ".csv")),
            row.names = FALSE)
  results_by_celltype[[ct]] <- res_df
  cat(sprintf("%-28s %d donors, %d significant genes\n", ct, ncol(pseudobulk), sum(res_df$significant, na.rm = TRUE)))
}

all_results <- bind_rows(results_by_celltype)
write.csv(all_results, file.path(out_dir, "03_de_all_celltypes.csv"), row.names = FALSE)

# --- where do the bulk RNA-seq consensus genes (stage 05) come from? --------
# The 67-gene, >=3/4-cohort consensus set from the bulk pipeline was found
# without knowing which liver cell type drives each gene. Once this stage
# has actually been run, that answer follows directly by joining on gene
# symbol - showing e.g. whether TREM2/SPP1/FASN/COL1A1 (the WGCNA turquoise
# hub genes from stage 04) are Kupffer-cell-driven, hepatocyte-driven, or
# both.

cat("Saved within-cell-type DE results to", out_dir, "\n")
