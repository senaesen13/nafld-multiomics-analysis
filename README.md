# NAFLD Multi-Omics Analysis

Transcriptomic characterisation of Non-Alcoholic Fatty Liver Disease (NAFLD) across
three independent bulk RNA-seq cohorts and one single-cell RNA-seq dataset, with
cross-cohort validation and cell-type resolution.

---

## Research Question

What are the reproducible transcriptional signatures of NAFLD across independent human
liver cohorts, and which genes are robustly dysregulated regardless of fibrosis stage,
platform, or cohort composition?

---

## Datasets

| Folder | GEO Accession | Type | Samples (Normal / NAFLD) | Role |
|---|---|---|---|---|
| `01-NAFLD-discovery-cohort-GSE162694/` | GSE162694 | Bulk RNA-seq | 31 / 112 | Discovery |
| `02-NAFLD-validation-cohort-GSE135251/` | GSE135251 | Bulk RNA-seq | 10 / 206 | Primary validation |
| `03-NAFLD-second-validation-GSE130970/` | GSE130970 | Bulk RNA-seq | 4 / 74 | Second validation |
| `04-NAFLD-cross-cohort-comparison/` | — | Analysis | — | Pairwise overlap + pathway concordance |
| `05-NAFLD-scRNA-seq-GSE136103/` | GSE136103 | scRNA-seq | 5 donors | Cell-type resolution |
| `06-NAFLD-coexpression-analysis/` | — | Analysis | — | Coexpression modules from 139 shared genes |

---

## Key Findings

### Gene-Level Replication (Discovery → Validation)

139 genes are significantly dysregulated in both GSE162694 (Dataset 1) and GSE135251
(Dataset 2). This overlap is 3.5× above chance (Fisher's OR = 5.11, p = 9.5e-43). 73.4%
of overlap genes change in the same direction in both datasets (binomial p = 1.6e-08).

**TREM2** (lipid-associated macrophage marker) is upregulated in NAFLD in both datasets:
log2FC +2.48 (Dataset 1) and +2.52 (Dataset 2). **SPP1** is similarly upregulated: log2FC
+1.87 (Dataset 1) and +1.41 (Dataset 2). Both genes are individually significant in both
cohorts — this only holds for the Dataset 1 vs Dataset 2 pairing; Dataset 3 lacks the
statistical power to confirm either gene (n = 4 normal controls).

### Cross-Cohort Pathway Concordance (GSEA)

GSEA was run using a signed significance ranking metric: **sign(lfc_apeglm) × −log10(pvalue_mle)**,
which weights genes by both direction and statistical confidence.

Two Hallmark gene sets are activated in NAFLD across all three cohorts and all three
pairwise comparisons: **MYOGENESIS** and **APICAL_JUNCTION**. These are the most
fibrosis-stage-independent pathway signals in this study.

Inflammatory pathways (TNF, IL-17, NF-kB) are dataset-dependent — likely reflecting
fibrosis-stage composition differences between cohorts rather than a core NAFLD effect.

### Single-Cell Resolution (GSE136103)

scRNA-seq of 5 NAFLD donors confirms TREM2+ lipid-associated macrophages (LAMs) as a
distinct, expanded macrophage subpopulation in NAFLD liver. SPP1 is a marker of a
separate macrophage cluster. GPNMB co-localises with TREM2+ LAMs.

---

## Why Dataset 1 vs Dataset 2 is the Primary Comparison

| Comparison | Overlap | Fisher's p | Direction concordance |
|---|---|---|---|
| A: Dataset 3 vs Dataset 1 | 38 | 8.7e-21 | 94.7% |
| B: Dataset 3 vs Dataset 2 | 52 | 1.6e-16 | 84.6% |
| **C: Dataset 1 vs Dataset 2** | **139** | **9.5e-43** | **73.4%** |

Comparisons A and B are limited by Dataset 3's n = 4 normal controls, which yields only
180 significant genes. The 139-gene overlap from Comparison C is nearly 4× larger than
Comparison A's 38 genes and provides a workable input for downstream GEM and coexpression
network analysis. These 139 shared genes are used as input for the next analysis phase.

---

## Folder Guide

```
01-NAFLD-discovery-cohort-GSE162694/
    scripts/deseq2_analysis.R       DESeq2 pipeline (GSE162694)
    scripts/gsea_analysis.R         GSEA: KEGG + Hallmark
    results/                        DESeq2 results, significant genes, GSEA CSVs
    plots/                          PCA, volcano, GSEA dotplots/ridgeplots
    NOTES.md                        Methods, results, interpretation

02-NAFLD-validation-cohort-GSE135251/
    scripts/validation_gse135251.R  DESeq2 pipeline + Day1 vs Day2 overlap
    scripts/gsea_analysis.R         GSEA: KEGG + Hallmark
    results/                        DESeq2 results, overlap summary, GSEA CSVs
    plots/                          PCA, volcano, LFC scatter, GSEA plots
    NOTES.md                        Methods, results, 139 shared genes, TREM2/SPP1

03-NAFLD-second-validation-GSE130970/
    scripts/deseq2_analysis.R       DESeq2 pipeline (GSE130970, Entrez IDs)
    scripts/gsea_analysis.R         GSEA: KEGG + Hallmark
    results/                        DESeq2 results, significant genes (180), GSEA CSVs
    plots/                          PCA, volcano, GSEA plots
    NOTES.md                        Methods, n=4 caveat, GSEA results

04-NAFLD-cross-cohort-comparison/
    scripts/pairwise_comparison.R   Comparisons A + B (C reuses Day2 numbers)
    results/pairwise_comparison.md  Full report: all three pairwise comparisons
    results/pairwise_summary.csv    One row per comparison: Fisher's p, concordance
    results/lfc_corr_*.png          LFC scatter plots for each pairing
    NOTES.md                        Summary tables, why Comparison C was chosen

05-NAFLD-scRNA-seq-GSE136103/
    scripts/                        Seurat pipeline: QC → clustering → annotation
    plots/                          UMAP, dotplots, feature plots (TREM2, SPP1, GPNMB)
    results/                        Cell-type composition, within-cluster DE
    NOTES.md                        Methods, cell-type annotation, key findings

06-NAFLD-coexpression-analysis/
    scripts/coexpression_analysis.R VST normalisation → Pearson correlation → hierarchical clustering
    results/vst_matrix_139genes.csv VST expression of 139 genes × 143 samples
    results/correlation_matrix.csv  139 × 139 Pearson correlation matrix
    results/module_assignments.csv  Gene → module mapping (k = 7)
    plots/coexpression_heatmap.png  Correlation heatmap ordered by module
    NOTES.md                        Methods, M1/M2 module characterisation, key gene report
```

---

## Tools and Key Packages

- **DESeq2** — differential expression (MLE + apeglm LFC shrinkage)
- **clusterProfiler** — GSEA (gseKEGG + GSEA with MSigDB Hallmark)
- **msigdbr** — MSigDB Hallmark gene sets
- **org.Hs.eg.db / biomaRt** — gene ID mapping
- **Seurat** — scRNA-seq processing and clustering
- **ggplot2 / ggrepel / enrichplot** — visualisation
