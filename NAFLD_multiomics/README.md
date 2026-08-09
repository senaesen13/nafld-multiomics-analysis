# NAFLD Multi-Omics

An end-to-end, reproducible pipeline for studying NAFLD/MASLD (metabolic
dysfunction-associated steatotic liver disease) across four data modalities:
bulk RNA-seq, single-cell RNA-seq, spatial transcriptomics, and genome-scale
metabolic modelling.

## Directory structure

```
NAFLD_multiomics/
├── data/                    raw and lightly-processed input data
│   ├── RNAseq/                two active bulk liver RNA-seq cohorts
│   ├── single_cell/             GSE136103 (human liver scRNA-seq)
│   ├── spatial/                   Vu et al. 2025 Visium spatial transcriptomics
│   ├── metabolomics/                two independent NAFLD metabolomics cohorts
│   └── GEM/                           Human-GEM 2.0 metabolic network topology
├── code/                    analysis pipeline, organised by stage
│   ├── 00_download/           one script per dataset, run once
│   ├── 01_preprocessing/       gene filtering + expression matrix assembly
│   ├── 02_differential_expression/  DESeq2
│   ├── 03_enrichment/           GSEA + ORA (Hallmark, KEGG, GO:BP)
│   ├── 04_coexpression/          WGCNA + cross-cohort module comparison
│   ├── 05_cross_cohort/           overlap/concordance/meta-analysis across cohorts
│   ├── 06_drug_repositioning/      Connectivity Map query preparation
│   ├── 07_metabolic_modelling/      Reporter Metabolite analysis
│   ├── 08_single_cell/                built, not yet run (see below)
│   ├── 09_spatial/                     built, not yet run (see below)
│   ├── 10_metabolomics/                  Ji2022 + MTBLS174 metabolomics
│   └── utils/                              shared functions used by every stage
├── results/                 pipeline output, mirrors the data/ subfolders
├── backup_data/             archived cohorts, not deleted (see below)
└── document/                 methods write-ups and figures for reporting
```

Every stage after `01_preprocessing` reads its input from `results/`, not
`data/` - `data/` is raw material, `results/` is everything the pipeline
derives from it. This means any stage can be re-run in isolation as long as
the stage before it has already been run once.

## The two active bulk RNA-seq cohorts

| Dataset | Design | Samples used | Source |
|---|---|---|---|
| GSE162694 | Discovery cohort | 31 Normal / 112 NAFLD | Govaere et al. |
| GSE135251 | Validation cohort | 10 Normal / 206 NAFLD | GEO |

Both use the same significance threshold throughout the pipeline (adjusted
p < 0.05 and |log2 fold-change| > 1; see `code/utils/analysis_helpers.R`),
so gene counts and overlap statistics are directly comparable.

**History:** the project originally carried 4 cohorts. GSE130970 was dropped
early on (control arm of 4-6 samples, too small to reliably power
differential expression). GSE126848 and E-MTAB-12807 were added to replace
it. The project was then narrowed to GSE135251+GSE126848 (largest raw DEG
overlap, 178 genes), then swapped to the pairing above instead:
GSE162694+GSE135251 has slightly less raw overlap (138 genes) but is the
strongest pairing on every *cross-cohort consistency* metric - genome-wide
LFC correlation (0.38 vs -0.06), direction concordance (73.2% vs 65.2%),
and Fisher's exact test significance (9.2e-54 vs 4.4e-30) - and is Sena's
own original discovery+validation pairing. GSE126848 and E-MTAB-12807 are
archived, not deleted, to `backup_data/` (same internal folder structure -
`data/RNAseq/<id>/`, `results/RNAseq/<id>/`,
`code/00_download/download_<id>.R`). To restore either: move its three
folders back out of `backup_data/` and add its entry back to
`RNASEQ_COHORTS` in `code/utils/paths.R`. See
`backup_data/results/RNAseq/cross_cohort_4dataset_archive/` for the full
original 6-pair comparison both decisions were based on, and
`backup_data/results/RNAseq/cross_cohort_GSE135251_GSE126848_archive/` for
the intermediate pairing's results.

Worth knowing: of the two active cohorts, **both were part of the original
student project** this repository grew out of
(`bioinformatics-learning/`) - GSE162694 was her "Day 1 discovery cohort"
(the dataset her whole cross-cohort framework was originally built around)
and GSE135251 her "Day 2 validation cohort". The now-archived GSE126848,
by contrast, was never used in her original work - it was found and added
later this session, specifically to fix an underpowered control arm
elsewhere in the project.

## Running the pipeline

Each RNA-seq stage takes a `--dataset` flag and must be run once per active
cohort (`GSE162694`, `GSE135251`), except stage 00 (one script per dataset)
and stage 05 (compares all active cohorts at once, no flag needed):

```bash
# once per dataset, e.g.:
Rscript code/00_download/download_GSE162694.R

# then, once per dataset, in order:
Rscript code/01_preprocessing/run_preprocessing.R       --dataset GSE162694
Rscript code/02_differential_expression/run_deseq2.R    --dataset GSE162694
Rscript code/03_enrichment/run_enrichment.R              --dataset GSE162694
Rscript code/04_coexpression/run_wgcna.R                  --dataset GSE162694
Rscript code/06_drug_repositioning/prepare_cmap_query.R     --dataset GSE162694
Rscript code/07_metabolic_modelling/run_reporter_metabolites.R --dataset GSE162694

# once, after every active dataset has been through stage 02:
Rscript code/05_cross_cohort/run_cross_cohort.R
# then, treating the combined consensus signature as its own pseudo-dataset:
Rscript code/03_enrichment/run_enrichment.R --dataset cross_cohort
Rscript code/06_drug_repositioning/prepare_cmap_query.R --dataset cross_cohort
Rscript code/07_metabolic_modelling/run_reporter_metabolites.R --dataset cross_cohort

# module preservation between the two cohorts (needs stage 04 run on both first):
Rscript code/04_coexpression/run_module_comparison.R
```

Both active cohorts, plus the combined consensus, have already been run
through every stage; `results/` contains real output, not placeholders.

## Cross-cohort headline result

The Stouffer-weighted meta-analysis combining both cohorts (see the comment
in `code/05_cross_cohort/run_cross_cohort.R` for why this is more principled
than a simple "significant in both" vote) finds **381 significant genes**
(330 up, 51 down) - `results/RNAseq/cross_cohort/meta_analysis_results.csv`.
That consensus signature feeds enrichment (178 significant GSEA gene sets,
13 ORA terms), reporter metabolites (top hit: 4Fe4S iron-sulfur cluster,
p=5.5e-06 - relevant to mitochondrial electron transport chain function),
and the CMap drug-repositioning query (150 up / 51 down genes).

Unlike the intermediate GSE135251+GSE126848 pairing this project passed
through, this pairing has *strong* cross-cohort consistency: genome-wide
log2FC correlation of 0.38 (vs -0.06) and 73.2% direction concordance among
shared significant genes (vs 65.2%) - see "History" above for the full
comparison. The trade-off runs the other way from before: fewer down-
regulated genes survive combination here (51, vs 125 in the archived
pairing) because GSE162694 alone only had 15 down-regulated significant
genes to begin with (see `results/RNAseq/GSE162694/deseq2_significant_genes.csv`) -
a real asymmetry in that cohort's power to detect down-regulation, not an
artefact of the meta-analysis.

## Single-cell and spatial - pipelines built, not run

`code/08_single_cell/` (3 stages) and `code/09_spatial/` (4 stages) are
complete, real analysis code - not stubs - written against the real
downloaded data/metadata but deliberately not executed. Both use one shared
liver cell-type marker reference (`code/utils/liver_cell_markers.R`) so a
cluster label means the same thing in either modality.

- Single-cell: GSE136103 is not NAFLD-specific (only 7/24 libraries are
  NAFLD-caused; others are alcohol/hemochromatosis/PBC-related cirrhosis
  and PBMC blood samples) - stage 1 restricts to the correct 15-library/
  7-donor subset from the real metadata fields. Stage 3 uses pseudobulk
  DE (sum counts per donor per cell type, not per-cell testing) to avoid
  pseudo-replication. Needs `scDblFinder` installed before it can run
  (not yet installed here).
- Spatial: `data/spatial/Vu2025_Visium/` requires a manual download step
  first (UQ eSpace rejects scripted downloads - see
  `code/00_download/download_spatial_Vu2025.R`). Deconvolution uses real
  RCTD (`spacexr`, confirmed installed) rather than the `AddModuleScore`
  fallback used in the predecessor project. Two integration points
  (RCTD's exact input object shape; a METAFlux internal error hit during
  a package smoke-test) are flagged in-script as needing verification
  before first run, rather than assumed correct.

## Metabolomics

`code/10_metabolomics/run_ji2022_metabolomics.R` re-groups the Ji et al.
2022 plasma metabolomics cohort (published as Control n=25 / NAFL n=42 /
NASH n=19 group-level summary statistics, not raw per-patient data) into
this project's Control-vs-NAFLD convention by pooling NAFL+NASH with the
statistically correct combined-subgroup formula (not just re-labelling),
then re-deriving significance via Welch's t-test from the pooled summary
stats. Result: 14/79 metabolites significant, including alpha-ketoglutarate
(TCA cycle) - this is the same metabolite the reporter-metabolite analysis
on the original *4-cohort* consensus independently flagged as its top hit.
That specific cross-modality match is a property of the archived 4-cohort
consensus, not the current one - the current GSE162694+GSE135251 consensus's
top reporter-metabolite hit is a 4Fe4S iron-sulfur cluster instead (see
`results/RNAseq/cross_cohort/reporter_metabolites.csv`; the original
4-cohort finding is preserved at
`backup_data/results/RNAseq/cross_cohort_4dataset_archive/reporter_metabolites.csv`).
Checked directly: in the current consensus, alpha-ketoglutarate ranks
105th of 2224 metabolites tested (Z=1.84, padj=0.69 - not significant).
The cross-modality match with Ji2022 was real for the archived 4-cohort
consensus but does not hold up in the current 2-cohort one - don't cite it
as a finding of the active pipeline.

`code/10_metabolomics/run_mtbls174_metabolomics.R` covers the second
metabolomics cohort. As agreed, kept as its own within-disease High-vs-Low
steatosis severity analysis, not forced into the Control-vs-NAFLD
framework - every subject in MTBLS174 is already NAFLD (no healthy control
arm at all), so "Low steatosis" is the mild end of the disease spectrum,
not a stand-in for "healthy". Method: limma moderated t-test on log2
NMR intensities (18 samples, 16 metabolites; Low n=8 vs High n=3,
Medium n=7 excluded from the test) - reproduces the original analysis in
`bioinformatics-learning/week4-day5-nafld-metabolomics/` exactly (verified
decimal-for-decimal), with one bug fixed along the way: that script's
steatosis-range parser intended "<5" to be treated as 0 (per its own
comment) but `gsub("<", "0", "<5")` actually produces the string "05" -> 5,
not 0 - fixed to match the stated intent. Only one sample in this cohort is
affected, and it doesn't change that sample's Low/Medium/High bucket
either way. Result: 0/16 metabolites significant (padj<0.05, |logFC|>0.5)
- an honest null result given only 3 samples in the High-steatosis group,
not a bug (matches what the original, smaller-sample-size-aware analysis
also found).

## GEM data

`data/GEM/human_gem_topology_network.csv` is the generic (non-context-
specific) Human-GEM 2.0 metabolite-gene network (Robinson et al. 2020),
used as-is by stage 07. Context-specific model extraction via tINIT
(RAVEN/MATLAB, both confirmed installed and runnable on this machine) was
considered and deliberately deferred: the mandatory model open/close
workflow on file references `getINITModel2`, but the installed RAVEN
2.10.3 only exposes `getINITModel` (no "2") and a newer `ftINIT` - whether
the same sequencing rules apply was not verified, so no tINIT code has
been written rather than guessing.
