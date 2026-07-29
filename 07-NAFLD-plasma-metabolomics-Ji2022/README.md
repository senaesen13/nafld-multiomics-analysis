# Plasma Metabolomics of NAFLD Progression — Ji et al. 2022

## Research Question

Which plasma metabolites are significantly elevated (or depleted) as liver disease
progresses from healthy controls through non-alcoholic fatty liver (NAFL) to
non-alcoholic steatohepatitis (NASH), and do those differences survive multiple-testing
correction across a panel of 79 measured metabolites?

---

## Data Source

**Paper:** Ji S, et al. "Plasma Metabolomics Provides Insight into Non-Alcoholic Fatty Liver
Disease Pathogenesis." *Biomedicines* 2022, 10(7):1669.  
**DOI:** [10.3390/biomedicines10071669](https://doi.org/10.3390/biomedicines10071669)

**Cohort:** 86 fasting plasma samples — Control (n = 25), NAFL (n = 42), NASH (n = 19).  
**Metabolite panel:** 79 plasma metabolites measured by GC-MS and LC-MS (33 amino acids,
4 kynurenine-pathway metabolites, 4 nucleosides, 18 organic acids, 20 fatty acids).

**Data file:** `data/Table_S1_Metabolite_Levels.csv` — the paper's Supplementary Table 1,
containing group means, standard deviations, and pre-calculated Kruskal-Wallis p-values
and BH-FDR q-values for all 79 metabolites. Individual patient-level measurements were
not released with the paper.

---

## Methods

1. **Data inspection.** Confirmed the CSV structure before analysis: 79 rows (one per
   metabolite, not one per patient), 16 columns including group means, SDs, and pre-computed
   statistics.

2. **Quality control.** Checked for missing values (none), verified group sizes against
   the paper, and identified the known aspartic acid outlier (~382.4 µg/uL in one NASH
   patient). The paper retained this sample for the group comparisons (Kruskal-Wallis is
   rank-based and robust to outliers) but excluded it before building machine learning
   classifiers; after exclusion, the combined NAFLD mean and median aspartic acid were
   25.4 and 16.9 µg/uL respectively.

3. **Univariate statistics.** BH-FDR correction was independently re-applied to the
   provided Kruskal-Wallis p-values using `statsmodels.multipletests` (method = `fdr_bh`)
   as a cross-check against the paper's reported Qval_FDR column. Pairwise group
   comparisons (Control vs NAFL, Control vs NASH, NAFL vs NASH) were taken directly from
   the paper's table.

4. **Trend assessment.** A Jonckheere-Terpstra ordered-groups test cannot be run without
   individual patient data. As a proxy, Spearman correlation between group order
   (Control = 0, NAFL = 1, NASH = 2) and group means was computed to characterise
   monotone direction.

5. **Visualisation.** Mean ± SD bar plots for the six significant metabolites; a
   pseudo-Z-score heatmap across all 79 metabolites (Z-scores computed from the three
   group means per metabolite); a volcano plot (log₂ fold change NASH/Control vs
   −log₁₀ FDR q-value).

6. **Multivariate (from published paper).** The paper ran PLS-DA on log-transformed,
   autoscaled individual-level data. Their reported model performance: R² = 0.753,
   Q² = 0.341, classification accuracy = 53.3%. Of the 79 metabolites, 23 had VIP
   scores above 1.0; the top five were glutamic acid, myristoleic acid,
   α-ketoglutaric acid, 3-hydroxypropionic acid, and tyrosine.

---

## Key Result

**Six metabolites passed BH-FDR correction (q < 0.05) across all 79 tested**, matching
the paper's reported set exactly:

| Metabolite | KW p | FDR q | NASH / Control ratio | Direction |
|---|---|---|---|---|
| Glutamic acid | < 0.001 | < 0.001 | 2.4× | ↑ Control < NAFL < NASH |
| α-Ketoglutaric acid | < 0.001 | 0.004 | 2.2× | ↑ Control < NAFL < NASH |
| Myristoleic acid | < 0.001 | 0.004 | 2.5× | ↑ Control < NAFL < NASH |
| Tyrosine | 0.001 | 0.015 | 1.6× | ↑ Control < NAFL < NASH |
| Kynurenic acid | 0.001 | 0.015 | 1.7× | ↑ Control < NAFL (plateau in NASH) |
| Palmitoleic acid | 0.004 | 0.048 | 2.2× | ↑ Control < NAFL < NASH |

All six are elevated with disease progression; none decrease. Five show a strictly
monotone increase from Control to NAFL to NASH (Spearman ρ = +1.0); kynurenic acid
rises Control to NAFL then plateaus (ρ = +0.5).

**Biological note:** Glutamic acid and α-ketoglutaric acid are not independent signals.
Glutamate dehydrogenase converts glutamate directly into α-ketoglutarate, which enters
the TCA cycle. Both rising together across disease stages reflects increased flux through
this single enzymatic step — one connected metabolic event, not two separate findings.

---

## Connection to the Broader Project

This analysis provides a second, independent line of evidence for the same
NAFL → NASH disease-progression structure identified in the transcriptomic arms of
this project. The bulk RNA-seq and scRNA-seq modules found TREM2, SPP1, and COL1A1 as
robust gene-level markers of NAFLD severity; this module identifies the corresponding
metabolite-level markers in plasma from the same disease spectrum.

---

## Limitations

This analysis uses the publicly released summary statistics table (group means, SDs,
pre-computed p-values), not individual patient measurements. As a result:

- The six significant metabolites and their FDR q-values are validated against the paper's
  own reported statistics, but were not independently computed from raw data.
- Individual-level classification metrics — sensitivity, specificity, ROC curves, and the
  53.3% PLS-DA accuracy — cannot be independently reproduced or verified from this file.
- The heatmap and trend analysis are computed from group means, not from 86 individual
  Z-scores as in the paper's Figure 1A.
- Shapiro-Wilk normality testing, formal Jonckheere-Terpstra trend tests, and independent
  PLS-DA computation would require per-patient data, which would need to be requested from
  the corresponding author.

---

## File Guide

```
data/
    Table_S1_Metabolite_Levels.csv         Group means, SDs, and pre-computed statistics
                                            for all 79 metabolites (paper's Table S1)

scripts/
    ji2022_metabolomics_analysis.py        Full analysis pipeline: Steps 1–8 with console
                                            output documenting QC, statistics, and what
                                            each step found

plots/
    significant_metabolites_bar.png        Mean ± SD bar plots for the 6 FDR-significant
                                            metabolites across the three groups
    heatmap_79metabolites_zscore.png       Pseudo-Z-score heatmap, all 79 metabolites,
                                            colour-coded by metabolite category
    volcano_NASH_vs_Control.png            log₂FC vs −log₁₀(FDR q-value) for all 79

results/
    significant_metabolites_FDR05.csv      The 6 FDR-significant metabolites with statistics
    all_metabolites_statistics.csv         Full table: KW p-values, FDR (paper + recomputed),
                                            Spearman trend ρ, log₂FC
    all_metabolites_ranked.csv             All 79 metabolites ranked by FDR q-value
    trend_direction_proxy.csv              Spearman ρ for all 79 metabolites
    sig_metabolites_for_metaboanalyst.txt  6 metabolite names for MetaboAnalyst upload
```

---

## Tools

- **Python 3** — pandas, scipy, statsmodels, matplotlib, seaborn
- **MetaboAnalyst** (external) — pathway enrichment; input file provided in `results/`
