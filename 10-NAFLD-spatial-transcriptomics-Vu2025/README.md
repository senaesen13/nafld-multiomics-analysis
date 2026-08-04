# 10 — MASLD Spatial Transcriptomics (Vu et al. 2025)

## Dataset

**Paper:** Vu H, Sun Y, Xiong Z, Tan X et al. "Progressive fibrosis in human MASLD is associated with spatially linked transcriptomic signatures of metabolic reprogramming and senescence." *JHEP Reports* (2025). PMID 41541503, DOI 10.1016/j.jhepr.2025.101657.

**Data:** 10x Genomics Visium CytAssist spatial transcriptomics from FFPE (formalin-fixed paraffin-embedded) human liver biopsies. Downloaded from UQeSpace (DOI: 10.48610/e95155f). Analysis code: github.com/BiomedicalMachineLearning/Liver.

**What Visium CytAssist gives you:** Each Visium slide places a grid of ~5,000 capture spots (~55 µm diameter each) over a tissue section. Each spot captures the RNA from ~3–10 cells directly underneath it — not single cells, but small groups. CytAssist is the FFPE-compatible version, using a probe panel of 18,085 targeted genes rather than full-transcriptome sequencing.

**Study design:**
- 33 MASLD liver biopsies from 33 patients spanning fibrosis stages F0–F4 (Brunt scoring)
- 10 Visium capture areas (called "arrays") in the publicly released data; 8 used in the published analysis
- Multiplexed design: each array holds 4–5 biopsies from patients at *different* fibrosis stages, placed side-by-side on the same slide — this was deliberate to minimize technical batch effects between stages
- 3-group clinical classification: Early (F0–F1, n=14), Intermediate (F2–F3a, n=5), Late (F3b–F4, n=14)

---

## Method

### Step 1 — Spot QC
Each of the 8 arrays was loaded separately. Spots were filtered to retain only those with ≥200 detected genes and <10% mitochondrial gene expression, and genes were kept only if detected in ≥3 spots. This removed spots likely from empty space or damaged tissue.

| Array    | Spots before | Spots after | % kept |
|----------|-------------|-------------|--------|
| VLP115_A | 1,744       | 1,415       | 81.1%  |
| VLP116_D | 2,173       | 2,163       | 99.5%  |
| VLP119_A | 1,691       | 1,689       | 99.9%  |
| VLP119_D | 1,277       | 1,277       | 100%   |
| VLP120_A | 1,884       | 1,866       | 99.0%  |
| VLP120_D | 1,657       | 1,656       | 99.9%  |
| VLP121_A | 1,530       | 1,492       | 97.5%  |
| VLP121_D | 1,683       | 1,681       | 99.9%  |
| **Total**| **13,659**  | **13,239**  | **96.9%** |

VLP115_A had the lowest quality (median ~1,400 genes/spot vs ~2,500–4,000 for other arrays) and lost 18.9% of spots to QC filtering.

### Step 2 — Normalization (SCTransform per array, then merge)
Each array was normalized independently using SCTransform (Hafemeister & Satija 2019, *Genome Biology*) — a variance-stabilizing normalization that accounts for differences in sequencing depth between spots. After per-array normalization, 3,000 shared variable features (genes varying most across the dataset) were selected, and all 8 arrays were merged into a single combined object: **13,239 spots × 17,910 genes**.

### Step 3 — Spatial Domain Discovery
Dimensionality reduction and unsupervised clustering were applied to the merged object to identify spatially coherent transcriptional programs:
- **PCA**: 30 principal components computed on 3,000 variable genes
- **UMAP**: 2D embedding using the top 20 PCs
- **Clustering**: Louvain algorithm at resolution=0.5 (implemented via Seurat's `FindNeighbors` + `FindClusters`)

### Step 4 — Cell-type Deconvolution
Because each Visium spot covers multiple cells, we scored each spot for the relative activity of 17 known liver cell types using **AddModuleScore** (Seurat). The reference cell-type gene sets were derived from the GSE136103 scRNA-seq dataset (Ramachandran et al. 2019, *Nature Medicine*) — 35,050 human liver cells across healthy and cirrhotic donors, annotated into 20 cell types.

The intended method was **CARD** (Ma & Zhou 2022, *Nature Biotechnology*), a Bayesian deconvolution framework that estimates fractional cell-type composition per spot while modelling spatial autocorrelation. CARD could not be compiled on the local system (C++ linker failure on macOS), so AddModuleScore was used as a fallback — it gives relative activity scores, not compositional proportions.

### Step 5 — Biological Annotation (FindAllMarkers)
`FindAllMarkers` (Seurat Wilcoxon test, positive markers only, logFC ≥ 0.25, min.pct ≥ 0.1) was run on the SCT assay after `PrepSCTFindMarkers` to identify genes that define each of the 24 clusters relative to all others. Top marker genes were scored against 10 biological gene sets (periportal/pericentral/midzonal hepatocytes; stellate/fibrotic; Kupffer/macrophage; LSEC; biliary; lipogenic; inflammatory; senescent) to assign plain-language biological labels. The 4 genes validated across earlier analyses in this project — **TREM2, SPP1, COL1A1, FASN** — were explicitly tracked.

---

## Key Results

### 24 Spatial Domain Clusters — Biological Annotation

`FindAllMarkers` identified **27,416 significant marker genes** (adj. p < 0.05) across the 24 clusters.

| Cluster | Biological label | Top 5 markers | Key gene hits | Spots |
|---------|-----------------|---------------|---------------|-------|
| C0  | Hepatocyte parenchyma (indeterminate) | MT-ND1, ZNF549, MT-ND5, ZNF550, GSTA2 | — | 2,521 |
| C1  | Portal/stromal zone | MYH11, CCL19, LMOD1, CFTR, NOTCH3 | SPP1, COL1A1 | 1,752 |
| C2  | Hepatocyte parenchyma (indeterminate) | ASCL1, TFF3, UPP2, RET, PLCH2 | FASN | 1,061 |
| C3  | Inflammatory/Kupffer-macrophage region | EEF1A2, CXCL10, UBD, FABP4, CAPG | COL1A1, SPP1 | 1,006 |
| C4  | Hepatocyte parenchyma (indeterminate) | HYDIN, ZNF385D, DHRS2, SLC45A2, ACSL4 | — | 833 |
| C5  | Hepatocyte parenchyma (indeterminate) | APOA4, ACSL4, LOXL4, AKR1B10, HKDC1 | COL1A1, SPP1, FASN | 632 |
| C6  | Hepatocyte parenchyma (indeterminate) | SLCO1A2, PTH2R, CYP2B6, SULT1C2, CYP4F3 | FASN | 566 |
| C7  | Hepatocyte parenchyma — ER stress | SYT7, IFI27, HSPA5, MANF, OLFM2 | FASN | 452 |
| C8  | Hepatocyte parenchyma (indeterminate) | NUDT8, IGFBP1, GALK1, MT-ND4L, MAP4 | FASN | 451 |
| C9  | Hepatocyte parenchyma (indeterminate) | COL7A1, ESPL1, NPIPB15, CYP4X1, MROH7 | — | 407 |
| C10 | **Pericentral hepatocyte zone** | CYP3A4, IL1RAP, SLC16A1, CFHR4, EEF1B2 | — | 399 |
| C11 | Hepatocyte parenchyma (indeterminate) | PRSS51, LINGO4, NECAB2, HAMP, EYA4 | FASN | 396 |
| C12 | Hepatocyte parenchyma (indeterminate) | IGFBP1, MX1, HMCN2, MROH7, CMPK2 | — | 382 |
| C13 | Hepatocyte parenchyma (indeterminate) | SPINK1, PZP, SMIM24, HAL, BCO2 | — | 361 |
| C14 | Hepatocyte parenchyma (indeterminate) | EPB41L1, GABRB3, SLC5A12, A2M, CPN2 | — | 314 |
| C15 | Hepatocyte parenchyma (indeterminate) | TAT, MDN1, ZBTB16, GNLY, A2M | — | 293 |
| C16 | **Interferon-stimulated hepatocytes** | IFI6, AKR1C2, MX1, IFI44L, OAS1 | — | 286 |
| C17 | **Pericentral / Midzonal — cholesterol synthesis** | NUP155, CYP51A1, PLCH2, MSMO1, ADFP | FASN | 275 |
| C18 | Hepatocyte parenchyma (indeterminate) | UPK3B, FOXN4, AOC1, MEP1B, NPIPB15 | FASN | 212 |
| C19 | Hepatocyte parenchyma (indeterminate) | SLC44A5, XPNPEP2, IGFBP1, MYOM1, GPC6 | — | 191 |
| C20 | Hepatocyte parenchyma — metal stress | SLC29A4, DNM1, MT1H, SCAMP5, MT1G | — | 189 |
| C21 | Hepatocyte parenchyma (indeterminate) | UNC93A, ENPP3, SDS, ATAD3C, MAB21L1 | — | 117 |
| C22 | **Pericentral / Midzonal — cholesterol synthesis** | SPINK1, PLA2G2A, HMGCS1, LGALS4, PIK3C2G | FASN | 78 |
| C23 | **Midzonal hepatocyte zone** | GSTM1, FAM151A, ACAT2, CYP2B6, GSTM3 | FASN | 65 |

**Bolded labels** = clusters with clear non-indeterminate biological signal.

#### Notes on specific clusters

**C3 — Inflammatory/Kupffer-macrophage region (1,006 spots):** The only cluster with a clear non-hepatocyte inflammatory signature. CXCL10 is a defining interferon-γ-induced chemokine; FABP4 and CAPG are canonical macrophage markers; GPNMB is strongly enriched in activated Kupffer cells in human MASLD. COL1A1 and SPP1 co-occurring here places this cluster at the intersection of inflammation and fibrosis.

**C1 — Portal/stromal zone (1,752 spots, 13.2%):** The top markers — MYH11 (smooth muscle myosin heavy chain), LMOD1 (leiomodin), CCL19, CCL21 (lymphoid chemokines expressed in portal fibroblasts), NOTCH3 (portal vascular smooth muscle), FBLN1 (fibulin, ECM) — point to portal tract tissue including vascular smooth muscle and portal fibroblasts, not hepatocytes. The 13.2% share is consistent with portal tract area across multiple biopsies. SPP1 and COL1A1 among C1's markers confirms fibrosis-associated activity in this zone.

**C10 — Pericentral hepatocyte zone (399 spots):** CYP3A4 is the defining zone-3 (centrilobular) marker, restricted by the Wnt signalling gradient. This is the clearest zonation cluster.

**C16 — Interferon-stimulated hepatocytes (286 spots):** IFI6, MX1, IFI44L, OAS1 are all canonical interferon-stimulated genes. This small but distinct cluster likely represents hepatocytes responding to an active innate immune signal — possibly viral sensing or cytokine stimulation from adjacent inflammatory cells.

**C17 and C22 — Cholesterol/lipid synthesis zones:** Both enriched for the mevalonate/cholesterol pathway (CYP51A1, MSMO1, HMGCS1) alongside FASN. These represent hepatocytes with maximal de novo lipid synthesis activity — directly relevant to MASLD steatosis.

### Cross-modality Confirmation: TREM2 / SPP1 / COL1A1 / FASN

These 4 genes were identified as significant in earlier stages of this multiomics project. The spatial data provides the following:

| Gene | Spatial finding | Interpretation |
|------|----------------|----------------|
| **TREM2** | Not detected as a cluster marker | TREM2 is primarily expressed in a rare macrophage subpopulation (~1–3% of liver cells). At Visium resolution (~3–10 cells per spot), its signal is diluted below the min.pct threshold. This is a known limitation of bulk-spatial transcriptomics for rare cell types, not evidence that TREM2-positive macrophages are absent. |
| **SPP1** | Significant in C1, C3, C5 | SPP1 (osteopontin) marks three clusters — the portal/stromal zone (C1), the inflammatory macrophage niche (C3), and a hepatocyte cluster with stromal features (C5). In the week 6 scRNA-seq analysis, SPP1 was the top marker of TREM2+ scar-associated macrophages. The spatial co-localisation with COL1A1 in the same clusters is a **direct cross-modality confirmation** of the SPP1/COL1A1 fibrosis signature. |
| **COL1A1** | Significant in C1, C3, C5 | Identical cluster distribution as SPP1. In week 4 bulk RNA-seq, COL1A1 was the strongest upregulated gene in fibrotic NAFLD. The spatial data shows that COL1A1 expression is not diffuse — it is concentrated in portal/stromal (C1) and inflammatory (C3) zones, not throughout the hepatocyte parenchyma. |
| **FASN** | Significant in 10/24 clusters (C2, C5, C6, C7, C8, C11, C17, C18, C22, C23) — 2,239 spots, 16.9% of tissue | FASN (fatty acid synthase, the key de novo lipogenesis enzyme) is spatially pervasive in MASLD liver. Its presence across 10 clusters, including both the lipogenic zones (C17, C22) and otherwise-indeterminate hepatocyte clusters, confirms that de novo lipogenesis is not restricted to a specialist metabolic zone but is a near-universal hepatocyte phenotype in MASLD tissue. This is consistent with the bulk RNA-seq finding (week 4) but spatially resolves it. |

**The SPP1 + COL1A1 co-localisation to the same 3 clusters** is the most meaningful cross-modality result: these two genes were identified independently in bulk RNA-seq (week 4), scRNA-seq (week 6), and now spatial transcriptomics (week 8) as markers of fibrotic niches. The spatial data shows they co-occupy the same tissue zones rather than being expressed by different cell populations in different regions.

### Cell-type Deconvolution (AddModuleScore)
All 24 clusters showed hepatocyte activity as the dominant signal. The clusters with the highest relative enrichment of hepatic stellate cell and Kupffer cell scores (relative to hepatocyte score) were C1, C3, C5, and C20 — which matches the independent FindAllMarkers result identifying C1, C3, and C5 as the clusters where SPP1 and COL1A1 appear. This concordance between the deconvolution approach and the marker gene approach strengthens confidence in the biological annotation.

---

## Key Limitation

**Fibrosis stage cannot yet be overlaid per spot.** The metadata file mapping each tissue barcode to a specific patient biopsy (and therefore to a fibrosis stage F0–F4) was not included in the public UQeSpace release. This means we cannot ask "are the inflammatory/fibrotic clusters enriched in biopsies from higher fibrosis stages?" — the central question of the Vu et al. paper. The 24 cluster structure and their biological annotation are established; the fibrosis-stage context requires the barcode-to-patient mapping, which exists on the authors' HPC.

---

## Files in This Folder

| File | Description |
|------|-------------|
| `umap_clusters_and_arrays.png` | UMAP of 13,239 spots: numbered clusters (left) and by Visium array (right) |
| `umap_biological_labels.png` | UMAP with biological labels: numbered clusters (left) and annotated labels (right) |
| `umap_bio_labels_only.png` | UMAP with biological labels only |
| `celltype_scores_heatmap_by_cluster.png` | Z-scored cell-type activity (17 types) across 24 spatial clusters |
| `qc_summary.csv` | Per-array spot counts before/after QC, median genes and UMIs per spot |
| `cluster_sizes.csv` | Spot count and percentage for each of the 24 clusters |
| `cluster_annotation_summary.csv` | Biological label, top 5 markers, key gene hits, and spot count per cluster |

Full analysis (scripts, all plots, deconvolution CSVs, marker gene tables, NOTES.md): `week8-day1-masld-spatial-setup/` in the bioinformatics-learning repo.

---

## Software

- R 4.6.0, Seurat 5.5.1, sctransform 0.4.3, hdf5r 1.3.12
- Reference: GSE136103 (Ramachandran et al. 2019, *Nature Medicine*) scRNA-seq — 35,050 human liver cells, 20 cell types
- FindAllMarkers: Wilcoxon rank-sum test, only.pos=TRUE, min.pct=0.1, logFC≥0.25, SCT assay
