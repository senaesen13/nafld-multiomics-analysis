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

---

## Key Results

### 24 Spatial Domain Clusters
Unsupervised clustering of 13,239 spots identified **24 transcriptionally distinct spatial domains**. The 5 largest clusters account for ~54% of all spots:

| Cluster | Spots | Share |
|---------|-------|-------|
| C0      | 2,521 | 19.0% |
| C1      | 1,752 | 13.2% |
| C2      | 1,061 | 8.0%  |
| C3      | 1,006 | 7.6%  |
| C4      | 833   | 6.3%  |
| C5–C23  | 6,066 | 45.8% |

These clusters do not represent distinct cell types — in liver parenchyma, ~70–80% of cells (and ~90%+ of tissue volume) are hepatocytes. Instead, the clusters represent **different transcriptional states of the hepatocyte parenchyma**: metabolic zonation (periportal vs centrilobular), fibrosis-associated activation programs, inflammatory niche microenvironments, and biliary/portal vs lobular spatial zones. This is consistent with Vu et al.'s own findings, where spatial domains corresponded to distinct metabolic and cellular senescence signatures.

The UMAP plot (`umap_clusters_and_arrays.png`) shows both the 24-cluster structure and coloring by array — confirming that clusters mix spots from multiple arrays rather than reflecting technical batch effects.

### Cell-type Deconvolution
All 24 clusters show hepatocyte activity as the dominant signal (confirmed by AddModuleScore), which is **biologically expected** for liver tissue. Non-hepatocyte signals — Kupffer cells (resident macrophages), hepatic stellate cells (fibrosis drivers), LSEC (liver sinusoidal endothelial cells), Monocytes — are all present at lower levels and vary between clusters. The heatmap (`celltype_scores_heatmap_by_cluster.png`) shows the relative enrichment of 17 cell types across the 24 spatial clusters.

Three cell types from the reference (B cells, CD4+ T cells, Naive T cells) could not be scored because fewer than 5 of their top 50 marker genes are targeted by the Visium probe panel — a known limitation of probe-based spatial transcriptomics for immune cell detection.

---

## Key Limitation

**Fibrosis stage cannot yet be overlaid per spot.** The Visium data files (expression matrices + spatial images) are publicly available, but the metadata file that maps each tissue barcode to a specific patient biopsy — and therefore to a fibrosis stage (F0–F4) — was not included in the public UQeSpace release. This mapping file exists on the authors' HPC at the University of Queensland.

What we do have: patient-level fibrosis staging for all 33 biopsies (derived from the paper's GitHub analysis scripts — see `week8-day1-masld-spatial-setup/data/metadata/patient_fibrosis_staging.csv` in the bioinformatics-learning repo). What we cannot yet do: link individual spots on the slide to their patient of origin, and therefore to their fibrosis stage.

Without this mapping, we cannot ask "do spots in fibrotic tissue show a different transcriptional state than spots in non-fibrotic tissue?" — which is the central question of the Vu et al. paper. The 24 clusters provide spatial domain structure, but their biological interpretation in the context of fibrosis progression requires the barcode-to-patient mapping.

---

## Files in This Folder

| File | Description |
|------|-------------|
| `umap_clusters_and_arrays.png` | UMAP of 13,239 spots: left panel = 24 clusters, right panel = by Visium array |
| `celltype_scores_heatmap_by_cluster.png` | Z-scored cell-type activity (17 types) across 24 spatial clusters |
| `qc_summary.csv` | Per-array spot counts before/after QC, median genes and UMIs per spot |
| `cluster_sizes.csv` | Spot count and percentage for each of the 24 clusters |

Full analysis (scripts, all plots, deconvolution CSVs, NOTES.md): `week8-day1-masld-spatial-setup/` in the bioinformatics-learning repo.

---

## Software

- R 4.6.0, Seurat 5.5.1, sctransform 0.4.3, hdf5r 1.3.12
- Reference: GSE136103 (Ramachandran et al. 2019, *Nature Medicine*) scRNA-seq — 35,050 human liver cells, 20 cell types
