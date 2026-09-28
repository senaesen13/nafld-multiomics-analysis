# Archive

Older versions of analyses, kept for transparency and comparison. Nothing here is deleted.

Rule: when a module is improved, the previous outputs are moved here with git mv
(so history is kept) into a folder named <module>_<YYYY-MM>_<version>, and a note
below explains what changed and why.

## 09-drug-repositioning_2026-07_first-signature
- Moved on: 2026-09-28
- Query: LINCS CMap L1000, submitted 2026-07-30 as NAFLD_repositioning
- Input: genes from the GSE130970 DESeq2 results (padj < 0.05, |log2FC| > 0.5, top 150 by effect size); 150 up, 51 down submitted; 99 up / 80 down matched L1000 landmarks
- Why replaced: GSE130970 has only 4 controls and is no longer part of the main pipeline; the new query uses the cross-cohort consensus signature
- Known issues: the original README said 150 down genes but the .grp file has 51; the old plot script used simulated placeholder Tau values, so the old candidates table should be treated as preliminary
