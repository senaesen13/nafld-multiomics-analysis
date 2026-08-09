#!/usr/bin/env Rscript
# Download E-MTAB-12807 - Grandt et al. 2025, healthy/NAFLD/cirrhosis liver
# bulk RNA-seq (ArrayExpress/BioStudies, not GEO).
#
# This series is deposited as per-sample Salmon quant.sf files (transcript-
# level quantification) rather than a ready-made gene-count matrix, so this
# script also aggregates transcripts to genes with tximport before saving.
# Run from the project root: Rscript code/00_download/download_EMTAB12807.R

suppressPackageStartupMessages({
  library(tximport)
  library(biomaRt)
  library(dplyr)
})
source(here::here("code", "utils", "paths.R"))

out_dir   <- dataset_path("RNAseq", "EMTAB12807")
quant_dir <- file.path(out_dir, "quant_files")
dir.create(quant_dir, recursive = TRUE, showWarnings = FALSE)

base_url <- "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/807/E-MTAB-12807/Files"

# --- sample metadata (SDRF) -------------------------------------------------
sdrf_path <- file.path(out_dir, "E-MTAB-12807.sdrf.txt")
if (!file.exists(sdrf_path)) {
  download.file(file.path(base_url, "E-MTAB-12807.sdrf.txt"), sdrf_path, quiet = FALSE)
}
sdrf <- read.delim(sdrf_path, check.names = FALSE, stringsAsFactors = FALSE)

# The SDRF has one row per FASTQ read (R1 and R2), so two rows per sample.
# One row per sample is enough to know its condition and which quant.sf file
# belongs to it.
meta <- sdrf |>
  distinct(`Source Name`, .keep_all = TRUE) |>
  transmute(
    sample_id  = `Source Name`,
    quant_file = `Derived Array Data File`,
    disease    = `Characteristics[disease]`,
    timing     = `Factor Value[clinical information]`,
    age        = `Characteristics[age]`,
    sex        = `Characteristics[sex]`
  )
write.csv(meta, file.path(out_dir, "meta_emtab12807.csv"), row.names = FALSE)
cat("Disease breakdown:\n")
print(table(meta$disease))
cat("Timing breakdown:\n")
print(table(meta$timing))

# --- per-sample Salmon quant.sf files ---------------------------------------
for (i in seq_len(nrow(meta))) {
  dest <- file.path(quant_dir, meta$quant_file[i])
  if (!file.exists(dest)) {
    download.file(file.path(base_url, meta$quant_file[i]), dest, mode = "wb", quiet = TRUE)
  }
}
cat("Downloaded", length(list.files(quant_dir)), "quant.sf files\n")

# --- transcript -> gene map --------------------------------------------------
# quant.sf transcript IDs carry a version suffix (ENST00000456328.2) that
# needs stripping before mapping. org.Hs.eg.db's ENSEMBLTRANS table was
# tried first and resolves a gene symbol for only ~10% of these transcripts
# (it is curated conservatively); querying Ensembl's own BioMart directly
# resolves ~87%, so that is what is used here. This step is slow (a few
# minutes, one BioMart query for ~230k transcripts) and only needs to run
# once - the result is cached to disk.
tx2gene_path <- file.path(out_dir, "tx2gene.csv")
one_quant    <- read.delim(file.path(quant_dir, meta$quant_file[1]))
tx_ids       <- sub("\\.[0-9]+$", "", one_quant$Name)

if (!file.exists(tx2gene_path)) {
  mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  bm <- getBM(
    attributes = c("ensembl_transcript_id", "hgnc_symbol"),
    filters = "ensembl_transcript_id", values = unique(tx_ids), mart = mart
  ) |> filter(hgnc_symbol != "")
  write.csv(bm, tx2gene_path, row.names = FALSE)
}
tx2gene_map <- read.csv(tx2gene_path)

# tximport needs a tx2gene table indexed by the *exact* transcript ID as it
# appears in the quant.sf files (i.e. with the version suffix restored).
tx2gene_full <- data.frame(tx_id = one_quant$Name, ensembl_transcript_id = tx_ids) |>
  inner_join(tx2gene_map, by = "ensembl_transcript_id") |>
  transmute(tx_id, gene_symbol = hgnc_symbol)

# --- aggregate to gene level with tximport -----------------------------------
quant_paths <- file.path(quant_dir, meta$quant_file)
names(quant_paths) <- meta$sample_id
txi <- tximport(quant_paths, type = "salmon", tx2gene = tx2gene_full,
                 countsFromAbundance = "lengthScaledTPM")

write.csv(txi$counts, file.path(out_dir, "EMTAB12807_gene_counts.csv"))
saveRDS(txi, file.path(out_dir, "EMTAB12807_txi.rds"))

cat(sprintf("EMTAB12807: %d genes x %d samples, matrix at %s\n",
            nrow(txi$counts), ncol(txi$counts),
            file.path(out_dir, "EMTAB12807_gene_counts.csv")))
