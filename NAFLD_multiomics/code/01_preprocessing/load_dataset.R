# Dataset loading for the four bulk RNA-seq cohorts.
#
# Each GEO/ArrayExpress series deposits its data in a different shape (one
# assembled matrix vs. per-sample files; Ensembl IDs vs. gene symbols; a
# metadata column that lines up with the count matrix column names vs. one
# that doesn't). Rather than writing four near-duplicate analysis scripts,
# every difference is handled once, here, behind a single interface:
#   load_counts(dataset_id)   -> genes x samples count matrix
#   load_conditions(dataset_id) -> data.frame(sample_id, condition)
# Every later stage of the pipeline (02 onward) only ever calls these two
# functions, so it does not need to know which dataset it is looking at.

suppressPackageStartupMessages(library(dplyr))
source(here::here("code", "utils", "paths.R"))

load_counts <- function(dataset_id) {
  dir <- dataset_path("RNAseq", dataset_id)

  if (dataset_id == "GSE162694") {
    read.csv(file.path(dir, "GSE162694_raw_counts.csv.gz"), row.names = 1, check.names = FALSE) |> as.matrix()

  } else if (dataset_id == "GSE135251") {
    read.csv(file.path(dir, "GSE135251_count_matrix.csv"), row.names = 1, check.names = FALSE) |> as.matrix()

  } else if (dataset_id == "GSE126848") {
    m <- read.delim(file.path(dir, "GSE126848_Gene_counts_raw.txt.gz"), row.names = 1, check.names = FALSE) |> as.matrix()
    # Column names use inconsistent zero-padding ("0869" vs "2683"); the
    # metadata join key is normalised the same way in load_conditions(), so
    # both sides must match here too.
    colnames(m) <- sub("^0+", "", colnames(m))
    m

  } else if (dataset_id == "EMTAB12807") {
    read.csv(file.path(dir, "EMTAB12807_gene_counts.csv"), row.names = 1, check.names = FALSE) |> as.matrix()

  } else {
    stop("Unknown dataset_id: ", dataset_id)
  }
}

load_conditions <- function(dataset_id) {
  dir <- dataset_path("RNAseq", dataset_id)

  if (dataset_id == "GSE162694") {
    meta <- read.csv(file.path(dir, "meta_gse162694.csv"), row.names = 1, check.names = FALSE)
    # The GEO title field encodes the sample ID as its second space-separated
    # token (e.g. "nash1_F0 548nash1" -> "548nash1"), which is what the count
    # matrix columns are actually named - geo_accession (GSM...) is not.
    sample_id <- sub(".* ", "", meta$title)
    fibrosis  <- sub("fibrosis stage: ", "", meta$"characteristics_ch1.3")
    condition <- ifelse(fibrosis == "normal liver histology", "Normal", "NAFLD")
    data.frame(sample_id = sample_id, condition = condition)

  } else if (dataset_id == "GSE135251") {
    meta <- read.csv(file.path(dir, "meta_gse135251.csv"), row.names = 1, check.names = FALSE)
    condition <- ifelse(sub("disease: ", "", meta$"disease:ch1") == "Control", "Normal", "NAFLD")
    data.frame(sample_id = meta$geo_accession, condition = condition)

  } else if (dataset_id == "GSE126848") {
    meta <- read.csv(file.path(dir, "meta_gse126848.csv"), row.names = 1, check.names = FALSE)
    # The count matrix column names and the "description" metadata field
    # both hold the numeric sample ID, but with inconsistent zero-padding
    # (count column "0869" vs description "869") - strip leading zeros from
    # both sides before matching, or every "08xx"-numbered sample silently
    # fails to join.
    sample_id <- sub("^0+", "", trimws(as.character(meta$description)))
    disease   <- sub("disease: ", "", meta$"characteristics_ch1.2")
    # "healthy" and "obese" are both non-NAFLD control arms in this study
    # design (Suppli et al. 2019) - "obese" means metabolically-healthy
    # obese without steatosis, not a disease group.
    condition <- ifelse(disease %in% c("healthy", "obese"), "Normal", "NAFLD")
    data.frame(sample_id = sample_id, condition = condition)

  } else if (dataset_id == "EMTAB12807") {
    meta <- read.csv(file.path(dir, "meta_emtab12807.csv"), check.names = FALSE)
    meta <- meta[meta$disease != "excluded", ]
    condition <- case_when(
      meta$disease == "normal" ~ "Normal",
      meta$disease == "non-alcoholic fatty liver disease" ~ "NAFLD",
      meta$disease == "cirrhosis of liver" ~ "Cirrhosis",
      TRUE ~ NA_character_
    )
    data.frame(sample_id = meta$sample_id, condition = condition)

  } else {
    stop("Unknown dataset_id: ", dataset_id)
  }
}

#' Load counts + conditions for a dataset, already aligned to each other and
#' restricted to Normal/NAFLD (drops Cirrhosis or any sample missing a call).
load_aligned <- function(dataset_id, keep_conditions = c("Normal", "NAFLD")) {
  counts <- load_counts(dataset_id)
  conds  <- load_conditions(dataset_id) |> filter(.data$condition %in% keep_conditions)

  common <- intersect(colnames(counts), conds$sample_id)
  if (length(common) == 0) {
    stop("load_aligned(", dataset_id, "): no overlap between count matrix ",
         "columns and condition sample_id values - check the join key.")
  }

  conds <- conds[match(common, conds$sample_id), ]
  list(
    counts    = counts[, common, drop = FALSE],
    condition = factor(conds$condition, levels = keep_conditions)
  )
}
