#!/usr/bin/env Rscript
# Spatial stage 4: per-spot metabolic flux estimation with METAFlux.
#
# NOT YET RUN in this project - see the note at the top of
# 01_qc_and_normalize.R. Depends on stage 2's clustered/annotated object.
#
# Usage: Rscript code/09_spatial/04_metabolic_flux.R
#
# This is the mechanistic counterpart to stage 07's reporter metabolite
# analysis: reporter metabolites ask "is this metabolite's enzymatic
# neighbourhood transcriptionally hot" (a statistical/topological question,
# answered per bulk cohort); this stage asks "what is the predicted flux
# through the metabolic network at each spatial location" (a mechanistic
# question, using Human-GEM + FBA, resolved per spot). Running both and
# checking whether they agree - do the reporter metabolites flagged in
# stage 07 also show altered predicted flux here, and does that flux
# co-localise with the fibrotic/inflammatory spatial domains found in stage
# 02 - is a stronger claim than either analysis alone.
#
# Note: this is the actual implementation of "COMPASS/METAFlux per spatial
# cluster", which is listed only as a placeholder synthetic-data template in
# Workshops/Module_3_Spatial_Omics_Analysis/Code/04_spatial_metabolic_flux.R
# (that script simulates a fake 10x10 tumour grid - it is a structural
# example, not runnable on this project's real liver data). This script
# uses the real, installed METAFlux R package instead.
#
# API note: a development smoke-test of this exact call sequence
# (calculate_avg_exp -> calculate_reaction_score -> compute_sc_flux) on
# METAFlux's own bundled sc_test_example hit an internal package error
# ("invalid class 'LogMap' object: Duplicate rownames not allowed", inside
# calculate_avg_exp's bootstrap resampling step) - this looks like a
# Bioconductor S4Vectors version-compatibility issue in the installed
# METAFlux build, not a mistake in how the functions are called here (the
# function names and argument order below are confirmed correct via
# `args()` against the installed package). Resolve this - likely by
# checking S4Vectors/METAFlux version compatibility, or updating METAFlux -
# before this stage is first run.

suppressPackageStartupMessages({
  library(Seurat)
  library(METAFlux)
})
source(here::here("code", "utils", "paths.R"))

spatial_dir <- results_path("spatial", "Vu2025_Visium")
obj <- readRDS(file.path(spatial_dir, "02_clustered_annotated.rds"))
data(human_blood)  # bundled blood/plasma metabolite availability constraints -
                    # the standard medium assumption for any human tissue FBA,
                    # and particularly apt for liver given its central role in
                    # processing blood-borne metabolites

# calculate_avg_exp treats each level of the identity column as a
# pseudo-cell-type and bootstrap-resamples within it - here that identity is
# spatial cluster (from stage 02), so "cell type" in METAFlux's own
# terminology below actually means "spatial domain" for this analysis.
avg_exp    <- calculate_avg_exp(obj, "seurat_clusters", n_bootstrap = 100, seed = 42)
rxn_scores <- calculate_reaction_score(avg_exp)

cluster_fraction <- table(obj$seurat_clusters) / ncol(obj)
flux <- compute_sc_flux(
  num_cell  = ncol(obj),
  fraction  = cluster_fraction,
  fluxscore = rxn_scores,
  medium    = human_blood
)

write.csv(as.data.frame(flux), file.path(spatial_dir, "04_metabolic_flux_by_cluster.csv"))
cat(sprintf("Computed flux for %d reactions across %d spatial domains\n", nrow(flux), ncol(flux)))
