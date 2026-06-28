#!/usr/bin/env Rscript

# run_complete_pipeline.R
# Public reproducibility pipeline for the HCC ferroptosis signature analysis.
# Usage: Rscript scripts_final/run_complete_pipeline.R

options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

message("================================================================================")
message("HCC ferroptosis signature reproducibility pipeline")
message("================================================================================")
message("Start time: ", Sys.time())

if (!dir.exists("data/processed")) {
  stop("Run this script from the repository root.")
}

message("\n[Stage 1] Differential expression and HCC-context ferroptosis set...")
source("scripts_final/00c_prepare_hcc_ferro_genes.R")
source("scripts_final/02b_multi_cohort_DEG.R")

message("\n[Stage 2] Prognostic modeling and validation...")
source("scripts_final/02c_prognostic_model.R")
source("scripts_final/02d_nomogram_calibration_dca.R")
source("scripts_final/02e_external_validation.R")
source("scripts_final/02h_model_diagnostics.R")
source("scripts_final/02i_external_calibration_brier.R")
source("scripts_final/02j_incremental_value_tcga.R")
source("scripts_final/02k_spline_risk_effect.R")
source("scripts_final/02m_tcga_multiomics_characterization.R")

message("\n[Stage 3] Immune-marker summaries...")
source("scripts_final/03a_immune_infiltration.R")
source("scripts_final/03b_immune_checkpoint.R")

message("\n[Stage 4] Additional validation and benchmark summaries...")
source("scripts_final/02n_model_validation_sensitivity.R")
source("scripts_final/02o_published_signature_benchmark.R")

message("\n[Stage 5] Figures...")
source("scripts_final/07_final_publication_figures.R")
source("scripts_final/07b_fig1_workflow.R")

message("\n================================================================================")
message("Reproducibility pipeline completed.")
message("Results are in 'results/'. Figures are in 'plots/'.")
message("================================================================================")
