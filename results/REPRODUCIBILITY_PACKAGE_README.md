# Derived Results

This directory contains derived result tables for the HCC ferroptosis prognostic-signature reproducibility package.

Key files include:

- `prognostic_model_coef.csv`: fixed nine-gene Cox model coefficients.
- `risk_score_data.csv`: discovery-cohort risk scores and survival data.
- `external_validation_stats.csv`: external validation performance summary.
- `external_validation_brier.csv`: external prediction-error summaries.
- `external_calibration_slope_intercept_3y.csv`, `external_dca_3y.csv`, and `fixed_threshold_net_benefit_3y.csv`: external 3-year calibration and decision-curve summaries.
- `lasso_selection_stability.csv`: bootstrap feature-selection stability.
- `tcga_proliferation_baseline_benchmark.csv`, `tcga_risk_adjusted_for_proliferation.csv`, and `tcga_risk_proliferation_correlation.csv`: TCGA proliferation-baseline and adjustment summaries.
- `published_signature_sources.csv`, `published_signature_gene_coverage.csv`, `published_signature_benchmark.csv`, and `benchmarks/published_signature_head_to_head.tsv`: published-signature and proliferation-baseline head-to-head benchmarking.
- `benchmarks/prediction_eval.tsv`, `benchmarks/method_benchmark.tsv`, and `effect_sizes/claim_effects.tsv`: compact reporting tables for prediction evaluation, method benchmarking, and claim-level effect sizes.
- `model_diagnostics_ph.csv`, `model_diagnostics_rmst.csv`, and `model_diagnostics_random_signature.csv`: model diagnostic summaries.
- `immune_risk_correlation.csv`, `immune_risk_group_diff.csv`, `checkpoint_risk_correlation.csv`, and `checkpoint_risk_diff.csv`: exploratory immune/checkpoint summaries.
- `immune_risk_group_diff_fdr.csv`, `checkpoint_risk_diff_fdr.csv`, and `tcga_multiomics_keygene_mutation_assoc_fdr.csv`: FDR-adjusted exploratory summaries.
- `dysfunction_exclusion_group_summary.csv` and `dysfunction_exclusion_surrogate.csv`: exploratory dysfunction/exclusion marker-score summaries.
- `tcga_multiomics_burden.csv`, `tcga_multiomics_keygene_mutation_assoc.csv`, and `tcga_subtype_assoc.csv`: TCGA multi-omics summaries.
