# Statistical Decision Rules

- Overall survival is the primary endpoint.
- Discovery modeling uses GSE14520 and a ferroptosis-reference plus tumor/non-tumor differential-expression candidate set.
- The final risk score is a fixed-coefficient nine-gene LASSO-Cox model.
- External cohorts are evaluated with fixed model coefficients after gene-symbol harmonization, including PKM2-to-PKM harmonization where needed.
- TCGA-LIHC, GSE76427, and ICGC-LIRI-JP/HCCDB18 have complete nine-gene coverage; GSE10143-HCC and GSE27150 are lower-coverage sensitivity cohorts.
- Risk groups are median-split for visualization; continuous risk-score analyses are used for model effect estimates.
- Time-dependent AUC is reported at 1, 3, and 5 years when follow-up and event counts support stable estimation.
- Calibration and Brier score summaries are reported as model-performance diagnostics.
- Decision curve analysis uses fitted Cox baseline-hazard predictions and landmark outcomes, not hardcoded baseline-risk constants.
- Immune analyses are marker-gene signature score analyses and checkpoint-expression summaries; they are exploratory and bulk-transcriptome based.
- TCGA multi-omics analyses are exploratory associations with mutation, copy-number, and subtype features.
