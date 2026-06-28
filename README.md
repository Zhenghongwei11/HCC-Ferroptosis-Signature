# HCC Ferroptosis Prognostic Signature

This repository provides the public reproducibility package for a hepatocellular carcinoma (HCC) ferroptosis-related prognostic signature study. The package contains the R scripts, processed public-data objects, reference files, result tables, and exported figures needed to rebuild the current analysis tables and figure set.

## Study scope

The current analysis develops and evaluates a nine-gene ferroptosis-related risk score for overall survival stratification in HCC. It includes discovery analysis in GSE14520, external evaluation in TCGA-LIHC, ICGC-LIRI-JP/HCCDB18, and GEO cohorts, calibration and prediction-error diagnostics, published-signature benchmarking, proliferation-baseline analyses, marker-gene immune signature analyses, and TCGA mutation/copy-number/subtype characterization.

This repository does not include manuscript submission files, cover letters, internal planning files, or retired modules.

## Repository structure

- `scripts_final/`: R analysis scripts used for the current pipeline.
- `scripts/`: shell entrypoints for reproduction and verification.
- `data/processed/`: processed public-data objects required to rebuild the current result tables and figures.
- `data/references/`: ferroptosis and immune reference files.
- `results/`: derived tables supporting the analysis and figures.
- `plots/publication/`: final main figures.
- `plots/supplementary/`: supplementary figures and diagnostics.
- `docs/`: data manifest, figure provenance, and statistical decision rules.

## Quick reproduction

Install R and the required R packages listed in `environment.yml`, then run:

```bash
bash scripts/verify_release.sh
bash scripts/rebuild_results_and_figures.sh
```

To run the one-click public reproduction path, which verifies the release and rebuilds the current result tables and figure set from the included processed public-data objects:

```bash
bash scripts/reproduce_one_click.sh
```

The one-click public path rewrites files in `results/`, `plots/publication/`, and `plots/supplementary/`. Source accessions and local file roles are documented in `docs/DATA_MANIFEST.tsv`.

## Runtime expectation

The full public reproduction path can take substantially longer than a figure-only refresh because it recomputes survival modeling, external validation, diagnostics, immune marker scores, TCGA multi-omics summaries, validation sensitivity analyses, and published-signature benchmarks.

## Data availability

All datasets used here are public. See `docs/DATA_MANIFEST.tsv` for accession/source information and intended analysis role.

## Citation

If using this repository, cite the corresponding GitHub release and Zenodo DOI for the exact archived version.
