# HCC Ferroptosis Prognostic Signature

This repository provides the public reproducibility package for a hepatocellular carcinoma (HCC) ferroptosis-related prognostic signature study. The package contains the R scripts, derived result tables, minimal figure-rebuild inputs, and exported figures needed to reproduce the current analysis tables and figure set.

## Study scope

The current analysis develops and evaluates a nine-gene ferroptosis-related risk score for overall survival stratification in HCC. It includes discovery analysis in GSE14520, external evaluation in TCGA-LIHC, ICGC-LIRI-JP/HCCDB18, and GEO cohorts, calibration-aware model diagnostics, marker-gene immune signature analyses, and TCGA mutation/copy-number/subtype characterization.

This repository does not include manuscript submission files, cover letters, internal planning files, or retired modules.

## Repository structure

- `scripts_final/`: R analysis scripts used for the current pipeline.
- `scripts/`: shell entrypoints for reproduction and verification.
- `data/processed/`: minimal processed public-data objects required to rebuild the current publication figures.
- `data/references/`: ferroptosis and immune reference files.
- `results/`: derived tables supporting the analysis and figures.
- `plots/publication/`: final main figures.
- `plots/supplementary/`: supplementary figures and diagnostics.
- `docs/`: data manifest, figure provenance, and statistical decision rules.

## Quick reproduction

Install R and the required R packages listed in `environment.yml`, then run:

```bash
bash scripts/verify_release.sh
bash scripts/rebuild_figures_only.sh
```

To run the one-click public reproduction path, which verifies the release and rebuilds the current figure set from included derived outputs:

```bash
bash scripts/reproduce_one_click.sh
```

The one-click public path rewrites files in `plots/publication/` and `plots/supplementary/`. Full end-to-end reanalysis from raw public repositories can be performed with the scripts in `scripts_final/` after downloading the public datasets listed in `docs/DATA_MANIFEST.tsv`.

## Runtime expectation

Rebuilding figures from existing outputs usually takes minutes on a laptop. A full pipeline run can take substantially longer because it recomputes survival modeling, external validation, diagnostics, immune marker scores, and TCGA multi-omics summaries.

## Data availability

All datasets used here are public. See `docs/DATA_MANIFEST.tsv` for accession/source information and intended analysis role.

## Citation

If using this repository, cite the corresponding GitHub release and Zenodo DOI for the exact archived version.
