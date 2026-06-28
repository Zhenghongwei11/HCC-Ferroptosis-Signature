# Compute Plan

The public release is designed to rebuild the current result tables and figure set from included processed public-data objects and reference files.

- Recommended memory: 8 GB RAM or higher.
- Recommended disk: 2 GB free space for outputs and temporary files.
- Main public path: `bash scripts/reproduce_one_click.sh` verifies the package and rebuilds result tables and figures from the included processed public-data objects.
- Equivalent explicit path: `bash scripts/verify_release.sh` followed by `bash scripts/rebuild_results_and_figures.sh`.
- Raw public source accessions are listed in `docs/DATA_MANIFEST.tsv`; the release includes processed public-data objects to make the archived reproduction path deterministic.
- Runtime depends on hardware and R package versions. On a laptop, the complete public path may take from tens of minutes to several hours because it recomputes modeling, validation, diagnostics, multi-omics summaries, sensitivity analyses, and benchmarking.
