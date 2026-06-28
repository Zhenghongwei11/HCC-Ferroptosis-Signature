#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
Rscript - <<'RS'
files <- list.files('scripts_final', pattern='[.]R$', full.names=TRUE)
for (f in files) parse(f)
cat('R parse OK:', length(files), 'files\n')
RS
test -f docs/DATA_MANIFEST.tsv
test -f docs/FIGURE_PROVENANCE.tsv
test -f docs/TABLE_PROVENANCE.tsv
test -f scripts/rebuild_results_and_figures.sh
test -f data/processed/GSE14520_expr.rds
test -f data/processed/TCGA_LIHC_expr.rds
test -f data/processed/HCCDB18_expr_symbol_hcc.rds
test -f results/prognostic_model_coef.csv
test -f results/published_signature_benchmark.csv
test -f plots/publication/Figure1_workflow.png
test -f plots/publication/Figure5_immune_analysis.png
echo 'Release verification OK'
