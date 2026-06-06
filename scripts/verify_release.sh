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
test -f results/prognostic_model_coef.csv
test -f plots/publication/Figure1_workflow.png
test -f plots/publication/Figure5_immune_analysis.png
echo 'Release verification OK'
