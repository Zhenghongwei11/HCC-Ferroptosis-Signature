#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
Rscript scripts_final/07_final_publication_figures.R
Rscript scripts_final/07b_fig1_workflow.R
