#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

Rscript scripts_final/run_complete_pipeline.R
