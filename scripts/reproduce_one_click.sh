#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
bash scripts/verify_release.sh
bash scripts/rebuild_results_and_figures.sh
