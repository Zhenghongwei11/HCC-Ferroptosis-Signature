#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
bash scripts/verify_release.sh
bash scripts/rebuild_figures_only.sh
