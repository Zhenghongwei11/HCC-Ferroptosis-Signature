# Compute Plan

The public release is optimized for reviewer-facing reproduction of the current derived tables and figure set.

- Recommended memory: 8 GB RAM or higher for figure rebuilding.
- Recommended disk: 1 GB free space for outputs and temporary files.
- Main public path: `bash scripts/reproduce_one_click.sh` verifies the package and rebuilds figures from included derived outputs and minimal figure-rebuild inputs.
- Full end-to-end reanalysis from raw repositories requires downloading public datasets listed in `docs/DATA_MANIFEST.tsv` and may require 16 GB RAM or higher.
