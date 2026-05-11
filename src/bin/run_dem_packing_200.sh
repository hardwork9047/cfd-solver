#!/usr/bin/env bash
set -euo pipefail

# DEM-only gravity packing example for 200 particles.
export PYTHONPATH="${PYTHONPATH:-}:src"
python3 src/runners/run_dem_packing.py \
  --config configs/dem_packing/cases/packing_200_particles.json \
  "$@"
