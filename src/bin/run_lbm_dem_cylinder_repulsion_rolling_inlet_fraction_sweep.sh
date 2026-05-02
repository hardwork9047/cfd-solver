#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/../.."

export PYTHONPATH="${PWD}/src${PYTHONPATH:+:${PYTHONPATH}}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/matplotlib}"
mkdir -p "$MPLCONFIGDIR"

# Filtration-style cylinder runs:
# particles are uniformly supplied from the left inlet, right-outflow particles
# are deleted, and results are saved under:
#   src/results/run_lbm_dem/cylinder_with_repulsion_rolling_left_inlet/phi_XXpct/
for phi in 0.05 0.10 0.20 0.40; do
  python -u src/demos/run_lbm_dem.py \
    --cylinder \
    --particle-repulsion \
    --rolling-friction \
    --particle-source left-inlet \
    --particle-volume-fraction "$phi" \
    "$@"
done
