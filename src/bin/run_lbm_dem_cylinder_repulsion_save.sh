#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/../.."

export PYTHONPATH="${PWD}/src${PYTHONPATH:+:${PYTHONPATH}}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/matplotlib}"
mkdir -p "$MPLCONFIGDIR"

# Fixed-cylinder LBM-DEM run with Hamaker-like particle-particle repulsion.
# Default cylinder: x=45, y=35, radius=6.
# Results are saved under src/results/run_lbm_dem/cylinder_with_repulsion_rolling/.
python src/demos/run_lbm_dem.py --cylinder --particle-repulsion --rolling-friction "$@"
