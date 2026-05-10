#!/usr/bin/env bash
set -euo pipefail

# Membrane-pore reference run:
# - transverse y boundary: periodic unit cell
# - streamwise x boundary: pressure inlet/outlet
# - fixed cylinders create the pore-scale pressure loss

poetry run python src/demos/run_lbm_dem.py \
  --nx 100 \
  --ny 50 \
  --total-steps 3000 \
  --snapshot-every 300 \
  --warmup-steps 100 \
  --y-boundary periodic \
  --streamwise-boundary pressure \
  --pressure-drop 1e-5 \
  --rho-out 1.0 \
  --reynolds-number 20 \
  --particle-source left-inlet \
  --n-particles 120 \
  --particle-radius 1.0 \
  --radius-variation 0.05 \
  --cylinder-spec 34 13 4.0 \
  --cylinder-spec 34 37 4.0 \
  --cylinder-spec 54 25 4.5 \
  --cylinder-spec 74 13 3.8 \
  --cylinder-spec 74 37 3.8 \
  --particle-fluid-coupling point_force \
  --fluid-accelerator auto \
  --compute-accelerator auto \
  --snapshot-storage none \
  --paraview-output \
  --paraview-every 0 \
  --max-stable-speed 0.5 \
  --max-stable-pressure 2.0 \
  --result-tag membrane_pressure_periodic_reference \
  "$@"
