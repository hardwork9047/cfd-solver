#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/../.."

export PYTHONPATH="${PWD}/src${PYTHONPATH:+:${PYTHONPATH}}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/matplotlib-cfd-solver}"
mkdir -p "$MPLCONFIGDIR"

TAG="four_cylinders_attr2x_phi40_muR3x"
LOG_DIR="src/results/run_lbm_dem/four_cylinder_cases"
LOG_FILE="${LOG_DIR}/${TAG}.log"
mkdir -p "$LOG_DIR"

cmd=(
  python -u src/demos/run_lbm_dem.py
  --cylinder-spec 75 25 6
  --cylinder-spec 75 45 6
  --cylinder-spec 105 25 6
  --cylinder-spec 105 45 6
  --particle-source left-inlet
  --particle-volume-fraction 0.40
  --particle-attraction
  --attraction-strength 0.002
  --rolling-friction
  --rolling-friction-coeff 0.15
  --total-steps 60000
  --snapshot-every 120
  --paraview-every 10
  --result-tag "$TAG"
)

echo "=== START ${TAG} ===" | tee "$LOG_FILE"
printf "Command:" | tee -a "$LOG_FILE"
printf " %q" "${cmd[@]}" | tee -a "$LOG_FILE"
printf "\n" | tee -a "$LOG_FILE"

"${cmd[@]}" >> "$LOG_FILE" 2>&1

echo "=== OK ${TAG} ===" | tee -a "$LOG_FILE"
echo "Log: ${LOG_FILE}"
