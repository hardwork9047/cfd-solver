#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/../.."

export PYTHONPATH="${PWD}/src${PYTHONPATH:+:${PYTHONPATH}}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/matplotlib-cfd-solver}"
mkdir -p "$MPLCONFIGDIR"

usage() {
  cat <<'USAGE'
Usage:
  src/bin/run_lbm_dem_solver_method_examples.sh [mode] [extra run_lbm_dem.py args...]

Modes:
  bgk_hertz   Fluid: LBM BGK + Guo forcing, DEM: Hertz contact (default)
  trt_hertz   Fluid: LBM TRT + Guo forcing, DEM: Hertz contact
  bgk_linear  Fluid: LBM BGK + Guo forcing, DEM: linear normal contact
  trt_linear  Fluid: LBM TRT + Guo forcing, DEM: linear normal contact
  trt_ibm     Fluid: LBM TRT + Guo forcing, DEM: Hertz, immersed boundary
  all         Run all five method combinations sequentially

Examples:
  src/bin/run_lbm_dem_solver_method_examples.sh bgk_hertz
  src/bin/run_lbm_dem_solver_method_examples.sh trt_ibm --ibm-stiffness 0.8
  src/bin/run_lbm_dem_solver_method_examples.sh trt_linear --total-steps 3000 --no-video
  src/bin/run_lbm_dem_solver_method_examples.sh all --particle-volume-fraction 0.10

Outputs:
  Results are saved under src/results/run_lbm_dem/... with timestamped folders.
  The selected methods are recorded in each run's metadata.json.
USAGE
}

mode="${1:-bgk_hertz}"
if [[ "$mode" == "-h" || "$mode" == "--help" ]]; then
  usage
  exit 0
fi
shift || true

run_case() {
  local tag="$1"
  local fluid_method="$2"
  local particle_method="$3"
  local coupling_method="$4"
  shift 4

  local log_dir="src/results/run_lbm_dem/solver_method_examples"
  local log_file="${log_dir}/${tag}.log"
  mkdir -p "$log_dir"

  local cmd=(
    python -u src/demos/run_lbm_dem.py
    --fluid-method "$fluid_method"
    --particle-method "$particle_method"
    --particle-fluid-coupling "$coupling_method"
    --cylinder
    --particle-source left-inlet
    --particle-volume-fraction 0.05
    --rolling-friction
    --total-steps 2000
    --snapshot-every 100
    --warmup-steps 200
    --paraview-every 5
    --result-tag "$tag"
    "$@"
  )

  echo "=== START ${tag} ===" | tee "$log_file"
  printf "Command:" | tee -a "$log_file"
  printf " %q" "${cmd[@]}" | tee -a "$log_file"
  printf "\n" | tee -a "$log_file"

  "${cmd[@]}" >> "$log_file" 2>&1

  echo "=== OK ${tag} ===" | tee -a "$log_file"
  echo "Log: ${log_file}"
}

case "$mode" in
  bgk_hertz)
    run_case "method_bgk_hertz" "lbm-bgk-guo" "dem-hertz" "point_force" "$@"
    ;;
  trt_hertz)
    run_case "method_trt_hertz" "lbm-trt-guo" "dem-hertz" "point_force" "$@"
    ;;
  bgk_linear)
    run_case "method_bgk_linear" "lbm-bgk-guo" "dem-linear" "point_force" "$@"
    ;;
  trt_linear)
    run_case "method_trt_linear" "lbm-trt-guo" "dem-linear" "point_force" "$@"
    ;;
  trt_ibm)
    run_case "method_trt_hertz_immersed_boundary" "lbm-trt-guo" "dem-hertz" "immersed_boundary" "$@"
    ;;
  all)
    run_case "method_bgk_hertz" "lbm-bgk-guo" "dem-hertz" "point_force" "$@"
    run_case "method_trt_hertz" "lbm-trt-guo" "dem-hertz" "point_force" "$@"
    run_case "method_bgk_linear" "lbm-bgk-guo" "dem-linear" "point_force" "$@"
    run_case "method_trt_linear" "lbm-trt-guo" "dem-linear" "point_force" "$@"
    run_case "method_trt_hertz_immersed_boundary" "lbm-trt-guo" "dem-hertz" "immersed_boundary" "$@"
    ;;
  *)
    echo "Unknown mode: ${mode}" >&2
    usage >&2
    exit 2
    ;;
esac
