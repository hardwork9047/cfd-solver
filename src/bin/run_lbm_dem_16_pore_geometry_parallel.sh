#!/usr/bin/env bash
set -euo pipefail
trap '' HUP

cd "$(dirname "$0")/../.."

export PYTHONPATH="${PWD}/src${PYTHONPATH:+:${PYTHONPATH}}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/matplotlib-cfd-solver}"
mkdir -p "$MPLCONFIGDIR"

NX="${NX:-540}"
NY="${NY:-210}"
PARTICLE_RADIUS="${PARTICLE_RADIUS:-9}"
PHI="${PHI:-0.40}"
ATTRACTION_STRENGTH="${ATTRACTION_STRENGTH:-0.002}"
ROLLING_FRICTION_COEFF="${ROLLING_FRICTION_COEFF:-0.15}"
TOTAL_STEPS="${TOTAL_STEPS:-600000}"
SNAPSHOT_EVERY="${SNAPSHOT_EVERY:-1200}"
PARAVIEW_EVERY="${PARAVIEW_EVERY:-10}"
MAX_WORKERS="${MAX_WORKERS:-16}"
LOG_ROOT="src/results/run_lbm_dem/pore_geometry_16_parallel"
STATUS_FILE="${LOG_ROOT}/status.tsv"
SUMMARY_FILE="${LOG_ROOT}/skipped_conditions.txt"
LOCK_DIR="${LOG_ROOT}/locks"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --max-workers)
      MAX_WORKERS="$2"
      shift 2
      ;;
    --total-steps)
      TOTAL_STEPS="$2"
      shift 2
      ;;
    --snapshot-every)
      SNAPSHOT_EVERY="$2"
      shift 2
      ;;
    --paraview-every)
      PARAVIEW_EVERY="$2"
      shift 2
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 2
      ;;
  esac
done

if [[ "$MAX_WORKERS" -lt 1 ]]; then
  echo "--max-workers must be >= 1" >&2
  exit 2
fi

mkdir -p "$LOG_ROOT" "$LOCK_DIR"
printf "condition\tnx\tny\tparticle_radius\tphi\tattraction_strength\trolling_friction_coeff\tstatus\texit_code\tlog_file\n" > "$STATUS_FILE"
: > "$SUMMARY_FILE"

append_status() {
  local line="$1"
  (
    flock 9
    printf "%s\n" "$line" >> "$STATUS_FILE"
  ) 9>"${LOCK_DIR}/status.lock"
}

append_skipped() {
  local tag="$1"
  (
    flock 9
    printf "%s\n" "$tag" >> "$SUMMARY_FILE"
  ) 9>"${LOCK_DIR}/skipped.lock"
}

run_case() {
  local tag="$1"
  shift

  local log_file="${LOG_ROOT}/${tag}.log"
  local cmd=(
    python -u src/runners/run_lbm_dem.py
    --nx "$NX"
    --ny "$NY"
    --particle-radius "$PARTICLE_RADIUS"
    --particle-source left-inlet
    --particle-volume-fraction "$PHI"
    --particle-attraction
    --attraction-strength "$ATTRACTION_STRENGTH"
    --rolling-friction
    --rolling-friction-coeff "$ROLLING_FRICTION_COEFF"
    --total-steps "$TOTAL_STEPS"
    --snapshot-every "$SNAPSHOT_EVERY"
    --paraview-every "$PARAVIEW_EVERY"
    --result-tag "$tag"
    "$@"
  )

  echo "=== START ${tag} ===" | tee "$log_file"
  printf "Command:" | tee -a "$log_file"
  printf " %q" "${cmd[@]}" | tee -a "$log_file"
  printf "\n" | tee -a "$log_file"

  set +e
  "${cmd[@]}" >> "$log_file" 2>&1
  local exit_code=$?
  set -e

  if [[ "$exit_code" -eq 0 ]]; then
    append_status "${tag}	${NX}	${NY}	${PARTICLE_RADIUS}	${PHI}	${ATTRACTION_STRENGTH}	${ROLLING_FRICTION_COEFF}	ok	0	${log_file}"
    echo "=== OK ${tag} ===" | tee -a "$log_file"
  else
    append_status "${tag}	${NX}	${NY}	${PARTICLE_RADIUS}	${PHI}	${ATTRACTION_STRENGTH}	${ROLLING_FRICTION_COEFF}	failed	${exit_code}	${log_file}"
    append_skipped "$tag"
    echo "=== FAILED ${tag} exit=${exit_code} ===" | tee -a "$log_file"
  fi
}

launch_case() {
  run_case "$@" &
  while [[ "$(jobs -pr | wc -l)" -ge "$MAX_WORKERS" ]]; do
    wait -n || true
  done
}

echo "Launching 16 pore-geometry cases with max_workers=${MAX_WORKERS}"
echo "Grid=${NX}x${NY}  phi=${PHI}  total_steps=${TOTAL_STEPS}  paraview_every=${PARAVIEW_EVERY}"

launch_case pore16_01_2cyl_vertical_r18 \
  --cylinder-spec 270 75 18 \
  --cylinder-spec 270 135 18

launch_case pore16_02_3cyl_centerline_r18 \
  --cylinder-spec 210 105 18 \
  --cylinder-spec 270 105 18 \
  --cylinder-spec 330 105 18

launch_case pore16_03_4cyl_square_r18 \
  --cylinder-spec 225 75 18 \
  --cylinder-spec 225 135 18 \
  --cylinder-spec 315 75 18 \
  --cylinder-spec 315 135 18

launch_case pore16_04_4cyl_square_r21 \
  --cylinder-spec 225 75 21 \
  --cylinder-spec 225 135 21 \
  --cylinder-spec 315 75 21 \
  --cylinder-spec 315 135 21

launch_case pore16_05_4cyl_square_r15 \
  --cylinder-spec 225 75 15 \
  --cylinder-spec 225 135 15 \
  --cylinder-spec 315 75 15 \
  --cylinder-spec 315 135 15

launch_case pore16_06_5cyl_cross_r16 \
  --cylinder-spec 210 75 16 \
  --cylinder-spec 210 135 16 \
  --cylinder-spec 270 105 16 \
  --cylinder-spec 330 75 16 \
  --cylinder-spec 330 135 16

launch_case pore16_07_6cyl_grid_r18 \
  --cylinder-spec 180 75 18 \
  --cylinder-spec 180 135 18 \
  --cylinder-spec 270 75 18 \
  --cylinder-spec 270 135 18 \
  --cylinder-spec 360 75 18 \
  --cylinder-spec 360 135 18

launch_case pore16_08_6cyl_grid_r15 \
  --cylinder-spec 180 75 15 \
  --cylinder-spec 180 135 15 \
  --cylinder-spec 270 75 15 \
  --cylinder-spec 270 135 15 \
  --cylinder-spec 360 75 15 \
  --cylinder-spec 360 135 15

launch_case pore16_09_6cyl_staggered_r12 \
  --cylinder-spec 170 70 12 \
  --cylinder-spec 170 140 12 \
  --cylinder-spec 270 105 12 \
  --cylinder-spec 370 70 12 \
  --cylinder-spec 370 140 12 \
  --cylinder-spec 430 105 12

launch_case pore16_10_8cyl_grid_r12 \
  --cylinder-spec 150 75 12 \
  --cylinder-spec 150 135 12 \
  --cylinder-spec 230 75 12 \
  --cylinder-spec 230 135 12 \
  --cylinder-spec 310 75 12 \
  --cylinder-spec 310 135 12 \
  --cylinder-spec 390 75 12 \
  --cylinder-spec 390 135 12

launch_case pore16_11_8cyl_grid_r15 \
  --cylinder-spec 150 75 15 \
  --cylinder-spec 150 135 15 \
  --cylinder-spec 230 75 15 \
  --cylinder-spec 230 135 15 \
  --cylinder-spec 310 75 15 \
  --cylinder-spec 310 135 15 \
  --cylinder-spec 390 75 15 \
  --cylinder-spec 390 135 15

launch_case pore16_12_9cyl_grid_r12 \
  --cylinder-spec 180 60 12 \
  --cylinder-spec 180 105 12 \
  --cylinder-spec 180 150 12 \
  --cylinder-spec 270 60 12 \
  --cylinder-spec 270 105 12 \
  --cylinder-spec 270 150 12 \
  --cylinder-spec 360 60 12 \
  --cylinder-spec 360 105 12 \
  --cylinder-spec 360 150 12

launch_case pore16_13_10cyl_staggered_r11 \
  --cylinder-spec 150 70 11 \
  --cylinder-spec 150 140 11 \
  --cylinder-spec 220 105 11 \
  --cylinder-spec 270 70 11 \
  --cylinder-spec 270 140 11 \
  --cylinder-spec 320 105 11 \
  --cylinder-spec 390 70 11 \
  --cylinder-spec 390 140 11 \
  --cylinder-spec 440 90 11 \
  --cylinder-spec 440 120 11

launch_case pore16_14_12cyl_grid_r10 \
  --cylinder-spec 150 60 10 \
  --cylinder-spec 150 105 10 \
  --cylinder-spec 150 150 10 \
  --cylinder-spec 230 60 10 \
  --cylinder-spec 230 105 10 \
  --cylinder-spec 230 150 10 \
  --cylinder-spec 310 60 10 \
  --cylinder-spec 310 105 10 \
  --cylinder-spec 310 150 10 \
  --cylinder-spec 390 60 10 \
  --cylinder-spec 390 105 10 \
  --cylinder-spec 390 150 10

launch_case pore16_15_12cyl_mixed_r10_r14 \
  --cylinder-spec 150 60 10 \
  --cylinder-spec 150 105 14 \
  --cylinder-spec 150 150 10 \
  --cylinder-spec 230 60 14 \
  --cylinder-spec 230 105 10 \
  --cylinder-spec 230 150 14 \
  --cylinder-spec 310 60 10 \
  --cylinder-spec 310 105 14 \
  --cylinder-spec 310 150 10 \
  --cylinder-spec 390 60 14 \
  --cylinder-spec 390 105 10 \
  --cylinder-spec 390 150 14

launch_case pore16_16_16cyl_dense_r8 \
  --cylinder-spec 150 50 8 \
  --cylinder-spec 150 85 8 \
  --cylinder-spec 150 125 8 \
  --cylinder-spec 150 160 8 \
  --cylinder-spec 240 50 8 \
  --cylinder-spec 240 85 8 \
  --cylinder-spec 240 125 8 \
  --cylinder-spec 240 160 8 \
  --cylinder-spec 330 50 8 \
  --cylinder-spec 330 85 8 \
  --cylinder-spec 330 125 8 \
  --cylinder-spec 330 160 8 \
  --cylinder-spec 420 50 8 \
  --cylinder-spec 420 85 8 \
  --cylinder-spec 420 125 8 \
  --cylinder-spec 420 160 8

wait || true

if [[ ! -s "$SUMMARY_FILE" ]]; then
  echo "No skipped or failed conditions." > "$SUMMARY_FILE"
fi

echo "16-condition pore geometry sweep complete."
echo "Status:  $STATUS_FILE"
echo "Skipped: $SUMMARY_FILE"
