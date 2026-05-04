#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/../.."

export PYTHONPATH="${PWD}/src${PYTHONPATH:+:${PYTHONPATH}}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/matplotlib-cfd-solver}"
mkdir -p "$MPLCONFIGDIR"

NX=540
NY=210
PARTICLE_RADIUS=9
PHI=0.40
ATTRACTION_STRENGTH=0.002
ROLLING_FRICTION_COEFF=0.15
TOTAL_STEPS=600000
SNAPSHOT_EVERY=1200
PARAVIEW_EVERY=10
MAX_WORKERS="${MAX_WORKERS:-4}"
LOG_ROOT="src/results/run_lbm_dem/highres_pore_geometry_sweep"
STATUS_FILE="${LOG_ROOT}/status.tsv"
SUMMARY_FILE="${LOG_ROOT}/skipped_conditions.txt"
PREVIOUS_TAG="four_cylinders_attr2x_phi40_muR3x"
PREVIOUS_LOG="src/results/run_lbm_dem/four_cylinder_cases/${PREVIOUS_TAG}.log"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --max-workers)
      MAX_WORKERS="$2"
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

mkdir -p "$LOG_ROOT"
printf "condition\tnx\tny\tparticle_radius\tphi\tattraction_strength\trolling_friction_coeff\tstatus\texit_code\tlog_file\n" > "$STATUS_FILE"
: > "$SUMMARY_FILE"

wait_for_previous_run() {
  while [[ ! -f "$PREVIOUS_LOG" ]] || ! grep -q "=== OK ${PREVIOUS_TAG} ===" "$PREVIOUS_LOG"; do
    echo "Waiting for previous run (${PREVIOUS_TAG}) to finish..."
    if [[ -f "$PREVIOUS_LOG" ]]; then
      tail -n 5 "$PREVIOUS_LOG" || true
    else
      echo "Previous log does not exist yet: $PREVIOUS_LOG"
    fi
    sleep 60
  done
}

run_case() {
  local tag="$1"
  shift

  local log_file="${LOG_ROOT}/${tag}.log"
  local cmd=(
    python -u src/demos/run_lbm_dem.py
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
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\tok\t0\t%s\n" \
      "$tag" "$NX" "$NY" "$PARTICLE_RADIUS" "$PHI" "$ATTRACTION_STRENGTH" \
      "$ROLLING_FRICTION_COEFF" "$log_file" >> "$STATUS_FILE"
    echo "=== OK ${tag} ===" | tee -a "$log_file"
  else
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\tfailed\t%s\t%s\n" \
      "$tag" "$NX" "$NY" "$PARTICLE_RADIUS" "$PHI" "$ATTRACTION_STRENGTH" \
      "$ROLLING_FRICTION_COEFF" "$exit_code" "$log_file" >> "$STATUS_FILE"
    printf "%s\n" "$tag" >> "$SUMMARY_FILE"
    echo "=== FAILED ${tag} exit=${exit_code} ===" | tee -a "$log_file"
  fi
}

launch_case() {
  run_case "$@" &
  while [[ "$(jobs -pr | wc -l)" -ge "$MAX_WORKERS" ]]; do
    wait -n || true
  done
}

wait_for_previous_run
echo "Launching high-resolution pore geometry sweep with max_workers=${MAX_WORKERS}"

launch_case highres_2cyl_r18_attr2x_phi40_muR3x \
  --cylinder-spec 270 75 18 \
  --cylinder-spec 270 135 18

launch_case highres_3cyl_r18_attr2x_phi40_muR3x \
  --cylinder-spec 210 105 18 \
  --cylinder-spec 270 105 18 \
  --cylinder-spec 330 105 18

launch_case highres_4cyl_r18_attr2x_phi40_muR3x \
  --cylinder-spec 225 75 18 \
  --cylinder-spec 225 135 18 \
  --cylinder-spec 315 75 18 \
  --cylinder-spec 315 135 18

launch_case highres_4cyl_r21_attr2x_phi40_muR3x \
  --cylinder-spec 225 75 21 \
  --cylinder-spec 225 135 21 \
  --cylinder-spec 315 75 21 \
  --cylinder-spec 315 135 21

launch_case highres_4cyl_r15_attr2x_phi40_muR3x \
  --cylinder-spec 225 75 15 \
  --cylinder-spec 225 135 15 \
  --cylinder-spec 315 75 15 \
  --cylinder-spec 315 135 15

launch_case highres_6cyl_r18_attr2x_phi40_muR3x \
  --cylinder-spec 180 75 18 \
  --cylinder-spec 180 135 18 \
  --cylinder-spec 270 75 18 \
  --cylinder-spec 270 135 18 \
  --cylinder-spec 360 75 18 \
  --cylinder-spec 360 135 18

launch_case highres_6cyl_r15_attr2x_phi40_muR3x \
  --cylinder-spec 180 75 15 \
  --cylinder-spec 180 135 15 \
  --cylinder-spec 270 75 15 \
  --cylinder-spec 270 135 15 \
  --cylinder-spec 360 75 15 \
  --cylinder-spec 360 135 15

launch_case highres_8cyl_r12_attr2x_phi40_muR3x \
  --cylinder-spec 150 75 12 \
  --cylinder-spec 150 135 12 \
  --cylinder-spec 230 75 12 \
  --cylinder-spec 230 135 12 \
  --cylinder-spec 310 75 12 \
  --cylinder-spec 310 135 12 \
  --cylinder-spec 390 75 12 \
  --cylinder-spec 390 135 12

launch_case highres_9cyl_r12_attr2x_phi40_muR3x \
  --cylinder-spec 180 60 12 \
  --cylinder-spec 180 105 12 \
  --cylinder-spec 180 150 12 \
  --cylinder-spec 270 60 12 \
  --cylinder-spec 270 105 12 \
  --cylinder-spec 270 150 12 \
  --cylinder-spec 360 60 12 \
  --cylinder-spec 360 105 12 \
  --cylinder-spec 360 150 12

wait

if [[ ! -s "$SUMMARY_FILE" ]]; then
  echo "No skipped or failed conditions." > "$SUMMARY_FILE"
fi

echo "High-resolution pore geometry sweep complete."
echo "Status:  $STATUS_FILE"
echo "Skipped: $SUMMARY_FILE"
