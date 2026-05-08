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
TOTAL_STEPS="${TOTAL_STEPS:-600000}"
SNAPSHOT_EVERY="${SNAPSHOT_EVERY:-1200}"
PARAVIEW_EVERY="${PARAVIEW_EVERY:-0}"
MAX_WORKERS="${MAX_WORKERS:-8}"
NO_VIDEO="${NO_VIDEO:-1}"
POST_ANALYZE="${POST_ANALYZE:-1}"

LOG_ROOT="src/results/run_lbm_dem/design_sweep_100_parallel"
STATUS_FILE="${LOG_ROOT}/status.tsv"
MANIFEST_FILE="${LOG_ROOT}/condition_manifest.tsv"
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
    --with-video)
      NO_VIDEO=0
      shift
      ;;
    --no-post-analyze)
      POST_ANALYZE=0
      shift
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
printf "condition\tgeometry_id\tgeometry_label\tcylinder_count\tmean_cylinder_radius\tnx\tny\tparticle_radius\treynolds_number\tphi\tsurface_force\tattraction_strength\trolling_friction\trolling_friction_coeff\tstatus\texit_code\tlog_file\n" > "$STATUS_FILE"
printf "condition\tgeometry_id\tgeometry_label\tcylinder_count\tmean_cylinder_radius\treynolds_number\tphi\tsurface_force\tattraction_strength\trolling_friction\trolling_friction_coeff\n" > "$MANIFEST_FILE"
: > "$SUMMARY_FILE"

geometry_args() {
  local geometry_id="$1"
  GEOMETRY_LABEL=""
  CYLINDER_COUNT=""
  MEAN_CYLINDER_RADIUS=""
  CYL_ARGS=()

  case "$geometry_id" in
    g01)
      GEOMETRY_LABEL="2cyl_vertical_r18"
      CYLINDER_COUNT=2
      MEAN_CYLINDER_RADIUS=18
      CYL_ARGS=(--cylinder-spec 270 75 18 --cylinder-spec 270 135 18)
      ;;
    g02)
      GEOMETRY_LABEL="4cyl_square_r15"
      CYLINDER_COUNT=4
      MEAN_CYLINDER_RADIUS=15
      CYL_ARGS=(--cylinder-spec 225 75 15 --cylinder-spec 225 135 15 --cylinder-spec 315 75 15 --cylinder-spec 315 135 15)
      ;;
    g03)
      GEOMETRY_LABEL="4cyl_square_r21"
      CYLINDER_COUNT=4
      MEAN_CYLINDER_RADIUS=21
      CYL_ARGS=(--cylinder-spec 225 75 21 --cylinder-spec 225 135 21 --cylinder-spec 315 75 21 --cylinder-spec 315 135 21)
      ;;
    g04)
      GEOMETRY_LABEL="5cyl_cross_r16"
      CYLINDER_COUNT=5
      MEAN_CYLINDER_RADIUS=16
      CYL_ARGS=(--cylinder-spec 210 75 16 --cylinder-spec 210 135 16 --cylinder-spec 270 105 16 --cylinder-spec 330 75 16 --cylinder-spec 330 135 16)
      ;;
    g05)
      GEOMETRY_LABEL="6cyl_grid_r15"
      CYLINDER_COUNT=6
      MEAN_CYLINDER_RADIUS=15
      CYL_ARGS=(--cylinder-spec 180 75 15 --cylinder-spec 180 135 15 --cylinder-spec 270 75 15 --cylinder-spec 270 135 15 --cylinder-spec 360 75 15 --cylinder-spec 360 135 15)
      ;;
    g06)
      GEOMETRY_LABEL="6cyl_staggered_r12"
      CYLINDER_COUNT=6
      MEAN_CYLINDER_RADIUS=12
      CYL_ARGS=(--cylinder-spec 170 70 12 --cylinder-spec 170 140 12 --cylinder-spec 270 105 12 --cylinder-spec 370 70 12 --cylinder-spec 370 140 12 --cylinder-spec 430 105 12)
      ;;
    g07)
      GEOMETRY_LABEL="8cyl_grid_r12"
      CYLINDER_COUNT=8
      MEAN_CYLINDER_RADIUS=12
      CYL_ARGS=(--cylinder-spec 150 75 12 --cylinder-spec 150 135 12 --cylinder-spec 230 75 12 --cylinder-spec 230 135 12 --cylinder-spec 310 75 12 --cylinder-spec 310 135 12 --cylinder-spec 390 75 12 --cylinder-spec 390 135 12)
      ;;
    g08)
      GEOMETRY_LABEL="9cyl_grid_r12"
      CYLINDER_COUNT=9
      MEAN_CYLINDER_RADIUS=12
      CYL_ARGS=(--cylinder-spec 180 60 12 --cylinder-spec 180 105 12 --cylinder-spec 180 150 12 --cylinder-spec 270 60 12 --cylinder-spec 270 105 12 --cylinder-spec 270 150 12 --cylinder-spec 360 60 12 --cylinder-spec 360 105 12 --cylinder-spec 360 150 12)
      ;;
    g09)
      GEOMETRY_LABEL="12cyl_mixed_r10_r14"
      CYLINDER_COUNT=12
      MEAN_CYLINDER_RADIUS=12
      CYL_ARGS=(--cylinder-spec 150 60 10 --cylinder-spec 150 105 14 --cylinder-spec 150 150 10 --cylinder-spec 230 60 14 --cylinder-spec 230 105 10 --cylinder-spec 230 150 14 --cylinder-spec 310 60 10 --cylinder-spec 310 105 14 --cylinder-spec 310 150 10 --cylinder-spec 390 60 14 --cylinder-spec 390 105 10 --cylinder-spec 390 150 14)
      ;;
    g10)
      GEOMETRY_LABEL="16cyl_dense_r8"
      CYLINDER_COUNT=16
      MEAN_CYLINDER_RADIUS=8
      CYL_ARGS=(--cylinder-spec 150 50 8 --cylinder-spec 150 85 8 --cylinder-spec 150 125 8 --cylinder-spec 150 160 8 --cylinder-spec 240 50 8 --cylinder-spec 240 85 8 --cylinder-spec 240 125 8 --cylinder-spec 240 160 8 --cylinder-spec 330 50 8 --cylinder-spec 330 85 8 --cylinder-spec 330 125 8 --cylinder-spec 330 160 8 --cylinder-spec 420 50 8 --cylinder-spec 420 85 8 --cylinder-spec 420 125 8 --cylinder-spec 420 160 8)
      ;;
    *)
      echo "Unknown geometry_id: $geometry_id" >&2
      exit 2
      ;;
  esac
}

append_status() {
  local line="$1"
  (
    flock 9
    printf "%s\n" "$line" >> "$STATUS_FILE"
  ) 9>"${LOCK_DIR}/status.lock"
}

append_manifest() {
  local line="$1"
  (
    flock 9
    printf "%s\n" "$line" >> "$MANIFEST_FILE"
  ) 9>"${LOCK_DIR}/manifest.lock"
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
  local geometry_id="$2"
  local re="$3"
  local phi="$4"
  local surface_force="$5"
  local attraction_strength="$6"
  local rolling="$7"
  local rolling_coeff="$8"

  geometry_args "$geometry_id"
  local log_file="${LOG_ROOT}/${tag}.log"
  local cmd=(
    python -u src/demos/run_lbm_dem.py
    --nx "$NX"
    --ny "$NY"
    --particle-radius "$PARTICLE_RADIUS"
    --particle-source left-inlet
    --particle-volume-fraction "$phi"
    --reynolds-number "$re"
    --total-steps "$TOTAL_STEPS"
    --snapshot-every "$SNAPSHOT_EVERY"
    --paraview-every "$PARAVIEW_EVERY"
    --result-tag "$tag"
  )

  if [[ "$surface_force" == "attraction" ]]; then
    cmd+=(--particle-attraction --attraction-strength "$attraction_strength")
  elif [[ "$surface_force" != "none" ]]; then
    echo "Unknown surface_force: $surface_force" >&2
    exit 2
  fi

  if [[ "$rolling" == "rolling" ]]; then
    cmd+=(--rolling-friction --rolling-friction-coeff "$rolling_coeff")
  else
    cmd+=(--no-rolling-friction)
  fi

  if [[ "$NO_VIDEO" -eq 1 ]]; then
    cmd+=(--no-video)
  fi

  cmd+=("${CYL_ARGS[@]}")

  echo "=== START ${tag} ===" | tee "$log_file"
  printf "Command:" | tee -a "$log_file"
  printf " %q" "${cmd[@]}" | tee -a "$log_file"
  printf "\n" | tee -a "$log_file"

  set +e
  "${cmd[@]}" >> "$log_file" 2>&1
  local exit_code=$?
  set -e

  if [[ "$exit_code" -eq 0 ]]; then
    append_status "${tag}	${geometry_id}	${GEOMETRY_LABEL}	${CYLINDER_COUNT}	${MEAN_CYLINDER_RADIUS}	${NX}	${NY}	${PARTICLE_RADIUS}	${re}	${phi}	${surface_force}	${attraction_strength}	${rolling}	${rolling_coeff}	ok	0	${log_file}"
    echo "=== OK ${tag} ===" | tee -a "$log_file"
  else
    append_status "${tag}	${geometry_id}	${GEOMETRY_LABEL}	${CYLINDER_COUNT}	${MEAN_CYLINDER_RADIUS}	${NX}	${NY}	${PARTICLE_RADIUS}	${re}	${phi}	${surface_force}	${attraction_strength}	${rolling}	${rolling_coeff}	failed	${exit_code}	${log_file}"
    append_skipped "$tag"
    echo "=== FAILED ${tag} exit=${exit_code} ===" | tee -a "$log_file"
  fi
}

launch_case() {
  local tag="$1"
  local geometry_id="$2"
  local re="$3"
  local phi="$4"
  local surface_force="$5"
  local attraction_strength="$6"
  local rolling="$7"
  local rolling_coeff="$8"

  geometry_args "$geometry_id"
  append_manifest "${tag}	${geometry_id}	${GEOMETRY_LABEL}	${CYLINDER_COUNT}	${MEAN_CYLINDER_RADIUS}	${re}	${phi}	${surface_force}	${attraction_strength}	${rolling}	${rolling_coeff}"

  run_case "$tag" "$geometry_id" "$re" "$phi" "$surface_force" "$attraction_strength" "$rolling" "$rolling_coeff" &
  while [[ "$(jobs -pr | wc -l)" -ge "$MAX_WORKERS" ]]; do
    wait -n || true
  done
}

GEOMETRIES=(g01 g02 g03 g04 g05 g06 g07 g08 g09 g10)
RE_VALUES=(50 100)
PHYSICS_CASES=(
  "p01_phi05_none_free|0.05|none|0|free|0"
  "p02_phi10_attr1x_free|0.10|attraction|0.001|free|0"
  "p03_phi20_attr05x_roll3x|0.20|attraction|0.0005|rolling|0.15"
  "p04_phi20_attr1x_roll1x|0.20|attraction|0.001|rolling|0.05"
  "p05_phi40_attr2x_roll3x|0.40|attraction|0.002|rolling|0.15"
)

echo "Launching 100 LBM-DEM design cases with max_workers=${MAX_WORKERS}"
echo "Grid=${NX}x${NY}  particle_radius=${PARTICLE_RADIUS}  total_steps=${TOTAL_STEPS}"
echo "Video output: $([[ "$NO_VIDEO" -eq 1 ]] && echo disabled || echo enabled)"

case_index=0
for geometry_id in "${GEOMETRIES[@]}"; do
  for re in "${RE_VALUES[@]}"; do
    for physics_case in "${PHYSICS_CASES[@]}"; do
      IFS="|" read -r physics_label phi surface_force attraction_strength rolling rolling_coeff <<< "$physics_case"
      case_index=$((case_index + 1))
      tag="$(printf "design100_%03d_%s_Re%s_%s" "$case_index" "$geometry_id" "$re" "$physics_label")"
      launch_case "$tag" "$geometry_id" "$re" "$phi" "$surface_force" "$attraction_strength" "$rolling" "$rolling_coeff"
    done
  done
done

wait || true

if [[ ! -s "$SUMMARY_FILE" ]]; then
  echo "No skipped or failed conditions." > "$SUMMARY_FILE"
fi

echo "100-condition design sweep complete."
echo "Status:   $STATUS_FILE"
echo "Manifest: $MANIFEST_FILE"
echo "Skipped:  $SUMMARY_FILE"

if [[ "$POST_ANALYZE" -eq 1 ]]; then
  python src/bin/analyze_lbm_dem_design_sweeps.py
fi
