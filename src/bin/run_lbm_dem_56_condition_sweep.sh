#!/usr/bin/env bash
set -u

cd "$(dirname "$0")/../.."

export PYTHONPATH="${PWD}/src${PYTHONPATH:+:${PYTHONPATH}}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/matplotlib}"
mkdir -p "$MPLCONFIGDIR"

TOTAL_STEPS=60000
SNAPSHOT_EVERY=120
SWEEP_PREFIX="t3x"

LOG_DIR="src/results/run_lbm_dem/sweep_56_conditions_3x"
mkdir -p "$LOG_DIR"
STATUS_FILE="$LOG_DIR/status.tsv"
SUMMARY_FILE="$LOG_DIR/skipped_conditions.txt"

printf "condition\tphi\tsurface_mode\tstrength_multiplier\trolling\tstatus\texit_code\tlog_file\n" > "$STATUS_FILE"
: > "$SUMMARY_FILE"

run_condition() {
  local phi="$1"
  local surface_mode="$2"
  local strength_multiplier="$3"
  local rolling="$4"
  local tag="$5"
  local strength="$7"

  local log_file="$LOG_DIR/${tag}.log"
  local cmd=(python -u src/demos/run_lbm_dem.py
    --cylinder
    --particle-source left-inlet
    --particle-volume-fraction "$phi"
    --total-steps "$TOTAL_STEPS"
    --snapshot-every "$SNAPSHOT_EVERY"
    --result-tag "$tag"
  )

  case "$surface_mode" in
    attraction)
      cmd+=(--particle-attraction --attraction-strength "$strength")
      ;;
    repulsion)
      cmd+=(--particle-repulsion --repulsion-strength "$strength")
      ;;
    none)
      ;;
    *)
      echo "Unknown surface mode: $surface_mode" | tee "$log_file"
      printf "%s\t%s\t%s\t%s\t%s\tfailed\t2\t%s\n" "$tag" "$phi" "$surface_mode" "$strength_multiplier" "$rolling" "$log_file" >> "$STATUS_FILE"
      printf "%s\t%s\t%s\t%s\t%s\n" "$tag" "$phi" "$surface_mode" "$strength_multiplier" "$rolling" >> "$SUMMARY_FILE"
      return 2
      ;;
  esac

  if [[ "$rolling" == "rolling" ]]; then
    cmd+=(--rolling-friction)
  else
    cmd+=(--no-rolling-friction)
  fi

  echo "=== START $tag ===" | tee "$log_file"
  printf "Command:" | tee -a "$log_file"
  printf " %q" "${cmd[@]}" | tee -a "$log_file"
  printf "\n" | tee -a "$log_file"

  "${cmd[@]}" >> "$log_file" 2>&1
  local exit_code=$?
  if [[ "$exit_code" -eq 0 ]]; then
    printf "%s\t%s\t%s\t%s\t%s\tok\t0\t%s\n" "$tag" "$phi" "$surface_mode" "$strength_multiplier" "$rolling" "$log_file" >> "$STATUS_FILE"
    echo "=== OK $tag ===" | tee -a "$log_file"
  else
    printf "%s\t%s\t%s\t%s\t%s\tfailed\t%s\t%s\n" "$tag" "$phi" "$surface_mode" "$strength_multiplier" "$rolling" "$exit_code" "$log_file" >> "$STATUS_FILE"
    printf "%s\t%s\t%s\t%s\t%s\n" "$tag" "$phi" "$surface_mode" "$strength_multiplier" "$rolling" >> "$SUMMARY_FILE"
    echo "=== FAILED $tag exit=$exit_code ===" | tee -a "$log_file"
  fi
  return 0
}

for phi in 0.05 0.10 0.20 0.40; do
  case "$phi" in
    0.05) phi_tag="05pct" ;;
    0.10) phi_tag="10pct" ;;
    0.20) phi_tag="20pct" ;;
    0.40) phi_tag="40pct" ;;
  esac
  for rolling in rolling free_roll; do
    run_condition "$phi" none 0 "$rolling" "${SWEEP_PREFIX}_surface_none_${rolling}_${phi_tag}" none 0.0
    run_condition "$phi" attraction 1.0 "$rolling" "${SWEEP_PREFIX}_attraction_1x_${rolling}_${phi_tag}" attraction 0.001
    run_condition "$phi" attraction 2.0 "$rolling" "${SWEEP_PREFIX}_attraction_2x_${rolling}_${phi_tag}" attraction 0.002
    run_condition "$phi" attraction 0.5 "$rolling" "${SWEEP_PREFIX}_attraction_0p5x_${rolling}_${phi_tag}" attraction 0.0005
    run_condition "$phi" repulsion 1.0 "$rolling" "${SWEEP_PREFIX}_repulsion_1x_${rolling}_${phi_tag}" repulsion 0.001
    run_condition "$phi" repulsion 2.0 "$rolling" "${SWEEP_PREFIX}_repulsion_2x_${rolling}_${phi_tag}" repulsion 0.002
    run_condition "$phi" repulsion 0.5 "$rolling" "${SWEEP_PREFIX}_repulsion_0p5x_${rolling}_${phi_tag}" repulsion 0.0005
  done
done

if [[ ! -s "$SUMMARY_FILE" ]]; then
  echo "No skipped or failed conditions." > "$SUMMARY_FILE"
fi

echo "Sweep complete."
echo "Status:  $STATUS_FILE"
echo "Skipped: $SUMMARY_FILE"
