#!/usr/bin/env bash
set -u

cd "$(dirname "$0")/../.."

export PYTHONPATH="${PWD}/src${PYTHONPATH:+:${PYTHONPATH}}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/matplotlib}"
mkdir -p "$MPLCONFIGDIR"

TOTAL_STEPS=60000
SNAPSHOT_EVERY=120
BASE_ROLLING_FRICTION_COEFF=0.05

LOG_ROOT="src/results/run_lbm_dem/sweep_56_conditions_rolling_friction_2x_3x"
mkdir -p "$LOG_ROOT"
STATUS_FILE="$LOG_ROOT/status.tsv"
SUMMARY_FILE="$LOG_ROOT/skipped_conditions.txt"

printf "condition\tphi\tsurface_mode\tstrength_multiplier\trolling\trolling_friction_multiplier\tstatus\texit_code\tlog_file\n" > "$STATUS_FILE"
: > "$SUMMARY_FILE"

rolling_coeff_for_multiplier() {
  case "$1" in
    2) echo "0.10" ;;
    3) echo "0.15" ;;
    *)
      awk -v base="$BASE_ROLLING_FRICTION_COEFF" -v mult="$1" 'BEGIN { printf "%.8g\n", base * mult }'
      ;;
  esac
}

run_condition() {
  local phi="$1"
  local surface_mode="$2"
  local strength_multiplier="$3"
  local rolling="$4"
  local rolling_multiplier="$5"
  local tag="$6"
  local strength="$8"

  local coeff
  coeff="$(rolling_coeff_for_multiplier "$rolling_multiplier")"

  local log_dir="$LOG_ROOT/rolling_friction_${rolling_multiplier}x"
  mkdir -p "$log_dir"
  local log_file="$log_dir/${tag}.log"
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
      printf "%s\t%s\t%s\t%s\t%s\t%s\tfailed\t2\t%s\n" "$tag" "$phi" "$surface_mode" "$strength_multiplier" "$rolling" "$rolling_multiplier" "$log_file" >> "$STATUS_FILE"
      printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$tag" "$phi" "$surface_mode" "$strength_multiplier" "$rolling" "$rolling_multiplier" >> "$SUMMARY_FILE"
      return 2
      ;;
  esac

  if [[ "$rolling" == "rolling" ]]; then
    cmd+=(--rolling-friction --rolling-friction-coeff "$coeff")
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
    printf "%s\t%s\t%s\t%s\t%s\t%s\tok\t0\t%s\n" "$tag" "$phi" "$surface_mode" "$strength_multiplier" "$rolling" "$rolling_multiplier" "$log_file" >> "$STATUS_FILE"
    echo "=== OK $tag ===" | tee -a "$log_file"
  else
    printf "%s\t%s\t%s\t%s\t%s\t%s\tfailed\t%s\t%s\n" "$tag" "$phi" "$surface_mode" "$strength_multiplier" "$rolling" "$rolling_multiplier" "$exit_code" "$log_file" >> "$STATUS_FILE"
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$tag" "$phi" "$surface_mode" "$strength_multiplier" "$rolling" "$rolling_multiplier" >> "$SUMMARY_FILE"
    echo "=== FAILED $tag exit=$exit_code ===" | tee -a "$log_file"
  fi
  return 0
}

for rolling_multiplier in 2 3; do
  sweep_prefix="t3x_muR${rolling_multiplier}x"
  for phi in 0.05 0.10 0.20 0.40; do
    case "$phi" in
      0.05) phi_tag="05pct" ;;
      0.10) phi_tag="10pct" ;;
      0.20) phi_tag="20pct" ;;
      0.40) phi_tag="40pct" ;;
    esac
    for rolling in rolling free_roll; do
      run_condition "$phi" none 0 "$rolling" "$rolling_multiplier" "${sweep_prefix}_surface_none_${rolling}_${phi_tag}" none 0.0
      run_condition "$phi" attraction 1.0 "$rolling" "$rolling_multiplier" "${sweep_prefix}_attraction_1x_${rolling}_${phi_tag}" attraction 0.001
      run_condition "$phi" attraction 2.0 "$rolling" "$rolling_multiplier" "${sweep_prefix}_attraction_2x_${rolling}_${phi_tag}" attraction 0.002
      run_condition "$phi" attraction 0.5 "$rolling" "$rolling_multiplier" "${sweep_prefix}_attraction_0p5x_${rolling}_${phi_tag}" attraction 0.0005
      run_condition "$phi" repulsion 1.0 "$rolling" "$rolling_multiplier" "${sweep_prefix}_repulsion_1x_${rolling}_${phi_tag}" repulsion 0.001
      run_condition "$phi" repulsion 2.0 "$rolling" "$rolling_multiplier" "${sweep_prefix}_repulsion_2x_${rolling}_${phi_tag}" repulsion 0.002
      run_condition "$phi" repulsion 0.5 "$rolling" "$rolling_multiplier" "${sweep_prefix}_repulsion_0p5x_${rolling}_${phi_tag}" repulsion 0.0005
    done
  done
done

if [[ ! -s "$SUMMARY_FILE" ]]; then
  echo "No skipped or failed conditions." > "$SUMMARY_FILE"
fi

echo "Sweep complete."
echo "Status:  $STATUS_FILE"
echo "Skipped: $SUMMARY_FILE"
