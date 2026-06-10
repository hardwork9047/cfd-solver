#!/usr/bin/env bash
set -euo pipefail

stage="${1:-stage1}"

case "$stage" in
  stage1)
    configs=(
      "configs/lbm_dem/cases/fouling_3d_sweep_20260609_b0_baseline.json"
      "configs/lbm_dem/cases/fouling_3d_sweep_20260609_a2_attr_0010.json"
      "configs/lbm_dem/cases/fouling_3d_sweep_20260609_a3_attr_0005.json"
    )
    ;;
  stage2)
    configs=(
      "configs/lbm_dem/cases/fouling_3d_sweep_20260609_c1_phi_010.json"
      "configs/lbm_dem/cases/fouling_3d_sweep_20260609_c2_phi_020.json"
      "configs/lbm_dem/cases/fouling_3d_sweep_20260609_f1_roll_005.json"
      "configs/lbm_dem/cases/fouling_3d_sweep_20260609_f2_roll_030.json"
      "configs/lbm_dem/cases/fouling_3d_sweep_20260609_g1_gap_wide.json"
      "configs/lbm_dem/cases/fouling_3d_sweep_20260609_g2_gap_narrow.json"
    )
    ;;
  stage3)
    configs=(
      "configs/lbm_dem/cases/fouling_3d_sweep_20260609_r1_rep_0005.json"
      "configs/lbm_dem/cases/fouling_3d_sweep_20260609_r2_rep_0010.json"
      "configs/lbm_dem/cases/fouling_3d_sweep_20260609_a1_attr_0015.json"
    )
    ;;
  all)
    configs=(
      "configs/lbm_dem/cases/fouling_3d_sweep_20260609_b0_baseline.json"
      "configs/lbm_dem/cases/fouling_3d_sweep_20260609_a2_attr_0010.json"
      "configs/lbm_dem/cases/fouling_3d_sweep_20260609_a3_attr_0005.json"
      "configs/lbm_dem/cases/fouling_3d_sweep_20260609_c1_phi_010.json"
      "configs/lbm_dem/cases/fouling_3d_sweep_20260609_c2_phi_020.json"
      "configs/lbm_dem/cases/fouling_3d_sweep_20260609_f1_roll_005.json"
      "configs/lbm_dem/cases/fouling_3d_sweep_20260609_f2_roll_030.json"
      "configs/lbm_dem/cases/fouling_3d_sweep_20260609_g1_gap_wide.json"
      "configs/lbm_dem/cases/fouling_3d_sweep_20260609_g2_gap_narrow.json"
      "configs/lbm_dem/cases/fouling_3d_sweep_20260609_r1_rep_0005.json"
      "configs/lbm_dem/cases/fouling_3d_sweep_20260609_r2_rep_0010.json"
      "configs/lbm_dem/cases/fouling_3d_sweep_20260609_a1_attr_0015.json"
    )
    ;;
  *)
    echo "Usage: $0 [stage1|stage2|stage3|all]" >&2
    exit 2
    ;;
esac

for cfg in "${configs[@]}"; do
  echo "==> Running $cfg"
  PYTHONPATH=src python -m src.runners.run_lbm_dem --config "$cfg"
done
