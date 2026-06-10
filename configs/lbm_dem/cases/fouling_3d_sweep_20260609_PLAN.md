# 3D Fouling Sweep Plan (2026-06-09)

Base reference: `simulate_20260601_095407`

Purpose:
- Reproduce the fully fouling baseline.
- Find the attraction-strength threshold between progressive fouling and near-complete blockage.
- Separate the effects of concentration, rolling friction, and membrane throat width.
- Check whether repulsion is strong enough to suppress cake growth.

Recommended execution order:
1. `configs/lbm_dem/cases/fouling_3d_sweep_20260609_b0_baseline.json`
2. `configs/lbm_dem/cases/fouling_3d_sweep_20260609_a2_attr_0010.json`
3. `configs/lbm_dem/cases/fouling_3d_sweep_20260609_a3_attr_0005.json`
4. `configs/lbm_dem/cases/fouling_3d_sweep_20260609_c1_phi_010.json`
5. `configs/lbm_dem/cases/fouling_3d_sweep_20260609_c2_phi_020.json`
6. `configs/lbm_dem/cases/fouling_3d_sweep_20260609_f1_roll_005.json`
7. `configs/lbm_dem/cases/fouling_3d_sweep_20260609_f2_roll_030.json`
8. `configs/lbm_dem/cases/fouling_3d_sweep_20260609_g1_gap_wide.json`
9. `configs/lbm_dem/cases/fouling_3d_sweep_20260609_g2_gap_narrow.json`
10. `configs/lbm_dem/cases/fouling_3d_sweep_20260609_r1_rep_0005.json`
11. `configs/lbm_dem/cases/fouling_3d_sweep_20260609_r2_rep_0010.json`
12. `configs/lbm_dem/cases/fouling_3d_sweep_20260609_a1_attr_0015.json`

Case summary:

| Case | Focus | Key change from baseline |
| --- | --- | --- |
| `b0_baseline` | Reproduction | `A=0.002`, `phi=0.15`, `mu_r=0.15`, baseline geometry |
| `a1_attr_0015` | Attraction | `A=0.0015` |
| `a2_attr_0010` | Attraction | `A=0.0010` |
| `a3_attr_0005` | Attraction | `A=0.0005` |
| `r1_rep_0005` | Repulsion | switch to repulsion, `R=0.0005` |
| `r2_rep_0010` | Repulsion | switch to repulsion, `R=0.0010` |
| `c1_phi_010` | Concentration | `phi=0.10` with `A=0.0010` |
| `c2_phi_020` | Concentration | `phi=0.20` with `A=0.0010` |
| `f1_roll_005` | Rolling friction | `mu_r=0.05` with `A=0.0010` |
| `f2_roll_030` | Rolling friction | `mu_r=0.30` with `A=0.0010` |
| `g1_gap_wide` | Geometry | widen throat: x = `19, 29` |
| `g2_gap_narrow` | Geometry | narrow throat: x = `21, 27` |

One-by-one run example:

```bash
python -m src.runners.run_lbm_dem --config configs/lbm_dem/cases/fouling_3d_sweep_20260609_a2_attr_0010.json
```

Batch run example:

```bash
for cfg in \
  configs/lbm_dem/cases/fouling_3d_sweep_20260609_b0_baseline.json \
  configs/lbm_dem/cases/fouling_3d_sweep_20260609_a2_attr_0010.json \
  configs/lbm_dem/cases/fouling_3d_sweep_20260609_a3_attr_0005.json \
  configs/lbm_dem/cases/fouling_3d_sweep_20260609_c1_phi_010.json \
  configs/lbm_dem/cases/fouling_3d_sweep_20260609_c2_phi_020.json \
  configs/lbm_dem/cases/fouling_3d_sweep_20260609_f1_roll_005.json \
  configs/lbm_dem/cases/fouling_3d_sweep_20260609_f2_roll_030.json \
  configs/lbm_dem/cases/fouling_3d_sweep_20260609_g1_gap_wide.json \
  configs/lbm_dem/cases/fouling_3d_sweep_20260609_g2_gap_narrow.json \
  configs/lbm_dem/cases/fouling_3d_sweep_20260609_r1_rep_0005.json \
  configs/lbm_dem/cases/fouling_3d_sweep_20260609_r2_rep_0010.json \
  configs/lbm_dem/cases/fouling_3d_sweep_20260609_a1_attr_0015.json
do
  python -m src.runners.run_lbm_dem --config "$cfg"
done
```

Staged helper script:

```bash
chmod +x bin/run-fouling-3d-sweep-20260609.sh
bin/run-fouling-3d-sweep-20260609.sh stage1
```

Primary comparison metrics to track:
- `generated_particles`
- `removed_particles`
- final `particle_count`
- `u_max` drop from first to last snapshot
- time to visible bridge formation in the pore throat
