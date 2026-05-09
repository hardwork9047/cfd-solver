# LBM-DEM Fouling Workflow

This note records the current practical workflow for developing the membrane
fouling model.

## 1. Particle-fluid coupling modes

`run_lbm_dem.py` supports two coupling modes.

- `point_force`: the previous two-way Stokes-drag back reaction.  Particles
  exchange momentum with the fluid but do not geometrically block lattice nodes.
- `solid_boundary`: active DEM particles are overlaid on the LBM `solid` mask.
  This makes deposited particles act as moving no-slip obstacles for the fluid
  solver.  It is still a 2-D approximation, but it is the preferred mode when
  testing pore blockage and flux decline.

Example:

```bash
python src/demos/run_lbm_dem.py \
  --particle-fluid-coupling solid_boundary \
  --fluid-method lbm-trt-guo \
  --particle-method dem-hertz
```

## 2. Verification

Use the production solver path for verification:

```bash
src/bin/verify_lbm_dem_fluid.py
```

The current suite checks:

- plane Poiseuille velocity-profile shape,
- flux balance around a fixed cylinder,
- target maximum-velocity control,
- flux reduction when a DEM particle is imposed as a solid boundary.

These are not final validation of membrane fouling, but they catch the most
important numerical failures before running large sweeps.

## 3. Standard fouling metrics

Each `run_lbm_dem.py` run writes `analysis/time_series.csv` and
`analysis/time_series.npz`.  The standardized fouling metrics include:

- `normalized_permeate_flux`,
- `pressure_drop_ratio`,
- `fouling_resistance_index`,
- `passed_particle_ratio`,
- `retained_particle_ratio`,
- `particle_area_fraction`,
- `dynamic_particle_solid_fraction`,
- `total_contacts`,
- `particle_particle_contacts`,
- `cylinder_contacts`.

For screening, lower `fouling_resistance_index`, lower retained-particle ratio,
and higher normalized permeate flux are favorable.  Visual inspection and
ParaView fields are still needed before treating a case as physically meaningful.

## 4. Config-file sweeps

Use JSON configs to keep parameter sweeps reproducible:

```bash
src/bin/run_lbm_dem_config_sweep.py src/bin/lbm_dem_sweep_example.json --dry-run
src/bin/run_lbm_dem_config_sweep.py src/bin/lbm_dem_sweep_example.json
```

The runner writes logs and `status.tsv` under
`src/results/run_lbm_dem/config_sweeps/<timestamp>/`.
