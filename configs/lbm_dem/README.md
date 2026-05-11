# LBM-DEM Configuration Layout

Run conditions are managed as JSON files instead of one-off shell scripts.

```text
configs/lbm_dem/
  templates/   reusable calculation-type defaults
  cases/       named single-run studies
  sweeps/      multi-case or factor-sweep definitions
  geometries/  reusable cylinder layouts for copy/reference
  materials/   reusable particle/surface-interaction settings for copy/reference
```

Single cases are launched with:

```bash
poetry run python src/runners/run_lbm_dem.py \
  --config configs/lbm_dem/cases/fouling_four_cylinder_supply.json
```

Sweeps are launched with:

```bash
poetry run python src/bin/run_lbm_dem_config_sweep.py \
  configs/lbm_dem/sweeps/fouling_screening_example.json --dry-run
```

Config files may use research-friendly sections such as `domain`, `flow`,
`particles`, `physics`, `runtime`, `outputs`, and `stability`. The loader
flattens these sections into runner arguments. A case can reuse a template with
`extends`.

