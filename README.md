# Particulate Flow Simulator

This repository contains a two-dimensional particulate-flow simulation toolkit
for studying particle transport, deposition, pore-scale blockage, and fouling
around configurable solid obstacles.

The active code path is centered on `src/particulate_flow/` and the production runner
`src/runners/run_lbm_dem.py`.

## Active Capabilities

| Component | Path | Purpose |
|---|---|---|
| Coupled solver | `src/particulate_flow/lbm_dem.py` | LBM-DEM time integration |
| Fast solver facade | `src/particulate_flow/fast_solver.py` | Cached production solver path |
| DEM solver | `src/particulate_flow/dem_solver.py` | Particle contact, rolling friction, surface forces |
| Geometry layer | `src/particulate_flow/geometry.py` | Cylinder pore layouts and masks |
| Config loader | `src/particulate_flow/simulation_config.py` | JSON-based run configuration |
| Fluid verification | `src/particulate_flow/fluid_verification.py` | Shared-solver verification problems |
| Cylinder flow runner | `src/particulate_flow/cylinder_flow.py` | Fluid-only cylinder calculations through the same solver path |
| Main runner | `src/runners/run_lbm_dem.py` | Simulation execution and output orchestration |

## Standard Run

```bash
poetry run python src/runners/run_lbm_dem.py \
  --config configs/lbm_dem/cases/fouling_four_cylinder_supply.json
```

The recommended membrane-pore boundary setting is:

```bash
--y-boundary periodic --streamwise-boundary pressure
```

This represents a periodically repeated transverse unit cell with a pressure
inlet on the left, a pressure outlet on the right, and fixed cylinders creating
the pore-scale hydraulic resistance.

## Outputs

Runs write timestamped result folders under:

```text
src/results/run_lbm_dem/
```

Typical outputs include:

- reproducibility metadata
- time-series CSV/NPZ data
- pressure and velocity fields
- particle and cylinder ParaView VTK files
- snapshots and videos when enabled
- summary JSON/Markdown reports

## Documentation

The current model is documented under:

```text
docs/fouling_model/
```

Key files:

- `SYSTEM_ARCHITECTURE.md`
- `PROCESS_FLOW.md`
- `ALGORITHMS.md`
- `NUMERICAL_WORKFLOW.md`
- `LBM_DEM_COUPLING_LIMITATIONS.md`

## Project Structure

```text
configs/lbm_dem/      Reusable JSON run configurations
docs/fouling_model/   Architecture, algorithms, limitations, and workflow notes
paper/                Manuscript drafts
src/bin/              Reproducible run, sweep, plotting, and benchmark scripts
src/particulate_flow/      Active LBM-DEM solver package
src/runners/            Main executable runner
src/results/          Generated outputs grouped by program name
tests/                pytest suite for the active solver path
```

## Development

```bash
poetry install
poetry run python -m pytest tests
```

The repository intentionally keeps simulation outputs file-based and
self-contained. Result directories can be moved or archived without requiring a
database.
