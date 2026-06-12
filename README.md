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

## Multi-Machine Sweeps

Manually staged case sweeps can be split across several machines, using git as
the coordination channel. Each machine pushes lightweight result artifacts
(summary JSON, time-series CSV — never VTK/video) to `sweep_results/`, so
`./bin/sweep status` on any machine shows cross-machine progress.

```bash
# On any machine, after git clone/pull:
./bin/sweep list                       # registered sweeps and case aliases
./bin/sweep status                     # who finished what, where
./bin/sweep run fouling_3d_20260609 a1 a2 a3   # run this machine's assigned cases
```

`run` handles everything: `git pull`, `poetry install`, sequential execution,
result harvesting, and `git push` after each case. Already-completed cases are
skipped (override with `--force`). Case assignment is manual — decide per
machine which aliases to run, and use `status` to avoid double assignment.

To register a new sweep, add a manifest JSON listing case configs with short
aliases under `configs/lbm_dem/sweeps/manifests/` (see `smoke_test.json` for
the template). Case configs must set `outputs.result_tag` in the leaf JSON.

### Web UI

A lightweight web dashboard wraps the same sweep machinery so you can monitor
progress, start/stop runs, and find result paths without typing commands.

```bash
poetry install --with web   # one-time: pulls FastAPI/uvicorn (optional group)
./bin/sweep-web             # serves http://127.0.0.1:8765
```

What it does:

- **Status table** — every case with a color-coded status badge, machine, and
  timestamp. The view is multi-machine (read from the git-synced `sweep_results/`).
- **Live progress** — for runs launched on *this* machine, a progress bar derived
  from snapshot-file count (updates each `snapshot_every` interval; the solver is
  untouched).
- **Run / Stop** — launches `./bin/sweep run` in the background and can terminate
  it. Control is local to the machine running the UI.
- **Results** — result directory path (with copy/reveal buttons), key metrics from
  `summary.json`, and a time-series sparkline.

Host/port override with `SWEEP_WEB_HOST` / `SWEEP_WEB_PORT`. Runtime logs and the
local process registry live in the gitignored `.sweep_web/` directory.

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
