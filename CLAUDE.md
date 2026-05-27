# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Install dependencies
poetry install

# Run all tests
poetry run python -m pytest tests

# Run tests excluding slow ones
poetry run python -m pytest tests -m "not slow"

# Run a single test file
poetry run python -m pytest tests/test_lbm_dem.py

# Lint
poetry run ruff check src tests

# Format
poetry run black src tests

# Run a simulation
poetry run python src/runners/run_lbm_dem.py \
  --config configs/lbm_dem/cases/fouling_four_cylinder_supply.json

# Run DEM packing
poetry run python src/runners/run_dem_packing.py \
  --config configs/dem_packing/<config>.json

# Run a parameter sweep
poetry run python src/runners/run_lbm_dem_sweep.py \
  --sweep configs/lbm_dem/sweeps/fouling_screening_example.json

# Analyse sweep results
poetry run python src/tools/analyze_lbm_dem_design_sweeps.py <results_dir>

# Benchmark accelerator backends
poetry run python src/tools/benchmark_lbm_accelerators.py
```

## Architecture

This is a 2D particulate-flow simulation toolkit coupling a **Lattice-Boltzmann Method (LBM)** fluid solver with a **Discrete Element Method (DEM)** particle solver, designed for membrane fouling analysis.

### Solver Hierarchy

```
FastLBMDEM (src/particulate_flow/fast_solver.py)
  └── LBMDEMSolver (src/particulate_flow/lbm_dem.py)   ← core 9,500-line solver
        ├── D2Q9 LBM fluid (BGK or TRT collision + Guo forcing)
        ├── DEMSolver (src/particulate_flow/dem_solver.py)
        └── PoreGeometry (src/particulate_flow/geometry.py)
```

`FastLBMDEM` adds per-timestep caching of macroscopic fields (ρ, ux, uy) on top of `LBMDEMSolver`. The runners in `src/runners/` instantiate `FastLBMDEM` directly.

### Two-Way Coupling

- **Fluid → Particles:** Stokes drag from interpolated fluid velocity at particle centers
- **Particles → Fluid:** Newton 3rd law body force distributed back to lattice nodes

### Configuration System

JSON configs in `configs/lbm_dem/cases/` support `"extends"` inheritance. `SimulationConfig` (`src/particulate_flow/io/config.py`) loads and flattens sections (domain, flow, solver, accelerator, numerics, particles, physics, runtime, outputs, stability) and binds them to CLI argument parsers. CLI flags override config values.

Config directory layout:
```
configs/
  lbm_dem/
    cases/       — runnable case configs (may use extends)
    geometries/  — geometry fragments (cylinder arrays)
    materials/   — material property fragments
    templates/   — base templates for common setups
    sweeps/      — parameter sweep definitions
  dem_packing/
    cases/       — DEM-only packing cases
```

### Numba Acceleration

Hotspots in `dem_solver.py` and the LBM collision operator have optional Numba JIT paths. They fall back gracefully to pure NumPy. Controlled at runtime via `--fluid-accelerator` and `--compute-accelerator` flags.

### Outputs

Results land in `src/results/run_lbm_dem/YYYYMMDD_HHMMSS/`:
- `run_status.json`, `summary.json/md` — metadata and metrics
- `*.csv` / `*.npz` — timeseries (pressure, velocity, particle positions)
- `fields_*.vtk`, `particles_*.vtk` — ParaView visualization files
- `snapshot_*.png`, `*.mp4` — visual outputs (ffmpeg required for video)

### Key Docs

Detailed design notes live in `docs/fouling_model/`:
- `SYSTEM_ARCHITECTURE.md` — high-level design
- `ALGORITHMS.md` — numerical methods
- `LBM_DEM_COUPLING_LIMITATIONS.md` — known constraints
