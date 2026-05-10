# LBM-DEM Process Flow

This document describes the current processing flow for membrane-pore fouling
simulations.  It focuses on the standard production path using JSON
configuration files, periodic transverse boundaries, pressure inlet/outlet, and
fixed-cylinder pore geometry.

## Standard Entry Point

The recommended command path is:

```bash
poetry run python src/demos/run_lbm_dem.py \
  --config configs/lbm_dem/membrane_pressure_periodic_smoke.json
```

Additional CLI arguments may be added after `--config`.  Explicit CLI arguments
override values loaded from the JSON file.

## High-Level Flow

```text
User command
  -> JSON config in configs/lbm_dem/
  -> SimulationConfig
  -> src/demos/run_lbm_dem.py
  -> PoreGeometry
  -> FastLBMDEM / LBMDEMSolver
  -> DEMSolver
  -> time stepping
  -> metrics, stability checks, snapshots
  -> result directory under src/results/run_lbm_dem/
```

## Step-By-Step Flow

### 1. Case Definition

Reusable cases live in:

```text
configs/lbm_dem/
```

The preferred geometry format is:

```json
{
  "geometry": {
    "cylinders": [
      {"x": 34, "y": 13, "radius": 4.0},
      {"x": 54, "y": 25, "radius": 4.5}
    ]
  }
}
```

The older CLI-style `cylinder_spec` format is still accepted for compatibility.

### 2. Config Loading

`src/cfd_dem_lbm/simulation_config.py` loads the JSON file through
`SimulationConfig`.

Responsibilities:

- normalize dashed CLI-style keys to argparse destination names
- translate `geometry.cylinders` into `cylinder_spec` for runner compatibility
- preserve the source config for later reproducibility output
- write the final effective run arguments to `config.json`

### 3. Runner Setup

`src/demos/run_lbm_dem.py` is the current production runner.

It handles:

- CLI parsing
- JSON config defaults
- output directory creation
- particle count and source-fraction calculation
- solver construction
- warmup and main time stepping
- stability checks
- scalar metrics
- figures, ParaView files, CSV/NPZ, and metadata

The runner is intentionally being migrated toward smaller modules.  For now, it
remains the orchestration point.

### 4. Geometry Construction

`src/cfd_dem_lbm/geometry.py` owns pore geometry through `PoreGeometry`.

Responsibilities:

- validate fixed cylinders
- convert cylinders into solver tuples when needed
- create fixed-cylinder solid masks
- handle periodic-y cylinder masks
- estimate available water area
- define global and pore-local pressure sampling sections
- write cylinder geometry to VTK for ParaView

`run_lbm_dem.py` builds:

```text
GEOMETRY = PoreGeometry.from_cylinders(...)
```

and passes it into:

```text
FastLBMDEM(..., geometry=GEOMETRY)
```

### 5. Coupled Solver Initialization

`FastLBMDEM` extends `LBMDEMSolver`.

`LBMDEMSolver` currently acts as the coupled time-integration facade.  It owns:

- LBM distribution functions
- pressure inlet/outlet and transverse boundary modes
- fixed solid mask from `PoreGeometry`
- particle state arrays
- particle-fluid coupling mode selection
- Numba/NumPy backend selection
- shared `DEMSolver`

The fixed pore geometry is stored as:

```text
sim.geometry
sim.cylinders
sim.fixed_solid
```

### 6. Particle Solver

`src/cfd_dem_lbm/dem_solver.py` owns DEM contact mechanics.

It handles:

- particle-particle contact
- particle-wall contact when `y_boundary=wall`
- particle-cylinder contact
- rolling and sliding friction
- particle-particle attraction or repulsion
- particle-cylinder attraction or repulsion

The solver accesses the active coupled simulation state but keeps DEM load
calculation separated from the runner.

### 7. Boundary And Coupling Modes

For membrane-pore calculations, the standard boundary mode is:

```bash
--y-boundary periodic --streamwise-boundary pressure
```

Meaning:

- top and bottom are periodic, not no-slip
- left boundary is pressure inlet
- right boundary is pressure outlet
- fixed cylinders generate pore-scale pressure loss

Supported particle-fluid coupling modes:

| Mode | Role |
|---|---|
| `point_force` | Stokes drag and distributed reaction force |
| `solid_boundary` | particles are overlaid as moving solid cells |
| `immersed_boundary` | particles use IBM marker forces |

### 8. Time-Stepping Loop

The runner performs:

```text
warmup:
  sim.advance(warmup_steps)

main loop:
  sim.advance(snapshot_every)
  get rho, ux, uy
  compute pressure and speed
  check stability
  compute metrics
  store snapshot or final fields
  optionally write ParaView data
```

Stability checks currently include:

- non-finite `rho`, `ux`, `uy`, pressure, particle position, velocity, force
- fluid speed limit
- particle speed limit
- pressure magnitude limit

### 9. Metrics

The runner currently computes:

- inlet and outlet flux
- global pressure drop
- pore-local pressure drop
- normalized permeate flux
- pressure-drop ratio
- local pressure-drop ratio
- particle generated, active, pending, and passed counts
- particle-particle contact count
- wall contact count
- cylinder contact count
- effective inlet particle fraction
- estimated dimensionless groups

These calculations are a good candidate for future extraction into
`metrics.py`.

### 10. Output Directory

Each run writes to a timestamped directory under:

```text
src/results/run_lbm_dem/
```

The path includes geometry mode, interaction mode, rolling mode, particle source,
optional result tag, and timestamp.

Typical files:

| File | Purpose |
|---|---|
| `config.json` | source config and final effective arguments |
| `metadata.json` | solver and run metadata |
| `environment.json` | Python, platform, git, and accelerator details |
| `git_commit.txt` | commit hash used for the run |
| `run_status.json` | running, completed, or failed state |
| `analysis/time_series.csv` | scalar time series |
| `analysis/time_series.npz` | numeric time series arrays |
| `analysis/summary.json` | final comparison metrics |
| `analysis/summary.md` | human-readable summary |
| `lbm_dem_final.png` | final static visualization |
| `paraview/*.vtk`, `*.pvd` | ParaView fluid, particle, marker, and cylinder data |
| `lbm_dem_simulation.mp4` | optional movie |

### 11. Failure Handling

When instability is detected, the runner writes:

```text
run_status.json
```

with a failed status and the reason.  For unexpected exceptions, the global
exception hook also attempts to persist a failed `run_status.json` before Python
prints the traceback.

### 12. Verification Path

The test suite currently checks:

- geometry mask generation
- periodic-y particle wrapping
- periodic-y particle-fluid interpolation and force spreading
- IBM marker wrapping
- pressure inlet/outlet with periodic-y boundary
- fluid method combinations
- particle method combinations
- coupling method combinations
- Numpy and Numba execution paths

Use:

```bash
poetry run python -m pytest tests
```

## Current Architecture Boundary

The current system is partly layered but not fully separated yet.

Already separated:

- `SimulationConfig`
- `PoreGeometry`
- `DEMSolver`
- `FastLBMDEM`

Still mixed inside `run_lbm_dem.py`:

- result writing
- metric calculations
- stability monitoring
- plotting and animation
- ParaView export orchestration

Still mixed inside `lbm_dem.py`:

- LBM core
- boundary conditions
- coupling models
- IBM helpers
- Numba kernels

## Recommended Next Refactors

1. Extract metrics from `run_lbm_dem.py` into `src/cfd_dem_lbm/metrics.py`.
2. Extract output writing into `src/cfd_dem_lbm/result_writer.py`.
3. Extract stability checks into `src/cfd_dem_lbm/stability.py`.
4. Extract boundary condition logic from `lbm_dem.py` into `boundaries.py`.
5. Extract coupling implementations into `coupling.py`.
6. Move run/sweep shell scripts out of `src/bin` into a top-level `scripts/`
   hierarchy when the command surface stabilizes.
