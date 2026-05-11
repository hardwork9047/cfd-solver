# Current System Architecture

This document summarizes the current architecture of the LBM-DEM membrane fouling simulation system.

For the execution sequence from config file to result artifacts, see
`docs/fouling_model/PROCESS_FLOW.md`.
For the numerical algorithms used inside the solver, see
`docs/fouling_model/ALGORITHMS.md`.

## Purpose

The current system is designed to evaluate how suspended particles are transported, retained, attached, repelled, or passed through idealized membrane-pore geometries represented by fixed cylinders. The main target is systematic membrane fouling analysis using repeatable simulations, parameter sweeps, time-series metrics, and visualization outputs.

## High-Level Flow

```text
User command or sweep script
  -> src/runners/run_lbm_dem.py
  -> cfd_dem_lbm.FastLBMDEM
  -> LBM fluid solver + DEM particle solver + particle-fluid coupling
  -> src/results/run_lbm_dem/.../run_TIMESTAMP/
  -> metadata, run_status, analysis CSV/NPZ, summary JSON/MD, images, videos, ParaView data
  -> plotting / summarization / benchmark scripts
```

The system is currently command-line based. Web UI or job-queue execution can be added on top by generating command arguments or config files and monitoring `run_status.json` and `summary.json`.

## Main Directories

| Path | Role |
|---|---|
| `src/cfd_dem_lbm/` | Production LBM-DEM solver implementation |
| `src/runners/run_lbm_dem.py` | Main simulation runner for fouling calculations |
| `src/bin/` | Reproducible scripts for runs, sweeps, verification, benchmarks, and plotting |
| `src/results/` | Simulation, verification, benchmark, and plotting outputs |
| `tests/` | Unit and regression tests |
| `docs/fouling_model/` | Design notes, limitations, workflow, and architecture docs |
| `paper/` | Paper draft and manuscript material |

## Core Solver Components

### Fluid Solver

The fluid path is implemented in `src/cfd_dem_lbm/lbm_dem.py` and exposed through `FastLBMDEM`.

Supported fluid methods:

| Option | Meaning |
|---|---|
| `lbm-bgk-guo` | BGK LBM with Guo forcing |
| `lbm-trt-guo` | TRT LBM with Guo forcing |

Supported fluid accelerators:

| Option | Meaning |
|---|---|
| `numpy` | Use NumPy implementation |
| `numba` | Require Numba-compiled implementation |
| `auto` | Use Numba when available, otherwise fall back to NumPy |

Flow control modes:

| Option | Meaning |
|---|---|
| `fixed-pressure` / `constant-pressure` | Fixed body-force style pressure-driven flow |
| `target-max-velocity` | Adjust drive force to target the maximum velocity |
| `constant-flux` | Alias for the current target max velocity control proxy |

The Reynolds number in `run_lbm_dem.py` is defined using particle diameter as representative length and maximum flow velocity as representative velocity.

## Run Configuration And Reproducibility

Simulation cases can be defined either with CLI arguments or with JSON config
files.  The preferred research workflow is to keep reusable cases under:

```text
configs/lbm_dem/
```

For example:

```bash
poetry run python src/runners/run_lbm_dem.py \
  --config configs/lbm_dem/membrane_pressure_periodic_smoke.json
```

Explicit CLI arguments override values loaded from the config file.  Every run
writes a small reproducibility bundle into its result directory:

| File | Purpose |
|---|---|
| `config.json` | Source config plus final effective arguments |
| `metadata.json` | Solver, physics, output, and command metadata |
| `environment.json` | Python, platform, git, and accelerator information |
| `git_commit.txt` | Commit hash used for the run |
| `run_status.json` | Running/completed/failed state for sweep recovery |

The current architecture intentionally keeps run tracking file-based rather than
using SQLite.  This keeps each result folder self-contained and easy to move,
archive, or compare on another machine.

## Maintainable Module Layout

The production path is being moved toward a layered structure:

```text
configs/lbm_dem/*.json
  -> SimulationConfig
  -> run_lbm_dem.py / Runner
  -> FastLBMDEM / LBMDEMSolver
  -> DEMSolver, PoreGeometry, coupling and boundary helpers
  -> Result artifacts under src/results/run_lbm_dem/
```

The intended responsibility split is:

| Layer | Current module | Responsibility |
|---|---|---|
| Case definition | `configs/lbm_dem/` | Reusable geometry, physics, output, and backend settings |
| Config loading | `src/cfd_dem_lbm/simulation_config.py` | Convert JSON config files into runner arguments |
| Geometry | `src/cfd_dem_lbm/geometry.py` | Cylinder definitions, pore masks, pressure probe sections, cylinder VTK |
| Runner | `src/runners/run_lbm_dem.py` | CLI, simulation execution, output orchestration |
| Coupled solver | `src/cfd_dem_lbm/lbm_dem.py` | LBM-DEM time integration facade |
| Particle solver | `src/cfd_dem_lbm/dem_solver.py` | DEM contact, wall/cylinder loads, surface interactions |
| Fast path | `src/cfd_dem_lbm/fast_solver.py` | Cached solver variant used by production runs |

New pore layouts should be implemented in `geometry.py` and referenced from
JSON configs.  Legacy `--cylinder-spec X Y R` remains supported, but the
preferred config style is:

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

Boundary modes:

| Option | Meaning |
|---|---|
| `--y-boundary wall` | No-slip top and bottom walls |
| `--y-boundary periodic` | Periodic transverse boundary for membrane-unit-cell calculations |
| `--streamwise-boundary periodic-force` | Periodic x streaming with body-force pressure-gradient approximation |
| `--streamwise-boundary pressure` | Left pressure inlet and right pressure outlet using density boundary conditions |

For the membrane-pore fouling model, the standard boundary setting is:

```bash
--y-boundary periodic --streamwise-boundary pressure
```

This represents a periodically repeated transverse unit cell.  The top and
bottom boundaries are not no-slip walls in this mode.  The pressure drop is
imposed between the left inlet and right outlet, and the fixed cylinders inside
the domain create the pore-scale hydraulic resistance and local pressure loss.

For pressure inlet/outlet runs, the pressure difference is specified with:

```bash
--pressure-drop 1e-4 --rho-out 1.0
```

Internally, `p = cs^2 rho`, so the inlet density is computed from the requested pressure drop. This mode is intended for filtration-style calculations where left and right should represent upstream and downstream reservoirs. The older `periodic-force` mode remains useful for fast periodic-channel verification and screening, but it should not be used as the main interpretation mode for pore-scale pressure-loss studies.

## DEM Particle Solver

Particle mechanics are handled by the shared DEM solver in:

```text
src/cfd_dem_lbm/dem_solver.py
```

Supported particle contact methods:

| Option | Meaning |
|---|---|
| `dem-hertz` | Hertzian normal contact |
| `dem-linear` | Linear normal contact |

Supported particle search methods:

| Option | Meaning |
|---|---|
| `cell_list` | Cell-list neighbor search for particle pairs |
| `all_pairs` | Brute-force all-pair search |

Supported particle interactions:

| Feature | Current role |
|---|---|
| Particle-particle contact | DEM contact force and damping |
| Particle-wall contact | Boundary collision and friction |
| Particle-cylinder contact | Contact with fixed membrane-pore cylinders |
| Particle-particle Hamaker attraction/repulsion | Enabled by `--particle-attraction` or `--particle-repulsion` |
| Particle-cylinder Hamaker attraction/repulsion | Uses the same attraction/repulsion mode and strength parameters |
| Rolling friction | Enabled by `--rolling-friction` |

## Particle-Fluid Coupling

The coupling mode is selected with:

```bash
--particle-fluid-coupling point_force
--particle-fluid-coupling immersed_boundary
--particle-fluid-coupling solid_boundary
```

| Coupling | Meaning | Use case |
|---|---|---|
| `point_force` | Stokes drag and point-force feedback without geometric blockage | Fast screening and dilute particle effects |
| `immersed_boundary` | IBM marker/direct-forcing representation of particle surfaces | Smoother particle-fluid boundary approximation |
| `solid_boundary` | Particles are overlaid onto the LBM solid mask | Strong blockage and clogging upper-bound behavior |

The IBM marker spacing is controlled by:

```bash
--ibm-marker-spacing
```

The current marker count is approximately based on particle circumference divided by marker spacing.

## Membrane-Pore Geometry

Fixed cylinders represent pore obstacles or membrane microstructure. They can be specified as:

```bash
--cylinder
--cyl-x 45 --cyl-y 35 --cyl-r 6
```

or repeated explicit cylinders:

```bash
--cylinder-spec X Y R --cylinder-spec X Y R ...
```

Each cylinder is treated as a fixed solid boundary in the LBM field and as a DEM interaction object for contact and Hamaker-type surface forces.

## Particle Supply

Two particle supply modes exist:

| Option | Meaning |
|---|---|
| `--particle-source initial` | Place particles at initialization |
| `--particle-source left-inlet` | Gradually inject particles from the left inlet |

For left-inlet supply, the solver tracks:

| Quantity | Meaning |
|---|---|
| `last_inlet_flow_rate` | Latest inlet flow rate |
| `cumulative_inlet_flow_area` | Integrated inlet flow area over time |
| `injected_particle_area` | Total area of injected particles |
| `inlet_particle_area_budget` | Remaining particle area budget |
| `effective_inlet_particle_fraction` | `injected_particle_area / cumulative_inlet_flow_area` |
| `inlet_particle_fraction_error` | Difference from requested source fraction |

This makes it possible to verify whether the injected particle concentration matches the requested volume or area fraction.

## Main Runner

The primary entry point is:

```bash
poetry run python src/runners/run_lbm_dem.py [options]
```

This runner handles:

- command-line configuration
- output directory creation with timestamps
- solver initialization
- warmup
- main time stepping
- instability detection
- snapshot generation
- ParaView export
- final image and optional MP4 output
- time-series analysis output
- final summary output
- run status output

## Output Structure

Each run is saved under:

```text
src/results/run_lbm_dem/<case_tags>/run_YYYYMMDD_HHMMSS_xxxxxx/
```

Important files:

| File | Role |
|---|---|
| `metadata.json` | Reproducibility metadata, command, git commit, solver configuration |
| `run_status.json` | Current or final status: `running`, `completed`, `failed`, `failed_unstable` |
| `analysis/time_series.csv` | Scalar time-series metrics |
| `analysis/time_series.npz` | Same time-series data in NumPy format |
| `analysis/summary.json` | Final and extrema metrics for comparisons and dataset generation |
| `analysis/summary.md` | Human-readable summary |
| `lbm_dem_final.png` | Final snapshot image |
| `lbm_dem_simulation.mp4` | Optional movie |
| `paraview/*.vtk`, `*.pvd` | Fluid, particles, cylinders, and optional IBM markers for ParaView |

## Analysis Metrics

The main time-series output includes:

| Metric group | Examples |
|---|---|
| Flow | inlet flux, outlet flux, permeate flux, normalized permeate flux |
| Pressure | inlet pressure, outlet pressure, global pressure drop, local pore pressure drop |
| Particles | active, generated, pending, passed particles |
| Inlet supply | cumulative inlet flow, injected particle area, effective inlet particle fraction |
| Contacts | particle-particle, wall, cylinder, total contacts |
| Fouling | retained particle ratio, passed particle ratio, fouling resistance index |
| Dimensionless groups | observed Reynolds number, particle Reynolds number, Stokes estimate |
| Solid/cake approximation | dynamic solid fraction, porous resistance fraction |

Global pressure drop is measured between near-inlet and near-outlet sections. Local pore pressure drop is measured around the fixed cylinder or cylinder group.

## Stability Detection

`run_lbm_dem.py` now checks for:

- NaN or Inf in fluid fields
- NaN or Inf in particle positions, velocities, or forces
- excessive fluid speed
- excessive particle speed
- excessive pressure magnitude

Controls:

```bash
--max-stable-speed
--max-stable-pressure
```

If instability is detected, the run is stopped and `run_status.json` records `failed_unstable` with the reason.

## Verification

The main verification entry point is:

```bash
poetry run python src/bin/verify_lbm_dem_fluid.py
```

It verifies the production solver path `cfd_dem_lbm.FastLBMDEM`, including:

- plane Poiseuille profile
- streamwise flux balance
- target max velocity control
- solid-boundary particle flux reduction
- near-clogged solid-boundary leakage suppression
- immersed-boundary particle load
- porous resistance flux reduction
- particle-cylinder Hamaker attraction

Verification reports are written under:

```text
src/results/verify_lbm_dem_fluid/
```

## Benchmarks

Important benchmark scripts:

| Script | Purpose |
|---|---|
| `src/bin/benchmark_lbm_accelerators.py` | Fluid-only NumPy vs Numba benchmark |
| `src/bin/benchmark_numba_components.py` | Component-level DEM, solid mask, IBM benchmark |
| `src/bin/benchmark_lbm_dem_coupling_backends.py` | Coupling mode and backend comparison |
| `src/bin/benchmark_ibm_marker_spacing.py` | IBM marker spacing convergence benchmark |

The current practical recommendation is to use Numba-capable paths for production sweeps, typically:

```bash
--fluid-accelerator auto --compute-accelerator auto
```

or, when explicitly fixing the fluid path:

```bash
--fluid-accelerator numba --compute-accelerator auto
```

## Plotting and Dataset Generation

Fouling metrics can be plotted from existing run directories:

```bash
poetry run python src/bin/plot_lbm_dem_fouling_metrics.py src/results/run_lbm_dem
```

This generates:

- metric time-series plots
- final pressure-drop vs retention scatter plot
- `fouling_metrics_dataset.csv`
- plotting report JSON

This script is the current bridge from raw simulation outputs toward surrogate modeling datasets.

## Sweep Execution

Sweep scripts and config-driven sweeps are in `src/bin/`, including:

- `run_lbm_dem_config_sweep.py`
- `run_lbm_dem_16_pore_geometry_parallel.sh`
- `run_lbm_dem_100_design_sweep_parallel.sh`

The current sweep model is script-based. Failed or unstable runs should be identified by reading each run's `run_status.json`.

## Web Application Boundary

A future Web app should not run heavy simulations inside the web request process. The recommended boundary is:

```text
Web UI
  -> create job configuration
  -> launch worker/subprocess
  -> worker runs src/runners/run_lbm_dem.py
  -> results are written under src/results/
  -> Web UI reads run_status.json, summary.json, time_series.csv, images, and videos
```

The existing output files already support this architecture:

- `run_status.json` for job state
- `metadata.json` for reproducibility
- `analysis/summary.json` for result cards and tables
- `analysis/time_series.csv` for plots
- `lbm_dem_final.png` and MP4 for visual inspection
- ParaView files for detailed post-processing

## Current Limitations

- The solver is currently 2-D.
- Particle-fluid coupling modes have different physical meanings and should not be compared without noting the coupling model.
- `constant-flux` is currently a target-max-velocity style proxy, not a strict integral flux controller.
- IBM does not create hard blocked cells; it applies marker/direct-forcing constraints.
- `solid_boundary` is better for hard blockage checks but can be more stair-stepped geometrically.
- Porous resistance is a simplified Brinkman-like local drag approximation, not a full cake-layer model.
- Web app, database-backed job queue, and resumable scheduler are not yet implemented.

## Recommended Next Architecture Step

The next architectural improvement should be a lightweight job layer:

1. Define a run configuration schema.
2. Add a small job launcher that writes one job directory per condition.
3. Store job state by reading and updating `run_status.json`.
4. Let a Web UI or CLI submit jobs through the same interface.
5. Keep all numerical results in the existing `src/results/` structure.

This would make the future Web app mostly a control and visualization layer, while preserving the current reproducible command-line workflow.
