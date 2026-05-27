# CHANGELOG

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [0.6.0] - 2026-05-27

### Added

- `src/tools/`: 8 analysis and benchmark scripts moved from the abolished `src/bin/`, see PR #6.
- `src/runners/run_lbm_dem_sweep.py`: sweep runner promoted from `src/bin/`.

### Changed

- `src/bin/` abolished — all scripts redistributed to `src/runners/` (sweep runner) and `src/tools/` (analysis/benchmark).
- `configs/lbm_dem/cases/membrane_pressure_periodic_smoke.json` moved from `configs/lbm_dem/` root.
- `CLAUDE.md`: updated paths, added config directory layout reference, added sweep/analysis command examples.

---

## [0.5.0] - 2026-05-27

### Added

- `particulate_flow/builder.py`: `build_lbm_dem_solver()` and `build_dem_packing_solver()` centralise solver construction previously inline in runners, see PR #5.
- `SimulationConfig` now accepts `solver` and `accelerator` JSON sections. `accelerator.fluid` maps to `fluid_accelerator`; `accelerator.compute` maps to `compute_accelerator`. The `numerics` section remains supported for backward compatibility.

### Changed

- `run_lbm_dem.py` and `run_dem_packing.py` are now thin entry-points: solver construction is delegated to the builder functions.
- Template configs (`fouling_supply`, `cylinder_flow`, `dem_settling_pack`) updated to use `solver`/`accelerator` sections.

---

## [0.4.0] - 2026-05-27

### Changed

- Restructured `particulate_flow` into subpackages: `lbm/`, `dem/`, `ibm/`, `geometry/`, `io/` for token-efficient AI navigation, see PR #4.
- Extracted D2Q9 constants and Numba kernels into dedicated modules (`lbm/constants.py`, `lbm/kernels.py`, `dem/kernels.py`, `ibm/kernels.py`).
- Moved visualization helpers to `io/visualization.py`; `plot_fields` / `plot_particles` remain importable from the top-level package.
- Migrated tests to mirror the subpackage structure (`tests/lbm/`, `tests/dem/`, `tests/ibm/`, `tests/geometry/`, `tests/io/`).
- Added `pytest.ini_options` with `importmode=importlib` to avoid stdlib name collisions in test subdirectories.

---

## [0.3.0] - 2026-03-03

### Added
- **CylinderFlow solver** (`src/cfd/cylinder.py`): Incompressible flow past a
  circular cylinder using the fractional-step projection method. Supports
  detection of Kármán vortex shedding at Re > 47.
- **2nd-order linear-upwind advection** (`advection_scheme="upwind2"`): O(Δx²)
  one-sided stencils with near-zero numerical viscosity, enabling Re_eff ≈
  Re_physical on coarser grids. Successfully reproduces Kármán shedding at
  Re = 100 (St ≈ 0.10) and Re = 200 (St ≈ 0.13).
- **D2Q9 Lattice-Boltzmann solver** (`cfd_lbm` package): BGK collision operator,
  Zou-He moving-lid boundary condition, suitable for lid-driven cavity and
  Poiseuille channel benchmarks.
- **DEM particle system** (`dem` package): Discrete Element Method with Hertzian
  normal contact, velocity-Verlet time integration, and spatial-grid collision
  detection for O(N) performance.
- **Coupled LBM-DEM solver** (`particulate_flow` package): Bilinear interpolation
  (fluid → particles) and 4-point force distribution (particles → fluid) for
  particle-laden flow simulation.
- **GitHub Actions CI** (`.github/workflows/ci.yml`): Automated tests on Python
  3.11, 3.12, 3.13 matrix; lint and format checks; package build verification.
- **SECURITY.md**: Vulnerability disclosure policy and supported versions.
- **CODE_OF_CONDUCT.md**: Contributor Covenant 2.1.
- **Unit tests for `cfd_lbm`** (`tests/test_lbm.py`): 13 tests covering
  equilibrium distribution, mass conservation, boundary conditions, and solver
  stability.
- **Unit tests for `dem`** (`tests/test_dem.py`): 11 tests covering particle
  initialization, gravity, contact forces, and domain confinement.

### Changed
- **CavityFlow projection method**: Rewrote `step()` as a proper fractional-step
  method (predict u* → solve ∇²p = +ρ/dt·∇·u* → correct u). Fixed RHS sign
  (was negative, now positive) and SOR update sign (was +=, now -=).
- **CavityFlow**: Added `advection_scheme` parameter (`"upwind"` / `"upwind2"`).
- **CylinderFlow**: Unified advection dispatch via `_advect()` — eliminates
  intermediate `phi_new` allocations in `step()`.
- **pyproject.toml**: Now exports all 4 packages (`cfd`, `cfd_lbm`, `dem`,
  `particulate_flow`).
- **CONTRIBUTING.md**: Fixed build-backend reference (hatchling → poetry-core),
  updated project structure to reflect all 4 packages, corrected coverage command.

### Fixed
- **CylinderFlow SOR sign bug**: Poisson solver was updating pressure in the
  wrong direction.
- **CylinderFlow Poisson RHS sign**: RHS was −ρ/dt·∇·u (wrong), now +ρ/dt·∇·u*.
- **docs/VIBES.md**: Removed accidentally committed credential; replaced with
  environment-variable guidance.

---

## [0.2.0] - 2026-02-20

### Added
- **Non-Newtonian Solver**: New `PowerLawPlanePoiseuille` class for
  shear-thinning and shear-thickening fluids.
- **Animation Workflow**: `make_animation_frames.py` script for producing
  frame sequences suitable for ffmpeg-based video generation.
- **Stabilized Navier-Stokes**: Re-implemented `CavityFlow` with upwind
  differencing, dynamic time-stepping, and velocity clipping for improved
  numerical stability.
- **Project-wide logging**: Standardized `logging` usage across all modules
  with structured `%(levelname)s: %(message)s` format.

### Changed
- **SOR Upgrade**: Upgraded Poiseuille solvers from Gauss-Seidel to SOR
  (Successive Over-Relaxation) for faster convergence.
- **Documentation**: Updated `README.md`, `docs/TUTORIAL.md` to reflect new
  capabilities.

---

## [0.1.0] - 2024-02-17

### Added
- **Initial release**
- `PlanePoiseuille`: Analytical and numerical solution for 2D channel flow.
- `CircularPoiseuille`: Hagen-Poiseuille flow in a circular pipe.
- `CavityFlow`: Lid-driven cavity flow solver.
- Comprehensive test suite covering analytical solutions and boundary conditions.
- matplotlib-based visualization for all solvers.
