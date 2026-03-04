# CHANGELOG

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
- **Coupled LBM-DEM solver** (`cfd_dem_lbm` package): Bilinear interpolation
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
  `cfd_dem_lbm`).
- **CONTRIBUTING.md**: Fixed build-backend reference (hatchling → poetry-core),
  updated project structure to reflect all 4 packages, corrected coverage command.

### Fixed
- **CylinderFlow SOR sign bug**: Poisson solver was updating pressure in the
  wrong direction.
- **CylinderFlow Poisson RHS sign**: RHS was −ρ/dt·∇·u (wrong), now +ρ/dt·∇·u*.
- **VIBES.md**: Removed accidentally committed credential; replaced with
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
- **Documentation**: Updated `README.md`, `TUTORIAL.md` to reflect new
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
