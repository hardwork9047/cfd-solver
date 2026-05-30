# CHANGELOG

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [0.17.0] - 2026-05-29

### Changed

- Extracted the 2D DEM contact-force laws into a fluid-independent `ContactLaw` dataclass (`particulate_flow/dem/contact_law.py`), see PR #30: the clamp-to-zero Hertz/linear normal force, Coulomb-limited tangential force, and Coulomb-limited rolling torque now live in a small class holding only the 8 material parameters, so the contact physics can be built and tested without constructing a fluid solver. `DEMSolver` delegates its three contact-law methods to a `ContactLaw` built from the coupled solver's params; the formulas are unchanged, so `FastLBMDEM` production results are bit-identical. This is the first step of the larger God-class split (issue #26).

### Notes

- Pure structural refactor — no physics or numerical change. `ContactLaw` is the 2D counterpart of the scalar contact-law helpers on the 3D `DEM3D` and shares the same clamp/Coulomb concepts. Decoupling the rest of `DEMSolver` (pair search, drag, IBM) from the fluid solver remains future work.

---

## [0.16.0] - 2026-05-29

### Added

- `/simulate` skill 3D support (issue #16), see PR #29: the skill now recognises 3D intent ("3D" / `nz`) and generates a config that extends `fouling_supply_3d.json` with `domain.dimensions=3` + `nz`, uses the 3D cylinder geometry (`four_cylinder_3d.json`), shows `dimensions`/`nz` in the confirmation table, and warns when a 2D-only template/geometry is requested under 3D intent. 2D default behaviour is unchanged. `.claude/skills/simulate/SKILL.md` only (no Python change).

### Notes

- Completes the 3D membrane-fouling effort: the full stack (#14/#17/#18/#19/#20 solver, #15 config/runner, #16 skill) now lets a natural-language request drive a 3D simulation end-to-end, mirroring the 2D operability.

---

## [0.15.0] - 2026-05-29

### Added

- 3D membrane-fouling config + runner plumbing (issue #15), see PR #28: `run_lbm_dem.py --config <3d_case>` now dispatches to a dedicated 3D path (`particulate_flow/runner3d.py`) when `domain.dimensions=3`, building `LBMDEMSolver3D` and writing `analysis/time_series.csv`+`.npz`, `summary.json`, `run_status.json`, and ParaView 3D VTK (structured fluid field + particle point cloud + `.pvd` time series). The 2D pipeline is untouched.
- `configs/lbm_dem/templates/fouling_supply_3d.json` (3D fouling template), `configs/lbm_dem/geometries/four_cylinder_3d.json` (z-aligned cylinder fragment, addable via `extends`), and `configs/lbm_dem/cases/fouling_3d_smoke.json` (tiny end-to-end case).

### Notes

- Makes the merged 3D solver stack (#14/#17/#18/#19/#20) runnable from a single command. Remaining: `/simulate` skill 3D support (#16). The 3D path writes ParaView VTK rather than the 2D matplotlib `.mp4`; 3D matplotlib animation is out of scope.

---

## [0.14.0] - 2026-05-29

### Added

- 3D left-inlet particle injection/removal for `LBMDEMSolver3D` (`src/particulate_flow/lbm3d.py`), see PR #25: `particle_source="left_inlet"` (pressure mode) injects particles at the inlet plane over time per a volume-flux budget (`+= source_volume_fraction · inlet_flux`, sphere cost `(4/3)π r³`) at flux-weighted, non-overlapping (y, z) points, and removes particles once they pass the outlet. The injected solid-volume fraction tracks `source_volume_fraction`.
- `DEM3D.add_particles` / `DEM3D.remove_particles` for a dynamic particle population (grow/shrink all state arrays consistently).

### Notes

- Final slice of the 3D membrane-fouling effort (issue #20); the 3D stack (#14 pressure flow, #18 DEM contact, #17 IBM coupling, #19 cylinder obstacles, #20 inlet injection) is now complete. Retires the `particle_source` `NotImplementedError` stub from #14. Per-particle radius distributions and Numba acceleration of the 3D kernels remain future work.

---

## [0.13.0] - 2026-05-29

### Added

- 3D fixed finite z-aligned cylinder obstacles for `LBMDEMSolver3D` (`src/particulate_flow/lbm3d.py`), see PR #24: `cylinders` build a solid mask (`(x-cx)² + (y-cy)² ≤ r²` within the z-extent) with halfway bounce-back for no-slip, mirroring the 2D cylinder support. IBM force spreading skips solid nodes, and the same cylinders are forwarded to the internal `DEM3D` so particles collide with the obstacles. Accepts `(cx, cy, r)` or `(cx, cy, r, z_lo, z_hi)`.

### Notes

- Fourth slice of the 3D membrane-fouling effort (issue #19). Retires the `cylinders` `NotImplementedError` stub from #14. Remaining: inlet particle injection (#20). Numba acceleration and arbitrary cylinder orientation are out of scope.

---

## [0.12.0] - 2026-05-29

### Added

- 3D immersed-boundary particle-fluid coupling for `LBMDEMSolver3D` (`src/particulate_flow/lbm3d.py`), see PR #23: D3Q15 Guo body forcing, spherical Fibonacci-lattice surface markers, trilinear velocity interpolation and force spreading, and a coupled `advance()` that drives an internal `DEM3D` (issue #18) via the IBM reaction force. Enabled with `n_particles > 0` + `particle_fluid_coupling="immersed_boundary"` + `particle_positions`.
- `DEM3D.step()` gains `external_forces` / `external_torques` arguments so a coupling layer can inject per-particle fluid loads while keeping `DEM3D` fluid-agnostic.

### Notes

- Third slice of the 3D membrane-fouling effort (issue #17). Completes two-way fluid↔particle coupling. Remaining: 3D fixed-cylinder fluid masks (#19) and inlet particle injection (#20). Numba acceleration of the 3D IBM kernels is deferred.

---

## [0.11.0] - 2026-05-29

### Added

- Fresh 3D DEM contact physics `DEM3D` (`src/particulate_flow/dem/contact3d.py`): sphere-sphere, sphere-wall, and sphere-(z-aligned finite cylinder) Hertz/linear normal contact with 3-vector angular velocity and torque, Coulomb-limited tangential friction along the in-plane slip direction, Coulomb-limited rolling resistance, buoyancy-corrected gravity, and semi-implicit Euler sub-stepped integration, see PR #22. The 2D `DEMSolver` is unchanged.

### Notes

- Second slice of the 3D membrane-fouling effort (issue #18). `DEM3D` is fluid-free; wiring it into `LBMDEMSolver3D.advance()` with drag + IBM back-reaction is issue #17. Numba acceleration of the 3D kernels is deferred.

---

## [0.10.0] - 2026-05-29

### Added

- 3D pressure-driven flow for `LBMDEMSolver3D`: selectable `streamwise_boundary="pressure"` applies a D3Q15 Zou-He density inlet/outlet on x (y, z periodic), mirroring the 2D contract `rho_in = rho_out + pressure_drop / cs²`, see PR #14. Default stays `"periodic"` (Lees-Edwards path unchanged).
- `build_lbm_dem_solver` forwards `streamwise_boundary`, `pressure_drop`, `rho_out` to the 3D solver.

### Notes

- First vertical slice of the 3D membrane-fouling effort (issue #14). The 3D particle stack — IBM coupling (#17), DEM contact (#18), fixed obstacles (#19), inlet injection (#20) — is not yet implemented; `LBMDEMSolver3D` raises `NotImplementedError` if particles, cylinders, or an inlet source are requested.

---

## [0.9.0] - 2026-05-28

### Added

- 3D LBM solver `LBMDEMSolver3D` using D3Q15 lattice (`src/particulate_flow/lbm3d.py`), with Lees-Edwards shear BC (x/z periodic, y LE), see PR #9.
- `build_lbm_dem_solver` dispatches to 3D solver when `dimensions=3` in config/CLI.
- New CLI args `--dimensions` (2 or 3) and `--nz` for 3D grid depth.
- Sample config `configs/lbm_dem/cases/lees_edwards_shear_flow_3d.json`.

---

## [0.8.0] - 2026-05-28

### Added

- DEM surface roughness (`physics.surface_roughness` / `--surface-roughness`): extends particle-particle contact threshold by h_r lattice units while keeping geometry-based surface-force gaps unchanged, see PR #8.
- Apparent viscosity evaluation (`runtime.viscosity_eval.enabled`): `ViscosityEvaluator` in `particulate_flow/rheology.py` computes η^H (hydrodynamic) and η^C (collisional) contributions per interval, writes `viscosity_timeseries.csv`, and records time-averaged η_s in `summary.json`, see PR #8.

---

## [0.7.0] - 2026-05-28

### Added

- Lees-Edwards (LE) boundary conditions (`y_boundary="lees_edwards"`) for wall-less shear flow, see PR #7.
- `isp` particle-fluid coupling mode (interpolated Stokes point-force), complementing existing `point_force` and `immersed_boundary`.
- `le_shear_rate`, `le_shear_axis`, `le_boundary_axis`, `le_interpolation_order` parameters on `LBMDEMSolver`.
- `solver.lees_edwards` JSON config sub-section; `enabled: true` automatically sets `y_boundary` and `streamwise_boundary`.
- Example config `configs/lbm_dem/cases/lees_edwards_shear_flow.json`.

### Changed

- `Y_BOUNDARIES` extended with `"lees_edwards"`; `PARTICLE_FLUID_COUPLINGS` extended with `"isp"`.
- `_flatten_sections` now expands `solver.lees_edwards` sub-dict into flat `le_*` argparse keys.

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
