# LBM-DEM Fouling Workflow

This note records the current practical workflow for developing the membrane
fouling model.

## 1. Solver method options

`run_lbm_dem.py` exposes three method axes that can be included in JSON
sweeps.

### Fluid solver method

Use `--fluid-method`.

- `lbm-bgk-guo`: single-relaxation-time BGK LBM with Guo forcing.  This is the
  default and the simplest reference method.  It is useful for baseline runs,
  quick comparisons, and continuity with earlier results.  Its weakness is that
  the single relaxation parameter controls all non-conserved modes, so it can be
  less robust near solid boundaries, narrow gaps, and high blockage.
- `lbm-trt-guo`: two-relaxation-time LBM with Guo forcing.  The even and odd
  distribution components relax separately, which usually gives better boundary
  behavior than BGK for bounce-back solids.  This is the preferred candidate for
  pore-blockage studies, especially with `solid_boundary` coupling.  It is a
  little more expensive and should be compared against BGK during validation.

### Particle contact method

Use `--particle-method`.

- `dem-hertz`: nonlinear Hertz normal contact with damping.  This is the
  default and is closer to elastic sphere/disc contact.  Contact stiffness grows
  with overlap, so it is suitable when the qualitative packing and collision
  response matter.  It can require smaller DEM substeps or careful stiffness
  tuning in dense fouling layers.
- `dem-linear`: linear spring normal contact with damping.  This is simpler and
  often easier to reason about when debugging or performing sensitivity sweeps.
  It is less physically specific than Hertz contact, but can be useful for
  isolating whether a result comes from geometry, adhesion/repulsion, or the
  nonlinear contact law.

### Particle-fluid coupling method

`run_lbm_dem.py` supports two coupling modes.

- `point_force`: the previous two-way Stokes-drag back reaction.  Particles
  exchange momentum with the fluid but do not geometrically block lattice nodes.
  This is computationally light and useful for legacy comparison, dilute
  transport, and sensitivity studies where complete pore blockage is not the
  central question.  It should not be interpreted as a fully blocked membrane
  pore even if particles visually fill a gap.
- `solid_boundary`: active DEM particles are overlaid on the LBM `solid` mask.
  This makes deposited particles act as moving no-slip obstacles for the fluid
  solver.  It is still a 2-D approximation, but it is the preferred mode when
  testing pore blockage and flux decline.  It is more physical for fouling
  resistance, but also more sensitive to lattice resolution, particle radius,
  and bounce-back boundary accuracy.
- `immersed_boundary`: circular DEM particle surfaces are represented by
  Lagrangian marker points.  The solver interpolates fluid velocity to each
  marker, compares it with the local particle surface velocity, spreads a
  direct-forcing correction back to the lattice, and applies the opposite force
  and torque to the particle.  This avoids stair-step solid masks and represents
  moving particle surfaces more naturally than `solid_boundary`.  It is still a
  penalty/direct-forcing IBM, so `--ibm-stiffness`, `--ibm-marker-spacing`, and
  grid resolution must be checked with verification cases.

### Recommended combinations

- Baseline/legacy comparison:
  `--fluid-method lbm-bgk-guo --particle-method dem-hertz --particle-fluid-coupling point_force`
- Pore-blockage screening:
  `--fluid-method lbm-trt-guo --particle-method dem-hertz --particle-fluid-coupling solid_boundary`
- Moving-particle hydrodynamic coupling:
  `--fluid-method lbm-trt-guo --particle-method dem-hertz --particle-fluid-coupling immersed_boundary`
- Debug/sensitivity runs:
  `--fluid-method lbm-bgk-guo --particle-method dem-linear --particle-fluid-coupling solid_boundary`

For publication-quality interpretation, compare at least BGK vs TRT and point
force vs solid-boundary coupling on a small set of representative geometries.

Example:

```bash
python src/demos/run_lbm_dem.py \
  --particle-fluid-coupling solid_boundary \
  --fluid-method lbm-trt-guo \
  --particle-method dem-hertz
```

### Acceleration and neighbour-search options

Use `--fluid-accelerator` to choose the LBM execution backend.

- `numpy`: stable default path.  Use this for reference runs and when comparing
  against older results.
- `numba`: compiled LBM collision/streaming/bounce-back path when the optional
  `numba` package is installed.  If numba is not available, the solver records
  that the NumPy path was used.
- `auto`: use numba when available, otherwise NumPy.

Use `--particle-search` to choose the DEM pair search.

- `cell_list`: production default.  Particles are binned into cells sized by
  the maximum contact/surface-force interaction distance, so only neighbouring
  cells are checked.
- `all_pairs`: exhaustive `O(N^2)` search for debugging small cases.  This is
  useful when verifying that cell-list results match a brute-force reference.

For large fouling sweeps, prefer:

```bash
python src/demos/run_lbm_dem.py \
  --fluid-method lbm-trt-guo \
  --fluid-accelerator auto \
  --compute-accelerator auto \
  --particle-search cell_list \
  --output-profile analysis
```

To measure the backend speed on the local machine, run:

```bash
poetry run python src/bin/benchmark_lbm_accelerators.py
```

On the 180 x 70 fixed-cylinder benchmark run on 2026-05-09, the compiled Numba
backend reached 1093 median steps/s versus 315 median steps/s for the NumPy
path after compilation warmup, a 3.47x speedup.  The first Numba call includes
JIT compilation overhead, so short runs should be interpreted separately from
long production runs.

Use `--compute-accelerator` to choose the non-LBM numerical backend for DEM
boundary loads, particle-solid-mask updates, and immersed-boundary marker
bookkeeping.  To measure these components independently:

```bash
poetry run python src/bin/benchmark_numba_components.py
```

On the 300-particle, 180 x 70 benchmark run on 2026-05-09, the Numba backend
gave these median component speedups:

- particle pair loads: 14.93x
- DEM wall/cylinder boundary loads: 1.82x
- solid-boundary particle mask update: 240.09x
- immersed-boundary marker bookkeeping: 3.18x

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
- hydrodynamic load generation when a particle is represented by the immersed
  boundary method.

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
src/bin/run_lbm_dem_config_sweep.py src/bin/lbm_dem_factor_sweep_example.json --jobs 4
src/bin/run_lbm_dem_config_sweep.py src/bin/lbm_dem_factor_sweep_example.json --dry-run
```

The runner writes logs and `status.tsv` under
`src/results/run_lbm_dem/config_sweeps/<timestamp>/`.

Two config styles are supported:

- `cases`: explicitly list named cases.
- `factors`: define lists of values and let the runner expand the Cartesian
  product into cases.
- `--jobs N`: run up to `N` cases in parallel.  Use this for independent
  sweeps, choosing `N` so memory and CPU usage remain comfortable.

After a sweep, summarize all available run directories with:

```bash
src/bin/summarize_lbm_dem_results.py --root src/results/run_lbm_dem
```

## 5. Flow Conditions

Use `--flow-condition` to make the intended boundary/control condition explicit.

- `fixed-pressure`: keeps the body force fixed.  This is closest to a fixed
  pressure-gradient calculation.  If fouling blocks the channel, the permeate
  flux can decline.
- `target-max-velocity`: adjusts the body force to keep the maximum fluid speed
  near `--u-max`.  This is useful for fair geometry comparisons at a controlled
  velocity scale, but the changing drive force must be interpreted as an
  imposed control action.

The older `--flow-control` / `--no-flow-control` flags are still accepted for
backward compatibility, but `--flow-condition` is preferred for new sweeps.

## 6. Performance Notes

The current execution path includes several speed-oriented choices:

- `solid_boundary` updates only local particle bounding boxes instead of
  rebuilding particle solids from a full-grid mesh each step.
- DEM particle-particle interactions use a cell-list neighbour search by
  default; `--particle-search all_pairs` remains available only for small
  debugging/reference cases.
- Stokes drag and particle-fluid interpolation are batched across particles.
- `immersed_boundary` computes marker interpolation and marker forces in
  batched arrays before spreading them back to the lattice.
- LBM streaming uses precomputed periodic source indices instead of repeated
  `np.roll` calls.
- `--fluid-accelerator numba` or `--fluid-accelerator auto` can use a compiled
  LBM step when numba is installed.
- `--output-profile analysis` disables movie generation and writes only
  analysis outputs plus final ParaView data without storing intermediate
  snapshot NPZ files; `--output-profile minimal` skips video, ParaView, and
  intermediate snapshots for fast parameter screening.
- Config sweeps can run independent cases concurrently with `--jobs`.

For large sweeps, reduce I/O first: increase `--snapshot-every`, set
`--paraview-every 0` for final-only VTK, and use `--no-video` unless movies are
needed for every condition.

## 7. Solver API Direction

The production solver now accepts explicit numerical-method controls:

- `fluid_method`
- `fluid_accelerator`
- `particle_method`
- `particle_search`
- `particle_fluid_coupling`
- `spatial_dimension`

Only `spatial_dimension=2` is implemented today, but keeping the dimension as
an explicit constructor argument avoids hiding 2-D assumptions in scripts.  A
future 3-D implementation should preserve this API and swap in D3Q19/D3Q27 LBM,
3-D DEM contact geometry, and 3-D particle-fluid coupling behind the same
method-selection surface.
