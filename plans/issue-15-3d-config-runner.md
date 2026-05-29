# Plan: Issue #15 — 3D Membrane-Fouling Config & Runner Plumbing

## Verbatim Goal (from issue #15)

> 3D膜ファウリングシミュレーションをコマンド1行で起動できるよう、テンプレート
> JSON・ジオメトリJSON・ランナーの後処理を整備する。

### Acceptance Scenarios
1. `fouling_supply_3d.json` テンプレートを extends して最小限の差分で3Dケースを定義できる
2. `run_lbm_dem.py --config <3d_case>.json` が正常終了し結果ディレクトリを生成する
3. 結果ディレクトリに time_series.csv・summary.json・ParaView用VTKファイルが出力される
4. 3Dシリンダー配置ジオメトリJSONを extends に追加するだけで障害物が有効になる

### Constraints
- Package: particulate_flow / configs
- #14 (3D core) — merged; the full 3D stack (#17/#18/#19/#20) is also merged
- 2D ケースの JSON 設定・出力フォーマットと後方互換を保つ
- `run_lbm_dem.py` の既存 CLI フラグを変更しない

## Scope decision

The runner (`run_lbm_dem.py`, 1,761 lines) is a monolithic **2D** fouling pipeline:
the main loop, fouling-specific analysis, matplotlib animation, and VTK writers all
assume the 2D `FastLBMDEM` (`get_fields()` → 2D, particle circles, 2D structured VTK).
`--dimensions 3` exists and `build_lbm_dem_solver` dispatches to `LBMDEMSolver3D`, but
**no 3D output path exists** — the existing `lees_edwards_shear_flow_3d.json` was a
solver-level test (#9), never run through this runner end-to-end.

Threading 3D through every 2D helper would be invasive and risky.  Instead:

### In scope (this PR)
- **A dedicated 3D execution branch** in `run_lbm_dem.py`: when `args.dimensions == 3`,
  run a streamlined loop (build `LBMDEMSolver3D` via the existing builder → advance in
  snapshot chunks → collect dimension-agnostic scalar rows → write
  `analysis/time_series.csv`+`.npz`, `summary.json`, ParaView 3D VTK, `run_status.json`,
  `metadata.json`) and **return before** the 2D matplotlib/fouling machinery.  The 2D
  path is left byte-for-byte unchanged.
- **3D VTK writers**: fluid as a `STRUCTURED_POINTS` scalar/vector field over
  (nx,ny,nz); particles as `POLYDATA` points with radius/velocity; a `.pvd` time series
  (reuse the existing `_write_pvd`, `_vtk_float`).
- **`configs/lbm_dem/templates/fouling_supply_3d.json`**: the 3D analogue of
  `fouling_supply.json` (`domain.dimensions=3`, `nz`, pressure flow, `left_inlet`,
  `immersed_boundary`).
- **`configs/lbm_dem/geometries/` 3D cylinder fragment**: a z-aligned finite-cylinder
  array that a case can `extends` to enable obstacles.
- **A runnable 3D case** under `configs/lbm_dem/cases/` (small grid) for the smoke test.
- Tests: a fast smoke test that the 3D runner path produces the required artifacts;
  config-load tests for the new templates; a 2D-regression guard.

### Out of scope (other issues / future)
- `/simulate` skill 3D support → #16.
- 3D matplotlib animation / mp4 (ParaView VTK is the 3D visualization path; the 2D
  `.mp4` animation is not ported).
- Refactoring the 2D runner into shared dimension-agnostic helpers (that is the #26
  God-class-split track); here the 3D branch only *reuses* already-agnostic helpers.
- Numba for the 3D solver.

## Approach

### Runner 3D branch
Add `run_3d(args, sim)` invoked right after the solver is built when
`dimensions == 3`.  It:
1. warms up (optional) then advances in `snapshot_every` chunks to `total_steps`,
2. each chunk records a scalar row: step, `|u|_max`, mean ρ, particle count, mean/max
   particle speed, injected/removed counts (from `LBMDEMSolver3D` / its `dem`),
3. writes the standard artifacts via the **existing dimension-agnostic** helpers
   (`_write_analysis_outputs`, `_write_run_status`, `_write_metadata`) plus new 3D VTK,
4. prints a concise completion summary and `return`s (no matplotlib).

### 3D VTK
- Fluid: `STRUCTURED_POINTS` (DIMENSIONS nx ny nz), POINT_DATA with `rho` scalar and
  `velocity` vector; mark solid cells via a `solid` scalar.
- Particles: `POLYDATA` points with `radius` and `velocity` arrays.
- `.pvd` index over exported frames (reuse `_write_pvd`).

### Config templates
- `fouling_supply_3d.json`: small-ish 3D domain, `dimensions:3`, `nz`, pressure BC,
  `left_inlet`, `immersed_boundary`, sphere radius, source fraction.
- `geometries/four_cylinder_3d.json` (or similar): z-aligned cylinder specs valid for
  the template's domain, addable via `extends`.

## Alternatives considered & rejected
- **Thread 3D through every 2D helper**: rejected — invasive, high regression risk on
  the production 2D path; the issue only needs the 3D *artifacts*, not the 2D
  matplotlib pipeline in 3D.
- **A separate `run_lbm_dem_3d.py` script**: rejected — the issue says
  `run_lbm_dem.py --config <3d>` must work (scenario 2), so the dispatch must live in
  the existing entrypoint.
- **Port the mp4 animation to 3D**: rejected — 3D visualization is ParaView's job;
  a 3D matplotlib animation is low-value and heavy.

## Assumptions about other modules (docstring-grounded)
- `build_lbm_dem_solver` "Dispatches to the 3-D solver when `args.dimensions == 3`"
  (builder docstring) — the 3D branch builds the solver via the same call.
- `LBMDEMSolver3D` exposes `advance`, `get_fields` (→ rho,ux,uy,uz), `solid`, `dem`
  (with `pos/vel/radii/n_p`), and inlet counters (`injected_particle_volume`,
  `removed_particles`) per its docstring — the 3D rows/VTK read these.
- `SimulationConfig` "loads and flattens sections … binds them to CLI argument
  parsers" with `extends` inheritance (CLAUDE.md) — the new templates use `extends`.
- `_write_analysis_outputs`/`_write_pvd`/`_vtk_float` are dimension-agnostic (write
  dict rows / generic VTK scaffolding).

## Testing plan (acceptance → test) — `tests/test_runner_3d.py` + config tests
- **Scenario 1** (template extends): a minimal 3D case that `extends`
  `fouling_supply_3d.json` loads via `SimulationConfig` and yields `dimensions==3`,
  `nz>1`, pressure flow.
- **Scenario 2/3** (runner end-to-end + artifacts): invoke the 3D runner branch on a
  tiny grid for a few steps (subprocess or direct `run_3d` call) and assert the result
  dir contains `analysis/time_series.csv`, `summary.json`, and a `paraview/` dir with a
  fluid `.vtk` and `.pvd`; assert it exits 0 / returns without error.
- **Scenario 4** (cylinder geometry via extends): a case extending the 3D template +
  the 3D cylinder geometry builds a solver whose `solid.any()` is True.
- **Regression**: a 2D config still runs the 2D path unchanged (existing
  `test_runner_refactor.py` stays green); the full suite passes.
