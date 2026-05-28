# Issue-7: Phase A — Lees-Edwards BCs and iSP fluid-particle coupling (2D)

## Goal (verbatim from issue)

Lees-Edwards (LE) 境界条件と iSP 流体-粒子結合を2Dソルバーに追加する。
完了後は、JSON config に `solver.lees_edwards` セクションを記述するだけで
壁なしのずり流れ懸濁液シミュレーションを実行できる。
既存の `immersed_boundary` / `point_force` モードは引き続き動作する。

## Approach

### 1. Constants & enumerations

Add `"lees_edwards"` to `Y_BOUNDARIES` in `src/particulate_flow/lbm/constants.py`.
Add `"isp"` to `PARTICLE_FLUID_COUPLINGS`.

### 2. Lees-Edwards streaming (LBM side)

Lees-Edwards BCs impose a shear rate γ̇ between the top and bottom periodic images.
During the streaming step, populations that cross the y boundary at y=0 or y=ny-1
are shifted in x by ±Δx_LE = γ̇ · ny · t (accumulated fractional shift tracked in
`self._le_shift`). A sub-lattice shift is implemented via first-order interpolation
(or third-order with `interpolation_order=3`).

Implementation:
- New attributes on `LBMDEMSolver.__init__`: `le_shear_rate`, `le_shear_axis`,
  `le_boundary_axis`, `le_interpolation_order`, `_le_shift` (float, accumulated).
- After the normal streaming step, apply the LE correction to populations that
  crossed the y boundary (directions 2,4,7,8 in D2Q9 numbering cross y=0 or y=ny-1).
- The correction shifts these rows in x by the current accumulated displacement,
  using `np.roll` for integer part + linear/cubic interpolation for fractional part.

### 3. iSP fluid-particle coupling

"iSP" (interpolated Stokes point-force) coupling: like `point_force` but the
interpolated velocity at the particle centre uses the smoothed (background) field
excluding the particle's own contribution. For Phase A this is identical to
`point_force` except:
- Uses third-order (cubic) velocity interpolation (controlled by `le_interpolation_order`).
- Adds the local background velocity correction to the Stokes drag.

This is the minimal iSP needed for correctness at φ~0.01; the full regularised
iSP is Phase B scope.

### 4. DEM particle wrapping

Particles that cross the LE y boundary must have their x-position incremented by
±`le_shift` (the same accumulated shift) to stay in the co-moving frame. Implemented
in `_dem_substep` when `y_boundary == "lees_edwards"`.

### 5. Config system

`SECTION_KEYS` already includes `solver`. No new section key needed.
`_flatten_sections` will pass `lees_edwards` sub-keys (enabled, shear_rate,
shear_axis, boundary_axis, interpolation_order) through as flat keys prefixed
by `le_*` after remapping in `_flatten_sections`.

`SimulationConfig.from_json` flattens all solver keys; the runner passes them to
`build_lbm_dem_solver` via argparse defaults, which reads them with `getattr`.

### 6. Builder & runner

`build_lbm_dem_solver` picks up the new `le_*` args and sets:
- `y_boundary="lees_edwards"`
- `streamwise_boundary="periodic_force"` (wall-less shear uses body-force drive)
- `le_shear_rate`, `le_shear_axis`, `le_boundary_axis`, `le_interpolation_order`

### 7. Test config

Add `configs/lbm_dem/cases/lees_edwards_shear_flow.json` as a minimal example.

## Alternatives considered

- **Full spectral interpolation** for LE streaming: overkill for Phase A; linear/cubic
  is sufficient for convergence check at < 1% L2 error.
- **GPU-accelerated LE kernel**: deferred to Phase C.
- **Separate LE solver class**: rejected in favour of conditional branching inside
  existing `LBMDEMSolver` to avoid duplicating ~400 lines of DEM logic.

## Out of scope (Phase A)

- Surface roughness / viscosity evaluation (Phase B)
- 3D extension (Phase C)
- Full regularised iSP with Rotne-Prager corrections (Phase B)

## Assumptions

- `D2Q9` direction numbering matches `src/particulate_flow/lbm/constants.py`.
- Streaming precomputed indices in `_stream_x_src` / `_stream_y_src` do not need
  to change; LE correction is applied as a post-stream fix to the crossing rows.

## Acceptance scenario test mapping

| Scenario | Test |
|----------|------|
| Linear velocity profile u_x = γ̇·y, L2 < 1% | `test_le_linear_shear_profile` |
| shear_axis / boundary_axis configurability | `test_le_axis_configuration` |
| interpolation_order 1 and 3 | `test_le_interpolation_order_toggle` |
| iSP coupling mode selectable | `test_isp_coupling_selectable` |
| particle wrap across LE boundary | `test_le_particle_wrapping` |
| existing fouling/cylinder regression | existing `test_integration.py` |
