# Plan: Issue #14 — 3D LBM-DEM Core Solver Extension

## Verbatim Goal (from issue #14)

> `LBMDEMSolver3D`（lbm3d.py）に圧力駆動流・IBM粒子結合・DEM粒子運動・固定障害物・
> 左inlet粒子注入を追加し、2Dの `FastLBMDEM` と同等の膜ファウリング計算を
> 3D（D3Q15）で実行できるようにする。

### Acceptance Scenarios (issue text)
1. `dimensions: 3` を指定したコンフィグで圧力駆動流（x方向inlet/outlet）が動作する
2. 3D空間で粒子が流体drag・重力・DEM接触力を受けて運動する
3. 3D固定シリンダー障害物（有限長円柱）が固体境界として機能する
4. 左inletから粒子が注入され、右outletで削除される
5. `build_lbm_dem_solver(args)` が `dimensions=3` のとき新しい3Dソルバーを返す
6. 既存の2Dソルバー（`FastLBMDEM`）の挙動が変わらない

## Scope Decision (agreed with user: "vertical slice end-to-end")

Issue #14 bundles **five** large subsystems that would together rival the 9,500-line
2D solver. The existing 2D `DEMSolver` is **not** reusable for 3D (it is hardwired to
2D contact: `compute_contact_stress_xy`, `tangent_from_normal`, and reads 2D arrays).
Delivering all five in one PR would produce an unreviewable diff with high physics-bug
risk.

**This PR delivers a complete, reviewable vertical slice:**

### In scope (this PR)
- **3D pressure-driven flow** — Zou-He pressure inlet/outlet on the x-boundaries of
  `LBMDEMSolver3D`, mirroring the 2D contract `rho_in = rho_out + pressure_drop / cs²`.
  Selectable via a `streamwise_boundary="pressure"` mode; default stays `"periodic"`
  so Lees-Edwards behaviour is untouched.
- **Dimension dispatch parity** — `build_lbm_dem_solver(args)` already dispatches to 3D
  when `dimensions==3`; extend `_build_lbm3d_solver` to forward `streamwise_boundary`,
  `pressure_drop`, `rho_out` so a pressure-driven 3D config actually runs.
- **Honest stubs for the particle stack** — IBM coupling, DEM contact, fixed obstacles,
  and inlet injection are exposed as constructor parameters that raise
  `NotImplementedError` with a message pointing at the follow-up issue. This makes the
  remaining surface explicit in code rather than silently absent.
- **Tests** for acceptance scenarios 1, 5, 6 (the fluid-side + regression scenarios),
  plus tests asserting the particle-stack stubs raise informatively (guarding scenarios
  2–4 against accidental silent no-ops).

### Out of scope (follow-up issues, split from #14)
- **#17** — 3D IBM particle–fluid coupling (scenario 2, fluid↔particle force exchange)
- **#18** — 3D DEM contact physics: sphere–sphere, sphere–wall, sphere–cylinder (scenarios 2–3)
- **#19** — 3D fixed finite-cylinder solid masks (scenario 3)
- **#20** — 3D left-inlet particle injection/removal (scenario 4)
- GPU / TRT-3D (already out of scope per issue Constraints)

These four follow-ups are filed and linked from the PR. #14 stays open as the tracking
epic; the PR will note the slice + children rather than auto-closing #14.

## Approach

### 1. Pressure BC in `LBMDEMSolver3D`
- Add constructor params: `streamwise_boundary: str = "periodic"`, `pressure_drop: float = 0.0`,
  `rho_out: float = 1.0`. Compute `rho_in = rho_out + pressure_drop / CS2_3`.
- Streaming stays periodic (np.roll); after streaming, when
  `streamwise_boundary == "pressure"`, apply Zou-He density BC on x=0 (inlet, ρ=rho_in)
  and x=nx-1 (outlet, ρ=rho_out). Implement a D3Q15 Zou-He that fixes density and the
  transverse momentum to zero at the boundary planes, solving for the unknown incoming
  populations. Keep y/z periodic and skip the LE correction when in pressure mode
  (mutually exclusive: pressure flow has le_shear_rate == 0).
- Validation hook: a pressure-driven duct with periodic y,z should reach a uniform
  bulk velocity in x with ρ decreasing linearly from inlet to outlet.

### 2. Builder forwarding
- `_build_lbm3d_solver`: forward `streamwise_boundary`, `pressure_drop`, `rho_out`
  from args (with getattr defaults matching current behaviour). No change to the 2D path.

### 3. Stubs
- Constructor accepts (and stores) `particle_fluid_coupling`, `n_particles`,
  `cylinders`, `particle_source`. If any requests the not-yet-built particle stack
  (e.g. `n_particles > 0`, non-empty `cylinders`, `particle_source != "none"`),
  raise `NotImplementedError` pointing at the relevant follow-up (#17 IBM, #18 DEM,
  #19 obstacles, #20 inlet).

## Alternatives considered & rejected
- **All-in-one PR**: rejected — unreviewable, high risk (user agreed).
- **Reuse 2D DEMSolver in 3D**: rejected — it is 2D-hardwired (verified via docstrings:
  `compute_contact_stress_xy`, `tangent_from_normal`). 3D needs fresh contact physics.
- **Wall (bounce-back) instead of Zou-He for x**: rejected — issue explicitly says
  pressure-driven inlet/outlet; Zou-He matches the 2D contract.

## Assumptions about other modules (docstring-grounded)
- 2D `LBMDEMSolver` sets `rho_in = rho_out + pressure_drop / CS2` and applies "Zou-He
  pressure inlet/outlet on the left/right boundaries" (lbm_dem.py:269, :894 docstrings).
  3D mirrors this with `CS2_3 = 1/3`.
- `build_lbm_dem_solver` "Dispatches to the 3-D solver when `args.dimensions == 3`"
  (builder.py docstring) — dispatch exists; only forwarding is missing.
- `LBMDEMSolver3D.advance()` / `get_fields()` signatures must stay unchanged (issue
  Constraints) — pressure BC is added inside `advance`'s per-step loop, no signature change.

## Testing plan (acceptance → test)
- **Scenario 1** (3D pressure flow): `tests/test_lbm3d_pressure.py` — build a 3D solver
  with `streamwise_boundary="pressure"`, `pressure_drop>0`, run to quasi-steady, assert
  (a) positive bulk ux, (b) ρ monotonically decreasing inlet→outlet, (c) finite/stable.
- **Scenario 5** (dispatch): assert `build_lbm_dem_solver` with `dimensions=3,
  streamwise_boundary="pressure"` returns an `LBMDEMSolver3D` whose `rho_in>rho_out`.
- **Scenario 6** (2D regression): existing suite must stay green; add an explicit test
  that `dimensions=2` still returns `FastLBMDEM`.
- **Scenarios 2–4** (stubs): assert that requesting particles/cylinders/inlet on the 3D
  solver raises `NotImplementedError` mentioning the follow-up issue — so the gap is
  loud, not silent.
- Pre-existing LE tests (`tests/test_lbm3d.py`) must remain green (default periodic mode).
