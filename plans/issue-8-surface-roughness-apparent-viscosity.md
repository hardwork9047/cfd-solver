# Plan: Issue 8 — DEM Surface Roughness and Apparent Viscosity Evaluation

## Goal (verbatim from issue)

DEM 接触モデルに表面粗さ h_r を追加し、見かけ粘度 η_s の計算・出力機能を実装する。
完了後は `physics.surface_roughness` で粒子間の有効接触距離を調整でき、
η_s のタイムシリーズと定常平均を自動出力できる。

## Approach

### 1. Surface roughness (`physics.surface_roughness`)

The particle-particle contact threshold in `DEMSolver._apply_particle_pair_loads` is
currently `min_dist = radii[i] + radii[j]`. Adding roughness means contact triggers when
the centre-to-centre distance `dist < radii[i] + radii[j] + h_r`, where `h_r` is stored
on the `LBMDEMSolver` instance.

- Add `surface_roughness: float = 0.0` parameter to `LBMDEMSolver.__init__` (stored as
  `self.surface_roughness`).
- In `DEMSolver._apply_particle_pair_loads` add `h_r = sim.surface_roughness` and replace
  `min_dist = radii[i] + radii[j]` with `min_dist = radii[i] + radii[j] + h_r` for the
  contact test. Surface-force computation keeps `phys_min_dist = radii[i] + radii[j]` so
  attraction/repulsion gaps are unchanged.
- Add `--surface-roughness` CLI flag to the runner (default 0.0).
- Map `physics.surface_roughness` in the config system (it naturally passes through
  `_flatten_sections` already; just add the CLI arg and wire it).

### 2. Apparent viscosity (`runtime.viscosity_eval`)

The shear viscosity of a suspension under Lees-Edwards shear is:

  η_s = -⟨σ_xy⟩ / γ̇

where σ_xy is the off-diagonal stress, γ̇ = `le_shear_rate`. Two contributions:

- **Hydrodynamic (η^H):** from the fluid stress tensor, summed over all lattice nodes.
  σ_xy^H ≈ -ρ ν (∂ux/∂y + ∂uy/∂x) estimated from the non-equilibrium part of f.
  In LBM with BGK: σ_xy = -(1 - τ/2) Σ_α e_αx e_αy (f_α - f_α^eq).
  Available quantities: the f array lives in `lbm_dem.py`.
- **Collisional (η^C):** impulse-based. For each particle pair collision:
  η^C contribution per step = (force × lever arm) integrated over time / (V * γ̇).

For the first implementation, compute η^H from the fluid stress alone (the dominant term
in dense suspensions); η^C will be accumulated from particle-particle contact forces.

Implementation location: a new `ViscosityEvaluator` helper class in
`src/particulate_flow/rheology.py`.

- Constructor: `ViscosityEvaluator(sim, *, start_step, viscosity_interval, average_steps,
  out_dir)`.
- `record(step, contact_stress_xy)` — called each step (zero-cost when disabled).
- `flush_csv(step)` — writes a row to `viscosity_timeseries.csv` every `viscosity_interval`
  steps.
- `finalize()` → dict with `eta_H_mean`, `eta_C_mean`, `eta_s_mean` over `average_steps`.

In `lbm_dem.py` (step loop) accumulate particle-particle contact forces × separation
vector for η^C. Compute fluid-side η^H using the existing f array.

Runner wires up config keys from `runtime.viscosity_eval.*`:
- `enabled` (bool, default false)
- `start_step` (int, default 0)
- `viscosity_interval` (int, default 100)
- `average_steps` (int, default 1000)

`finalize()` result is merged into `summary.json` under `"apparent_viscosity"`.

### 3. Config plumbing

`runtime.viscosity_eval` is a nested dict; it gets flattened by `_flatten_sections` into
flat keys `viscosity_eval_enabled`, etc. Alternatively, expand it similarly to
`lees_edwards` using a `_expand_runtime_section` helper. The latter is cleaner.

We'll add a `_VISCOSITY_EVAL_KEY_MAP` and expand in `_flatten_sections`.

## Alternatives Considered

- **Modify `min_dist` for both contact AND surface forces:** rejected; surface gap for
  attraction/repulsion should remain geometry-based (touching spheres), not roughness-shifted.
- **In-place η^H computation using macroscopic gradients:** less accurate; using f directly
  avoids storing extra arrays.

## Out of Scope

- 3D suspension; wall contributions to viscosity; non-uniform roughness per particle.
- Numba path for particle contact stress accumulation (pure NumPy is fine here).

## Assumptions

- `le_shear_rate` is accessible on `sim`; confirmed in lbm_dem.py:317.
- f array is accessible as `self.f` in `LBMDEMSolver`; confirmed.
- `viscosity_eval` only produces meaningful output when LE BC is active, but the code does
  not gate on it — the user is responsible (per issue constraint).

## Test Plan

1. **Regression (h_r = 0.0):** run a 2-particle collision with no roughness, verify
   contact threshold matches `r_i + r_j`.
2. **Roughness contact (h_r = 0.05):** verify contact triggers at `dist < 2a + h_r`.
3. **viscosity_eval disabled:** simulate step — verify no extra arrays allocated and no CSV
   created.
4. **viscosity_eval enabled:** run ≥ `start_step + viscosity_interval` steps, verify CSV
   exists with expected columns and ≥ 1 row; verify `finalize()` returns dict with expected
   keys.
5. **average_steps averaging:** verify `eta_s_mean` is the mean of rows in the average
   window.
6. **summary.json integration:** verify `apparent_viscosity` key in summary.json when
   viscosity_eval enabled.
