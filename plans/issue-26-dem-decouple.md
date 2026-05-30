# Plan: Issue #26 — Decouple 2D DEM Contact Law from the Fluid Solver

## Verbatim Goal (from issue #26)

> 2D の DEM 接触物理を `LBMDEMSolver`/`coupled_solver` への依存から切り離し、3D の
> `DEM3D`（流体非依存）と対をなす独立した接触ソルバーに整理する。これにより接触
> 物理が流体ソルバーを構築せずに単体でテスト・再利用できる。

### Acceptance Scenarios
1. 2D の DEM 接触（法線 Hertz・接線摩擦・転がり抵抗）が流体ソルバーを構築せずに
   単体で実行・テストできる
2. 既存の `FastLBMDEM` の本番計算結果が変わらない（数値リグレッションなし）
3. `DEM3D` と 2D 接触ソルバーが同じ接触則の概念（clamp 付き法線・Coulomb 接線・
   転がり抵抗）を共有する形に近づく

### Constraints
- Package: particulate_flow
- `FastLBMDEM` の挙動を壊さない
- 公開インターフェース（`advance`, `get_fields`, `dimensionless_groups`）を維持
- 接触則の物理は変えない（純粋な構造リファクタ、数値結果は不変）
- スコープは DEM 接触の切り出しのみ（境界条件・IBM・粒子注入・統計の分割は別 issue）

## Scope decision

`DEMSolver` reads ~30 distinct `sim.X` attributes — material params, particle state
arrays, geometry, AND callbacks (`_periodic_y_delta`, `_particle_pair_candidates`,
`_particle_drag_forces`, `ibm_forces_p`).  Fully decoupling `compute_loads` would
require reworking the pair-search/drag/IBM wiring — that is the broader God-class
split, explicitly **out of scope** here.

The cleanly extractable, fluid-independent core is the **contact-law magnitudes**:
`normal_contact_magnitude`, `tangential_force_magnitude`, `rolling_resistance_torque`.
These are pure functions of (contact state) + (7 scalar material params:
`k_n, damping, sliding_friction, tangential_damping, rolling_friction,
rolling_friction_coeff, rolling_damping`) + `contact_model`.  They are exactly the 2D
analogue of the `DEM3D._normal_magnitude` / `_tangential_magnitude` / `_rolling_torque`
methods.

### In scope (this PR)
- **New module `dem/contact_law.py`** with a `ContactLaw` dataclass holding the 8
  contact parameters and three methods (`normal_magnitude`, `tangential_magnitude`,
  `rolling_torque`) — **no fluid-solver reference**, constructible and testable
  standalone.  The formulas are copied **verbatim** from `DEMSolver` (clamp-to-0
  normal Hertz/linear, Coulomb-limited tangential, Coulomb-limited rolling).
- **`DEMSolver` delegates** its three contact-law methods to a `ContactLaw` built from
  `self.sim`'s params.  The public method signatures stay identical (back-compat); the
  bodies become thin wrappers.  No change to `compute_loads` orchestration or any other
  `sim.X` access — only the three pure methods are redirected.
- **A `ContactLaw.from_solver(sim)` constructor** that reads the 8 params off any object
  exposing them (the solver), so `DEMSolver` builds it in one line.
- **Tests** (`tests/test_contact_law.py`): the three laws run and are correct
  **without constructing any fluid solver** (scenario 1); they match the 3D `DEM3D`
  concepts (scenario 3); and a numerical-equivalence guard that `DEMSolver`'s methods
  return identical values to a `ContactLaw` with the same params (scenario 2 support).
- **Regression**: the full existing suite stays green (scenario 2).

### Out of scope (other issues)
- Decoupling `compute_loads`, pair search, drag, IBM from `sim` (broader God-class
  split — future issues).
- Unifying 2D and 3D into one dimension-agnostic contact module (that is the larger
  dimension-agnostic effort; here we only make the 2D core fluid-free and *conceptually*
  aligned with `DEM3D`).
- Any change to `DEM3D`.

## Approach

`ContactLaw` (frozen dataclass):
```
fields: contact_model, k_n, damping, sliding_friction, tangential_damping,
        rolling_friction, rolling_friction_coeff, rolling_damping
normal_magnitude(overlap, v_n, mass)        # == DEMSolver.normal_contact_magnitude
tangential_magnitude(v_t, normal_force, mass)  # == DEMSolver.tangential_force_magnitude
rolling_torque(omega, normal_force, radius, mass)  # == DEMSolver.rolling_resistance_torque
from_solver(sim) -> ContactLaw              # read the 8 params off sim
```
`DEMSolver.__init__` builds `self._contact_law = ContactLaw.from_solver(coupled_solver)`
and the three methods delegate.  Formulas are unchanged → numbers are bit-identical.

## Alternatives considered & rejected
- **Full `compute_loads` decoupling**: rejected — pulls in pair-search/drag/IBM rewiring,
  which is the broader split and out of scope (high regression risk in one PR).
- **Move the 3 methods onto a shared base with `DEM3D`**: rejected for now — 2D uses
  scalar omega/torque, 3D uses vectors; unifying needs the dimension-agnostic rework.
  This PR aligns them *conceptually* (same `ContactLaw` formulas) without forcing that.
- **Make `DEMSolver` itself fluid-free**: rejected — it legitimately needs particle
  state + pair search + drag from the solver; that is the God-class split, not this slice.

## Assumptions about other modules (docstring-grounded)
- `DEMSolver` "reads particle and geometry arrays from the coupled simulation object
  and asks that object for fluid drag" (solver docstring) — only the 3 pure contact-law
  methods are redirected; the array/drag reads are untouched.
- `DEM3D._normal_magnitude/_tangential_magnitude/_rolling_torque` (contact3d.py
  docstrings) define the same clamp/Coulomb laws — `ContactLaw` mirrors them in 2D form.
- `FastLBMDEM` inherits `DEMSolver`'s contact path indirectly via `_dem_*`; the delegated
  methods return identical values, so production results are unchanged.

## Testing plan (acceptance → test) — `tests/test_contact_law.py`
- **Scenario 1** (standalone): build `ContactLaw(...)` directly (no solver) and assert
  normal Hertz (`k_n·overlap^1.5`) clamps at 0 under strong damping; tangential is
  Coulomb-bounded `±μ·f_n`; rolling torque is bounded `±rolling_friction_coeff·f_n·r`;
  all return 0 when `rolling_friction=False` or `normal_force<=0`.
- **Scenario 2** (equivalence / no regression): for a range of inputs, assert a
  `DEMSolver` instance's `normal_contact_magnitude` / `tangential_force_magnitude` /
  `rolling_resistance_torque` equal a `ContactLaw` with the same params — proving the
  delegation is behaviour-preserving.  Plus the full suite stays green.
- **Scenario 3** (concept parity): assert `ContactLaw` and `DEM3D` agree on the scalar
  contact-law results for an equivalent 1-D (head-on / tangential) configuration
  (normal magnitude and Coulomb bounds match), documenting the shared concept.
- A lightweight construction test that `ContactLaw.from_solver` reads the 8 params.
