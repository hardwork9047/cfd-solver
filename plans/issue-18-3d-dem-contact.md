# Plan: Issue #18 — 3D DEM Contact Physics

## Verbatim Goal (from issue #18)

> Fresh 3D DEM contact physics for `LBMDEMSolver3D` (the 2D `DEMSolver` is
> 2D-hardwired: `compute_contact_stress_xy`, `tangent_from_normal`).
> Sphere-sphere, sphere-wall, and sphere-cylinder normal+tangential Hertz
> contact with gravity and time integration in 3D.

### Acceptance Scenarios
1. Two approaching spheres in 3D repel via Hertz normal contact
2. A sphere settles under gravity onto a wall and rests
3. Tangential friction + rolling resistance act in 3D

### Constraints
- Package: particulate_flow
- Depends on #14 (3D pressure BC) being merged — **DONE** (PR #21 squashed to main)
- Do not regress 2D DEM

## Scope decision

This issue is the **contact-physics core only** — pure DEM (particles + walls +
cylinders + gravity + integration) in 3D. It is intentionally decoupled from the
fluid: fluid→particle drag and particle→fluid IBM back-reaction are issue **#17**.
So this PR delivers a standalone, fluid-free 3D DEM stepper that the IBM work (#17)
will later wire into `LBMDEMSolver3D.advance()`.

### In scope (this PR)
- New module `src/particulate_flow/dem/contact3d.py` with a `DEM3D` class that owns
  3D particle state (pos, vel, omega as 3-vectors, radii, masses) and computes
  + sphere-sphere Hertz normal contact (f_n = k_n · overlap^1.5) with velocity damping
  + tangential Coulomb-limited friction along the contact-plane slip direction
  + rolling resistance as a torque vector opposing relative angular velocity
  + gravity (buoyancy-corrected, matching 2D `density_ratio` convention)
  + sphere-wall contact for the y=0 / y=ny-1 planes (configurable)
  + sphere-(finite z-axis-aligned) cylinder contact
  + velocity-Verlet / semi-implicit Euler time integration with sub-stepping
- A thin opt-in hook on `LBMDEMSolver3D` (e.g. `attach_dem(...)`) that constructs a
  `DEM3D` and exposes `step_dem(n_substeps)` WITHOUT touching the fluid path. The
  existing `NotImplementedError` particle-stub in `__init__` stays for the fluid-coupled
  config route (that path is #17); `attach_dem` is the explicit, fluid-free entry.
- Tests for the three acceptance scenarios + 2D-regression guard.

### Out of scope (other issues)
- Fluid drag + IBM back-reaction coupling → #17
- Pressure-driven advection of particles, inlet injection → #20
- 3D solid masks for fixed cylinders in the FLUID (#19); here cylinders act only as
  DEM collision bodies, not fluid obstacles.
- Numba acceleration of the 3D kernels (NumPy first; optimise later if needed).

## Approach

### Physics (mirrors 2D `DEMSolver`, lifted to 3D vectors)
- **Normal:** `n = (p_j - p_i)/|p_j - p_i|`; `overlap = (r_i + r_j) - dist`; if > 0,
  `f_n = max(k_n·overlap^1.5 - damping·v_n·sqrt(k_n·m_eff), 0)` along `n`.
  (Clamp to ≥0 to avoid artificial tension, exactly as 2D.)
- **Tangential:** relative surface velocity
  `v_rel = (v_j - v_i) + (omega_j × (-r_j n) - omega_i × (r_i n))`;
  tangential component `v_t = v_rel - (v_rel·n)n`; direction `t = v_t/|v_t|`;
  `f_t = clip(-tangential_damping·sqrt(k_n·m)·|v_t|, -μ·f_n, μ·f_n)·t`.
  Apply `±f_t` to the pair and torque `τ_i += (r_i n) × f_t_on_i`.
- **Rolling resistance:** `τ_roll = clip(-rolling_damping·sqrt(k_n·m)·r²·ω_rel,
  ±rolling_friction_coeff·f_n·r)` opposing the relative angular velocity vector.
- **Wall (y-planes):** treat as an infinite plane; overlap with radius; same normal
  + tangential + rolling formulas with the wall normal ±ŷ and wall velocity 0.
- **Cylinder (z-aligned, finite):** radial distance in the x-y plane to the cylinder
  axis; collide when within `r_cyl + r_p` and within the cylinder's z-extent.
- **Gravity:** `f -= m·g·(1 - 1/density_ratio)·ĝ` (ĝ default -y), matching 2D.
- **Integration:** semi-implicit Euler with `dem_substeps`, as 2D uses.

### Reuse without regression
- I will NOT modify the 2D `DEMSolver`. The 3D code is a separate module. The 2D
  scalar-torque formulas are special cases of the 3D cross-product ones, but sharing
  would force a 2D rewrite and risk regression — rejected (see Alternatives).
- Constants/material params (`k_n`, `damping`, `sliding_friction`,
  `rolling_friction_coeff`, `rolling_damping`, `density_ratio`) reuse the same names
  and semantics as 2D so configs stay consistent.

## Alternatives considered & rejected
- **Generalise `DEMSolver` to N-D**: rejected — it is pervasively 2D (scalar `omega_p`,
  `tangent_from_normal`, `forces[:,1]`, numba kernels with 2-column arrays). Refactoring
  risks regressing the production 2D path; the issue explicitly says "fresh 3D".
- **Couple to fluid in this PR**: rejected — that is #17; keeping DEM fluid-free makes
  this PR testable in isolation and reviewable.
- **Numba from the start**: rejected — premature; correctness first, NumPy is fine for
  test sizes.

## Assumptions about other modules (docstring-grounded)
- 2D `DEMSolver.normal_contact_magnitude` docstring: "Normal contact force with
  damping, clamped to avoid artificial tension" and uses `k_n·overlap^1.5` for
  dem-hertz — 3D mirrors this clamp + exponent.
- 2D `tangential_force_magnitude`: "Coulomb-limited tangential force opposing slip" —
  3D uses the same Coulomb clip `±μ·f_n`.
- 2D `rolling_resistance_torque`: "Coulomb-limited rolling resistance torque opposing
  angular velocity" — 3D uses the same limit `rolling_friction_coeff·f_n·r`.
- `LBMDEMSolver3D` exposes `nx, ny, nz`; the `__init__` particle stub raises
  `NotImplementedError` for the fluid-coupled config route (added in #14) — `attach_dem`
  is a distinct, explicit, fluid-free path that does not go through that stub.

## Testing plan (acceptance → test) — `tests/test_dem3d.py`
- **Scenario 1** (two spheres repel): place two overlapping spheres at rest on the
  x-axis, step DEM, assert they separate (centre distance increases; x-velocities point
  apart; no NaN). Also assert NO force when separated beyond contact.
- **Scenario 2** (settle on wall): one sphere above the y=0 wall under gravity; step
  many sub-steps; assert it comes to rest at ≈ radius above the wall (penetration small,
  velocity → ~0), and does not tunnel through.
- **Scenario 3** (friction + rolling): a sphere with tangential surface velocity against
  the wall; assert a tangential force opposes slip (Coulomb-bounded) and a rolling
  resistance torque opposes spin; assert both vanish when normal force is zero.
- **Energy/sanity**: a head-on equal-mass elastic-ish collision conserves momentum to
  tolerance; fields stay finite.
- **2D regression**: existing `tests/` (incl. `test_lbm_dem`, dem suites) stay green;
  add an explicit assertion that importing `contact3d` does not perturb 2D `DEMSolver`.
- **Cylinder**: a sphere driven toward a z-aligned cylinder is deflected (no penetration
  inside `r_cyl + r_p`).
