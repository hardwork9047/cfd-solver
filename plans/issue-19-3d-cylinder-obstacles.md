# Plan: Issue #19 — 3D Fixed Finite-Cylinder Solid Obstacles

## Verbatim Goal (from issue #19)

> 3D fixed finite-length cylinder obstacles as solid masks in `LBMDEMSolver3D`,
> with bounce-back at solid faces, mirroring the 2D cylinder support.

### Acceptance Scenarios
1. A finite cylinder placed in 3D pressure flow blocks fluid (no-slip) and deflects it
2. Particles collide with the fixed cylinder (needs 3D DEM sibling)

### Constraints
- Package: particulate_flow
- Depends on #14 (3D pressure BC) being merged — **DONE**

## Scope decision

Two halves, both already partly in place:
- **Fluid side (new here):** build a 3D solid mask for z-aligned finite cylinders
  and apply halfway bounce-back on solid nodes after streaming, mirroring the 2D
  `f[i, solid] = f[opposite[i], solid]`.  This is the piece #14 stubbed with
  `NotImplementedError` for `cylinders`.
- **Particle side (already done in #18):** `DEM3D` already collides spheres with
  z-aligned finite cylinders.  This PR wires the same cylinder list into the
  internal `DEM3D` so particles and fluid see the same obstacles.

### In scope (this PR)
- Lift the `cylinders` `NotImplementedError` guard in `LBMDEMSolver3D.__init__`.
- Accept `cylinders` as `(cx, cy, r)` or `(cx, cy, r, z_lo, z_hi)`; build a 3D
  boolean `self.solid` mask `(x-cx)² + (y-cy)² ≤ r²` within `[z_lo, z_hi]`
  (full z by default).
- Halfway bounce-back on `self.solid` after streaming, before the BC step, in
  `advance()`.
- Make the IBM spreading skip solid nodes (`active = ~solid`) so reaction force is
  not deposited inside obstacles, mirroring 2D `_distribute_lagrangian_forces`.
- Pass the same `cylinders` into the internal `DEM3D` so particles collide with them.
- Tests for both acceptance scenarios + regression (no cylinders ⇒ unchanged).

### Out of scope (other issues)
- Inlet particle injection → #20.
- Numba acceleration.
- Arbitrary-orientation cylinders (only z-aligned, matching `DEM3D` and the 2D model
  which is effectively z-extruded discs).

## Approach

### Solid mask
`_build_cylinder_solid(cylinders)` → bool array (nx, ny, nz).  For each cylinder,
mark cells whose x-y radial distance ≤ r and whose z is within `[z_lo, z_hi]`.
Stored as `self.solid`; `self.solid.any()` gates the bounce-back cost.

### Bounce-back
After `_stream()`, when `self.solid.any()`:
`f_tmp = f.copy(); for i: f[i, solid] = f_tmp[OPPOSITE3[i], solid]`.
This enforces no-slip on the obstacle surface (halfway bounce-back), identical in
form to the 2D path.  Applied before the pressure/LE boundary step (same order as 2D).

### IBM interaction with solids
In `_spread_forces_3d`, mask out solid target nodes so spread force lands only on
fluid cells (mirrors 2D `active = ~self.solid[ix, iy]`).  Interpolation reads
velocity at markers on particle surfaces (in fluid), so no change needed there; the
solid cells carry the post-bounce-back populations whose velocity is ~0 (no-slip),
which is the physically correct thing to interpolate near an obstacle anyway.

### DEM coupling
Forward `cylinders` to the internal `DEM3D(..., cylinders=...)` so spheres collide
with the obstacles (the contact physics already exists from #18).

## Alternatives considered & rejected
- **Interpolated (curved) bounce-back:** rejected — the 2D model uses simple halfway
  bounce-back on a staircased mask; mirroring it keeps 2D/3D consistent and the issue
  says "mirroring the 2D cylinder support".
- **Reuse `PoreGeometry`:** rejected — `PoreGeometry.cylinder_solid_mask` is 2D
  `(nx, ny)`; 3D needs the z-extent, so a small local builder is clearer than bending
  the 2D geometry class.
- **Arbitrary cylinder axis:** rejected — out of scope; `DEM3D` cylinders are
  z-aligned, so the fluid mask matches.

## Assumptions about other modules (docstring-grounded)
- 2D bounce-back (`lbm_dem.py` _lbm_step): `f[i, solid] = f_tmp[OPPOSITE[i], solid]`
  after streaming — 3D mirrors with `OPPOSITE3`.
- `DEM3D` docstring: cylinders are "z-aligned finite cylinders as (cx, cy, r) or
  (cx, cy, r, z_min, z_max); collide in the x-y plane within the z-extent" — the
  fluid mask uses the identical convention so fluid and particles agree.
- `LBMDEMSolver3D.advance/get_fields` signatures unchanged (issue #14 constraint).

## Testing plan (acceptance → test) — `tests/test_lbm3d_cylinder.py`
- **Scenario 1a** (no-slip): a cylinder in 3D pressure flow → velocity inside the
  solid mask ≈ 0; assert `|u|` on solid cells ~0 to tolerance after running.
- **Scenario 1b** (deflection): the flow speeds up around the cylinder / is diverted —
  assert the transverse velocity (uy) near the cylinder flanks is nonzero, and the
  centreline x-velocity just downstream is reduced vs the no-cylinder case.
- **Scenario 2** (particle collision): a particle advected toward a fixed cylinder
  does not penetrate it (centre stays ≥ r_cyl + r_p from the axis within the z-extent).
- **Mask geometry**: `_build_cylinder_solid` marks the right cells (inside radius &
  z-extent true; outside false; finite z respected).
- **Regression**: with no cylinders, `self.solid` is all-False and results match the
  pre-#19 behaviour (no bounce-back path taken); full suite stays green.
