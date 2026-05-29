# Plan: Issue #17 — 3D IBM Particle-Fluid Coupling

## Verbatim Goal (from issue #17)

> 3D Immersed Boundary Method coupling for `LBMDEMSolver3D`: distribute particle
> reaction force to lattice nodes and interpolate fluid velocity at particle
> markers, mirroring the 2D `_ibm_particle_reaction` path.

### Acceptance Scenarios
1. A particle in a 3D pressure-driven flow experiences fluid drag
2. Newton's 3rd-law back-reaction is applied to the fluid body force
3. 2D IBM behaviour unchanged

### Constraints
- Package: particulate_flow
- Depends on #14 (3D pressure BC) being merged — **DONE**
- Requires 3D DEM contact (#18) for full particle motion — **DONE** (`DEM3D`)

## Scope decision

This is the **integration capstone** for the 3D fouling slice: it wires the
fluid-free `DEM3D` (#18) into the `LBMDEMSolver3D` fluid step (#14) with two-way
coupling.  The coupling is the classic direct-forcing IBM lifted to 3D, mirroring
the 2D `_apply_immersed_boundary_forces` path.

### In scope (this PR)
- **3D Guo-style body force in the LBM step.** The current 3D BGK collision has no
  forcing term.  Add a `Fx, Fy, Fz` Eulerian body-force field and a D3Q15 Guo
  forcing source term so spread IBM reaction forces actually drive the fluid.
- **3D Lagrangian markers.** Distribute markers over each sphere's surface
  (spherical Fibonacci/lat-long lattice) with surface velocity
  `ub = v + ω × r`.
- **Spread + interpolate (3D trilinear).** Interpolate fluid velocity at marker
  points (trilinear) and spread marker forces back to the 8 surrounding lattice
  nodes, with periodic y, z and the x in/outlet handled like the 2D code handles
  its boundaries.
- **Direct-forcing law:** `f_marker = stiffness · (ub - uf) · ds`; particle
  reaction `-Σ f_marker` (force) and `-Σ r × f_marker` (torque vector).
- **A coupling driver on `LBMDEMSolver3D`.** Replace the construction path: when
  `particle_fluid_coupling="immersed_boundary"` with particles, build an internal
  `DEM3D` and run a coupled `advance()`:
  1. compute macroscopic fields,
  2. IBM: markers → interp uf → marker forces → particle reaction + spread to `F`,
  3. LBM collide (with Guo forcing) → stream → boundary BC,
  4. DEM sub-step using the IBM reaction force as the external fluid load.
- **`DEM3D` external-force hook.** Add an optional `external_forces`/
  `external_torques` argument to `DEM3D.step()` (per-particle (n,3)) so the coupling
  layer can inject the IBM reaction without `DEM3D` knowing about fluid.  This is the
  small, clean seam #18 deliberately deferred.
- Tests for the 3 acceptance scenarios + regression.

### Out of scope (other issues)
- Pressure-driven inlet particle injection / removal → #20 (here particles are
  placed initially and tracked; no inlet feed).
- 3D fixed-cylinder FLUID solid masks → #19 (cylinders still act only as DEM bodies).
- Numba acceleration of the 3D IBM kernels (NumPy first).
- iSP / point-force coupling modes in 3D (only `immersed_boundary` here).

## Approach detail

### Guo forcing for D3Q15
Velocity is shifted by half a force: `u = (Σ c_i f_i + F/2)/ρ`.  Collision adds
`S_i = (1 - 1/(2τ)) w_i [ (c_i-u)/cs² + (c_i·u)/cs⁴ c_i ] · F`.  This is the
standard Guo term; in 3D `c_i` and `u` are 3-vectors.

### Markers on a sphere
For radius r, choose `n ≈ max(12, ceil(4π r² / ds²))` markers via a spherical
Fibonacci lattice (even area coverage); each marker area `ds² ≈ 4π r² / n`.  The
spread/interp weight uses the trilinear (1 - |dx|)(1 - |dy|)(1 - |dz|) stencil.

### Coupled advance
Mirror 2D ordering: IBM force is computed from the pre-collision field and added to
`F`; the LBM collision consumes `F`; the same per-particle reaction drives the DEM
sub-steps.  DEM contact (walls/cylinders/pairs) comes from `DEM3D.compute_loads`,
with the IBM reaction added as `external_forces`.

## Alternatives considered & rejected
- **Body-force via simple velocity shift (no Guo):** rejected — Guo is the
  consistent second-order forcing already used in 2D (`_guo_forcing`); matching it
  keeps the physics comparable across 2D/3D.
- **Make `DEM3D` aware of the fluid:** rejected — breaks #18's clean fluid-free
  separation. The `external_forces` seam keeps `DEM3D` agnostic.
- **One marker at the centre (point force):** rejected — issue asks for IBM markers
  distributed over the surface, mirroring 2D.

## Assumptions about other modules (docstring-grounded)
- 2D `_apply_immersed_boundary_forces` docstring: "Direct-forcing IBM for circular
  DEM particles" — 3D mirrors with spherical markers.
- 2D marker force `stiffness·(ub-uf)·ds` and reaction `-f`, torque `r×f`
  (lbm_dem.py:1343-1359) — 3D uses the same law with 3-vectors.
- `DEM3D.step` (contact3d.py docstring): "semi-implicit Euler … dt_sub = 1/dem_substeps"
  — I extend it with optional external loads without changing that integration.
- `LBMDEMSolver3D.advance/get_fields` signatures must stay unchanged (issue #14
  constraint) — coupling lives inside `advance`'s loop; no signature change.

## Testing plan (acceptance → test) — `tests/test_lbm3d_ibm.py`
- **Scenario 1** (drag): place one neutrally/heavier sphere in a 3D pressure-driven
  flow; run; assert the particle accelerates downstream (vx > 0, trends toward the
  local fluid velocity), and that with zero pressure drop and quiescent fluid it does
  not spontaneously move.
- **Scenario 2** (back-reaction): assert that after an IBM step the fluid body-force
  field `F` is nonzero where markers sit and that total IBM reaction on particles
  equals `-Σ` spread force (momentum book-keeping: |Σ particle reaction + Σ spread|
  ≈ 0 to tolerance).
- **Scenario 3** (2D unchanged): full existing suite stays green; explicit assertion
  that the 2D `_apply_immersed_boundary_forces` and `DEMSolver` are untouched.
- **DEM3D external-force hook**: a standalone `DEM3D.step(external_forces=F)` moves a
  particle per the injected force, and omitting it reproduces the prior behaviour
  (regression for #18 tests).
- **Stability/finiteness**: coupled run stays finite and subsonic over many steps.
