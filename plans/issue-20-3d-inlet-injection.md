# Plan: Issue #20 — 3D Left-Inlet Particle Injection/Removal

## Verbatim Goal (from issue #20)

> 3D `left_inlet` particle source for `LBMDEMSolver3D`: inject particles across
> the inlet plane per the flow-rate budget and delete them after the outlet,
> mirroring the 2D left-inlet logic.

### Acceptance Scenarios
1. Particles appear at the inlet plane over time and advect downstream
2. Particles past the outlet are removed
3. Injected volume fraction tracks the configured source fraction

### Constraints
- Package: particulate_flow
- Depends on #14 (3D pressure BC), #17 (IBM), #18 (DEM) — all **merged**

## Scope decision

The final slice of the 3D fouling stack.  It makes the particle population
*dynamic*: particles are injected at the inlet (x≈0) plane over time according to
an area(volume)-flux budget, and removed once they pass the outlet (x=nx-1).  The
fluid/IBM/DEM machinery from #14/#17/#18/#19 is reused unchanged per step; only
the particle array grows/shrinks between steps.

### In scope (this PR)
- **Dynamic particle arrays in `DEM3D`.** Add `add_particles(pos, vel, radii)` and
  `remove_particles(mask)` that grow/shrink `pos, vel, omega, radii, masses,
  inertias, n_p` consistently.  This is the seam the fixed-count `DEM3D` lacks.
- **`particle_source="left_inlet"` in `LBMDEMSolver3D`.** Lift the #14 stub.  With
  a source fraction `phi`, each step:
  1. accumulate budget `+= phi · inlet_flow_rate` (3D inlet flux = Σ positive ux
     over the inlet y-z plane, mirroring 2D `_left_boundary_flow_rate`),
  2. while budget ≥ a particle's volume, place it just inside the inlet at a
     flux-weighted (y, z) sample with no overlap, and decrement the budget,
  3. after the DEM step, delete particles whose sphere has fully passed x=nx-1.
- **Volume budget in 3D.** Use sphere volume `(4/3)π r³` (3D analogue of the 2D
  disc-area budget) so the injected solid-volume fraction tracks `phi`.
- **IBM array sync.** The coupled `advance()` already reads `self.dem.pos/vel/...`
  fresh each step and rebuilds markers per particle, so a changing `n_p` is handled
  automatically as long as `_ibm_reaction`/`_ibm_reaction_torque` are sized to the
  current `n_p` (they are, since `_apply_ibm` allocates `(n, 3)` each call).
- Tests for the 3 acceptance scenarios + DEM3D add/remove unit tests + regression.

### Out of scope (other issues / future)
- Per-particle radius distributions (use a single `particle_radius`, as the rest of
  the 3D path does; the 2D queue supports variation but that is not required here).
- Numba acceleration.
- Inlet injection in the periodic/Lees-Edwards mode (left-inlet only makes sense
  with the pressure streamwise boundary; guard against the combination).

## Approach

### DEM3D growth/shrink
`add_particles(pos, vel, radii)`: vstack/append onto every state array; new `omega`
rows are zero; recompute `masses`/`inertias` for the appended radii; bump `n_p`.
`remove_particles(mask)`: boolean-keep every array; update `n_p`.  Both are no-ops
on empty input.

### Inlet flux + placement (3D)
- `inlet_flow_rate = Σ max(ux[0, :, :], 0)` over the inlet plane (area flux).
- Budget accumulates `phi · inlet_flow_rate`; a sphere of radius r costs
  `(4/3)π r³`.
- Placement: x = r + 0.75; sample (y, z) weighted by the local inlet ux (fallback
  to a deterministic sweep before flow develops); reject if it overlaps an existing
  particle or a cylinder; seed velocity from the local inlet ux.
- Use a seeded RNG (`np.random.default_rng`) for reproducibility, like 2D.

### advance() integration
At the top of each step (pressure mode + left_inlet): feed the inlet (grow `dem`),
run the existing coupled IBM+LBM+DEM step, then delete outflow particles.  Keep the
fixed-population path unchanged when `particle_source != "left_inlet"`.

## Alternatives considered & rejected
- **Pre-place all particles at t=0 (initial source):** rejected — the issue wants
  *gradual* injection tracking the flow-rate budget.
- **A separate inlet manager class:** rejected — the state lives in `DEM3D`; adding
  two small array ops there plus a feed method on the solver is simpler and mirrors
  the 2D in-solver approach.
- **Per-particle radius queue (2D-style pending arrays):** rejected as unnecessary —
  the 3D path uses a single radius; a uniform-radius budget is sufficient and
  simpler. (Left as a noted future extension.)

## Assumptions about other modules (docstring-grounded)
- 2D `_try_feed_left_inlet_particles`: budget `+= phi · inlet_flow_rate`, place when
  budget ≥ particle area, delete when `pos_x - r > nx` (lbm_dem.py:683-756). 3D
  mirrors with sphere volume and the y-z inlet plane.
- 2D `_left_boundary_flow_rate`: Σ positive ux at x=0 over the cross-section. 3D sums
  over the y-z plane.
- `DEM3D` (contact3d.py docstring): owns `pos/vel/omega/radii/masses/inertias/n_p`;
  `compute_loads`/`step` read those each call — so add/remove between steps is safe.
- `LBMDEMSolver3D._apply_ibm` allocates `(n_p, 3)` reaction arrays each call, so a
  changing `n_p` needs no extra bookkeeping.

## Testing plan (acceptance → test) — `tests/test_lbm3d_inlet.py`
- **Scenario 1** (appear + advect): start with `n_particles=0`, `left_inlet`,
  pressure flow; after enough steps `dem.n_p > 0`, and the mean particle x increases
  over time (downstream advection); all finite.
- **Scenario 2** (removal): run long enough that early particles cross the outlet;
  assert the running count stops growing unbounded and that no particle persists with
  `x - r > nx` (outflow deleted); cumulative removed > 0.
- **Scenario 3** (phi tracking): with a higher `phi`, more solid volume is injected
  over the same number of steps than with a lower `phi` (monotone), and the injected
  volume ≈ `phi · cumulative_inlet_flux` within tolerance.
- **DEM3D add/remove unit tests**: `add_particles` grows all arrays consistently
  (shapes, masses/inertias recomputed, omega zero); `remove_particles(mask)` keeps
  the complement; round-trip leaves a consistent state; empty ops are no-ops.
- **Guard**: `left_inlet` with periodic streamwise boundary raises (only valid with
  pressure flow).
- **Regression**: fixed-population coupling (`n_particles>0`, no inlet) unchanged;
  full suite green.
