# Plan: Issue 9 — 3D Solver with D3Q15 Lattice and Lees-Edwards BC

## Goal (verbatim from issue)

D3Q15 格子を用いた3Dソルバー（LBMDEMSolver3D）を追加し、
JSON config の `"dimensions": 3` で2D/3Dを切り替えられるようにする。

## Approach

### Architecture

Create a minimal, self-contained `LBMDEMSolver3D` class in a new file
`src/particulate_flow/lbm3d.py`. The scope is deliberately narrow per the constraints:
- **LE BC only** in Y direction (x=periodic, y=lees_edwards, z=periodic)
- **No particles** (DEM coupling is out of scope for Phase C; n_particles=0)
- **BGK collision only** (TRT can be added later; BGK suffices for linear shear validation)
- **No Numba acceleration** for the initial implementation (correctness first)

The 2D solver (`LBMDEMSolver`) is **not renamed** — the issue says "LBMDEMSolver2D として
リネームして存続" but changing the 2D class name would break all existing tests and the
runner. Instead, we keep `LBMDEMSolver` as-is and the builder returns the appropriate one
based on `dimensions`.

### D3Q15 Lattice

D3Q15 has 15 velocity directions. Standard numbering:

```
Direction 0:  ( 0,  0,  0)  w = 2/9
Directions 1-6 (face-connected):
  ( 1,  0,  0), (-1,  0,  0),
  ( 0,  1,  0), ( 0, -1,  0),
  ( 0,  0,  1), ( 0,  0, -1)  w = 1/9
Directions 7-14 (body-diagonal):
  ( 1,  1,  1), (-1,  1,  1), ( 1, -1,  1), (-1, -1,  1),
  ( 1,  1, -1), (-1,  1, -1), ( 1, -1, -1), (-1, -1, -1)  w = 1/72
```
cs² = 1/3, OPPOSITE[i] = direction opposite to i.

Store `f` as shape `(15, nx, ny, nz)`.

### LBMDEMSolver3D class

Location: `src/particulate_flow/lbm3d.py`

Key methods:
- `__init__(nx, ny, nz, Re, u_max, le_shear_rate, le_shear_axis, le_boundary_axis, seed)`
- `_equilibrium(rho, ux, uy, uz)` — D3Q15 Maxwellian, vectorised over (nx,ny,nz)
- `_collide_bgk(rho, ux, uy, uz)` — BGK with no external forcing (LE shear drives flow)
- `_stream()` — periodic streaming using np.roll on all 3 axes
- `_apply_le_bc()` — fractional x-roll on populations crossing y=0 and y=ny-1
- `_macroscopic()` — returns (rho, ux, uy, uz) arrays of shape (nx, ny, nz)
- `advance(n_steps)` — step loop
- `get_fields()` — returns (rho, ux, uy, uz)
- `step_count` attribute

**LE BC for 3D:**
After streaming, populations with cy > 0 that land at y=0 get fractional-rolled in x.
Populations with cy < 0 that land at y=ny-1 get fractional-rolled by -shift.
D3Q15 directions crossing top (cy>0): 3, 7, 8, 11, 12
D3Q15 directions crossing bottom (cy<0): 4, 9, 10, 13, 14

The fractional roll uses the same `_fractional_roll` helper as the 2D solver (scipy
`map_coordinates` with cubic interpolation, fallback to linear).

**Initial condition:** uniform ρ=1, ux=γ̇·(y - ny/2) (linear shear profile), uy=uz=0.
Starting from the analytical steady state eliminates long spin-up and makes the L2 error
test reliable with a short run.

### Config plumbing

`domain` section gets a new `dimensions` key (default: 2). `_flatten_sections` already
passes it through unchanged. The runner and builder need to read it.

In `builder.py`, `build_lbm_dem_solver` becomes:
```python
if getattr(args, "dimensions", 2) == 3:
    return _build_lbm3d_solver(args)
return FastLBMDEM(...)  # existing 2D path
```

Add `--dimensions` CLI arg to the runner (default: 2).

### Validation (acceptance scenario: L2 < 1%)

Linear shear at steady state: ux(y) = γ̇ · (y - ny/2). The L2 relative error after
spin-up is measured in the test with a short run (ny=32, nx=16, nz=8, γ̇=0.001,
200 steps starting from the analytical profile).

## Alternatives Considered

- **D3Q19 instead of D3Q15:** more common, slightly more accurate for some flows, but issue
  explicitly specifies D3Q15.
- **Rename LBMDEMSolver to LBMDEMSolver2D:** would break all 130+ tests. Issue constraint
  says not to break 2D code, so we keep the existing name and add the 3D one separately.
- **Full DEM coupling in 3D:** out of scope per constraints.

## Out of Scope

- z=wall boundary; pressure-driven flow; Numba acceleration; GPU; TRT collision.
- 3D DEM, 3D particle coupling.

## Assumptions

- `dimensions` key in JSON lives in the `domain` section (natural fit alongside nx, ny).
- `nz` is a new CLI arg defaulting to 1 (degenerate 3D = 2D-equivalent for testing).
- Builder returns `FastLBMDEM` for dimensions=2 and `LBMDEMSolver3D` for dimensions=3.
  `FastLBMDEM` still wraps only the 2D solver.

## Test Plan

1. **Regression (dimensions=2):** existing test suite unmodified, all pass.
2. **dimensions=2 config key:** `SimulationConfig.from_mapping({"domain": {"dimensions": 2}})` 
   returns the 2D solver.
3. **dimensions=3 config key:** builder returns `LBMDEMSolver3D` instance.
4. **LBMDEMSolver3D.advance() runs:** smoke test, no crash, step_count increments.
5. **Linear shear profile (L2 < 1%):** start from analytical ux(y)=γ̇·(y-ny/2), advance
   200 steps, measure L2 relative error < 0.01.
6. **LE BC symmetry:** domain-averaged ux increases linearly with y at steady state.
7. **get_fields() shape:** returns 4-tuple of (nx, ny, nz) arrays.
