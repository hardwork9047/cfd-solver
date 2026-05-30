# Findings: issue-15-3d-config-runner

## Round 1

### Code Reviewer
- Blocking: 0
- Should Fix: 3
- Informational: 6
- Key Issues:
  - **Should Fix** — `_scalar_row` called `sim.get_fields()` twice (full O(nx·ny·nz)
    recompute) — once discarding rho, once for `rho_mean`. **Fixed**: capture rho once.
  - **Should Fix** — `streamwise_boundary="periodic-force"` (a valid 2D CLI value) would
    crash deep in `LBMDEMSolver3D` with a bare `ValueError`. **Fixed**: `build_3d_solver`
    maps `periodic-force → periodic` before constructing the solver.
  - **Should Fix** — tests passed `output_root` but `run_3d` ignored it, writing into the
    repo's `results/` tree on every run. **Fixed**: `run_3d` honors an optional
    `output_root` attribute; added an assertion that results land under `tmp_path`.
  - Informational addressed: `metadata.json` was promised in the plan but not written —
    **added** a metadata writer (grid, Re, pressure, source, radii, cylinders). Test now
    asserts its presence.
  - Informational left as-is (by design / out of scope): `build_3d_solver` ignores
    Lees-Edwards args (fouling/pressure scope); 3-tuple cylinders only (no z_lo/z_hi from
    config yet); the 2D runner's pre-existing ruff debt (264 lines, identical on main);
    VTK axis ordering verified correct by the reviewer.
- Judgment: No blocking issues — the reviewer confirmed the VTK ordering, .pvd format,
  dispatch placement, and all four acceptance scenarios are sound. The three Should-Fix
  items were all real: a redundant full-field recompute, a genuine crash path for a
  copied-from-2D config value, and test pollution of the repo tree. All fixed.

### Doc Parrot
- Divergences Found: 0
- Details:
  - `build_3d_solver` / `run_3d`: full Args blocks, updated to reflect the fixes
    (`periodic-force` normalization, `output_root`, `metadata.json`); prose matches code.
  - `_scalar_row`, `_write_pvd`, `_vtk_float`: intentional one-line docstrings for trivial
    private helpers; prose fully describes behaviour, params self-evident — judged OK,
    consistent with the codebase style.
  - VTK writers (`_write_fluid_vtk_3d`, `_write_particles_vtk_3d`): Args documented; prose
    matches the STRUCTURED_POINTS / POLYDATA implementation.
- Judgment: No prose/impl divergence; the public-API docstrings were kept in sync with the
  Round-1 fixes.
