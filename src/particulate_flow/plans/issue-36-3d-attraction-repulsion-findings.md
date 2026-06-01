# Findings: issue-36-3d-attraction-repulsion

## Round 1

### Code Reviewer
- Blocking: 1
- Should Fix: 2
- Informational: 4
- Key Issues:
  - **Blocking**: `LBMDEMSolver3D` docstring missing 10 new Args entries (`sliding_friction`, `rolling_friction_coeff`, and 8 Hamaker params). Fixed: added all entries to the Args section.
  - **Should Fix 1**: `normal_i` was defined twice inside `_sphere_sphere_loads` (once in the contact block, once unconditionally after). Fixed: hoisted to a single definition before the contact block.
  - **Should Fix 2**: No test for repulsion zero beyond cutoff (asymmetry with attraction test). Fixed: added `test_repulsion_zero_beyond_cutoff`.
- Judgment: Reviewer caught real issues. The Blocking item (missing docstring) was the plan's own stated requirement that was missed. Both Should Fix items were genuine gaps.

### Doc Parrot
- Divergences Found: 0
- Details:
  - `DEM3D` class: OK — new Hamaker params documented correctly in Args section; mutual exclusivity ValueError is communicated.
  - `_sphere_sphere_loads()`: OK — "even without contact" comment is accurate (Hamaker fires outside contact range); Hamaker also fires during contact (DLVO-style superposition), which is intentional and mirrors 2D implementation. Not surprising enough to flag.
  - `_cylinder_loads()`: OK — same reasoning as above.
- Judgment: No doc-code divergences. All changed callables are correctly documented.
