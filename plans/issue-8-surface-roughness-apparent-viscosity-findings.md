# Findings: issue-8-surface-roughness-apparent-viscosity

## Round 1

### Code Reviewer
- Blocking: 3
- Should Fix: 4
- Informational: 3
- Key Issues:
  - η^C always zero — collisional stress accumulation not wired (fixed: added `compute_contact_stress_xy()` to DEMSolver)
  - Solid node bias in `_fluid_stress_xy` — `np.mean` over all nodes including bounce-back solid cells (fixed: masked to `~sim.solid`)
  - `record()` per-frame mismatch — modulo check skipped flush points when called coarsely (fixed: replaced with `_next_flush_step` threshold logic)
  - Summary.json integration untested; `workflow_state.md` committed to branch
- Judgment: Reviewer caught three real correctness defects; all were fixed before Round 2.

### Doc Parrot
- Divergences Found: 3
- Details:
  - `_fluid_stress_xy` Args section incorrectly listed `rho`, `ux`, `uy` as required attributes (they are computed internally from `f`) → fixed
  - `record()` docstring did not mention internal call to `sim.dem_solver.compute_contact_stress_xy()` → fixed with Note
  - `compute_contact_stress_xy()` formula description said "contact force" without clarifying normal-only (tangential excluded) → fixed
- Judgment: Parrot caught one wrong docstring and two real gaps that would mislead integrators; all fixed.

## Round 2 (post-fix verification)

### Code Reviewer
- All three blocking items resolved.
- Blocking: 0
- Should Fix: 0 (remaining items accepted: `workflow_state.md` gitignore is out-of-scope for this PR; summary.json test deferred)
- Informational: workflow_state.md in diff (cosmetic)

### Doc Parrot
- Divergences Found: 0
- All three docstring fixes applied and verified.
