# Findings: issue-18-3d-dem-contact

## Round 1

### Code Reviewer
- Blocking: 1
- Should Fix: 3
- Informational: 5
- Key Issues:
  - **Blocking** — `attach_dem`/`step_dem` hook on `LBMDEMSolver3D` was listed
    in-scope in the plan but not delivered. **Resolved by documenting a scope trim**:
    deferred to #17 because a fluid-free DEM hook has no real use (the only reason to
    embed a DEM stepper in the fluid solver is to exchange drag + IBM back-reaction,
    which is #17). `DEM3D` is fully standalone and unit-tested; adding the hook now
    would be dead surface area #17 would rework. Plan updated with rationale.
  - **Should Fix** — `DEM3D` not exported from `dem/__init__.py`. Fixed: added import
    and `__all__` entry.
  - **Should Fix** — Scenario-3 test asserted a rolling-resistance torque but `omega=0`
    at call time, so `_rolling_torque` returned zeros; the assertion passed via the
    tangential-friction torque instead. Fixed: corrected the comment and added
    `test_rolling_resistance_opposes_spin` (nonzero omega, asserts z-torque opposes
    spin) and `test_rolling_resistance_vanishes_without_contact`.
  - **Should Fix** — `_tangential_magnitude` docstring said `v_t` is "signed" but the
    caller always passes a non-negative norm. Fixed the docstring.
  - Informational: arithmetic-mean effective mass (matches 2D); wall planes at bare
    integer coords vs 2D's half-cell offset (noted for #17); no `density_ratio>0` /
    `dem_substeps>0` validation (matches 2D, low risk); rolling uses absolute omega
    (matches 2D, plan wording softened). Left as-is with notes for #17.
- Judgment: Real and valuable. The blocking item was a genuine plan/delivery mismatch;
  resolving it as a documented scope trim is the honest call. The scenario-3 test gap
  was a true coverage hole the reviewer correctly caught — now closed with a test that
  actually exercises `_rolling_torque`. Physics verified correct across all paths.

### Doc Parrot
- Divergences Found: 0 (after Round 1 fixes)
- Details:
  - `__init__`: all params documented in the class docstring (standard convention).
  - `_normal_magnitude`, `_rolling_torque`, `_add_contact`, `compute_loads`, `step`:
    prose matches implementation.
  - `_tangential_magnitude`: the one real divergence ("signed" vs always-non-negative)
    was fixed; now matches.
- Judgment: The "signed" wording was a genuine (if minor) divergence that would mislead
  the next agent; once fixed, prose and code agree for every changed callable.
