# Findings: session-2026-05-28-simulate-skill

## Round 1

### Code Reviewer
- Blocking: 4
- Should Fix: 3
- Informational: 4
- Key Issues:
  - `result_tag` must not contain `/` or `\` — skill lacked this constraint
  - `four_cylinder_staggered.json` coordinates are fixed to 180×70 domain; no warning for other domain sizes
  - `adhesive_rolling_particles.json` material fragment not mentioned for attraction scenarios
  - Re stability threshold was vague (no number given)
  - No Bash timeout specified for runner invocation (default 120 s kills long runs)
  - Result directory path omitted `result_tag` suffix
  - `paraview_every` / `snapshot_storage` interaction not documented
  - Timestamp generation not instructed (agent may use wrong year)
- Judgment: Reviewer caught real and actionable issues — all Blocking items were fixed before proceeding.

### Doc Parrot
- Divergences Found: 2
- Details:
  - `_fractional_roll_x()` (lbm3d.py:228): scipy unavailable → silently falls back to linear interpolation even when `le_interpolation_order >= 3`. Not mentioned in docstring.
  - `_apply_le_bc()` (lbm3d.py:263): `le_shear_axis` parameter is stored but `cx = C3[i, 0]` is hardcoded — non-default `le_shear_axis` values are silently ignored. Docstring does not warn of this constraint.
- Judgment: Both findings are in pre-existing code not touched by this session's changes. The `_apply_le_bc` finding is a real latent bug worth filing as a separate issue. Not blocking this PR.
