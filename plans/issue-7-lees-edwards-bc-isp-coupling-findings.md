# Findings: issue-7-lees-edwards-bc-isp-coupling

## Round 1

### Code Reviewer
- Blocking: 3
- Should Fix: 4
- Informational: 3
- Key Issues:
  - `_le_shift` updated after `_lbm_step` but before DEM — timing inconsistency (fluid and particles saw different shifts)
  - `_fractional_roll_x` cubic path used `np.interp` (linear) not actual cubic interpolation
  - iSP coupling was identical to point_force with no behavioural difference
  - Particle wrap unit tests were structurally broken (called `_dem_substep` directly without ever setting `_le_shift`)
- Judgment: Reviewer caught real issues — all three blocking items were genuine bugs fixed in this session.

### Doc Parrot
- Divergences Found: 2
- Details:
  - `_expand_solver_section` — docstring omitted the boundary-injection side-effect when `enabled: true`
  - `_interp_velocity_cubic` — "periodic" in docstring was misleading; also missing scipy-fallback note
- Judgment: Parrot caught real gaps; both docstrings updated.
