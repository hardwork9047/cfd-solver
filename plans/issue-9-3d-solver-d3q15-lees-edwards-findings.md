# Findings: issue-9-3d-solver-d3q15-lees-edwards

## Round 1

### Code Reviewer
- Blocking: 1
- Should Fix: 3
- Informational: 5
- Key Issues:
  - Linear fallback in `_fractional_roll_x` broken for negative shifts: `int()` truncates toward zero giving negative `frac`, causing extrapolation instead of interpolation. Fixed with `math.floor`. (Blocking)
  - `_le_shift` updated before collide/stream, causing shift to lead by one step. Fixed by moving increment after streaming and before `_apply_le_bc`. (Should Fix)
  - No test for linear fallback path when scipy absent. (Should Fix — deferred; scipy is a project dependency)
  - `_fractional_roll_x` docstring says "positive = shift right" but implementation shifts left. (Should Fix → fixed as docstring)
- Judgment: Reviewer caught one real correctness defect (math.floor fix) and one step-ordering bug; both fixed. D3Q15 constants and LE BC physics verified correct.

### Doc Parrot
- Divergences Found: 2
- Details:
  - `_apply_le_bc` docstring said "Top-crossers receive +dv" but code applies -dv (correct sign physically; docstring was wrong) → fixed
  - `_fractional_roll_x` docstring said "positive = shift right" but np.roll(-shift) shifts left → fixed
- Judgment: Both were real wrong-sign errors in docstrings. Physics confirmed correct by shear profile test (L2≈0).

## Round 2 (post-fix verification)

### Code Reviewer
- Blocking: 0
- All blocking and should-fix items resolved.

### Doc Parrot
- Divergences Found: 0
- Both docstring corrections applied.
