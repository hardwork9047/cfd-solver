# Findings: issue-33-numba-le-correction

## Round 1

### Code Reviewer
- Blocking: 0
- Should Fix: 2
- Informational: 4
- Key Issues:
  - `test_le_numba_linear_shear_profile` and `test_le_numpy_numba_field_agreement` did not assert `sim.uses_numba_lbm is True` — if the accelerator silently fell back to numpy both tests would validate the numpy path twice
  - Tolerance of `1e-4` in `test_le_numpy_numba_field_agreement` may be flaky due to floating-point accumulation order differences between the two paths; loosened to `1e-3`
- Judgment: Reviewer caught real correctness gaps in the test guards; both applied.

### Doc Parrot
- Divergences Found: 0
- Details: none — `_le_sim_numba` has a minimal but accurate docstring adequate for a private test helper; `test_*` functions skipped per parrot rules
- Judgment: No action required.
