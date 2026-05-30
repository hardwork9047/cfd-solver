# Findings: issue-27-analytical-benchmarks

## Round 1

### Code Reviewer
- Blocking: 0
- Should Fix: 2
- Informational: 4
- Key Issues:
  - `test_poiseuille_2d_umax_matches_body_force` used 15% tolerance — twice the documented ~7% offset; tightened to 12% (covers documented offset + noise, rejects 8–14% regressions)
  - `test_poiseuille_2d_grid_convergence` assertion `u_max_errors[2] <= u_max_errors[0] * 1.1` allowed 10% regression; replaced with strict monotone-decrease assertions at each refinement step (16→24 and 24→32), matching measured data (14.4% → 9.3% → 6.8%)
  - Docstring updated to accurately describe the assertion as monotone-decrease (not "2nd-order convergence")
- Judgment: Reviewer caught real correctness gaps in assertion semantics; both applied.

### Doc Parrot
- Divergences Found: 0
- Details: `_poiseuille_2d` and `_parabola_fit_l2` are private test helpers with minimal but accurate docstrings; test functions skipped per parrot rules
- Judgment: No action required.
