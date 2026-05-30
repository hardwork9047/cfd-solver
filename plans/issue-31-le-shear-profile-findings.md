# Findings: issue-31-le-shear-profile

## Round 1

### Code Reviewer
- Blocking: 0
- Should Fix: 3
- Informational: 6
- Key Issues:
  - **Should Fix** — `init_analytical` was not in the `LBMDEMSolver` class docstring
    (public param other agents must discover). **Fixed**: added a full entry mirroring
    the 3D solver, including the `le_shear_rate==0` no-effect note.
  - **Should Fix** — `tests/conftest.py` comment still said `--import-mode=importlib`
    is "set in pyproject.toml" after the migration. **Fixed**: updated to `pytest.ini`.
  - **Should Fix** — `init_analytical=True` with `le_shear_rate==0` is silently ignored
    (the `and le_shear_rate != 0.0` guard), diverging from 3D. **Fixed by documenting**
    the no-effect case in the docstring (chose documentation over raising, since the
    profile would be all-zeros anyway and raising would be a surprising behaviour change).
  - Informational acted on: CHANGELOG date was 2026-05-30 → corrected to 2026-05-29.
  - Informational surfaced (not fixed, out of scope): the **Numba LBM path
    (`_lbm_step_numba`) returns early and never calls `_apply_le_streaming_correction`**,
    so the LE velocity boost is silently absent when `fluid_accelerator="numba"`. This
    is a pre-existing bug that this PR makes more impactful (NumPy now correct, Numba
    still wrong). The test pins `fluid_accelerator="numpy"`. Worth a follow-up issue.
- Judgment: No blocking issues — the reviewer independently confirmed the velocity-boost
  physics (`Δf_i = w_i·ρ·c_ix·(±dv)/cs²`, top −dv / bottom +dv, mass-conserving) exactly
  mirrors the proven 3D `_apply_le_bc`, and that `init_analytical` is bit-identical to 3D
  and truly opt-in. The diagnosis (three layered causes + the pytest.ini header) was
  validated. High-value review; the Numba-skip catch is the most important follow-up.

### Doc Parrot
- Divergences Found: 0
- Details:
  - `_apply_le_streaming_correction`: docstring describes both the x-shift and the new
    velocity boost (formula + sign convention) — matches the implementation.
  - `LBMDEMSolver.__init__` `init_analytical`: class docstring matches the
    `init_analytical and le_shear_rate != 0` guard, including the no-effect case.
- Judgment: prose and implementation agree for every changed callable.

## Follow-up worth filing
- **Numba LE correction gap**: `_lbm_step_numba` skips `_apply_le_streaming_correction`,
  so `fluid_accelerator="numba"` silently omits the LE velocity boost (and the x-shift).
  Pre-existing, but now a NumPy/Numba physics asymmetry for Lees-Edwards runs.
