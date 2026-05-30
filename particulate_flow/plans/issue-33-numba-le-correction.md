# Plan: issue-33 — Numba LBM path missing Lees-Edwards correction

## Goal (from issue)

When `fluid_accelerator="numba"` (or auto with numba available), the
Lees-Edwards boundary correction (x-shift + velocity-jump boost) must be
applied identically to the NumPy path, so both accelerators produce the same
linear shear profile `u_x = γ̇·y` with L2 error < 1%.

## Root Cause

`LBMDEMSolver._lbm_step` (`lbm_dem.py:1057`) has an early-return branch for
the numba path that calls `_lbm_step_numba` and then `return`s without ever
calling `_apply_le_streaming_correction`. The numpy path calls it on line 1091.

## Approach

**One-line fix:** remove the early `return` after the numba path and fall
through to the same post-streaming logic as the numpy path. The LE correction
is applied to `self.f` in-place, so there is no ordering conflict.

Concretely, `_lbm_step` becomes:

```python
def _lbm_step(self, rho, ux, uy):
    if self.uses_numba_lbm:
        method_id = 1 if self.fluid_method == "lbm-trt-guo" else 0
        self.f = _lbm_step_numba(...)
        # Fall through to shared post-streaming block
    else:
        # existing numpy collision + streaming code

    # Shared post-streaming for both paths:
    if self.y_boundary == "lees_edwards":
        self._apply_le_streaming_correction()
    # Bounce-back (numba kernel already does bounce-back internally, so we
    # must NOT duplicate it for the numba path — only LE correction is shared)
    if not self.uses_numba_lbm:
        f_tmp = self.f.copy()
        for i in range(Q):
            self.f[i, self.solid] = f_tmp[OPPOSITE[i], self.solid]
    self._apply_pressure_boundaries()
```

Wait — bounce-back is handled *inside* the numba kernel; pressure boundaries
are already called inside the numba branch. Therefore the refactor is:

1. Move `_apply_pressure_boundaries()` out of the numba branch.
2. Add `_apply_le_streaming_correction()` after the numba kernel (when
   `y_boundary == "lees_edwards"`).
3. The numpy branch stays unchanged (it already calls LE + bounce-back +
   pressure in order).

Simplest correct restructure:

```python
if self.uses_numba_lbm:
    self.f = _lbm_step_numba(...)
    if self.y_boundary == "lees_edwards":
        self._apply_le_streaming_correction()
    self._apply_pressure_boundaries()
    return
# numpy path unchanged below...
```

This is minimal: no restructuring of the numpy path, no new helpers, just
two lines added before the existing `return` in the numba branch.

## Alternatives Rejected

- Implementing LE correction inside the numba JIT kernel: requires passing
  many additional arguments (W, C, LE_TOP_CROSS_DIRS, etc.) into the kernel
  and complicates the hot loop. Overkill for a correctness fix.
- Merging all post-streaming logic: requires careful handling of the fact that
  bounce-back is already inside the numba kernel but external in numpy. The
  minimal patch is cleaner.

## Out of Scope

- 3D Lees-Edwards (handled separately).
- Numba kernel performance improvements.

## Assumptions

- `_apply_le_streaming_correction` reads/writes `self.f` and uses `self._le_shift`,
  `self.le_shear_rate`, `self.ny` — all of which are set before `_lbm_step` is
  called. No threading issues.
- `_apply_pressure_boundaries` is idempotent if called once.

## Test Plan

1. Add `test_le_numba_linear_shear_profile` in `tests/lbm/test_lees_edwards.py`
   — mirrors the existing `test_le_linear_shear_profile` but passes
   `fluid_accelerator="numba"`, marked `@pytest.mark.slow` and skipped if
   numba is not installed.
2. Add `test_le_numpy_numba_field_agreement` — runs same sim with both
   accelerators for 500 steps from the same seed, asserts steady-state `ux`
   fields agree within 1e-4 (verifying physical symmetry).
3. Both tests are gated with `pytest.importorskip("numba")` so they skip
   cleanly in environments without numba.
