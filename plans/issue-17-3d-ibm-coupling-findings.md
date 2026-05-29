# Findings: issue-17-3d-ibm-coupling

## Round 1

### Code Reviewer
- Blocking: 1
- Should Fix: 3
- Informational: 8
- Key Issues:
  - **Blocking** — stale velocity passed to `_collide_bgk` after IBM. `advance()`
    computed `_macroscopic()` once (with the *previous* step's F half-correction),
    then zeroed F and ran IBM, but still passed the old velocity to the collision —
    inconsistent with the Guo scheme and with the 2D path (which re-reads velocity
    after IBM). **Fixed**: compute the IBM exchange against the un-forced velocity,
    then re-fetch `_macroscopic()` so the collision sees the new force's
    half-correction. Verified the fix improved the tracer velocity ratio
    (0.972 → 0.983, closer to the ideal ~1.0).
  - **Should Fix** — dead `self.nu` attribute (`tau_to_nu()` stored but never read).
    **Fixed**: removed the assignment and the now-unused `tau_to_nu` method.
  - **Should Fix** — `test_reaction_balances_spread_force` ran the audit on a
    quiescent freshly-built solver, so both sides were trivially zero. **Fixed**:
    `advance(50)` first and assert the exchange is genuinely nonzero before checking
    the Newton balance.
  - **Should Fix** — no unit test for the Guo source-term moments. **Fixed**: added
    `TestGuoSourceMoments` asserting `Σ S_i = 0` (mass) and `Σ c_i S_i = (1-ω/2)F`
    (momentum) for a known state.
  - Informational (not actioned, by design / low value): per-step Fibonacci marker
    recompute (caching is a future perf win); `.any()`-scan early-exit in collide;
    `workflow_state.md` could be gitignored; `test_2d_ibm_untouched` is a light
    hasattr check; no rho divide-by-zero guard (matches 2D). Left as-is.
- Judgment: High-value review. The blocking item was a real Guo-consistency bug the
  integration tests masked (they only checked stability/finiteness, not scheme
  correctness); fixing it both corrected the physics and was confirmed by the new
  moment unit test. The trivial-zero test was a genuine false-confidence hole.

### Doc Parrot
- Divergences Found: 0
- Details:
  - `_collide_bgk`, `_macroscopic`, `_sphere_markers`, `_interp_velocity_3d`,
    `_spread_forces_3d`, `_apply_ibm`, `_ibm_force_audit`, `advance`, `DEM3D.step`:
    every changed callable's prose matches its implementation; all parameters
    documented (class docstring for `__init__`).
  - The `advance` docstring was updated alongside the blocking fix so it still
    describes the actual ordering (IBM → re-read velocity → collide → stream → BC →
    DEM step).
- Judgment: No divergences — the one risk (an `advance` docstring going stale after
  the ordering fix) was caught and the prose updated to match.
