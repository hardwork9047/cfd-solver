# Findings: issue-19-3d-cylinder-obstacles

## Round 1

### Code Reviewer
- Blocking: 0
- Should Fix: 3
- Informational: 7
- Key Issues:
  - **Should Fix** — `_ibm_force_audit` docstring claimed Newton's 3rd law holds
    unconditionally, but the solid-skip in `_spread_forces_3d` breaks it when a
    particle's markers overlap a cylinder. **Fixed**: qualified the docstring (balance
    holds without obstacles; the obstacle absorbs the dropped momentum otherwise).
  - **Should Fix** — no test covered the IBM audit with cylinders present. **Fixed**:
    added `TestIBMSolidInteraction` — a particle touching a cylinder shows a bounded
    spread/reaction imbalance, while the no-cylinder control still balances to ~0.
  - **Should Fix** — latent defect: a cylinder touching the pressure inlet/outlet
    plane would have its bounce-back silently overwritten by the Zou-He BC. **Fixed**:
    `_build_cylinder_solid` now raises `ValueError` in pressure mode if a cylinder's
    x-extent reaches x=0 or x=nx-1; added a test for the guard.
  - Informational (not actioned): `np.where` solid-skip is less efficient than 2D
    boolean indexing (Numba deferred); `has_solid`/`solid.any()` could be cached at
    construction; single-point centreline test; generous particle-penetration
    tolerance; z_hi-inclusive convention (now documented). Left as-is.
- Judgment: Strong review. No blocking issues — the core bounce-back, ordering, and
  OPPOSITE3 were all verified correct. The three Should-Fix items were real: the
  docstring/test gap around the *intentional* momentum imbalance near obstacles was
  the most valuable catch (a future agent would otherwise read it as a conservation
  bug), and the inlet-plane guard turns a silent latent defect into a loud error.

### Doc Parrot
- Divergences Found: 0
- Details:
  - `_build_cylinder_solid`: prose (radial + z-extent mask, tuple conventions, new
    `Raises` for pressure-plane overlap) matches the implementation.
  - `_apply_bounce_back`: halfway BB + no-op guard match.
  - `_spread_forces_3d`: updated prose about skipping solid nodes matches the
    `np.where(solid, 0, contrib)` guard.
  - `_ibm_force_audit`: the corrected docstring now matches the obstacle-aware
    behaviour.
- Judgment: No divergences after the Should-Fix docstring correction; the one real
  risk (the audit docstring) was exactly what the reviewer flagged and is now aligned.
