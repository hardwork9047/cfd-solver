# Findings: issue-14-3d-fouling-core-solver

## Round 1

### Code Reviewer
- Blocking: 1
- Should Fix: 4
- Informational: 5
- Key Issues:
  - **Blocking** — `_zou_he_reconstruct` declared+documented two parameters
    (`known_neg_dirs`, `sign`) that the body never reads: misleading public API.
    Fixed by removing both; algorithm drives entirely from `OPPOSITE3[i]` and the
    sign already lives in the sign of `u_normal`.
  - Docstring of `_apply_pressure_bc` overclaimed "uy=uz=0 at both planes". Only
    `rho` and `ux` are enforced exactly. Fixed the docstring.
  - `n_particles` guard `if n_particles and n_particles > 0` let `-1` slip through.
    Simplified to `if n_particles > 0`.
  - Missing test for the `particle_fluid_coupling` stub. Added
    `test_coupling_raises_not_implemented`.
  - `streamwise_boundary.replace("-","_")` in builder was a no-op for valid values
    and mangled invalid ones before the constructor's ValueError. Removed.
  - Informational: redundant write-backs (`f[:,0,:,:]=plane`) on views — cleaned up;
    generous stability/monotonicity test thresholds — left as smoke bounds for the
    slice (noted for a future tightening pass).
- Judgment: Real and useful. The blocking item was a genuine API-honesty problem
  that would mislead the follow-up IBM/DEM agents; the physics itself was confirmed
  correct (direction groups, exact ux derivation, exact mass conservation).

### Doc Parrot
- Divergences Found: 0 (after Round 1 fixes)
- Details:
  - `__init__`: all params (incl. new pressure + stub params) documented; behaviour
    (guards, `rho_in` formula, mutual exclusion) matches prose.
  - `_apply_pressure_bc`: corrected docstring now matches numerics (rho/ux exact,
    uy/uz biased toward zero, not pinned — verified ~1e-18).
  - `_zou_he_reconstruct`: signature matches docstring after dead-param removal;
    `f_i = feq_i + (f_opp - feq_opp)` rule described accurately.
  - `advance`: periodic→LE vs pressure→Zou-He branching matches.
  - builder `_build_lbm3d_solver`: forwarding of pressure params documented.
- Judgment: The docstring overclaim caught by the reviewer was the only real
  divergence; once fixed, prose and implementation agree for every changed callable.
