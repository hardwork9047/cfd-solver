# Findings: issue-26-dem-decouple

## Round 1

### Code Reviewer
- Blocking: 0
- Should Fix: 3
- Informational: 4
- Key Issues:
  - **Should Fix** — `ContactLaw.from_solver` read `getattr(sim, "contact_model",
    "dem-hertz")`, but `LBMDEMSolver` stores the model as `particle_method`, so
    `from_solver` silently returned `"dem-hertz"` for any `dem-linear` solver. Not a
    production bug (DEMSolver uses the direct constructor), but `from_solver` is public.
    **Fixed**: fall back to `particle_method` (`contact_model or particle_method or
    'dem-hertz'`); verified it now reads `dem-linear` correctly.
  - **Should Fix** — `from_solver` docstring said "contact_model-less material
    attributes" then described reading `contact_model` — contradictory. **Fixed**:
    rewrote to describe the `contact_model`/`particle_method` fallback and added a Note
    that DEMSolver builds the law directly (not via `from_solver`).
  - **Should Fix** — the equivalence test only covered `dem-hertz`; the `dem-linear`
    delegation through DEMSolver was untested. **Fixed**: parametrized
    `test_normal_matches` over both models (and pass `particle_method` to the solver).
  - Informational left as-is: three unused private contact methods on `LBMDEMSolver`
    (`_normal_contact_magnitude` etc.) predate this PR — noted for a future cleanup;
    `TestThreeDConceptParity` only checks normal (tangential/rolling differ 2D↔3D by
    design); `workflow_state.md` housekeeping.
- Judgment: No blocking issues — the reviewer verified the three formulas are
  byte-identical to the original, the `contact_model` routing is correct (DEMSolver's own
  arg, not the coupled solver's), and param caching is safe (params set once before
  DEMSolver is built). The behaviour-preserving claim holds. The `from_solver` linear-model
  fallback was a real latent bug in the new public API; fixed before any caller adopts it.

### Doc Parrot
- Divergences Found: 0 (after Round 1 fixes)
- Details:
  - `ContactLaw.normal_magnitude/tangential_magnitude/rolling_torque`: prose (clamp/
    Coulomb laws, zero-guards) matches the implementation, which is verbatim-equal to the
    original DEMSolver formulas.
  - `ContactLaw.from_solver`: the corrected docstring matches the new
    contact_model/particle_method fallback and the DEMSolver note.
  - `DEMSolver` delegating methods: docstrings state "delegates to ContactLaw; behaviour
    unchanged" — accurate.
- Judgment: The one real divergence (from_solver's contradictory wording) was exactly the
  reviewer's Should-Fix and is now aligned; all changed callables match.
