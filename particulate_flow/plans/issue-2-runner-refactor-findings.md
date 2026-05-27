# Findings: issue-2-runner-refactor

## Round 1

### Code Reviewer
- Blocking: 2
- Should Fix: 4
- Informational: 3
- Key Issues:
  - Logic duplication: runner recomputed N_PARTICLES/SOURCE_VOLUME_FRACTION/GEOMETRY independently from builder → synced from sim after construction
  - TestBuildDEMPackingSolver supplied wrong arg keys (e_n, e_t, k_t…) → fixed with correct keys + forwarding assertions
  - Unused FLUID_METHODS import → replaced with meaningful FLUID_METHODS + FLUID_ACCELERATORS validation
  - _poiseuille_flow_rate duplicated between runner and builder → runner copy kept (needed for pre-sim logging), sim-based sync fixes divergence risk
  - Tests used cwd-relative Path("src/...") → fixed with _REPO_ROOT anchor
  - args.particle_source comparison used hyphenated form after normalisation → fixed to use normalised particle_source variable
- Judgment: Reviewer caught real issues — logic duplication was a genuine correctness risk; test schema mismatch would have silently passed with wrong defaults tested.

### Doc Parrot
- Divergences Found: 0
- Details: build_lbm_dem_solver and build_dem_packing_solver both OK — raises ValueError for invalid fluid_method is undocumented but acceptable for runner-targeted API; required args are implied by "parsed namespace after config defaults applied"
- Judgment: No real divergences found; parrot confirmed docstrings adequate for intended callers.
