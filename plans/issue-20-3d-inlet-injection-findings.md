# Findings: issue-20-3d-inlet-injection

## Round 1

### Code Reviewer
- Blocking: 2
- Should Fix: 4
- Informational: 5
- Key Issues:
  - **Blocking** — `_can_place_inlet` only checked particle-particle overlap, not
    cylinder bodies (the 2D analogue checks both). With #19 cylinders co-present at
    the inlet (the intended fouling use case), a sphere could be injected inside a
    cylinder → instant force spike. **Fixed**: added cylinder-overlap rejection
    (radial distance within the z-extent) and updated the docstring.
  - **Blocking** — `_feed_inlet_particles` called `_macroscopic()` twice per step
    (once via `_inlet_flow_rate`, once for `ux_inlet`). **Fixed**: compute the
    macroscopic field once and pass `ux_inlet` into `_inlet_flow_rate`.
  - **Should Fix** — `advance()` docstring omitted the inlet feed / outflow removal
    steps. **Fixed**: docstring now lists the full per-step sequence.
  - **Should Fix** — `left_inlet` with `source_volume_fraction=None` silently injected
    nothing. **Fixed**: construction now raises `ValueError` requiring a positive
    `source_volume_fraction`; added a test.
  - **Should Fix** — the `left_inlet` + non-`immersed_boundary` guard was untested.
    **Fixed**: added `test_left_inlet_requires_immersed_boundary_coupling`.
  - **Should Fix** — `remove_particles` did not validate `mask.size == n_p` (silent
    state corruption on misuse). **Fixed**: raises `ValueError` on size mismatch;
    added a test plus an add-then-remove round-trip test.
  - Informational: `_remove_outflow_particles` docstring said "x=nx-1" but the
    condition is `pos_x - r > nx` (fixed wording); per-step string checks for mode
    (negligible overhead); test uses 400 steps vs the 1500-step physics check; the
    `_le_shift` hardcodes `ny` (pre-existing, unrelated). Left as-is.
- Judgment: Strong, high-value review. Both blocking items were real: the cylinder
  overlap is a genuine fouling-scenario bug (injecting into a membrane pore), and the
  double `_macroscopic` was a real inefficiency. The Should-Fix guards/tests closed
  actual silent-failure and misuse paths. Core add/remove and budget physics verified
  correct.

### Doc Parrot
- Divergences Found: 0 (after Round 1 fixes)
- Details:
  - `add_particles`, `remove_particles`, `_inlet_flow_rate`, `_sample_inlet_point`,
    `_feed_inlet_particles`, `_can_place_inlet`, `_remove_outflow_particles`,
    `advance`, and the class `__init__` docstring all match their implementations
    after the fixes; every parameter is documented.
  - The two docstring corrections the reviewer flagged (`_can_place_inlet` omitting
    cylinders; `_remove_outflow_particles` "x=nx-1" vs `> nx`; `advance` missing
    inlet/removal steps) are now aligned.
- Judgment: The real divergences were exactly those the reviewer caught; once fixed,
  prose and code agree for every changed callable.
