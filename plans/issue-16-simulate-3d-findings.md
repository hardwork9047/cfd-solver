# Findings: issue-16-simulate-3d

## Round 1

### Code Reviewer
- Blocking: 0
- Should Fix: 4
- Informational: 6
- Key Issues:
  - **Should Fix** — step 5 (結果サマリー) still listed the 2D output tree
    (`summary.md`, `fields_*.vtk`, `snapshot_*.png`, root-level `*.csv`), which does not
    match what the 3D runner (`runner3d.py`) writes; an agent following step 5 for a 3D
    run would look for nonexistent files. **Fixed**: added a "3D の場合は出力ツリーが
    異なる" override describing the actual 3D tree (`paraview/`, `analysis/`,
    `metadata.json`, no png/mp4).
  - **Should Fix** — ambiguous `nz` wording ("set dimensions and nz" then "nz optional").
    **Fixed**: rephrased to "set `dimensions:3`; only set `nz` when the user gives a
    depth (else template default 24)".
  - **Should Fix** — `test_2d_default_behaviour_documented` checked `"fouling_supply"`,
    which `fouling_supply_3d` also satisfies (no-op guard). **Fixed**: now asserts the
    3D-specific strings `"3D 以外"` + `"デフォルト挙動は不変"`.
  - **Should Fix** — `test_warns_about_2d_only_assets` first assert (`"2D" in text`) was
    trivially true. **Fixed**: now asserts `"2D 専用"` (the actual warning phrase).
  - Informational left as-is: no 3D NL example in step 1 (step-2 detection rule is
    clear); `fields_*.vtk` is also inaccurate for the 2D runner (pre-existing, out of
    scope); version-bump judgment; `workflow_state.md` housekeeping.
- Judgment: No blocking issues — the reviewer confirmed the referenced assets exist, the
  3D example JSON dispatches correctly through the config pipeline, and the 2D path is
  intact. The step-5 mismatch was the most valuable catch (an agent would otherwise hunt
  for 2D-only files after a 3D run); the two weak test assertions were real no-op guards.
  All four fixed.

### Doc Parrot
- Divergences Found: 0 (N/A)
- Details: This PR changes no production Python callables — the only functional file is
  `.claude/skills/simulate/SKILL.md` (a markdown prompt interpreted at runtime). The new
  `tests/test_simulate_skill_3d.py` functions are documentation-presence guards with no
  docstrings to validate. There is no docstring↔implementation surface for the parrot to
  check.
- Judgment: Not applicable to a docs-only skill change; recorded for completeness.
