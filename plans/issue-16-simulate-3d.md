# Plan: Issue #16 — /simulate Skill 3D Support

## Verbatim Goal (from issue #16)

> `/simulate` スキルが3Dパラメータ（`dimensions: 3`、`nz`）を自然言語から解釈し、
> 3D膜ファウリングシミュレーションを2Dと同じ操作感で実行できるようにする。

### Acceptance Scenarios
1. 「3D、nz=16、Re=10で走らせて」のような入力から3D用コンフィグが生成される
2. 確認テーブルに `nz`・`dimensions` が表示される
3. 3D未対応のテンプレートを指定した場合に警告が表示される
4. 2Dシミュレーションの挙動が変わらない

### Constraints
- `.claude/skills/simulate/SKILL.md` の変更のみ（**Python コード変更なし**）
- #15 完了後 — **merged**
- `fouling_supply_3d.json` テンプレートが存在することを前提 — **存在する** (#15 で追加)

## Scope decision

Pure prompt/skill edit: teach `/simulate` to recognise 3D intent and emit a config
that extends `fouling_supply_3d.json` with `domain.dimensions=3` + `nz`.  The runner
side already works (`run_lbm_dem.py --config <3d>` dispatches to `runner3d.run_3d`,
verified in #15).  No Python changes.

### In scope (this PR)
- **SKILL.md step 2 (config generation)**: add a 3D rule — when the user mentions 3D
  / `nz` / a depth, base the config on `../templates/fouling_supply_3d.json` (instead
  of the 2D `fouling_supply.json`), set `domain.dimensions: 3` and `domain.nz`, and use
  the 3D cylinder geometry (`../geometries/four_cylinder_3d.json`) when cylinders are
  requested.  Default to 2D when no 3D intent is expressed (scenario 4).
- **A 3D config example** block mirroring the existing 2D example.
- **SKILL.md step 3 (confirmation table)**: show `dimensions` and `nz` rows when 3D.
- **A warning rule (scenario 3)**: if the user pins a 2D-only template/geometry while
  asking for 3D (e.g. the 2D `four_cylinder_staggered.json`, which is sized for
  nx=180,ny=70 and has no z-extent), warn that it is 2D-only and use the 3D template/
  geometry instead.
- **Tests**: a docs-presence guard (`tests/test_simulate_skill_3d.py`) asserting the
  skill file documents the 3D template path, `dimensions`/`nz`, and the warning — so a
  future edit that drops 3D support is caught.

### Out of scope
- Any Python change (runner/solver/config) — all already done in #15 and the solver
  stack.
- The `/grill-me` confirmation-table redesign sitting in `stash@{0}` — unrelated to
  #16; left in its stash, not mixed into this PR.
- `/sweep`, `/analyze` skills' 3D support (separate skills, not in this issue).

## Approach

Edit `.claude/skills/simulate/SKILL.md`:
1. Step 2 rules: add a "3D の場合" bullet group — base template
   `fouling_supply_3d.json`, set `domain.dimensions:3` + `nz` (default 24 if the user
   says "3D" without a depth, matching the template), 3D cylinder geometry fragment.
2. Add a 3D generation example (Re=10, nz=16) next to the 2D one.
3. Step 2 param reference: note `domain.nz` / `domain.dimensions` and the 3D template.
4. Step 3 confirmation table: include `dimensions` / `nz` rows for 3D runs.
5. Notes: a warning rule for 2D-only template/geometry misuse under 3D intent.

## Alternatives considered & rejected
- **Auto-detect 3D from `nz` alone**: kept, but also trigger on explicit "3D"/"三次元"
  wording so `nz` is not required (scenario 1 says "3D、nz=16" — both signals).
- **Make 3D the default**: rejected — scenario 4 requires 2D behaviour unchanged when
  no 3D intent is expressed.
- **Incorporate the stashed `/grill-me` table**: rejected — out of scope; would
  conflate two unrelated changes in one PR.

## Assumptions about other modules (docstring-grounded)
- `fouling_supply_3d.json` exists with `domain.dimensions=3`, `nz`, pressure flow,
  `left-inlet`, `immersed_boundary` (added in #15; verified in repo).
- `four_cylinder_3d.json` geometry fragment exists (added in #15) for 3D obstacles.
- `run_lbm_dem.py --config <3d>` dispatches to the 3D path on `dimensions==3`
  (`run_lbm_dem.py` dispatch + `runner3d.run_3d`, #15) — the skill's step-4 run command
  is unchanged (same `poetry run python src/runners/run_lbm_dem.py --config ...`).

## Testing plan (acceptance → test) — `tests/test_simulate_skill_3d.py`
- **Scenario 1/2** (3D config + table): assert SKILL.md references
  `fouling_supply_3d.json`, mentions `dimensions` and `nz`, and contains a 3D example
  with `"dimensions": 3`.
- **Scenario 3** (warning): assert SKILL.md contains a warning about 2D-only
  templates/geometries under 3D intent.
- **Scenario 4** (2D unchanged): assert the 2D base template `fouling_supply.json`
  and the 2D example are still present (the 3D additions don't remove the 2D path).
- These are documentation-presence checks — the skill is interpreted by the agent at
  runtime, so the test guards that the instructions remain present, not runtime
  behaviour.
