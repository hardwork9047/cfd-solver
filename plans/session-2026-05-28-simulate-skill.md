> Note: This plan was generated retroactively from the completed changes.
> This is expected for interactive development sessions (/workon + /ship).

# Plan: /simulate skill — natural language membrane fouling simulation

## Goal

Add a Claude Code skill (`/simulate`) that lets users describe a membrane fouling simulation in plain Japanese (or English) and automatically generates a JSON config, shows it for approval, runs the simulation, and returns a results summary with file paths.

## Approach

- Implemented as a single `SKILL.md` file at `.claude/skills/simulate/SKILL.md`, consistent with all other skills in this repo.
- Claude Code itself (the current conversation) generates the config JSON — no separate API call required.
- The skill uses `configs/lbm_dem/templates/fouling_supply.json` as the base via `extends`, writing only the parameters the user mentions as a diff.
- Simulation execution uses `poetry run python src/runners/run_lbm_dem.py` with a 3600-second timeout.
- Results are surfaced as a text summary plus file path listing.

## Alternatives considered

- **Separate Claude API call for config generation**: rejected — the skill already runs inside a Claude Code conversation, making an extra API round-trip unnecessary.
- **Full automation (no confirmation)**: rejected — simulations are long-running and resource-intensive; user approval before execution is safer.
- **Template-based parameter extraction (regex/fixed form)**: rejected — LLM-based generation handles open-ended natural language better and requires no maintenance as parameters evolve.

## Scope

**Changed:**
- `.claude/skills/simulate/SKILL.md` — new skill file

**Explicitly left alone:**
- All Python source files, runners, configs, and tests — no code changes.

## Assumptions

- The `configs/lbm_dem/` directory structure (templates, geometries, materials) remains stable.
- `poetry run python src/runners/run_lbm_dem.py` is the correct entry point for all fouling simulations.
- Results land in `src/results/run_lbm_dem/<timestamp>_<result_tag>/` as documented in CLAUDE.md.
