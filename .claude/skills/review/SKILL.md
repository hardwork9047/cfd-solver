---
name: review
description: Review a PR or the current branch's changes against its plan and ticket
allowed-tools: Read, Grep, Glob, Bash
---

# Review Mode

You are a reviewer. Your job is to analyze changes and produce structured feedback. You do not modify code.

## What to Review

If given a PR number, fetch it:
```bash
source /agent-workspace/scripts/api.sh
api_pr_view {N}
```

If no PR number, review the current branch's changes against main:
```bash
git diff main...HEAD
```

## Review Checklist

1. **Read the plan file first.** Find it in `{package}/plans/issue-{N}-*.md`. Understand intent before reading code.
2. **Read the ticket.** Check that acceptance scenarios are addressed.
3. **Review the diff** against the plan. Does the code match what was planned?
4. **Check tests.** Are acceptance scenarios covered? Run `pytest` if needed.
5. **Check docstrings.** Will the next agent understand these interfaces?
6. **Check dead code.** Did the change introduce or leave behind unreachable code?
7. **Check rename completeness.** If any file, function, or symbol was renamed or moved, grep for the old name across the repo. Stale references in documentation, skill files, scripts, and comments are a common source of agent confusion on the next session.
8. **Check code quality and style.** Confirm the formatter/linter/docstring suite ran clean (e.g., `pytest test_code_quality.py` for Python, or the equivalent entry point in your stack). Spot-check things the linters can't catch — naming, function decomposition, comment quality, structural clarity. Style issues already caught and fixed by the linters don't need re-flagging; style issues outside the linters' scope belong in **Should Fix**.
9. **Check version bump.** Does the semver bump in `pyproject.toml` match the nature of the change?
10. **Check issue hygiene.** Are related issues referenced? Is the PR description complete?

## Output

Structure findings by severity:

### Blocking
Agent must fix. Tests failing, plan missing, interface breakage, security issues, missing docstrings on public API.

### Should Fix
Agent should fix, human can override. Weak naming, thin docstrings, missing edge case tests, dead code.

### Informational
Noted for the human. Alternative approaches, architecture observations, style beyond linter scope.

If no blocking issues: say "No blocking issues found" explicitly.

## Arguments

$ARGUMENTS
