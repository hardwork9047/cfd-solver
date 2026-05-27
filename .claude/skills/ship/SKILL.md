---
name: ship
description: Finalize interactive work — generate retroactive plan, self-review, push, PR with annotations
allowed-tools: Read, Write, Edit, Bash, Glob, Grep, Agent
---

# Ship Mode

You are finalizing an interactive development session. The developer has been working alongside
you and is ready to file a PR. Your job is to run the quality gates and create a proper PR.

## Arguments

$ARGUMENTS

Parse the arguments:
- **Issue number (optional):** if provided (e.g., `/ship 42`), link the PR to that issue.
- **No arguments:** file a PR without issue linkage.

## Detect Context

Before starting, determine:
1. **Which package** was modified — use `git diff --name-only upstream/main...HEAD` to identify changed files. If files span multiple packages, ask the developer which is primary.
2. **Whether an issue number was provided** in the arguments.
3. **Current branch name** — used to derive the PR title if needed.
4. **Whether a workflow state file exists** — check `{package}/plans/workflow_state.md`. If it exists, read it to understand what was in progress. This is informational context for generating the plan; delete it when filing the PR.

## Workflow

### 1. Pre-flight Checks

```bash
source poetry_env/bin/activate
cd {package} && pytest
```

All tests must pass. If they don't, stop and tell the developer what failed.

```bash
pre-commit run --all-files
```

Fix any formatting/lint issues automatically. If fixes require judgment, ask the developer.

### 2. Generate Retroactive Plan

Check if a plan file already exists for the current branch:
```bash
ls {package}/plans/
```

If no plan file exists, generate one from the changes.

Analyze the diff:
```bash
git diff upstream/main...HEAD
git log upstream/main..HEAD --oneline
```

Create a plan file:
- **If issue number provided:** `{package}/plans/issue-{N}-slug.md`
- **If no issue number:** `{package}/plans/session-{YYYY-MM-DD}-slug.md`

The slug is derived from the branch name or a summary of the changes.

The plan must contain:
- **Goal:** what was accomplished (derived from the diff and session context)
- **Approach:** what was done and why
- **Alternatives considered:** if any were discussed during the session
- **Scope:** what was changed, what was explicitly left alone
- **Assumptions about other packages** (with docstring references if applicable)

Mark the plan as retroactive:
```
> Note: This plan was generated retroactively from the completed changes.
> This is expected for interactive development sessions (/workon + /ship).
```

Commit the plan:
```bash
git add {package}/plans/
git commit -m "Add retroactive plan for {slug}"
```

### 3. Version and Changelog

Check if the developer already bumped the version in `pyproject.toml`. If not:
- Determine the appropriate semver bump from the nature of the changes
- Bump the version
- Update `CHANGELOG.md`:
  - With issue#: `v{version}: {summary}, see PR #{N}`
  - Without: `v{version}: {summary}`

Commit if not already committed.

### 4. Self-Review (MANDATORY)

Run the code-reviewer agent:
```
Use the code-reviewer agent to review my changes
```

Wait for findings. Fix all **Blocking** items before proceeding.

For **Should Fix** items: the developer is present. Present each finding and ask whether to fix it or proceed. The developer's judgment overrides Should Fix items.

**Do not proceed to step 4b until the reviewer has run and all blocking items are resolved.**

### 4b. Doc-Parrot

After the code-reviewer, run the doc-parrot skill to validate docstring–implementation alignment for every changed callable in the diff:

```
/doc-parrot
```

The parrot identifies every changed callable, extracts docstring prose, and asks you to derive a description and usage example from prose alone — then compare that against the implementation. For each callable, record a judgment: **Fix docstring**, **Fix code**, or **OK**. Act on any **Fix** judgments before proceeding. Present findings to the developer for any judgment calls.

The parrot is not a gate. It produces an artifact you act on with judgment.

### 4c. Write (or update) findings file

After both the code-reviewer and the doc-parrot have run, append a round section to the
findings file. If the file doesn't exist yet, create it with the file header first.
Each review iteration adds a new round — do not overwrite previous rounds.

- **If issue number provided:** `{package}/plans/issue-{N}-slug-findings.md`
- **If no issue number:** `{package}/plans/session-{YYYY-MM-DD}-slug-findings.md`

```markdown
# Findings: {issue-N-slug or session-date-slug}

## Round 1

### Code Reviewer
- Blocking: N
- Should Fix: N
- Informational: N
- Key Issues: {bullet list or "none"}
- Judgment: {one-line assessment — did the reviewer catch something real?}

### Doc Parrot
- Divergences Found: N
- Details: {bullet list of callable + what diverged, or "none"}
- Judgment: {one-line assessment — did the parrot catch something real?}
```

If the developer asks you to fix blocking items and re-run, append `## Round 2` (and so on)
with the same structure. Commit after each round. Do not push until the latest round
shows `Blocking: 0`.

### 5. Push and PR

If `{package}/plans/workflow_state.md` exists, delete it before pushing:
```bash
git rm --ignore-unmatch {package}/plans/workflow_state.md
git diff --cached --quiet || git commit -m "chore: remove workflow state (PR filed)"
```

```bash
git push origin HEAD
```

Create the PR and capture the PR number:
```bash
PR_URL=$(api_pr_create "{concise title}" "{body}" "main")
PR_NUMBER=$(echo "$PR_URL" | grep -oP '\d+$')
COMMIT_SHA=$(git rev-parse HEAD)
```

**PR body when issue number is provided:**
```
Closes #{N}

{summary of changes}

Plan: {package}/plans/issue-{N}-slug.md
```

**PR body when no issue number:**
```
{summary of changes}

Plan: {package}/plans/session-{date}-slug.md
```

Reference any related issues mentioned during the session, even if not the primary issue.

### 6. Post PR Annotations

Using the `annotation-data` JSON block from the code-reviewer's output (step 4 Self-Review), post targeted
annotations to help the human reviewer.

**Guided tour comment** — post a summary comment on the PR:
```bash
api_pr_comment $PR_NUMBER "## Guided Review Tour
### Boundary Violations
{list from annotation-data, or 'None'}

### Architectural Decisions (worth a second opinion)
{list with file:line and descriptions}

### Areas of Uncertainty
{list with file:line and descriptions}

### Mechanical Changes (safe to skim)
{list of mechanical files}

---
_Generated by the implementing agent's self-review._"
```

**Line-specific annotations** — for each item in the annotation data:
```bash
# Boundary violations:
api_pr_review_comment $PR_NUMBER "$COMMIT_SHA" "{path}" {line} "⚠️ **BOUNDARY VIOLATION**: {reason}"

# Architectural decisions:
api_pr_review_comment $PR_NUMBER "$COMMIT_SHA" "{path}" {line} "🏗️ **ARCHITECTURAL DECISION**: {description}"

# Uncertainties:
api_pr_review_comment $PR_NUMBER "$COMMIT_SHA" "{path}" {line} "❓ **UNCERTAINTY**: {description}"
```

**Rules:**
- If no items to flag, post a simplified tour: "All changes are routine. No items flagged for special attention."
- If a line comment fails (line not in diff hunk), skip it silently.
- Limit to 15 line-specific comments maximum.
- Target `tri-projects/tri_amr_code` (upstream, where the PR lives).

### 7. Issue Hygiene (only if issue number provided)

```bash
api_issue_comment {N} "Addressed in PR #${PR_NUMBER}.

{brief summary of what was done}"
```

If no issue number was provided, skip this step entirely.
