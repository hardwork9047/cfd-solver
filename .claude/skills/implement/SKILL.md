---
name: implement
description: Start an implementation task for a specific package and GitHub issue
allowed-tools: Read, Write, Edit, Bash, Glob, Grep, Agent
---

# Implementation Mode

You are an implementing agent. Your job is to pick up a GitHub issue and deliver a clean PR.

## Setup

Identify your target from the arguments:
- **Package:** which package are you working in?
- **Issue:** what GitHub issue number?

If either is missing, ask the user before proceeding.

### Resume Check

Before asking the plan-review question, check whether a prior session left a workflow state:

```bash
cat {package}/plans/workflow_state.md 2>/dev/null
git branch --show-current
git log --oneline -10
```

If the file exists and `issue:` matches your target issue, **validate it before trusting it**:

1. **Branch check** — does the current branch name match `feature/issue-{N}-*` or `bugfix/issue-{N}-*`? If you're on a branch for a *different* issue, the state file is stale — ask the user before proceeding.
2. **Step plausibility** — does the recorded `step:` make sense given recent commit messages? If the state says `step: 5 (write tests)` but git log shows commits that look like cleanup or review work, the file is behind reality — tell the user what you see and ask how to proceed.
3. **Package match** — does the `package:` field match where you're working? A state file from a different package in the same repo is always stale.

If any check fails, surface what you found:
> "I found a workflow_state.md for issue #{N} at step X, but [what's inconsistent]. Should I resume from that step, or start fresh?"

If all checks pass, **resume from the recorded step** — do not restart from step 1. Honor the recorded `plan-review-preference` and do not ask the question again.

### Plan Review Preference (skip if resuming)

**Ask this BEFORE starting any work.** The user may step away during the long tail of this skill, so you must set expectations upfront. Ask:

> "Before I start: when I finish writing the plan (step 4), do you want to review it before I commit and start coding, or should I proceed autonomously?"

Remember the answer and honour it when you reach step 4:
- **review**: pause at step 4, show where the plan file is, wait for their feedback, apply any edits they request, then commit.
- **proceed**: commit the plan at step 4 and keep going without prompting.

Do not ask this question later — the point is to lock it in before the user walks away.

### Workflow State Format

At each step transition, write `{package}/plans/workflow_state.md` and include it in the step's commit. This file survives context compaction and lets a resumed agent pick up exactly where the previous one stopped.

```
skill: implement
package: {package}
issue: {N}
plan-review-preference: {proceed|review}
step: {N} ({step name})
next: {one sentence — what to do first on resume}
```

Delete this file in the same commit as the PR is filed (step 9).

## Workflow

Follow this sequence exactly:

### 0. Write Initial State

Before doing anything else, write the workflow state so a resumed agent knows where you are:

```bash
mkdir -p {package}/plans
cat > {package}/plans/workflow_state.md << 'EOF'
skill: implement
package: {package}
issue: {N}
plan-review-preference: (pending — not yet asked)
step: 1 (sync and branch)
next: sync, branch, set up environment, read the issue
EOF
```

Update this file at every subsequent step transition.

### 1. Sync and Branch
```bash
git fetch upstream && git rebase upstream/main
git checkout -b {prefix}/issue-{N}-short-slug
```
Create a new branch from the latest upstream main. Name it after the issue.

**Branch prefix:** read the issue first (step 3), then choose:
- `feature/` — new functionality, enhancements
- `bugfix/` — fixing broken behaviour

### 2. Set Up Environment

All packages share a single virtualenv (`poetry_env`) at the repo root. Check if it exists and is activated:

```bash
# If poetry_env doesn't exist yet:
python3 -m venv poetry_env
source poetry_env/bin/activate
pip install --upgrade pip
pip install poetry

# Install in dependency order using poetry (matches CI exactly):
for pkg in wcs_tri_api_models tri_arrival_estimate_models tri_amr_common amr_optimal_routing rms_replay_tri tri_amr_sim tri_planner_api_server; do
  (cd $pkg && poetry install) || exit 1
done
pip install mypy types-PyYAML
```

If `poetry_env` already exists, just activate it:
```bash
source poetry_env/bin/activate
```

Also load the GitHub API helpers — available for all subsequent steps:
```bash
source /agent-workspace/scripts/api.sh
```

Verify by running tests on your target package:
```bash
cd {package}
pytest
```

Ensure existing tests pass **before** making any changes.

### 3. Understand the Issue
```bash
api_issue_view {N}
api_issue_list | jq '.[] | {number, title}'  # scan for related issues
```
Read it thoroughly. Note any acceptance scenarios.

### 4. Plan
Create `{package}/plans/issue-{N}-short-slug.md` with:
- Verbatim goal from the issue
- Your chosen approach and decomposition
- Alternatives you considered and rejected
- What's explicitly out of scope
- Assumptions about other packages (with docstring references)
- How you'll test each acceptance scenario

**Honour the plan-review preference captured in Setup.**
- If the user asked to review: tell them the plan is ready and the file path, and wait. Apply any revisions they request to the file, then commit.
- If the user said proceed: commit the plan and continue.

**Commit the plan before writing any code.** Include a `workflow_state.md` update in the plan commit:

```
step: 5 (write tests)
next: translate acceptance scenarios into failing tests, then commit
```

### 5. Write Tests First
Translate acceptance scenarios from the issue into failing tests. Commit them. Include a `workflow_state.md` update in the same commit:

```
step: 6 (implement)
next: write code until tests pass; run fast quality checks iteratively
```

### 6. Implement
Write code until tests pass. During the inner loop, run your new tests and the fast quality tier as two quick commands:
```bash
pytest {your-new-test-files}
pytest -m quality_fast tests/test_code_quality.py
```
Your tests catch regressions; `quality_fast` catches formatting and lint before they accumulate. The `-m` filter is global, so combining them into one command would silently drop your unmarked tests — keep them separate. Run plain `pytest` (all unit + integration + quality tests) before step 7.

**Never suppress quality check failures by excluding files.** If a quality check fails because of missing modules, broken imports, or pre-existing issues, fix the root cause. Adding files to exclude lists (like `IMPORT_EXCLUDE`) to make tests pass is not fixing — it's hiding. If the root cause is genuinely outside your scope, escalate to the human rather than masking the failure.

When all tests pass, update `workflow_state.md` and commit before moving to step 7:

```
step: 7 (clean up)
next: check dead code, bump version in pyproject.toml, update CHANGELOG.md
```

### 7. Clean Up
- Check for dead code in files you touched — clean it up
- Bump the version in `pyproject.toml` (semver)
- Update `CHANGELOG.md`: `v{version}: {summary}, see PR #{N}`

Include a `workflow_state.md` update in the cleanup commit:

```
step: 8 (self-review)
next: run code-reviewer agent, fix blocking items, run doc-parrot, write findings file
```

### 8. Self-Review (MANDATORY — DO NOT SKIP)

**You must run the code-reviewer agent before pushing. This step is not optional.** If you push without running the reviewer, the PR will be rejected.

Call the code-reviewer agent now:
```
Use the code-reviewer agent to review my changes
```

Wait for its findings. Fix all **Blocking** items before proceeding. **Should Fix** items should also be addressed unless there's a good reason not to.

**Do not proceed to step 8b until the reviewer has run and all blocking items are resolved.**

### 8b. Doc-Parrot (MANDATORY — DO NOT SKIP)

After the code-reviewer, run the doc-parrot skill to validate docstring–implementation alignment for every changed callable in the diff:

```
/doc-parrot
```

The parrot identifies every changed callable, extracts docstring prose, and asks you to derive a description and usage example from prose alone — then compare that against the implementation. For each callable, record a judgment: **Fix docstring**, **Fix code**, or **OK**. Act on any **Fix** judgments before proceeding.

The parrot is not a gate. It produces an artifact you act on with judgment.

### 8c. Write (or update) findings file

After both the code-reviewer and the doc-parrot have run, append a round section to
`{package}/plans/issue-{N}-slug-findings.md`. If the file doesn't exist yet, create it
with the file header first. Each review iteration adds a new round — do not overwrite
previous rounds.

```markdown
# Findings: issue-{N}-{slug}

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

If you fix blocking items and re-run steps 8 and 8b, append `## Round 2` (and so on)
with the same structure. Commit after each round. Do not push until the latest round
shows `Blocking: 0`.

Include a `workflow_state.md` update in the findings commit:

```
step: 9 (push and PR)
next: run pre-commit, push, create PR, delete workflow_state.md
```

### 9. Push and PR

```bash
pre-commit run --all-files
git rm {package}/plans/workflow_state.md
git commit -m "chore: remove workflow state (PR filed)"
git push origin HEAD
```

Create the PR and capture the PR number and commit SHA:
```bash
PR_URL=$(api_pr_create "{concise title}" "Closes #{N}

{summary}

Plan: {package}/plans/issue-{N}-slug.md" "main")
PR_NUMBER=$(echo "$PR_URL" | grep -oP '\d+$')
COMMIT_SHA=$(git rev-parse HEAD)
```
Reference all related issues in the PR body.

### 10. Post PR Annotations

Using the `annotation-data` JSON block from the code-reviewer's output (step 8 Self-Review), post targeted annotations to help the human reviewer.

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
- If the annotation-data has no items to flag, post a simplified tour: "All changes are routine. No items flagged for special attention."
- If a line comment fails (line not in diff hunk), skip it — don't let annotation failures block the PR.
- Limit to 15 line-specific comments maximum.
- All `api_*` functions target the upstream repo automatically (detected from the `upstream` git remote) — no extra flags needed.

### 11. Issue Hygiene
```bash
api_issue_comment {N} "Addressed in PR #${PR_NUMBER}. {brief summary}"
```
Comment on related issues the same way.

## Arguments

$ARGUMENTS
