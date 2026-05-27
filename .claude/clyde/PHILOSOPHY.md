# Paradigm Philosophy

A short account of the *why* behind the agent-driven development paradigm
implemented by the skills in `.claude/skills/` and the supporting files in
`.claude/clyde/`. The mechanics are documented in the skill files and in
`CLAUDE.md`. This document is about the ideas those mechanics exist to
enforce.

This paradigm is provisionally called **Junior Engineer Clyde** — a name
that's allowed to age out and be renamed later if something stickier
emerges. Throughout this doc we mostly just say "the paradigm". The name
matters less than the rules.

This file is deliberately implementation-agnostic: no specific package
names, no specific CI system, no specific issue tracker quirks. If you
find yourself wanting to edit this document to add repo-specific content,
stop — that content belongs in `CLAUDE.md` or a package README, not here.

## Core principle: the agent is a collaborator, not a tool

The skills treat the agent as a junior collaborator with full read access,
bounded write access, and a clear handoff protocol to a human reviewer.
The paradigm's job is to make that collaboration legible — so a human can
quickly see what the agent did, why, and where to push back.

Every rule below falls out of that principle.

## Plan before code

Every non-trivial change begins with a written plan: verbatim goal, chosen
approach, alternatives rejected, explicit scope boundaries, assumptions
about neighbouring code.

**Why.** The expensive part of reviewing an agent's PR is *not* reading
the diff — it's reconstructing what the agent was trying to do and
deciding whether that was the right thing. A plan front-loads that
reconstruction. A missing plan means the reviewer has to infer intent
from code, which is exactly the activity the paradigm exists to avoid.

The plan is a living document until the PR is filed. If the agent pivots
during implementation, the plan gets updated to reflect what actually
happened. Stale plans are worse than no plans.

## Fork-and-PR, always

Agents work on a fork (or a feature branch that behaves like one) and
open PRs against upstream. They never push directly to the shared trunk.

**Why.** Fork-and-PR gives three things for free:

1. **A review boundary.** Every change passes through a PR, which is
   where the annotations, self-review output, and CI results live.
2. **A blast-radius cap.** A mistaken agent run cannot corrupt upstream
   history — the worst case is a bad PR, which is cheap to close.
3. **Clean attribution.** The human who launched the agent is the commit
   author, with the agent added as a `Co-Authored-By` trailer. This
   keeps the audit trail honest without inflating the agent's
   contribution graph.

## Sandbox the agent — why the VM matters

Agents operate at a permission level that lets them do real work
unsupervised. That only stays safe if the *blast radius* of a misbehaving
agent is bounded. Two complementary boundaries:

- **The fork** caps damage to upstream code (covered above).
- **The VM** caps damage to your host machine.

A short-lived virtual machine — Vagrant, a container, a remote
ephemeral host — is what lets you run the agent with the
`--dangerously-skip-permissions` flag (or its equivalent) without
gambling your laptop. Inside the VM, the agent can install packages,
edit dotfiles, run subprocesses, and otherwise have a bad day. None of
it touches the host. When something goes wrong, you `vagrant destroy`
and start over.

This is paradigm-relevant, not infrastructure detail. The skills assume
the agent has wide latitude. That assumption is only true if the
sandbox is real. Adopters who skip the VM and run agents directly on
their development machine have re-introduced a problem the paradigm
takes for granted is solved.

A sample Vagrant configuration is shipped alongside this kit (see the
`vagrant/` directory) — Ubuntu 24.04, the system packages you need,
shared folder back to the project, ready in one `vagrant up`. Adopt it,
fork it, or replace it with your own equivalent. Just don't skip the
sandbox.

## Package-scoped agents

An implementing agent works in **one package at a time**. Reading across
packages is free; modifying across packages requires explicit permission
or a coordinating agent that merges both sides atomically.

**Why.** Most agent mistakes are scope mistakes — "fixing" code in a
neighbouring package that had a reason to be the way it was. A hard
scope boundary turns those mistakes into asks ("can I touch this other
package?") rather than surprises discovered at review time.

The invariant worth preserving: **the trunk never sees a partial
cross-package change.** Either both sides land together or neither
lands.

## Tests are the paradigm's enforcement gate

Tests aren't a checkbox the agent ticks before pushing — they're how the
paradigm proves it did the work it claimed to. Three rules together
make tests load-bearing:

1. **Tests come from the ticket, not from the implementation.** The
   ticket states acceptance scenarios in plain language ("when X, then
   Y is true"). The agent's first concrete commit after the plan is the
   tests that encode those scenarios — *before* any production code.
   Writing tests after code biases them toward what was built; writing
   them from the ticket biases them toward what was asked for.
2. **The agent runs the full test suite — including style and quality
   checks — and passes them before the human ever sees the PR.** The
   self-review step is gated on this. A red CI on a freshly-opened PR
   means the self-review didn't happen, which means the PR shouldn't
   have been opened.
3. **Existing tests are inviolable without sign-off.** When integration
   tests start failing, that is a signal to investigate, not a
   signal to update the tests. Behavior changes that intentionally
   invalidate a test require an explicit decision in the PR — not a
   silent edit.

The reviewer can trust that "everything green" actually means
something. Without that trust, the human reviewer becomes the test
runner, which is the role the paradigm exists to avoid.

## Same code, two contexts

Quality checks and review logic live in the repo as code. The agent
runs them locally as a fast feedback loop; CI runs the *same* code as
an enforcement gate.

**Why.** If the local feedback loop and the CI gate are different code,
the agent gets surprised at PR time. Surprise costs reviewer attention.
Identity between the two means the agent never pushes a PR that CI will
reject — it already ran the same checks five seconds before pushing.

This shape applies to multiple things:

- **Style and lint checks** invoked via a single test entry point
  (e.g., `pytest test_code_quality.py` for Python — wrapping
  black/ruff/pylint/pydocstyle/darglint as test cases). Local feedback
  loop, CI enforcement gate, same code.
- **The reviewer agent.** Lives in the repo as a skill the agent
  invokes against its own diff during self-review. CI can invoke the
  same skill on every PR.
- **Integration tests.** The agent runs them locally during
  implementation; CI runs them as a PR acceptance gate. Same suite,
  two triggers.

A pre-commit hook can layer on top, but the paradigm doesn't depend on
it. Pre-commit failures should reproduce the agent's local-test
output, not introduce new constraints the agent didn't see coming.

## Coverage and dead code, by default

Coverage and dead-code detection aren't aspirational reviewer concerns
— they're part of the same enforcement gate as style and lint:

- **Coverage as a floor.** A per-package coverage threshold runs as
  part of the test suite. Below the floor → tests fail → agent
  doesn't ship. The floor exists to catch obviously untested code,
  not to enforce a culture of writing tests for trivial getters.
  Branch coverage matters more than line coverage for catching real
  bugs.
- **Dead code is found and deleted.** A tool like `vulture` runs
  alongside the other quality checks. When an agent touches a file,
  it cleans up dead code in that file as part of the same PR — not a
  separate ticket. Deletion is cheap; git has the history. The bar
  for keeping unused code is high: a concrete plan to use it in the
  current release.

Both belong in the same `test_code_quality.py`-style entry point as
the formatters and linters. One pytest invocation, one truth, no
out-of-band tooling.

## Done means done

Acceptance criteria must be **satisfied, not perfected.** When the
ticket's scenarios pass, the agent files the PR. It does not gold-plate,
refactor surrounding code, or hold the branch open for continuous
improvement.

**Why.** Agents left running without a clear stop condition drift. They
optimize the wrong thing, introduce unrelated changes, and produce PRs
that are hard to review because the scope is unclear. A tight "done"
rule is the discipline that keeps PRs reviewable.

Concretely: a bug fix doesn't need surrounding cleanup. A new endpoint
doesn't need refactoring the handler it lives beside. Three similar
lines is better than a premature abstraction.

## Self-review is mandatory

Before pushing, the agent runs a code-reviewer pass on its own diff.
Blocking findings must be fixed. Should-Fix findings are addressed
unless there's a good reason not to.

**Why.** Self-review catches three classes of problem before the human
sees them:

1. **Obvious mistakes** — dead imports, wrong exception types, failing
   tests the agent hadn't rerun. Cheap to fix, expensive to explain in
   PR comments.
2. **Scope violations** — edits that wandered into a neighbouring
   package or an unrelated module. Easier to drop at self-review than
   at PR review.
3. **Uncertainty surfacing** — places where the agent guessed at an
   interface. Self-review flags these so the human reviewer knows
   exactly where to focus.

A PR pushed without self-review is rejected. The paradigm treats this
as non-negotiable — not because self-review catches everything, but
because skipping it signals that the agent is not taking the review
contract seriously.

## PR annotations are the guided tour

When the PR is filed, the agent posts a summary comment and a set of
line-specific annotations derived from the self-review. The annotations
fall into three categories:

- **Boundary violations** — changes outside the agent's declared scope.
- **Architectural decisions** — choices that deserve a second opinion.
- **Uncertainties** — places where the agent guessed and wants
  confirmation.

**Why.** A large PR is unreviewable without signposting. The agent
knows which lines matter and which are mechanical; the human doesn't,
unless the agent tells them. Annotations are the handoff: "here is
where your attention actually adds value." Mechanical changes go
unannotated and are safe to skim.

The annotation protocol is what makes large agent PRs *practically*
reviewable. Without it, reviewers either rubber-stamp or die of
exhaustion.

## Docstrings are contracts, not decoration

Every public interface has a complete, correct docstring. Docstrings
get enforced mechanically (style, correctness, argument
documentation). A stale docstring is treated as a bug, not a minor
issue.

**Why.** Agents read docstrings to understand interfaces they didn't
write. If the docstring is wrong, the next agent will build on top of
a false premise. The cost of a wrong docstring is not the comment;
it's every downstream decision made from it.

This is why the paradigm enforces docstring correctness automatically
— it is not a style preference, it's a tooling guarantee for the next
agent.

## Escalate instead of spinning

When an agent gets stuck — tests won't pass, approach isn't working,
ticket was bad — it stops, comments on the issue with what it tried,
applies a `needs-human-review` label, and leaves the branch alone.

**Why.** Stuck agents don't fail gracefully by default — they thrash,
producing commits that make the problem harder to diagnose. The
escalation protocol turns "stuck" into a structured handoff: the
branch state, the agent's notes, and the issue label together form a
clear picture for a human to pick up.

The mess, intact, is diagnostic information. Cleaning it up destroys
the signal.

## Language honesty: Python-shaped tooling, language-agnostic principles

The principles in this document are language-agnostic. The *tooling*
references in the skills, the templates, and the CLAUDE.md sections
that talk about quality enforcement are Python-specific:
`python3 -m venv`, `pip install -e`, `pytest`, `pyproject.toml`,
`darglint`, `pre-commit`, the `test_code_quality.py` template.

Adopters whose stack is not Python need to translate three things:

1. **Skill bodies** — replace tool invocations in
   `.claude/skills/{implement,workon,ship}/SKILL.md` with the
   equivalent for your language (e.g., `go test ./...`, `cargo test`,
   `npm test`).
2. **Quality gate** — `test_code_quality.py` is a Python pattern: a single
   pytest entry point that wraps your formatter, linter, and docstring/type
   tooling. See the corpus repo (`clyde/implementations/tri_amr/`) for a
   working example. Build the equivalent for your language.
3. **CLAUDE.md** — its "Quality Gates" and "Prerequisites" sections
   describe Python tooling by default. Replace with your stack.

The shape of the workflow — plan-before-code, fork-and-PR,
package-scoped, tests-as-gate, self-review, PR annotations — does not
change with language. Only the toolchain references do. If you find
yourself wanting to change a workflow rule because of your language,
you've probably mis-translated; ask first.

## What the paradigm is *not*

A short list of things the paradigm does not attempt, because naming
them sharpens what it *does* attempt:

- **Automated code generation.** The agent is slow and deliberate, not
  a throughput engine. If the goal is to produce 100 PRs a day, this
  is the wrong shape.
- **A substitute for architecture work.** Agents implement within an
  existing architecture. They do not design new subsystems from
  scratch — that is human work that the agent's plans build on.
- **Protection against malicious input.** The paradigm assumes the
  tickets and repo state are good-faith. Adversarial scenarios need
  additional controls the paradigm does not supply.
- **A one-shot bootstrap for a new team.** The paradigm rewards
  investment in tickets, docstrings, and plan hygiene. A team that
  skips those will not see the benefits.

## Adopting the paradigm

The paradigm is synthesized for new repos via the `/bootstrap` skill in the
Clyde corpus repo. Bootstrap interviews you about your repo shape, auth model,
and toolchain, then writes a tailored set of files — skills, `CLAUDE.md`,
API wrapper, optional Vagrant config — directly into the target repo. See the
corpus repo for the concrete checklist and working examples.

This document is about the *why*. The corpus is about the *how*.
