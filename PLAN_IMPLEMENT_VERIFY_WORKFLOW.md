# Plan / Implement / Verify — A Three-Phase Claude Code Workflow

A portable export of the three-phase workflow used in this project. Drop it
into another repo, swap the project-specific bits, and you have the same
loop: planner → implementer → reviewer, each a fresh Claude Code session.

## What this workflow gives you

A pipeline that decouples three responsibilities that the same Claude tends
to conflate when given a complex multi-step task:

1. **Planner** (`/plan`) — reads long historical context, writes a focused
   task packet (`ACTIVE_STEP.md`) with everything the next coder needs.
   Never writes production code.
2. **Implementer** (`/implement`) — reads only the packet (not the long
   plan), writes code + tests, appends completion notes and takeaways.
   Never re-plans the next step.
3. **Reviewer** (`/verify`) — independent second pass. Reads the diff and
   the spec; writes an `accepted` or `rejected` verdict. Never edits code
   or rewrites the packet.

Each phase is a separate `claude` session. The phases pass context through
files, not through one long conversation, which keeps each call's context
window focused on its single job.

## When this workflow earns its keep

- **Multi-step engineering work** where a single Claude run drifts: scope
  creeps, takeaways get forgotten, the next session has to re-derive
  context from scratch.
- **Validation-heavy work** (numerical sims, ML training pipelines,
  compiler passes) where each step has explicit pass/fail criteria that
  must be quoted into the packet, not paraphrased.
- **Long-running projects** where you'll come back next week and need to
  pick up where you left off. The packet is the durable handoff.

When NOT to use it: one-off scripts, quick bug fixes, exploratory
prototyping. The packet overhead doesn't pay off below ~3 days of work.

## Files in the workflow

```
.claude/commands/plan.md         # /plan slash command (template below)
.claude/commands/implement.md    # /implement slash command (template below)
.claude/commands/verify.md       # /verify slash command (template below)
CLAUDE.md                        # Project entry point: rules, commands, context loading order
ROADMAP.md                       # Short status + active step pointer
ACTIVE_STEP.md                   # Current task packet (the only "code" Claude reads at /implement)
ARCHIVE_DONE.md                  # Completed-step history (append-only)
<your-long-plan>.md              # Full historical plan / equations / references (planner reads slices only)
```

The split between `CLAUDE.md` (operational, short) and your long-plan
file (detailed, never loaded by default) is the key trick: the long plan
stays comprehensive without poisoning every session's context budget.

## The loop

```
/plan       -> ACTIVE_STEP.md is fresh and ready
/implement  -> code + tests + completion notes appended to ACTIVE_STEP.md
/verify     -> Review findings appended; verdict: accepted | rejected
  if accepted: /plan   (archives takeaways, preps the next packet)
  if rejected: /plan   (revises ACTIVE_STEP.md based on review findings, then /implement again)
```

Each invocation is a separate `claude` session. The slash commands are
defined in `.claude/commands/*.md` and Claude picks them up automatically.

## Why three phases (not one, not two)

- **One phase** (just /implement) drifts: scope creeps as Claude finds
  "while I'm here" cleanups, the takeaways never get written down, the
  next session re-derives context.
- **Two phases** (plan + implement) catches scope drift but misses
  correctness bugs. The implementer is invested in the work; an
  independent reviewer catches things the implementer rationalized.
- **Three phases** is the minimum that keeps each session focused on a
  single role. The reviewer doesn't have to be "fair" — they re-derive
  from the spec.

The reviewer phase is the highest-leverage one. The planner has the
hardest job. The implementer is the most fun.

---

# Drop-in templates

Copy the four blocks below into a new project. Rename `MY_PROJECT_LONG_PLAN.md`
and the validation-reference filenames to match your repo. The slash commands
are largely project-agnostic; the only thing you need to swap is the
"long sources" list in `/plan` and the build/run commands in `/implement`.

## 1. `.claude/commands/plan.md`

```markdown
---
description: Plan the next coding step. Rewrites ACTIVE_STEP.md and refreshes ROADMAP.md from the long plan.
---

# /plan — Planner Phase

You are the **planner**. Your job is to refresh the focused task packet so a
later `/implement` call can code the next step without opening the long
historical plan.

You do **not** write production code in this phase. You only read context and
rewrite planning files.

## Inputs you should read, in order

1. `CLAUDE.md` (already loaded — project rules and architecture notes).
2. `ROADMAP.md` — current status, active milestone, next milestones.
3. `ACTIVE_STEP.md` — the current packet. Check its status first:
   - If status is `completed`, you are starting a fresh step. Read the
     `Claude completion notes` and `Implementation takeaways` so they can be
     archived.
   - If status is `ready for implementation` or similar, the step has not run
     yet — confirm with the user whether they want to re-plan or skip planning.
   - If status is `rejected by review` (or there is a `Review findings`
     section), revise the packet to address the reviewer's concerns; do not
     simply move on to the next step.
4. `ARCHIVE_DONE.md` — only if you need historical detail to summarize a
   just-completed step or to disambiguate context for the next step.
5. **Long sources, read only the slice you need:**
   - `<MY_PROJECT_LONG_PLAN>.md` — full historical plan and checklist. Read
     the section for the upcoming step.
   - `<spec/equations/references>` — read only the relevant section.
6. Inspect current code/tests that the next step will touch — just enough to
   make `ACTIVE_STEP.md` concrete (file paths, existing patterns to mirror).

## What to produce

### 1. If the previous step is completed, archive it first

- Move the durable summary and `Implementation takeaways` from
  `ACTIVE_STEP.md` into `ARCHIVE_DONE.md` (append; do not rewrite history).
- Update the relevant completed-step section in the long plan with a compact
  `Implementation summary` of **20 lines or fewer**. Summarize Claude's longer
  notes; do not paste raw logs, verbose diffs, or transient debugging detail.
- Update the `## Current Status` section of `ROADMAP.md` to list the
  just-completed step.

### 2. Rewrite `ACTIVE_STEP.md` for the next coding step

`ACTIVE_STEP.md` must contain everything the implementer needs without
opening the long plan. Required sections:

- `# Active Step: <N> - <short title>`
- `Status: ready for implementation`
- `## Context` — one short paragraph on why this step exists and what came
  before. Carry forward any caveats from prior steps that constrain this one.
- `## Sources` — exact pointers (file:section) into the long plan or
  reference files. Quote constants, equations, and pass criteria directly
  here so the implementer does not need to open the long files.
- `## Goal` — numbered deliverables (files to add/edit, functions to
  implement, regressions to add).
- `## Guardrails` — what *not* to touch, what must remain default-off, what
  earlier validation thresholds must remain green, etc.
- `## Commands` — exact build, run, and test commands for this step
  (mirror the patterns in `CLAUDE.md`).

Keep it tight: equations, constants, file paths, validation criteria, and
guardrails go in; verbose prose stays out.

### 3. Refresh `ROADMAP.md`

- `## Current Status` lists the just-completed step (if any).
- `Active step:` names the new step and points at `ACTIVE_STEP.md`.
- `Next milestones after the active step:` lists the next few in one-line
  form.
- Preserve and update `## Known Stale/Important Notes` — do not delete
  caveats from prior steps unless the user explicitly says they are resolved.

## Guardrails for the planner

- Do **not** edit production code, tests, or `src/`. This phase only touches
  `ROADMAP.md`, `ACTIVE_STEP.md`, `ARCHIVE_DONE.md`, and the long plan.
- Do **not** rewrite the long plan unless the user explicitly asks.
  Append-only updates to the matching step's `Implementation summary` are
  fine.
- Preserve warnings about stale references, changed assumptions, input names,
  and validation caveats.
- Remove transient debugging detail from the new `ACTIVE_STEP.md`; keep only
  what the implementer needs.
- If you are unsure which step to plan next, ask the user before rewriting
  anything.

## End of phase

Report:
- Which step was archived (if any), with a one-line summary.
- Which step is now active, with its title and the one-line goal.
- Any unresolved caveats the implementer must respect.
```

## 2. `.claude/commands/implement.md`

```markdown
---
description: Implement the active step. Reads ACTIVE_STEP.md and writes code, tests, and completion notes.
---

# /implement — Implementer Phase

You are the **implementer**. Your job is to code the active step described in
`ACTIVE_STEP.md` and produce a runnable, tested result.

## Inputs you should read, in order

1. `CLAUDE.md` (already loaded — project rules, architecture, commands, and
   do-not-touch list).
2. `ROADMAP.md` — confirm the active step matches what `ACTIVE_STEP.md` says.
3. `ACTIVE_STEP.md` — your task packet. This is the source of truth for
   scope, sources, goal, guardrails, and commands.
4. The specific source/test files named in `ACTIVE_STEP.md`. Read existing
   patterns nearby and mirror them.

Open the long sources **only** if the active task explicitly calls for a
constant or equation that is not already inlined in `ACTIVE_STEP.md`. If you
find yourself reaching for them often, the packet is under-specified —
stop and ask the user to run `/plan` again.

## What to produce

1. **Code** that satisfies every numbered deliverable in `ACTIVE_STEP.md`.
   - Mirror existing patterns.
   - Default-off any new behavior unless the step explicitly turns it on.
   - Do not touch the do-not-touch paths listed in `CLAUDE.md`.
   - Do not regress completed validation tests.

2. **Tests / regressions** as required by the step. Mirror existing test
   patterns in the repo.

3. **Build and run** the regressions named in `ACTIVE_STEP.md`. If a
   regression cannot run on this machine, say so explicitly — do not claim
   success on un-run tests.

4. **Update `ACTIVE_STEP.md` at the end of the phase**:
   - Change `Status:` to `completed` (or `blocked` if you genuinely cannot
     finish; explain).
   - Append `## Claude completion notes` with: files changed, tests run, test
     result (pass/fail/skipped), and any tests not run with the reason.
   - Append `## Implementation takeaways` with: gotchas, design choices,
     input parameter names, changed assumptions, caveats. Write these for the
     reviewer and the next planner — they are the durable record of this
     step.

Do **not** rewrite the next `ACTIVE_STEP.md`. That is the planner's job. Do
**not** archive yourself in `ARCHIVE_DONE.md` or the long plan. The planner
moves the takeaways later.

## Guardrails for the implementer

- Stay inside the scope of `ACTIVE_STEP.md`. If a needed change is outside
  scope, stop and ask — do not silently expand the step.
- Do not invent validation constants or reference data. If a constant is
  missing from `ACTIVE_STEP.md`, ask, or pull it from the source files
  explicitly named in the packet.
- Keep new behavior default-off unless the step explicitly enables it.
- Preserve completed regression behavior. If you must touch a shared file,
  re-run the affected regressions.
- Do not amend prior commits or force-push. Create new commits only if the
  user asks.
- Run UI/feature checks where the step requires them. Type-check and tests
  verify code correctness, not feature correctness — if you cannot exercise
  the feature, say so explicitly.

## Build and run

Use the commands block in `ACTIVE_STEP.md` first. Fall back to the
project-level commands in `CLAUDE.md`.

## End of phase

Report:
- Files changed.
- Tests run and their pass/fail status.
- Any tests not run, with the reason.
- A pointer to the `## Claude completion notes` and
  `## Implementation takeaways` sections you just appended to
  `ACTIVE_STEP.md`.
```

## 3. `.claude/commands/verify.md`

```markdown
---
description: Review the implementer's work against ACTIVE_STEP.md. Flags issues only; does not rewrite the packet.
---

# /verify — Reviewer Phase

You are the **reviewer**. Your job is to verify that the implementer's
output actually satisfies the active step. You are a secondary correctness
check — independent of the implementer.

You do **not** rewrite `ACTIVE_STEP.md` for the next step, and you do
**not** revise the active packet when you find problems. If review fails,
flag the issues; the user will run `/plan` again to revise.

## Inputs you should read, in order

1. `CLAUDE.md` (already loaded).
2. `ACTIVE_STEP.md` — the spec, plus the implementer's `Claude completion
   notes` and `Implementation takeaways`.
3. `ROADMAP.md` — to understand which prior validations must still hold.
4. The diff and the touched code/tests. Inspect them directly; do not
   trust the implementer's narrative.
5. `ARCHIVE_DONE.md` — only if you need historical context for regressions
   that must remain green.

## What to check

Walk through these in order and write each conclusion down:

1. **Scope match.** Does the implementation address every numbered
   deliverable in `ACTIVE_STEP.md`? Anything skipped, hand-waved, or
   silently expanded?
2. **Correctness.** Read the actual code. Does it implement what the
   `Sources` section says (equations, constants, BC patterns, pass
   criteria)?
3. **Guardrails honored.** Default-off behavior is actually default-off?
   No edits to the do-not-touch paths in `CLAUDE.md`? No unscoped refactors?
4. **Tests really ran.** Inspect test output or rerun the regression named
   in `ACTIVE_STEP.md`. The completion notes claim something — verify it.
5. **No regressions on completed steps.** If shared files changed,
   sample-run the affected prior regressions or read the diff carefully
   enough to be sure.
6. **Caveats preserved.** Anything `ROADMAP.md` flags under
   `## Known Stale/Important Notes` that this step might have touched —
   still respected?
7. **Takeaways are durable.** The `Implementation takeaways` section should
   capture gotchas, design choices, input names, and caveats. Are they
   useful enough that the next planner can write a good packet from them?

## Verdict

Write a short `## Review findings` section at the **bottom of
`ACTIVE_STEP.md`** with one of these verdicts:

### `Verdict: accepted`

- Bullet what was checked and how (which tests rerun, which files
  inspected). Be specific — file paths and a sentence each.
- The user will run `/plan` next; the planner will archive the step.

### `Verdict: rejected`

- Bullet each problem with a one-line statement, a pointer to file:line or
  test name, and what the implementer or planner needs to do.
- Group findings as **must-fix** (blocks acceptance) vs **nit** (acceptable
  if the user decides to ship). Do not block on nits.
- Update the top of `ACTIVE_STEP.md` to set
  `Status: rejected by review` so the next `/plan` call sees it.
- Do **not** rewrite the `Goal`, `Sources`, or `Guardrails` sections — that
  is the planner's job.

## Guardrails for the reviewer

- You are the second pair of eyes. Re-derive your verdict from the diff
  and the spec; do not just paraphrase the implementer's notes.
- Do not edit production code, tests, or `src/`. The only file you write
  to is `ACTIVE_STEP.md`, and only the `Status:` line and a new
  `## Review findings` section.
- Do not move work into `ARCHIVE_DONE.md` or the long plan. That belongs to
  the next `/plan` call after acceptance.
- If you are unsure whether something is in scope, ask the user before
  rejecting.

## End of phase

Report:
- One-line verdict (`accepted` or `rejected`).
- One-line summary of the most important finding (positive or negative).
- The exact next step the user should take:
  - on accepted: "run `/plan` to archive this step and prep the next packet."
  - on rejected: "run `/plan` to revise `ACTIVE_STEP.md` based on the review
    findings."
```

## 4. `CLAUDE.md` template (project root)

```markdown
# CLAUDE.md

This file is the small entry point for Claude. Keep it short. The long
project plan lives in `<MY_PROJECT_LONG_PLAN>.md`, but Claude should not
read that file by default.

## Three-Phase Workflow

This project is driven by three separate Claude invocations. Each one is a
slash command defined under `.claude/commands/`. The same Claude does all
three roles, but each phase is a fresh call with a focused responsibility.

1. **`/plan`** — Planner. Reads the long sources, archives any completed
   step, and rewrites `ACTIVE_STEP.md` so the next coding step is fully
   specified. Does not write production code.
2. **`/implement`** — Implementer. Reads `ACTIVE_STEP.md` and writes the
   code, tests, and regressions. Updates `ACTIVE_STEP.md` with completion
   notes and takeaways at the end. Does not rewrite the next packet.
3. **`/verify`** — Reviewer. Verifies the implementer's output against
   `ACTIVE_STEP.md`. Writes a `## Review findings` section with an
   `accepted` or `rejected` verdict. On rejection, the user re-runs `/plan`
   to revise the packet. The reviewer does not edit code or rewrite the
   packet itself.

Run them in separate Claude sessions. Each command's spec lives in
`.claude/commands/<phase>.md`.

The typical loop:

```
/plan       -> ACTIVE_STEP.md is fresh and ready
/implement  -> code + tests + completion notes
/verify     -> accepted | rejected
  if accepted: /plan (archives + preps next step)
  if rejected: /plan (revises ACTIVE_STEP.md, then /implement again)
```

## Context Loading Order

After opening this file, read only these additional files by default:

1. `ROADMAP.md` — current status, active milestone, and next milestones.
2. `ACTIVE_STEP.md` — the detailed task packet for the active coding step.

Then use the project rules, commands, gotchas, and guardrails below.

Open the long/reference files only when the active phase explicitly
requires them or the user asks:

- `<MY_PROJECT_LONG_PLAN>.md` — full historical plan and checklist. The
  planner phase reads slices of this; the implementer should not need it.
- `<other-references>` — domain references, equations, validation data.
- `ARCHIVE_DONE.md` — completed-step history. Read only for debugging
  regressions or when the planner is updating project context.

## Project Goal

<one paragraph: what we are building, for whom, and the success criterion.>

## Architecture Notes

<short list of project conventions: file layout, patterns to mirror,
language/framework gotchas, anti-patterns.>

## Lessons Learned

<numbered list of project-specific gotchas that bit a past session.
Update sparingly; this is the operational memory.>

## Commands

```bash
# Build
<your build command>

# Run / test
<your test command>
```

## Do Not Touch

- <vendored deps>
- <build artifacts>
- <auto-generated files>
- <licenses, config files Claude shouldn't rewrite>

## Context Budget Rules

- Keep `CLAUDE.md` short and operational.
- Keep `ROADMAP.md` short enough to read every session.
- Put detailed completed history in `ARCHIVE_DONE.md`, not `ROADMAP.md`.
- Keep each completed-step summary in `<MY_PROJECT_LONG_PLAN>.md` to 20
  lines or fewer.
- Put only the current coding packet in `ACTIVE_STEP.md`.
- Include equations, constants, file paths, validation criteria, and
  guardrails in `ACTIVE_STEP.md` when the implementer needs them.
- Prefer summaries and exact pointers over copying large blocks from the
  long plan.
```

## 5. `ROADMAP.md` template

```markdown
# Roadmap

## Current Status

- [completed step N]: <one-line summary>
- [completed step N-1]: <one-line summary>

Active step: **Step <N+1> — <title>** (see `ACTIVE_STEP.md`).

Next milestones after the active step:
- Step <N+2>: <one-line goal>
- Step <N+3>: <one-line goal>

## Known Stale/Important Notes

- <caveat from a prior step that still constrains future work>
- <validation threshold that must stay green>
- <deferred-but-not-forgotten issue>
```

## 6. `ACTIVE_STEP.md` template

```markdown
# Active Step: <N> — <short title>

Status: ready for implementation

## Context

<one short paragraph: why this step exists, what came before, what caveats
from earlier steps still apply.>

## Sources

<exact pointers into the long plan + reference files. Quote the constants,
equations, and pass criteria directly here. The implementer should not need
to open the long files.>

## Goal

1. <numbered deliverable: file to add or edit, function to implement>
2. <numbered deliverable: test or regression to add>
3. <numbered deliverable: build/run verification>

## Guardrails

- <what NOT to touch>
- <what must remain default-off>
- <earlier validation thresholds that must stay green>

## Commands

```bash
# Build
<exact command>

# Run / test
<exact command>
```

## Expected outcomes

- <if-pass: what success looks like>
- <if-fail-mode-1: what diagnostic to check>
- <if-fail-mode-2: what diagnostic to check>
```

## 7. `ARCHIVE_DONE.md` template

```markdown
# Archive of completed steps

Append-only. Each entry is the durable summary + takeaways from the
implementer phase, archived by `/plan` after `/verify` accepted the step.

---

## Step <N> — <short title>  (completed YYYY-MM-DD)

### Summary

<one paragraph: what this step delivered and the headline result.>

### Files changed

- `path/to/file.ext` — <one-line description of the change>

### Implementation takeaways

1. <gotcha or design choice worth carrying forward>
2. <input parameter name, default value, why>
3. <caveat for future planner: what this step does NOT cover>
```

---

# Setup checklist (new project)

1. Create the four files: `CLAUDE.md`, `ROADMAP.md`, `ACTIVE_STEP.md`,
   `ARCHIVE_DONE.md`.
2. Create `.claude/commands/` and drop the three slash command files in.
3. Write the long plan (`<MY_PROJECT_LONG_PLAN>.md`) — the full
   step-by-step roadmap with equations/constants/references. This is the
   only place verbose project context lives.
4. Fill in `CLAUDE.md`'s project goal, architecture notes, build commands,
   and do-not-touch list.
5. Write the first `ACTIVE_STEP.md` by hand. After that, `/plan` writes it.
6. Run `/plan` in one Claude session, `/implement` in another,
   `/verify` in a third. Iterate.

# Tips from running this workflow in practice

- **Quote, don't paraphrase, in `ACTIVE_STEP.md`'s `## Sources`.** If the
  pass criterion is "rel err ≤ 8%", paste that string in the packet.
  Paraphrased criteria drift between phases.
- **The reviewer should re-run at least one regression.** Not all of them
  — the implementer already did — but enough to confirm the test output
  is real and not invented.
- **`Implementation takeaways` is the most valuable artifact.** Future
  planners read it instead of the diff. If a takeaway is a bullet list of
  "I edited X.cpp", it's not durable. Good takeaways read like:
  "Parser key `confining.p` is parsed BEFORE the model-branch early
  return so command-line overrides work under both branches. Default-zero
  short-circuit is mandatory for byte-identity — verified by X test."
- **Don't skip `/verify` because it feels like overhead.** Most of the
  bugs caught by this workflow get caught in `/verify`, not `/implement`.
- **If `/verify` rejects, `/plan` revises — do not have `/implement`
  patch directly.** The planner is the one who reframes scope; the
  implementer treats the packet as a contract.
- **One step at a time.** Multi-step packets defeat the purpose. If a
  step is too big to fit in `ACTIVE_STEP.md` cleanly, split it.
- **The long plan is append-only inside its step sections, but the
  planner is free to add new step sections.** Don't let the long plan
  become an archive of stale ideas — prune it deliberately, but never
  silently.

# What you'll need to swap when porting

- The long-plan filename in `/plan`'s "Inputs" section.
- The build/run commands in `CLAUDE.md` and `/implement`.
- The do-not-touch list in `CLAUDE.md`.
- The "validation reference" file list in `/plan` (replace with your
  spec/equation/data sources).
- The mention of MPI ranks / test-pattern conventions in `/implement` —
  swap for your project's test command pattern.

The phase structure itself is project-agnostic. Don't change it on the
first attempt; let it run a few cycles before you decide what's missing.
