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

## Where the pre-implementation check lives — NOT here

⚠⚠ This spec used to describe running `/verify` on an unimplemented packet. **The loop
does not do that** (`/plan` → `/implement` → `/verify`), and a check that depends on a
phase nobody runs is not a check. The ledger-against-code read now lives in
**`/implement`'s step zero**, which is where it already executes and where it has
already stopped a packet. The planner's own re-read is scoped to **internal
consistency only** and says so.

⇒ **Three readers, three distinct jobs, no overlap:** the planner checks the packet
against itself; the implementer checks the packet against the code and stops on a failed
row; the reviewer checks the outputs against the packet, independently. If any two of
those start describing the same job, the duplication these specs exist to remove has
come back.

## Inputs you should read, in order

1. `CLAUDE.md` (already loaded).
2. `ACTIVE_STEP.md` — the spec, plus the implementer's `Claude completion
   notes` and `Implementation takeaways`.
3. `ROADMAP.md` — to understand which prior validations must still hold.
4. The diff and the touched code/tests. Inspect them directly; do not
   trust the implementer's narrative.
5. `docs/project/ARCHIVE_DONE.md` — only if you need historical context for regressions
   that must remain green.

## What to check

Walk through these in order and write each conclusion down:

1. **Scope match.** Does the implementation address every numbered
   deliverable in `ACTIVE_STEP.md`? Anything skipped, hand-waved, or
   silently expanded?
2. **Correctness.** Read the actual code. Does it implement what the
   `Sources` section says (equations, constants, BC patterns, pass
   criteria)? ⚠⚠ Reading without a deliverable is not a check — the
   reading is written down as the `### Diff review` table below, one row
   per hunk.
3. **Guardrails honored.** Default-off behavior is actually default-off?
   No edits to `ext/`, `bin/`, `obj/`, `build/`, `compile_commands.json`,
   `configure`, or `LICENSE`? No unscoped refactors?
4. **Tests really ran.** Check the existing test outputs (output dirs,
   logs, `thermo.dat`, plotfile timestamps vs the binary) against the
   numbers in the completion notes. Do **not** rerun tests or rebuild; see
   the rerun rule under Guardrails.
5. **No regressions on completed steps.** If shared files (operators,
   integrator base, microstructure paths) changed, read the diff carefully
   enough to be sure, and check the implementer's recorded witness results.
6. **Caveats preserved.** Anything `ROADMAP.md` flags under
   `## Known Stale/Important Notes` that this step might have touched —
   still respected?
7. **Takeaways are durable.** The `Implementation takeaways` section should
   capture gotchas, design choices, input names, and caveats. Are they
   useful enough that the next planner can write a good packet from them?

## Diff review — the code itself, as quality control

⚠⚠ Why this exists: three defects in ONE `SurfaceCellFlux` call (Pa vs MPa; the Newton
length; the coefficient cap) and a hard-coded gas temperature in a scoring script all
passed review, because the reviewer's only deliverable was the criteria table (log
2026-10-03 §2–§3). The diff is read **with a checklist and an output**, not "for bugs".

**Scope:** the diff, plus every call site and sibling the touched code mirrors. Test and
scoring scripts are code. No reruns. Time-boxed — this is a second pair of eyes, not a
second implementation.

**One row per hunk, or per new/changed call, comparison, conditional or test**, tagged
**bug / risk / nit / ok**, against these classes — each of which has shipped here:

1. **Duplication.** For every new function or block, name the existing code that computes
   the same thing, or write `none`. A copy of existing arithmetic or of a sibling call
   (the four lateral flux sites were four copies) needs a stated reason; a recurring copy
   becomes a shared helper or a lock test, not a fifth copy.
2. **Consistency across consumers.** A changed call's arguments against its siblings
   (kernel / march / ledger / criterion); a conditional that exists elsewhere must be the
   **same** conditional.
3. **Units at every comparison and table lookup** (MPa·√m vs Pa·√m, m vs mm, K vs °C).
4. **Scope predicate vs data source** — a value read only where it is zero by scope.
5. **Hard-coded quantities the run logs**, in source **and** in scripts; parse or invert
   them from the run instead.
6. **Floating-point equality, or `>=` at an exact boundary**; bars below the logging
   precision (`%g` is six digits).
7. **New keys and fields have readers.** Each new key is read on the path the packet
   claims and is default-off behind the stated conjunction; each new field has a consumer
   and a logged column. Inert keys have shipped three times.
8. **Accounting.** Every new removal or sink path has an **unconditional** ledger term.
9. **Removed, relaxed or deleted checks** — list every one with the deliverable that
   authorised it. An unlisted relaxation is a must-fix.
10. **Every new test is RED on the pre-fix binary.** Name the archived pre-fix output it
    fails on. A test that passes on both binaries tests nothing (the clamped cavity
    assertion passed by construction for two steps).
11. **Comments the diff makes false** — flag any comment in a touched hunk that no longer
    describes the code, and any comment quoting the bug the code fixes.
12. **Self-limited quantities used as discriminators; thresholds carried from a different
    experiment or population.**
13. **Guards on every new division or logarithm** (`ell_eff` reaches 3 km on a flat
    column), and steps-per-second from the run log against the previous leg.

**Tags:** a **bug** touching a gating row or a standing rule is **must-fix** ⇒
`rejected`. A **risk** is recorded in the deferred list **by name**. A **nit** ships.
Where a class is already covered by an automated lock, the row says `covered by <test>`
instead of re-reading it.

## Verdict

Write a `## Review findings` section at the **bottom of `ACTIVE_STEP.md`**. Its body is
**the packet's `## Criteria` table again**, with two columns added: the **value you
derived independently** and the **file you read it from**. Not the implementer's number
re-typed — your own reduction of the output. A row you could not re-derive is a row that
is **unbacked**, and ⚠ before calling it that, check any scratch path the packet's
commands name. ⚠ A cited artefact that is not in the repository is **unbacked**, whatever
it says.

Then the **`### Diff review` table** from the section above — one row per hunk, tagged.

- **`accepted`** means **every gating row reproduces AND no diff-review row is tagged
  `bug` on a gating row or standing rule.**
- **`rejected`** names the rows that do not, and the `bug` rows.
- **must-fix** is a gating row, a standing rule, or a `bug` against one; everything else
  is a nit or a named risk.
- ⚠ Confirm the completion notes **quote nothing the packet marked unquotable**, checked
  against **both** `ROADMAP.md`'s `## Standing rules` (programme invariants) and
  `## Known Stale/Important Notes` (step caveats).

Then one of these verdicts:

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
- **Do not rerun tests, simulations, or builds by default.** `/implement`
  already ran them, so rerunning wastes compute. Rerun one only when it is
  really necessary: when you have a concrete reason to think the
  implementer's results missed an important detail. Examples are output
  that is missing or older than the code, numbers that contradict the
  notes, or a required check that was never run. Rerun just that test, and
  say in the findings why you did.
- Do not edit production code, tests, or `src/`. The only file you write
  to is `ACTIVE_STEP.md`, and only the `Status:` line and a new
  `## Review findings` section.
- Do not move work into `docs/project/ARCHIVE_DONE.md` or
  `docs/project/in-main-tex-you-will-quizzical-treasure.md`. That belongs to the next
  `/plan` call after acceptance.
- If you are unsure whether something is in scope, ask the user before
  rejecting.

## End of phase

Report:
- One-line verdict (`accepted` or `rejected`).
- One-line summary of the most important finding (positive or negative).
- The diff-review tally: how many rows, how many `bug` / `risk` / `nit`, and the one
  duplication or consistency finding most worth the planner's attention.
- The exact next step the user should take:
  - on accepted: "run `/plan` to archive this step and prep the next packet."
  - on rejected: "run `/plan` to revise `ACTIVE_STEP.md` based on the review
    findings."
