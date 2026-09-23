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

1. `CLAUDE.md` (already loaded — project rules and ALAMO architecture notes).
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
4. `docs/project/ARCHIVE_DONE.md` — only if you need historical detail to summarize a
   just-completed step or to disambiguate context for the next step.
5. **Long sources, read only the slice you need:**
   - `docs/project/in-main-tex-you-will-quizzical-treasure.md` — full historical plan and
     checklist. Read the section for the upcoming step.
   - `main.tex` — equations and physics. Read only the relevant section.
   - The relevant validation `.txt` file (Zhang/Oglesby, Hu, Kant/von Rohr) if
     and only if the next step is that validation.
6. Inspect current code/tests that the next step will touch — just enough to
   make `ACTIVE_STEP.md` concrete (file paths, existing patterns to mirror).

## What to produce

### 1. If the previous step is completed, archive it first

- Move the durable summary and `Implementation takeaways` from
  `ACTIVE_STEP.md` into `docs/project/ARCHIVE_DONE.md` (append; do not rewrite history).
- Update the relevant completed-step section in
  `docs/project/in-main-tex-you-will-quizzical-treasure.md` with a compact
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
- `## Sources` — exact pointers (file:section) into `main.tex`, the long
  plan, or validation `.txt` files. Quote constants, equations, and pass
  criteria directly here so the implementer does not need to open the long
  files.
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

## Verification before writing (mandatory)

Do these against the code and the frozen studies **before** drafting; carry 1–4
in the packet. **Guardrails are assertions too — a wrong one suppresses the check
that would have caught it, so never write one that forbids a cheap verification.**

1. **Reachability ledger** — one row per control/test pair: the exact key diff,
   the code site reading each key, and the guard condition there. A named run
   resolves to keys read from its `run.py`, never to what its name implies.
2. **Legs as literal key lines**, not prose — "at default X" must be unwriteable.
3. **Sources** — every constant, measurement and literature number carries a
   location you opened *this session*. Take constants from where they are
   **consumed**, not declared.
4. **Refuting evidence** — list the on-disk runs that would **refute** each
   prediction. If one exists, the measurement governs and the conflict is the
   finding.
5. **"What would make this inert?"** — enumerate the non-physics explanations for
   a null and give each a witness diagnostic.
6. **Context carries only what is needed to execute**; argument and history go to
   the daily log and `ARCHIVE_DONE.md`.

## Guardrails for the planner

- Do **not** edit production code, tests, or `src/`. This phase only touches
  `ROADMAP.md`, `ACTIVE_STEP.md`, `docs/project/ARCHIVE_DONE.md`, and
  `docs/project/in-main-tex-you-will-quizzical-treasure.md`.
- Do **not** rewrite the long `docs/project/in-main-tex-you-will-quizzical-treasure.md`
  unless the user explicitly asks. Append-only updates to the matching step's
  `Implementation summary` are fine.
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
