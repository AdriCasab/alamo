# AGENTS.md

> **Historical document.** This file describes the previous Codex+Claude
> split workflow and is no longer active. The live workflow is documented in
> `CLAUDE.md` and the slash commands under `.claude/commands/`
> (`/plan`, `/implement`, `/review`). Keep this file for reference only.

Codex restart instructions for this thermal-spallation project.

## Role Split

- Claude is the coding worker. Claude should start from `CLAUDE.md`,
  `ROADMAP.md`, and `ACTIVE_STEP.md`, then code the active step.
- Codex is the context curator and plan rewriter. Codex should read the long
  sources only when preparing or refreshing the focused context for Claude.

## On Every Codex Restart

Read these first:

1. `CLAUDE.md`
2. `ROADMAP.md`
3. `ACTIVE_STEP.md`
4. `ARCHIVE_DONE.md` if the user asks for status/history, if
   `ACTIVE_STEP.md` says the step is completed, or if a completed step needs to
   be summarized.

Then inspect current code/tests only as needed for the user's request.

## When Preparing A New Claude Task

1. Read the relevant section of
   `in-main-tex-you-will-quizzical-treasure.md`.
2. Read the relevant section of `main.tex` for equations/physics.
3. Read only the validation `.txt` file for the active validation step, if any.
4. Inspect the current code/tests that the next step will touch.
5. Rewrite `ACTIVE_STEP.md` so Claude has the necessary context without opening
   the long plan by default.
6. Update `ROADMAP.md` with the active milestone and next few milestones.
7. Keep `ARCHIVE_DONE.md` as the detailed completed-work log.

Do not rewrite the long `in-main-tex-you-will-quizzical-treasure.md` unless the
user explicitly asks. It is the full historical plan, not the startup context.

## When Claude Finishes A Coding Step

Claude is expected to mark `ACTIVE_STEP.md` completed and add
`Claude completion notes` plus `Implementation takeaways`. When Codex is called
after that:

1. Read those completion notes first.
2. Inspect the diff, current code/tests, and any test output the user provides.
3. Review Claude's implementation as a secondary correctness check before
   accepting the step as completed. Verify that the requested behavior,
   validation criteria, tests, guardrails, and documented caveats were actually
   satisfied; inspect relevant code paths and run or review tests as needed.
4. If the work is incomplete, incorrect, insufficiently tested, or fails an
   important guardrail, tell the user clearly what is missing or suspect. Do
   not archive the step as complete, do not rewrite the roadmap as if it is
   done, and do not generate the next `ACTIVE_STEP.md` until the gap is
   resolved or the user explicitly accepts the risk.
5. If the work passes Codex's review, move the durable summary and takeaways
   into `ARCHIVE_DONE.md`.
6. Update the relevant completed-step section in
   `in-main-tex-you-will-quizzical-treasure.md` with a compact
   `Implementation summary` of 20 lines or fewer. Summarize Claude's longer
   notes; do not paste raw logs, verbose diffs, or transient debugging detail.
7. Update `ROADMAP.md` so the completed step and next active step are clear.
8. Generate a fresh `ACTIVE_STEP.md` for the next coding step.
9. Preserve warnings about stale references, changed assumptions, input names,
   and validation caveats.
10. Remove transient detail from the new `ACTIVE_STEP.md`; keep only what Claude
   needs for the next coding step.

## Context Budget Rules

- Keep `CLAUDE.md` short and operational.
- Keep `ROADMAP.md` short enough to read every session.
- Put detailed completed history in `ARCHIVE_DONE.md`, not `ROADMAP.md`.
- Keep each completed-step summary in
  `in-main-tex-you-will-quizzical-treasure.md` to 20 lines or fewer.
- Put only the current coding packet in `ACTIVE_STEP.md`.
- Include equations, constants, file paths, validation criteria, and guardrails
  in `ACTIVE_STEP.md` when Claude needs them.
- Prefer summaries and exact pointers over copying large blocks from `main.tex`
  or the validation references.

## Project Guardrails

- Do not edit `ext/`, `bin/`, `obj/`, `build/`, `compile_commands.json`,
  `configure`, or `LICENSE`.
- Keep long physics/reference files on demand, not in the default context path.
- Preserve completed regression behavior unless the active step explicitly
  changes it.
- Use existing ALAMO patterns and tests as the implementation guide.
- Run test and validation simulations with 4 MPI ranks by default on this
  machine, unless the test metadata or user request explicitly requires a
  different rank count.
