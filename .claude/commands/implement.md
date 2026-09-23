---
description: Implement the active step. Reads ACTIVE_STEP.md and writes code, tests, and completion notes.
---

# /implement — Implementer Phase

You are the **implementer**. Your job is to code the active step described in
`ACTIVE_STEP.md` and produce a runnable, tested result.

## Step zero — before writing code

Verify the packet's **reachability ledger** and **key table**, and diff the
harness's printed commands against the table. **Stop on any failed row: that is a
planner bug, not yours.** Re-validate any frozen baseline whose keys the identity
set does not exercise by re-running its cheapest leg and comparing byte for byte.

**Smoke before campaign.** Run every leg briefly at the coarse mesh first: the
mechanism's witness diagnostic must be non-zero on test legs and zero on
controls. A grep proves reachability in principle; a smoke run proves it
happened, and catches what no grep can.

## Inputs you should read, in order

1. `CLAUDE.md` (already loaded — project rules, ALAMO architecture, commands,
   and do-not-touch list).
2. `ROADMAP.md` — confirm the active step matches what `ACTIVE_STEP.md` says.
3. `ACTIVE_STEP.md` — your task packet. This is the source of truth for
   scope, sources, goal, guardrails, and commands.
4. The specific source/test files named in `ACTIVE_STEP.md`. Read existing
   ALAMO patterns nearby and mirror them.

Open the long sources (`docs/project/in-main-tex-you-will-quizzical-treasure.md`,
`main.tex`, validation `.txt` files) **only** if the active task explicitly
calls for a constant or equation that is not already inlined in
`ACTIVE_STEP.md`. If you find yourself reaching for them often, the packet is
under-specified — stop and ask the user to run `/plan` again.

## What to produce

1. **Code** that satisfies every numbered deliverable in `ACTIVE_STEP.md`.
   - Mirror existing ALAMO patterns (integrator inheritance, ParmParse,
     `Numeric/` reuse, IC/BC `select_default`, microstructure paths, etc.).
   - Default-off any new behavior unless the step explicitly turns it on.
   - Do not touch `ext/`, `bin/`, `obj/`, `build/`, `compile_commands.json`,
     `configure`, or `LICENSE`.
   - Do not regress completed validation tests.

2. **Tests / regressions** as required by the step. Mirror existing
   `tests/MMWSpalling/<name>/` patterns: `input`, Python `test`, optional
   `reference/`. Use 4 MPI ranks by default unless the step says otherwise.

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
**not** archive yourself in `docs/project/ARCHIVE_DONE.md` or
`docs/project/in-main-tex-you-will-quizzical-treasure.md`. The planner moves the takeaways
later.

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

**Always hold the laptop awake for anything longer than a couple of minutes.**
Simulations here are not checkpointed, so a host sleep loses the run outright
and it has to restart from t = 0. Wrap every study launch in `caffeinate -dis`,
not just the ones whose packet happens to spell it out:

```bash
caffeinate -dis $PY <study>/run.py --cases ... --jobs N
# already launched without it? attach one to the running harness:
caffeinate -dis -w <harness pid> &
```

Prefer harnesses that assert this themselves (spawn
`caffeinate -dis -w <own pid>` at the top of the run loop), so it cannot be
forgotten at the call site.

Use the commands block in `ACTIVE_STEP.md` first. The standard fallbacks
from `CLAUDE.md`:

```bash
# Build (from repo root)
EIGEN=$PWD/ext \
  CPLUS_INCLUDE_PATH=/opt/homebrew/include \
  LIBRARY_PATH=/opt/homebrew/lib:/opt/homebrew/Cellar/gcc/15.2.0_1/lib/gcc/current \
  make -j8

# Run a testcase or validation simulation (4 ranks by default)
mpirun --oversubscribe --bind-to none -np 4 \
  bin/mmwspalling-3d-g++ tests/MMWSpalling/<step>/input

# Regression check
python3 tests/MMWSpalling/<step>/test tests/MMWSpalling/<step>/output
# Or, if the test needs numpy/yt on this machine:
/Users/tzetze20/Desktop/code/.venv/bin/python tests/MMWSpalling/<step>/test
```

## End of phase

Report:
- Files changed.
- Tests run and their pass/fail status.
- Any tests not run, with the reason.
- A pointer to the `## Claude completion notes` and
  `## Implementation takeaways` sections you just appended to
  `ACTIVE_STEP.md`.
