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

⚠⚠ **Step zero is the project's only outside read of the packet.** The planner's own
re-read is a self-check and cannot catch a planner blind spot. If a ledger row does not
resolve against the code — including the **conjunction** of conditions a deliverable
requires — **stop and report it. That is a planner bug, not yours** (log 2026-10-02 §1).

**Run order and gating.** Run the packet's `## Criteria` rows **in its stated run
order**. After a **gating** row fails, launch **no further long legs or sweeps**: report
that row as the result. Diagnostics that fall out of runs already made are still
reported.

**Tests per edit, and the one gate.** Pick the two or three tests covering the sites in
the packet's ledger from the table in
`tests/MMWSpalling/studies/tree_sweep_0922/witness/WHICH_TESTS.md`
— that is a ~2-minute check per build. The **full tree gate runs ONCE
per packet, on the binary that ships.** Its claim is a property of that binary: with
every new key off, it reproduces the frozen tree. ⚠⚠ **A second full gate needs a stated
discriminating reason, and "attribution" is not one** — bisecting with targeted tests
after a failure is cheaper than pre-paying a 10-minute sweep before it. ⚠ The harness
memoises by binary sha; **do not pass `--force`** unless the binary is unchanged and you
mean to re-run anyway.

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
   - ⚠ **Any change to a file the identity gate hashes — including an added event-log
     column — is either keyed off by default or declared as a re-pin request in the
     completion notes. Never both silently.**
   - ⚠⚠ **Every new removal regime owes its own ledger term.** An omitted term takes a
     cell's enthalpy out of the books unrecorded; the energy invariant is what finds it,
     and the per-cell arithmetic is what turns the suspicion into a diagnosis.
   - ⚠⚠ **A fix claimed to be bookkeeping-only must be PROVEN so**, not asserted: keep
     the pre-fix run and show the rebuilt binary reproduces its event log byte for byte.
     On a toolchain that is not byte-reproducible, "it only touches a reported scalar" is
     an argument, not evidence.
   - ⚠⚠ **`--check` cannot re-validate an anchor after a source change** — it re-reduces
     output written by the old binary. Re-run the cheapest leg and compare bytes.
   - ⚠ **Archive a harness's output directory before re-running it** if the harness
     `rmtree`s its own legs and their hashes live in an earlier manifest; on a
     non-byte-reproducible toolchain that provenance is otherwise unrecoverable.
   - ⚠ **Do not demand bit-exactness from a text diagnostic.** `thermo.dat` is written
     with `%g` — six significant digits — so compare the printed token to
     `"%g" % want`, not the full double.
   - ⚠ **Non-vacuity:** an assertion that probes a quantity which is zero **by
     construction** on both legs will "match" two zeros and prove nothing. Probe
     something the mechanism actually moves.
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
   - Append `## Claude completion notes`: files changed, and **the packet's
     `## Criteria` table filled in** — measured value, pass/fail, and the output path,
     one row per criterion. **No free-form prose for what the table already carries.**
     Tests not run are rows with the reason in the measured column.
   - ⚠ **Do not quote anything the packet marked unquotable.**
   - Append `## Implementation takeaways` with **only what the diff and the outputs
     cannot show**: gotchas, design choices, input parameter names, changed assumptions,
     caveats, and anything a reader would otherwise have to re-derive. ⚠ Do not restate
     what the filled criteria table already carries. Write these for the reviewer and
     the next planner — they are the durable record of this step.

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
