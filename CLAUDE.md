# CLAUDE.md

This file is the small entry point for Claude. Keep it short. The long project
plan remains in `docs/project/in-main-tex-you-will-quizzical-treasure.md`, but
Claude should not read that file by default.

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
`.claude/commands/<phase>.md` and contains the full phase instructions.

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

1. `ROADMAP.md` — **read `## Standing rules` first**: programme invariants that are
   never archived (unit conventions, scoring traps, coupled calibrations). Then current
   status, the active step, and next milestones.
2. `ACTIVE_STEP.md` — the detailed task packet for the active coding step.

Then use the project rules, commands, ALAMO gotchas, and guardrails below.

Open the long/reference files only when the active phase explicitly requires
them or the user asks:

- `docs/project/in-main-tex-you-will-quizzical-treasure.md` — full historical
  plan and checklist. The planner phase reads slices of this; the implementer
  should not need it.
- `main.tex` — full physics proposal and equations. Prefer the relevant
  excerpts already distilled into `ACTIVE_STEP.md`.
- `tests/MMWSpalling/validation/zhang/zhang_oglesby_mmw_validation_summary.txt`
  — Zhang/Oglesby MMW thermal validation reference.
- `tests/MMWSpalling/validation/hu/hu_thermal_spallation_validation_tests_augmented_vv_framework.txt`
  — Hu prescribed-temperature, thermoelastic, breakage, and LRST validation data.
- `tests/MMWSpalling/validation/kant/kant_von_rohr_validation_test_for_mmwspalling_plan.txt`
  — Kant/von Rohr flame-spallation threshold validation data.
- `docs/paper/paper.md` — ALAMO architecture paper; use only for ALAMO design
  background.
- `docs/project/ARCHIVE_DONE.md` — completed-step history. Read only for
  debugging regressions or when the planner is updating project context.
- `docs/project/AGENTS.md` — historical record of the prior Codex+Claude split
  workflow. Kept for reference; not active.

## Project Goal

Implement a unified thermomechanical millimetre-wave drilling model in ALAMO:
enthalpy heat conduction, Gaussian/Beer-Lambert MMW source, surface losses,
Voronoi mineral microstructure, thermoelastic stress, Drucker-Prager damage,
grain-boundary/cohesive effects, material removal by spallation or vaporisation,
and staged validation tests.

## ALAMO Architecture Notes

- Integrators in `src/Integrator/` are the top-level physics drivers. Mirror the
  existing multi-physics pattern in `src/Integrator/ThermoElastic.H`.
- Use virtual multiple inheritance for coupled integrators. Do not simplify the
  inheritance graph to avoid plumbing problems.
- Models in `src/Model/Solid/` and `src/Model/Interface/` are templated
  constitutive classes used by mechanics integrators.
- Operators in `src/Operator/` wrap AMReX linear operators. Mechanics uses
  `Operator::Elastic`; diffusion-like solves use the implicit/cell operators.
- `Solver::Nonlocal::Linear` wraps `amrex::MLMG`.
- IC/BC classes use `pp.select_default<Type1, Type2, ...>(...)` and the
  `IO::ParmParse` family of macros.
- Reuse `src/Numeric/` utilities such as `Stencil`, `Gradient`, `Laplacian`, and
  `Interpolate::CellToNodeAverage` instead of hand-rolling stencils.
- Tests live under `tests/<Name>/` with an `input`, a Python `test`, and optional
  `reference/` data. Mirror existing `tests/MMWSpalling/*` patterns.
- The build is header-driven: `.H` files under `src/` are auto-discovered; each
  top-level `.cc` under `src/` becomes a binary in `bin/`.

## Lessons Learned

1. Mirror existing integrator inheritance exactly. Earlier attempts to simplify
   base classes caused hard-to-debug vtable/thunk crashes during AMR setup.
2. Reuse `HeatConduction::temp_mf` as the temperature field. It already has BC,
   ghost-cell, FillPatch, and refinement plumbing.
3. `Parse` is static and `pp.queryclass<Foo>("prefix", value)` strict-checks
   unused keys. When an input fails, inspect each base class `Parse` before
   guessing at parameter names.
4. The Zhang/Oglesby validation reference file corrects an old beam-waist
   mistake: use `omega_0 = 0.02 m` for that case, not the stale `0.20 m` note in
   older plan text.
5. The implicit thermal path currently targets material thermal validation and
   single-level runs. Revisit compatibility before coupling it to later fracture
   mechanics steps.

## Commands

```bash
# Build from /Users/tzetze20/amr_tools/alamo.
# Required env on this machine: eigen3 lives in ext/, libpng in homebrew,
# gfortran from homebrew gcc.
EIGEN=$PWD/ext \
  CPLUS_INCLUDE_PATH=/opt/homebrew/include \
  LIBRARY_PATH=/opt/homebrew/lib:/opt/homebrew/Cellar/gcc/15.2.0_1/lib/gcc/current \
  make -j8

# Run a testcase or validation simulation.
# Use 4 MPI ranks by default on this machine unless directed otherwise.
mpirun --oversubscribe --bind-to none -np 4 \
  bin/mmwspalling-3d-g++ tests/MMWSpalling/<step>/input

# Regression check.
python3 tests/MMWSpalling/<step>/test tests/MMWSpalling/<step>/output
# Or, if the test needs numpy/yt on this machine:
/Users/tzetze20/Desktop/code/.venv/bin/python tests/MMWSpalling/<step>/test

# Remote Control (drive phases from the iPhone app / claude.ai/code).
# Server rooted in this repo; -c resumes the last session. Header has setup.
docs/project/tools/rc.sh
```

## Do Not Touch

- `ext/` — vendored AMReX and Eigen sources.
- `bin/`, `obj/`, `build/` — build artefacts managed by `make`.
- `compile_commands.json` — auto-generated.
- `configure` — generated configuration script.
- `LICENSE` — leave as-is.
- Existing unrelated integrators unless the active step explicitly requires a
  small, justified integration point.

## Context Budget Rules

- Keep `CLAUDE.md` short and operational.
- **`ROADMAP.md`: budget the hand-written sections, not the total.** The file is read at
  startup by every phase, so the cost is context; the danger is the notes turning into a
  second archive. A single total cap measured the wrong thing: the **generated**
  `## Current Status` list grows one line per archived step and carries no decisions, so a
  total cap shrank the hand-written part every time a step was archived. Rule (set
  2026-10-07, replacing the 300 → 400 total cap):
  - **Generated list excluded.** The block from the `## Current Status` heading to the line
    before `Active step:` is written by `docs/project/tools/roadmap_status.py` and is not
    counted. Never hand-edit it.
  - **Hand-written content ≤ 350 lines.** Everything else in the file. Check:
    `awk '/^## Current Status/{g=1} /^Active step:/{g=0} !g{n++} END{print n}' ROADMAP.md`
    (287 on 2026-10-07).
  - **Section budgets inside that:** `## Standing rules` **10–20 bullets, never archived** (raised 15 → 20 by the user on 2026-10-12:
    the block had sat at 15/15 while two new programme invariants earned a slot)
    (unit conventions, scoring traps, coupled calibrations); `## Known Stale/Important Notes`
    **≤ 40 bullets of ≤ 3 lines**, each a pointer to the evidence; `Active step:` plus
    `Next milestones` together **≤ 100 lines**; `## Validation Story` ≤ 20 lines.
  - **Total alarm at 600 lines** (`wc -l`), so a runaway is still caught by a one-line check.
  - ⚠ A line cap is gamed by long lines. If lines run past ~200 characters, count words
    instead: the hand-written part should stay under **5,000 words**.
  - ⚠ The cap is what matters: without it the standing block becomes the new notes section
    and the notes section becomes the new archive. A check on a label ("exactly one
    `Active step:`") passes while the section under it grows (log 2026-10-02 §1).
- **`ACTIVE_STEP.md` under 500 lines**, carrying only what the implementer executes.
  Argument and history go to `Claude_markdowns/<date>.md` and `ARCHIVE_DONE.md`.
  ⚠ Keep the number that is the *specification*; move the number that is *provenance*.
  - ⚠ Raised 300 -> 500 by the user on 2026-10-07: compressing the packet was costing more
    planning time than the bloat cost, and the cap also has to hold the implementer's
    `## Claude completion notes` and the reviewer's `## Review findings`, which append to the
    same file and are only archived at the next `/plan`. 300 bound the *planner's* packet
    alone and was breached by the phases that follow it (P2-2 reached 444).
- **`ARCHIVE_DONE.md` carries a heading index at the top**, regenerated on append.
  Pointers into a 14,000-line file with no index are worse than the bloat they replace.
- Put detailed completed history in `docs/project/ARCHIVE_DONE.md`, not
  `ROADMAP.md`.
- Keep each completed-step summary in
  `docs/project/in-main-tex-you-will-quizzical-treasure.md` to 20 lines or fewer.
- Put only the current coding packet in `ACTIVE_STEP.md`.
- Include equations, constants, file paths, validation criteria, and guardrails
  in `ACTIVE_STEP.md` when the implementer needs them.
- Prefer summaries and exact pointers over copying large blocks from `main.tex`
  or the validation references.
