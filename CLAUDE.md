# CLAUDE.md

This file is the small entry point for Claude. Keep it short. The long project
plan remains in `in-main-tex-you-will-quizzical-treasure.md`, but Claude should
not read that file by default.

## Context Loading Order

After opening this file, read only these additional files by default:

1. `ROADMAP.md` - current status, active milestone, and next milestones.
2. `ACTIVE_STEP.md` - the detailed task packet for the next coding step.

Then use the project rules, commands, ALAMO gotchas, and guardrails below.

Open the long/reference files only when the active task explicitly requires
them or the user asks:

- `in-main-tex-you-will-quizzical-treasure.md` - full historical plan and
  checklist. Treat as the long source of truth, not default context.
- `main.tex` - full physics proposal and equations. Prefer the relevant
  excerpts already distilled into `ACTIVE_STEP.md`.
- `zhang_oglesby_mmw_validation_summary.txt` - Zhang/Oglesby MMW thermal
  validation reference.
- `hu_thermal_spallation_validation_tests_augmented_vv_framework.txt` - Hu
  prescribed-temperature, thermoelastic, breakage, and LRST validation data.
- `kant_von_rohr_validation_test_for_mmwspalling_plan.txt` - Kant/von Rohr
  flame-spallation threshold validation data.
- `docs/paper/paper.md` - ALAMO architecture paper; use only for ALAMO design
  background.

## Project Goal

Implement a unified thermomechanical millimetre-wave drilling model in ALAMO:
enthalpy heat conduction, Gaussian/Beer-Lambert MMW source, surface losses,
Voronoi mineral microstructure, thermoelastic stress, Drucker-Prager damage,
grain-boundary/cohesive effects, material removal by spallation or vaporisation,
and staged validation tests.

Claude does the coding. Codex periodically rewrites `ACTIVE_STEP.md`, updates
`ROADMAP.md`, and archives completed work so Claude does not need to ingest the
full project history on every session.

## End-of-Task Handoff

When you complete the active coding step, update `ACTIVE_STEP.md` before ending
the session:

- Mark the step status as completed.
- Add a short `Claude completion notes` section with the files changed, tests
  run, test result, and any tests not run.
- Add a short `Implementation takeaways` section with gotchas, design choices,
  input names, changed assumptions, or caveats that Codex should preserve in
  `ARCHIVE_DONE.md`, `ROADMAP.md`, or the next active task.
- Do not rewrite the next step yourself unless the user asks. Codex will read
  those notes later and regenerate the planning files.

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
```

## Do Not Touch

- `ext/` - vendored AMReX and Eigen sources.
- `bin/`, `obj/`, `build/` - build artefacts managed by `make`.
- `compile_commands.json` - auto-generated.
- `configure` - generated configuration script.
- `LICENSE` - leave as-is.
- Existing unrelated integrators unless the active step explicitly requires a
  small, justified integration point.
