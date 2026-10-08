# `unit/support_rule`, retired by P2-2a (2026-10-07)

`spall.support_rule` and its two companion keys were **retired on supersession**.
This directory holds the unit test that exercised them, **moved and not deleted**,
because `../MANIFEST.md` hashes its eight legs by path (standing rule 9: never
`rm` what a manifest cites).

## What the rule was

It admitted a column with >= 2 exposed sides that the ordinary criterion
refused, because a thin tower heated from all sides has no cold interior and so
no skin to shed. It borrowed the criterion's own threshold but pinned the drive
to the self-similar `1 - e^{-c_f} = 0.3297` through `drive_frac_force`.

## Why it went, and why the guard could not simply be kept

Its **80 um** flaw was matched to that reduced drive **and to no other**: against
full `absolute` stress the same flaw fires about **240 K lower** (525 K against
the criterion's 768 K). A parse guard refused `support_rule = 1` under
`driving_stress = absolute` for exactly that reason.

⚠⚠ That guard was a **calibration** guard, not a code dependency.
`drive_frac_force` overrode the fraction *before* the reference-depth branch and
worked under either driving form. So deleting `differential` and the guard would
have **broken nothing and aborted nothing** — `support_rule = 1` would have kept
parsing and started firing 240 K early, **in silence**. Retiring the rule
outright was the alternative to re-fitting 80 um to full stress, which would
have been a new calibration.

## The ground: supersession, not failure

A column with exposed sides is what a **side-face** criterion attacks, with
physics rather than a neighbour count, and lateral removal **ships**
(`spall.lateral_drive / lateral_removal / lateral_subcell / lateral_thermal`).
P1e's staircase correction removes what *makes* towers but has no removal path;
the lateral criterion is the removal path.

⚠ **For a future planner:** there is now no geometric removal path for a
free-standing hot column, **by design**. If towers recur the answer is the
lateral criterion, not a revived neighbour count.

## What is frozen here

These legs **cannot be regenerated** — the keys no longer exist and the binary
that produced them is gone (the toolchain is not byte-reproducible). The outputs
and their hashes in `../MANIFEST.md` are the only record.

Three analyses read `../output_fh` and are likewise **frozen and no longer
regenerable**:

- `studies/p1c_isotropy/analyze.py`
- `studies/p1c_isotropy/roundness.py`
- `studies/p1c_isotropy/attribution.py`
- `studies/p1e_staircase/analyze.py`

⚠ `test` and `input` in this directory are kept **as artefacts**. `input` still
sets `spall.driving_stress = differential`, and the legs still pass
`spall.support_rule` on the command line, so running either against a current
binary leaves those keys **unconsumed** — which `unit/lateral_flux_lock`
asserts deliberately, and which is the evidence that the keys are gone rather
than parsed and ignored.
