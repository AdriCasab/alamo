# D2z-0 CRITERION — hashed before the first scored leg

`spall.removal_increment = full_cell` on the frozen arbiter pair. Written from
`ACTIVE_STEP.md`, with every reference value reproduced independently first
(see `PREFLIGHT.md`). Hashed with **no scored-leg output on disk**; the smoke
runs that preceded it were deleted, and their result is recorded in
`PREFLIGHT.md`.

Metric: `col_rate` and the metadata reader are **imported unchanged** from
`studies/d2j0b_plane/analyze.py`, the same estimator the frozen `A_f02`
baseline was scored with. Frozen studies are read by path and never re-run.

Legs: `H2` (A_f02's argv, key absent), `F2` (+ full_cell, 2 mm),
`F1` (+ full_cell, 1 mm). Blocks are 50–100 / 100–150 / 150–200 s throughout.

---

## (a) SCORED — the rate is essentially unchanged

Disk-mean rate over **r < 25 mm**, `F2` against frozen `A_f02` and `F1`
against frozen `A_f02_1mm`, per block.

**Bar: the change lies in +0.2 % to +0.6 % at BOTH meshes**, point prediction
**+0.36 %** — one timestep of delay per cell removed.

⚠ **NOT r < 30 mm.** Reproduced here: the same frozen pair reads
**−0.2 / +3.0 / −1.8 %** at r < 25 mm and **−6.4 / −8.0 / −14.0 %** at
r < 30 mm, because that disc includes the 25–30 mm rim stall. A 5 % bar on
r < 30 mm is failed by the configuration it exists to protect.

Outside the band, the mode is doing something other than removing the delay.
Report the r < 30 mm figures alongside so the rim is never silently scored.

## (b) SCORED — the mesh gap does not move

`F1` vs `F2` at r < 25 mm, against the `f = 0.2` reference
**−0.2 / +3.0 / −1.8 %** (reproduced in `PREFLIGHT.md`).

**Bar: the change is < 0.2 percentage points in every block.**

## (c) SCORED — plumbing

1. `firings_per_void` collapses to **exactly 1.000** on `F2` and `F1`.
2. Exactly **one `removal_events.csv` regime-1 row per voided cell**, i.e.
   `n_voided == 1` on every such row.
3. `removal_increment_cells` reads **exactly 1.0** on every firing row.

⚠ **Baseline value.** `ACTIVE_STEP.md` states the baseline as **1.977**. That
is **Meier `A`'s** value; reproduced here it is 1.977 on Meier `A` but
**1.957 on the arbiter `A_f02`** and 1.962 on `A_f02_1mm`. The arbiter is what
this packet runs, so the honest statement of (c) is **1.957 → 1.000** at 2 mm
and **1.962 → 1.000** at 1 mm. The bar is on the destination (exactly 1.000),
which is unaffected. The disagreement is reported, not substituted.

## (d) SCORED — identity

1. `H2` is **byte-identical** to frozen `d2l_scan/output/A_f02` — `thermo.dat`,
   `_removal_events.csv`, and every `Level_0/Cell_D_*` — over the interval both
   completed. Both run to `stop_time = 250.0`, so the interval is the whole run;
   if either truncates, compare only the common prefix and say so.
2. Tree identity **2639/2639** with the key unset.

⚠ The 2639-file gate **cannot** stand in for (d1): no input in that set turns
spalling on with an event log at the arbiter's configuration, so `A_f02` is not
covered by it. `H2` is what re-validates the frozen baseline against this
binary.

---

## Pre-registered predictions — contiguous, no dead zones

* **P1** rate change in **+0.2…+0.6 %** at both meshes (point +0.36 %) held;
  outside refuted.
* **P2** mesh-gap change **< 0.2 pts** in every block held; larger refuted.
* **P3** `firings_per_void` → **exactly 1.000** held; anything else refuted.
* **P4** `H2` byte-identical to frozen `A_f02` held; any difference refuted.

**Refuting evidence already on disk:** the frozen `A_f02` pair. If `full_cell`
moves the rate by more than ~1 %, the mode is not merely removing a one-dt
delay — check `n_voided` and the phi tie **first**, before reaching for any
physical explanation.

## What would make this inert

| failure | witness | bar |
|---|---|---|
| key never reaches the scan | `firings_per_void` | 1.957 → 1.000 |
| phi tie misses or doubles | in-code assertion + `n_voided` | == 1 on every firing |
| default perturbed | `H2` vs frozen `A_f02` | byte-identical |
| new columns break frozen headers | tree sweep | 2639/2639 with the key unset |
