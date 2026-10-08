# D2z-2 PRE-FLIGHT — item 0, read-only, before any leg

## Reachability ledger — 12 of 12 cited line numbers EXACT

Unlike D2z (where D2z-0/D2z-1 edits had drifted every citation by ~+11), every
site the packet names is exact in the working tree:

| site | packet | verified |
|---|---|---|
| `PairDelta` | `:6477` | ✓ `static Set::Scalar PairDelta(` |
| `StepCellWeight` | `:6508` | ✓ |
| `AreaOfDelta` | `:6518` | ✓ |
| `StepConductance` | `:6527` | ✓ |
| `side_face_area_law` default | `:6819` | ✓ `query_default(..., std::string("step"))` |
| `side_face_factor` default | `:6748` | ✓ `query_default(..., 1.0)` |
| `delta_sum` accumulation | `:3315` | ✓ `delta_sum += dsum;` |
| kernel `g` | `:1687` | ✓ `kf *= StepConductance(...)` |
| `WallStepDiagnostic` | `:2109` | ✓ |
| diagnostic `g` | `:2178` | ✓ `g = StepConductance(...)` |
| `AssertStepSplit` | `:2057` | ✓ |
| `AssertStepSplit` call | `:1465` | ✓ |

## Deliverable 4 — the derivation, reproduced before it was implemented

```
       a  n  delta  D2z form  derived   packet
    1.29  2    2.3     1.31    1.026     1.03
    0.65  2    2.3     1.62    1.066     1.07
    0.65  3    3.2     1.83    1.020     1.02
    0.65  2    1.8     1.44    0.948     0.95
```

All four rows reproduce. Three further properties, each checked independently and
now asserted in `unit/side_face_step_height/test` (checks j2–j6):

* **`n = 1` is exact** — `max |derived − D2z| = 0.0` over a 7 × 8 grid of
  `(a, delta)`. Nothing D2z measured inside r < 25 mm can move.
* **`a → 0` degenerates to `delta / n`**, not `min(delta, 1)`. The `a ≈ 0` branch
  is fixed accordingly.
* **The never-fired fallback is exactly 1 for any `a`** (`phi = 1` on both
  columns ⇒ `delta = n`), with no special case. The D2z form returned 1.67 at
  `n = 3` and 2.90 at `n = 7`.

**Reconciliation of the packet's two amplification figures.** The packet quotes a
bound of "1.85x at 2 mm and 3.09x at 1 mm" and "1.67 at n = 3, 2.90 at n = 7",
while D2z's own PREFLIGHT derived `a = 1.29` at 2 mm, where `1/(1 − e^{−a})` is
1.380. There is no contradiction: each figure is evaluated at its own mesh's
measured `a` (1.85 ⇒ `a = 0.778`; 3.09 ⇒ `a = 0.391`), and `a = dz²/(dt·alpha)`
varies with the local flip interval. Both reproduce exactly at the `a` they
imply. Recorded so the disagreement is not rediscovered as an error.

## ⚠ `g` is NOT bounded by 1 — and the packet's framing needs this correction

`delta = n + phi_S − phi_T` with both phases in [0, 1], so `delta` spans
`[n − 1, n + 1]` and `g > 1` whenever `phi_S > phi_T`. That is physical: the
shorter column is due to fire and the taller has just fired, so the riser really
is momentarily taller than the `n` cells the grid resolved.

The derivation's bound is `g ≤ g(n+1, n, a) ≤ (n+1)/n`:

| n | derived bound at a = 0.391 | D2z bound |
|---|---|---|
| 1 | 1.676 | 1.676 |
| 2 | 1.273 | 2.134 |
| 3 | 1.145 | 2.443 |
| 7 | 1.022 | 2.955 |

**The derived bound FALLS with `n` where the D2z bound ROSE.** That, not "g ≤ 1",
is the correct statement of the repair, and it is what criterion (j) records.
An early draft of the unit test asserted `g ≤ 1` and failed against correct code.

## ⚠ Criterion (d) compares two DIFFERENT means

`wall_step_q_eff_in / wall_step_q_in` is **flux-weighted** (`Σ q g / Σ q`);
`step_g_mean` is **face-weighted** (`Σ g / N`). They coincide only when `g` and
`q` are uncorrelated, and they are strongly correlated. `q_raw` is also
**negative** wherever the neighbour is hotter, so the flux-weighted mean is not
even bracketed by the range of `g`. Measured on the unit test's 4 × 4 grid:
flux-weighted 0.132–1.015 against face-weighted 0.333–1.020, **max difference
0.383**.

This is recorded in `CRITERION.md` before launch. (d) is scored **exactly as
written**; the exact identity `wall_step_q_in == wall_step_q` with the gate off
is reported beside it as the check that actually has teeth.

## Binaries — fingerprinted, because the comparison IS the binary

| form | sha256 | legs |
|---|---|---|
| derived | `2dde94489142ce84597102a171ffcfde1f88d41324783250c49437a97bda75b6` | `R2d` `R1d` `SW` `SWs` `S0` |
| D2z | `f9085c65ca28959aa650d2248a0f9b2277afedd93bf56032ab2846e42c15ef5b` | `R2` `R2w` `SW0` |

The D2z-form binary was produced by temporarily restoring the old normaliser,
building, copying the binary to the session scratchpad, and **restoring the
source immediately** (verified: `grep "SW0 BUILD ONLY"` returns 0 lines). No
second conductance form is alive behind a switch. `run.py` re-hashes the binary
before every leg and **refuses to launch** on a mismatch — a leg run against the
wrong binary is not a wrong number, it is an unattributable one.

Both builds happened **before the first leg**, so no rebuild occurs mid-sweep.

## Deliverable 1 — a decision the packet's wording does not settle

The packet says `delta_sum` should accumulate `PairDelta` **"once per stepped
face"**. The branch runs for *every exposed cell* of a riser, so a per-face
accumulation weights an `n`-cell pair `n` times, while `delta_offline.py` weights
every **pair** equally — and criterion (b) compares the two directly to 10 %.
Accumulating per face would reintroduce the population mismatch this packet
exists to remove.

**Implemented once per PAIR**, at the taller column's top exposed cell
(`k == kT`). The per-face count is preserved as the new `delta_faces` column, so
nothing is lost and `delta_faces / delta_n` is the mean riser height in cells.
Recorded here and in `CRITERION.md` because it is a deviation from the packet's
literal wording in favour of its stated purpose.

## Smoke — the witnesses before the campaign is trusted

`unit/side_face_step_height` (30 checks, all pass) exercises every new path on a
4 × 4 × 60 grid before any study leg:

* the eight new columns are **absent** with the keys off and **present** with
  them on (k, k4c: 38 shared columns bit-identical between the gate pair);
* `step_diag_radius` **removes faces** (n_in 4 vs n 24) and, with the key unset,
  the `_in` accumulators are **bit-identical** to the ungated ones;
* the gate is **diagnostics-only** — `GATE` and `GATE0` agree byte for byte on
  the removal log and all 12 `Cell_D` files;
* `delta_mean` is now unbounded above (0.0000–10.3833) where
  `StepCellWeight`'s per-cell share was pinned near 1.
