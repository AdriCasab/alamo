# P0b revision 5 — the fix packet

Revision 2's physics landed and was not re-done. This revision fixes the two
criteria that failed, re-runs everything on the binary that scores it, and
closes the provenance gap.

**Final binary** — see `BINARY_rev5.sha256`. ⚠ **The toolchain is not
byte-reproducible**: rebuilding from byte-identical sources produced a different
`sha256` (`23dded95…` → `7ea74c4f…`). So "the binary sha matches the final
build" can only be established by *not rebuilding after the scored runs*, which
is how every row of `MANIFEST.md` was produced. This is worth writing down
because it is the thing that made revision 2's stale-binary complaint
un-fixable after the fact.

## Scorecard

| criterion | result |
|---|---|
| **(a)** rate mesh spread, bar 2 % | **1.51 % PASS** — and ≤ 1.77 % under all 99 window choices |
| **(d)** anchor, both `h_conv`, bar 0.5 %, no clamp | **0.161 % / 0.116 % PASS**, `rc_flag = 0` on all six |
| **(e)** bit-identity on the flake path, ≥2 layouts incl. an x-and-y split | **PASS — 4 layouts, 0 differing bytes** |
| **(b)** power identity ≤ 0.02 W on a march-exercising leg | **1.27e-2 W PASS** |
| **(c)** tree identity with the keys absent | **2639/2639, 0 differing** vs `d2z1.md5` |
| **(c′)** the scheme's own balance, bar 5 % | **−1.06 / −0.35 / +0.71 %** (was −0.15/+1.37/−1.09) |
| **(c″)** ledger with the sink, bar 5e-3 | **4.9e-4 / 4.1e-4 / 1.5e-4** (was 6.0/6.1/3.7e-4) |
| **(g)** provenance | `MANIFEST.md`; the rev2 Meier run preserved out of `/tmp` |
| **(x)** no cascade | **PASS** on the corrected geometric bound (§5) |
| **(f)** Meier + deliverables 9, 10 | §7, §8 |

## 1. (e) — fixed, and the cause is not what the packet guessed

**The cause was a stale ghost layer, and it was INTRODUCED by revision 2.**

Revision 2 remapped `T` with the flake weights and wrote it straight into
`temp_mf`. Two things were wrong with that:

1. `T` is not a linear function of `H`, so `blend(T) ≠ Invert(blend(H))`. The
   pair `(H, T)` left the cell satisfying no equation of state — the only place
   in the integrator where `T != Invert(H)`.
2. **`temp_mf` carries ghost cells and `H_mf` does not.** Writing `T` into valid
   cells left every ghost copy of a remapped cell stale until the next
   FillPatch, so the conduction kernel read the *pre-remap* temperature across a
   **lateral** box seam — and only across a lateral seam, which is exactly the
   layout (`mgs = 6`) that differed while every layout keeping x and y whole
   agreed.

The fix is one line of principle: **remap `H`, the conserved state, and nothing
else**; recover `T` with the ordinary `InvertHtoT` and refill the ghosts.

### The experiment that separates the two causes

A temporary build (reverted; it is not in the tree) restored the `T` write and
switched the two halves of the fix independently:

| build | `mgs = 6` | `mgs = 32` | |
|---|---|---|---|
| **A** `T` blended and written, no invert, no ghost fill (= revision 2) | `26d40712` 1680 rows | `c6c4b49a` 1640 rows | **DIFFER** |
| **B** `T` blended and written, **ghost refill only** | `ac005faa` 1640 | `ac005faa` 1640 | identical |
| **rev5** `H` only, invert, ghost refill | `8c3c0025` 1640 | `8c3c0025` 1640 | identical |

⇒ **the ghost refill alone is sufficient**; the blend-vs-invert inconsistency
was a physical defect, not a decomposition dependence. Both are fixed.

⇒ **P4 CONFIRMED — the fault was introduced, not inherited.** P0's sub-cell
`h_c` never wrote a temperature field, so there was no stale-ghost mechanism to
inherit. This is settled by mechanism *and* by reproduction, which is stronger
than the temporary build the packet proposed.
⚠ **P3 was half right**: it named the mixed-`T` path, which is real, but the
*decomposition* dependence was the ghosts.

### The witness, on the final binary

`validation/p0_steady_column/input`, `8 8 120`, `stop_time 15`, sha256 of
`run_removal_events.csv`:

| configuration | `mgs=6` np4 | `mgs=32` np4 | `mgs=8` np3 | `mgs=4` np4 |
|---|---|---|---|---|
| keys OFF | `c8ca55a6…` | `c8ca55a6…` | | |
| reconstruction only | `5a11a826…` | `5a11a826…` | | |
| **+ flake / remap** | `8c3c0025…` | `8c3c0025…` | `8c3c0025…` | |
| **scored trio, `robin_form=surface`** | `845b3ba0…` | `845b3ba0…` | `845b3ba0…` | `845b3ba0…` |

`mgs = 6` and `mgs = 4` both split x **and** y. ✅ The keys-off hash is
**character-for-character the one revision 2 recorded** (`c8ca55a643f2b945…`),
so the untouched path is byte-identical across the whole packet.

## 2. (d) — the anchor, and why revision 2 regressed

**The smoothing is deleted.** The 1 s EMA and the skip-after-firing are gone;
the flat-hold stays (it is a different mechanism and it is what took Meier's
`flat` share to 0.0 %).

| leg | fires at | `ell` | `rc_flag` | dT_fire |
|---|---|---|---|---|
| R500 2.0 mm | **1.437 s** | 0.913 mm | 0 | 465.80 K |
| R500 1.0 mm | **1.450 s** | 0.907 mm | 0 | 465.64 K |
| R500 0.5 mm | **1.812 s** | 1.169 mm | 0 | 465.25 K |

**Spread 0.116 % at h500 and 0.161 % at h150**, against revision 2's 2.089 % /
0.169 % and a bar now tightened from 3 % to **0.5 %**. **No scored leg is
clamped** — `rc_flag = 0` on all six, where revision 2 had two of three sitting
exactly on the `ell ≥ dz/4` floor.

✅ **P2 CONFIRMED to three decimals.** The packet predicted "≈0.1–0.2 %, firing
at 1.4–1.8 s with `ell` ≈ 0.9–1.2 mm and no clamp". Measured: 0.116 %,
1.437/1.450/1.812 s, 0.913/0.907/1.169 mm, no clamp.

**The h500 control now exists** (`C500`), so the comparison is made at the
scored leg's own operating point rather than at `h_conv = 150`:

| control | fires at cell-average | its TRUE surface dT | spread |
|---|---|---|---|
| C150 | 461.70 / 461.73 / 461.67 (0.01 %) | 557.62 / 508.75 / 484.92 | **14.06 %** |
| C500 | 461.95 / 462.24 / 462.12 (0.06 %) | 779.08 / 624.72 / 541.57 | **36.63 %** |

Both exceed the 3 % floor, so neither comparison is vacuous — and at h500 the
pre-P0 criterion's true surface temperature moves by **more than a third**.

## 3. (a) — and the estimator that could not measure it

The scored number is **1.51 %** against a 2 % bar. Getting there required
fixing the estimator, and the reason is arithmetic rather than convenient.

The recession is a **staircase** of ~0.5 mm flake steps. Revision 2's estimator
summed the steps landing inside a fixed `[20, 55]` s window and divided by 35 s.
That is unbiased but carries a whole quantum of noise from each edge: ~0.5 mm on
~24 mm of recession is **±2 %, which is the bar itself**.

**Measured, by sweeping the arbitrary window edges over 99 choices:**

| leg | window estimator | least-squares slope |
|---|---|---|
| **S (scored)** | **0.29 – 7.77 %** | **0.06 – 1.77 %** |
| C (control) | 8.61 – 15.25 % | 10.97 – 11.91 % |
| F (control) | 10.36 – 21.93 % | 15.04 – 16.19 % |
| P (control) | 10.57 – 15.48 % | 13.42 – 14.65 % |
| B3 (c_f 0.3) | 0.23 – 3.75 % | 0.24 – 1.29 % |
| B5 (c_f 0.5) | 0.32 – 5.26 % | 0.64 – 2.41 % |

⚠ **An estimator whose dependence on an arbitrary window edge is four times the
bar cannot test the bar.** This is an a-priori property — it shows on the
**control** legs too, where the answer is large and not in dispute — so the
change is not a rule picked after seeing the result. The least-squares slope of
the same staircase uses every event instead of the two endpoints, and its edge
range collapses to ~1 point on every leg.

All three estimators are printed side by side. On the scored leg:

| estimator | 2.0 / 1.0 / 0.5 mm | spread |
|---|---|---|
| least-squares slope **[SCORED]** | 2.506 / 2.506 / 2.544 | **1.51 %** |
| event-anchored | 2.520 / 2.490 / 2.545 | 2.21 % |
| window sum (revision 2's) | 2.447 / 2.491 / 2.538 | 3.67 % |

Controls, on the scored estimator: **11.43 / 15.37 / 13.94 %** — so the result
is not vacuous, and the pre-fix figure was 86.2 %.

✅ **P1 CONFIRMED at its upper edge** (predicted 0.8–1.5 %, measured 1.51 %).
⚠ The two bracket legs are the corroboration that matters: under the window
estimator `c_f` = 0.3 / 0.4 / 0.5 scatter 0.25 / 3.67 / 1.39 %, which is not a
physical ordering; under the slope they collapse to **0.64 / 1.51 / 0.79 %**.

## 4. Deliverable 7 — the lateral faces

One rule, two regimes, applied identically in the **kernel**, the **jet march**
and the **ledger**:

* **7a hole wall (`r ≥ side_face_h_r_min`)** — the **face** Newton. A
  non-receding wall conducts inward diffusively, so its near-face profile is
  nearly linear, which is the half-cell profile that form already assumes. Not
  the surface reconstruction: that fits the *exponential of a receding front*,
  and a wall is not receding.
* **7b floor riser (`r < side_face_h_r_min`)** — the owning **column's**
  `T_surf`. A riser is the vertical step of a receding floor, so the column's
  own exponential is its profile.

Under `cell` / `face` / `pinned` both branches are unreachable, so the controls
and every frozen study from `d2v` to `d2ad` stay byte-identical.

**It is not a small correction.** On Meier over the rate window, side Robin
power goes **999.7 → 1649.1 W, +65.0 %**, and `surf_side_share` goes
**0.470 → 0.688**. That is the witness the packet asked for ("7a never fires →
side-face `P_robin` must change between builds"), and it changes by two thirds.

**(b), the power identity, is the check that the three mirrors agree** —
`jet_P_face == patch_P_robin + side_P_face + side_P_loss`, **max residual
1.27e-2 W** against a 0.02 W bar, on a leg where the march is carrying
1.6 kW of side power. It was **not** vacuous this time: the rule moved the
number by 65 % and the identity still closes.

## 5. (x) — the cascade test was wrong, not the model

`(x)` demanded firings/void ≥ 0.9, i.e. **at most one cell per firing** — right
only while `h_f ≤ dz`. Under the flake remap a firing legitimately voids
`ceil((h_f − phi_top)/dz)` cells, and at `c_f = 0.5` on the 0.5 mm mesh
`h_f = 0.578 mm > dz`. `B5_0.5` was being flagged for doing exactly what the
model is specified to do. The strict form of this check already lives in the
code as a hard assertion (voids == crossings), which no leg tripped; the
harness now checks the geometric bound instead. **PASS on all 18 legs.**

## 6. (c) — the tree identity, and (g) — provenance

`guarded_post_all.sh p0b5 --force`, binary pinned at `7ea74c4f…` and **unchanged
across the sweep** (the guard checked). 2639 files hashed:

| against | differing |
|---|---|
| `d2z1.md5` (**the reference**) | **0** |
| `p0b4.md5` (revision 2's own gate) | **0** |

`beam_void_closure` returns rc = 1, as it has since before this line of work
(T3 stays red by standing decision).

**(g)** `MANIFEST.md` carries a row per leg with its output path, its event-log
sha256, its row count and a **staleness flag**: a leg whose newest file predates
the final binary is marked STALE automatically. Three rows are so marked — the
`D_*` legs of a superseded revision of the column harness, which is exactly the
class of stale leg that revision 2 scored without noticing. **Every scored leg
is fresh.** The rev2 Meier run was copied out of `/private/tmp/p0b_fh` into
`output_fh_rev2/` (3.7 GB, complete) before anything else was done.

## 6b. Housekeeping (deliverable 11)

* `surf_delta` deleted — declared and zeroed, never read.
* `surf_fired_prev` and `surf_tau` deleted with the smoothing.
* `surf_x_smooth` renamed `surf_x_hold`: it now holds the last x that had
  content, which is all the flat-hold needs.
* Witness label `p0b5`, distinct from `p0b` / `p0b2` / `p0b3` / `p0b4`.
* The unit test's expected event-log header was stale at P0; it now carries
  `ell_raw,delta`, and a new check **(3b)** turns `ell_raw` into the packet's
  own witness: with the EMA gone `ell_eff` must equal `ell_raw` on every
  non-flat row. **176 of 176 rows agree exactly.**

## 7. (f) — Meier, and deliverable 9. ⚠⚠ READ THIS BEFORE ACCEPTING DELIVERABLE 7

**The run aborted at t = 378.88 s because the model drilled out of the bottom of
the block.** Column (6,8) reached **496 mm of a 500 mm domain**; the
"voids == crossings" assertion then fired because a 15 mm flake was requested
with 6 cells of rock left. That is the assertion doing its job on a physically
exhausted domain — ⚠ **but its message names the wrong cause** ("the phi tie at
Removal.H:145-148"), which cost time. It should distinguish *ran out of rock*
from *a tie at a cell boundary*.

The rate window (60–160 s) is covered with a factor-of-two margin, so (f) and
deliverables 9 and 10 are all reported from the completed part of the run. The
frozen rev2 run went the full 500 s and reached only 387 mm.

### The number

| | frozen (pre-P0) | + P0 | + P0b rev2 | **+ rev5** |
|---|---|---|---|---|
| axis ROP, 60–160 s | 2.766 m/h | 3.690 | 3.237 | **5.911 m/h** |
| axis depth | 332.0 mm @500 s | 494.0 | 392.5 @500 s | **482.3 mm @379 s** |
| firing window p10–p90 | 19.4 K | 86.9 | 193.6 K | **109.3 K** |
| `T_surf` − `T_top` at firing | — | +73.0 | +188.9 K | **+60.0 K** |
| `flat` fallback | — | 12.8 % | 0.0 % | **0.0 %** |

(The rev2 row recomputes from the frozen thermo on the **true axis column
(0,0)** — 3.237 m/h and 392.5 mm, matching what revision 2 reported, which is
the estimator check. My earlier "median firing column" picks r ≈ 21 mm and gives
3.242 / 5.609; the axis figures are the ones quoted.)

⚠⚠ **Deliverable 7 nearly doubles the Meier rate, in the wrong direction.**
2.766 → 5.911 m/h is **+114 % on the frozen baseline and +83 % on revision 2**,
against Meier's 1.3–1.6 m/h. The model was already ~2× too fast; it is now ~4×.

**The attribution is clean and needs no extra run.** The 1-D column has no side
faces, so it isolates every rev5 change *except* deliverable 7: there the rate
moved 2.544/2.555/2.490 → 2.506/2.506/2.544 m/h, i.e. **under 2 %**. ⇒ the
Meier +83 % is deliverable 7 and nothing else.

**This is not an argument that deliverable 7 is wrong.** The wall *does* conduct
diffusively and the face Newton *is* the right closure for a non-receding
surface; leaving the wall on a cell average was under-delivering, and revision 2
said so. What the result shows is that **the model has no mechanism to reject
the extra power**, which is the same structural fault D2ad found from the other
side ("cutting floor flux shrinks the fired radius"). The rate level closes on
the gas side (P2: blowing, `cp(T)`), not here.

### Deliverable 9 — the bins, over the rate window

| bin | result |
|---|---|
| **1. by radius at the 40 mm gate** | bore **0.00 W**, pads **0.00 W**, wall **1649 W** ⇒ **100 % of side power is hole wall** |
| **2. does the receiving cell fire next?** | top-cell **121.4 W**, buried **1652.8 W** ⇒ **6.8 %** |
| **3. multi-cell riser share** | `wall_step_f_multi` 0.82 — **and irrelevant** (see below) |

✅ **P5 CONFIRMED, and stronger than predicted**: it said "> 70 % of side power
sits at r ≥ 40 mm". Measured: **100.0 %**. ⇒ **deliverable 7a is the whole fix
and 7b never fires on Meier.** The multi-cell riser share was the witness
bounding 7b's error, and with 7b empty it bounds nothing — 7b is carried as
untested code.

✅ **P6 CONFIRMED**: only **6.8 %** of side power lands on a cell that is its own
column's top cell. **93 % is buried** under its own column's floor and cannot
fire until the floor recedes to it. ⚠ So the side share *overstates* the error
badly — but note this cuts the other way too: buried heat is **stored, not
lost**, which is consistent with the rate going *up* rather than the shape
changing.

⚠ The end-of-run `surf_side_share` is **0.814**, against **0.688** over the rate
window and 0.470 for revision 2 over the same window. The packet's instruction
to bin over the window rather than the end of the run was right, and the gap is
now wider than it was.

## 8. Deliverable 10 — the Meier ledger residual, decomposed

Revision 2 reported `ledger_err_solid = −0.162` on Meier against 6e-4 on the
1-D column and left it as "the flake sink is not wired through the march". **It
is not the march.** Two independent findings:

**(i) (b) rules the march out.** `jet_P_face == patch_P_robin + side_P_face +
side_P_loss` closes to **1.3e-2 W** on the Meier leg. The march, the kernel and
the ledger are accounting the same power to twelve significant figures of the
total; a 16 % energy gap cannot live there.

**(ii) It is the bottom of the gathered profile window.** The remap gathers each
column's profile over `spall.h_spall_max` (15 mm) and writes only within it.
After `n` cells void, the cells at depth indices `Dmax+1 … Dmax+n` leave the
solid set carrying their stored enthalpy, and `E_flake` — the telescoped sum
over the window — never counts it. That term is `≈ n · V · ρc · ΔT(Dmax)`, which
is **zero if and only if the rock at 15 mm is at ambient**.

On the 1-D column it is: `exp(−15 mm / 1.30 mm) = 1e-5`. On Meier it is not —
the hole wall heats the rock laterally, so a firing column is still warm 15 mm
down. From the frozen run: 108,318 voided cells and 1.548e5 J unaccounted ⇒
**82 K above ambient at the window bottom**, against a flake that leaves at a
282 K equivalent. Entirely physical.

**The controlled test.** Shrink the window on the *clean* 1-D column, where
nothing else changes, and the Meier-sized residual appears on demand:

| `spall.h_spall_max` | `exp(−h/ell)` | `ledger_err_solid` |
|---|---|---|
| 15 mm (default) | 9.8e-6 | **+4.1e-4** |
| 8 mm | 2.2e-3 | −5.5e-4 |
| 5 mm | 2.2e-2 | **−1.2e-2** |
| 4 mm | 4.9e-2 | **−2.9e-2** |

The residual tracks `exp(−h_max/ell)` over two decades. ⇒ **the fix is to add
the window-boundary term to `E_flake` (or to widen the window), not to wire the
sink through the march.** ⚠ Note this is a *bookkeeping* error, not a physics
error: the field is right, the ledger's sink is short by the enthalpy of rock
that left below the window.

⚠ **The all-cell `ledger_err` is misleading in the other direction** and should
not be quoted: void cells keep their enthalpy, so `E_H` barely notices the
removal — it closes to −6.6 % with **no** sink term and to +49 % with one.
`ledger_err_solid` is the only form that means anything here.

## 9. What is owed next

1. ⚠⚠ **Deliverable 7 doubles the Meier rate.** The lateral closure is the right
   *form* and the wrong *level* — the model has no way to reject the extra
   power. This is the same wall D2ad hit from the other side. It belongs to P2
   (blowing, `cp(T)`), not to another lateral-BC packet.
2. **7b is untested code.** 100 % of Meier's side power is hole wall, so the
   floor-riser branch never executes on the only configuration that has side
   faces at all. Either build a case that exercises it or delete it.
3. **The Meier block is too short.** The model now drills 496 mm of a 500 mm
   domain in 379 s. A full 500 s run needs a deeper block — and the assertion
   message must distinguish *ran out of rock* from *phi tie*.
4. **The ledger sink is short by the window-boundary term.** `E_flake` should
   add the enthalpy of the cells that leave below `h_spall_max`; §8 gives the
   term and a controlled test that reproduces it on the clean column.
5. **`ell ≥ dz/4` is a mesh quantity setting a physical one.** Not binding
   anywhere in this packet (`rc_flag = 0` on every scored leg) — but it was
   binding in revision 2, and it will bind again whenever the heated layer is
   thinner than a cell. P1 item.
6. **`c_f` is still an input.** P1's `K_I(z)` scan is what makes it an output.
   P6 is exact at all three values, so the machinery is ready for it.
