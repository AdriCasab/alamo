# P0b revision 7 — bound the flake, re-run Meier, say what the cap is not

Revision 5's fixes stand and were not re-done. This revision adds one empirical
bound, one diagnostic split, one ledger term, one message, and re-measures Meier.

**Final binary** — `BINARY_rev7.sha256`. Nothing was rebuilt after the scored runs.

## ⚠⚠ What the cap is, and what it is not

`spall.flake_h_max` bounds the flake **thickness** and nothing else.

It exists because on a side-heated Meier column the two-cell temperature ratio
goes nearly flat — 1.3 where a front receding at 5.9 m/h demands 117 — so
`ell_eff` reads 7.5 mm against a physical `κ/v` of 0.42 mm and `h_f = c_f·ell_eff`
follows it. **A flat ratio means the gradient is unresolved**, not that the flake
is thick.

**It does not fix the physics.** The firing criterion is **gradient-blind by
construction**: a Tada edge crack over a 4–20 µm flaw sees an effectively uniform
stress, so a flat slab at the firing temperature fires exactly as a steep skin
does. Nothing in this model says a preheated slab should spall less. The real fix
is a gradient-aware criterion — buckling or energy release with a flake thickness
— and that is a **P1 design item**. ⚠ P1's planned `K_I(z)` scan does **not** fix
it either: run on a flat profile it returns a thick slab too. The cheaper P1 route
worth testing first is a **grain-scale flaw** (0.5–1 mm ≈ `ell`), which makes the
*existing* `K_I` integral gradient-sensitive with no new law.

**This key buys the Meier comparison; it does not earn it.**

### What it is applied to

`h_f` only. `x` and `T_surf` stay raw, because on a genuinely flat profile the
correct surface temperature **is** the top-cell average and the raw
reconstruction returns exactly that (`fac → 1` as `x → 0`). Holding a steeper `x`
would fire the column hundreds of degrees early — that is the h500 clamp failure
revision 5 removed. Revision 6's deliverable 1a was withdrawn for this reason and
nothing here reinstates it.

### It is not `spall.h_spall_max`

That key is a **search depth** — how far down the model looks for a flake at all
— and revision 5 measured it acting as the thickness bound on 13.5 % of Meier's
firings. They are now separate quantities with separate keys, and the parser
**refuses** `flake_h_max >= h_spall_max`, because a bound at or above the search
depth is inert and hides that the search depth is doing the bounding.

### Empirical, and declared so

1.5 mm on the Meier legs, from Wilkinson's measured 0.1–1 mm flake band and
Rauenzahn Fig 4.12. Default **0 = inert**, so no frozen leg and no control sees it.

## Two reporting defects found on the way

⚠⚠ **The event log's `rc_flag` never carried the cap bits.** It encoded
flat / `ell` / `T_surf` only, so revision 2's "the `h_f` cap binds on 85 % of
columns but **0 % of firings**" was reading a bit that is not written. Bits 8
(search depth) and 16 (flake bound) are now emitted where the firings are. The
column diagnostic `recon_hf_capped` was always right; the per-firing one was not.

⚠⚠ **`surface_patch.side_face_T` already exists** (D2ab's side gas temperature,
`= live` on every Meier leg). The new switch was first written under that name
and aborted at parse on the first smoke run. It is `side_face_closure`.
**A grep for the key would have found it; only the smoke run did.**

## The revision-5 criteria: re-run, not re-derived

Every one reproduces **digit for digit** on the revision-7 binary, with the flake
bound **carried on the column legs** so (a) is measured in the shipped
configuration rather than in a configuration nothing ships.

| | rev 5 | rev 7 | |
|---|---|---|---|
| **(a)** lsq slope, 2/1/0.5 mm | 2.506 / 2.506 / 2.544, **1.51 %** | 2.506 / 2.506 / 2.544, **1.51 %** | PASS, bar 2 % |
| controls C / F / P | 11.43 / 15.37 / 13.94 % | 11.43 / 15.37 / 13.94 % | not vacuous |
| **(d)** anchor h150 / h500 | 0.161 % / 0.116 % | 0.161 % / 0.116 % | PASS, bar 0.5 % |
| clamped legs | none | none (`rc_flag = 0` × 6) | PASS |
| controls' true surface | 14.06 % / 36.63 % | 14.06 % / 36.63 % | > 3 % |
| **(e)** 4 layouts, flake path | `8c3c0025…` | `8c3c0025…` | 0 differing |
| **(c′)** | −1.06 / −0.35 / +0.71 % | −1.06 / −0.35 / +0.71 % | PASS |
| **(c″)** | 4.89 / 4.14 / 1.51e-4 | 5.06 / 4.16 / 1.52e-4 | PASS, bar 5e-3 |
| **P6** | exact at 4 dp | exact at 4 dp | |
| unit test | 7 checks | 7 checks, 6 guards | PASS |

⚠ (c″) moves in the fourth digit — that is **deliverable 5's boundary term**,
which is ~1e-5 of `E_src` on a column whose rock is at ambient 15 mm down. It is
supposed to be invisible here; that it is, is the check.

✅ **P7 CONFIRMED**: (a) is untouched by the cap. Not by argument — the flake
bound was **switched on for every column leg** and the (e) hashes are still
character-for-character revision 5's, so it never once bound.

## (i) — the cap discriminates

| | binds | `h_f` range |
|---|---|---|
| 1-D column, 2.0 mm | **0.00 %** | 0.518–0.696 mm |
| 1-D column, 0.5 mm | **0.00 %** | 0.449–0.568 mm |
| 1-D column, all 18 legs | **0.00 %** (worst leg) | |
| Meier (45 s smoke, t ≥ 30 s) | **15.9 %** | p50 0.454, p90 1.500 mm |

The column's flakes are 0.45–0.70 mm against a 1.5 mm bound, so it cannot touch
them; on Meier it catches the unresolved-gradient tail.

## (f) — Meier with the bound. P1 and P3 confirmed, P2 and P6 refuted

All four runs reduced with **one estimator**, which reproduces the packet's
quoted revision-2 and revision-5 figures exactly (volume 2.216e-6 / 2.712e-6,
power 1180.3 + 999.7 / 808.3 + 1649.1 W, standoff 134.4 / 205.2 mm) — so the
comparison is like for like.

| 60–160 s window | rev 2 | rev 5 | **rev 7 follow** | rev 7 cell |
|---|---|---|---|---|
| axis ROP [m/h] | 3.237 | 5.911 | **3.235** | 3.232 |
| removed volume [m³/s] | 2.216e-6 | 2.712e-6 | **1.793e-6** | 1.841e-6 |
| `h_f` p50 / p90 [mm] | 0.602 / 2.596 | 3.015 / 15.000 | **0.709 / 1.500** | 0.687 / 1.500 |
| bound by `flake_h_max` | — | — | **30.6 %** | 30.6 % |
| bound by `h_spall_max` | 0.0 % | 0.0 % | **4.5 %** | 4.8 % |
| floor power [W] | 1180.3 | 808.3 | **1173.2** | 1159.7 |
| side power [W] | 999.7 | 1649.1 | **906.6** | 989.5 |
| min standoff [mm] | 134.4 | 205.2 | **133.1** | 133.6 |
| `ledger_err_solid` | −0.162 | −0.269 | **+1.0e-4** | +9.9e-5 |
| power identity max [W] | 1.37e-2 | 1.27e-2 | **9.9e-3** | 1.37e-2 |

✅ **P1 CONFIRMED** — the rate returns to **3.235 m/h** (predicted 3.0–3.5) and
the standoff to **133.1 mm** (predicted ≈130). ⇒ the thickness *was* the cause,
and withdrawing "the wall closure doubles the rate" was right.
✅ **P3 CONFIRMED** — floor power recovers **808.3 → 1173.2 W** (predicted ≈1180).
✅ **P4 CONFIRMED** — binds 30.6 % on Meier, 0.00 % on the column.
✅ **P5 CONFIRMED** — at the live 40 mm gate the split is **371.6 W inside /
610.8 W outside = 38 / 62 %**, not 100 / 0. The dead-gate bins still read
`side_P_bore = side_P_pads = 0.00 W`, caught in the act.
✅ **Deliverable 5 closed the Meier ledger**: `ledger_err_solid` **−0.162 →
+1.0e-4**, the same order as the clean 1-D column. The window-boundary term was
the whole residual, as revision 5's controlled experiment predicted.

### ⚠⚠ P2 is REFUTED, and the sign is the opposite of revision 5's claim

At the **same flake bound**, switching the lateral closure off changes almost
nothing: **volume −2.6 %, axis ROP +0.1 %, delivered power −3.2 %.** Not the
+10–20 % predicted, and not remotely the doubling revision 5 claimed.

⚠ **The face Newton delivers LESS side power than the cell average** — 906.6 W
against 989.5 W — because it solves for a face temperature *above* the cell
average, which cuts the driving `h·(T_gas − T_f)`. Revision 5's "+65 % side
power" was a property of two different trajectories, not of the closure. **My
revision-5 attribution was wrong in exactly the way the planner said, and this
run reverses its sign.** The withdrawal stands and is now measured.

### ⚠ P6 is REFUTED — removal is not purely energy-limited

Removed volume is **not** invariant: 2.216e-6 → 2.712e-6 → **1.793e-6 m³/s**.
From revision 5 the delivered power fell 15 % (2457 → 2080 W) while the volume
fell **34 %**. How removal is parcelled does change how much rock leaves, so the
"energy-limited, barely moves" framing does not hold at this bound.

## ⚠⚠ (h) FAILS — and it cannot be fixed by choosing a different bound

**Capped fraction 30.6 %, against a < 20 % bar.** (`h_f` p50 = 0.709 mm passes
the 0.3–1.0 mm band.)

This is not a mis-set bound. The unbounded `c_f·ell_eff` distribution over the
window is **bimodal**: p10 0.426, **p50 0.709**, p75 2.242, p90 4.394, p99 15.000 mm.
No value inside the declared empirical band gets under the bar:

| bound | binds on |
|---|---|
| 1.0 mm | 38.6 % |
| **1.5 mm (used)** | **30.6 %** |
| 2.0 mm | 27.1 % |
| 3.0 mm | 17.0 % |
| 5.0 mm | 8.4 % |

Only 3–5 mm clears 20 %, and that is three to five times Wilkinson's measured
0.1–1 mm band — choosing it to pass the criterion would be tuning, so it was not
done.

**And the grip tightens as the hole deepens:**

| t [s] | 30–60 | 60–90 | 90–120 | 120–160 | 160–200 |
|---|---|---|---|---|---|
| binds | 18.7 % | 24.3 % | 31.6 % | 36.3 % | **38.6 %** |
| `ell_eff` p50 [mm] | 1.21 | 1.50 | 1.81 | 2.10 | **2.33** |

⇒ **About a third of Meier's firings happen on columns whose gradient the 2 mm
mesh does not resolve, and that fraction grows monotonically with depth.** By the
criterion's own logic — "if it binds on half the firings the cap *is* the
thickness law" — the cap is partway there and heading the wrong way. That is the
result (h) was written to detect, and it is the strongest argument yet for the
**P1 gradient-aware criterion**: the stopgap does not converge, it degrades.

## (c) identity and (g) provenance

`guarded_post_all.sh p0b7 --force`, binary pinned and **unchanged across the
sweep** (the guard checked). 2639 files:

| against | differing |
|---|---|
| `d2z1.md5` (**the reference**) | **0** |
| `p0b5.md5` (revision 5's gate) | **0** |

So `spall.flake_h_max`, `surface_patch.side_face_closure` and
`surface_patch.side_power_interval` are all genuinely default-off, and
`cell` / `face` / `pinned` are byte-identical. `beam_void_closure` returns
rc = 1, red by standing decision.

**(g)** every scored leg has a `MANIFEST.md` row with its output path, event-log
sha256, row count and a staleness flag against the final binary.

## (b) the power identity

**9.9e-3 W** on the Meier leg (bar 0.02), with the march carrying ~0.9 kW of
side power, and **1.37e-2 W** on the `cell`-closure control. The identity spans
kernel ↔ march ↔ ledger and both new side bins were added to the ledger's own
accumulation, so it is also the check that the re-bin did not double-count: the
two independent splits of the same side power (by radius, and by whether the
receiving cell is its column's top cell) agree to **2.4e-5 W** and reproduce
`surf_side_share` to four decimals.

## Deliverables

| | |
|---|---|
| **1** flake bound | `spall.flake_h_max`, default 0, `h_f` only, parser refuses ≥ `h_spall_max` |
| **2** Meier re-run + live-gate re-bin | `surf_side_P_rin` / `surf_side_P_rout`, plus a per-column `_side_power.csv` (`surface_patch.side_power_interval`) — no plotfile carries a side flux |
| **3** deeper block | 0.6 m by extending **downward** (`prob_lo z = −0.1`), so the surface, the standoff and `dz` are untouched: verified as exactly 100 mm of added rock, first firing 0.288 vs 0.304 s |
| **4** assertion message | now compares the flake's demand with the cells available above the floor and says "RAN OUT OF ROCK … this is NOT the phi tie" |
| **5** window-boundary term | in `E_flake` at void time in PASS 3; Meier `ledger_err_solid` −0.162 → **+1.0e-4** |
| **6** closure on/off | `surface_patch.side_face_closure = follow \| cell`, one binary, one key |

## The full 500 s run

`rc = 0`, all 31250 steps, **no abort** — where revision 5 died at 379 s having
used 496 mm of a 500 mm block.

| | |
|---|---|
| axis depth at 500 s | **396.9 mm** (rev 2: 392.5 mm at 500 s; rev 5: 482.3 mm at 379 s, aborted) |
| deepest column | 398 mm of the **600 mm** block — 200 mm of headroom |
| firings | 200,836 |
| `ledger_err_solid` at end | **+5.1e-5** (rev 2 −0.162, rev 5 −0.269) |
| `surf_side_share` at end | 0.592 |

⚠ The depth at 500 s is back within 1 % of revision 2's, which is the
independent confirmation of P1: the same rock, the same time, the same hole.

## What is owed next

1. ⚠⚠ **The gradient-aware criterion (P1).** (h) fails at 30.6 % and the bind
   fraction **grows with depth** (18.7 → 38.6 %). The stopgap does not converge,
   it degrades, and no bound inside the measured 0.1–1 mm flake band clears the
   bar. Test the **grain-scale flaw** route first — a 0.5–1 mm flaw makes the
   *existing* `K_I` integral gradient-sensitive with no new law.
2. **`c_f·ell_eff` is bimodal on Meier** (p50 0.709, p75 2.24, p90 4.39 mm). The
   resolved half is physical; the tail is an unresolved gradient. A criterion
   that could tell them apart would not need a bound at all.
3. ⚠ **The wall closure is worth ≈ nothing at this bound** (−2.6 % volume,
   +0.1 % rate). It is still the more defensible closure, but it is not where
   the rate lives — **P2 blowing + `cp(T)`** is.
4. **7b (the floor-riser branch) is still untested.** The live-gate split gives
   38 % of side power inside 40 mm, so unlike revision 5's dead-gate reading
   there *is* a population to exercise it — but nothing yet checks it.
5. **Removal is not energy-limited** (P6 refuted): volume fell 34 % from
   revision 5 while power fell 15 %. How removal is parcelled changes how much
   rock leaves, which is worth understanding before the next calibration.
