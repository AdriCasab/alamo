# P1f Leg B — RESULT on the corrected binary: every gate green

Leg B is `output_rev7sc`'s argv **verbatim** + `spall.lateral_drive=1` +
`spall.lateral_removal=1` + a new `plot_file`. The diff is **+3 / −1** keys and nothing
else (`output_legB/cmd.txt`), asserted at build time: `side_face_staircase=1` present
exactly once, and **no** `spall.driving_stress` and **no** flaw key on a control's argv.
Binary sha in `output_legB/BINARY.sha256`. The pre-fix run is kept as
`output_legB.void_unitbug/` (3.6 GB) — the bug's only witness.

## The one change that mattered

`Removal.H:624` compared `K_I` in **Pa·√m** against `K_Ic` straight from `KIcNasseri` in
**MPa·√m** — a **1e6×** score inflation. Deliverable 9 converts both branches.

| | control | **Leg B corrected** | Leg B VOID (bug) |
|---|---|---|---|
| events | 194,771 | **194,674** | 395,023 |
| voided cells | 88,628 | **87,716 (0.990×)** | 289,017 (3.261×) |
| regime 1 | 88,628 | 86,606 | 64,658 |
| **regime 6** | — | **1,110** | **224,359** |
| interior-cavity cells at 500 s | 0 | **9 (0.010 %)** | 72,925 (25.2 %) |
| `ledger_err_solid`, max t > 50 s | 1.453e-4 | **1.453e-4** | 8.771e-3 |
| regime-6 radius, max | — | **62.8 mm** | 166.9 mm |
| regime-6 `riser_depth`, max | — | **6.0 mm** | 374 mm |

Regime 6 is **202× smaller**, the volume is back to **99.0 %** of the control, and the
cavity population is **8,100× smaller**.

## Scored

| # | bar | measured | verdict |
|---|---|---|---|
| **(u)** | ≤ 1 % | logged `q_lat` **1.000634** vs hand **1.000634** = **0.000 %**; `K_Ic_lat` min **7.929e5** Pa·√m | **PASS** |
| **(v)** | monotone, one crossing, fires 0/0/>0 | `lat_qmax` **0.5202 < 0.8190 < 1.0006**, 1 crossing; **0 / 0 / 48** | **PASS** |
| **(v2)** | ≤ 5 % | worst **0.73 %** (0.01 / 0.01 / 0.73) | **PASS** |
| **(w)** | all three | **0 < 48 < 80**; hot `lat_qmax` 1.0006 ≥ 1; **warm 0.8190 ∈ (0.5, 1)** | **PASS** |
| **(b2)** | abs > 0, diff = 0, byte-identical off | **48 / 0**; scores straddle 1.0 (**0.5004 < 1 ≤ 1.0006**); lateral-off pair **byte-identical**, 322,758 B | **PASS** |
| **(a)** | 0 of 2639 | **0** vs `p1f`, `p1e` **and** `p1b`; binary pinned across 447 s | **PASS** |
| **(j)** | byte-identical | **840 files** — Kant 357 + `p0_steady_column` 483 | **PASS** |
| **(d)** | 500.0 s, rc = 0 | `TIME = 500`, rc 0, finalized, **0** aborts | **PASS** |
| **(e)** | ≤ 5e-3 | last **+5.3391e-05**, max **1.453e-04**; control +5.2700e-05 / 1.453e-04 | **PASS** |
| **(m2)** | assert the scan never reads a buried void; counts | **0 violations** on 23 events over the 4 columns that carry one; **9** cells (0.010 %) | **PASS** |
| **(h4)** | exactly 0 at r < 40 mm | **0** of 1,110 (min radius exactly 40.0; p50 50.7, max 62.8) | **PASS** |
| (h1) | report; unusable if monotonic to the edge | `riser_depth` **p50 0.0, p95 0.0, max 6.0 mm**; all five windows concentrated in 0–5 mm | **bounded** |
| (h2) | ≤ Ø110 mm at ≥ 150 mm | worst **98.0 mm** at 150–200; mouth 50–60 reported separately at **100.0** | **PASS** |
| (x) | report | **0.990×** the control | **report** |
| (l) | report | regime 6 lives only in 40–60 (1,052 cells) and 60–90 (58); nothing beyond 90 mm | **report** |
| (p) | report | Leg B vs control by band: **−1.41 / −0.17 / 0.00 / −1.10 %** | **report** |
| (k) | report | cap binds on **99.8 %** of lateral firings (vertical: 61.9 %) | **report** |
| **(r)** | report deficit + factor vs 41.42 % | see below | **report — informative miss** |
| (c2) | report | field gaps **1.3 / 5.3 / 4.7 %** over `ell/dz` 0.52 / 1.22 / 3.11 — **non-monotone** | **miss** |
| (c3) | report | field **+4.72 / +3.79 / +6.11 %** — a band, no mesh trend | **report** |
| (k3), (P6) | report | **not scored** — see Takeaways | — |
| (f), (g), (i) | — | deferred to P1g; the control is already round at 1.015–1.027 | — |

## (r) — and the miss is the informative one

Firing rate per column, census-normalised, refusing any sector under 5 firings:

| band | regime | axis/col | diag/col | deficit | factor vs 41.42 % |
|---|---|---|---|---|---|
| 40–60 | **6** | 3.5000 | 1.2672 | **+63.79 %** | **0.6×** |
| 60–90 | 6 | 0.0705 | 0.0610 | +13.41 % | 3.1× |
| 40–60 | 1 | 93.6615 | 91.6107 | +2.19 % | 18.9× |
| 60–90 | 1 | 0.3691 | 0.3220 | +12.76 % | 3.2× |

⚠⚠ Regime 6's 40–60 mm deficit is **above** the per-column 41.42 % it replaces. Per the
packet this is the **informative** outcome and **indicts P1e's flux isotropy, not the
lateral criterion**: under the `absolute` form `h_f` never enters σ, so the criterion is
**`ell`-blind** and the face temperature — hence the deficit — is set by the
staircase-corrected **per-area flux**. ⚠ The tracked-front stopping rule is **not armed
here** and is **not** hung on (r); it stays **(f)**'s call on Leg A in P1g.

⚠⚠ **(r)'s drive half is NOT SCOREABLE on this leg.** The field reduction at r ≥ 40 mm
returns axis **59.5 mm** / diagonal **86.7 mm** — 30–43 cells. A two-cell log read at
those lengths is reading a nearly flat interior field, not a decay length, and it fails
the same containment test the (c2) sweep has to filter for. The +45.7 % it prints is
**reported as degenerate, not as a result**.

## Pre-registered predictions

| | prediction | outcome |
|---|---|---|
| **P1** | sparse but non-zero; ignition at the fresh 40–44 mm corner; bounded; volume near the control's | **CONFIRMED on all four.** 1,110 firings, p05 radius 43.6 mm, `riser_depth` max 6.0 mm, volume 0.990× |
| **P2** | (d,e) and (m2) clear once the firings are physical | **CONFIRMED.** (e) is now **indistinguishable from the control** (1.453e-4 on both) ⇒ the ledger term was never a second defect |
| **P3** | firing starts at 40–44 mm, not 44–50 | **CONFIRMED** — min exactly 40.0, p05 43.6 |
| **P4** | (r) shows a reduction | **REFUTED in 40–60 mm** (63.79 %, factor 0.6×), held in 60–90 (3.1×). Per the packet this indicts P1e, and no bar is set from it |
| **P5** | volume near 88,628, not 3.26× | **CONFIRMED** — 87,716 |
| **P7** | the bore moves little | **CONFIRMED** — ≤ 1.41 % in every band |

## What this says for P1g / P1h

1. ⚠⚠ **P1h stays a SUCCESSOR.** My earlier "P1h is a prerequisite for P1g" rested on
   72,925 cavity cells that were the bug's. The real figure is **9**, all within 1–3
   cells of the wall top, and the depth scan reads none of them. **Withdrawn, as the
   packet required, and now measured rather than assumed.**
2. ⚠⚠ **"The predicate must confine to a step corner" is also withdrawn** — and the
   corrected mode confines itself: `riser_depth` max **6.0 mm**, radius max **62.8 mm**.
3. **The `absolute` form makes this leg `ell`-blind**, so P1f cannot test the azimuthal
   thesis at all. (r)'s 40–60 mm miss points at **P1e's per-area flux**, and the
   `differential` form's own test belongs to P1g with Leg A.
4. ⚠ **The cap binds on 99.8 % of lateral firings** against the vertical mode's 61.9 %.
   Under `absolute` that does not reach the drive, but it will the moment P1g runs
   `differential` — where the cap is the only route by which an `ell` anisotropy reaches
   the stress at all.
5. ⚠ **No Meier rate, volume, width or standoff is quotable from P1 … P1f.** The bore
   numbers above are Leg-B-vs-control only, which is what (p) asks.
