# D2c hand predictions for the D2d Meier campaign

Written 2026-09-17 13:53 by `predict.py`, **before any D2d run**. The sha256 and mtime of this file are
recorded in `RESULTS.md` §0. Nothing here may be edited after a D2d run starts.

## Model

- **Steady state:** `d2b_feet_rop/feetmodel.py` with the D2c closures (`Cfg`: far law,
  momentum D_e = jet_De_ref·√(T_mix/T_ent), core length). The D2b geometry is unchanged:
  the ring at r_q = 38.97 mm at s = 50 mm; s linear from s_c (r ≤ 18.75 mm) to 50 mm;
  face at T_fire = 821 K; T_rec = T_gas(r_q) as a fixed point.
  Deep in the hole no free-surface column marches, because the nozzle plane is below
  the original surface, so dilution does not enter the steady state.
- **Transient** (`handmodel.transient`, dt 0.25 s):
  - the ring (feet) and the centre advance at their pinned closed-form rates;
  - outer columns on the 2 mm mesh grid march while s > 0, with ALAMO's per-bin
    free-surface dilution;
  - T_rec = the mouth inlet T (the exhaust T with no free-surface column);
  - the annulus band r_q–40 mm is cut to the feet depth when proud (clearance).
  It gives the time to a steady window (the score's rule applied to the smooth
  series: 150 s with |ds_c/dt| < 0.02 mm/s and 50 s sub-fits within ±10 %), the wall
  profile, and the whole-excavation and mechanical volume. It covers the full jet;
  shares and rates per quarter follow by /4.
- **Anchors per T_nozzle:** Martin h_ref = `walljet.h_anchor('M', T)`; jet_De_ref =
  `nozzle.de_ref(T)` (γ 1.3, J including the pressure thrust).

## Model check against D2b (not a prediction)

The transient model with the D2b closures (clamp, nozzle D, core 5, no free surface, no
clearance), compared with the D2b scored run B_JM_A2:

| quantity | hand | D2b sim | hand / sim − 1 |
|---|---|---|---|
| burner ROP 171–341 s [m/h] | 1.54 | 1.36 | +13 % |
| centre depth at 341 s [mm] | 291.20 | 262.00 | +11 % |
| s_c at 341 s [mm] | 203.10 | 191.33 | +6 % |
| mean T_rec 171–341 s [K] | 1502.82 | 1483.88 | +1 % |

The hand model runs fast and hot against the simulation. The refutation bands below
are **±27 % on rates**, i.e. twice the hand model's D2b ROP error of
13 %, rounded; and **±20 mm on Ø and stand-off** (5 mesh cells in Ø).

## P1. Centre equilibrium, time to a steady window, domain height

Steady state (full jet; face power per quarter):

| case | ring ROP [m/h] | centre s_c | T_rec | T_stag | φ | D_e | P_face / 4 |
|---|---|---|---|---|---|---|---|
| scored (n = 1) | 1.93 | 180 mm | 1733 K | 1797 K | 0.385 | 8.64 mm | 0.79 kW |
| n = 0.5 | 1.96 | 335 mm | 1743 K | 1776 K | 0.207 | 8.67 mm | 0.74 kW |
| clamp | 1.82 | none (deep-pit limit) | 1683 K | 1690 K | 0.034 | 8.52 mm | 1.03 kW |

Transient:

| case | steady window from | ROP in window | run end (window + 50 s) | nozzle-plane depth then | s_c then | domain needed | centre 40 mm above the bottom of 0.30 / 0.40 m at |
|---|---|---|---|---|---|---|---|
| scored (n = 1) | 750 s | 1.88 m/h | 950 s | 430 mm | 177 mm | 0.65 m | 325 s / 498 s |
| n = 0.5 | 2205 s | 1.94 m/h | 2405 s | 1195 mm | 318 mm | 1.55 m | 289 s / 433 s |
| clamp | never (to 3000 s) | — | — | — | — | — | 249 s / 354 s |

**Domain-height rule** (ACTIVE_STEP deliverable 8): use 0.30 m if nozzle-plane depth
at the end of the steady window + s_c + 40 mm ≤ 0.30 m, else 0.40 m. The run end
(window + 50 s steady-stop) is used as the end.

- **Scored run: the rule gives 0.40 m, but neither height fits.** The
  hand steady window opens at 750 s. By the end of the run the
  nozzle plane is 430 mm deep, so the domain
  must be ≥ 0.65 m. In a 0.40 m domain the centre reaches the bottom
  watchdog at 498 s, before any window.
  **The hand model therefore predicts that D2d run 1 at 0.40 m stops at the bottom
  with no steady window.** That decision (a taller domain, or scoring without a window)
  is D2d's; this packet does not trim or change anything.
- **What delays the window:** the s_c drift criterion. The ROP sub-fits are within
  ±10 % from early on, while |ds_c/dt| < 0.02 mm/s is first met at 750 s.
- **n = 0.5:** s_c ≈ 335 mm; the domain would have to be
  ≥ 1.55 m. **n = 0.5 does not fit** either height, so D2d reports
  n = 0.5 from the hand model only (run 3 is not run).
- **Clamp:** no centre equilibrium (s_c → the 2 m bisection bound). The
  deep-pit limit is 1.82 m/h, so run 2 ends at the bottom
  (354 s at 0.40 m).

## P2. Ring ROP and its sensitivity to n

| case | steady ring ROP | transient ROP, 150 s to the 0.40 m bottom | ROP in the hand steady window |
|---|---|---|---|
| scored (n = 1) | 1.933 m/h | 1.79 m/h | 1.88 m/h |
| n = 0.5 | 1.956 m/h | 1.74 m/h | 1.94 m/h |
| clamp | 1.823 m/h | 1.66 m/h | — |

Computed sensitivity (steady):
- n = 1 → 0.5 changes the ring ROP by +1.2 %;
- n = 1 → clamp (n = 0) by -5.7 %.

The ring hardly feels n. The far law moves the centre (s_c), and the ring feels it only
through T_rec. **The scored steady ring ROP, 1.93 m/h, is above Meier's band [1.04, 1.92].**

## P3. T_nozzle sweep (scored closures)

| T_nozzle | Martin h_ref | jet_De_ref | steady D_e (at T_rec) | steady ROP | transient ROP 150 s → 0.40 m bottom | steady T_rec |
|---|---|---|---|---|---|---|
| 1600 K | 1413 | 0.0037512240761077564 m | 8.43 mm | 1.34 m/h | 1.24 m/h | 1482 K |
| 1750 K | 1431 | 0.003647390490008811 m | 8.54 mm | 1.64 m/h | 1.51 m/h | 1607 K |
| 1900 K | 1449 | 0.0035555966839916735 m | 8.64 mm | 1.93 m/h | 1.79 m/h | 1733 K |

Steady ROP rises 0.147 %/K (of the 1600 K value) between 1600 and 1900 K.
Linear interpolation gives 1.3 m/h at T_nozzle ≈ below 1600 K and
1.6 m/h at ≈ 1732 K.
This is conditional on J-M, n = 1, free-surface dilution and core 8, and is never reused
as an input.

## P4. Hole: mouth, depth-mean, minimum, volume (scored)

Two windows:
- **(a) 150 s to the 0.40 m bottom stop.** This is what D2d run 1 will measure, and the
  P6 bands use it.
- **(b) The hand steady window.** It needs a ≥ 0.65 m domain.

The profile is taken at the window end, over the drilled depth (0 to the nozzle-plane
depth). Ø is in mm. The volume rate is the whole excavation of the full hole, in cm³/s
(score target 1.98–2.96).

| case | window | drilled depth | Ø at 20 mm | Ø at 25 mm | Ø at 50 mm | depth-mean Ø | min Ø | whole volume rate | burner ROP [m/h] | mechanical share |
|---|---|---|---|---|---|---|---|---|---|---|
| scored, (a) to the 0.40 m bottom | 150–498 s | 194 mm | 138 | 134 | 118 | 111 | 94 | 4.33 | 1.79 | 0.04 % |
| scored, (b) hand steady window | 750–900 s | 403 mm | 138 | 134 | 118 | 99 | 86 | 3.18 | 1.88 | 0.06 % |

Reading of window (a) against the targets:
- **Mouth:** Fig. 8.8 has 92 mm at 25 mm and 88 mm at 50 mm (±10 mm). The hand gives
  134 and 118 mm, so the **mouth criterion is predicted to fail**:
  free-surface dilution narrows the funnel but does not remove it. The D2b closures give
  the Ø in the P5 row.
- **Minimum Ø:** 94 mm. The hand floor is 2·r_q = 78 mm; the wall
  grid is 2 mm, so Ø moves in 4 mm steps.
- **Depth-mean Ø:** 111 mm against the 85–93 mm target.
- **Volume rate:** 4.33 cm³/s against 1.98–2.96.
- **Mechanical share:** 0.04 % against the 5 % finding threshold.
  This is a lower bound: the hand model has one smooth clearance band (r_q–40 mm),
  while the mesh has discrete columns proud of a nearest-rank pad.

## P5. Variants (same rows)

Hole metrics are for window (a), 150 s to the 0.40 m bottom stop.

| case | steady ROP | steady s_c | hand window from | domain needed | bottom stop (0.40 m) | ROP (a) | Ø 20 | Ø 25 | Ø 50 | depth-mean Ø | min Ø | volume rate (a) | mech share (a) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| scored (n = 1) | 1.93 | 180 mm | 750 s | 0.65 m | 498 s | 1.79 | 138 | 134 | 118 | 111 | 94 | 4.33 | 0.04 % |
| core 5 | 1.93 | 175 mm | 765 s | 0.64 m | 518 s | 1.75 | 142 | 138 | 122 | 113 | 94 | 4.35 | 0.04 % |
| Ricou (core 3.125) | 1.93 | 172 mm | 780 s | 0.64 m | 532 s | 1.72 | 142 | 138 | 122 | 114 | 94 | 4.37 | 0.04 % |
| D2b closures + n = 1 | 1.93 | 174 mm | 800 s | 0.64 m | 557 s | 1.71 | 190 | 174 | 134 | 125 | 94 | 4.70 | 0.03 % |
| aspect 0.5 | 1.93 | 180 mm | 755 s | 0.65 m | 502 s | 1.78 | 146 | 142 | 122 | 114 | 94 | 4.40 | 0.04 % |
| aspect 2 | 1.93 | 180 mm | 745 s | 0.65 m | 494 s | 1.80 | 130 | 126 | 114 | 108 | 90 | 4.24 | 0.04 % |
| exponent 0.8 | 1.93 | 180 mm | 750 s | 0.65 m | 498 s | 1.79 | 142 | 138 | 122 | 113 | 94 | 4.35 | 0.04 % |

## P6. What D2d would refute

Bands: ±27 % on rates, ±20 mm on Ø and stand-off (see the model check).

- **P1, scored:** refuted if run 1 finds a steady window before its centre reaches the
  0.40 m bottom, or if its s_c settles (drift < 0.02 mm/s for 150 s) outside
  160–200 mm.
- **P1, clamp (run 2):** refuted if run 2 shows a centre equilibrium, i.e. |ds_c/dt|
  < 0.02 mm/s for 150 s before the bottom.
- **P1, n = 0.5:** not tested by D2d (does not fit); refutable only with a ≥ 1.6 m domain.
- **P2:** refuted if run 1's burner ROP from 150 s to the bottom stop is outside
  1.31–2.27 m/h.
- **P2, clamp:** refuted if run 2's burner ROP from 150 s to its bottom stop is outside
  1.21–2.10 m/h (hand 1.66 m/h).
- **P3, 1600 K:** refuted if the ROP from 150 s to the bottom stop is outside
  0.91–1.57 m/h.
- **P3, 1750 K:** refuted if the ROP from 150 s to the bottom stop is outside
  1.11–1.91 m/h.
- **P3, ordering:** refuted if ROP is not increasing in T_nozzle across runs 4, 5 and 1.
- **P4 (window a):** refuted if run 1's Ø at 25 or 50 mm is within ±10 mm of Fig. 8.8, i.e.
  the mouth passes; the hand says 134 / 118 mm, ±20 mm.
  Also refuted if the depth-mean Ø is outside 91–131 mm,
  the mechanical share exceeds 5 % (hand 0.04 %), or the whole volume rate
  is outside 3.17–5.48 cm³/s.
- **P5, run 6 (D2b closures + n = 1):** refuted if its Ø at 25 mm is not wider than run 1's
  (hand: 174 vs 134 mm; the free-surface dilution is the funnel
  lever), or if its ROP from 150 s to the bottom differs from the hand value 1.71 m/h by more than ±27 %.
- **P5, core 5 / Ricou (optional runs):** refuted if their ROP ordering against run 1
  differs from the hand's (core 5 1.75, Ricou 1.72, scored 1.79 m/h).
- **P5, aspect / exponent:** hand only (no D2d run). Reported for the size of the
  mouth lever.

## dt for D2d (overshoot rule, 2 mm)

The hand transient's maximum T_stag is 1900 K (scored). run.py evaluates the
overshoot rule, ≤ 5 K per step, at q = h(0, s ≥ 50 mm)·(T_stag,max − 821) − loss.

