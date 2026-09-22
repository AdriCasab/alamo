# D2f pre-registration: the contact support rule (pads, no clip), q = 1.0 / 0.9 / 0.8

Written 2026-09-18 17:27:07 by `predict.py`, before any D2f run. **Frozen: never edit.**
Its sha256 and mtime are in RESULTS §0. These are diagnostics. No number
here or in the runs is scored against Meier's band.

## Hand model

- **Model:** `d2c_steady/handmodel.py` with the scored closures:
  - J-M at 1900 K;
  - far law power n = 1;
  - momentum D_e, core 8;
  - exhaust recirculation with entrained mass;
  - free surface with aspect 1, exponent 1.
- **The ring is at the radius the quantile selects in an axisymmetric hole,**
  R_Q(q) = √(r_i² + q(r_o² − r_i²)), set by rebinding `feetmodel.R_Q`.
- **The clip is off** (`transient(clear=False)`).
- **Steady:** `feetmodel.equilibria`.
- **Transient:** 0–300 s on the 2 mm outer grid, with the 0.36 m centre stop.
  The ROP is a fit of the feet depth over 150 s to the end.

| q | R_Q [mm] | steady ring ROP [m/h] | transient ROP 150–300 s [m/h] | s_c / T_rec at 300 s [mm / K] | steady s_c / T_rec | Ø 25 / 50 mm at 300 s | min Ø to the feet depth [mm] |
|---|---|---|---|---|---|---|---|
| 1.0 | 40.00 | 1.874 | 1.678 | 152.8 / 1604 | 184 / 1729 | 134 / 122 | 80.0 |
| 0.9 | 38.97 | 1.933 | 1.741 | 150.7 / 1614 | 180 / 1733 | 134 / 118 | 77.9 |
| 0.8 | 37.91 | 1.998 | 1.809 | 148.5 / 1624 | 175 / 1736 | 130 / 118 | 75.8 |

- **q-sensitivity (hand):**
  - transient **-0.066 m/h per +0.1 of q** (-3.8 % of the mean);
  - steady -0.062 m/h per +0.1 of q.
- **Spread** of the three transient ROPs about their mean (1.742 m/h):
  **±3.8 %**.
- **Gate (ii) predicted to HOLD** (limit ±5 %).
- **Direction against D2d R1** (1.65 m/h steady, q 0.9 + the clip, 0.70 m):
  - the hand model with the clip on, at q 0.9, gives a transient
    1.741 m/h, against 1.741 m/h with it off;
  - in the hand model the clip does not touch the ring rate. It only cuts
    the annulus band, cutting 0.5 cm³ by 300 s;
  - so the hand model predicts **no change** in ROP from removing the clip.
    D2e's P-a at 2 mm ran +1–2 % faster than B0 (the simulation, not the
    hand model).
- **The minimum Ø to the feet depth is 2 R_Q in the hand model by
  construction** (the wall is floored at the ring). It carries no
  information beyond R_Q.

## Limits (stated)

- **The hand model is axisymmetric with no column scatter.** It cannot
  predict the mesh behaviour of the per-pad maximum, which is set by the
  single slowest column. **Gate (i) has no hand prediction.** The expectation
  below comes from D2e's diagnosis, not from the hand model.
- **Its q-dependence is geometric only:** a larger R_Q gives a lower h and
  cooler gas at the ring. In the simulation, q also picks which columns'
  scatter sets the pad. The hand spread is therefore a lower bound on the
  simulated q-sensitivity.
- **The hole numbers are direction only** (packet: the wellhead, two phases,
  and Fig. 8.8's lower-bound widths).

## Expectations and the decision rule (the packet's, fixed before any run)

| id | expectation | criterion |
|---|---|---|
| E1 (gate i) | Q10 converges | Q10 burner ROP, 1 mm / 2 mm over 150–300 s, within **5 %**, and every 50 s block after 50 s within 5 %. No hand number: the maximum rides on the slowest column. D2e's P-a (q 0.9) converged to +0.1 %, and whether q 1.0 keeps that is the open question |
| E2 (gate ii) | q does not control the rate | the 2 mm burner ROP over 150–300 s for q 0.8, 0.9 and 1.0 is within **±5 %** of their mean. Hand: ±3.8 %, so predicted to hold |
| E3 | the ordering | the ROP falls as q rises (a larger R_Q gives a lower h and cooler gas). Hand: 1.809 > 1.741 > 1.678 m/h |
| E4 | the transient ROP, 2 mm, 150–300 s, per q | the hand value ±27 % (the D2c P6 rate band): q 1.0 1.68, q 0.9 1.74, q 0.8 1.81 m/h |
| E5 | s_c at 300 s, 2 mm, per q | the hand value ±20 mm (the D2c P6 length band): 153 / 151 / 149 mm |
| E6 | Q09_2mm is P-a's configuration | thermo and CSV rows to t ≤ 200 s byte-identical to `d2e_mesh/output/Pa_2mm` (harness) |

**Outcomes:**
- **(a) E1 and E2 hold:** D2g scores **q = 1.0**, the literal "highest rock
  under each foot". The jolt allowance is shown to be immaterial.
- **(b) E2 holds, E1 fails:** D2g scores **q = 0.9**, justified as the jolt
  allowance and shown converged by P-a (and by Q09_1mm to 300 s, if run).
  - The report covers the maximum's mesh failure and its mechanism, from
    `rim.py`.
  - If Q09_1mm was not run, the adoption is conditional on running it in
    D2g.
- **(c) E2 fails:** q controls the rate. The report gives the ROP per 0.1 of
  q and stops. D2g then needs a physical basis for the allowance. **q is not
  picked to hit the band.**
