# P2-11 research note — what should the heat on each rock face be? (2026-10-28, planner; no code, no runs)

Question: replace the hand-typed position maps (Höser's Goldstein map + flange lid; Meier's G4) with ONE rule that
computes each face's heat from the hole's own shape. Before planning, check the literature and our own data.

## 1. Sources read
- **Rauenzahn 1986 / Wilkinson 1989 MIT theses** (in-hole k-ε CFD, mock-hole copper calorimetry, steady hole radii):
  `validation/hoeser/rauenzahn_wilkinson_inhole_heat_notes.md` (OCR'd; every claim with a page number).
- **Zuckerman & Lior 2006** (Adv. Heat Transfer 39), correlation tables; **Goldstein, Behbahani & Heppelmann 1986**
  (IJHMT 29:1227) abstract; **Martin 1977**; **Lim & Lee 1997** and **Fenot et al. 2008** (concave targets);
  **Höser & Rudolf von Rohr 2018** (Exp. Therm. Fluid Sci. 91:166) — no radial profile (one axis sensor).
- **Our data:** P2-9's Höser leg P47n, re-analysed offline: `r3_p47n.py` → `R3.txt`.

## 2. Findings

**F1. Höser's map uses an AREA-AVERAGE correlation as a LOCAL h — a defect, independent of any new rule.**
Goldstein et al. 1986's `(533 + 44 (r/D)^1.285)` form is "a correlation … for the average heat transfer" (abstract);
Zuckerman & Lior tabulate it as Nu_avg. The local value is (1/2r)·d(r²·Nu_avg)/dr ⇒ g_loc = 533(533 + 0.3575u)/(533+u)²,
u = 44(r/D)^1.285. Half-width: average 7.0 D, local 4.35 D. At r = 36 mm local/average = 0.76. Meier's G4 is NOT
affected: its radial term is exactly the local form of Martin's average (converted 2026-09-30, "disc average retired" 10-06).

**F2. In-hole CFD and the only in-hole measurement both use DISTANCE ALONG THE SURFACE.** Rauenzahn fitted his
mock-hole calorimetry with Q(s) = Q_max exp(−s/s₀), "s … the distance along the curved face … not the radial
distance" (p.416), justified by wall-jet similarity on flat, concave and convex targets (Hrycak 1981). His CFD flux
plots run "along curved wellbore surface" to 15 R_noz, past the hole radius. ⇒ the exposure idea is the established
convention, not an invention.

**F3. The heat profile is much narrower near stagnation than today's map, with a fat tail along the walls.**
Rauenzahn's in-cavity CFD (Z/R_noz = 10): half-maximum at 2.2–3.4 R_noz, 1/3 at ~3.5 R_noz, then a slowly falling tail
(0.12–0.16 of peak at 15 R_noz); the cavity tail holds 2–5× the flat-plate value at equal distance. Today's Höser map
(average used as local) keeps 0.88 of peak to 3 R_noz and halves at 14 R_noz.

**F4. The profile half-width scales with standoff: ≈ 0.23 × H.** Two CFD points: 1.8–2.7 R_noz at Z/R = 10 (Rauenzahn)
and 13.8 R_noz at Z/R = 60 (Wilkinson) both fit b ≈ 0.23 H; copper fits bracket it (exp half-width 0.17–0.46 H).
Measured steady hole radii R_h ≈ 0.41–0.66 H (mean 0.53; CFD slope 0.43–0.46) ≈ 2b. In both theses the radius comes
from the SHAPE of Q(s)/Q_bottom through u_n = V cos θ; a uniform change of level leaves R_h unchanged (W p.231).
No minimum flux appears in either model.

**F5. Concave confinement raises the coefficient but cools the gas.** Lim & Lee 1997: concave hemispheres raise local
Nu everywhere (D/d 6–12). Fenot 2008: strong confinement stalls the adiabatic wall temperature and lowers Nu downstream.
Rauenzahn: wall-adjacent gas ~900 K vs a 1,940–2,830 K jet. ⇒ h along the path AND a gas temperature that depletes
along the same path (our enthalpy march, re-indexed to the path) are both needed.

**F6. Our data (R3, P47n shapes, equal total conductance).** Share of heat within 10 mm of the axis:

| rule | 15 s | 30 s | 60 s | 90 s | 120 s |
|---|---|---|---|---|---|
| today (average as local, at r, own standoff) | 0.165 | 0.084 | 0.038 | 0.026 | 0.025 |
| local Goldstein at r (F1 only) | 0.196 | 0.101 | 0.048 | 0.034 | 0.032 |
| local Goldstein on the path, standoff to stagnation | 0.216 | 0.125 | 0.070 | 0.055 | 0.050 |
| cavity shape, half-width 0.23 H, on the path | 0.267 | 0.160 | 0.089 | 0.065 | 0.057 |
| Rauenzahn cavity shape at its own R_noz scale | 0.279 | 0.173 | 0.104 | 0.087 | 0.080 |

Floor share (n_z > cos 45°) moves 0.63 → 0.77 at 15 s and 0.42 → 0.55–0.61 at 120 s. Unlike the path alone (15 s
unchanged), the profile width acts from the first second — where the model is already 39 % too shallow.
Wall cells fire at a calorimeter-equivalent h of 2.8–5.3 kW/m²K, ≫ Kant's minimum ~0.5: a minimum flux is not the lever.

## 3. The general rule this points to (for /plan to design; nothing here is implemented)

h(face) = h_stag(H) · g( s / b(H) ),  q = h · (T_g(s) − T_s)

- H = nozzle → stagnation-point distance (one per jet per step); s = shortest distance ALONG the solid surface from the
  stagnation point (face-list graph, with a staircase correction tested on planes/sphere/bore); hardware (flange,
  burner body, skirt) blocks the path.
- h_stag(H): measured (Höser calorimeter vs L) or from the jet closure (Meier). b(H) = 0.23 H (theses; to bracket).
- g: a peaked core with a wall-jet tail — the in-cavity CFD shape (Rauenzahn Fig 6.16/6.20) normalised by its
  half-width; bracket against the local flat-plate shape (F1) for the tail.
- T_g(s): the existing enthalpy march, indexed by s instead of r.
- Same rule for Höser and Meier: only h_stag, H, the nozzle and the hardware differ. Deletes the flange lid
  (`r>0.02 and s<0.0633`), Höser's average-as-local Goldstein, and Meier's G4 taper/plateau/cut.
- Test targets: Höser depth/diameter vs time (Fig 5.10; dL/dt = 1.8 dD/dt); steady Meier/Tester R_h ≈ 0.45–0.53 H
  (P2-12). The level S stays a stated input (uptake gap open).

## 4. What is NOT settled
1. b(H) = 0.23 H rests on two CFD points (supersonic, 1,940–2,830 K torches) — Höser's jet is near-sonic methane-air at
   1,866 K, L = 8 D, inside a Ø40 × 35 mm confining ring. Bracket b (e.g. 0.15–0.35 H) and report, never fit.
2. The tail shape: cavity CFD (fat) vs flat-plate local (thin) — bracket.
3. Rauenzahn's high-pressure jets (p_jet/p₀ > 9) trap a hot bubble and peak OFF axis: not Höser's regime; check Meier's.
4. Graph distance on a staircase overestimates s by up to ~8 % (26-neighbour): needs its own offline gate.
5. Höser's lab test is TRANSIENT (burner fixed, H grows 63 → 140 mm); R_h ≈ 0.45 H is a STEADY-drilling law — do not
   score Höser against it.
6. The cut F1 alone is a straight correction of a misapplied correlation; it could be its own small step or the first
   leg of P2-11 (attribution: F1 alone, then the path, then the width).

## 5. Double-check (planner, same day) — CORRECTIONS that supersede §2–§3 where they conflict
Figures re-read by eye (`rauenzahn_wilkinson_inhole_heat_notes.md` renders of Figs 6.16, 6.20); quote p.416 verified.

1. **F4 WITHDRAWN as a law.** In-cavity Fig 6.20 (Z/R = 10, low p) halves at ≈ 3.2 R_noz = 0.32 H, not 0.23 H; the
   1.8–2.7 R_noz matched before was the FLAT-PLATE run (Fig 6.12), and Rauenzahn's flat plates at Z/R = 5 and 10 have the
   SAME width (nozzle-scaled). Wilkinson's 13.8 R_noz at Z/R = 60 is vs RADIUS on a prescribed ellipse. No width law in H
   is established; at Z/D ≤ 10 the width is nozzle-scaled.
2. **Her own crater discriminates the width (zero compute, from P2-7 a2, review-accepted):** her 30 s bowl implies
   q(17.5 mm)/q(axis) = 0.7095–0.7632 (n_q 1.44–1.83). Goldstein AS USED (average) 0.813 — outside; Goldstein LOCAL (F1)
   **0.716 — inside**; Rauenzahn's in-cavity width (half at ~13 mm) ≈ 0.4 — far outside (cf. P2-5's σ = 0.48 D pencil
   bore). ⇒ for Höser the profile WIDTH should be Goldstein's local width; the narrow supersonic-torch CFD width does not
   transfer (his torches were inside their supersonic core). F6's two "cavity shape" rows are therefore too narrow for
   Höser; the "LOCAL Goldstein on path" row is the candidate (share within 10 mm ×1.3 at 15 s → ×2.0 at 120 s).
3. **Tail shape check:** Goldstein's local shape rescaled to its half-width gives 0.10 at 4.6 half-widths; the in-cavity
   CFD gives 0.14 there — the local flat-plate shape already carries a wall-jet tail close to the cavity's. One shape
   (Goldstein local, in s/D) is enough; the "fat cavity tail" claim (2–5×) rested on cross-figure normalisations and is
   not needed.
4. **Minimum flux — restated.** Not active under TODAY's map (walls 2.8–5.3 kW/m²K calorimeter-equivalent). Under the
   candidate rule a mouth wall at s ≈ 86 mm, H ≈ 118 mm gets ≈ 7,000 × 0.21 × 0.44 ≈ 650 W/m²K, i.e. NEAR Kant's ~500:
   the minimum flux may become live at late times. Report it; do not add it in the same step.
5. **Range of validity.** Goldstein: cold air, Re 61k–124k, H/D 6–12, r/D 0.5–32. Her Re ≈ 2e4 (planner estimate,
   1866 K, near-sonic, D 7.9 mm) is below range; H/D grows 8 → 18 by 120 s (beyond 12 the width should start to grow
   with jet spreading — unvalidated, bracket fixed vs growing). Meier H/D ≈ 12–13: at the edge.
6. **Simplest rule after the check:** h(face) = h_stag(H) · g_GL(s/D), q = h·(T_g(s) − T_s), with g_GL the LOCAL
   Goldstein shape, s the surface distance from the stagnation point, H the nozzle→stagnation distance, hardware as gas
   barriers, T_g the existing march indexed by s. No new constant. Same for Meier (replaces G4's local-Martin + taper).

## 6. Status of the design (2026-10-28, after the step-back with the user)
The user chose a BOUNDED heat-model check: P2-11 runs G1 only (input change). §5.6's path rule (X1) is a DRAFT for a
follow-up packet, written only if P2-11's e1 row says GO. Risks recorded at the step-back: shape tests the flame heat
field more than the rock; the uptake gap may carry spatial structure a uniform S hides; the width rests on one crater
number; the planner's last four mechanism hypotheses were refuted by the next step.
