# Wilkinson & Tester 1993 — "Experimental measurement of surface temperatures during flame-jet induced thermal spallation", Rock Mech. Rock Eng. 26, 29–62 — reading notes

**On disk:** `franck_material/SurfaceTemperature.pdf` (34 pp). Read in full 2026-10-14. Data: `wilkinson_1993_tables.csv`.

**Setup.** Barre (VT) and Westerly (RI) granite; (a) Rauenzahn's propane–oxygen field torch
(supersonic, T_jet ≈ 2,830 K, R_noz 3.2 mm), standoff SOD = Z/R_noz 32–96; (b) a 25 kW CO₂ laser with
a 4:1 annular beam spun at 4–6 krpm, 5–15.8 cm diameter, flux 0.1–3 MW/m²; IR scanner 20 Hz,
240 × 512 zones, calibrated ±1 % at 560–600 °C; flame interference ≤ 2 %. Heat flux by (i) the
transient 1-D rise T_s = T0 + 2Q√(αt/π)/k (Eq. 10–11), ±35 % at high flux from the 50 ms sampling,
or (ii) the drilling-rate energy balance Q = ρc_p V_dr (T_s − T0) using Rauenzahn's quarry rates.

**The central distinction (Fig. 15 vs Figs 11–14), which the programme had been missing:**
- **Global ONSET temperature** (first spall anywhere on the heated surface) is **"relatively
  insensitive to the applied heat flux level"**: 400–600 °C over 0.2–2.7 MW/m² for both granites,
  within ~30 % of the mechanical-Weibull prediction (m ≈ 20, σ0 ≈ 70 MPa·m^{3/20}).
- **STEADY-STATE spalling surface temperature** (area-weighted stagnation-region mean while spalling
  continues) **rises with flux**: Barre 500 → 930 °C over 0.6 → 2.9 MW/m² (Table 3), Westerly 480 →
  740; the laser series agree (Table 5). Fitted Weibull m = 7.7 (Barre, σ0 42), 12.8 (Westerly, σ0 51)
  — the "m ≈ 8" our flaw-size exponent discussion used.
- Their mechanism for the difference (Fig. 17): at high flux a partly buckled spall **overheats**
  while stress build-up competes with relief in neighbouring soft minerals before ejection; the
  onset is set by the weakest flaw and is flat. Sioux quartzite (93 % quartz, no soft phase) shows
  **no** flux dependence at all (≈ 355 °C, 0.7–2.3 MW/m²) — consistent with that picture.
- ⇒ **Our criterion fires at the ONSET temperature of each column (`T_surf` at firing), and that
  quantity is measured flat by Wilkinson (Barre/Westerly), Rauenzahn (propane torch, Fig. 9), Kant
  2016 (Gotthard/Bethel) and Höser (Grimsel): 0.2–7.4 MW/m², five rocks.** The 0.15–0.31 flux
  exponent that P2-2d chased is a property of the steady-state surface mean, i.e. of flakes
  overheating before they leave. The model's analogue is the time-mean floor temperature during
  drilling (`T_bar_flake`, the plotfile surface-cell mean), not the firing temperature. A future
  check: compare the model's steady floor-surface mean vs flux to Table 3/5, not its firing T.

**Other numbers kept:** Table 2 quarry rates (Rauenzahn 1986): Barre 6.8 → 2.7 m/h, hole radius
3.2 → 8.9 cm as SOD goes 20 → 68; Westerly 6.3 → 3.4 m/h. Stanton numbers 1–4e-3 at SOD 20–30
(Fig. 16). Spalls "lamellar, 10–20 grains in diameter, one to two grains thick". IR-measured
steady-state T_s 710–730 °C on Barre under the torch (Fig. 9).
