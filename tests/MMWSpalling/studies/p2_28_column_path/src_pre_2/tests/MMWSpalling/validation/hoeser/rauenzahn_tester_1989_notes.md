# Rauenzahn & Tester 1989 — "Rock failure mechanisms of flame-jet thermal spallation drilling — theory and experimental testing", IJRMMS 26, 381–399 — reading notes

**On disk:** `franck_material/RockFailureMechanism.pdf` (19 pp). Read in full 2026-10-14.

- **Theory (Weibull over the thermally stressed volume, Eqs 9–16):** median spalling condition
  0.693 = (σ_s/σ0)^m (πC_L²/4)∫x²e^{−ωx}dx with ω = m u_r/α; ⇒ Eq. 13–14:
  T_s − T0 = [(Q/ρc_p)³ ((1−ν)σ0/(β_r E))^m (2·0.693/(πC_L²)) (m/α_r)³]^{1/(m+3)} — flux exponent
  **3/(m+3)** (0.13 at m = 20; 0.28 at m = 7.7). Aspect ratio from buckling: C_L = π/(6(1+ν)β_r ΔT_s)
  ≈ **10** (Eq. 15–16). Predicted ΔT_s 400–550 °C for Q 1–10 MW/m², ρc_p 2.64e6, E 45 GPa, ν 0.25,
  β 1e-5, σ0 70 MPa, m 20, α 1e-6.
- **Spall thickness distribution (Eqs 17–24):** G₁(x) in ξ0 = m V_dr x/α_r; median ξ0 = 2.9, mode 2.0,
  number mean 3.35, volume-weighted 4.85 ⇒ for V_dr 6 m/h, m 20, α 1.2e-6: x_med **0.10 mm**, mode
  0.07, mean 0.12, volume-weighted 0.17 mm. (Our memory's "2.67 ℓ/m" is this family; here 2.9.)
  Field spalls: average thickness **0.4 mm** (Canadian Mines Bureau), less than a grain.
- **Propane-torch tests (Barre/Westerly, 2,000 fps video, copper-block calorimetry):** torch tip
  25 mm, standoff 15/20 cm, O₂ 3.2/5.1 L/s; flux profile fitted as **Q = Q_max exp(−r/r0)** with
  Q_max 0.68 / 0.89 / 0.92 / 1.19 MW/m² and **r0 = 90 / 79 / 87 / 75 mm** (Table 3; measured Nu
  44–76 vs Kataoka 47–74). Spallation temperature RISE **465/473, 435/446, 445/423, 472/455 K**
  (Barre/Westerly) — i.e. T_s ≈ 700–750 K — with "quite small, if any" flux dependence (Fig. 9), and
  ALL below the 573 °C quartz transition. Spalls lamellar, 10–20 grains across, 1–2 grains thick.
  ⚠ r0/standoff ≈ 0.45–0.5 is a 25 mm subsonic torch, not a 7.9 mm jet; do not transfer the width.
- **Laser tests (500 W CO₂, 5–15 mm beams):** 520–1,120 K, inconsistent — "extremely localized
  heating … invalidates the 1-D assumption"; single grains lifted rather than spalls ejected.
  A warning for the beam phase: a spot comparable to the grain size does not spall like a flame.
- **Field rates (Table 1):** Browning, Conway NH: 15.8 m/h, Ø 0.2–0.25 m, 3.4 MPa chamber, 34 m³/min
  air; Barre VT: 7.6 m/h, Ø 0.35–0.4, 0.86 MPa; Los Alamos Pedernal: 6–7 m/h, Ø 0.35–0.45, 0.76 MPa.
  Air flows ~50× Meier's; context only.
- Grain-scale analysis (Eq. 6–8, Fig. 7): mismatch stresses σ ≈ EΔβΔT; a 40 °C rise near a free
  surface gives ~10 MPa tension around a quartz inclusion ⇒ heating creates microcracks beyond the
  room-temperature population ("mechanically determined Weibull parameters are not directly
  applicable", echoed by Wilkinson).
