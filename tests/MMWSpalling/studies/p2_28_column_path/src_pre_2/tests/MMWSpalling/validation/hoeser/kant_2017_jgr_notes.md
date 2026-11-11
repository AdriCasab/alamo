# Kant, Rossi, Madonna, Höser, Rudolf von Rohr 2017 — "A theory on thermal spalling of rocks with a focus on thermal spallation drilling", JGR Solid Earth 122, 1805–1815 — reading notes

**On disk:** `franck_material/JGR Solid Earth - 2017 - Kant - ….pdf` (11 pp). Read in full 2026-10-14. The LEFM Sp criterion our `sp_weibull` path implements.

- **Criterion:** plane-strain edge crack of length a in a half-space under a convective BC (h, T_fl)
  and a pressure p parallel to the crack (normal to the surface); σ_xx = −EαΘ/(1−ν) − pν/(1−ν)
  (full lateral constraint); K_I by Tada's edge-crack weight function; **Sp = K_I/K_Ic ≥ 1**.
  Reduced form (Eq. 14): Sp_red = 2A*√(a³/π)·Eα/(K_Ic(1−ν)λ)·(T_fl − T_i)·h + 2·1.763·√(a/π)/K_Ic·pν/(1−ν),
  A* = 1.763. Sensitivity (Fig. 4): crack length a dominates; E, α, h·ΔT, 1/λ linear.
  ⚠ Our Step 18 confining term (+pν/(1−ν) on the compression-positive σ_xx) is exactly this p term.
- ⚠⚠ **ROCK IDENTITY — corrects my note of this morning:** "Central Aare granite … samples Ø 85 ×
  150 mm were taken from **Chüenzentennlen close to the Grimsel Pass**" (p. 1810). Höser's "Grimsel
  granite/granodiorite" comes from the quarry in **Chuenzentennlen** (thesis p. 41; 2018 paper).
  **Same quarry.** "Central Aare" is the massif, "Grimsel granodiorite" the local name; Kant's Table 1
  cites Nagra Grimsel-test-site reports for its properties. ⇒ Kant 2017's onset data, Höser's drilling
  data and our `Grimsel.H` properties are ONE rock source. The Meier pilot block (Grimsel granite,
  thesis Ch. 8) is the same family.
- **Setup:** the 2016 10 kW burner at ~1,400 °C with Höser's flange (Ø 46 mm opening, Fig. 5), pyrometer
  1.5 mm / 1,000 Hz; critical temperature = **failure-plane temperature θ_FP after the spall leaves**.
- **Table 1 (literature, room T):** E **30–40 GPa**, ν 0.19–0.33, K_Ic **1.5 MPa√m**, **a = 20 µm**
  ("the largest initially present cracks"), α 7.8–8.2e-6. ⚠ Our Kant deck uses a0 = 8 µm and the
  Grimsel property tables (E 53.3 GPa at 293 K per the plotfile header; Höser's 47.3) — three sources
  span E 30–53 GPa, a 1.3–1.8× spread in the predicted onset ΔT (Θ ∝ 1/E). Worth a line in the
  properties reference; not acted on here.
- **Results (Fig. 7, Table 2):** measured θ_FP at 0 / 27 / 48 MPa: **553–694 °C (mean ≈ 630 ± 40)** /
  ≈ 575 ± 30 / ≈ 560 ± 35 **[fig]**; model (a = 20 µm): 475 (390–560) / 450 (365–535) / 430 (340–510) °C.
  Pressure LOWERS the onset (−70 °C measured over 48 MPa; model −45). Measured sits above the model
  because the first spall rarely detaches inside the 1.5 mm pyrometer spot.
- Our Kant-deck onset 843.8 K = **570.6 °C** (0 MPa) is at the low edge of their measured 553–694 °C.
