# Höser, Meier, Patru, Kant, Rudolf von Rohr 2018 — "An experimental study of the influence of uniaxial load on flame jet drilling", IJRMMS 106, 311–318 — reading notes

**On disk:** `franck_material/UniaxialLoading.pdf` (8 pp). DOI 10.1016/j.ijrmms.2018.04.007. Read in full 2026-10-14.
Journal version of thesis Ch. 5.3.1–5.3.2 (same data). Adds the exact drilling-leg conditions:

- **Rock:** Grimsel Granodiorite, quarry Chuenzentennlen near the Grimsel Pass; ρ 2,706 kg/m³,
  c_p 1,025 J/kgK, E 47.3 GPa, ν 0.33, α 8.5e-6 /K, k 3.22 W/mK. Cylinders Ø 100 × 150 mm.
- **Flame:** 7.9 mm nozzle, 38 kW, λ 1.1, p_ch 3.25 bar (underexpanded supersonic free jet, eight
  shock diamonds), **T_nozzle = 1,866 K** (two-method mean), **axial distance 63.2 mm** (L = 8), "the
  smallest feasible distance". Spallation onset **< 0.5 s** after impingement.
- **Confinement (Fig. 4):** rock between two flanges 200 × 200 mm; upper flange **20 mm** thick with a
  **Ø 40 mm** opening; a **15 mm steel disc** between flange and rock ⇒ the jet enters a Ø 40 × 35 mm
  channel above the rock — the SAME geometry as the calorimeter's confining ring (Ø 40 × 35), so the
  calorimeter fluxes apply to the drilling leg by construction. Four M16 bolts; torque 0–180 Nm ⇒
  0–47 MPa (Fig. 7, linear). At torque 0 the flange's own weight gives ~9 kPa.
- ⚠ **"0 MPa" means the flange is ON but untorqued.** Without the flange the cylinder splits (mode I)
  within 10 s and spallation stops (Fig. 1). ⇒ a simulation of the 0 MPa leg has a Ø 40 mm steel
  opening above the rock and a free heated face, no axial load.
- **Results:** removal (Fig. 10, 30 s, ×3 repeats, 5 % std): **26 / 52 / 54 / 64 / 66 cm³/min at 0 /
  17.5 / 28.9 / 37.4 / 47 MPa**. Time series at 47 MPa (Fig. 8): 64 / 66 / 81 / 90 cm³/min at 15 /
  30 / 60 / 90 s. Cavity sections at 47 MPa (Fig. 9) **[fig]**: 15 s depth ≈ 20 mm, half-width ≈ 22;
  30 s ≈ 30 / 25; 60 s ≈ 45 / 28; 90 s ≈ 60 / 30; 120 s ≈ 80 / 32 mm; dL/dt = 1.8 dD/dt; the cavity
  undercuts the Ø 40 opening (top width 44–64 mm). Drilling stopped at ~120 s when the cavity reached
  the sample's diameter and the thin walls broke.
- **Her stress model (Meier et al. 2016 code):** uncoupled 2-D axisymmetric thermoelastic, constant
  properties, **h = 7 kW/m²K applied UNIFORMLY over 0 < r < 20 mm** with T_ref = 1,866 K (the
  calorimeter's stagnation value spread over the whole opening — an assumption, not a measurement),
  adiabatic elsewhere, flange load as a normal stress on r > 20 mm of the top face, no-slip bottom.
  Result at 0.01 s: a thin compressive skin (σ_rr to −300 MPa) over a tensile region (σ_rr, σ_θθ
  +5 to +12 MPa) whose tension GROWS with the axial load (Fig. 13, 15); the compressive skin barely
  changes. Their explanation of the ×2.5: tension below the skin loosens grain bonds. ⚠ Not an
  effect our criterion can express: Kant 2017's confining term (Step 18 in our code) is a Poisson
  compression under full lateral constraint and lowers the onset by ~40–70 °C, which moves a rate by
  ~10 %, not ×2.5.
- **Gaps NOT closed by this paper:** no measured radial flux profile (uniform-over-opening is her
  model's choice; the energy bound from 26 cm³/min says a ~1.5 D spot, and a uniform 7 MW/m² over the
  opening would deliver 9 kW to 0.6 kW of flakes — impossible by conduction, so the real profile is
  concentrated and/or the cavity's growth cuts the flux within seconds); no 0 MPa cavity section.
