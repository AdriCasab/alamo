# Kant & Rudolf von Rohr 2016 — "Minimal required boundary conditions for the thermal spallation process of granitic rocks", IJRMMS 84, 177–186 — reading notes

**Source on disk:** `franck_material/MinimalBoundaryConditions.pdf` (10 pages). DOI 10.1016/j.ijrmms.2015.12.009.
Read 2026-10-14 by the advisor agent, in full. This is the paper behind **Table 4** in
`validation/kant/kant_von_rohr_validation_test_for_mmwspalling_plan.txt` (h_fl 766 / 1009 / 873 / 747 /
633 W/m²K at SOD 50 / 80 / 110 / 130 / 150 mm), i.e. the number the P2-3 floor scale (×0.743) was derived from.

## 1. What the experiment was

- **Burner: a 10 kW methane–air laboratory burner (Pharos E-354)** with a multi-orifice head
  (Fig. 3: 10 × Ø 5 mm holes on a Ø 15 mm circle, head Ø 28 mm). Flows **V̇_CH4 = 0.3 Nm³/h,
  V̇_air = 2.7 Nm³/h** (≈ stoichiometric; ≈ 3.3 kW HHV from the methane flow; the burner's nominal
  rating is 10 kW). Flame temperature ≈ **1,400 °C** at the nozzle, falling with SOD (Table 4: 1,450 →
  1,100 °C over 50 → 150 mm). SOD 50–160 mm. **Unconfined** flame on a flat block.
- Rock: **Gotthard granite (Gurtnellen) and Bethel White granite (Vermont)**, 250 × 250 × 150 mm
  blocks, air dried, free surface. Table 2: Gotthard λ 3.27 W/mK, c_p 792, ρ 2,594, κ 1.6e-6, E 54 GPa,
  UCS 214 MPa; Bethel 2.62 / 784 / 2,566 / 1.3e-6 / 58 GPa / 152 MPa.
- Surface temperature: two **Optris OPTCTL3MH pyrometers**, 1.5 mm spot, 1,000 Hz, λ 2.3 µm, at 30°;
  emissivity measured (Gotthard 0.549 ± 12 %, Bethel 0.434 ± 17 %); total uncertainty 5–6 % of the
  reading (Fig. 8). Resolves single spalls: θ_SP (just before) and θ_FP (failure plane, just after),
  Δθ ≈ 100 °C (Fig. 9).
- Heat flux: the rock replaced by a **sensor plate of the same dimensions** carrying a Wuntronic FCR
  heat-flux sensor at the stagnation point, water-cooled, low-conductivity steel surface so that
  θ_s,sens = 390–550 °C ≈ the rock's spall temperature; flame T by thermocouples 10 mm above the
  stagnation point (Pitts correction). h_fl = q_sens/(θ_f − θ_s,sens); uncertainty ≈ 17 %.
- Purpose, in the authors' words: "the minimal required boundary conditions were defined as the
  conditions where at least a single spall could be detected … A successful, efficient application
  of this process must be conducted at significantly higher heat transfer coefficient."

## 2. Results

- **Table 1 (spall / no spall):** Gotthard spalled at SOD 50–130 mm and **not at 140–160**; Bethel
  spalled at 50–150 and not at 160.
- **Table 4:** q_sens 720 / 866 / 755 / 609 / 450 kW/m², θ_s,sens 511 / 552 / 495 / 435 / 389 °C,
  θ_f 1,450 / 1,410 / 1,360 / 1,250 / 1,100 °C, h_fl 766 / 1,009 / 873 / 747 / 633 W/m²K at
  50 / 80 / 110 / 130 / 150 mm (maximum at 80 mm = end of the potential core).
- **Spall temperatures (Fig. 11):** means ≈ 600 °C (θ_SP) and ≈ 550 °C (θ_FP) for both rocks, with
  no clear SOD trend (Gotthard 620 → 570 °C over 50 → 130 mm; Bethel 630 → 540 over 50 → 150).
  Individual spalls scatter 550–700 °C (Fig. 10).
- **Minimal conditions (Fig. 13, 15):** θ_SP,min ≈ **500 °C** and h_fl,min ≈ **500 W/m²K**; the
  no-spall cases sit at h ≈ 600–700 W/m²K with θ_f 1,100–1,250 °C, i.e. a surface flux at the spall
  temperature of q ≈ h(θ_f − θ_SP) ≈ **0.3–0.5 MW/m²**. Gotthard (UCS 214) needs more than Bethel.
- **Fig. 14** (their points + Wilkinson Barre/Westerly + Rauenzahn): "the spallation temperature
  increases linearly with heat flux" over 0.35–2.4 MW/m² (regression ≈ 530 → 950 °C), with the
  caveat that the high-flux points are other rocks and another method.

## 3. What it means for the programme (advisor, 2026-10-14)

1. ⚠⚠ **Table 4 is not a drilling flame.** It is a 3–10 kW multi-orifice laboratory burner run at the
   threshold of spallation on an unconfined block. **The P2-3 scale (floor ×0.743 to 985 W/m²K at
   85 mm) scaled Meier's 38 kW drilling flame to this threshold burner — a category error, mine
   (advisor, 2026-10-13).** Höser 2016/2018 measured the matching flame (7.9 mm nozzle, 38 kW,
   confined) at **5.0–7.1 kW/m²K** on the axis at L 8–11; Table 4 is 5–9× below that because it is
   a different burner, not because of confinement or sensor averaging. The S0/W1/B1L runs remain
   valid as runs; the premise under their floor level is withdrawn. Standing rule 17's sentence
   ("the plate's coefficient is a flow property") was true but applied to the wrong flow.
2. **Onset temperature: flux-independent for Alpine granites over a 16× flux range.** Kant:
   550–630 °C at 0.45–0.87 MW/m² (Gotthard/Bethel, pyrometer); Höser: 430–600 °C at 4.3–7.4 MW/m²
   (Grimsel, IR), no trend in either. Our 571–579 °C is inside both. Wilkinson's rising trend is
   Barre/Westerly by another method; Fig. 14's "linear" regression mixes rocks and methods. ⇒ the
   model's fixed-temperature onset is supported up to at least 7 MW/m² on this rock family.
3. ⚠⚠ **A MEASURED MINIMUM FLUX FOR SPALLATION: ≈ 0.3–0.5 MW/m² (h ≈ 500–700 W/m²K at a ~1,200 °C
   flame).** Below it the surface can exceed θ_SP,min and still not spall (Fig. 2/15: the operating
   range is bounded by h_min independently of temperature). This is the rule the model lacks: the
   P2-3/P2-4 bore wall spalled at 0.03–0.16 MW/m², 3–10× below it, and the hole's outer floor at
   0.1–0.6 MW/m² straddles it, while the pad ring (0.9–1.6 MW/m²) is above it. A flux cut-off at
   ~0.4 MW/m² would quiet the wall without touching the ring. ⚠ The cut-off must be on the surface
   FLUX (ℓ_q = λΔT/q ≈ 2.75 mm at 0.5 MW/m²), not on the reconstructed layer thickness ℓ_eff: the
   ring's ℓ_eff is 8.5 mm from side preheat while its q is ~1 MW/m² (q_ratio 0.2 in the event log).
   Note the existing cap `flake_h_max` 1.5 mm = c_f·ℓ with ℓ ≈ 3.75 mm is numerically the same
   order as Kant's threshold layer — the cap is the minimum-flux rule in disguise, applied to the
   wrong variable and as a clamp instead of a cut-off. Rock-dependent (Gotthard > Bethel); no Grimsel
   value exists (Höser's lowest drilled flux was 2.5 MW/m²).
