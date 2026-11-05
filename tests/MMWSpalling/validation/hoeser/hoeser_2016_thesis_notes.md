# Höser 2016 — "Flame-jet assisted drilling technology" (Diss. ETH No. 23896) — reading notes

**Source on disk:** `franck_material/hoser.pdf` (194 PDF pages; printed page p = PDF page p + 20).
DOI 10.3929/ethz-a-010825044. Author Dragana Höser, supervisor Ph. Rudolf von Rohr, co-examiner
K. Boulouchos, ETH Zürich 2016. Read 2026-10-14 by the advisor agent: front matter, Ch. 3 (setup),
Ch. 4 (surface temperature + heat-transfer study), Ch. 5 (drilling), Ch. 6 (field-tool test),
Ch. 7–8. Ch. 2 (state of the art) and the appendices were not read. Numbers marked **[fig]** were read
off a plotted figure by eye (±5 % of axis span); everything else is quoted from text or a table.

⚠ **The rock is Grimsel granite/granodiorite from the Chuenzentennlen quarry at the Grimsel Pass** (p. 41).
Kant 2017 (JGR) took its "Central Aare granite" from the SAME quarry (Chüenzentennlen, p. 1810): the two
names are the massif and the local rock, one source. My 10-14 morning note calling the SGW 2017 caption
"wrong" is withdrawn (log 2026-10-14 §4). Our `Numeric/Material/Grimsel.H` is this rock.

## 1. Setup (Ch. 3, pp. 27–43)

- Methane–air, water-cooled chamber (inner Ø 80 mm × 140 mm), λ = **1.1** for all experiments
  (gave the highest heat flux; stable range 0.9–1.3). Firing rate **23–54 kW** (air line 7 bar limits).
- Two convergent nozzles: **D1 = 7.9 mm**, **D2 = 12.1 mm**. Chamber pressure 1.54–3.25 bar. D1 jet
  subsonic at 23/28 kW, **underexpanded (shock cells) at 33/38 kW** (transition 31 kW, p_ch 2.25 bar);
  D2 subsonic at all rates. Exit velocity up to ~850 m/s, Re up to 35,000 (p. 31).
- Flame temperature at the nozzle exit (Fig. 3.8, mean of two thermocouple methods) **[fig]**:
  D1: 1,270 °C (23 kW), 1,400 (28), 1,480 (33), **1,600 °C (38 kW)**; D2: 1,730 (38), 1,800 (54).
  ⇒ our Meier `jet_T_nozzle = 1900 K` (1,627 °C) at 38 kW is consistent with D1/38 kW.
- Dimensionless standoff L = S/D_nozzle, tested 5–14. Rock samples: cylinders Ø 100 × 150 mm.
- Table 3.1 (literature values used in her stress model): k 3.22 W/mK, ν 0.33, E 47.3 GPa,
  ρ 2,706 kg/m³, c_p 1,025 J/kgK, α 8.5e-6 /K. UCS 116.9 MPa, tensile 9.54 MPa (p. 101, ref [79]).
- Rock screening (p. 41): Grimsel, Bethel (USA), Gotthard, Paradiso (IND) granites, a gabbro and
  Rorschach sandstone; granite spalled most, sandstone melted. Grimsel chosen.

## 2. Spallation onset temperature on Grimsel (Ch. 4.2, pp. 65–69) — DIRECT CHECK FOR US

- Emissivity of Grimsel granite measured: **0.82** (own method, uncertainty up to 10 %; p. 64–66).
- IR camera (Optris PI450, 80 Hz, 60° to the axis), unconfined cylinder moved under the flame at
  0.2 m/s, L = 8, D1, 23–38 kW, 3 repeats; the temperature one frame before the first spall left:
  **430–600 °C**, Fig. 4.11 **[fig]**: ≈500 (23 kW), ≈545 (28), ≈460 (33), ≈555 °C (38 kW), error bars
  ±25–45 °C (frame-to-frame rise up to 45 °C at 38 kW). **"Any trend was not recorded"** with
  firing rate (i.e. with heat flux 4.3 → 7.4 MW/m²). First spall within 5 s; spalls > 1.5 mm visible.
- ⇒ Our reconstructed firing temperature on Meier, **852 K = 579 °C** (P2-2f/P2-3), and the Kant-deck
  onset 843.8 K = 571 °C, both sit INSIDE her band, at its upper end. Her scatter (±50 °C over a
  ×1.7 flux range) cannot distinguish a flux exponent of 0 from Kant's 0.23–0.31 (which would give
  ≈ +60–90 °C across her range).

## 3. Heat-transfer study (Ch. 4.3, pp. 69–98) — THE CALORIMETER, confined like a cavity

- Two water-cooled calorimeters (CU: surface 70–200 °C; **PSI: surface 250–600 °C**, i.e. at rock
  spall temperatures), point sensing element on the axis, guard ring, and a **confining ring Ø 40 mm
  × 35 mm tall** on top "to simulate conditions present in drilling experiments" (p. 72). Heat flux
  by water calorimetry; **h = q_sen / (T_f − T_surf)** with T_f the nozzle-exit flame temperature
  (Eq. 4.16). Uncertainty 7–45 % on q, < 20 % above 4 MW/m² (p. 82).
- **h does not depend on the target surface temperature** (180–550 °C) for L ≥ 9 (p. 89–90, Fig. 4.22)
  — the fact behind ROADMAP standing rule 17.
- **Confinement raises the flux** vs the unconfined jet and removes the shock-cell oscillation
  (Fig. 4.19, p. 87): confined 38 kW D1: 9.0 (L 5) → 3.7 MW/m² (L 14); unconfined 38 kW: 5.5 → 2.5.
- PSI, confined, D1 (Fig. 4.25 / 4.26) **[fig]** — the dataset for the jet closure:

  | L | q, 23 kW | q, 28 kW | q, 33 kW | q, 38 kW [MW/m²] | h, 23 kW | h, 38 kW [kW/m²K] |
  |---|---|---|---|---|---|---|
  | 5 | 5.5 | 6.8 | 7.9 | 9.1 | 6.5 | 9.0 |
  | 8 | 4.2 | 5.3 | 6.5 | 7.4 | 4.9 | 7.1 |
  | 10 | 3.3 | 4.3 | 5.2 | 6.1 | 3.7 | 5.7 |
  | 11 | 2.9 | 3.8 | 4.7 | 5.4 | 3.4 | 5.0 |
  | 13 | 2.2 | 2.9 | 3.7 | 4.3 | 2.5 | 3.7 |
  | 14 | 1.9 | 2.5 | 3.2 | 3.8 | 2.1 | 3.2 |

  Flux is linear in firing rate at fixed L (Fig. 4.27): slopes 0.235 / 0.214 / 0.180 / 0.124 MW/m²
  per kW at L = 5 / 8 / 11 / 14.
- Nozzle diameter (Fig. 4.28): at 38 kW, D1 gives **2.35×** the flux of D2 at L = 5 (= the area
  ratio); equal nozzle heat flux q_f = Q_f/A_nozzle (474 MW/m²: D1/23 kW ≡ D2/54 kW) gives equal
  target flux at small L. Table 4.6: from (38 kW, L 5, D1), −39 % flux needs L +120 %, or Q_f −38 %,
  or D +38 %.
- Jet force on the target 10–13.4 kPa (Fig. 4.30) — negligible against UCS.
- ⚠ **Magnitude vs the numbers the programme has used:** at L ≈ 11 (S = 87 mm, Meier's median
  standoff) and 38 kW, D1 (the closest analogue of Meier's 7.5 mm nozzle at 38 kW): **q ≈ 5.4 MW/m²,
  h ≈ 5.0 kW/m²K at the stagnation point**. Our Meier map at the axis reads 1,326 W/m²K unscaled
  (CTRL8 axis flux 1.2 MW/m²) and 985 after the P2-3 scale to Kant & von Rohr 2016 Table 4
  (766–1,009 W/m²K). Höser's calorimeter is a POINT sensor on the axis inside a Ø 40 mm ring; the
  radial profile was not measured (her Comsol assumed a Gaussian). See §7.

## 4. Drilling experiments (Ch. 5, pp. 103–138) — all Grimsel, D1 unless stated, 30 s runs

- Sample held between two flanges with a **Ø 40 mm opening for the jet**; uniaxial load 0–47 MPa via
  four bolts (Fig. 5.4–5.6; 47 MPa ≈ 2–3 km depth). **Without load the cylinder split (mode I) within
  10 s** and spallation stopped (Fig. 5.1); with load the cavity grew until its diameter reached the
  sample's (≈120 s at 47 MPa). Removal rate from mass loss / ρ; 3 repeats; std 5 % at 30 s.
- **Time series** (38 kW, 47 MPa, L = 8 initially, burner FIXED; Fig. 5.9–5.11) **[fig]**:
  removal rate 64 (15 s), 66 (30 s), 81 (60 s), 90 cm³/min (90 s); cavity depth and diameter linear
  in time with dL/dt = 1.8·dD/dt; at 120 s depth ≈ 80 mm, top width ≈ 75 mm (Fig. 5.10 profiles at
  15/30/60/90/120 s); "drilling velocity" 8.0 (15 s), 6.3 (30), 5.3 (60), 4.1 (90), 3.9 cm/min (120 s)
  = 4.8 → 2.3 m/h, falling as the standoff opens. After 30 s the cavity undercuts the flange opening.
- **Load series** (38 kW, L 8, 30 s; Fig. 5.12) **[fig]**: 26 (0 MPa), 52 (17.5), 54 (28.9),
  64 (37.4), 66 cm³/min (47 MPa) — ×2.5. Attributed to tensile stress below the compressive skin.
- **Flux series A, firing rate at L = 8** (47 MPa, 30 s; Fig. 5.18) **[fig]**:
  23 kW 17 cm³/min @ 4.3 MW/m²; 28 kW 27 @ 5.2; 33 kW 33 @ 6.4; 38 kW 46 @ 7.4.
- **Flux series B, standoff at 38 kW** (47 MPa, 30 s; Fig. 5.19) **[fig]**:
  L 8: 46 @ 7.4; L 10: 35 @ 6.1; L 11: 23 @ 5.4; L 13: 21 @ 4.3 MW/m².
  Both series on one curve: removal ∝ q_sen^(≈1.5–1.8) over 4.3–7.4 MW/m² (my fit of the read-off
  points; +72 % flux ⇒ +120–170 % removal).
- **Nozzle series** (same nozzle heat flux 474 MW/m², Table 5.2, Fig. 5.25) **[fig]**:
  D2/54 kW: 15.4 @ 2.6, 18.0 @ 3.4, 19.4 cm³/min @ 4.0 MW/m²; D1/23 kW: 7.0 @ 2.5, 11.0 @ 3.4,
  14.5 @ 4.0. Fig. 5.26: after 90 s both ≈ 30 mm deep, D2 cavity wider (≈ 60 vs ≈ 50 mm).
- **Spalls** (Fig. 5.27–5.30): sieve mass fractions at L 8 — > 50 % below 0.5 mm, 20–25 % at
  0.5–1 mm, ≈ 10 % at 1–2.8 mm, a few % above; laser diffraction of the < 500 µm fraction peaks at
  **200–300 µm**; 3–4 mm flakes exist; high-speed video (10,000 fps): spalls break up while
  buckling. ⚠ The "average spall size 200 µm" that Kant SGW 2017 uses as the failure-plane distance
  Δx is this particle-size peak — a FRAGMENT size, not a measured flake thickness. Her stress model
  assumes a 250 µm × Ø 3 mm "future spall" pad (Fig. 5.16).
- Her 2-D thermoelastic model (Patru / Meier et al. Matlab, uncoupled, constant properties, Table
  3.1; h from the calorimeter, T_ref = 1,866 K; Fig. 5.13–5.24): a thin compressive skin with a
  tensile region beneath; compressive stress at the surface +50 % from 23 to 38 kW; load raises the
  tensile stresses below the skin and barely changes the compression.

## 5. The field-tool test (Ch. 6, pp. 139–143) — A SECOND MEIER-TYPE ANCHOR

- Burner designed by **Thierry Meier** at IPE: Ø 80 mm double-walled stainless tube, 0.8 m long,
  water-cooled, **6 mm de Laval nozzle**, **58 mm feet**; 38 kW; Grimsel block 0.5 × 0.7 × 0.7 m
  (0.5 m high). **Drilled through 0.5 m in 22 min ⇒ ≈ 1.4 m/h** (text says "around 1.5 m/h");
  **hole Ø ≈ 95 mm, straight, slightly varying** (Fig. 6.3). ⚠ "At the very beginning the burner could
  not smoothly be placed in the hole and remained at one position for some time (around 10 min),
  thus the cavity below the burner was becoming wider than required" — the sawn block (Fig. 6.3)
  shows a wider top ~0.1 m: the pre-crossing lip, in the experiment.
- Spalls collected through a cyclone.

## 6. Conclusions she draws (Ch. 7–8)

Heat flux (not flame temperature) is the lever; fluid temperatures above ~1,000 °C are not needed;
flux up to 4 MW/m² drilled Grimsel successfully; tensile stress below the skin matters and uniaxial
load helps; her model does not predict removal or hole formation (outlook: flaw-based removal).

## 7. What this means for the programme (advisor, 2026-10-14)

1. **Onset check passed:** Grimsel spalls at 430–600 °C by IR; our 571–579 °C is inside.
2. **A same-rock rate-vs-flux series exists** (§4, 10 points over 4.3–7.4 MW/m² plus the load and
   time series), lab-scale, 30 s, under a Ø 40 mm flange opening at 47 MPa. It is a target for the
   floor criterion in a configuration without a body in the hole, but it needs (i) a radial flux
   profile (only the axis was measured) and (ii) the load, which our `confining.p` does not represent
   (hers is axial on the annulus around a free heated face).
3. ⚠⚠ **The flux level conflicts with the scale P2-3 applied.** For the D1/38 kW flame — nozzle and
   power essentially Meier's — the calorimeter reads h ≈ 5 kW/m²K and q ≈ 5.4 MW/m² at Meier's
   standoff (L ≈ 11), against Kant 2016 Table 4's 0.77–1.0 kW/m²K that the floor map was scaled to,
   and against the map's own 1.3. The two ETH measurements differ by ~5× and both are "the same
   lab's burner". Until the Table 4 burner/firing rate/standoff definition is reconciled with
   Höser's (nozzle, kW, confinement, T_f convention), **neither number should be called "the
   measured coefficient"**. Höser's is the better-documented of the two and is a stagnation-point
   value inside a confining ring.
4. **Meier-type anchor #2:** 0.5 m through in 22 min (≈ 1.4 m/h), Ø ≈ 95 mm, with a wider mouth
   caused by the burner sitting above the hole for ~10 min. Same rate band and diameter as the pilot.
5. **Spall thickness is NOT measured here**; the 200 µm is a fragment-size peak.

## 8. The journal version of the heat-transfer study (read 2026-10-14)

Höser & Rudolf von Rohr, *Exp. Therm. Fluid Sci.* 91 (2018) 166–174, DOI 10.1016/j.expthermflusci.2017.10.014,
on disk as `franck_material/HeatTransferStudy.pdf`. Same experiment and figures as thesis Ch. 4.3;
adds Table 2 (T_nozzle 1,269 / 1,397 / 1,481 / 1,593 °C at 23 / 28 / 33 / 38 kW for D1; 1,729 at
38 kW and 1,801 °C at 54 kW for D2; p_ch 1.54 / 2.00 / 2.54 / 3.25 bar), sensing element Ø 9 mm (Fig. 3),
calorimetric water 2.0 kg/min, max flux 10.5 MW/m² (CU). **No radial profile and no Gaussian widths**
— gap 1 in `2026-10-14.md` §2 is NOT closed by this paper; the inverse route (cavity aspect ratio) or
the impinging-jet literature remains. Confined vs unconfined (Fig. 7): confinement raises the flux
×1.3–1.6 at L 8–14 and removes the shock-cell oscillation.
