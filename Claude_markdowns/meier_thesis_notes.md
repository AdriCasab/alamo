# Notes on Meier (2017), *Thermal spallation drilling* — what matters for the ALAMO Meier validation

**Source.** Thierry Meier, ETH Zürich Diss. No. 24021 (2017). Local copies:
`~/Downloads/ThierryMeier_PhD.pdf`,
`~/Desktop/researchDocs/thermalSpallation/raw/ContactlessDrilling_Meier_2017.pdf`.
**Printed page = PDF page − 27.** `pdftotext -layout` gives a clean text layer;
figures need `pdftoppm`. Burner drawings: thierrymeier.ch/documents/drillBabydrill_VT5_R4.pdf
(appendix A.2, p. 262). His FEM codes: FEM_2D_heatConduction.m / FEM_2D_thermoelastic.m (A.3).

Read on 2026-09-17: §2.4.1, §2.5, Ch. 3, §4.2, Ch. 6, Ch. 7, Ch. 8 (all), §9.1, Ch. 10, Ch. 11,
appendix. Written for the planner and for future reviews; page numbers are printed pages.

---

## 1. The pilot demonstration (Ch. 8, pp. 205–228) — the experiment we simulate

### Setup
- **Block:** Grimsel granite, 0.7 × 0.7 × 0.5 m, 660 kg, loaded symmetrically from the sides to
  SHmin = SHmax ≈ 1 MPa "to prevent axial splitting of the block rather than to support
  spallation" (p. 208). Wellhead glued to the top with mastic (p. 213).
- **Burner** ("contactless drilling tool", Fig. 8.3, p. 210): ≈ 850 mm long, 15 kg, hung from the
  hall crane on coil tubing, "keeping the coil essentially in tension using the drawworks. Its
  length ensures smooth and vertical progression of the drill, due to gravity" (p. 209).
  Chamber shield 5 mm tube; cooling mantle 1 mm tube; pressure wall 2.5 mm; internal chamber
  Ø 56 × 545 mm (1.34 L), rated 100 bar (p. 211). Axial swirler 62°, radial/axial velocity 1.88.
  Fuel injector: 12 holes Ø 1 mm at 30° (p. 212). Pure turbulent diffusion flame, no mixing
  spacers (p. 212). "The combustion chamber of the burner is also expressively oversized to
  ensure complete combustion of the fuel and prevent any high enthalpy oxidation in the annulus
  formed with the open hole" (p. 215).
- **Nozzle:** converging–diverging **Laval** nozzle (Zimmermann design); text gives outlet bore
  **7.1 mm** (p. 212); the drawings give **7.5 mm** (sheet 14/17); throat **6 mm** (p. 222, 225).
  ⇒ the jet is choked and leaves supersonic (our estimate: chamber ≈ 6 bar abs at this flow;
  the air line at 6 bar(g) "is limiting the combustion chamber pressure", p. 211; a manometer
  upstream of the non-return valve "allows to estimate the pressure in the combustion chamber",
  p. 208, but no value is reported).
- **Feet:** "burner feet are used to keep a distance between the nozzle and the exposed rock
  surface. In the present case, the feet height leads to a dimensionless stand off distance of
  7 — using the outlet diameter of the nozzle as reference — which generally leads to optimal
  heat transfer conditions" (p. 212). Drawings: the Ø 80 tube continues 50 mm below the nozzle
  flange, bore Ø 56, three exhaust slots ⇒ three arc pads on r = 28–40 mm; nozzle exit to foot
  tip 50 mm (sheet 1/17).
- **Flows:** CH₄ 2.47 kg/h, air 52 kg/h, λ = 1.2 (Fig. 8.7); 38 kW HHV (Tab. 8.2). Air from the
  building at 6 bar(g). **Cooling water 160 L/h** through the tool, K-type thermocouples
  upstream and downstream (p. 208) — **the temperature rise is not reported anywhere in Ch. 8.**
  (160 L/h = 186 W per kelvin of rise; 10–30 K would be 1.9–5.6 kW ⇒ 100–300 K off the nozzle
  gas temperature.)
- **Exhaust path:** phase 1: up the annulus between burner and wellhead tube, out the top;
  phase 2: a vacuum cleaner on the exhaust line reversed the flow in that gap, "no spall is
  dragged out of the system" (p. 219). Cyclone separates > 10 µm (p. 213).

### Procedure and the "relative temperature" (Fig. 8.7, p. 220)
- Fig. 8.7 plots T/(1163 °C), ṁ_CH₄/(2.47 kg/h), ṁ_air/(52 kg/h). "Relative" = normalised by the
  set point / maximum.
- The temperature is the **igniter thermocouple** (NiCr-80/20 custom sensor, Ch. 5), coil tip
  33 mm from the injector, "positioned in the central recirculation of the purely diffusion
  flame" (p. 212). "The sensor temperature is reported rather than the local fluid temperature,
  i.e. no correction is applied for radiation and convection losses" (p. 215). Calibrated in an
  oven to 1000 °C (Fig. 8.5, p. 217); ±15 °C systematic.
- Time history: ignition at ≈ 50 s with 0.11 kg/h CH₄; flows ramped stepwise over ≈ 10 min;
  at ≈ 1 kg/h CH₄ / 20.8 kg/h air "the temperature drops significantly, corresponding to the
  flame detachment from the fuel injection nozzle. Thereafter, the flame is lifted up" (p. 219).
  **The 1163 °C peak is the attached-flame reading during the ramp; during both drilling phases
  the sensor reads ≈ 0.4–0.5 × 1163 °C ≈ 500–580 °C** — it is upstream of the lifted flame and
  does not see the products. ⇒ **There is no measurement of the nozzle gas temperature.**
  1436 K must not be quoted as a measured or lower-bound gas temperature.
- "During drilling, the temperature log is relatively steady, despite the constant axial
  acceleration imposed to the burner by the drawworks, which is necessary to keep the drill
  free" (p. 219). "Important fluctuations … are all related to rapid variations in the air
  flowrate"; a sudden air-flow drop at t = 1250 s in phase 2 (p. 220). "The intermittent and
  rough axial displacement of the chamber while drilling did not influence significantly
  combustion" (p. 220).
- Drilling phases: 647 s (phase 1) and 736 s (phase 2), both after the flows were at set point.

### Results
- **ROP:** "determined by analyzing the burner position with respect to the top of the
  wellhead. The hook load (viz. crane scale) giving evidence whether the burner hangs or lies
  on the rock. During the demonstration, an average and constant penetration rate of 1.5 m/h
  has been reached" (p. 221). Position read from marks on the coil at 10 cm intervals (p. 215).
  Self-consistency: 1.5 m/h × 1383 s = 0.58 m > 0.50 m block ⇒ 1.3–1.6 m/h is the consistent
  band (our 09-15 analysis).
- **Volume:** "the total excavated volume (3.42 L) is estimated by filling the hole with water"
  (p. 221) ⇒ includes any mouth widening. 3.42 L / 0.5 m ⇒ mean Ø 93 mm.
- **TSE:** 15.4 J/mm³ = 38 kW × 1383 s / 3.42 L, HHV basis (Tab. 8.2, p. 222). Same order as
  Browning/Fleming field tests (8.7–60 J/mm³) despite 40× lower power ⇒ "it is likely that the
  penetration rate was limited by the power of the combustion and not by poor heat transfer"
  (p. 222).
- **Hole Ø 85 mm** (Tab. 8.2, D column, p. 222). "The hole shape is cylindrical and calibrated
  … only slightly underreamed allowing: the necessary clearance for the contactless burner to
  advance under the effect of gravity and, sufficient velocity of the exhaust stream to drag
  the cuttings out of the hole. At about 10 cm from the top, there is an underreamed section
  which was induced by holding the burner steady for several seconds instead of continuously
  advancing the drill" (p. 223). "The clearance to the hole is adjusted by changing the
  advancing rate of the drill" (p. 223). "As the diameter of the drill controls the smallest
  hole diameter, the diameter of the nozzle sets an upper limit to the surface heated by the
  concentrated jet and thus the maximal spall size" (p. 223).
- **Fig. 8.8 (p. 224):** sawn block, near-cylindrical hole top to bottom, local underream at
  ≈ 10 cm, thermally altered band beside the hole. **No collar funnel.** Fig. 8.9: optical
  televiewer log of the first 40 cm.
- **Central pit foreseen:** "increasing solely the combustion power using the same nozzle …
  is inefficient. This is because the advancing geometry of the hole would look like a half
  elongated spheroid" (p. 222). Recommends nozzle Ø = 10–30 % of drill Ø (pp. 223, 227).
- **Cuttings:** 82 % collected; sieved < 500 µm (72 % of mass) laser-diffracted; average
  < 100 µm, 90 % < 300 µm (Fig. 8.10, p. 226); smaller than the ≈ 2 mm grains ⇒ intra-crystalline
  fracture ⇒ thermal-stress spallation, pore pressure negligible (p. 225).
- Closing (p. 227): TSE "compares well with similar tests … performed in the field"; the lower
  ROP "is attributed to a lack of combustion power … and a poor utilization of the thermal
  energy"; "it is better to overestimate the length of the burner and under dimension the
  outer bore diameter of the nozzle".

## 2. What the thesis says about the spallation physics (§2.4.1, pp. 42–52)
- Mechanism (Fig. 2.14, p. 52): pre-existing flaws → compression in the heated layer makes
  them grow and combine by tensile fractures parallel to the compression → when crack length
  ≈ distance to the free surface the growth is unstable → buckling → flake ejection. Spall
  aspect ratio 10:1–15:1 (p. 48; also [186] on p. 223: diameter/thickness 8–15).
- Stress: σ = −Eα(T − T₀)/(1 − ν) in the heated layer (Gray, Thirumalai; Eqs. 2.5–2.6). Our
  criterion uses exactly this stress on a Weibull-drawn microcrack; the thesis's lineage is
  Rauenzahn/Dey–Kranz Weibull theory (Eq. 2.7), **m = 15–25 for granite** (Dey et al. 1985,
  p. 48). Our m = 20 sits in the middle.
- Fluxes: flaking needs "1–10 MW/m²" (p. 52, citing [75]); too low a flux gives a heat wave and
  axial splitting of unconfined cores (p. 52); Browning: heat should penetrate ≈ 3 mm for max
  ROP in granite; higher flux → "dusting" (p. 44). Thirumalai (quartzite, 1.7 MW/m²): gradients
  175–260 °C/mm, spall temperatures 150–255 °C (pp. 46–47).
- **Fig. 2.12 (p. 49), Rauenzahn's compilation:** measured surface-temperature rise at
  spallation for Barre/Westerly granite ≈ 300–900 °C over 0.6–1.8 MW/m², bounded by Weibull
  m = 15 (upper) and m = 25 (lower). Our ΔT_fire = 528 K at 0.4–1.0 MW/m² lies inside this band.
- Rauenzahn's rules for unconfined samples: heated area ≤ 10 % of the total area; heat flux
  high enough to reach spalling T before ≈ 10 % of the volume is heated (pp. 48–49). The block:
  hole area / top area ≈ 1.2 % ✓.
- Calaman: spallability ∝ diffusivity × expansion × grain size / compressive strength (p. 46).
- Incident radiation "accounts at most for 1 % of the overall incident heat flux" (p. 55).

## 3. Impingement heat-transfer measurements (Ch. 3, pp. 73–105)
- Convective calibration rig: process-air heaters (3 and 18 kW), **solid-stream nozzle
  Ø 5.33 mm**, air 4–25 Nm³/h at 400–700 °C, **SOD = 2, 4, 6** (12, 22, 32 mm), Vatell
  microsensors as reference (Tab. 3.2, p. 78).
- **Measured stagnation h = 0.4–1.6 kW/m² K** over that range (Fig. 3.12, p. 91; Fig. 3.14).
  Martin (1977) evaluated at one of his points (10 Nm³/h, 500 °C, SOD 4; Re ≈ 2.6 × 10⁴) gives
  ≈ 1.0 kW/m² K — same order. **Nowhere in the thesis is 5–10 kW/m² K measured**; those values
  enter only in Ch. 7 as an assumption borrowed from Potter Drilling.
- Stanton number falls with heater power because "there is more entrainment of the
  surrounding air owing to the increasing density and velocity difference between the hot
  air-jet and the cold surrounding air … cooling the hot air-jet intensively" (p. 79). Relevant
  to our free-surface funnel: a hot jet over a free surface dilutes fast.
- Radiative loss from a 400 °C sensor surface is 0.2 % of the convective flux (Tab. 3.4,
  p. 100); "the incident heat flux is fully absorbed".
- Invariant-h concept (h depends on fluid-side conditions only) used throughout (p. 88).
- Heat fluxes on the sensors reached ≈ 0.8 MW/m² at 400–440 °C surface (Fig. 3.18, Tab. 3.4).

## 4. His own numerical model (Ch. 7, pp. 187–204)
- 2-D axisymmetric FEM, uncoupled quasi-static thermoelasticity (coupling δ = 0.0012, p. 189),
  128 700 quads, Δt = 2 ms, ≤ 10 s, **no material removal**, T-independent properties
  (Tab. 7.1: k 1.5, ρ 2750, Cp 790, E 30 GPa, ν 0.3, α 8e-6).
- Geometry: lab core r = 42.5 mm, 150 mm long, in a confining assembly; heated patch
  **r < 20 mm** (Ø 40) with **Robin h = 10 kW/m² K, T_ref = 1000 K** — "rather conservative
  values to describe the powerful flame-jets commonly used for drilling [73 = Potter Drilling]"
  (p. 193). Adiabatic elsewhere.
- Findings: surface reaches 90 % of T_ref within 1 s (p. 197); a thin compressive skin with a
  **tension zone beneath**, whose magnitude grows linearly with the axial load (Fig. 7.8,
  p. 203); lateral friction with the confining ring adds to it; ⇒ spallation efficiency should
  **increase with depth** (pp. 201–204); unconfined cores split axially (mode I) — explains lab
  failures vs field successes.
- Verification by manufactured solutions, 2nd order (Fig. 7.2).
- **Outlook (p. 250):** "Formulate mathematically a failure criteria for the cleavage of a spall
  in order to analyze and model the rate of penetration of the drill. This requires
  additionally a moving boundary condition or an artifact, e.g. a binary cell property." Also:
  T-dependent properties; measure the impinging heat flux of the flames (never done).

## 5. Laboratory drilling attempts (Ch. 6, pp. 171–186)
- Hydrothermal (ethanol–water/O₂) flame-jets at 260 bar on Ø 83–85 × 150 mm cores; SOD 10–15
  mm; 86–132 kW; jet velocities 10–250 m/s depending on throat (2.5–10 mm).
- Unconfined: spallation initiates then stops; axial splitting; cooling water quenched the jet
  in pre-drilled cavities. With **axial confinement ≈ 73 MPa** (≈ 2.9 km overburden): a real
  cavity, ≈ 30 g in 10 min, no axial fracture (p. 184).
- "While the fluid pressure has relatively minor impact on spallation, the initial state of
  stress is very important" (p. 185). Mechanism classified as true thermal-stress spallation.
- Potter Drilling: up to 200 mL in 3 min with 41.3 MPa axial + 20.6 MPa radial load (p. 178).

## 6. Energy-balance method he used elsewhere (§4.2, pp. 119–123)
- For the pressure-vessel burner he estimates the jet temperature from a **mass and energy
  balance over the vessel including the cooling water** (Eq. 4.1, Fig. 4.11) and the jet
  velocity from ṁ/(ρA) (Eq. 4.2). This is the method to apply to the pilot burner for the
  nozzle temperature; the missing input is the 160 L/h cooling-water ΔT.

## 7. Context statements worth quoting
- Demonstration at ambient pressure with gaseous feeds "is more challenging than what is
  expected downhole" (p. 207); gaseous feeds are not for depth; liquid systems beyond ≈ 100 m.
- Lab tests before this "spalling only initiated … the process stopped after some time"; the
  open question was "if the process can be sustained or if the rate of penetration drops as a
  result of heating the rock volume [Silva]" (p. 207) — answered: constant 1.5 m/h (p. 221).
- Conclusions (p. 247): "the hole is underreamed and continuous heating is a sustainable
  drilling process."
- Techno-economics (Ch. 9–10): flame-jet coil-tubing rigs cannot drill casing shoes; best case
  for mining/blast holes, workovers, lateral slim holes.

---

## 8. Consequences for our project (corrections and confirmations, 2026-09-17)

**Corrections to things we had written**
1. **1436 K is not a gas temperature.** The reference sheet's "chamber T plateau ≈ 1163 °C" and
   our "1436 K lower bound" are wrong: the sensor is upstream of the lifted flame and read
   ≈ 550 °C during drilling. Nozzle gas T is unmeasured; 1900 K adiabatic is the upper bound;
   the 1436 K run is a generic low case. The report (§1 table, summary) now says so.
2. **Score the volume on the whole excavation**, not inside Ø 93: Meier filled the hole with
   water. D2b whole-top 5.68 cm³/s = 2.3× measured. My 09-16 scoring-addendum advice is
   withdrawn; RESULTS.md §0 addendum should be amended when D2b is archived.
3. **The mouth funnel is a real model error** (Fig. 8.8 shows none). Cause: free-surface wall
   jet never entrains ambient air; exhaust recirculation on from t = 0. Fix = near-surface
   ambient entrainment (Ch. 3, p. 79 supports strong dilution of hot free jets).
4. **Author: Thierry Meier.**
5. Burner "≈ 850 mm, 15 kg" (p. 209) vs drawing set "797 mm, ca. 14 kg" — immaterial, but quote
   the thesis figure when describing the tool.
6. "The burner is not pushed" is too strong: the drawworks applied intermittent axial pushes
   "to keep the drill free" and the operator set the advance rate; the hook load told whether it
   rested on the rock. The feet rule remains the right baseline; a stall = where the operator
   pushed. Report §2.6 reworded.

**Confirmations**
- Feet / 50 mm / Ø 80 / three slots: thesis p. 212 + drawings, consistent.
- 85 mm: Tab. 8.2. Volume 3.42 L by water filling. ROP constant across phases; exhaust routing
  and a mid-run air dip did not change it.
- m = 20 inside Dey's 15–25; ΔT_fire 528 K inside Rauenzahn's measured band; hole area 1.2 % of
  the top satisfies the 10 % rule; radiation ≤ 1 % of incident flux (we include it anyway).
- Meier's outlook asks for exactly the model we built (failure criterion + moving boundary).

**New items to carry**
- **Supersonic jet:** Laval nozzle, 6 mm throat, choked ⇒ Martin (subsonic) is being applied to
  an underexpanded supersonic impinging jet. Re is mass-flow based (unchanged); the impingement
  structure is not. Limitation, not a fix we can make from the rock side.
- **Cooling water 160 L/h**, ΔT not reported ⇒ nozzle-T energy balance needs the author or the
  design documents. Each 10 K = 1.9 kW ≈ 100 K off the gas.
- Property-temperature convention in Martin: ±7 % on h.
- Ch. 3's measured 0.4–1.6 kW/m² K (Re ~ 10⁴) is the only jet heat-transfer data in the thesis;
  it is of the order of the correlation value, not of the 5–10 kW/m² K guesses.
- Meier expected a "half elongated spheroid" face if power were raised at fixed nozzle
  (p. 222): the central pit our model produces is the effect he had in mind; the ring-limited
  ROP with a deep centre is physically plausible, its unbounded depth is our clamp artefact.
- The underream at ≈ 10 cm depth was caused by holding the burner still for seconds
  (p. 223): a direct datum for "hole widens when descent stops" — the model reproduces this
  behaviour qualitatively (columns keep firing until 50 mm behind the nozzle).

## 9. Quick-reference numbers (Ch. 8 unless noted)
| item | value | page |
|---|---|---|
| block | 0.7 × 0.7 × 0.5 m, 660 kg, 1 MPa side load | 208 |
| CH₄ / air / λ | 2.47 kg/h / 52 kg/h / 1.2 | Fig. 8.7 |
| power (HHV) | 38 kW | Tab. 8.2 |
| nozzle outlet / throat | 7.1 mm text, 7.5 mm drawing / 6 mm | 212, 222 |
| stand-off | 7 D ≈ 50 mm (feet) | 212 |
| burner OD at feet | 80 mm (drawing) | A.2 |
| chamber | Ø 56 × 545 mm, 1.34 L | 211 |
| air supply / cooling water | 6 bar(g) / 160 L/h | 208 |
| igniter TC | peak 1163 °C (attached flame); ≈ 550 °C while drilling; uncorrected | 215, 220 |
| drilling time | 647 + 736 = 1383 s | 220 |
| ROP | 1.5 m/h "average and constant" (consistent band 1.3–1.6) | 221 |
| volume | 3.42 L by water filling | 221 |
| TSE | 15.4 J/mm³ | Tab. 8.2 |
| hole Ø | 85 mm; underream at ≈ 10 cm from a pause | Tab. 8.2, 223 |
| cuttings | 82 % collected; avg < 100 µm; 90 % < 300 µm | 224–225 |
| granite | k 1.5, ρ 2750, Cp 790, E 30 GPa, ν 0.3, α 8e-6 | Tab. 7.1 |
| Weibull m (granite) | 15–25 | 48 |
| Ch. 7 model BC | h = 10 kW/m² K, T_ref = 1000 K, Ø 40 mm patch, ≤ 10 s, no removal | 193 |
| measured impingement h (air jets, Ch. 3) | 0.4–1.6 kW/m² K | Fig. 3.12 |
