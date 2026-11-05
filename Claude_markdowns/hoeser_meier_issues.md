# What is wrong with the Höser and Meier tests — a reference for choosing the research direction

Compiled 2026-10-24 from the daily logs (`2026-09-15.md` … `2026-10-14.md`; the 10-14 file also holds the
10-19 and 10-24 entries, §15–§20), the study result files, `ACTIVE_STEP.md` (P2-2e completion notes) and
`meier_thesis_notes.md`. **"Log 10-13 §4"** means `Claude_markdowns/2026-10-13.md`, section 4. Every row
gives its evidence; a number with no pointer was measured on 2026-10-24 and the source file is named.

Issue IDs (D1, S2, C3, G1 …) are stable so later notes can cite them. Updated 2026-10-24 (evening) with the
modelling principle (§0b), the typed-geometry issues (§3.9), Höser's energy bound (E8) and the source check of
her thesis and 2018 papers (log 10-14 §21).

---

## 0. The answer in one paragraph

**In every case where the flame's heat input was actually measured, real rock removes several times less
rock per joule of flame heat than the model does.** Höser measured essentially Meier's flame (same lab,
38 kW, 7.9 vs 7.5 mm nozzle) with a calorimeter and got **3.0 MW/m² at the model's own axis standoff, and
3.7–7.2 MW/m² across Meier's real 65–100 mm**. The model's Meier floor gets **0.77 MW/m²** on the axis (S0)
— 4–9× less — and **still drills 1.4× too fast**. Rauenzahn (1986) saw the same gap: a copper calorimeter in
a mock hole read 3–4× more heat than the drilling rate implied. **Höser's own data show it with no model at
all (E8):** at the face temperature her infrared camera saw, her rock could have taken in at most 1.4–1.6 kW
over 30 s, while her calorimeter implies 7.9–8.9 kW arriving — 5–6.5×. Each time, the model has been brought
into line by **turning the heat input down**: Martin's cold-air correlation (D2b), the Kant-burner scale
(P2-3), a narrow Gaussian spot (P2-5). That compensation is why the rate looks right and the shape is wrong
in both tests. The failure law passes every test whose heat input was measured. **The error is most likely
in how much heat a spalling rock face actually takes in.** Energy conservation fixes what happens to heat
once it is absorbed, so the difference is in the uptake, not in the conversion. **The second problem is
structural: the model does not compute the hole's radius and shape. Hand-set numbers decide them (§0b,
§3.9).**

---

## 0b. The modelling principle (user, 2026-10-24)

**One unified model that reproduces the flux law and computes the hole's radius and shape by itself.
Neither may be set by a defined parameter.**

| | allowed as an input | must come out of the model |
|---|---|---|
| Hardware | the rig's measured geometry: nozzle, flange opening, burner body, skirt and slots where they exist | which rock faces see which gas, and how much heat each one takes in |
| Flame | the measured heat transfer at the stagnation point (Höser's calorimeter) | how that heat spreads over the floor and walls of the cavity as its shape changes |
| Rock | measured properties, and material constants calibrated **once** on a lab test and never per case (flaw size, relief R, constraint depth) | — |
| Results | — | **the flux law**: onset temperature flat with flux; no spall below a minimum flux (Kant); removal rising faster than flux (Höser: ∝ q^1.4–1.8). **The geometry**: the hole's radius and shape against depth and time, and which faces count as floor and which as wall |

Today every deck breaks this somewhere: §3.9 lists the hand-set numbers that decide the hole's radius or
shape. Each is a defect to remove, not a parameter to tune. A hardware dimension is legitimate only as
hardware: Meier's Ø80 skirt may block gas because it is steel, but it may not decide where the rock wall is.

---

## 1. The two tests in brief

| | **Meier 2017 (field pilot)** | **Höser 2016 (lab)** |
|---|---|---|
| Rock | Grimsel granite, 0.7×0.7×0.5 m block, 1 MPa lateral load | Grimsel granite (same quarry), Ø100 × 150 mm cylinder, 0 or 47 MPa axial load |
| Flame | CH₄–air, 38 kW HHV, λ 1.2, Laval nozzle 7.5 mm | CH₄–air, 23–38 kW, nozzle D1 = 7.9 mm, T_nozzle 1,866 K at 38 kW |
| Geometry | burner **inside the hole**: Ø80 skirt with 3 slots resting on three pads at r 28–40 mm, Ø67 water-cooled body, sealed water-cooled wellhead at the mouth | fixed burner above a **Ø40 mm flange opening** (steel lid over the cavity) |
| Kinematics | crane-hung, operator-advanced, "contactless", hook load monitored | burner fixed; standoff opens as the cavity deepens |
| Measured | rate (1.30 m/h time-mean; 1.5 m/h quoted), hole Ø85 (Tab. 8.2), volume 3.42 L by water filling, Fig. 8.8 photo profile | heat flux on the axis (PSI calorimeter, Ø9 mm disc) vs standoff and power; removal vs flux, load, time; cavity sections at 47 MPa (30 s: a 47 × 30 mm bowl, wider than the Ø40 opening) |
| **Not measured** | gas temperature, wall temperature, heat flux anywhere, cooling-water ΔT, slot width, wellhead bore | radial flux profile (fitted to the sensing disc + guard ring, but never reported), cavity shape at 0 MPa |

Sources: `meier_thesis_notes.md`; log 09-18c §1; `validation/hoeser/hoeser_2016_thesis_notes.md`; log 10-14 §1–§4.

---

## 2. What is NOT wrong (do not re-open these)

- **The failure law on measured inputs.** Kant onset 843.76 K (= 570.6 °C) sits inside Kant 2017's measured
  553–694 °C on our own quarry (log 10-14 §4.5). Kant 2016's minimum flux is reproduced by P2-2e: the cut
  lands at 0.27–0.44 MW/m² against his 0.3–0.5 (P2-2e `D_REF.txt`; the flux worked out at each firing).
  Onset is flux-independent over 0.45–7.4 MW/m² for this rock family, and the model fires at 856–864 K at
  6× Meier's flux on Höser (P2-5 RESULTS c6).
- **Energy bookkeeping.** The ledger closes (`ledger_err_solid` ≤ 3e-3 on production legs). The price
  per removed volume, 1.19–1.33 GJ/m³ on every leg, is ρcΔT arithmetic at the firing temperature
  (log 10-12 §14; P2-5 RESULTS).
- **The long list of numerical defects already fixed.** None of them is the cause: E1 beam leak, the
  step-face cascade (D2j-0b/D2k/D2r-1), the mask edge (D2t), the mouth funnel as a heating BC (D2v), the
  mesh-dependent surface temperature (P0/P0b), the 1e6 unit bug in the lateral criterion (P1f), the shift
  leak (P1g-2c), the pad-clearance enthalpy leak (P2-2f), the `rc_ellq` reduction (P2-2d).

---

## 3. The issues, grouped by where they live

### 3.1 How much heat the flame delivers — the LEVEL

| ID | Issue | Evidence | Status |
|---|---|---|---|
| **D1** | **The model's axis flux is 4–9× below Höser's measurement of the matching flame.** S0 axis 0.77 MW/m² (r < 16 mm, 150–350 s) at the model's axis standoff 114 mm (L ≈ 14–15). Höser 38 kW on 850 K rock: 3.0 MW/m² at L 14 (**4.0×** like-for-like; Fig. 4.26 reads 3.0, the CSV's 3.2 is high), ≈3.7 at L 13, ≈5.0 at L 11, ≈7.2 at L 8. Meier's real nozzle-to-floor is 65–100 mm (K3; log 09-30a §1) ⇒ L ≈ 8–13 ⇒ 5–9×. ⚠ Her sensor face ran at 680–855 K (not cold), so no cold-sensor correction applies. | S0 `run_jet_profile.csv` + `thermo.dat` (`jet_s_c` 0.114 m); `hoeser_2016_heat_transfer.csv`; log 10-14 §7 | **Open — the central conflict** |
| D2 | The original map level (1,449 W/m²K = Martin 1977 at Meier's mass flow) is a cold, subsonic air-jet correlation applied to a supersonic flame. Its axis flux (CTRL8, 1.2 MW/m² at 85 mm) was already ~4× below Höser's at that standoff. The programme's first "untuned rate match" (D2b, 1.36 m/h) used it, and also rested on a numerical drain that was later removed (D2x: the un-throttled rate is 2.1 m/h). | log 09-16b §2; log 09-24 D2x §1; log 10-14 §1.4 | Open |
| D3 | P2-3 scaled the floor ×0.743 to Kant 2016 Table 4 — a 3–10 kW threshold burner, not a drilling flame. Category error. S0/W1/B1L stand as runs; their floor level does not. | log 10-14 §3.1 | Premise withdrawn |
| D4 | Nozzle gas temperature unmeasured. 1,900 K is used; equilibrium adiabatic is 2,040 K; the 160 L/h cooling water bounds it at 1,450–1,850 K; ΔT instrumented but never reported; the igniter thermocouple sits upstream of a lifted flame. | log 09-22a (Ch. 4 review); log 09-18c §3.5 | Open; needs Meier's data |
| D5 | The jet's axial decay has never been validated. Core length 8 D_e (only 3.8 D in `fixed` mode). `exhaust` mode props the floor gas at ~1,504 K whatever the decay; the same deck in `fixed` mode reads 750 K at 100 mm. Kant's flame decays 2–3.5× more slowly than the model's cold-entraining jet, because a lifted flame burns in flight. | log 10-13 §2–§3; log 10-14 §6 | Open |
| D6 | Delivery is limited by the coefficient, not by enthalpy: 75 % of the jet enthalpy leaves the floor unused at ~1,504 K. CTRL8 put 20 % of HHV into the rock faces, against Meier's implied 10–12 %. | log 10-13 §1 | Fact |
| D7 | Under `exhaust` recirculation, cutting the floor coefficient does **not** reduce delivered power, because recirculation feeds it back (P2-3: rock power +35 %, entrained mass +56 %). | log 10-13 §13 | Fact; constrains any level fix |

### 3.2 Where the heat lands — the SHAPE

| ID | Issue | Evidence | Status |
|---|---|---|---|
| **S1** | **Höser's radial profile was fitted (Gaussian on her sensing disc + guard ring, thesis Eq. 4.18–4.20) but never reported. The model's Gaussian (σ = 0.5 D) was picked because it matched her volume.** Her 30 s cavity (47 MPa) is a 47 × 30 mm bowl, and "spallation occurred at the whole rock surface exposed to the flame jet" (p. 117). Against Goldstein et al. 1986 (L/D 8, scaled to her axis value), the Gaussian gives 392 / 7 / 0.03 W/m²K at r = 10/15/20 mm where the literature gives 6,370 / 5,908 / 5,456. A realistic profile delivers 7.7 kW over Ø40 ⇒ ~116 cm³ in 30 s at the model's price, **9× her 13 cm³**. | log 10-14 §20.3; P2-5 RESULTS (ladder) | **Open — the Höser twin of D1** |
| S2 | Meier's map is flat inside r = 18.75 mm (Martin is invalid below r/D 2.5), so it has no stagnation peak. On S0 the floor flux *rises* outward, from 0.77 at the axis to 1.0 MW/m² at r 30 mm, because the axis sits deeper (pit lead) than the pad ring. | S0 `run_jet_profile.csv`; log 09-30 §3 | Open |
| S3 | **The bore diameter is an input, not a prediction.** S0's bore is 79.3 mm = 2 × `foot_r_outer` (the tool stance). The wall sits inside the map's typed 40→56 mm taper. D2x: "no f can place the wall at 43–46 mm without also multiplying the floor rate". | log 09-24 D2x §3; log 10-13 §4–§5; `ARCHIVE_DONE.md` :12628 | Fact |
| S4 | The free-jet law is used inside a capped cavity (Höser). The closure entrains 293 K air; T_stag falls 1,866 → 1,447 K as the pit deepens from 63 to 153 mm, and floor power falls 71 % in 30 s. Her cavity is flange-capped, so the gas inside is recirculating flame gas. | log 10-14 §20.1 | Open |
| S5 | **Contradiction on the radial shape.** The literature says in-hole flux is PEAKED (Rauenzahn/Wilkinson CFD: down to 1/3 by 3 nozzle radii). D2ad says Meier's two-sided window needs a nearly FLAT floor (ring ≥ 0.85–0.9 of the axis). This may be the support rule (K1) talking, not the flow. | log 09-29 §13; log 09-30 §5.5 | Unresolved |

### 3.3 The cavity: body in the hole, skirt, slots, wellhead (Meier only)

| ID | Issue | Evidence | Status |
|---|---|---|---|
| **C1** | **The skirt.** The bottom 50–58 mm of the burner is a Ø80 tube (bore 56) with three slots, resting on the rock. The pad ring (r 28–40 mm) sits under 12 mm of steel and is heated only through the slots. Over that band the rock wall faces steel across a 2.5–6.5 mm gap. **The model lets gas leave over the whole circumference and heats the pad ring with the map. The slot width is not recorded anywhere.** | `meier_thesis_notes.md` L36–41; log 09-22a (refrigerated-pipe review §2); log 09-24 D2x successor; log 09-28 D2aa review | **Open — never modelled** |
| C2 | One wall-gas temperature (`jet_T_exhaust` ≈ 1,504 K) is used for every side face at every height (`MMWSpalling.H:2481`). Real gas cools against the water-cooled Ø67 body. The wall march (P2-3) moved the mouth ≤ 5 mm; the body sink (P2-4) is ~3× too weak, because its coefficient falls from 142 to 47 W/m²K as the gap widens. | log 10-13 §4–§5, §12, §17 | Open |
| C3 | **Wall-coefficient estimates span 25×** and none is measured: laminar annulus 14–40 (log 09-22a; log 10-06 §2 above the nozzle); Gnielinski 47 below the nozzle (log 10-06 §2); typed 169.3; wall-jet near the corner 100–170 (log 10-06 §8); P2-1's derived composite, p50 372 (log 10-13 §4). | as listed | Open |
| C4 | **Wall heating is needed, yet it runs away.** S0 (wall off): 79 mm, the tool stance; Meier: 85. W1/B1L (wall on): 115–119 mm. User ruling 10-19: wall heating is physical, so something is missing. | log 10-13 §17; user, this session | **Open** |
| C5 | The mouth formed inside a sealed, water-cooled wellhead, over two phases, and in phase 2 a vacuum reversed the gap flow. The model has an open free surface at 293 K. The mouth is weak evidence either way. | log 09-18c §3.1–3.2 | Fact; mouth not scoreable |
| C6 | The ~6 mm mouth plateau comes from the free-surface 142 W/m²K plateau, which has no standoff decay, before the burner reaches the original surface. | log 10-13 §12.3 | Open, small |
| C7 | The block's 1 MPa lateral confinement is not modelled. | log 09-18c §3.6 | Small |

### 3.4 How the burner moves (Meier only)

| ID | Issue | Evidence | Status |
|---|---|---|---|
| **K1** | **Meier's burner hung from a crane.** The operator set the advance, the hook load showed hanging vs resting, and intermittent pushes "kept the drill free". In the model's feet mode the burner rests on the highest rock under its pads, **so the descent rate equals how fast the pad ring clears.** | log 09-18c §1, §3.4; log 09-28 "what physics is missing" | **Open — a modelling choice** |
| **K2** | **The pad ring is where the burner jams, every time.** D2o-0: switch jams the burner, 1.49 → 0.11 m/h. D2aa PD: jammed from t = 0. P1: 48 % of floor frozen, collision at 177 s. **P2-2e SD: pads silent at 36 s, centre at 260 s, burner never moves.** Four different physics changes, one failure location. | log 09-22a; log 09-28 D2aa; log 10-01 §11; SD event log (10-24) | **Open** |
| K3 | The rate is one-sided by kinematics (the operator runs at or below the contact limit), but the Fig. 8.8 pit lead of ≤ 5 cm makes the standoff two-sided (≤ 100 mm). Three independent routes agree. | log 09-28; log 09-29 §1; log 09-30a §1 | Fact |
| K4 | Prescribed and feet modes give different deliveries: feet pins the standoff at 50 mm, where map h is 2.4× its value at 200 mm. Numbers from the two modes are not one series. | log 10-12 §12 must-fix 1 | Fact |
| K5 | Under prescribed 1.5 m/h the standoff opens to 170–223 mm, which Fig. 8.8 excludes. Cutting the flux uniformly stalls the pad ring before the axis reaches ≤ 100 mm, so the window is empty (D2ad). | log 09-29 §7, §13 | Fact |

### 3.5 What happens at the spalling surface — HEAT UPTAKE AND REMOVAL

| ID | Issue | Evidence | Status |
|---|---|---|---|
| **E1** | **The surface is a perfect thermostat.** It fires the instant it reaches ~850 K and never overheats: on Höser the exposed surface sits **83–135 K BELOW** the firing T. Wilkinson measured the steady-state surface mean *rising* with flux (Barre 500 → 930 °C) and put it down to partly buckled spalls overheating before ejection. The model ejects instantly. ⚠ Her IR shows the face cycling ~600–870 K at onset (it drops ~100 °C after each spall; thesis p. 65 — unconfined, first 5 s). | P2-5 RESULTS c12; log 10-14 §4.2 | Open — cannot close E8 alone: that needs a face averaging ~1,400 K |
| E2 | The price per volume is fixed by arithmetic at 1.19–1.33 GJ/m³ on every leg, so whatever heat arrives is converted. The floor runs ~90 % efficient. | log 10-06 §1; P2-5 RESULTS | Fact |
| **E3** | **Rauenzahn's calorimeter.** A copper calorimeter in a mock hole read 3–6 MW/m²; his CFD matched it, and it was 3–4× above the flux implied by the drilling rate. His candidate causes: higher surface T during drilling, adhered spalls / stagnant gas, chips in the flow, losses ≤ 25 %. **This is the D1 gap, seen independently in 1986.** | log 09-30 §5.3; log 09-30a §1 | **Literature support for E** |
| E4 | Spall blowing (ejected spalls disturb the boundary layer): Wilkinson −35 %, film theory up to −50 % at Meier's feed. Never implemented. | log 09-30 §5.6; log 09-30a §2.4 | Not built. ⚠ E8's load comparison says it acts the wrong way (more spalling ⇒ MORE uptake) |
| E5 | Cuttings as a heat sink in the gas: Meier 6.5 g/s, model 12.3 g/s, 3–5 kW at realistic equilibration. The feedback has the right sign: faster drilling ⇒ more cuttings ⇒ less heat. | log 10-13 §1 (direction B) | Parked |
| E6 | The 1.5 mm flake cap binds on 95–99 % of wall firings and 95–100 % of pad-ring firings. Kant's notes: the cap is "the minimum-flux rule in disguise, applied to the wrong variable". | log 10-14 §5; `kant_2016_minimal_bc_notes.md` §3.3 | Open; touching it risks K2 |
| E7 | Höser's load raises removal ×2.5 (0 → 47 MPa) at the same flame. Unloaded rock wastes much of its heat, and the model's `confining.p` is a different load. | log 10-14 §2, §20.3 | Open |
| **E8** | **Her own data bound her rock's uptake.** At a face ≤ 870 K, the flakes plus the most that conduction (a stationary face — an upper bound) and radiation can carry over 30 s is **≤ 1.4–1.6 kW** at 0 MPa (12.9 cm³). Her calorimeter's h (7.0 kW/m²K, Ø9 disc mean) implies **7.9–8.9 kW** over Ø40 (Goldstein to uniform) ⇒ **5–6.5×**. It closes only if the face averaged 1,400–1,450 K. At 47 MPa (33 cm³) the bound is ~3.2 kW ⇒ 2.5–2.8×: **the same flame put ~2× more energy into loaded rock**, so how much heat the rock takes in depends on how it spalls (adhered or partly detached flakes shielding the face fit this; blowing does not) | scratchpad `hoeser_bound.py`; log 10-14 §21; P2-7 row a1 | **Open — the headline; P2-7 measures it** |

### 3.6 The failure law

| ID | Issue | Evidence | Status |
|---|---|---|---|
| F1 | A compressive stress is fed into a mode-I opening formula; R = 0.531 (UCS/σ_el) is the price paid. The 8 µm flaw compares with measured microcracks of 50–500 µm. R, a and K_Ic act as one knob on the level. | log 10-07 §4; log 10-14 §9 | Known; step B target |
| F2 | **Constraint depth (P2-2e) calibrated to 2.0 mm.** The minimum flux is reproduced (h150 never fires; a deep-absorbing beam gives no spall). Kant's fast standoffs move +47–63 K but stay inside Kant 2017's measured 553–694 °C. Höser loses 60 % under the Gaussian spot (S1), and still 27 % under a realistic profile (P2-7 c5; confounded, that leg also absorbed 24 % more). Meier jams at the pads (K2). Standing rule 1's flaw refit was skipped because the reference was cold at onset at 6.5 mm; at 2 mm it is not. | `ACTIVE_STEP.md` completion notes; `studies/p2_2e_cref/`; SD event log (10-24) | Completed, negative |
| F3 | The property tables clamp at 873.15 K. 22.6–34 % of Meier floor firings sit above the clamp, with a tail to 1,395 K (above melting); no melt limit. | log 10-12 §14; log 10-13 §9 | Parked |
| F4 | E is 30 / 47 / 53 GPa across three sources for one rock: a 1.3–1.8× spread on the predicted onset ΔT. | log 10-14 §4.5 | Open |
| F5 | The rock is one grain (`microstructure.mode = expression`), with no grain-scale heterogeneity. | log 10-07 §4.4 | Known |

### 3.7 The data themselves

**Meier.** The 1.30 m/h time-mean and the quoted 1.5 m/h disagree (log 09-18c §3.3). 3.42 L over 0.5 m
gives a mean Ø of 93 mm, not Tab. 8.2's 85; the self-consistent sets are (2.47 cm³/s, Ø93, 1.3 m/h) or
(2.47, Ø85, 1.57) (log 09-15b §1). Fig. 8.8 widths are lower bounds, since the saw cut need not pass
through the axis (log 09-18c §2). A pause underream sits at ~10 cm (standing rule 2). There were two
phases with opposite top-flow conditions. Cooling-water ΔT was recorded but never reported. The thermally
altered band cannot be measured from the PDF (log 09-22a §1).

**Höser.** Her calorimeter's sensing disc is Ø9 mm on the axis. She fitted a Gaussian to the disc and the
guard ring around it (thesis Eq. 4.18–4.20), but never reported its parameters, so no radial profile is
published. There are no cavity sections at 0 MPa. Without the flange the 0 MPa cylinder splits within 10 s,
and 30 s is the longest 0 MPa run that survives (p. 117). Her 2018 model assumed h = 7 kW/m²K spread
uniformly over Ø40. On her calorimeter that would deliver ~9 kW, while her rock could take in only ~1.5 kW
(E8). That is the uptake gap, not evidence for a narrow spot (S1). Her samples are Ø100 × **150 mm** (P2-5
used 100 mm, P2-2e 200 mm). Her "200 µm" is a fragment size, not a flake thickness (log 10-14 §1).

### 3.8 Measurement traps that have already misled a reading

- **The Höser lateral criterion was switched off by scope (G1).** The inherited W1 key `side_face_h_r_min = 0.036`
  limits sideways breaking to r ≥ 36 mm, and her whole cavity lies inside 20 mm, so P2-5's "no undercut"
  (c13) was vacuous. With the scope removed, the walls get the typed 169.3 because the derived wall law needs
  Meier's burner body (G2). That gives ~0.11 MW/m², below Kant's minimum, so they still do not fire
  (log 10-14 §20.2).
- **Rows that pass on a stalled run.** P2-2e c8/c9/c11 "passed" on SD because nothing was removed. P1's
  frozen census cannot see a floor that froze all at once (P2-2e takeaway 2).
- **Whole-run vs windowed references.** S0's wall share is 16.9 % in 150–350 s, not the whole-run 12.1 %
  (P2-2e takeaway 3).
- **Reduction slips.** Cap share "0 %" was really 95–99 % (wrong cap, wrong column; log 10-14 §5).
  Quarter-domain ratio: 0.309 instead of 1.237 GJ/m³ (standing rule 3). Per-step vs cumulative columns
  (standing rule 7).
- **Mesh.** At 2 mm the constraint depth is one cell and reads 25 K warm (P2-2e c3). Every 2 mm rate is
  provisional (standing rule 12; the 1 mm control is parked).

### 3.9 Hand-set numbers that decide the hole's radius or shape (each breaks the principle in §0b)

| ID | Hand-set number | What it decides | Why it is wrong | What should replace it |
|---|---|---|---|---|
| **G1** | `side_face_h_r_min` = 36 mm: Meier's Ø80 skirt radius, lowered from 40 so the modelled 39.77 mm bore counts as wall (`ARCHIVE_DONE.md` L19963) | Which side faces are "floor steps" and which are "walls", and two consequences tied together. Floor steps get the jet's heat but may **not** break sideways (`Removal.H:673`). Walls may break sideways but get the wall coefficient (`MMWSpalling.H:2397-2415`) | It is one tool's radius, and Höser's deck inherited it by copying. Her whole cavity lies inside 36 mm, so its walls can never break sideways. The split exists for a real numerical reason: the grid turns a smooth floor into a staircase. Unheated steps drained heat (D2j-0b), and steps breaking sideways would remove whole cells instead of thin flakes. But a fixed radius is the wrong way to tell a step from a wall | Decide "floor step or wall" from the hole's own shape (step height or slope; the code already measures step heights, D2z), and separate "which heat" from "may it break sideways" |
| **G2** | The derived wall law's built-in Ø67 body, `WJ_R_SHIELD = 0.0335` (`MMWSpalling.H:5629`) | The wall's heat. The law measures the gap to Meier's burner body and runs only when the hole radius exceeds 33.5 mm + 2 cells (`:5957`); otherwise it falls back to the typed 169.3 W/m²K | Höser's cavity (radius ≤ 20 mm, no body in it) can never use it, so every wall gets 169.3 (`side_h_used_n` = 0 on the r_min = 0 diagnostic; log 10-14 §20.2) | A wall heat law driven by the cavity that exists in the run: with a body in the hole when the deck has one, without it when it does not |
| **G3** | Side faces take the heat of their column's **top**, not of where they are | Whether rock beside the cavity can be heated at all. Under Höser's flange the column top is shaded (h = 0.001), so its side face facing the cavity is shaded too, even when the cavity beside it is 30 mm deep and full of hot gas | Her cavity widens under the flange: Ø39 at 15 s, **Ø47 at 30 s**, ~Ø59 at 60 s, ~Ø72–75 at 120 s (Fig. 5.10, 47 MPa, digitised 10-24). In the model it can never pass Ø40. That is ≤ 15 % of her volume at 30 s, and most of it by 120 s | **Face-local heat:** each side face's heat evaluated at its own position, so the flange shades only the rock it covers. First version: the jet's map at the face's radius and depth, run as a bracket, because no one has measured the heat on a cavity wall |
| G4 | Meier's floor map: literal radii inside the `h_expr` string — the 40→56 mm taper, the 142 W/m²K plateau, the 60→76 mm cut (`MMWSpalling.H:10035`; S0 `argv.sh`) | Where the floor's heating stops, and therefore where the bore wall sits | The bore wall sits inside the typed taper (S3; `ARCHIVE_DONE.md` :12628). The string has no dx, so mesh refinement cannot move it | A heat distribution from a flow model or a measured profile, not a hand-drawn taper |
| G5 | `foot_r_outer` (Ø80 feet) in feet mode, with the 0.9-quantile support rule | The bore: S0's 79.3 mm = 2 × `foot_r_outer`, the tool stance (D2m: the cut tracks the stance to 2 mm) | The tool's diameter becomes the hole's diameter by construction (S3, K1) | Burner kinematics that do not pin the hole to the tool (direction E); the skirt as a gas barrier, not as the wall |
| G6 | `spall.flake_h_max` = 1.5 mm | The flake thickness on 95–99 % of wall firings and 95–100 % of pad-ring firings (E6) | A clamp standing in for the minimum-flux rule (Kant notes §3.3), acting exactly where the hole's width is decided | A computed flake thickness and a physical minimum-flux rule (the constraint depth's job, F2) |
| G7 | `side_face_factor_above` (0 or 1) with the typed `side_face_h` = 169.3 | Whether, and how hard, the wall above the nozzle plane is heated. 0 gives the tool-stance bore (S0, 79 mm); 1 with the wall law gives the runaway (115–119 mm) | A switch, not physics: the result jumps between two hand-set states (C4) | The face-local heat of G3, with a gas temperature that follows the gas up the hole (C2) |
| G8 | Staircase risers at r < `side_face_h_r_min` inherit their column's full map heat (`MMWSpalling.H:2374-2381`) | How much heat a sloping floor takes in. A rim that lags exposes riser area at the full stagnation heat, takes in more and catches up | P2-7 c2: the Goldstein profile puts 9 % less on the floor map, yet risers take 22.7 vs 12.9 kJ and the totals match (127.0 vs 127.5 kJ, `p2_7_hoeser_heat/DIAG.txt`), so the floor stays flat (0.89 vs her 0.61). At 10 % heat the bowl still forms (0.64, `whatif_uptake/`), so G8 flattens but does not forbid it | Riser heat from the face's own exposure, with a staircase never taking in more than the smooth surface it stands for (a mesh-convergence gate) |

---

## 4. The six conflicts that matter

| # | Conflict | One side | Other side | Where it points |
|---|---|---|---|---|
| **X1** | **Measured heat vs needed heat** | Höser: matching flame at 3.3–7.2 MW/m² (D1); Rauenzahn 3–4× (E3) | model drills Meier 1.4× too fast at 0.77 MW/m²; Höser deck would over-remove 9× (S1); her own rock could take in ≤ 1/5 of the calorimeter's heat (E8) | **heat uptake at a spalling face (E8, E3)** — not the face temperature (her IR), not the conversion (energy conservation) |
| X2 | Wall heating physical vs wall runaway | S0 wall-off = 79 mm tool stance (C4) | wall-on 115–119 mm; coefficient unknown to 25× (C3) | the cavity flow: skirt, slots, gas temperature (C1–C2) |
| X3 | Kant's minimum flux vs Meier's descent | the law reproduces Kant (F2) | Meier's pad ring dies at 36 s (K2) | pad-ring heating (C1) and/or the support rule (K1) |
| X4 | Peaked vs flat floor | literature: in-hole flux is peaked (S5) | Meier's window needs a flat floor (S5) | the support rule may be forcing the "flat" requirement |
| X5 | Rate right, shape wrong — in both tests | Meier rate within 1.4×; Höser volume 1.09× | Meier bore = tool stance (S3); Höser 19×90 mm pencil vs a 47×30 mm bowl (S1) | X1: delivery was turned down to fix the rate, so it lands in the wrong places |
| **X6** | **The principle (§0b) vs every deck** | the hole's radius and shape must be computed | today hand-set numbers decide them: a skirt radius picks floor vs wall (G1), Meier's body is built into the wall law (G2), the flange shades the cavity's walls (G3), a typed taper places Meier's wall (G4), the feet set the bore (G5) | **computed geometry (direction F)** |

**X1 drives X5, and probably X2 and X3.** If the heat input were set to measured values, the model would
over-remove. The size of that over-removal is the size of the missing surface physics. The pad ring
(X3) would then get several times its present heat and might stop jamming, but the burner would descend
far too fast unless the uptake is also fixed. **X6 is independent of X1:** even with the heat right, the hole's
radius would still be set by G1–G7, so the principle needs its own work.

---

## 5. Research directions this suggests

| | Direction | Resolves | Cost | Risk / what it cannot do |
|---|---|---|---|---|
| **A** | **Measured delivery, both decks.** Put Höser's measured axis level and a literature radial profile (Martin/Goldstein, the same law on both decks) into Meier and Höser. Record the over-removal as the target. | Sizes X1 on both tests with the same inputs; tests whether the pad ring still jams (X3) | Inputs only; minutes per leg | Rates will be far too high by design. Must be read as a measurement of the gap, not a result |
| **B** | **Spalling-surface energy balance.** Flake overheating before ejection (E1), spall blowing (E4), cuttings sink (E5) — each a mechanism that scales with removal rate | X1 directly; Wilkinson's rising surface T; Rauenzahn's 3–4× | Source work, days per mechanism | Each needs its own measurement to avoid becoming a fitted factor; Wilkinson's data are Barre/Westerly, not Grimsel |
| **C** | **Cavity CFD** (2-D axisymmetric RANS: jet, Ø67 body, Ø80 skirt with slots as a porous ring, Ø85–122 hole) | h(r) on the floor, h(z) on the wall, pad-ring heating, recirculating gas T (S1–S5, C1–C3) | 1–2 days per geometry; external tool; your group has the tooling (log 10-06 §8 route 2) | Slot width unknown (bracket it); no in-hole data to validate against except Höser's axis |
| **D** | **Carry on the agreed failure-law sequence** (softening → α(T) → regime map) on Kant, Höser and the beam, with Meier as a regression | The beam application, which needs it whatever else happens | As planned | Meier and Höser cannot validate it until X1 is settled |
| **E** | **Burner kinematics.** Replace the feet support rule with prescribed feed plus a contact check (Meier's crane), so the pad ring stops being the bottleneck by construction | K1/K2; separates "the ring can't clear" from "the burner can't move" | Inputs (`nozzle_descent = prescribed` exists) | Prescribed feed opened the standoff to 200 mm before (K5); that is X1 again |
| **F** | **Computed geometry** (§0b). Face-local heat for side faces (G3); "floor step or wall" decided by the hole's shape and separated from sideways breaking (G1); a wall heat law that follows the cavity in the run (G2); the skirt as a gas barrier, not the wall (G5) | X6; Höser's undercut and time series (Ø39 → Ø75 over 120 s); Meier's bore without a typed taper (G4) | Source work, off by default, one step per item; tested first against Höser's Fig. 5.10 sequence | The heat on a cavity wall is unmeasured, so it must be bracketed. The staircase problems behind G1 (D2j-0b; whole-cell removal) must stay solved, so the shape rule needs a mesh-convergence gate |

**Recommendation (updated 2026-10-24).** **A is now P2-7** (`ACTIVE_STEP.md`), on Höser only: it measures
E8 from her data, then runs uniform and Goldstein profiles at her measured heat on her real 150 mm sample,
with nothing fitted. If P2-7 confirms that the gap is uptake, **B** follows, built around E8's load
comparison (which favours shielding by adhered flakes over blowing). **F** is required by the principle
whatever P2-7 shows, and Höser's 15–120 s cavity sequence can test it. **D** continues in parallel for the
beam. Meier is the weakest place to judge any of this: four unknowns (slots, gas T, kinematics, wellhead)
sit between the flame and the rock, and five of the seven G-items act there. **Höser should be the primary
development case**: same flame, no body, no feet, measured axis flux.

**Update and decision (2026-10-26).**
- **P2-7 accepted.** Her data alone show the uptake gap (4.7–6.5× at 0 MPa, E8); the model takes in 2.7–3.0× her maximum.
- **What-if** (`tests/MMWSpalling/studies/whatif_uptake/`): at ~10 % of her measured heat the model gives her volume, the
  flux law (1.31–1.65 vs her 1.44–1.83), a minimum flux, the bowl (0.64 vs 0.61) and a rising rate, with no other change.
  The flux law and the shape follow from the uptake LEVEL.
- **Literature** (`validation/hoeser/literature_uptake_gap_notes.md`): the gap is old and unsolved (Rauenzahn 1986:
  "remains unexplained"; Dey & Kranz: 15 % reached the rock). Gas-side effects reach ~2× at most each; no MMW data exist.
- **Heat budget** (`tests/MMWSpalling/studies/heat_budget/`): fast drilling (Höser 47 MPa, Rauenzahn) is explained by
  Wilkinson's measured hot surface plus a modest gas side. **Slow drilling (Höser 0 MPa, Meier) is not explained by
  anything measured.** Stuck or overheated flakes (direction B) have no support at slow speeds, which overturns the
  10-24 recommendation of B.
- **DECISION (user):** heat uptake at a spalling face stays an **OPEN PROBLEM on the task list**. Until it is solved,
  each test's uptake level is a stated input taken from its own measured removal, and everything else stays a
  prediction. **Next: F (computed geometry), G1–G8, with Meier's bore Ø(z) as the shape test.** What would settle the
  open problem: Grimsel's surface temperature during slow drilling; a flux measurement on a spalling face; or C, a
  hot-jet flow simulation of Höser's ring and cavity.

---

## 6. Questions only new data can answer

- **Meier (one email):** cooling-water ΔT (bounds the gas temperature); slot width (VT5 drawing set, sheet
  15/17 — not yet read); wellhead bore and height; depth at the end of phase 1; whether 1.5 m/h is a
  within-phase slope.
- **Höser (emailed by the user; no answer as of 2026-10-26):** the k and a of her Gaussian fit (thesis Eq. 4.18–4.20), the only radial-width information her
  rig produced; any radial flux measurement; a 0 MPa cavity section.
- **Anyone:** a calorimeter reading on a spalling surface (Rauenzahn-style), which would split X1 between
  "less heat arrives" and "the heat that arrives is used less well".
