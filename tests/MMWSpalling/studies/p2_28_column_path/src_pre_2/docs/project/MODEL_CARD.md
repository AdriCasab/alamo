# Model card — the MMWSpalling thermal-spallation model

**What this document is for.** It answers one question: *what would we need to know to run this model on a test
nobody has done yet?* Every number the model uses is listed with where it comes from. Numbers fitted to a test's
own drilling result cannot be used in a prediction, and this card is how we find them.

**Written** 2026-10-30, against the binary built from today's source.
- The model card itself was built on sha256 `9278aa02d8036a74…`, with X1's path mode removed (P2-12).
- The user then restored X1 the same day, default-off. Binary `0a96a095ad30f45e…`; every default output is unchanged
  (`studies/p2_12_retire_path/restore/`).
- **P2-13 (2026-10-31):** the rock and the flame model form are now two shared input files
  (`tests/MMWSpalling/validation/shared/rock_grimsel.inputs`, `model_form_flame.inputs`), and the Kant calibration
  deck runs in the production form (`validation/kant/grimsel_calibration/`). The banked constants hold there:
  onset **850.43 K**, inside the Grimsel band 826.15–873.15 K (§9). **Regenerate** the raw tables with `python3 docs/project/model_card/raw_keys.py` and
`python3 docs/project/model_card/runnable.py` (from the repo root). **Update this card** whenever a key is added,
removed, re-calibrated or moved between classes; a step that changes a number in a calibration or validation
deck owes a line here.

**Goal it is judged against** (user, 2026-10-30): an *academic*, generalisable model of spallation and
vaporisation that predicts *unseen* MMW tests. Flame drilling is the proving ground because flame data exist and
MMW removal data do not (yet). Höser and Meier are checks, never targets to tune to.

---

## 1. Summary: where the model stands

1. **The rock half is close to general; the heat-input half is not.**
   - Rock properties are measured (Grimsel reference).
   - The onset level is set by two constants calibrated once on Kant: an 8 µm flaw and a relief factor of 0.531.
   - The flame heat input, by contrast, is typed per test: Höser's calorimeter fit and Meier's hand-shaped map.
     Both carry a heat multiplier S fitted to that test's own removed volume.
2. **The model is currently two models.** Today's removal rule (the face list, P2-8/P2-9) runs flame drilling only.
   It refuses vaporisation, melt-aware spallation, Meier's burner-clearance rule, mesh refinement, and any heat
   source other than a flame map (`MMWSpalling.H` ~984–1021). The constraint depth was deleted from both paths at
   P2-16 (2026-11-03); its key now aborts. Those features live only
   on the older column path. **No MMW case can run on today's removal rule.**
3. **The calibration decks do not run the production model.** The Kant decks that set the flaw and the relief
   factor use:
   - a single fixed flaw;
   - the `cell` surface heat form;
   - no flake cap;
   - `flake_length = flux` (the key deleted in P2-26; it now aborts at any value but `twocell`);
   - 1 MPa confinement.

   The Höser and Meier decks use:
   - random (Weibull) flaws;
   - the `surface` form;
   - a 1.5 mm flake cap;
   - `twocell` flake length;
   - 0 or 1 MPa confinement.

   The calibrated settings themselves (Grimsel law, relief factor) exist only as command-line additions inside
   study folders. The shipped Kant input files still carry the pre-calibration constants.
4. **Only three flame configurations run on today's binary:** the Kant decks, today's Höser deck (G1), and the
   Meier base inputs. The last quotable Meier deck (P2-2f) and both Höser 0 MPa decks (P2-7, the uptake what-if)
   abort on keys P2-10 deleted. **Today's model has never been run on Meier or on Höser 0 MPa.** (§5)
5. **Every Höser depth, shape and growth comparison so far used her 47 MPa loaded series, but the deck has no
   load.** The code cannot represent her load mechanism. S = 0.1636 absorbs both the heat-uptake gap and the
   missing load effect. Only her 0 MPa test matches the model's physics, and it measured a single 30 s volume.
6. **One unexplained factor dominates prediction: the heat-uptake gap.** Rock spalling under a flame takes in 3–6.5×
   less heat than a calorimeter implies. It is measured by three independent groups and explained by nobody.
   Until it is explained, any rate prediction carries one global factor calibrated on one test. Whether MMW will
   show it is unknown: beam data would settle it.
7. **Validity is narrow in places that matter.**
   - The property tables are flat above 873 K, and 22–34 % of Meier firings sit above that.
   - The jet correlations are valid only to H/D 12, which deep holes pass after about 45 s.
   - The flat floor drifts 10.4 % between 2, 1 and 0.5 mm meshes.
   - (Numerics: §10, with the MMW part in §8.)
8. **The MMW side is validated thermally only, on one case, with three surface settings tuned away from the
   reference.**
   - Zhang/Oglesby runs R 0.2, ε 0.8 and h 51; its reference says 0, 0.7/0.3 and 56.
   - It runs on a path that cannot remove rock.
   - It was last run 2026-05-30 and is not in the tree gate.
   - The absorption coefficient (300 m⁻¹) has no source on disk.
   - Temperature-dependent absorption is not energy-conserving.
   - The surface reconstruction the criterion fires on returns "flat" when the beam's temperature peak is below the
     surface.
   - Spallation has never been compared with data with the beam on. (§8)

---

## 2. The model in one page

| piece | what it does | where |
|---|---|---|
| Heat conduction | enthalpy form; per-cell properties from the phase law (`grimsel` or `constant`) | `MMWSpalling.H` `AdvanceMicrostructure` (~1744), properties `UpdatePropertyFields` (~4212) |
| Flame heat input | Robin flux q = h(T_gas − T_s) on every exposed face; h from the deck's `h_expr`; T_gas from an enthalpy march outward from the stagnation point (`jet_closure = enthalpy`) | `BuildFlameColumns`, `JetEnthalpyMarch` (~2412), `FaceListHeat`/`FaceH` |
| MMW heat input | Gaussian beam, Beer–Lambert volumetric deposition | `Numeric::MMWBeam`; §8 |
| Surface losses | radiation (ε, to T_amb) and natural convection on every face | `SurfaceCellFlux`, `FrFaceFlux` |
| Surface temperature | reconstructed from a two-cell exponential profile (P0) | `ReconstructSurfaceT` / `ReconFromX` |
| Spall criterion | Sp = K_I/K_Ic ≥ 1; K_I from the stress along a flaw of size a (Tada edge crack); stress σ = R·E(T_s)·∫α dθ/(1−ν) (1-D confined), plus pν/(1−ν) for confinement; flaws per cell from a Weibull draw (or one fixed size) | `Removal.H` (`DepthScanH`), `Grimsel.H` (stress), `FaceRemoval.H` R2 |
| Removal (face list) | one line per surface cell into the rock; the flake thickness is c_f·ℓ, capped; the solid fraction drops and the heat is remapped | `FaceRemoval.H` R0–R3 |
| Removal (column path) | the older top-down path: vapour front, bit, melt-aware, feet clearance | `Removal.H` |
| Burner support | Meier: the burner rests on its feet at the highest rock under its pads | `FeetDescent` (~2217) |

**Single shared rule:** the stress in every site comes from `Grimsel.H`, one function (`unit/lateral_flux_lock`
locks this).

---

## 3. What each experiment measured

| test | rock | what was measured | uncertainty | what was NOT measured |
|---|---|---|---|---|
| **Kant 2017** (JGR 122) | Central Aare = Grimsel quarry | onset (failure-plane) temperature 553–694 °C at 0 MPa; ≈575 °C at 27 MPa; ≈560 °C at 48 MPa (onset falls with confinement) | pyrometer 1.5 mm spot; first spall often missed ⇒ upper anchor | the heat-transfer coefficient (h) |
| **Kant & von Rohr 2016** (IJRMMS 84) | Gotthard, Bethel (other rocks) | spall / no spall vs standoff (Gotthard: spalls ≤ 130 mm, not ≥ 140); h_fl and T_flame on a sensor plate (747 W/m²K, 1250 °C at 130 mm); minimum flux 0.3–0.5 MW/m²; spall T ≈ 550–600 °C | h ±17 %; a 3–10 kW threshold burner (NOT a drilling flame, SR 17) | Grimsel |
| **Höser 2016 thesis + 2018 papers** | Grimsel | **calorimeter** h and q vs standoff (L 5–14) and power (23–38 kW), confined Ø40 ring; IR onset 430–600 °C, flat with flux; **drilling**: load series 26 (0 MPa) → 66 cm³/min (47 MPa) at 30 s; time series 64/66/81/90 cm³/min and cavity profiles to 120 s at 47 MPa; flux series A (power) and B (standoff) at 47 MPa; nozzle series (D2/54 kW vs D1/23 kW); fragment sizes | q 7–45 % (< 20 % above 4 MW/m²); removal ±5 % (3 repeats) | the radial heat profile (point sensor on the axis); flake thickness; shapes at 0 MPa |
| **Höser thesis Ch. 6: field-tool test** | Grimsel block 0.5 × 0.7 × 0.7 m | Meier-type burner (6 mm Laval nozzle, 58 mm feet, 38 kW): 0.5 m in 22 min ≈ 1.4 m/h, hole Ø ≈ 95 mm, wider top from a 10 min pause | text-level | heat transfer; shape detail. **Never used by this project** |
| **Meier 2017 pilot** (thesis Ch. 8) | Grimsel block 0.7 × 0.7 × 0.5 m | advance 1.28–1.5 m/h (operator-set, but the floor led the burner by ≤ 5 cm, so effectively two-sided); Ø 85 ± 5 mm over 150–450 mm (Tab. 8.2, Fig. 8.8); volume 3.42 L (includes mouth and pause); TSE 15.4 J/mm³; particle sizes; 1 MPa lateral pads | profile read from a greyscale figure (±5 mm) | nozzle gas temperature (≤ ~1900–2040 K); any heat-transfer coefficient (only cold-air calibration jets, 0.4–1.6 kW/m²K) |
| **Hu 2018/2019** | other granite, sandstone | lowest spall surface temperature 773–821 K (granite) under a prescribed surface temperature ramp; onset times 26–37 s | — | flaw sizes; our decks run a different failure law here (`damage_law`) |
| **Zhang 2023 / Oglesby 2014** (MMW) | granite | surface temperature under the beam, 137 GHz radiometer, piecewise power to 3.2 kW | emissivity-corrected (factor 0.7) | any removal: **thermal only** |
| Rauenzahn 1986 / Wilkinson 1989, 1993 | other rocks (Barre, Westerly, …) | flame spallation rates; surface temperature vs speed; in-hole copper calorimetry (Q(s) = Q_max e^(−s/s0), s0 ≈ 9–16 R_noz) | fits poor by the authors' account | Grimsel |

Notes in `tests/MMWSpalling/validation/{kant,hoeser,meier,hu,zhang}/`.

---

## 4. Proposed roles of the tests

⚠ **Proposed 2026-10-30, pending the user's confirmation.** The principle: each number is set once on one test,
then frozen, and a different test is predicted.

| role | test | what it sets or scores |
|---|---|---|
| calibrate the rock (done) | Kant 2017 onset; Kant 2016 min-flux bracket | the flaw (8 µm), the relief factor (0.531); the minimum-flux mechanism is the buckling span (row 28; the constraint depth was retired at P2-16) |
| heat-source inputs (measured without rock) | Höser's calorimeter; Meier's nozzle, flow and geometry | heat level and profile. Never fitted to removal |
| calibrate the unknown (new) | **Höser 0 MPa, 30 s volume** | the one global uptake factor |
| blind prediction | **Meier**; **Höser field-tool test** | rate first; diameter with a declared uncertainty |
| report only | Höser 47 MPa series | the load ratio (2.5×), a known gap until load is modelled |
| rock-only checks | Hu; Kant 2017 pressure series | the criterion's form at a prescribed temperature; the confinement term |

---

## 5. Which decks run on today's binary

From `docs/project/model_card/RUNNABLE.txt` (each deck started for one step; regenerated at P2-16 on sha `38c86238`; the only change from P2-13's `0a96a095` is the two P2-2e Kant legs).

| deck | result | why |
|---|---|---|
| **Kant `grimsel_calibration/input` (P2-13)** | runs | **the calibration deck**: includes both shared files; production form, face list |
| Kant `input_model`, `input_weibull` (shipped, legacy) | runs | — but these carry the pre-calibration constants (§6 row 1); not re-pointed (their retirement needs the user's permission) |
| Kant P2-2e legs (V80, G2.0_130) | **aborts** | they set `spallation.sp.constraint_depth`, retired in P2-16 (records; not edited) |
| Höser G1 (P2-11) | runs | today's Höser deck |
| Höser `input_lab` | template | the harness supplies `amr.n_cell`; since P2-13 it includes both shared files (the rock and model form no longer come from argv); G1's argv over it reproduces G1 byte for byte to 15 s (P2-13 a3) |
| Höser 0 MPa: P2-7 HG, what-if s100 | **aborts** | 13 keys deleted in P2-10 (`spall.lateral_*`, `side_face_*`) |
| Meier P2-2f (last quotable), P2-2e SD | **aborts** | 14 / 16 deleted keys; SD also sets `side_face_h_scale`, which nothing reads |
| Meier `input_feet`, `input_feet_d2c`, `input_drilling` (shipped) | runs | column path; no wall heat since P2-10, so not the P2-2f model. Since P2-27 the face list carries the pad clearance (regime 4, `FrFootClearance`): `input_feet_face` is `input_feet` re-keyed onto it (input_feet's own constants, NOT the rock file -- item 5), the face-list Meier feet deck; the three stay column-path records |
| Hu granite 2 (prescribed patch, 638.22 + 3.74 t K) | re-expressed at P2-27 (a RECORD, `validation/hu/hu_spall_onset/test_lefm`) | Grimsel's criterion under Hu's ramp on site A: first firing **107.9 s at T_B1 1042 K** vs Hu's 37.0 s / 773 K (pre-registered 50-66 s / 826-885 K: OUTSIDE). The side is the PATH: with a prescribed surface temperature site A takes K_Ic at the boundary T but the stress from the 3.125 mm top CELL's average, 790.7 K at onset -- 250 K below the boundary. Re-expressed on the face list at H (decision 8), where the reconstruction should remove the lag. `hu_end_to_end`, the damage-law onset decks and `diag/` retired (`tests/MMWSpalling/retired_validation/hu_damage_law/`) |
| Zhang (MMW) | runs | thermal only |
| **Oglesby `validation/oglesby/input` (P2-23)** | runs | **the first beam deck on the face list**: includes both shared files with the surface patch OFF (`surface_patch.enabled = 0`): the list is built and every listed face is billed at h = 0 (its losses only); test #19 (12.88 mm spot, the measured ramp) on GRIMSEL's law and rule. The stress is the 1-D confined estimate × R (the elastic operator LEFT OUT by the user, 2026-11-06: it reads 0.26–0.65 of the exact stress at ℓ/dx ≤ 2, and a 20 mm layer under a 13 mm spot is not a 1-D slab at any mesh). **P2-24: the melt is a state** — the whole Fig. 34 trace (158 × 10 s, 1,580 s of heating, 220 s of cooling), the rock file's melt and vapour set, the beam's melt keys, the molten / lid exclusions and the face list's vapour rule (`vapor.enabled = 1`, the deck's one allowed override of the model form) |
| Hu onset, end-to-end | runs | `damage_law` failure model |

⚠ A deck rebuilt from a run's `metadata` file is **reconstruction, not replay** (SR 8): `metadata` also lists
every default that day's binary filled in.

---

## 6. Ledger A — the rock and the failure law

**Classes:**
- **MEAS-MAT**: a measured rock property.
- **MEAS-RIG**: measured or set on the experiment, without rock.
- **LIT**: a published value, possibly for another rock.
- **CAL-ONCE**: fitted once on a named test, then frozen.
- **TUNED**: set from the test's own removal, rate or shape; not allowed in a prediction.
- **CHOICE**: a value or switch we chose, with no measured source.
- **NUM**: numerical.

The decks: **K** = Kant calibration decks; **H** = Höser G1; **M** = Meier P2-2f; **Hu** = Hu; **Z** = Zhang.

| # | setting | value(s) | live where | class | source / evidence | consequence for prediction |
|---|---|---|---|---|---|---|
| 1 | `microstructure.phase0.law` | `grimsel` (H, M, K calibration deck and study legs); `constant` (legacy shipped K inputs, Hu) | all | CHOICE (form) | P2-2 | ⚠ the shipped Kant inputs are not the calibrated model |
| 2 | Grimsel laws: c_p(T), λ(T), ρ 2650, E₀ 53.3 GPa, ν 0.30, K_Ic,0 1.5 MPa√m, Table 3 ratios α/E/K_Ic to 650 °C | — | under `law = grimsel` | MEAS-MAT | `src/Numeric/Material/Grimsel.H`; `validation/grimsel_properties/grimsel_granite_properties_reference.md` | ⚠ α above 300 °C is an estimate (±30 %); **E/E₀ and K_Ic are flat above 600 °C (873 K) and α above 650 °C**, the tables' end. 22.6 % (m 13) / 34.2 % (m 8) of Meier firings sit above 873 K. E₀ sources differ: 53.3 (ours), 47.3 (Höser), 30–40 GPa (Kant) |
| 3 | `phase0.kappa/beta/E/mu/rho/Cp` (1.5 / 8e-6 / 30 GPa / … / 2750 / 790: Meier Table 7.1) | — | **inert under `grimsel`** (the law rewrites the cell fields each step) including E and μ since P2-14a: the confinement term pν/(1−ν) takes `grimsel.nu` under the law (`Removal.H` `ElasticConstants`); until then it read E/(2μ) − 1, equal by value only (P2-13 step zero, DR-1); live only under `constant` (shipped K inputs) | MEAS-MAT (Meier) | `MMWSpalling.H` ~4212, ~6978 | required keys that carry no meaning in production decks: a reader would assume they are used |
| 4 | `hc.heat.alpha` | 6.9e-7 (K, H, M); 1.3e-6 (Z); 1.14e-6 (Hu) | **inert in every deck listed**: read only by `AdvanceConstAlpha` (:1432, the no-material path) and by the `damage_law` spall-depth formula when spall is on (`Removal.H:293`; Hu runs spall off) | — | `MMWSpalling.H:141, 1432` | required but dead in these decks |
| 5 | melting and vaporisation (`phase0.T_m`, `fus_width`, `L_m`, `T_vap_lo/hi`, `L_v`) | **since P2-24 in `rock_grimsel.inputs`**: three named temperatures (P2-27) -- **T_m 1510.5 K** (the band CENTRE), **solidus 1410.5 K** = T_m − fus_width/2, **liquidus 1610.5 K** = T_m + fus_width/2 (fus_width 200 K); L_m 3.4e5 J/kg, vapour 3233–3503 K, L_v 4.8e6 J/kg; were 1e30 / 0 | every Grimsel deck: the phase table (`Table::Build`) integrates ρc_p plus the boxcars, so H(T) below 1410.5 K is unchanged (no flame gate leg passes 818 K, P2-24 z0); the liquid fraction f_L(T) feeds the beam's melt keys (rows 31–32) and the face list's molten / lid exclusions | LIT (Zhang, another rock: a Barre-like granite; no Grimsel melt data on disk) | `src/Numeric/Material/Zhang.H:59–64`; P2-24 | ⚠ the Grimsel λ(T), c_p(T) are EXTRAPOLATED through the melt (Zhang's κ cap at 1273 K is the recorded alternative). The face list carries vapour since P2-24: a LISTED cell whose H passes H_vap_hi leaves once (regime 2, `ledger_E_vapor`), a buried one is counted, not voided; `vapor.enabled` stays 0 in the flame model form. A partially molten cell (0 < f_L < 1) still carries K_Ic (declared). `spall.melt_aware` (P1's skeleton f(Λ_L) on the column path) is **deleted, P2-26** (census C3, user Q6): it never saw Λ_L > 0 on a gate deck (hottest cell 659 K vs **those gate decks' own** solidus 1448 K = their deck-local T_m 1473 − fus_width 50/2 -- not the rock file's 1410.5 K) |
| 6 | `alpha_beta.*` (quartz transition strain) | disabled | — | CHOICE | — | the 573 °C jump is in Table 3's α instead |
| 7 | `spallation.model` | `sp_weibull` (K, H, M); `damage_law` (Hu, Z) | — | CHOICE (form) | — | **two failure laws in the code.** Hu's check exercises the one we do not use for prediction |
| 8 | flaw size: `spallation.sp.a0` (fixed) / `weibull.a0_ig` (Weibull scale) | **8 µm** in both roles (`rock_grimsel.inputs`; `input_lab`'s inert 20 µm gone at P2-13) | criterion | **CAL-ONCE** | **P2-13: holds in the production form** (Kant calibration deck 850.43 K, inside 826–873 K; d1 not run). P0 (2026-09-30): "re-anchor the flaw against Kant at three meshes", 4 → 8 µm, inside Kant's binding 390–560 K rise (`ARCHIVE_DONE.md` ~L13967); held through P2-2d | ⚠ Kant 2017 states a = 20 µm (thin sections). The flaw and the drive form are one calibration (SR 1) |
| 9 | `weibull.enabled` | 1 (K calibration deck, H, M; `model_form_flame.inputs`); 0 (legacy K inputs; CAL_det control) | — | CHOICE (form) | P2-13 | ~~calibrated with one fixed flaw, applied with random flaws~~ **closed at P2-13**: the calibration deck draws Weibull flaws; the draw moves the onset +5.8 K over the fixed flaw (c3); the Weibull onset also moves ~3.7 % with mesh by construction (ROADMAP note) |
| 10 | `weibull.m` | **13** (H, M); 20 (K weibull input); 15 (default) | Weibull decks | LIT (other rock) | Westerly granite's strength modulus via Wilkinson (`Claude_markdowns/2026-09-30a.md` ~L250) | ⚠ the result is modulus-specific: P2-2f's firing T passes at m 13 and misses at m 8 (`ARCHIVE_DONE.md` ~L19529). No Grimsel value |
| 11 | `weibull.V0` with `size_volume = cell` | 1e-9 (H, 1 mm); 8e-9 (M, 2 mm) | — | CHOICE | V0 = V_cell in both ⇒ the size effect is neutralised (user choice 09-16) | flaws do not shrink with volume; the mesh dependence is accepted instead. ⚠ P2-27: the Weibull SIZE CORRECTION (V_cell/V0)^(1/m) is therefore OFF on every production deck by `check_shared.py` rule (iv); the member comment said the opposite until P2-27 (it now states what the code does). NOT deleted: switching it on is a physics decision for the mesh-check runs (§11) -- the 0.78 mm twin fired 2.09x the events with it off; on, flaws at 0.78 mm would be 8^(-1/13) ≈ 0.85x ⇒ K_I 0.92x. The one deck with V0 != V_cell: `kant/sp_kant_onset_meiercfg/input_weibull` (8e-9 on 1 mm, a legacy record) |
| 12 | `weibull.a0_gb` | 40 µm | **inert** (the rock is grain-free in every input; SR 13) | — | — | — |
| 13 | `spallation.sp.stress_relief` R | **0.531294** (`rock_grimsel.inputs`: K calibration deck, H; M, K study legs); 1.0 default (legacy K inputs) | criterion | **CAL-ONCE** (held at P2-13) | R = UCS/σ(850 K), UCS 169 ± 37 MPa (Grimsel reference), 850 K an onset target ⇒ Kant fires 843.757 K (P2-2d; `studies/p2_2d_kic/derive_onset.py`, `LEVEL.txt`) | ⚠ P2-2b: R sits ~25 % below its own physical band, "absorbing something that is not relief". With the flaw it forms one level knob |
| 14 | `spall.stress_law` / `stress_source` | `grimsel` / `local_thermoelastic` | criterion | CHOICE (form) | `Grimsel.H` (one stress formula) | 1-D confined stress from temperature; mechanics off (`el.type = disable`) ⇒ **an applied load (Höser 47 MPa) cannot be represented** |
| 15 | `spallation.sp.kic.law` | `grimsel_table` (H, M, K legs); `nasseri_table` (shipped K inputs) | criterion | MEAS-MAT | Grimsel Table 3 | ⚠ shipped K inputs use Westerly (Nasseri) |
| 16 | `confining.p` | 1 MPa (K, M); 0 (H) | criterion | M: MEAS-RIG (Meier's pads); K: CHOICE; H: — | Meier §2 (1 MPa lateral) | ⚠ Kant's binding case is p = 0 (the deck's 1 MPa moves onset ~1 K). Höser's 47 MPa run as p = 0 |
| 17 | ~~`spallation.sp.constraint_depth` d_ref~~ **RETIRED P2-16 (2026-11-03)** | — (the key aborts, naming P2-16) | deleted from both paths | — | P2-2e (2.0 mm on Gotthard, never reviewed); P2-14a carried it to the face list | retired on FORM: its onset climbs toward the minimum flux for every d (P2-14a c1: −122 K from h 500 to 2000 at d = 2 mm), against Kant 2016 Fig. 11; the buckling span (row 28) is the one minimum-flux mechanism kept. Study argvs that set d are records and no longer run (`studies/p2_16_span_rule/`, `unit/retired_keys`) |
| 18 | `spall.surface_T` | `reconstructed` | all | CHOICE (form) | P0 | firing temperature mesh-independent to 0.1–0.4 % (Kant cube) |
| 19 | `spall.removal_energy` | `flake` (all decks including `model_form_flame.inputs`); `cell` (legacy K inputs) | — | CHOICE (form) | P0b | the face list requires `flake`; calibration and production now agree (P2-13) |
| 20 | `spall.flake_c_f` | 0.4 | flake thickness h = c_f·ℓ; ℓ is the two-cell fit's e-fold length on the main branch and, **since P2-19, the measured e-fold depth ℓ_p of the line's own temperature rise on the face list's fallback branch** (`ProfileEFoldDepth`, `Shared.H`; was the tangent depth k·ΔT/q_net, undefined when the face loses heat) | CHOICE (literature-informed) | `MMWSpalling.H` ~5280: 0.3–0.5 from Wilkinson's 0.1–1 mm thicknesses and Rauenzahn at ℓ 0.8–1.6 mm; weak lever (±5 %). c_f stays: it carries the gradient physics the minimum flux rests on (`DESIGN_2026-11-05_onset_and_thickness.md` §4) | ⚠ on a constant-flux (erfc) profile the e-fold depth is 0.857× the tangent depth (P2-19 z0e), so the fallback's flakes are thinner than before. Every firing now also logs ℓ_p, the peak depth, the peak's e-fold and the energy thin limit h_E = K_Ic²(1+ν)/σ_s² (logged, never applied) |
| 21 | `spall.flake_h_max` | **1.5 mm** (H, M, K calibration deck: `model_form_flame.inputs`); 0 = off (legacy K inputs) | **NUM since P2-17: each face-list line's gather/remap window** (`FaceRemoval.H:368–380`: `L.Lg = Dmax + 2 + ⌈min(cap, h_spall_max)/(μ dx)⌉`); as a flake bound **DEAD under B** | NUM (was TUNED stopgap) | P0b rev 7, set on Meier runs; failed its own check (30.6 % of firings at the cap vs a 20 % bar). P2-17: B and the cap act on the same normal flake μ·c_f·ℓ, so the cap can bind a passing line only above T_s*(L) with C_L(T_s*) = L/1.5 mm: 1272.7 K at 12.0 mm, 1151.9 K at 13.2 mm (894 K at R1's 17.77 mm band edge) (`studies/p2_17_buckle_on/Z1.txt`, e0) | ⚠ NOT a no-op to switch off: 0 sends each line from 19 to 32 cells at 1 mm and moves the deep remap; decouple the window from the cap, then retire it (with the column path) |
| 22 | `spall.flake_length` | **DELETED P2-26** (census B3; was `twocell` on H, M and `flux` on the retired K_dref decks) | — | — | P1g-3 rejected on provenance | the two-cell length is the only one; the key is retired at `twocell` (any other value aborts naming P2-26), with `spall.flake_fraction` (`cf`), `weibull.size_volume` (`cell`) and `weibull.V_ref` (`*`) |
| 23 | `spall.flake_coherence_length` | **DELETED P2-26** (census C2; was 4 mm on K and the Meier pilot decks, 2 mm on the Hu extension decks, 0 on H, M) | — | — | — | the column path's needle cap; retired at `0` (any other value aborts). It was bypassed wherever `spall.surface_normal = 1` (deleted with it), so the two gate decks did not move through it; the out-of-gate decks that set it are records (`studies/p2_26_closed_switches/AB.txt`) |
| 24 | `spall.h_spall_max` | 15 mm | scan depth | NUM / CHOICE | — | binds only when the flake cap is off |
| 25 | `spall.sample`, `xi_s` 8, `C_h` 1, `shallow_depth` 6 mm, `damage_threshold` 0.99 | — | column path only (not read by `FaceRemoval.H`) | CHOICE (defaults) | — | inert under the face list |
| 26 | `damage.*`, `cohesive.*`, `bit.*`, `el.*` | off in flame decks; Hu runs a Drucker–Prager damage law (φ 30°, cohesion 20 MPa, …) and a static elastic solve | Hu only | Hu: MEAS-MAT (Hu's table) + CHOICE | Hu reference file | — |
| 27 | `surface_patch.robin_form` | `surface` (K calibration deck, H, M); `cell` (legacy K inputs) | the face's heat balance | CHOICE (form) | P0b | closed at P2-13 for the calibration deck; was: ROADMAP note "the Kant anchor runs the wrong BC"; ≈19 of a 21.6 % firing-time spread is this artefact |
| 28 | `spall.buckle_span` L_span (the buckling thickness condition, P2-15) | **12.0 mm in `rock_grimsel.inputs` (production, P2-17)**; 0 in the onset-only calibration deck (a declared override); 11.0 mm in the P2-15 study legs; 13.2 mm the P2-17 alternative | face list only (the column path refuses it); a fired line is withdrawn unless C_L·μ·c_f·ℓ ≤ L_span, C_L = π/√(6(1+ν)ε_th) | DERIVED-ON-KANT (Gotthard thermal; by construction on 130/140) | Rauenzahn & Tester 1989 Eqs. 15–16; `studies/p2_15_buckle/` (SPAN.txt, RESULTS.txt) | ⚠ **coupled to `c_f`** (one calibration, SR 1). ⚠ **Bethel refutes a rock-independent span under the conductivity-only transfer**: Bethel 150 and Gotthard 140 sit at the same ratio (1.090 / 1.092), so both rocks fit only a 0.02 mm window. On Grimsel at 11.0 mm, h 500 (0.41 MW/m²) never fires. Not in production. **Grimsel's span, two routes (P2-16, `studies/p2_16_span_rule/RESULTS.txt`):** (a) Gotthard carried over (conductivity-only transfer): **[9.915, 12.010) mm**; (b) the one surviving rule, **span ∝ 1/UCS** (fits both of Kant's rocks with a 10.0 % Gotthard window, L_G* 10.42 mm): **13.20 mm** (UCS 169 ± 37 MPa: 10.83–16.90 mm; with the window's edges too 10.30–17.77 mm). One span for every rock fits only a 0.15 % window (a coincidence), span ∝ 1/UCS² is excluded. ⚠ DERIVED CANDIDATE from two brackets and one fixed-exponent rule, not a law. The calibration deck at h 500 needs 12.84 mm at its earliest firing, so the two routes put Grimsel's minimum flux on opposite sides of that deck's conventional h 500. ⚠ **R1's advantage over one span rests on Bethel's onset equalling Gotthard's:** a 30 K lower Bethel onset (Kant 2016 Fig. 11/13 reads ≈ 30 K [fig]) moves Bethel's window −2.6 / −3.5 % (C_L·ℓ does not cancel: the skin shortens like √t, not like ΔT_s), which opens one span's window to 2.76 %; ≈ 21 K suffices for 2 % (P2-16 c1, a sensitivity only). **P2-17: 12.0 mm DECLARED** — one span for both Kant granites, inside one-span's window under both onset readings (Δ 0 [11.992, 12.010), Δ −30 K [11.683, 12.010) mm), chosen on parsimony; 13.2 mm (∝ 1/UCS, band 10.30–17.77) runs as the alternative. At 12.0 mm Grimsel's minimum net flux is **0.432 MW/m² (h_min 535.5 W/m²K), above the onset deck's h 500 (0.401 MW/m²)**: with B on that deck would not spall, which is why it stays onset-only (13.2 mm: 0.391 MW/m², below) (`studies/p2_17_buckle_on/Z2.txt`) |
| 29 | `beam.eps_r`, `beam.tan_delta` (P2-20) | **5.02, 0.01065** in `rock_grimsel.inputs` (all flame decks; inert at P0 = 0); 1.0 and a fixture string in the beam fixtures | the beam: α = 2π√ε′ tan δ/λ (POWER), R = ((√ε′−1)/(√ε′+1))² | MEAS-MAT (Barre granite, ANOTHER rock) | DOE DE-EE0005504 p. 44/52: n 2.24, 30 % single-pass absorption in 25 mm at 28 GHz | ⚠ flat with temperature (user, 2026-11-05); to be replaced by a Grimsel measurement |
| 30 | `beam.tan_delta_T` (P2-21) | **Hartlieb's 17 pairs** (`validation/shared/loss_table_hartlieb.inputs`, included by `rock_grimsel.inputs` and the Zhang deck): tan δ(T)/tan δ(22 °C) = 1.00 at 22 °C … 2.01 at 500 °C … 5.58 at 699 °C … 8.88 at 995 °C, T in K, linear, clamped | the beam: α(T) = α × r(T), deposited as the Beer–Lambert integral of α(T) down each column's solid path (the column optical-depth sum, `ColumnOpticalDepth`); built only when the table has ≥ 2 pairs AND the beam deposits | LIT (another granite, 2.467 GHz heating branch, cavity perturbation) | Hartlieb / Microwave Properties North, Zenodo 4835664 (`Claude_markdowns/research_notes/…/microwave_absorption.md:34–58`) | ⚠ the TRANSFER of the 2.467 GHz shape to 28 GHz is ASSUMED (ε″ ∝ f^−0.41 in the source: conduction-type loss); cooling branch (9× lower), the 675–875 °C transformation and ε′(T) not modelled; never fitted (P2-21 tests it on #19) |
| 31 | `beam.melt_factor` m (P2-24) | **4.0** in `rock_grimsel.inputs`; default 1.0 (an exact no-op) | the beam: α = α_table(T)·(1 + (m − 1)·f_L(T)) per cell in the column optical-depth sum (REQUIRED: the sum path; a deck with m ≠ 1 and a flat or absent table aborts) | LIT (Zhang eq. 12: α_λ 14 Np/m solid, 56 liquid ⇒ ×4), ADOPTED by the user 2026-11-06 | `MMWBeam.H`, `ColumnOpticalDepth` | ⚠ the DOE's "about two orders of magnitude" (DE-EE0005504 L2168–2171) is their MODEL's assumption, the recorded alternative (×100, leg M100). ⚠ m multiplies the TABLE value (124 m⁻¹ above 1268 K): ×4 ⇒ τ = 2 at 4.0 mm, ×100 ⇒ 0.16 mm (P2-24 z0) |
| 32 | `beam.R_melt` (P2-24) | **0.3** in `rock_grimsel.inputs`; default = the computed R (an exact no-op) | the beam: per column R_col = R + (R_melt − R)·f_L(T_top), at P_inc, P_trans and every cell's (1 − R) prefactor; the radiometer record's ε₁₃₇ = (1 − R_melt)·f_L + (1 − R)·(1 − f_L) | LIT (DOE L2172–2174: "0.3 after melting as observed") | P2-24 | the solid's R stays computed from ε′ (0.147) |
| 33 | `microstructure.phase0.grimsel.lambda.cap_theta` (the P2-24 addendum) | default **−1 (off)** in every deck; 1000 (°C) on the addendum's LAM leg only | `Grimsel.H` `lambda()`: θ held at the cap above it (λ 1.28 W/mK above 1000 °C instead of the law's 0.6 W/mK at 3000 K) | CHOICE (a declared sensitivity, never adopted) | `studies/p2_24_melt_vapour/RESULTS_SENS.txt` (M0R reproduces M0 on every non-timer column at the default) | gated at P2-25 (binary `c052c348`: a0 0/2639, a1 IDENTICAL at the default). The LAM leg EXCLUDED λ × 2 as the cause of the post-melt surface miss |

---

## 7. Ledger B — the flame heat source

Full per-key table with code lines: `docs/project/model_card/B_flame_source.md` (2026-10-30). Its key claims were
spot-checked against the files (7272.3958 = 7100/0.976294; the exponent fit points; `jet_De_ref`; the dead `T_ext`
and `jet_D`). The rows that decide predictions:

| setting | Höser (G1) | Meier (P2-2f) | class | note |
|---|---|---|---|---|
| level h₀ | 7272.3958 = her calorimeter 7,100 W/m²K (Ø9 mm disc mean, sensor-referenced to T_f = 1,866 K) / 0.976294 | **1449** = Martin (1977) local h at r = 2.5D, 50 mm, with air properties, measured flow, an assumed 1900 K nozzle and the model's own 821 K firing temperature; P2-2e's deck uses 1075.9 = 1449 × 0.7425, a Kant rescale (SR 17 category error) | H: MEAS-RIG; M: LIT + CHOICE | ⚠ G1 switched to Goldstein's local shape without renormalising: S = 1 gives 6,992, not 7,100 (−1.5 %, absorbed in S) |
| heat multiplier S | **0.1636** | — (1.0) | **TUNED** | S\* rule to her 33.0 cm³ at 30 s, **47 MPa**, deck unloaded ⇒ S = uptake gap × missing load |
| radial profile | Goldstein et al. 1986, local form | Martin 1977, local form, clamped at r ≥ 2.5D and H/D 2–12 | LIT | ⚠ two different correlations; both outside their range (hot, under-expanded or supersonic jets; Goldstein Re 61k–124k, H/D 6–12) |
| standoff law | (63.2/H)^1.3185 | Martin clamp, then 1/H beyond 90 mm | H: CAL-ONCE on her calorimeter at L 10–14 (79–110.6 mm), h₀ pinned at L 8; M: CHOICE | ⚠ extrapolated inside a deep hole; **double-counted past 112.35 mm** with the gas-temperature decay (`Claude_markdowns/2026-10-30.md` §1) |
| `jet_core_length` | 30 (× D_e = 112.35 mm) | 8 | H: CAL-ONCE ("K30", split unmeasured; K20 moves volume 27.9 %); M: LIT | — |
| `jet_T_ent_mode` | `fixed` | `exhaust` | CHOICE | ⚠ this switch alone moves Meier's floor gas between 1,742–1,804 K and 720–831 K. Never validated; opposite in the two decks |
| `jet_T_nozzle` | 1866 K | 1900 K | H: MEAS-RIG; M: CHOICE (unmeasured; equilibrium ≈ 2040 K, hardware bound 1450–1850 K) | — |
| `jet_mdot` | 3.4e-3 (from 38 kW / λ 1.1 set-points, ÷4) | 3.78e-3 (measured flows, ÷4) | MEAS-RIG (H derived) | quarter domain |
| `jet_cp` | 1250 | 1250 | CHOICE (no source) | ±8 % sensitivity |
| `jet_De_ref` | 3.745 mm = 7.9 × **Meier's** D_e/D ratio | 3.556 mm (isentropic Laval at an assumed 1900 K) | H: CHOICE; M: LIT (derived) | — |
| hardware | flange lid r > 20 mm (Ø40 opening) as h = 0.001 | feet 28/40 mm, 50 mm standoff, 3 pads; `foot_pad_quantile` 0.9; `foot_body_clearance` 1 | MEAS-RIG (+ NUM workarounds) | ⚠ Meier's bore = 2·`foot_r_outer` by construction |
| Meier map shape | — | taper 40→56 mm, plateau 142 W/m²K, cut 60→76 mm | **CHOICE, shape-setting** | the bore wall sits inside the taper; 142 came from a model run's own geometry |
| losses | ε 0.8, h 10 W/m²K to 293 K on every face | same | ε LIT (Höser measured 0.82); h CHOICE | no view factor; acts on in-hole faces |
| Kant decks | h 500 / 1673 K (onset: 500 is Kant 2016's minimum, 1673 Kant 2017's "≈ 1400 °C"); h 747 / 1523 K (d_ref, retired P2-16: Kant 2016 Table 4 at 130 mm) | — | CHOICE / MEAS-RIG | losses off in both |
| `surface_patch.face_heat` | `position` (default: G1) | — | CHOICE (form) | `path` = X1: heat at the surface-path distance from stagnation, with flange `barrier_expr`. **Default-off, kept as the physically preferred structure** (§11). It refuses feet, the free surface and non-surface Robin forms (the pad shield and above-plane h it also refused are deleted, P2-26) |
| dead keys | `surface_patch.T_ext`, `face`, `type`, `phase<i>.name`, `phase<i>.alpha_attenuation`, `grimsel.T_ref` **DELETED P2-18 (2026-11-05; the first three tolerated at their old default)**; `jet_D` (inert under momentum decay); `pinned_idle_cycles`, `radius` (inert under `surface`). **DELETED P2-26** (census B3, C2, C3: closed switches): `surface_patch.pad_shield_open`, `jet_above_h`, `jet_T_above`, `jet_T_ent_mode = rec_fixed` + `jet_T_rec_fixed` (strict scope: a setter aborts as unused); `spall.flake_length`, `spall.flake_fraction`, `weibull.size_volume`, `weibull.V_ref`, `spall.surface_normal`, `spall.flake_coherence_length`, `spall.melt_aware`, `melt.lambda_crit` (retired at their old defaults: the default runs, any other value aborts) | — | — | remove the rest |

**For Meier specifically:**
- **Measured by Meier:** the mass flow, the nozzle bore (7.5 mm, drawing), the 50 mm standoff, the feet (Ø80/Ø56)
  and the 3 pads.
- **From Höser or Kant:** nothing in the P2-2f deck's heat source.
- **Set by us:** the 1900 K nozzle temperature, c_p, the closure switches, and the map's taper, plateau and cut.

---

## 8. Ledger C — the MMW beam and the high-temperature physics

Full per-key tables with code lines: `docs/project/model_card/C_beam_numerics.md` (2026-10-30). Its key claims were
checked against the files:
- `alpha_attenuation` was never read (parse deleted P2-18);
- Zhang's reference values;
- Zhang is not in the tree witness list (`witness/wit.sh`);
- `RawSurfaceX` returns 0 when T₁ ≥ T_top (`MMWSpalling.H` ~8750).

**The only MMW validation is thermal: Zhang/Oglesby.**
- `validation/zhang/zhang_oglesby` runs on the global-`material` path, where spallation and vapour cannot be
  switched on.
- It compares the beam-centre surface temperature with Zhang 2023's digitised Fig. 7.
- It passed at Step 4b (peak 2871 K vs 2850 K) TUNED to Zhang's *simulation* curve (α 300, R 0.2, ε 0.8, h 51).
- **P2-20 (2026-11-06) re-ran it untuned as a recorded prediction** (`studies/p2_20_beam_deposition/RESULTS.txt`):
  28 GHz, Barre ε′ 5.02 / tan δ 0.01065 ⇒ α 14.0, R 0.147; ε 0.7/0.3, h 56; 150 s. Closure 6.5e-7, absorbed
  fraction 0.2954 (gated). Raw surface-centre T 468 / 549 / 601 K at 50 / 100 / 150 s (ω₀ 20 mm) and 565 / 670 / 734 K
  with Oglesby's 14.1 mm spot: 0.47–0.61 of the radiometer's uncorrected signal. ⚠ #28c is a glassed, re-melted disc
  (heats > 10× faster), so this shortfall measures the glass, NOT the virgin loss tangent. Zhang's own inputs (α 14,
  R 0) read 496 K at 50 s against Zhang's published 989 K: Zhang's curve is not reproducible from the paper's inputs.
  Virgin test #19 (2 kW, 50 mm): the axial PEAK reaches T_m at 160 / 318 s for ω₀ 15 / 20 mm and not by 400 s at 25 mm
  (DOE ≈ 300 ± 50 s), but 8–11 mm BELOW the surface against the DOE's 2–4 mm. **Not in the 2639 tree gate.**

| setting | Zhang deck | reference / source | class | note |
|---|---|---|---|---|
| power schedule 800–3200 W, ω₀ 0.02 m, zero standoff | as reference | summary :37–68 | MEAS-RIG | ω₀ corrected 0.20 → 0.02 against the reference; Zhang 2023 is not on disk |
| `beam.lambda` **10.71 mm (28 GHz)** since P2-20; was 2.19 mm ("137 GHz") | as the DOE report | 28 GHz heating, 137 GHz the **radiometer** (DOE DE-EE0005504) | MEAS-RIG | sets α through tan δ; z_R 0.117 m ≫ block |
| **`beam.eps_r` 5.02, `beam.tan_delta` 0.01065** (P2-20; the typed `beam.alpha_expr` 300 m⁻¹ is RETIRED, it aborts) | α = 2π√ε′ tan δ/λ = **14.0 m⁻¹** (POWER) | DOE p. 44/52: Barre, n 2.24, 30 % single-pass absorption in 25 mm at 28 GHz (= Zhang's α 14) | MEAS-MAT (Barre, another rock) | the COLD value since P2-21 (the user lifted 'flat with T', 2026-11-06); one declared pair per rock |
| **`beam.tan_delta_T`** Hartlieb's ratio table (P2-21) | α(T) = 14.0 × r(T): 28 m⁻¹ at 500 °C, 78 at 700 °C, 124 at 995 °C | none at 28 GHz (the only T-resolved granite loss on disk is 2.467 GHz) | LIT (another rock AND frequency) | the transfer to 28 GHz is the assumption P2-21 scores on #19's melt depth, onset and transmission |
| `beam.omega0` on test #19 **12.88 mm** (P2-22, used by P2-23) | ω = 12.7·√(1 + (7.9/47.3)²) mm | HE11 1/e² diameter 25.4 mm at the 32.5 mm aperture (DOE L978), standoff 7.9 mm (Appendix Table 2), z_R 47.3 mm | DERIVED | one spot, no bracket; the in-rock spread uses z_R × n (P2-22) |
| ~~`phase.T_m` 1510.5 K on the beam decks (P2-21, P2-23)~~ **left the #19 deck at P2-24** | `bc_z_melt` now reads the rock PHASE's T_m (1510.5 K, the rock file) | the Zhang material's solidus | LIT (as P2-21) | the melt is a state on the Grimsel path since P2-24 (§6 rows 5, 31, 32) |
| `beam.R` reflectivity (RETIRED as a key, P2-20) | **computed 0.1465** = ((√ε′−1)/(√ε′+1))² | the DOE's Fresnel 0.14 (Zhang applies none) | DERIVED | was typed 0.2 (TUNED) |
| `losses.epsilon` | **0.7 below T_m, 0.3 above** (P2-20) | 0.7 / 0.3 | as reference | was 0.8 (TUNED); the test no longer multiplies T by ε: it scores the raw T and the brightness (1−R)·T |
| `losses.h_conv` | **56** (P2-20) | 56 (itself "tuned" by Zhang) | as reference | was 51 (TUNED, commit d9f098b96); on #19 (P2-23): the DOE chamber's forward air purge, LIT |
| material (`material.type = zhang`) | ρ(T), κ(T), c_p(T), melt 1410–1610 K, L_m 3.4e5, vapour 3233–3503 K, L_v 4.8e6 | Zhang 2023 Table 1 via the summary | LIT | κ/ρ cap at 1273 K is "adopted"; **cannot run together with removal** |

**What the beam code does and does not do:**
- **Deposition is Beer–Lambert down each column**, energy-exact for constant α (`unit/beam` 0.00 %,
  `beam_void_closure` ≤ 1e-9).
- **α is computed from one dielectric description (P2-20), cold value plus an optional ratio table (P2-21).**
  `beam.eps_r` and `beam.tan_delta` (the rock file: Barre's 5.02 / 0.01065 until a Grimsel measurement) give the
  cold α = 2π√ε′·(tan δ/λ) and R = ((√ε′−1)/(√ε′+1))² with the deck's λ (`MMWBeam.H`; the association is locked by
  `unit/retired_keys` R6). `beam.tan_delta_T` (Hartlieb's ratios, LIT) makes α depend on T; the deposition is then
  the Beer–Lambert integral of α(T_old) down each column's solid path (`ColumnOpticalDepth`: one gathered
  optical-depth vector per step, built only when the table varies AND the beam is on, so every flat or beam-off
  deck is byte-identical). The per-column pairs telescope exactly (`beam_col_err_max` ≤ 1e-12); the closure
  keeps a residual of (d/z_R)² × the centre-to-edge contrast in α, because the Gaussian's Pi is evaluated at each
  cell's depth and a column sum cannot follow the power a spreading beam moves sideways (≈ 1 % near onset;
  `studies/p2_21_loss_table/step0`). The typed `beam.alpha_expr` and `beam.R` abort naming P2-20. Beam fixtures carry ε′ = 1
  (R = 0) and a tan δ string that reproduces their old α bit for bit (`studies/p2_20_beam_deposition/Z0.txt`).
- **P2-21 on test #19 (virgin Barre, 50 mm, at Fig. 34's MEASURED incident power 1.46 → ~2 kW; the user's
  decision, not 2 kW):** with Hartlieb's table the melt INITIATES 1.95 / 1.95 / 2.73 mm deep at ω 13 / 16.5 / 20 mm
  (DOE 2–4 mm; 0.78 mm cells), but onset comes at 52 / 84 / 122 s against ≈ 295 s, and the TOTAL transmission at onset
  is 0.85 / 0.76 / 0.72 of its beam-on value (the central column's 0.28 / 0.23 / 0.20), where Fig. 34's pick-off
  never falls below 0.79 and stands at 1.05 just before melt. The flat coefficient onsets at 384 s, 10.6 mm deep. No
  declared scaling (m = 1, 2, 3.5, 8.88) lands depth and transmission together. ⚠ The closure residual is +2–7 % of
  P_inc at onset (b1r): the column sum's limit under a Gaussian that spreads across the block. Records only; the
  table is not adopted as a 28 GHz property (`studies/p2_21_loss_table/RESULTS.txt`).
- **P2-24 on test #19 through its 26 minutes (the melt a state; Grimsel's law, Zhang's melt set, ×4 / R_melt 0.3, the
  face list's vapour rule; 1.56 mm; `studies/p2_24_melt_vapour/RESULTS.txt`).** Onset-aligned (model 43.2 s, DOE 285 s
  after beam-on), all REPORTS:
  - transmission per beam power falls to ×0.148 of its onset value in **198 s** over the whole block and **31 s** on
    the central column (DOE ×0.148 in 70 s, pick-off geometry not modelled); ×100 is SLOWER over the block (467 s):
    the top molten cell takes the whole beam, the pool spreads less and the unmelted annulus keeps transmitting;
  - absorbed fraction 0.66 at +70 s, 0.69 at ≈ 2 kW (DOE #28c ≈ 0.70); melt spot 34.9 mm at beam-off (DOE ≈ 38 mm);
  - the synthetic radiometer reads LOW after melt: 513 / 1126 °C at +20 / +100 s (DOE 1195 / 1649), the plateau 1449 /
    1623 °C at 2.0 / 3.3 kW (DOE 1774 / 2389, a +35 % rise against the model's +12 %); pre-melt 248 °C (DOE 153 mean,
    272 max); after beam-off 1009 / 699 / 327 °C at +15 / +65 / +215 s (DOE 1011 / 627 / 349);
  - ⚠ **vaporisation at 92 s and 1.57 kW**, the first buried count (92.0 s) one step before the first surface void; 324
    surface cells by beam-off, 6.6 % of the delivered energy booked as vapour (`ledger_err_solid` ≤ 2e-13 throughout).
    The reading: the ×4 absorption puts the power ~4 mm below a radiating face, and the Grimsel λ(T) extrapolated to
    3000 K (0.6 W/mK) carries little of it up. Against the DOE's 1.3 g on #28 at ≤ 2.5 kW this is a finding against the
    model's high-temperature conduction and surface balance (P2-24 P6), not a knob. ×100 never vaporises (absorbed in
    the radiating top cell); the flat-α leg (M0F) holds a buried core past H_vap_hi from 324 s, its peak AT THE PHASE
    TABLE'S 5000 K END from 396 s (a numerical ceiling: T there is clamped, H is not), and vents 1,300 cells in 7 s at
    1,250 s; the 0.78 mm leg peaks at 4,783 K inside 170 s;
  - after onset the crack rule still removes **one-cell flakes of MUSH** (0.47 ≤ f_L < 1) on the vapour pit's walls (60
    cells in M0, 312 in M0F, none before the first vapour): the declared choice "a partially molten cell carries K_Ic";
    the alternative (exclude f_L > 0) would remove none of them.
  - **ADDENDUM 2026-11-07 (`studies/p2_24_melt_vapour/RESULTS_SENS.txt`): the two recorded alternatives are EXCLUDED
    as the cause of the cold surface.** λ held at 1.28 W/mK above 1000 °C (new key `grimsel.lambda.cap_theta`,
    default off): surface unchanged (1445 / 1630 °C), vapour halved (1.74 g), but the pool deepens 3.9 → 7.0 mm and the
    cooling match is lost (812 / 427 vs DOE 627 / 349 °C); ε 0.5: +2–6 % on the surface, +21 % vapour, cooling worse.
    PARKED by the user (2026-11-07) as a stated limit: §10 row "melting / vaporisation", §11 item 6, §12 item 11.
- **P2-25 the regime map (no new physics; the #19 deck with a flat loss, 7 deposition depths × 6 absorbed fluxes, scored
  on the FIRST event; `studies/p2_25_regime_map/RESULTS.txt`).** Spall (≥ 20 removals before any melt) exactly where
  P1 put it, 1/α ≤ 1 mm at q ≥ 1 MW/m², plus 2 mm at q ≥ 4; melt-first or blocked everywhere at ≥ 4 mm (none removes,
  16 decided cells); q 0.5 blocked at every depth (the absorbed threshold is 0.48, so 0.5 is 1.04× it). Specific
  energy 1.8–24 GJ/m³ absorbed, falling with flux (1.8 at 8 MW/m², 4.4 at 2, 16 at 1); at q ≤ 2 the rate decays as
  the hole deepens and one cell (1 mm, q 1) stops at 41 s. K_Ic ×0.5 / ×2 and the spot ω 6.4 / 19.3 mm do not move the
  depth boundary; the 0.78 mm twin turns the 2 mm / q 4 cell melt-first. ⚠ 10 cells (not ≤ 6) ran out before firing:
  the 1-D stop clock is fast for deep absorption (lateral spread).
- **The deposition closure is computed on every advance path** (constant-α, material explicit and implicit,
  microstructure) and written to `thermo.dat`; the Zhang deck's test gates it (≤ 1e-3) on a flat deck, and on a
  table deck gates the per-column `beam_col_err_max` (≤ 1e-12) and reports the closure (P2-21).
- The per-phase `microstructure.phase0.alpha_attenuation` (2000 in Höser and Meier) was required by the parser and
  never read — deleted P2-18, with its deck lines.
- **Propagation is vertical (−z) only.** There is no oblique incidence, no wall absorption, no waveguiding in the
  bore, and no vapour or plasma shielding.
- **Vapour on the face list (P2-24) is a state, not a rate:** a LISTED solid cell (touching the void) whose
  enthalpy has passed its phase table's H_vap_hi leaves in the face list's one write (regime 2, the events CSV, the
  ledger term `ledger_E_vapor` = H·V booked once); an interior cell past H_vap_hi is counted (`fr_n_vapor_int`), not
  voided. The enthalpy came through deposition and conduction once and the face's losses are the face list's own
  (`unit/face_vapor`: the void at t_v, ledger_err ≤ 1e-15 across it). Melt and vapour leave no melt pool behind: no
  melt flow or evacuation is modelled.
- **The column-path vapour front (Step 11/17) is untouched and unreachable under the face list** (the dispatch
  returns first). It is tested only with synthetic properties, ignores radiation at 3300 K, uses the Gaussian formula
  under a uniform profile and, by the audit's reading, spends beam power the enthalpy field also receives with no
  ledger term — the reason the face list's rule replaced it; its deletion goes with the column path's.
- **The P0 surface reconstruction returns "flat" whenever the cell below is hotter than the surface cell.** That is
  the normal state under volumetric heating, so the firing criterion is undefined there. This is ROADMAP's
  "neither breaking rule is formulated for volumetric heating".
- **There are beam-on removal runs, but none is compared with data:**
  - `unit/beam_void_closure` (energy closure only);
  - S1 (09-15, pre-P0 criterion);
  - P2-2e B2000/B100 (1-D columns at 5 MW/m²; with d_ref = 2 mm and 1/α = 10 mm the column **never fires**;
    unreviewed);
  - the Hu-MMW extensions (09-15, "mostly vapour").

**What an unseen MMW test would need, and what we have:**

| needed | have |
|---|---|
| incident power, spot, frequency | deck keys (MEAS-RIG) |
| absorption α(T, f) of the rock | one literature α, no data on disk, no conservative T-dependence |
| reflectivity, emissivity vs T | tuned constants |
| ρ, c_p, λ to 1500–3500 K | Grimsel measured to ~500 °C, extrapolated with no clamp; Zhang's law cannot run with removal |
| melt and vapour properties, melt removal | Zhang's set (LIT, another rock) in the rock file since P2-24, with the melt's absorption ×4 and R 0.3 (LIT); no melt flow or evacuation |
| E, α, K_Ic above 650 °C, softening | flat tables; softening not implemented (DIRECTION step B) |
| a failure criterion for an interior temperature peak | none (above) |
| removal-rate data to compare with | none on disk ("MMW ablation-rate data (needed for D)", ROADMAP) |

---

## 9. Where the calibration and production models disagree

**Closed at P2-13 (2026-10-31) for the onset level.** The rock and the model form are two shared files,
`tests/MMWSpalling/validation/shared/rock_grimsel.inputs` and `model_form_flame.inputs`; every key carries its
class and source. The Kant 2017 calibration deck `validation/kant/grimsel_calibration/input` and the Höser deck
`validation/hoeser/input_lab` both include them, and `validation/shared/check_shared.py` checks that no deck
overrides a shared value undeclared and that the run resolved every one.

Results (`validation/kant/grimsel_calibration/RESULTS.txt`): per column, the first top-face non-fallback firing.

| row | measured | against |
|---|---|---|
| level, CAL (Weibull, 1 mm) | **850.43 K** | Grimsel onset band 826.15–873.15 K (Kant 2017 ∩ Höser IR): **inside ⇒ the banked constants stand** |
| fixed flaw, CAL_det | 844.63 K | banked 843.757 K: +0.88 K (the production form alone barely moves it) |
| Weibull draw | +5.80 K | offline twin +6.3 K |
| flux ×3.99 (h 500 → 2000) | +2.55 K (exponent +0.003) | data show no trend; P2-13 predicted +10 to +40 K |
| 27 / 48 MPa | −26.5 / −38.4 K | Kant measured ≈ −55 / −70 K; his model −25 / −45 K |
| mesh 2 / 1 / 0.5 mm | 853.33 / 850.43 / 853.03 K | spread 2.9 K (0.34 %) |
| all-firing median | 851.78 K | G1 851.55, P2-2f Meier 852.61 K |
| firings above 873.15 K (table clamp) | 20.7 % | — |

**The calibration deck is ONSET-ONLY (declared at P2-16, 2026-11-03).** It calibrates the onset rule A with the
buckling condition B off by definition, and nothing else. B is a pure withdrawal: a firing that passes it is identical
to the span-off firing (P2-15 DR-3, 576/576), so the onset measured there is the production onset wherever production
fires. Its h 500 W/m²K is a convention (Kant 2017 states no standoff, coefficient or flux), so the deck carries no
evidence about B or the minimum flux, and P2-15's DR-1 (the Gotthard span blocks the deck) is not a conflict between
two calibrations. With B off the onset is flat in flux (847.7 / 850.4 / 851.4 / 853.0 K at h 150 / 500 / 1000 / 2000).

**Still different, and why:**

| item | calibration deck | production |
|---|---|---|
| heat source | uniform flame, T 1673 K, h 500 W/m²K | Höser: Goldstein map; Meier: feet + jet |
| confinement | 0 (Kant's binding case) | 0 (Höser) / 1 MPa (Meier) |
| buckling condition B | off by definition: `spall.buckle_span = 0.0` after the FILE lines, declared by `# check_shared: allow spall.buckle_span` (P2-17; `check_shared.py` accepts exactly that override) | **on, 12.0 mm** (the rock file, P2-17) |
| removal path | face list | face list (Höser); column path (shipped Meier inputs, not yet re-pointed) |
| legacy Kant inputs (`sp_kant_onset*`) | `constant` law, Nasseri, R 1, `cell` | not re-pointed; their retirement needs the user's permission |

---

## 10. Validity envelope

| quantity | range with evidence | evidence | outside it |
|---|---|---|---|
| surface temperature for properties | ≤ 873 K (E, K_Ic), ≤ 923 K (α); α above 300 °C is an estimate | `Grimsel.H` Table 3 | tables are flat; a declared band, not a curve (ROADMAP) |
| firing temperature vs mesh | 0.1–0.4 % across 2, 1 and 0.5 mm (fixed flaw, Kant cube) | P0 | Weibull flaws add ~3.7 % by construction |
| flat-floor recession vs mesh | **10.4 % drift** across 2, 1 and 0.5 mm | P2-9 review DR-21 | compare shapes and rates at the same mesh only |
| thermal skin | Höser at 1 mm: ℓ ≈ 0.85–0.95 mm, about one cell (inverted from P2-7's top-cell centre 214–227 K below the face, `ARCHIVE_DONE.md` ~L22240). Meier and flat floor: 1.6–1.9 mm, under-resolved at 2 mm (ROADMAP) | as cited | rests on the surface reconstruction (clamped ℓ ≥ dz/4). ⚠ An earlier "0.4–0.8 mm (P2-7)" was the advisor's estimate, not in the repo |
| jet correlations | H/D ≤ 12 (Goldstein 6–12; Martin 2–12); subsonic cold jets | literature | both rigs are hot, under-expanded or supersonic; deep holes pass H/D 12 after ~45 s |
| Höser standoff law | 63–111 mm, open target in a Ø40 ring | P2-5 `MAP.txt` | extrapolated in deep holes; double count past 112 mm |
| onset vs flux | flat, 0.2–7.4 MW/m², five rocks (Kant, Höser, Wilkinson) | ROADMAP note | the steady-state mean surface T does rise with flux (different observable) |
| minimum flux | B off: none (the CAL deck fires at 0.12 MW/m²). With the buckling span (the one minimum-flux mechanism since P2-16): Gotthard's span 11.0 mm puts Grimsel's between 0.41 (h 500, blocked) and 0.82 MW/m² (h 1000, fires), above Kant's 0.3–0.5 for his rocks (P2-15); span ∝ 1/UCS (13.2 mm) would put it below h 500 (the deck needs 12.84 mm there; zero-compute, P2-16) | P2-15; P2-16 | span not in production; the onset where it fires is unchanged (B is a pure withdrawal) |
| onset vs flux at d > 0 | retired with d (P2-16); with B the onset where it fires is the B-off onset (flat, +2.6 K over ×4) | P2-14a; P2-15 DR-3 | — |
| applied load | none (1-D confined stress; Poisson term only) | §6 row 14 | Höser's 47 MPa data cannot be represented |
| melting / vaporisation | **conduction-only melt, tested on #19 (P2-24 + addendum):** absorbed fraction, melt-spot size and the cooling curve land; the post-melt surface reads **300–700 °C cold** with a power response of +12 % vs +35 %, and the interior reaches the boiling shelf (≈ 3 g voided where none was seen; ≈ 10× #28's loss). Not a λ or ε question (addendum). | §8 P2-24 record; `studies/p2_24_melt_vapour/RESULTS_SENS.txt` | quote NO post-melt surface temperature or vapour mass without this caveat; the spall result does not depend on it |
| volumetric (beam) heating | no criterion formulated for an interior temperature peak (ROADMAP note) | — | MMW |
| regime boundary (beam) | **P2-25 map** (#19 deck, flat loss, ω 12.88 mm, 1.56 mm, first event): spall only at 1/α ≤ 1 mm for absorbed centre flux q ≥ 1 MW/m², and at 2 mm for q ≥ 4; NOTHING removed in the 16 decided cells at 1/α ≥ 4 mm (9 blocked, 7 melt first; 8 cells undecided, C); q 0.5 (1.04× the 0.48 MW/m² absorbed threshold) fires and is blocked at every depth. At 1.56 mm the central column's e-fold `bc_ell_p` at the first firing separates the outcomes: ≤ 3.34 mm in every spall cell, ≥ 3.37 mm in every fired non-spall cell (map, K and W cuts; the 0.78 mm twin's melt-first 2 mm / q 4 cell reads 3.12 mm) | `studies/p2_25_regime_map/RESULTS.txt` | brackets: K_Ic ×0.5 / ×2 moves the onset T (428–1208 K) but NOT the depth boundary; ω 6.44 / 19.32 mm the same outcome at 1 and 4 mm; the 0.78 mm twin FLIPS the 2 mm / q 4 cell to melt-first (the boundary at 2 mm is a mesh statement) and removes 38 % less volume at 1 mm |

| time step | ≤ ~6 K overshoot per step: 16 ms at 2 mm (D2c ladder, one ladder, D2d-era deck, pre-P0 criterion); "8 ms at 1 mm" derived from the rule, not laddered | `studies/d2c_steady/RESULTS.md` ~611–634 | ⚠ Höser at her full measured heat (P2-7 HG): ≈ 19–27 K per step at 8 ms (audit arithmetic, not run). There is no code stability check |
| mesh, Höser cavity | 2/1/0.5 mm at 30 s on the 0 MPa deck: volume within 5.8 %, Ø within 2.6 mm | P2-9 c3 | **no mesh twin for G1 (120 s) or P2-7** |
| mesh, Meier | **none current**; 2 mm only (1 mm parked by the user, SR 12; older 1 mm checks failed by −28.5 % and −13.9 %) | `ARCHIVE_DONE.md` ~L10203, ~L10553 | every Meier rate is provisional at 2 mm (P2-2f) |
| decomposition | byte-identical across 1–4 ranks and box sizes 8–60 | P0, P2-2e, P2-9 | rebuilt binaries are not byte-reproducible (SR 9) |
| MMW (Zhang) numerics | implicit, dt 2 s (26× the explicit limit), 0.78 mm: **no dt or mesh study**; the latent-heat linearisation is not enthalpy-conservative across the fusion band (audit reading) | audit §2.2 | — |
| boundary conditions | Neumann-0 on every face of every deck; Höser and Meier are quarter domains (×4) | RAW_KEYS | the ledger is exact only with Neumann-0; floors must stay above the adiabatic base |

---

## 11. Open problems that limit prediction

1. **The heat-uptake gap** (3–6.5×, flame; cause unknown): `validation/hoeser/literature_uptake_gap_notes.md`,
   `studies/heat_budget/`. Beam data would show whether it is rock-side (MMW too) or gas-side (flame only).
   - **P2-17 ladder (Höser 0 MPa, face list, 30 s, target 13.0 cm³; `studies/p2_17_buckle_on/RESULTS.txt`):**
     S\*(O, B off) = 0.0994, S\*(N, B on at 12.0 mm) = 0.1014, S\*(A, 13.2 mm) = 0.1003. **B trims the gap on S
     from 10.06× to 9.87× (S\*(N)/S\*(O) = 1.019)**, and the absorbed energy at S\* by the same 1.019
     (`ledger_E_src` ×4: 24.71 → 25.18 kJ). The span choice is immaterial for this calibration (A/N 0.989). B is
     nearly inert at Höser's uptake level: it withdraws 2.9 % of firings at S 0.10 and 0.2 % at S 0.1636 (the 47 MPa deck's S\*).
     ⚠ The ratio 1.019 sits at the pre-written (a)/(b) boundary 1.02; either way the gap stays open (≈ 10×).
2. **Heat delivery to the walls of a deep hole** (the bounded check, closed 2026-10-29): free-jet rules bracket
   Höser (G1 too wide, X1 too narrow). Width predictions carry this uncertainty.
3. **The load effect** (Höser ×2.5 at 47 MPa): no mechanism in the model. It is relevant to deep MMW drilling.
4. **Shape is a strong feedback loop:** offline estimates on a fixed shape give the direction of a change, never
   its size (P2-11, P2-11b).
5. **Which flame heat rule is the more realistic basis** (advisor's judgement, 2026-10-30, adopted by the user:
   *realistic over matching*):
   - **X1 is more realistic.** It sets the heat by distance along the surface from the stagnation point, which
     follows the gas. Hrycak's wall-jet similarity and Rauenzahn's in-hole calorimeter fit both use this coordinate.
   - **G1 is less realistic.** It maps heat by position in space, so a wall near the mouth is heated as if the jet
     struck it.
   - **X1's starvation of the walls most likely comes from two open-plate pieces both rules borrowed:**
     - a standoff factor on every face, whereas Rauenzahn's in-hole peak heat was roughly independent of standoff
       and Höser's rate rises with depth;
     - a flat-plate decay along the wall, whereas Rauenzahn's CFD shows the in-hole fall-off 2–4.5× slower.
   - **Future wall-heat work starts from X1 with confined-hole versions of these two.** G1 stays the Höser deck only
     as the incumbent.
6. **Heat transport out of the melt (parked 2026-11-07).** The conduction-only melt holds the beam's energy a few mm
   below a radiating face: surface too cold, interior boiling, one cause. Excluded: λ ×2 (isotropic conduction cannot
   warm the surface without breaking the cooling curve), ε, the absorption factor (×4 and ×100 both cold), h_conv, mesh.
   Candidates that fit every constraint, in order: (i) the INSTRUMENT — the 137 GHz radiometer integrates emission over
   the melt's transparency depth and weights the centre of its view, our synthetic one reads the top cell's area mean
   (zero compute from the kept plotfiles); (ii) radiative transfer through the semi-transparent melt (λ_rad ∝ T³/κ_R,
   large only in the hot core, switches off on cooling; needs the melt's IR opacity); (iii) bubbles venting a buried
   superheat and continuous SURFACE evaporation (the physical mass-loss route; Zhang's boiling shelf is a model value).
   Pool flow is weak in granite melt (viscosity). `Claude_markdowns/2026-11-07.md`.

7. **The Weibull size correction (P2-27, recorded, not decided).** `vol_factor = (V_cell/V0)^(1/m)` makes the flaw
   scale shrink with the cell so a finer mesh does not break earlier just by drawing more often; every production deck
   sets V0 = V_cell (check_shared rule iv), which switches it OFF (§6 row 11). Decide on the mesh-check runs whether
   to fix V0 (on) or keep it per-cell (off); it moves K_I by 8^(-1/(2m)) per halving of the cell.
---

## 12. Before a blind prediction: the to-do list

The per-test and inconsistent items above, in the order they block a prediction:

1. **One removal path.** Bring vaporisation, melt, the beam and Meier's feet onto the face
   list, or decide which path is the model. ✅ P2-27: Meier's feet with the pad clearance run on the face list
   (`FrFootClearance`; `unit/face_feet`, the face gate's FEET_S / FEET_0, `input_feet_face`); what is left is the column
   path's deletion (H) and the Meier deck's rebuild on the rock file (item 5).
   - P2-23 ported the beam (the face list without a flame); P2-24 made the melt a state on it (the latent heat in the
     table the run inverts, the beam's melt absorption and reflection on f_L, no crack scan in or over a liquid) and
     moved vaporisation onto it (a listed cell past H_vap_hi leaves once, booked). Melt flow and evacuation are not
     modelled.
   - P2-19 made the face list's flake depth defined on ANY profile (the e-fold depth of the line's own rise, readable
     under an interior peak) and retired the "no heated skin" refusal (a face losing heat no longer refuses its
     line). That is the thickness prerequisite of the beam port; the port itself, and the choice between the
     face e-fold ℓ_p and the peak's e-fold ell_pk under a beam, are still to do (both are logged).
   - P2-15 added a buckling thickness condition (face list only, default off); P2-16 retired the constraint depth
     (deleted from both paths), so the buckling span is the one minimum-flux mechanism. Its span is a derived
     candidate (§6 row 28): one span for every rock fails Bethel, span ∝ 1/UCS fits both of Kant's rocks.
     It is switched on after the flake-rule check and BEFORE the uptake factor is fitted (ROADMAP).
2. ~~**One calibration deck in the production form**~~ **Done at P2-13** for the onset level
   (`validation/kant/grimsel_calibration/`; 850.43 K, constants unchanged). The deck is onset-only (§9, P2-16): the
   minimum flux is the buckling span's, scored on Kant 2016's brackets (§6 row 28), and d_ref is retired. Still to re-check, on the merged
   path (item 1). The original item:
   - Kant onset and minimum flux run with the production settings (surface heat form, Weibull flaws, the flake
     rule, removal energy, losses);
   - or production moved to the calibration form.

   Then re-check 843.757 K there.
3. **The flame heat source from rig inputs, by one rule for every flame test.** Its parts:
   - level and profile from the calorimeter where one exists, otherwise from one correlation;
   - one gas-temperature convention, which removes the double count;
   - one entrainment mode;
   - no hand-shaped taper.
4. **One global uptake factor**, calibrated once on Höser 0 MPa and frozen.
5. **Rebuild the Meier deck on today's code** with no retired keys and no Meier-specific tuned numbers, then
   predict it. Add the Höser field-tool test as a second blind case.
6. **The flake cap (1.5 mm) and m = 13** are, respectively, a stopgap and another rock's value. State them as
   such, or replace them with measured or derived values.
   - ⚠ Any change to the flake rule (`c_f`, the cap) re-opens the buckling span: the bracket fixes only
     L_span / c_f (P2-15).
7. **Dead and required-but-inert keys:**
   - `surface_patch.T_ext` (deleted P2-18, tolerated at any value);
   - ✅ P2-27: the six `phase0` constants under `grimsel` (`kappa`, `rho`, `Cp` refused under `law = grimsel`; `beta`,
     `E`, `mu` refused when `spall.stress_law = grimsel` too) and `A_D` (optional, 0.0 only) are struck; the t = 0
     property fields come from the law. STAY, with why: `hc.heat.alpha` (HeatConduction's shared required parse; read
     by the constant-alpha path), `alpha` (the homogeneous / Zhang path, D), `el.model1.{E,nu}` (Mechanics'
     `query_exactly<2>` even under `el.type = disable`, shared), `hc.bc.temp.*` (parser-optional, but the default is
     Dirichlet 0 K -- required for the answer, not by the parser);
   - `sp.a0` under Weibull;
   - `a0_gb`;
   - `jet_D` under momentum decay.
   - ✅ P2-26 deleted the closed study switches (census B3, C2, C3; twelve keys, `studies/p2_26_closed_switches/`):
     the pad shield, the above-plane pair, `rec_fixed`, the flux flake length and fraction, the flux stressed volume,
     `surface_normal`, the coherence cap and `melt_aware`. The top-level ones are retired at their old defaults
     (`unit/retired_keys` R7, R8); `side_power_interval` is now read under every `jet_closure`, and the flake-cap
     check reads the deck's `h_spall_max` (B1).

   Remove them, or stop requiring them, so no deck carries a number that does nothing (retire, don't stack).
8. **Load:** use Höser 47 MPa only for the load ratio until applied stress is modelled.
9. **Properties above 873 K:** carry a declared band.
10. **Mesh:** the 10.4 % floor drift; quote every rate and shape at one mesh. Höser G1 and Meier have no current
    mesh twin.
11. **For MMW, before any removal run:**
    - ~~check the vapour front's energy use~~ answered at P2-24 by replacement: the face list's vapour rule spends the
      enthalpy once and books it (`ledger_E_vapor`); the column-path front is unreachable under the face list;
    - ~~make α(T) deposition conservative~~ done at P2-21: α(T) is a ratio table on the cold value (a rock property,
      LIT) deposited through the column optical-depth sum (exact per column; the closure's lateral residual is
      reported); the 2.467 GHz shape's transfer to 28 GHz is scored on #19 (`studies/p2_21_loss_table/`);
    - give the failure criterion a meaning for an interior temperature peak; ⚠ P2-23 (the BEAM PORT) ran the one
      rule on a beam deck for the first time (#19 on Grimsel, 1-D stress × R, no operator): under the 28 GHz beam
      every fired line takes the fallback with NO depth inside the 15 mm window (the layer e-folds at ~19–21 mm
      on the full column) and B withdraws it -- zero removals before melt at 1.56 and 0.78 mm, the lower bound
      bk_lb ≥ 3.75 even if the layer were the window; the positive control (α 2000, 150 W, flame-like) removes
      under the same rule but only at B's MINIMUM-FLUX edge (6 events, bk 0.9994–1.0000, ℓ 2.4–2.8 mm = the
      conduction length, not the 0.5 mm absorption depth). The interior-peak criterion itself is unchanged: the
      onset still reads the 1-D stress, and which depth B should use under an interior peak stays open
      (`studies/p2_23_beam_port/RESULTS.txt`);
    - ~~restate Zhang with the reference's R, ε and h~~ done at P2-20: 28 GHz, the Barre pair (R computed), ε 0.7/0.3,
      h 56, scored raw and in brightness T as a recorded prediction (#28c is a glassed, re-melted disc; the virgin
      record is #19). It is still not in the 2639 tree gate;
    - delete the never-read `alpha_attenuation`;
    - find absorption and removal-rate data;
    - **the melt's heat transport (parked 2026-11-07, §11 item 6):** first the two instrument checks offline (depth- and
      centre-weighted brightness from the kept plotfiles), then a radiative-transfer packet scored as a bounded check if
      the gap survives; report vapour as MASS not cells; replace the boiling-shelf void with a surface evaporation flux
      when melt removal is the subject; gate or revert the addendum key `grimsel.lambda.cap_theta` (binary c052c348 is
      ungated); fix `face_gate.py`'s NaN comparator before the next face packet (P2-24 DR-1).

---

## 13. Records behind this card

- Raw tables:
  - `docs/project/model_card/RAW_KEYS.tsv` (358 keys × 9 decks, from each run's `metadata`);
  - `DECKS.tsv`;
  - `RUNNABLE.txt`.
- Audits:
  - `docs/project/model_card/B_flame_source.md` (flame source, per key, with code lines);
  - `docs/project/model_card/C_beam_numerics.md` (MMW beam, high-temperature physics, numerics, per key, with code
    lines).
- Shared files and the calibration deck (P2-13): `tests/MMWSpalling/validation/shared/`,
  `tests/MMWSpalling/validation/kant/grimsel_calibration/` (`RESULTS.txt`, `A0`–`A3.txt`, `UNITS_p2_13_vs_p2_12r.txt`,
  `legs/<leg>/argv.sh` + `BINARY.txt`).
- Daily logs: `Claude_markdowns/2026-10-30.md` §1–§4.
