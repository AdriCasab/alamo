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
   - `flake_length = flux`;
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
| Meier `input_feet`, `input_feet_d2c`, `input_drilling` (shipped) | runs | column path; no wall heat since P2-10, so not the P2-2f model. The face list cannot run them while `foot_body_clearance` is on |
| Zhang (MMW) | runs | thermal only |
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
| 5 | melting and vaporisation (`phase0.T_m`, `L_m`, `L_v`, `T_vap_*`) | 1e30 / 0 in every flame deck | — | CHOICE | — | **melting and vaporisation are switched off in every flame validation** ⇒ untested at the high end. The face list refuses `vapor.enabled` and `spall.melt_aware` |
| 6 | `alpha_beta.*` (quartz transition strain) | disabled | — | CHOICE | — | the 573 °C jump is in Table 3's α instead |
| 7 | `spallation.model` | `sp_weibull` (K, H, M); `damage_law` (Hu, Z) | — | CHOICE (form) | — | **two failure laws in the code.** Hu's check exercises the one we do not use for prediction |
| 8 | flaw size: `spallation.sp.a0` (fixed) / `weibull.a0_ig` (Weibull scale) | **8 µm** in both roles (`rock_grimsel.inputs`; `input_lab`'s inert 20 µm gone at P2-13) | criterion | **CAL-ONCE** | **P2-13: holds in the production form** (Kant calibration deck 850.43 K, inside 826–873 K; d1 not run). P0 (2026-09-30): "re-anchor the flaw against Kant at three meshes", 4 → 8 µm, inside Kant's binding 390–560 K rise (`ARCHIVE_DONE.md` ~L13967); held through P2-2d | ⚠ Kant 2017 states a = 20 µm (thin sections). The flaw and the drive form are one calibration (SR 1) |
| 9 | `weibull.enabled` | 1 (K calibration deck, H, M; `model_form_flame.inputs`); 0 (legacy K inputs; CAL_det control) | — | CHOICE (form) | P2-13 | ~~calibrated with one fixed flaw, applied with random flaws~~ **closed at P2-13**: the calibration deck draws Weibull flaws; the draw moves the onset +5.8 K over the fixed flaw (c3); the Weibull onset also moves ~3.7 % with mesh by construction (ROADMAP note) |
| 10 | `weibull.m` | **13** (H, M); 20 (K weibull input); 15 (default) | Weibull decks | LIT (other rock) | Westerly granite's strength modulus via Wilkinson (`Claude_markdowns/2026-09-30a.md` ~L250) | ⚠ the result is modulus-specific: P2-2f's firing T passes at m 13 and misses at m 8 (`ARCHIVE_DONE.md` ~L19529). No Grimsel value |
| 11 | `weibull.V0` with `size_volume = cell` | 1e-9 (H, 1 mm); 8e-9 (M, 2 mm) | — | CHOICE | V0 = V_cell in both ⇒ the size effect is neutralised (user choice 09-16) | flaws do not shrink with volume; the mesh dependence is accepted instead |
| 12 | `weibull.a0_gb` | 40 µm | **inert** (the rock is grain-free in every input; SR 13) | — | — | — |
| 13 | `spallation.sp.stress_relief` R | **0.531294** (`rock_grimsel.inputs`: K calibration deck, H; M, K study legs); 1.0 default (legacy K inputs) | criterion | **CAL-ONCE** (held at P2-13) | R = UCS/σ(850 K), UCS 169 ± 37 MPa (Grimsel reference), 850 K an onset target ⇒ Kant fires 843.757 K (P2-2d; `studies/p2_2d_kic/derive_onset.py`, `LEVEL.txt`) | ⚠ P2-2b: R sits ~25 % below its own physical band, "absorbing something that is not relief". With the flaw it forms one level knob |
| 14 | `spall.stress_law` / `stress_source` | `grimsel` / `local_thermoelastic` | criterion | CHOICE (form) | `Grimsel.H` (one stress formula) | 1-D confined stress from temperature; mechanics off (`el.type = disable`) ⇒ **an applied load (Höser 47 MPa) cannot be represented** |
| 15 | `spallation.sp.kic.law` | `grimsel_table` (H, M, K legs); `nasseri_table` (shipped K inputs) | criterion | MEAS-MAT | Grimsel Table 3 | ⚠ shipped K inputs use Westerly (Nasseri) |
| 16 | `confining.p` | 1 MPa (K, M); 0 (H) | criterion | M: MEAS-RIG (Meier's pads); K: CHOICE; H: — | Meier §2 (1 MPa lateral) | ⚠ Kant's binding case is p = 0 (the deck's 1 MPa moves onset ~1 K). Höser's 47 MPa run as p = 0 |
| 17 | ~~`spallation.sp.constraint_depth` d_ref~~ **RETIRED P2-16 (2026-11-03)** | — (the key aborts, naming P2-16) | deleted from both paths | — | P2-2e (2.0 mm on Gotthard, never reviewed); P2-14a carried it to the face list | retired on FORM: its onset climbs toward the minimum flux for every d (P2-14a c1: −122 K from h 500 to 2000 at d = 2 mm), against Kant 2016 Fig. 11; the buckling span (row 28) is the one minimum-flux mechanism kept. Study argvs that set d are records and no longer run (`studies/p2_16_span_rule/`, `unit/retired_keys`) |
| 18 | `spall.surface_T` | `reconstructed` | all | CHOICE (form) | P0 | firing temperature mesh-independent to 0.1–0.4 % (Kant cube) |
| 19 | `spall.removal_energy` | `flake` (all decks including `model_form_flame.inputs`); `cell` (legacy K inputs) | — | CHOICE (form) | P0b | the face list requires `flake`; calibration and production now agree (P2-13) |
| 20 | `spall.flake_c_f` | 0.4 | flake thickness h = c_f·ℓ | CHOICE (literature-informed) | `MMWSpalling.H` ~5280: 0.3–0.5 from Wilkinson's 0.1–1 mm thicknesses and Rauenzahn at ℓ 0.8–1.6 mm; weak lever (±5 %) | — |
| 21 | `spall.flake_h_max` | **1.5 mm** (H, M, K calibration deck: `model_form_flame.inputs`); 0 = off (legacy K inputs) | flake cap | **TUNED (stopgap)** | P0b rev 7, set on Meier runs; failed its own check (30.6 % of firings at the cap vs a 20 % bar); above Wilkinson's 0.1–1 mm | ⚠ calibration runs without it, production with it |
| 22 | `spall.flake_length` | `twocell` (H, M); `flux` (K_dref) | — | CHOICE (form) | P1g-3 | ⚠ the face list refuses `flux` |
| 23 | `spall.flake_coherence_length` | 4 mm (K); 0 (H, M) | column path | CHOICE | — | ⚠ the face list refuses > 0 |
| 24 | `spall.h_spall_max` | 15 mm | scan depth | NUM / CHOICE | — | binds only when the flake cap is off |
| 25 | `spall.sample`, `xi_s` 8, `C_h` 1, `shallow_depth` 6 mm, `damage_threshold` 0.99 | — | column path only (not read by `FaceRemoval.H`) | CHOICE (defaults) | — | inert under the face list |
| 26 | `damage.*`, `cohesive.*`, `bit.*`, `el.*` | off in flame decks; Hu runs a Drucker–Prager damage law (φ 30°, cohesion 20 MPa, …) and a static elastic solve | Hu only | Hu: MEAS-MAT (Hu's table) + CHOICE | Hu reference file | — |
| 27 | `surface_patch.robin_form` | `surface` (K calibration deck, H, M); `cell` (legacy K inputs) | the face's heat balance | CHOICE (form) | P0b | closed at P2-13 for the calibration deck; was: ROADMAP note "the Kant anchor runs the wrong BC"; ≈19 of a 21.6 % firing-time spread is this artefact |
| 28 | `spall.buckle_span` L_span (the buckling thickness condition, P2-15) | **0 = off** in every shared input; 11.0 mm in the P2-15 study legs | face list only (the column path refuses it); a fired line is withdrawn unless C_L·μ·c_f·ℓ ≤ L_span, C_L = π/√(6(1+ν)ε_th) | DERIVED-ON-KANT (Gotthard thermal; by construction on 130/140) | Rauenzahn & Tester 1989 Eqs. 15–16; `studies/p2_15_buckle/` (SPAN.txt, RESULTS.txt) | ⚠ **coupled to `c_f`** (one calibration, SR 1). ⚠ **Bethel refutes a rock-independent span under the conductivity-only transfer**: Bethel 150 and Gotthard 140 sit at the same ratio (1.090 / 1.092), so both rocks fit only a 0.02 mm window. On Grimsel at 11.0 mm, h 500 (0.41 MW/m²) never fires. Not in production. **Grimsel's span, two routes (P2-16, `studies/p2_16_span_rule/RESULTS.txt`):** (a) Gotthard carried over (conductivity-only transfer): **[9.915, 12.010) mm**; (b) the one surviving rule, **span ∝ 1/UCS** (fits both of Kant's rocks with a 10.0 % Gotthard window, L_G* 10.42 mm): **13.20 mm** (UCS 169 ± 37 MPa: 10.83–16.90 mm; with the window's edges too 10.30–17.77 mm). One span for every rock fits only a 0.15 % window (a coincidence), span ∝ 1/UCS² is excluded. ⚠ DERIVED CANDIDATE from two brackets and one fixed-exponent rule, not a law. The calibration deck at h 500 needs 12.84 mm at its earliest firing, so the two routes put Grimsel's minimum flux on opposite sides of that deck's conventional h 500. ⚠ **R1's advantage over one span rests on Bethel's onset equalling Gotthard's:** a 30 K lower Bethel onset (Kant 2016 Fig. 11/13 reads ≈ 30 K [fig]) moves Bethel's window −2.6 / −3.5 % (C_L·ℓ does not cancel: the skin shortens like √t, not like ΔT_s), which opens one span's window to 2.76 %; ≈ 21 K suffices for 2 % (P2-16 c1, a sensitivity only) |

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
| `surface_patch.face_heat` | `position` (default: G1) | — | CHOICE (form) | `path` = X1: heat at the surface-path distance from stagnation, with flange `barrier_expr`. **Default-off, kept as the physically preferred structure** (§11). It refuses feet, the pad shield, the free surface, above-plane h and non-surface Robin forms |
| dead keys | `surface_patch.T_ext` (no reader); `jet_D` (inert under momentum decay); `pinned_idle_cycles`, `pad_shield_open`, `radius` (inert under `surface`) | — | — | remove |

**For Meier specifically:**
- **Measured by Meier:** the mass flow, the nozzle bore (7.5 mm, drawing), the 50 mm standoff, the feet (Ø80/Ø56)
  and the 3 pads.
- **From Höser or Kant:** nothing in the P2-2f deck's heat source.
- **Set by us:** the 1900 K nozzle temperature, c_p, the closure switches, and the map's taper, plateau and cut.

---

## 8. Ledger C — the MMW beam and the high-temperature physics

Full per-key tables with code lines: `docs/project/model_card/C_beam_numerics.md` (2026-10-30). Its key claims were
checked against the files:
- `alpha_attenuation` is never read;
- Zhang's reference values;
- Zhang is not in the tree witness list (`witness/wit.sh`);
- `RawSurfaceX` returns 0 when T₁ ≥ T_top (`MMWSpalling.H` ~8750).

**The only MMW validation is thermal: Zhang/Oglesby.**
- `validation/zhang/zhang_oglesby` runs on the global-`material` path, where spallation and vapour cannot be
  switched on.
- It compares the beam-centre surface temperature with Zhang 2023's digitised Fig. 7.
- It passed at Step 4b: peak 2871 K vs 2850 K. That was tuned to Zhang's *simulation* curve (`ARCHIVE_DONE.md`
  ~L246).
- Last run 2026-05-30, before every later source change. **Not in the 2639 tree gate.**

| setting | Zhang deck | reference / source | class | note |
|---|---|---|---|---|
| power schedule 800–3200 W, ω₀ 0.02 m, zero standoff | as reference | summary :37–68 | MEAS-RIG | ω₀ corrected 0.20 → 0.02 against the reference; Zhang 2023 is not on disk |
| `beam.lambda` 2.19 mm ("137 GHz") | — | the reference names 137 GHz as the **radiometer** frequency; the source frequency is not stated | CHOICE | inert here (Rayleigh length ≫ block) |
| `beam.alpha_expr` 300 m⁻¹ (constant) | — | `main.tex` citing Woskov 2009 (not on disk); the reference gives no α | LIT (secondary) | the one MMW absorption number; no rock-specific or T-dependent data |
| `beam.R` reflectivity | **0.2** | **0** | **TUNED** | degenerate with P and α in the absorbed power |
| `losses.epsilon` | **0.8** | **0.7 below T_m, 0.3 above** | **TUNED** | the test also scores ε·T with this ε, so the raw peak is ≈ 3590 K, above the vapour window, with no removal |
| `losses.h_conv` | **51** | **56** (itself "tuned" by Zhang) | **TUNED** | changed 56 → 51 in commit d9f098b96 |
| material (`material.type = zhang`) | ρ(T), κ(T), c_p(T), melt 1410–1610 K, L_m 3.4e5, vapour 3233–3503 K, L_v 4.8e6 | Zhang 2023 Table 1 via the summary | LIT | κ/ρ cap at 1273 K is "adopted"; **cannot run together with removal** |

**What the beam code does and does not do:**
- **Deposition is Beer–Lambert down each column**, energy-exact for constant α (`unit/beam` 0.00 %,
  `beam_void_closure` ≤ 1e-9).
- **α(T) is not energy-conserving.** `BeamSource` uses the cell's own α without the cumulative optical depth; the
  code says so (`MMWSpalling.H` ~4423–4427). ROADMAP's and DIRECTION's pointer to `:6166` is stale.
- **Absorption is one global `beam.alpha_expr`.** The per-phase `microstructure.phase0.alpha_attenuation` (2000 in
  Höser and Meier) is required by the parser and never read.
- **Propagation is vertical (−z) only.** There is no oblique incidence, no wall absorption, no waveguiding in the
  bore, and no vapour or plasma shielding.
- **The vapour front (Step 11/17) is tested only with synthetic properties.**
  - It ignores radiation losses at 3300 K (≈ 5 MW/m²).
  - It uses the Gaussian formula even under a uniform profile.
  - By the audit's reading, it may spend beam power that the enthalpy field also receives, with no ledger term.
    **Unverified**; to check before any MMW removal run.
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
| melt and vapour properties, melt removal | scalars defaulted off in every deck; no melt evacuation |
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
| buckling condition B | off by definition (onset-only deck) | off (span 0) until the step that switches it on |
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
| melting / vaporisation | **untested** in any validation with removal | §6 row 5 | — |
| volumetric (beam) heating | no criterion formulated for an interior temperature peak (ROADMAP note) | — | MMW |

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

---

## 12. Before a blind prediction: the to-do list

The per-test and inconsistent items above, in the order they block a prediction:

1. **One removal path.** Bring vaporisation, melt, the beam and Meier's feet onto the face
   list, or decide which path is the model. Today no MMW case and no Meier-with-clearance case can run on the
   current removal rule.
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
   - `surface_patch.T_ext`;
   - `hc.heat.alpha` and the `phase0` constants under `grimsel` (E and μ no longer set ν: since P2-14a the confinement term takes `grimsel.nu`);
   - `sp.a0` under Weibull;
   - `a0_gb`;
   - `jet_D` under momentum decay.

   Remove them, or stop requiring them, so no deck carries a number that does nothing (retire, don't stack).
8. **Load:** use Höser 47 MPa only for the load ratio until applied stress is modelled.
9. **Properties above 873 K:** carry a declared band.
10. **Mesh:** the 10.4 % floor drift; quote every rate and shape at one mesh. Höser G1 and Meier have no current
    mesh twin.
11. **For MMW, before any removal run:**
    - check the vapour front's energy use (possible double use of beam power);
    - make α(T) deposition conservative (cumulative optical depth);
    - give the failure criterion a meaning for an interior temperature peak;
    - put Zhang back in the regression set and restate it with the reference's R, ε and h (or state why not);
    - delete the never-read `alpha_attenuation`;
    - find absorption and removal-rate data.

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
