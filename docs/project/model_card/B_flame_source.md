# Model card, part B: the flame heat source and surface boundary keys

Scope: every `surface_patch.*` and `losses.*` key in the flame decks, and every literal inside the Höser and Meier
`surface_patch.h_expr` strings. Compiled 2026-10-30 against **today's working tree** (`src/Integrator/MMWSpalling.H`,
9,148 lines, modified 2026-10-30 18:00, which includes P2-12's removal of the path mode). No repository file was edited
and no binary was run.

**Decks** (`docs/project/model_card/DECKS.tsv`): KO = KANT_onset, KD = KANT_dref, H47 = HOSER_47 (G1), H0 = HOSER_0 (HG),
H0S = HOSER_0_S (HG with h × 0.10), M = MEIER (P2-2f m13), Md = MEIER_dref (P2-2e SD).

**Where the deck values come from.** They were read from the base input file plus `argv.sh`, and cross-checked against
`RAW_KEYS.tsv`.
- ⚠ `RAW_KEYS.tsv` is built from `metadata`, which is `pp.dumpTable` (`src/IO/WriteMetaData.cpp:58`). It therefore
  also contains every key that the binary of the day added through `query_default`. Examples: `side_face_closure`,
  `side_face_area_exact`, `side_face_factor`, `side_face_h_max`, `side_face_r_min`, `jet_body*`, `jet_wall_mdot`, and
  `spall.lateral_* = 0` on KO.
- Those keys are not deck settings, and a deck rebuilt from `metadata` is reconstruction, not replay (SR 8).
- In the `overridden` column, `surface_patch.enabled ... :0` means the base input file sets 0 and argv sets 1.

**Parse sites.** Line numbers are `MMWSpalling.H` unless another file is given. `SurfacePatch::Parse` runs at
7672–8106 and is called through the strict `pp.queryclass<SurfacePatch>("surface_patch", …)` at 178. The losses are
parsed at 153–161.

---

## A. `surface_patch.*` keys

| key | KO | KD | H47 | H0 | H0S | M | Md | where read / consumed | LIVE? | CLASS | provenance | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| enabled | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 7675 (if 0, Parse returns at once) | LIVE all | FORM | `input_lab:164` and `input_drilling:167` set 0; argv sets 1 | — |
| mode | convective_flame | = | = | = | = | = | = | 7694; flux in `BuildFlameColumns` 2056 and `SurfaceCellFlux` 3628 | LIVE all | FORM | Step 15a | — |
| type | prescribed_temperature (default) | = | = | = | = | = | = | 7690; checked only under `prescribed_T` | INERT all | FORM (legacy) | 7688 comment | — |
| face | zhi | = | = | = | = | = | = | 7764 (only `zhi` allowed) | INERT | FORM | — | — |
| face_list | – (0) | – | **1** | – | – | – | – | 7681; guards 949–1060; `FaceListHeat` 3559, `FaceH` 3382 | LIVE H47 only | FORM | P2-8/P2-9 (`ARCHIVE_DONE.md` L22369, L23207) | Under the face list `foot_body_clearance` aborts (993–994), so M/Md cannot simply turn it on |
| robin_form | cell | cell | surface | surface | surface | surface (argv `pinned`, then `surface`; last wins) | surface | 7712; `SurfaceCellFlux` 3674–3687 (surface), ~3735 (cell) | LIVE all | FORM | P0b (L14048) | ⚠ ROADMAP note *"THE KANT ANCHOR RUNS THE WRONG BC"*: the rock-law anchor runs `cell` while every flame deck runs `surface` |
| pinned_idle_cycles | – | – | 1e9 | 1e9 | 1e9 | 1e9 | 1e9 | 7743, read only to satisfy the strict parser; consumer `PinRule` 3083 runs under `pinned` only; "INERT" message at 758 | INERT all | NUM | C1/D2i (memory: idle clock off) | — |
| pad_shield_open | – (−1) | = | = | = | = | = | = | 7754–7755; 2140, 2178, 3396 | INERT (sentinel) | FORM | D2aa | Never set on a scored deck |
| x0, y0 | 0.012 | 0.012 | 0 | 0 | 0 | 0 | 0 | 7768–7769; r at 2091, 2453, 3384 | KO/KD INERT (uniform h); H/M LIVE | FORM (geometry) | Jet axis on the domain corner: a quarter domain, ×4 (SR 3) | — |
| radius | 0.05 | 0.05 | 0.08 | 0.08 | 0.08 | 0.2 | 0.2 | 7770; hard mask at 2096, 2453, 1950, 3122, 3450, 3860 (`FaceH` ignores it) | INERT all: it covers the domain (Kant 24 mm cube / 4 mm column; Höser corner 70.7 mm; Meier corner 169.7 mm) | NUM | `run_p2_5.py:110` comment; D2t | — |
| T_expr | – | – | – | – | – | – | – | 7784 | INERT | — | — | — |
| T_ext | 293.15 | 293.15 | – (0) | – | – | – | – | 7791; **no consumer anywhere in src** | INERT (dead key) | DEFAULT-UNJUSTIFIED | grep of `src/` | Can be deleted from the Kant inputs |
| h_conv | **500** (input 150; harness sets 500) | **747** (input 150; argv 747) | – | – | – | – | – | 8092; `BuildFlameColumns` 2063–2065 (uniform `flame_h`) | LIVE (KO, KD) | KO: DEFAULT-UNJUSTIFIED; KD: MEAS-RIG | KO: `input_model:37`, `test:16-20`. Kant 2017 gives no h_fl (`kant-2017-central-aare-validation-reference.md:69-71`, "the free knob"); 500 is **Kant 2016's** h_fl,min ≈ 500 ± 100 (`kant_von_rohr_…plan.txt:121`). KD: Kant & von Rohr 2016 Table 4, SOD 130 mm (`plan.txt:73`) | KD's value is indirect: a sensor plate replacing the rock, ±17 % (`plan.txt:76-78`). Onset T should be h-independent in theory, but the P1 gradient criterion makes it weakly dependent, so the calibration point matters. Combining it with `h_expr` aborts (8009) |
| T_flame | 1673.15 | 1523.15 | – | – | – | – | – | 8093; 2066; `SurfaceGasTempCap` 8871 | LIVE (KO, KD) | KO: MEAS-RIG (approximate); KD: MEAS-RIG | KO: Kant 2017 "≈ 1400 °C" (`kant-2017-…md:68`). KD: Table 4, SOD 130, 1250 °C (`plan.txt:73,84`) | The flame-temperature convention in Table 4 is not stated beyond "T_flame" |
| h_expr | – | – | G1 string | HG string | HG × 0.1 | G4 (1449) | G4 (1075.9) | 7799 / 8013–8016; `h_f` at 2113 (enthalpy branch, **s > 0 only**), 2151, `FaceH` 3390 | LIVE (H, M) | see Tables C, D | — | Under `enthalpy` it is the floor h (ROADMAP note) |
| jet_closure | none (default) | none | enthalpy | = | = | = | = | 7806 | LIVE | FORM | D2a2 (L5637) | — |
| jet_mdot | – | – | 3.4e-3 | = | = | 3.78264e-3 | = | 7819; 2500 (m = ṁ/φ), 2504 (mcp), 2664 (P_nozzle), 2678 | LIVE (enthalpy budget, radial T_gas decay, cap) | MEAS-RIG (H: derived; M: measured) | H: `run_p2_5.py:96-97`: CH₄ 0.684 g/s = 38 kW / 55.5 MJ/kg (HHV), air ×17.2×1.1 ⇒ 13.63 g/s, /4. M: (52 + 2.47 kg/h)/4 (`d2a_jet_face/jet.py:23`; Meier Fig. 8.7, `meier_thesis_notes.md:42`) | Quarter-domain share confirmed for both. ⚠ H is computed from set-points, not a measured flow; whether her 38 kW is HHV-based is unverified (LHV gives +11 %) |
| jet_T_nozzle | – | – | 1866 | 1866 | 1866 | 1900 | 1900 | 7820; 2498 (T_stag), 2665, 8872 | LIVE | H: MEAS-RIG; M: DEFAULT-UNJUSTIFIED | H: Höser 2018 IJRMMS (`hoeser_2018_uniaxial_load_notes.md:9`); Höser & von Rohr 2018 Table 2, 1,593 °C (`hoeser_2016_thesis_notes.md:154`). M: "unmeasured; 1900 adiabatic is an upper bound" (`validation/meier/README.md:108,139-141`) | ⚠ M: equilibrium T_ad at λ 1.2 is ≈ 2040 K, so 1900 is not adiabatic (`Claude_markdowns/2026-09-22a.md:103-105`); the hardware bound is 1450–1850 K (ibid. :114-117). Höser's 1,600 °C was used only as a plausibility check (`thesis_notes.md:24`) |
| jet_cp | – | – | 1250 | = | = | = | = | 7821; 2504, 2631, 2664 | LIVE | DEFAULT-UNJUSTIFIED | "T-averaged cp of the products (planner's value)" (`d2a2_walljet/walljet.py:32`); ±8 % sensitivity (L6133) | — |
| jet_D | – | – | 0.0079 | = | = | 0.0075 | = | 7822; **only consumer is 2493, under `jet_decay_diameter = nozzle`** | **INERT on every flame deck** (all run `momentum`) | MEAS-RIG | H: D1 nozzle (`thesis_notes.md:19`). M: VT5 drawing sheet 14/17 (reference errata 1; the thesis text says 7.1 mm) | D enters the result only through the literals inside `h_expr` |
| jet_T_ent | – | – | 293.15 | = | = | = | = | 7823; 2472 (march floor and reference), 2492 (D_e ∝ √(T_mix/T_ent)), above-plane `flame_Tg` 2146 | LIVE | DEFAULT-UNJUSTIFIED | Assumed 20 °C lab ambient; no rig value was opened | H: in `fixed` mode this is the entrained-gas temperature |
| jet_core_length | – | – | 30 | 30 | 30 | 8 | 8 | 7827–7828 (explicit value required under `momentum`, 7891); 2496–2497: φ = min(1, core·D_dec/s_c) | H: LIVE once s_c > 112.35 mm (G1 from ~90 s); M: LIVE | H: CAL-ONCE; M: LIT | H: P2-5 rev 2 "K30", φ = 1 out to 30 × 3.745 = 112.35 mm so the coefficient carries all of the calorimeter's standoff decay (`p2_5_hoeser/MAP.txt`; L21030-21045). M: cold supersonic core L/D_e ≈ 4.2 + 1.1 M², M ≈ 1.86 (L7441, L7528) | H: the gas/coefficient split is **unmeasured**. The K20 bracket moves volume by 27.9 % (L20990). ⚠ Today's log §1: past 112.35 mm the standoff decay is double-counted. M: the core is in D_e units, and under `exhaust` D_e ≈ 8 mm, so the core is ~64 mm |
| jet_stagnation | – | – | decay | = | = | = | = | 7830; 2495 | LIVE | FORM | D2a2/D2b recommendation (L6177) | Today's log names `nozzle` as the candidate fix on Höser |
| jet_entrained_mass | – | – | 1 | = | = | = | = | 7836; 2499–2500 | LIVE | FORM | D2a2/D2b | — |
| jet_T_ent_mode | – | – | fixed | fixed | fixed | exhaust | exhaust | 7844; T_mix at 2484–2487; recirculation 2668–2681 | LIVE | FORM | H: "by design: a ventilated box" (`input_lab:62`). M: blind-hole recirculation (D2b) | ⚠⚠ M: this switch alone sets the floor gas temperature, T_stag 1,742–1,804 K under `exhaust` against 720–831 K under `fixed` at 85–107 mm (ROADMAP note *"THE MEIER FLOOR GAS NEVER CAME FROM THE JET CLOSURE'S DECAY"*). Never validated |
| jet_decay_diameter | – | – | momentum | = | = | = | = | 7876; 2491–2493 | LIVE | FORM | D2c (L7437) | — |
| jet_De_ref | – | – | 0.003745 | = | = | 0.0035556 | = | 7889; 2492 | H: LIVE only through the product core × D_e = 112.35 mm; M: LIVE | H: DEFAULT-UNJUSTIFIED; M: LIT (derived) | H: 7.9 mm × **Meier's** D_e/D ratio 0.47408 (`run_p2_5.py:98`). M: `d2c_steady/nozzle.py` (isentropic Laval nozzle at the measured ṁ, **assumed** T₀ 1900 K, γ 1.3, R 290, 6 mm throat; `meier-2017-…md:67-73`) | H: a ratio for a Laval nozzle borrowed for her convergent, underexpanded nozzle |
| jet_free_surface | – | – | – (0) | 0 | 0 | 0 | 0 | 7900; 2432–2465, 2617–2633 | INERT | FORM | D2c (scored = 1 there; 0 here) | — |
| jet_bin_update | – | – | exponential | = | = | = | = | 7917; 2569–2584 | LIVE (numerics) | NUM | D2c | — |
| jet_negative_flux | – | – | count | = | = | = | = | 7923; 2574, 2591–2596 | LIVE | FORM | D2c (rock hotter than the gas heats the gas) | — |
| jet_above_h / jet_T_above | – | – | 0 / −1 (defaults) | = | = | = | = | 7932 / 7937–7938; 2112–2113, 2146, 3389 | INERT | FORM | D2v | Because these are off, the Meier `h_expr` `s ≤ 0` branch is inert (Table D) |
| jet_dr | – | – | 0 (default = max(dx, dy)) | = | = | = | = | 7950; 2424 | LIVE (bin width = mesh) | NUM | — | — |
| jet_profile_interval | – | – | 1.0 | 1.0 | 1.0 | 10 | 2.0 | 7953; 2550 | INERT (CSV log only) | DIAG | — | — |
| side_power_interval | – | – | 10 | = | = | = | = | 7956; 3808, 3944–3950 | INERT (CSV) | DIAG | — | — |
| nozzle_descent | prescribed (default) | = | prescribed | = | = | **feet** | feet | 7987; 2078–2080 | LIVE | H: MEAS-RIG (fixed burner); M: FORM | H: "burner FIXED" (`thesis_notes.md:84`). M: the burner rests on its feet (Meier p. 212) | M: the support rule is a model (memory `meier-burner-support-rule`) |
| nozzle_z0 | – | – | 0.2132 | = | = | – | – | 8068; 2080 | LIVE | MEAS-RIG | 0.15 m block = her Ø100 × 150 sample + 63.2 mm (L = 8) (`input_lab` header; P2-7 L22100) | — |
| nozzle_feed | – | – | 0.0 | = | = | – | – | 8069; 2080 | LIVE (standoff only opens) | MEAS-RIG | Fixed burner | — |
| nozzle_collision_radius | – | – | 0.004 | = | = | 0.028 | = | 8079; 2076–2077, abort guard 2097–2103 | INERT unless the nozzle reaches the rock (guard) | NUM | SR 6; M: = `foot_r_inner` (the feet default) | — |
| foot_r_inner / foot_r_outer | – | – | – | – | – | 0.028 / 0.040 | = | 8029–8030; `FeetDescent` 2217–2320 | LIVE | MEAS-RIG | VT5 drawing sheet 15/17: Ø80 tube, Ø56 bore (`meier_thesis_notes.md:36-41`) | ⚠ The bore is 2·`foot_r_outer` **by construction** (issues file G5; P2-2f bore 79.3 mm) |
| foot_standoff | – | – | – | – | – | 0.050 | = | 8031; 2307 | LIVE | MEAS-RIG | Drawing sheet 1/17; thesis p. 212 ("L = 7") | 6.67 D at D = 7.5 mm |
| foot_rule | – | – | – | – | – | pads | = | 8039 | LIVE | FORM | D2b | — |
| foot_npads | – | – | – | – | – | 3 | = | 8045; 2251 | LIVE | MEAS-RIG | Three exhaust slots ⇒ three arc pads (`meier_thesis_notes.md:40`) | — |
| foot_pad_quantile | – | – | – | – | – | 0.9 | = | 8048; 2278, 2297 | LIVE | NUM | D2f: a guard against one never-fired 2 mm column freezing the burner; q 0.8/0.9/1.0 changes ROP by ±2.9 % (L9602; memory) | Meier's physical rule is the maximum (q = 1) |
| nozzle_feed_max | – | – | – | – | – | inf (default) | = | 8051; 2316 | INERT | — | — | — |
| foot_body_clearance | – | – | – | – | – | 1 | 1 | 8058; arms `foot_clear_k` in `FeetDescent`; `Removal.H:383-501` (regime 4) | LIVE | FORM | D2c (L7714). D2e found it to be the 1 mm/2 mm mesh seed (L9304; memory: "no physical basis") | P2-2f books its enthalpy (`foot_mech_E`); pad clearing was 0.026 L of 0.59 L. Under `face_list` it aborts |

### A2. Retired keys still in the argv of accepted decks (P2-10 deleted their code)

All of these abort today: the retired list at 94–120 aborts on any key present, and the strict parser aborts on the
rest. They did shape the H0, H0S, M and Md results when those decks were run.

| key | H0 / H0S | M | Md | class when it was live | provenance |
|---|---|---|---|---|---|
| side_face_flux | 1 | 1 | 1 | FORM | D2k |
| side_face_h | 169.3 | 169.3 | 169.3 | DEFAULT-UNJUSTIFIED (typed) | P2-1 found the derived wall law 1.43× this (L18838) |
| side_face_h_law | derived ‹hash› | = | = | LIT-derived | P2-1 |
| side_face_h_r_min | 0.036 | 0.040 | 0.036 | DEFAULT-UNJUSTIFIED | Meier skirt radius, lowered from 40 so a 39.77 mm bore counts as wall (issues G1). Höser inherited it by copying |
| side_face_T / _area_law / _step_height / _staircase / _h_datum | live / wetted / 1 / 1 / bore | live / wetted / 1 / 1 / – | live / wetted / 1 / 1 / bore | FORM | D2q–D2z, P2-3 |
| side_face_factor_above | – | 0.0 | 0.0 | FORM (a switch; issues G7) | D2v/D2w |
| step_diag_radius | – | 0.050 | 0.050 | DIAG | — |
| jet_wall_march | – | – | 0 | FORM | P2-3 |
| side_face_h_scale | – | – | 1.0 | — | **Not on the retired list and no reader** ⇒ strict "specified but not used" abort |
| spall.lateral_drive / _removal / _subcell / _thermal | 1 / 1 / 1 / 1 | = | = | FORM | P1f/P1g (outside this scope, but on the abort list) |

---

## B. `losses.*` keys

Parsed at 153–161; consumed in `SurfaceCellFlux` 3678–3681 (`surface` form, at the reconstructed T_s) and ~3735
(`cell` form), the march `q_col` 2523, and `FrFaceFlux` at `FaceRemoval.H:149,154,487`. Guards are at 700, 720 and 746.

| key | KO | KD | H47 / H0 / H0S | M / Md | LIVE? | CLASS | provenance | notes |
|---|---|---|---|---|---|---|---|---|
| epsilon | 0.0 | 0.0 | 0.8 | 0.8 | Kant INERT; flame LIVE | Kant: FORM (losses off); flame: LIT | Consistent with Höser's measured Grimsel ε = 0.82 ± 10 % (`thesis_notes.md:30`). The value was inherited from `input_drilling` (2026-06) | ⚠ It radiates to 293 K from every heated face, deep in-hole faces included: no view factor and no radiation from the hot gas |
| epsilon_high | 0.0 | 0.0 | 0.8 | 0.8 | INERT (single ε) | NUM | — | — |
| T_eps_break | 1e30 | = | = | = | INERT | NUM | — | — |
| sigma_SB | 0.0 | 0.0 | 5.670374e-8 | = | Kant: switch off; flame LIVE | LIT (physical constant) | — | — |
| h_conv | 0.0 | 0.0 | 10.0 | 10.0 | Kant INERT; flame LIVE | DEFAULT-UNJUSTIFIED | "h_ext ≈ 10 W/m²K natural convection in lab air" (`meier-2017-…md:188`) | ⚠ It acts on every face, including faces under the jet, as a second heat path to 293 K alongside the flame Robin term |
| T_amb | 293.15 | 293.15 | 293.15 | 293.15 | Kant INERT; flame LIVE | DEFAULT-UNJUSTIFIED | Assumed 20 °C | (Zhang 297.5, Hu 280.15) |

---

## C. Literals in the Höser `h_expr`

Strings (argv):
- H47 (G1): `if((r>0.02)and(s<0.0633),0.001,1189.763953*533*(533+0.3575*44*pow(r/0.0079,1.285))/pow(533+44*pow(r/0.0079,1.285),2)*pow(0.0632/max(0.0632,s),1.3185305455))`
- H0 (HG): `if((r>0.02)and(s<0.0637),0.001,7272.3958*533/(533+44*pow(r/0.0079,1.285))*pow(0.0632/max(0.0632,s),1.3185305455))`
- H0S: the same as HG with `727.239580`.

| literal | value(s) | meaning | LIVE? | CLASS | provenance | notes |
|---|---|---|---|---|---|---|
| `r>0.02` | 0.02 m | Flange lid: Ø40 mm opening | LIVE (shades the unrecessed top beyond r = 20 mm) | MEAS-RIG | 20 mm flange + 15 mm disc, Ø40 opening (`hoeser_2018_uniaxial_load_notes.md:11-13`) | Geometry coded as h; issues file G3 (the lid shades risers too) |
| `s<0.0633` / `s<0.0637` | G1 0.0633; HG 0.0637 | Lid applies only to faces on the original surface | LIVE | NUM | G1: s0 + 0.1 mm (P2-8 row L1, `ARCHIVE_DONE.md` L22522; 0.0637 tied at 1 mm under face-centre evaluation). HG: s0 + dz/2 (P2-5 `MAP.txt` (iii)) | — |
| `0.001` | W/m²K | Lid value (positive floor) | LIVE | NUM | `MAP.txt` | — |
| C: `1189.763953` (H47) | = S·7272.3958, S = **0.1636** | Level | LIVE | **TUNED** | S\* rule: log-log interpolation of V(30 s) between S 0.14 and 0.20 to **33.0 cm³** (P2-11 b1, L24300; `p2_11_heat/SSTAR.txt`) | ⚠ 33.0 cm³ is her **47 MPa** volume, while `confining.p = 0`, so S absorbs the uptake gap **and** the unmodelled load (`Claude_markdowns/2026-10-30.md` §3) |
| C: `7272.3958` (H0) | S = 1 | Level at her measured h | LIVE | MEAS-RIG | h = 7.1 kW/m²K at L 8, 38 kW, D1 (`hoeser_2016_heat_transfer.csv`; Fig. 4.26), as a Ø9 mm sensing-disc mean ÷ the disc mean of the Goldstein **average** shape 0.976294 (P2-7 sources, L22072) | The thesis text gives 7.0 ± 0.8; P2-7 DR-3 (7,000 vs 7,100) is unresolved. The sensor is a point on the axis inside a Ø40 ring; h is sensor-referenced to T_f = 1,866 K |
| C: `727.239580` (H0S) | S = 0.10 | What-if level | LIVE | TUNED (what-if) | "S ~ 0.10 was chosen from her depth via the 1-D model" (`whatif_uptake/README.md` P3) | — |
| `533`, `44`, `1.285` | — | Radial shape | LIVE | LIT | Goldstein, Behbahani & Heppelmann 1986 (IJHMT 29:1227), via Zuckerman & Lior 2006: Nu_avg ∝ 1/(533 + 44(r/D)^1.285) (P2-7 sources; P2-11 sources L24215) | ⚠ Cold-air jets, Re 61k–124k, H/D 6–12. Hers: Re ≤ 35,000, underexpanded at 38 kW, 1,866 K (`thesis_notes.md:19-21`). Goldstein's own standoff factor (24 − \|L/D − 7.75\|) and its Re^0.76 are **not** used |
| `0.3575` (G1 only) | 1 − 1.285/2 | Local form (1/2r)·d(r²·g_avg)/dr | LIVE | LIT (derived identity) | P2-11 a0, independently re-derived (L24310) | HG/H0S use the area-average as a local value, which is the defect P2-11 corrected |
| `533*…/pow(…,2)` | — | Normalisation g(0) = 1 | — | math | — | ⚠ G1 kept HG's 7272.3958, which normalises the **average** shape. Over Ø9 the g_loc disc mean is 0.96149, so S = 1 under G1 gives 6,992 W/m²K, not 7,100 (−1.5 %, my arithmetic). It is absorbed in S |
| `0.0079` | m | Nozzle D1 | LIVE | MEAS-RIG | `thesis_notes.md:19` | — |
| `0.0632` (×2) | m | Reference standoff s0 = L 8 | LIVE | MEAS-RIG | 8 × 7.9 mm, her drilling standoff | — |
| `1.3185305455` | — | Standoff exponent | LIVE once s > 63.2 mm (the pit deepens) | CAL-ONCE | Fitted to her calorimeter at **L 10/11/13/14 (79.0–110.6 mm)** with h0 pinned at L 8, under K30; residuals −7.2/−6.7/+1.2/+6.1 % (`p2_5_hoeser/MAP.txt`; L21047-21056). The K20 alternative is 0.4288 | ⚠ Measured on an open target in a Ø40 ring. Inside a deep hole it is extrapolation (`2026-10-30.md` §2), and it is double-counted with the T_stag decay past s_c = 112.35 mm (ibid. §1) |

---

## D. Literals in the Meier `h_expr` ("G4")

Structure, read from the string and from `d2v_hmap/run.py:103-141`:
`h = [ if(s>0, h_imp·w1 + 142·(1−w1), 142) ]·w2 + 0.001·(1−w2)`, where
- `h_imp = A·G_loc(max(r, 2.5D), clamp(s, 2D, 12D)) / G_loc(2.5D, 50 mm) · (12D/max(12D, s))^1`,
- `G_loc = (D/2r)(r² + 2cDr − 1.1cD²)/(r + cD)²` with c = 0.1(H/D − 6), and
- `w1` and `w2` are raised cosines from 1 to 0 over r = 40→56 mm and 60→76 mm.

**G_loc is the local form of Martin (1977)'s area-averaged single-round-nozzle correlation**, (1/2r)·d(r²·h̄)/dr. I
checked this algebraically and it matches `d2a_jet_face/jet.py:10-15`.

| literal | value | meaning | LIVE? | CLASS | provenance | notes |
|---|---|---|---|---|---|---|
| `1449.0` (M) | W/m²K | Martin **local** h at r = 2.5D, H = 50 mm | LIVE | LIT (computed) | `walljet.martin(1900)`: film T (1900 + T_fire 821)/2 = 1360.5 K, Sutherland **air** properties, Re = 5.16e4 from the measured ṁ 0.01513 kg/s, D 7.5 mm. I reproduced **1448.85** (`walljet.py:42-52`; `jet.py:22-25`; L6136) | Inputs: MEAS-RIG (ṁ, D, SOD) plus the assumed T_nozzle and the model's own T_fire. ⚠ Laval nozzle ⇒ supersonic, while Martin is subsonic (reference errata 1). Meier's only measured impingement h is 0.4–1.6 kW/m²K (cold Ch. 3 calibration jets, errata 4) |
| `1075.9` (Md) | = 1449 × 0.742499 | Floor scaled to Kant | LIVE | CAL-ONCE (cross-rig) | Kant & von Rohr 2016 Table 4 interpolated at the m13 run's own median standoff 85.3 mm (984.97) over the map's value there (1326.5) (P2-3 L19674, L19735, L19952) | ⚠⚠ SR 17: a category error (a 3–10 kW threshold burner). P2-3 was REJECTED, yet Md still carries it |
| `0.0075` | m | Nozzle D | LIVE | MEAS-RIG | VT5 drawing sheet 14/17 | The thesis text says 7.1 |
| `0.01875` | 2.5D | Martin validity r/D ≥ 2.5 (plateau inside) | LIVE | LIT (validity clamp) | `walljet.py:65-72` | Martin was verified against Zuckerman & Lior 2006 (L5452) |
| `0.015`, `0.09` | 2D, 12D | Martin validity H/D 2–12 (clamp on s) | LIVE | LIT | ibid. | — |
| `2.0`, `1.1`, `0.1`, `6.0` | — | Martin G-function constants (local form) | LIVE | LIT | Martin 1977; `jet.py:10-15` | — |
| `0.19763872491145218` | — | G_loc(2.5D, H = 50 mm), the normaliser | LIVE | NUM | `walljet.py:61` (NORM); I re-evaluated it exactly | Built on the 50 mm SOD (MEAS-RIG) |
| `pow(0.09/max(0.09,s),1.0)` | n = 1 beyond 12D | Far-field standoff decay | LIVE when s > 90 mm (the median standoff is 85 mm) | FORM | D2c: the stagnation velocity gradient falls as 1/s², so h ∝ 1/s; n = 0.5 is the alternative (L7417-7423) | Described later as "hand-set" (L18608) |
| `0.04` | m | w1 start = skirt radius (`foot_r_outer`) | LIVE | MEAS-RIG | `d2v_hmap/run.py:86`; reviewer fix, `2026-09-23.md:316-318` | — |
| `0.056`, `0.016` | m | w1 end and width (8 cells at 2 mm) | LIVE | NUM — **shape-setting** | `run.py:30-31` ("never tuned") | ⚠ "The wall sits at r ≈ 44 mm, inside the 40–56 mm blend, so the hole radius is substantially set by the taper" (D2v, L12655). The string has no dx, so mesh refinement cannot move it (issues G4) |
| `142.0` (s > 0 blend) | W/m²K | Annulus plateau at r > 40 mm | LIVE | DEFAULT-UNJUSTIFIED | "Area expansion" of a 2 mm wall jet into a 12.3 mm annulus (×6.15), computed on the **D2i runaway run's** geometry (D2r-0, L11942). Band 81 (Dittus-Boelter) / 142 / 250 (`run.py:83`) | Declared and frozen, not fitted (`d2v_hmap/CRITERION.md:99`) |
| `142.0` (s ≤ 0 branch) | — | Above-plane wall | **INERT today** | — | `jet_above_h = 0` hard-zeroes h for s ≤ 0 (2112-2113, 3389) | The retired side-face path made it reachable in the accepted runs |
| `0.06`, `0.076`, `0.016` | m | w2: the "mouth" cut | LIVE | DEFAULT-UNJUSTIFIED — **shape-setting** | Chosen from the model's own LI0 mouth ring (r 50–76 mm), "floor 1e-3 beyond r ≈ 60 mm" (`2026-09-23.md:272,290`) | The issues file (G4) calls the taper/plateau/cut "hand-drawn" |
| `0.001` | W/m²K | Outside the mouth (sealed, water-cooled wellhead) | LIVE beyond 76 mm | NUM | `run.py:84` | — |
| `3.141592653589793`, `0.5`, `1.0` | — | Raised cosine | — | math | — | — |

---

## E. Would each deck run on today's source? (reasoned from source, not run)

- **Retired-key loop** (94–120): `pp.contains` on 21 keys, aborting on the first one present in list order. That is
  `spall.lateral_drive` for every deck that carries any.
- **Strict scope** (`ParmParse.H:959-972`): any `surface_patch.*` key with `m_count == 0` after `SurfacePatch::Parse`
  throws "specified but not used".

| deck | inputs | retired keys present | unread `surface_patch.*` (strict abort) | verdict |
|---|---|---|---|---|
| KO | `input_model` + harness (`test:170-180`; `P1_KEYS = []`) | none | none (constants path: enabled, mode, face, x0, y0, radius, h_conv, T_flame, T_ext are all read) | **Runs.** ⚠ A deck rebuilt from metadata would abort (`side_face_closure`, `side_face_flux`, `spall.lateral_*` are defaults from the 10/07 binary) |
| KD | `input_model` + `G2.0_130/argv.sh` | none | none | **Runs** (`constraint_depth`, `flake_length=flux` abort only under `face_list`, which is off) |
| H47 | `input_lab` + `G1/argv.sh` | none (the 13 were dropped; `DROPPED.txt`) | none | **Runs.** It ran on P2-10 C2. P2-12 removed only `face_heat`/`barrier_expr`, which G1 does not set |
| H0 | `input_lab` + `HG/argv.sh` | 13: spall.lateral_drive/removal/subcell/thermal, side_face_flux, _h, _T, _h_r_min, _h_law, _h_datum, _area_law, _step_height, _staircase | none once those are dropped | **Aborts.** With the 13 dropped it would run on the non-face-list path with **no side-face heat**, so not a replay |
| H0S | as H0 | the same 13 | none | **Aborts** (as H0) |
| M | `input_drilling` + `m13/argv.sh` | 14: spall.lateral_×4, side_face_flux, _h, _T, _h_r_min, _h_law, _factor_above, _area_law, _step_height, _staircase, step_diag_radius | none once dropped (the robin_form last-wins value is `surface`, and `pinned_idle_cycles` is read) | **Aborts.** Dropped, it runs without face-list or side-face heat, so it is not the accepted P2-2f model. Adding `face_list=1` aborts on `foot_body_clearance` (993) |
| Md | `input_drilling` + `SD/argv.sh` | 16: M's 14 plus `side_face_h_datum` and `jet_wall_march` | **`side_face_h_scale`** (no reader and not on the retired list) | **Aborts** twice over |

---

## (a) What blocks a BLIND flame prediction

1. **S, the heat-uptake factor**, is TUNED per test: 0.1636 to Höser's 47 MPa 30 s volume, and 0.10 in the what-if.
   Nothing physical sets it, and at S = 1 Höser 0 MPa removes about 7.8× her volume.
2. **The calorimeter level** h0 = 7,100 W/m²K (Ø9 disc, Ø40 ring, L 5–14) exists **only for Höser's rig**. Meier has no
   measured h: his level is a subsonic Martin correlation, evaluated with air properties and an assumed T_nozzle
   (1449), or a cross-rig Kant scale (1075.9).
3. **The standoff law** (exponent 1.3185 under K30) is CAL-ONCE on Höser's open-target calorimeter over 63–111 mm. There
   is no Meier counterpart (Meier uses Martin clamps plus n = 1), and inside a deep hole it is extrapolation and
   double-counted past 112.35 mm.
4. **`jet_core_length`, and the gas/coefficient decay split behind it**, are unmeasured: K30 vs K20 moves Höser's volume
   by 27.9 %, and Meier uses a different closure (8 D_e from cold supersonic data).
5. **`jet_T_ent_mode`** (`fixed` vs `exhaust`) alone sets Meier's floor gas at 1,742–1,804 K vs 720–831 K. It has never
   been validated, and the two decks use opposite modes.
6. **Meier's `jet_T_nozzle` = 1900 K** is unmeasured. It is not the adiabatic value (≈ 2040 K), and the hardware bound
   is 1450–1850 K.
7. **`jet_cp` = 1250** is the planner's value and has no source.
8. **The radial profile differs by deck**: Goldstein for Höser, Martin for Meier. Neither was checked against a
   measured radial profile (Höser's Gaussian k and a are unreported), and both are outside their correlation's
   Re/nozzle range (an underexpanded or supersonic hot jet).
9. **Meier's taper, plateau and cut** (40→56 mm, 142 W/m²K, 60→76 mm) are hand-set and decide where the bore wall sits.
   The 142 comes from a model run's geometry.
10. **Meier's `foot_pad_quantile` 0.9 and `foot_body_clearance` 1** are numerical workarounds that control rate and
    bore, and the bore equals 2·`foot_r_outer` by construction.
11. **Höser's `jet_De_ref`** is Meier's Laval-nozzle D_e/D ratio applied to a different nozzle. On Höser it is absorbed
    into the K30 product.
12. **The Kant anchor's operating point** (`h_conv` 500 on KO, borrowed from Kant 2016) and its `cell` Robin form differ
    from the `surface` form every flame deck runs.
13. **The losses** (ε 0.8 and h 10 to 293 K on every in-hole face) are generic defaults. They matter on low-h walls.
14. **Shape and growth data**: Höser's shape and growth targets are 47 MPa loaded data, and the model has no load
    (`confining.p = 0`).

## (b) Meier: where each heat-source input comes from

- **Measured or specified by Meier** (thesis Ch. 8 and the VT5 drawings):
  - ṁ: 2.47 kg/h CH₄ + 52 kg/h air, giving `jet_mdot`.
  - λ 1.2 and 38 kW HHV (context only; not inputs).
  - Nozzle bore Ø7.5 mm (drawing; the text says 7.1; 6 mm throat). Used for the literals 0.0075, 0.01875, 0.015 and
    0.09 and for `jet_D` (which is inert).
  - SOD 50 mm: `foot_standoff` and the 0.19764 normaliser.
  - Ø80 skirt with Ø56 bore: `foot_r_inner/outer` and the 0.04 taper start.
  - Three slots: `foot_npads` = 3.
  - The burner rests on its feet: `nozzle_descent = feet`.
- **From Höser or Kant**:
  - Nothing in the MEIER deck's heat source comes from Höser. P2-5b, the re-levelling to Höser, is PARKED. Höser's
    1,600 °C was used only to call 1900 K "consistent", and his ε 0.82 supports ε = 0.8.
  - MEIER_dref's 1075.9 is a Kant & von Rohr 2016 Table 4 rescale. SR 17 calls it a category error.
- **Set by us**:
  - Values: `jet_T_nozzle` 1900, `jet_cp` 1250, `jet_T_ent` 293.15.
  - Closure choices: `exhaust` mode, `decay`, entrained mass, momentum decay diameter, `jet_core_length` 8 (literature
    at a computed Mach number), and `jet_De_ref` (computed from an assumed T₀, γ and R).
  - The 1449 level: Martin evaluated at an assumed film T that uses the model's T_fire = 821 K, with air properties.
  - The map shape: the Martin clamps, far-law n = 1, the taper radii 56/60/76 mm, h_ann 142 and h_out 0.001.
  - Feet and surface settings: `foot_pad_quantile` 0.9, `foot_body_clearance`, `foot_rule`, `radius` 0.2, the losses,
    and `robin_form`.

## (c) Where the evidence contradicts or refines the context

1. **The exponent 1.3185** was fitted at L 10/11/13/14 (**79.0–110.6 mm**) with h0 pinned at L 8, not to L 8–13
   (63–103 mm) (`MAP.txt`). Today's log §1 repeats the "63–103" wording.
2. **`jet_mdot` 3.4e-3** is a quarter share, as stated, but it is **derived** from the 38 kW / λ 1.1 set-points with
   HHV 55.5 MJ/kg, not measured. The HHV basis of her 38 kW is unverified.
3. **`jet_De_ref` 3.745 mm** is not Höser's momentum diameter. It is 7.9 mm × Meier's D_e/D ratio (`run_p2_5.py:98`).
4. **7272.3958 normalises Goldstein's *average* shape.** G1 switched to g_loc without renormalising, so its S = 1 Ø9
   mean is 6,992, not 7,100. The thesis text says 7.0 ± 0.8 and the CSV 7.1 (P2-7 DR-3 is open).
5. **S\* matches her 47 MPa 30 s volume**, which makes S an uptake × load factor while load is not modelled.
6. **Meier's 1449 is not hand-drawn.** It is Martin's local h at (2.5D, 50 mm), film 1360.5 K, air properties, measured
   ṁ (reproduced as 1448.85). The hand-drawn parts are the taper, plateau and cut. 1075.9 = 1449 × 0.7425 (the Kant
   scale).
7. **Meier's "adiabatic" 1900 K** is ~8 % below the equilibrium T_ad of ≈ 2040 K. It is an assumed value.
8. **`surface_patch.jet_D` is inert on every flame deck** (momentum decay diameter). D enters only through the
   `h_expr` literals.
9. **The Meier map's `s ≤ 0` branch (142)** is inert in today's source, because `jet_above_h = 0`.
10. **Not every orphaned `side_face_*` key is on the retired list.** `side_face_h_scale` (Md), and the metadata-only
    `side_face_area_exact`, `_closure`, `_factor`, `_h_max` and `_r_min`, abort through the strict parser instead of
    through the retired-key message.
11. **`surface_patch.T_ext` has no consumer.** The Kant decks' 293.15 is dead.
12. **Kant**:
    - KD's 747 / 1523.15 is Kant 2016 Table 4 at SOD 130, measured indirectly with ±17 % uncertainty.
    - KO's 500 / 1673.15 was **not** measured on the Kant 2017 rig: 500 is Kant 2016's h_fl,min, and 1673.15 is Kant
      2017's "≈ 1400 °C".
    - The Kant decks run with losses **off** (ε 0, h 0). `losses.epsilon` 0.8 and `losses.h_conv` 10 belong to the
      flame decks.
13. **The Höser lid constant** is 0.0633 in G1 (s0 + 0.1 mm, P2-8 L1) but **0.0637** in HG/H0S.
14. **The accepted Meier deck (P2-2f)** and the HG deck cannot be rerun today. Their side-face and lateral keys are
    deleted, so their published numbers are not reproducible on the current binary.
