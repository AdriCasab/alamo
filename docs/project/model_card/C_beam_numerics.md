# Model card, section C: the MMW beam, MMW-side physics, numerics and validity

Compiled 2026-10-30 from today's working tree (`master`, dirty) with read-only access. No repository file was edited.
Sources I opened are cited as `file:line`. `M` = `src/Integrator/MMWSpalling.H`, `R` = `src/Integrator/MMWSpalling/Removal.H`,
`B` = `src/Numeric/MMWBeam.H`, `A` = `docs/project/ARCHIVE_DONE.md`, `RM` = `ROADMAP.md`. Deck values come from
`docs/project/model_card/RAW_KEYS.tsv`, which holds the resolved values from each run's `metadata` (argv overrides the input file).

CLASS legend: MEAS-MAT, MEAS-RIG, LIT, CAL-ONCE, TUNED, FORM, NUM, DEFAULT-UNJUSTIFIED, DIAG. I add **INERT** for a key that is
parsed but has no effect in that deck, or anywhere.

Deck order used below: KANT_onset | KANT_dref | HOSER_47 | HOSER_0 | HOSER_0_S | MEIER | MEIER_dref | MMW_ZHANG | HU_onset.

---------------------------------------------------------------------------------------------------------------------
## PART 1 — The MMW beam and the physics on the MMW side

### 1.0 Headline facts (each one verified)

| # | fact | evidence |
|---|---|---|
| F1 | The beam is **off** (`beam.P0 = 0`, no schedule) in all 8 flame and Hu decks. The only one of the nine decks with the beam on is MMW_ZHANG. | RAW_KEYS `beam.P0`. MEIER's input file has `beam.P0 = 11400` (`validation/meier/sp_meier_pilot/input_drilling:139`), but argv sets `beam.P0=0.0` (`studies/p2_2c_meier/ctrl8/m13/argv.sh`). Running that input file alone would switch on a legacy beam surrogate on top of the flame. |
| F2 | MMW_ZHANG is a **thermal-only** check. It runs on the global-`material` path, where removal is impossible (spall and vapor abort without `microstructure.nphases`, `M:615-635`). `material` together with `microstructure` aborts (`M:1082-1084`). Thermal stress is inert (`alpha = 0 0`, `el.interval = 1000000`). | the deck, plus `M:1414` dispatch |
| F3 | **Spallation has never run with the beam on against data.** Runs with the beam on and removal on do exist: unit/beam_void_closure, studies/s1_surface_resolution (09-15, old criterion), the P2-2e B2000/B100 columns (10-24, not reviewed) and the two extensions (last run 09-15). None of them is compared with a measurement. | §1.5; RM:309-313, RM:337 |
| F4 | `microstructure.phase0.alpha_attenuation` is **required by the parser but never read** (`M:335-336` is its only use). The beam's absorption is `beam.alpha_expr` alone, one global α(T), not per phase. | grep of `alpha_atten` over `src/` returns only `M:335-336` and the struct member `M:8125`. A deck comment agrees: `extensions/hu_spall_onset_mmwbeam/input_granite2:56`. |
| F5 | **α(T) deposition with energy closure is NOT done.** `BeamSource` uses the cell's own `a = alpha_f(T_old)` with `exp(-a·path)` (`M:4632, 4643-4646`), so the column sum does not telescope when α varies. The code says so at `M:4423-4427`: "a correct alpha(T) deposition needs the cumulative optical depth tau = sum_k a_k*dz (follow-on step, required before real alpha(T) MMW work)". No test anywhere uses a T-dependent `alpha_expr`. | ⚠ The pointers `RM:281` and `DIRECTION_2026-10-14…:103` (`MMWSpalling.H:6166(-6169)`) are **stale**. Line 6166 is now inside the Sp cluster labeller. |
| F6 | Today's removal rule (`surface_patch.face_list = 1`, the HOSER_47 deck) **aborts** with `vapor.enabled`, `spall.melt_aware`, `constraint_depth > 0`, any non-`sp_weibull` model, or a stress source that is not `local_thermoelastic` (`M:991, 995, 997-1007`, "that is the MMW step's"). It also requires a `convective_flame` patch with an `h_expr` (`M:952-955`). A beam-only MMW deck therefore cannot use today's removal rule. | `M:949-1020` |
| F7 | The surface temperature reconstruction that the criterion fires on assumes a profile that decays from the face. `RawSurfaceX` returns 0 ("flat") whenever `T_1 >= T_top` (`M:8750-8760`). Under volumetric absorption with an interior peak, the criterion falls back to the cell average or to a held x. RM:360-363 states the same thing: "neither breaking rule is formulated for volumetric heating". | `M:8730-8800` |

### 1.1 `beam.*` keys

Parse: `B:35-89`. Deposition kernel: `BeamSource` `M:4618-4647`. Incident intensity: `BeamIncidentIntensity` `M:4655-4686`.
Call sites: `AdvanceConstAlpha` `M:1457`, `AdvanceMaterial` `M:1577`, `AdvanceMaterialImplicit` `M:1653`, `AdvanceMicrostructure` `M:2005`.
All of them pass `dz_cell = DX[2]`, so deposition is the analytic cell integral. Only the microstructure path passes the void-aware path length (`SolidFacePath`, `M:4028`).
Closure counter: `BeamEnergyBalance` `M:4448`. Vapour use: `R:323-336`.

| key | MMW_ZHANG | other 8 decks | read at | live when | CLASS | provenance (opened) | notes |
|---|---|---|---|---|---|---|---|
| `beam.P0` | 0 (schedule used) | 0 (MEIER input 11400 → argv 0) | `B:39`; `power_at` `B:93-99` | P ≠ 0 at time t | MEAS-RIG | Zhang summary `validation/zhang/zhang_oglesby_mmw_validation_summary.txt:53-68` ("averaged incident powers… replacing the noisy experimental power trace") | `MMWBeam::active()` (`B:126`) tests `P0` and ignores the schedule. It is unused (grep), so this is harmless today. |
| `beam.schedule.t_end` / `.P` | 150…1750 s / 800,1800,2300,2800,1700,2800,1700,3200,0 W | – | `B:83-86`, `B:93-99` | if present (P0 ignored) | MEAS-RIG | summary :58-68 (copied exactly) | Summary :137-139 treats these as incident power at the rock (zero standoff). |
| `beam.omega0` | 0.02 m | 0.02 (inert) | `B:41`; Gaussian `M:4674-4677`; vapour `A_beam` default `M:1066`; vapour Q_local `R:1713-1734` | Gaussian profile, P ≠ 0 | MEAS-RIG/LIT | summary :44-51 ("earlier draft 0.20 m from Zhang Table 4 was a mistake; intended 0.02") | Zhang 2023 is **not on disk**, so the 0.02 vs 0.20 correction is unverified against the paper. The default of 1.0 m is DEFAULT-UNJUSTIFIED. |
| `beam.lambda` | 2.19e-3 m | 2.19e-3 (inert) | `B:43`; only in `zR = π ω0²/λ` at `M:4674` | Gaussian only | MEAS-RIG (assumed) | deck comment: "c / 137 GHz gyrotron" | ⚠ The summary names 137 GHz only as the **radiometer** frequency (summary :12). The source frequency is not stated in the reference. `zR = 0.574 m` ≫ the 25 mm block, so λ is effectively inert for Zhang (waist change ~0.1 %). λ does not set α: the code has no α(f) link. |
| `beam.R` | 0.2 | 0.0 | `B:52`; `(1-R)` in `M:4643/4646`; vapour `R:324` | P ≠ 0 | **TUNED** | summary :133-139 says **R = 0** and "no (1−R) factor is applied". A Step 4b record: "Tuned settings reproduce Zhang's simulation curve" (`A:246-249`) | Constant: no R(T), no R(angle), no R(melt). Degenerate with P and α in the absorbed power. |
| `beam.alpha_expr` | 300 m⁻¹ | 1.0 (Kant, Hu), 2000 (Höser, Meier), all inert | `B:73-77` (parser in T); evaluated per cell at `T_old`, `M:4632` | P ≠ 0 | LIT, secondary | deck comment → `main.tex:359-366` (α ∈ {10,30,100,300,1000} ↔ 28/45/94/140/300 GHz "for granitic rock", citing woskov2009mmw) | ⚠ The summary gives **no α** (grep of the summary: only "alpha_lambda", :135). Woskov 2009 is **not on disk**. A T-dependent expression is accepted but not energy-closed (F5). Default "10" is DEFAULT-UNJUSTIFIED. |
| `beam.profile` | gaussian (default) | gaussian / uniform | `B:60-64` | – | FORM | `B:54-59` (uniform was added as a flame-jet surrogate) | ⚠ The vapour rate under `sp_weibull` always uses the **Gaussian** waist formula, even when `profile = uniform` (`R:1713-1734`). |
| `beam.radius` | 0 | 0 / 0.02 / 0.042 | `B:66-68` | uniform only | FORM | – | Uniform requires > 0 (abort). |
| `beam.x0/y0/z0` | 0.05/0.05/0.025 | various (inert) | `B:45-47` | always (geometry) | FORM | summary :37-42 (beam centred, zero standoff) | `z0` is the beam origin. Deposition is zero above `z0` (`M:4663`). |
| — | — | — | — | — | FORM | `M:4655-4686`, `M:1955-2012` | **Beam geometry limits (by reading the code):** propagation is vertical (−z) only, column by column. Power enters each (i,j) column at its top solid cell. There is no oblique incidence, no wall absorption, no reflection or waveguiding in the bore, and no vapour or plasma shielding. Power that crosses the domain bottom leaves as `P_trans`. |

### 1.2 `material.type = zhang` (`src/Numeric/Material/Zhang.H`). All coefficients are parser defaults; MMW_ZHANG overrides none of them.

Base class: `Material::effective_rhocp` adds latent-heat boxcars (`src/Numeric/Material/Material.H:30-43`).
H↔T table: `src/Numeric/Material/Table.H:29-58`, range **200–5000 K, dT = 0.5 K**. T clamps silently outside that range (`Table.H:103-106, 122-123`).

| constant | value | Zhang.H line | CLASS | provenance | note |
|---|---|---|---|---|---|
| ρ(θ) = a·θ + b | a = −0.2877, b = 2884 kg/m³ | :39-40, :67-71 | LIT | summary :74-78 (Zhang 2023 Table 1/Eq. 10; paper not on disk) | θ = T − 273.15 |
| ρ cap | 2300 above T = 1273 K | :41-42 | LIT | summary :77-78 | **Discontinuous**: 2596.3 → 2300 kg/m³ at 1273 K (−11 %). |
| κ(θ) = 1/(a + bθ) | a = 0.3154, b = 3.8e-4 | :45-46, :73-77 | LIT | summary :80-82 | 3.08 W/mK at 297.5 K; 1.44 at the cap |
| κ cap temperature | 1273 K | :47 | DEFAULT-UNJUSTIFIED | summary :86-89: "Table 1 does not give the explicit threshold. We adopt 1273 K" | – |
| c_p low branch | −2.73e-4 θ² + 0.894 θ + 788 | :50-52, :79-84 | LIT | summary :95-96 | – |
| c_p high branch | 0.184 θ + 1249 | :53-54 | LIT | summary :98-99 | Continuous to within 0.5 J/kgK at the break |
| c_p break = T_m | 1510.5 K | :55, :59 | LIT-derived | summary :91-93 (midpoint of 1488–1533 K) | – |
| fus_width | 200 K | :60 | LIT | summary :105 | Boxcar L_m/200 K on 1410.5–1610.5 K |
| L_m | 3.4e5 J/kg | :61 | LIT | summary :104 | – |
| T_vap_lo / T_vap_hi | 3233 / 3503 K | :62-63 | LIT | summary :107 | – |
| L_v | 4.8e6 J/kg | :64 | LIT | summary :106 | Spread over 270 K |
| T_ref | 273.15 K | :58 | FORM | – | Datum where H = 0 |

### 1.3 Melt and vapour machinery

| key / item | value (decks) | read at | live when | CLASS | provenance | notes |
|---|---|---|---|---|---|---|
| `phase.T_ref/T_m/T_v/lh_m/lh_v` | 0 / 1e30 / 1e30 / 0 / 0 in all 9 | `M:145-149`; used only in the const-α branch, `M:6386-6392, 6433-6440, 6483-6510` | no `material` and no `microstructure` (unit/stefan) | INERT in all 9 decks | – | Sharp-plateau Stefan model in K units (lh = L/c_p). Witness: unit/stefan, front error 4.23 % (`unit/stefan/sweep_run.log`). |
| `microstructure.phase0.T_m, fus_width, L_m, T_vap_lo, T_vap_hi, L_v` | 1e30 / 0 / 0 / 1e30 / 1e30 / 0 in all 7 flame decks and HU. Not applicable to ZHANG. | `M:376-381`; copied into the phase table `M:6788-6820` (Grimsel too, `M:6809-6817`) | always (table) | DEFAULT. **DEFAULT-UNJUSTIFIED for MMW** | – | **No melt or vapour is possible in any flame deck.** A cell can heat to any T with no phase change. RM:308 reports m = 8 Meier columns reaching **1395 K**, inside the 1200–1450 K solidus band (DIRECTION §2). Extensions copy Zhang's scalars (`extensions/*/input_granite2:88-93`). Witness: `unit/latent_heat`. |
| `spall.melt_aware` | 0 in all | `M:5219-5220`; `R:1027-1033, 1375-1379, 1479` (column depth scan only) | sp_weibull column path | INERT | P1 record `A:2758-2823` | Scales the driving stress by f(Λ_L) = max(0, 1 − Λ_L/Λ_crit) (`M:8362`). Not applied at the Sp_field onset site. **Aborts under `face_list`** (`M:995`). P1 takeaway 1: with α constant a spalling surface never reaches T_m ("spall thermostat"). |
| `melt.lambda_crit` | 0.4 | `M:5221-5226` | melt_aware | LIT (unverified) | "rigidity percolation" (`A:2764`); no source on disk | INERT |
| `vapor.enabled` | 0 in all 9 | `M:5516-5518` | – | INERT | Step 11 `A:999-1066` | **The vapour front is not exercised by any representative deck.** It aborts under `face_list` (`M:991`). |
| `vapor.A_beam` | default 0.5·π·ω0² | `M:1059-1069`, `M:5520-5527` | vapor on, damage_law branch | FORM | `M:1061-1064` | Gaussian effective area |
| `vapor.T0` | default `losses.T_amb` | `M:1071-1072`, `M:5528-5534` | vapor on | FORM | – | – |
| Vapour rate law (Step 11/17) | v = (1−R)·Q_local/(ρ·L_tot), with L_tot = c_p(T_vap_lo − T0) + L_m + L_v. Q_local = 2P/(πω0²)·exp(−2r²/ω0²) at the waist. Fires when H_top ≥ H(T_vap_lo) | `R:1688-1755`; spall has priority under sp_weibull, `R:1985-1990` | vapor on | FORM | `A:1020-1022` | **By reading the code:** (i) c_p is the top-cell value, constant over the whole span (`R:1692-1696` flags it); (ii) radiation and conduction losses are not in the balance (εσT⁴ ≈ 5 MW/m² at 3300 K); (iii) the uniform profile and beam spreading are ignored; (iv) `BeamSource` keeps depositing the same power into H while the rate formula also spends it, and no regime-2 ledger term was found (grep). This looks like double use of the beam power. **Unverified by any test**: the witnesses are synthetic (unit/regime_low_high_power uses ρ = c_p = 1 and α = 0; unit/sp_v_n_regime uses T_vap_lo 700 K and L_v = 0). |
| Melt evacuation (purge or shear) | – | – | – | not modelled | `validation/meier/TODO.md:50-56` (P2, not done) | Melt stays in place as a spall-proof blanket. |

### 1.4 `losses.*` in the MMW deck

Parse: `M:153-161`. Applied explicitly on the top-row `surface_mf` at the **cell** temperature (implicit path `M:1657-1661`; explicit material path `M:1582`).

| key | MMW_ZHANG | other decks | CLASS | provenance | notes |
|---|---|---|---|---|---|
| `losses.epsilon` | 0.8 | 0 (Kant, Hu); 0.8 (Höser, Meier) | **TUNED** (Zhang) | summary :121-126 gives **0.7 for T ≤ T_m and 0.3 above** | ⚠ The validation **metric** also uses this ε. The test reports ε(T)·T_raw (`zhang_oglesby/test:119-126`), while the summary's radiometer correction factor is 0.7 (:13-14). |
| `losses.epsilon_high` | 0.8 | = epsilon | TUNED | as above | – |
| `losses.T_eps_break` | 1510.5 K | 1e30 | LIT-derived, INERT | = T_m | No effect because ε = ε_high. |
| `losses.h_conv` | 51 W/m²K | 0 / 10 | **TUNED** | summary :128-131 gives **56**, itself "tuned rather than a free-convection value (5–25)". Changed 56 → 51 in commit d9f098b96 (2026-06-10); the last run's metadata already shows 51. | – |
| `losses.T_amb` | 297.5 K | 293.15 / 280.15 | MEAS-RIG | summary :34-35 | – |
| `losses.sigma_SB` | 5.670374419e-8 | 0 or 5.670374e-8 | FORM (physical constant) | – | – |

### 1.5 What MMW validation exists

| case | beam | removal | compared to data? | last run (binary) | in the 2639 tree gate? | result |
|---|---|---|---|---|---|---|
| `validation/zhang/zhang_oglesby` | ON: Gaussian, schedule, α 300, R 0.2 | **OFF** (material path) | Yes, digitised Zhang 2023 Fig. 7 (simulation and Oglesby experiment), on ε·T at the beam centre | **2026-05-30**, `04e11e99c-dirty` (output metadata) | **No** (`studies/tree_sweep_0922/witness/wit.sh:10-14`) | Step 4b PASS: "peak corrected surface temperature 2871 K vs Zhang's 2850 K" (`A:246-249`). Bars: peak 15 %, plateau 250 K, cooling 200 K (`test:8-15`). ⚠ The raw model T is the plotted value / 0.8, so the peak is ≈ **3590 K**, above T_vap_hi = 3503 K, with no removal. The plotfiles are deleted; only the PNG and metadata remain. It has not been re-run since the BeamSource integration fix, the Table perf change and the later source changes. |
| `unit/beam` | P0 1 W, α 10, Gaussian | off | analytic | outputs 10-29 | No | Beer–Lambert integral 0.00 %; hot spot on axis (`unit/beam/sweep_run.log`) |
| `unit/equilibrium` | flat beam | off | analytic radiation equilibrium | 10-28 | No | 0.00 % |
| `unit/beam_void_closure` | uniform 2.06 MW/m², α 2000, 2 mm | **ON** (sp_weibull) | energy closure only | 10-30 | **Yes** | T1 \|closure_err\| ≤ 1e-9; audit closure 1.00000 (`sweep_run.log`) |
| `studies/s1_surface_resolution` (09-15) | 1-D, 2.06 and 15 MW/m², α 100–32000 | ON (pre-P0 cell-average criterion, constant properties) | No | 09-15 | smoke only | ROP converged at fine dz; firing T 655–667 K; "no melt anywhere incl. M-100" (`RESULTS.md`, `A:3297-3302`) |
| `studies/p2_2e_cref` B2000 / B100 | 1-D 1 mm, 5 MW/m², α 2000 / 100, Grimsel, P2-2d criterion | ON | No (P2-2e **not reviewed**, `A:21572-21578`) | 10-24 | No | key 0 fires at 847.95 K (B2000) and 841.60 K (B100). With d_ref = 2 mm: 871.62 K and **no fire** (69.5 % negative brackets) (`CHECKS.txt:26-35`) |
| `extensions/hu_spall_onset_mmwbeam(_lefm)` | 10 kW, ω0 20 mm, α 100, R 0.2, Zhang losses and melt/vapour scalars | ON (spall + vapour) | No ("model extension, NOT a validation", `A:2662-2663`) | 2026-09-15 (`a463edd9c-dirty`) | No | "Drilling is mostly vapor (T_max 4233 K > T_vap_lo 3233 K)" (`A:2668`). Pre-P0, pre-P2-2. |

**Verdict:** the MMW validation is thermal only, on one case, with three tuned surface parameters (R, ε, h). It was tuned to Zhang's *simulation* curve, the absorption coefficient is borrowed from a secondary table, and it is not regression-protected. Spallation, melt and vapour under the beam have no data comparison at all. RM:309-313 and RM:337 are confirmed.

### 1.6 What an unseen MMW test would need, and what the code accepts today

| quantity | needed for a blind prediction | what the code accepts today | gap |
|---|---|---|---|
| Incident power P(t), spot, standoff/waist, frequency | MEAS-RIG | `beam.P0`/schedule, ω0, λ, z0 | Waveguide-to-rock geometry is reduced to a Gaussian waist at `z0`. There is no bore waveguiding. |
| Absorption α(T, f) (dielectric ε′, ε″ or tan δ vs T) | measured for the rock at the gyrotron frequency | one global `alpha_expr(T)`. Per-phase `alpha_attenuation` is inert (F4). | T-dependence is not energy-conserving (F5). No α(f) link. No data on disk (Woskov, Zhang, Oglesby all absent). |
| Reflectivity R(T, angle, melt) | measured or Fresnel from ε′ | constant `beam.R` | Tuned to 0.2 on Zhang (reference says 0) |
| Surface emissivity ε(T) | measured | two-level step function | Tuned to 0.8 on Zhang (reference 0.7/0.3) |
| Above-surface convection / purge gas | rig | one constant `h_conv` | Tuned 56 → 51 |
| ρ, c_p, λ up to 1500–3500 K | MEAS-MAT | Grimsel law (measured 25–500 °C) extrapolated without a clamp (`M` UpdatePropertyFields `g.lambda(Tc)`, `g.cp(Tc)`). The Zhang law runs to 1273 K but cannot be used with removal (F2). | No composable high-T law on the removal path (Step 17b deferred, `A:2352-2360`, `A:17132-17141`) |
| Melt: solidus/liquidus, L_m, melt removal | MEAS-MAT | per-phase scalars, defaulted off in every deck; `melt_aware` skeleton | No purge or shear removal; melt_aware is not carried by `face_list` |
| Vapour: T_vap, L_v, ablation-rate data | MEAS-MAT plus data | Step 11/17 front | Untested against data (RM:328, "MMW ablation-rate data (needed for D)" absent); possible double use of beam power (§1.3) |
| E(T), β(T), K_Ic(T) above 650 °C; softening s(T) | MEAS-MAT | Grimsel Table 3 flat ≥ 600 °C (`src/Numeric/Material/Grimsel.H:78-85, 280-290`); s(T) not implemented (DIRECTION step B) | Above-clamp behaviour is parked (RM:296-308) |
| A failure law valid for interior-peaked heating | FORM | P0 reconstruction plus the absolute or d_ref criterion | F7; RM:360-363 |

---------------------------------------------------------------------------------------------------------------------
## PART 2 — Numerics and validity

### 2.1 Numerical setup per deck (resolved values; `[in: …]` is the input-file value that argv overrode)

| deck | cells (mm) | domain (mm) | symmetry | n_cell | dt (s) | stop (s) | solver / path | max_level | ledger | follow_mask | IC (K) | BC (T) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| KANT_onset | 1.0 cubic | 24×24×24 | full; patch centred, r = 50 mm covers the block, so effectively 1-D | 24³ | 0.0125 | 45 [in: 90] | explicit, microstructure (law constant) | 0 | 0 | 1 | 293.15 | Neumann-0 ×6 |
| KANT_dref | 1.0 | 4×4×150 | 1-D column | 4×4×150 [in: 24³] | 0.0125 | 200 [in: 90] | explicit, microstructure (Grimsel) | 0 | 1 | 1 | 293.15 | N0 ×6 |
| HOSER_47 | 1.0 | 50×50×150 | **quarter** (patch at x0 = y0 = 0) | 50×50×150 | 0.008 | 120 [in: 30] | explicit (Grimsel), face_list = 1 | 0 | 1 | 1 | 293.15 | N0 ×6 |
| HOSER_0, HOSER_0_S | 1.0 | 50×50×150 | quarter | 50×50×150 | 0.008 | 30 | explicit (Grimsel) | 0 | 1 | 1 | 293.15 | N0 ×6 |
| MEIER, MEIER_dref | **2.0** | 120×120×600 (z −100…500) | quarter | 60×60×300 [in: 70×70×60] | 0.016 [in: 0.001] | 500 [in: 56] | explicit (Grimsel) | 0 | 1 | 1 | 293.15 | N0 ×6 |
| MMW_ZHANG | 0.781 cubic | 100×100×25 | full (beam centred) | 128×128×32 | **2.0** | 1750 | **implicit**, material path | 0 | n/a (needs microstructure, `M:661-664`) | n/a | 297.5 | Neumann, hard-coded in MLMG (`M:1687-1696`) |
| HU_onset | 3.125 | 100×100×100 | full (patch r 15 mm) | 32³ | 0.8 | 45 | explicit, damage_law, el static | 0 | 0 | 0 | 280.15 | N0 ×6 |

### 2.2 Numerics keys

| key | where read | live | CLASS | evidence the value is adequate | notes |
|---|---|---|---|---|---|
| `timestep` | ALAMO Integrator base | always | NUM | **D2c dt ladder**: R7b_2mm at 4/8/16/32/64/128 ms. 8 and 16 ms PASS (ROP −0.1 / −1.0 %, t_bottom +0.4 / +1.2 %). 32 ms FAILS (ROP +4.4 %, min Ø −2.4 mm), and the failure mode is the flake growing 1.01 → 1.11 mm. Rule adopted: overshoot q·dt/(ρc_p·dz) ≲ 6 K/step, giving **16 ms at 2 mm and 8 ms at 1 mm** (`studies/d2c_steady/RESULTS.md:611-634`; rule code `run.py:91-95`). | ⚠ The ladder was run once, on the D2d-era deck (pinned closure, pre-P0 criterion, constant properties) and at 2 mm only. "8 ms at 1 mm" is derived from the rule, not laddered. **No code CFL guard** exists (grep `stab/courant/dt_max` in `M`, `R`, `HeatConduction.H`: none). Margins are in §2.4. |
| `amr.n_cell` + `geometry.prob_lo/hi` | Integrator | always | NUM | §2.3 | Face list needs cubic cells (`M:3417`) |
| `amr.max_level` (+ `regrid_int` 1e6, `blocking_factor` 2, `grid_eff`, `n_error_buf`) | – | – | NUM | 0 everywhere. AMR is parked by user instruction (RM:89-90, SR 12). | If ever turned on: `P_now = owner ? … : 0` (`M:1782`) drops beam deposition outside the finest band. `CheckSurfaceCovered` ties the band to `h_spall_max` (`M:4163`), not to 1/α, so nothing guards beam + AMR. The ledger aborts at max_level > 0 (`M:665-667`). `RefillTempGhosts` uses FillBoundary, not FillPatch (RM:349-350). |
| `amr.max_grid_size` | – | – | NUM | Decomposition identity: P0 (e) gives removal_events.csv bit-identical across 4 (max_grid_size, ranks) layouts (`A:13964-13965`, `A:14075`). P2-2e a3: event logs byte-identical for np 1/2/4 and mgs 8/25/50 on Kant and Höser (`studies/p2_2e_cref/GATE_a3.txt`). P2-9 a3(iv): np1/np4, mgs 8/16 (`studies/p2_9_face_removal/A3_A6_CHECK.txt:20-23`). D2e: mgs 60 agrees except ledger_err at 4.4e-16 (`A:8970`). | SR 16: a zero pattern that follows the domain split is a decomposition bug first. SR 9: the toolchain is not byte-reproducible across rebuilds. |
| `solver` | `M:248-251`; dispatch `M:1411-1416` | implicit **only on the global-material path**. Every spallation deck is explicit forward Euler. | NUM | none for Zhang | The comment at `zhang_oglesby/input:23-26` is stale ("Forward Euler… 0.1 s"). |
| `implicit.tol_rel/abs/max_iter/verbose` | `M:252-255`, solve `M:1714` | Zhang only | NUM | none | ⚠ By reading the code: A = ρc_p,eff(T_old)/dt (`M:1648`) is a single linearisation, and H is reset from T_new through the table (`M:1720-1726`). A step that crosses the edge of a 200 K fusion boxcar is not enthalpy-conservative. There is no ledger on this path, and dt = 2 s has no ladder. Source and losses are explicit (`M:243-245`). MLMG non-convergence behaviour is unverified. |
| `surface.follow_mask` | `M:646-664` | removal on | FORM | P0 takeaway 10: **mandatory** with removal and a surface patch, otherwise heating stops at the first void (`A:14036-14041`) | Not supported with prescribed_T (Hu = 0) |
| `surface.refine_depth/height` | `M:675-686`, tagging `M:4688-4710` | max_level > 0 | NUM | INERT (0) | – |
| `energy_ledger.enabled` | `M:648, 661-667` | microstructure, single level | DIAG | "exact only with Neumann-0" (S1 goal, `A:3336`, `A:3711`). All decks are Neumann-0. Whole-domain closure ~1e-14 on S1 (`studies/s1_surface_resolution/RESULTS.md`) | SR 7(ii): `ledger_err_solid` is relative. Without `removal_energy = flake` it reads ≈0.98 on firing columns by construction (`studies/p2_2e_cref/CHECKS.txt:66-71`). |
| `hc.bc.temp.*` | `M:125` (HeatConduction parse); domain ghosts feed the stencil (`M:1888-1897`) | explicit paths | FORM | Neumann-0 on all six faces in all 9 decks. Quarter decks: x = 0, y = 0 are symmetry planes, x/y high and the bottom are adiabatic. | The base is adiabatic and a leg can drill into it (RM:433-435). Meier m13 margin is 71 mm at 500 s (`A:19414`, `A:19507`). Höser block is 150 mm; P2-5 DR-8 asks for 200 mm (`A:20993`). |
| `hc.ic.*` | HeatConduction | t = 0 | MEAS-RIG/FORM | = T_amb of the case. It is also the stored stress-free reference T_init (`M:1088-1095`). | – |
| `hc.heat.alpha` | `M:141` | const-α path (`M:1432`); damage_law h_spall (`R:293`) | INERT on sp_weibull and material paths | Zhang deck comment :83-84 | – |
| `stop_time` | – | – | FORM | Zhang 1750 s = full schedule; Meier 500 s; Höser 30 or 120 s; Kant 45 or 200 s; Hu 45 s | – |
| `amr.thermo.plot_int` (6 on Höser and Meier), `amr.plot_int` | – | – | DIAG | – | SR 7(i): never sum a per-step thermo column |

### 2.3 Convergence evidence

| study | what varied | result | pointer | does it cover today's deck? |
|---|---|---|---|---|
| P0 Kant cube | dz 2/1/0.5 mm | firing T_surf spread **0.161 % / 0.116 %** (h 150 / 500). Bar ±0.5 % (`TOL_MESH = 0.005`, `validation/kant/sp_kant_onset_meiercfg/test:83`) | `A:13960`; P0b (d) `A:14075` | Firing T only, P0-era criterion. KANT_onset today is a 10-07 run carrying unused `differential` keys. Not re-scored after P2-2. |
| P0b 1-D column | 2/1/0.5 mm | rate spread **1.51 %** (flake bound on) | `A:14072-14073` | Column only |
| P2-2e c3, Kant SOD 80, d_ref | 2/1/0.5 mm, dt 12.5 ms fixed | 915.8 / 890.4 / 888.1 K. Spread 4.58 % overall; **1 ↔ 0.5 mm 0.38 %**. 2 mm is +25 K because d_ref is one cell there. | `studies/p2_2e_cref/CHECKS.txt:21-25, 73-82` | KANT_dref (1 mm): yes. MEIER_dref (2 mm): under-resolved. |
| P2-9 c4, flat floor T00 (h 700, T_gas 1866 K) | 2/1/0.5 mm (dt 8/8/**4** ms) | 6.681 / 6.132 / 6.027e-4 m/s: **10.4 %** overall, 1.7 % from 1 to 0.5 mm, the same on both keys | `studies/p2_9_face_removal/A3_A6_CHECK.txt:89-100`; `A:23559`; RM:345-346 | Mesh and dt changed together at 0.5 mm |
| P2-9 c3, Höser 0 MPa M0, 30 s | 2/1/0.5 mm | V(f): 2 mm +5.8 %, 0.5 mm +5.4 % vs 1 mm. Max Ø +1.2 / +2.6 mm, axis 0.0 mm | `studies/p2_9_face_removal/MESH.txt:15-16` | 30 s and 0 MPa only. **No mesh twin for HOSER_47 (G1, 120 s) or HOSER_0 (P2-7 HG)** (`studies/p2_11_heat/legs`, `p2_7_hoeser_heat/legs`). |
| P2-5 rev 2, Höser G05 | 1 vs 2 mm (dt 8 vs 16 ms) | volume 1.50 % | `A:20975-20978` | Pre-face-list deck (side-face path, deleted in P2-10) |
| S1, 1-D beam and Robin | dz 2 → 0.0625 mm | beam ROP converged from 2 mm (F1 9.40 vs 9.37 m/h); Robin not converged at 2 mm (+42 % / −82 %) | `A:3297-3302`; `studies/s1_surface_resolution/RESULTS.md` | Old criterion |
| D2g-1 / D2i, Meier | 2 vs 1 mm | **FAILED**: −28.5 % and −13.9 % (matched 150–57x s) | `A:10203-10215`, `A:10553-10565` | Older deck. Since then the 1 mm run is parked (SR 12), and P2-2f c4 says "no absolute 2 mm rate is quotable" (`A:19369`). **MEIER has no current mesh partner.** |
| D2d dt check | 16/8/4 ms | ±0.5 % (gap seed is not time stepping) | `A:8287+`; memory `mmwspalling-d2d-mesh-diagnosis.md` (background) | Old deck |
| Hu onset | 40³ → 32³ | granite onset 37.0 → 37.6 s | `A:1127-1137` | Two points only. A_D is calibrated ("v1"). |
| Zhang | – | **none** (mesh or dt) | – | – |
| Weibull legs | mesh | onset moves ~3.7 % by construction (per-cell draws, V0 = V_cell at 2 mm) | RM:401-403 | Höser and Meier (weibull on) |

### 2.4 Explicit-step margins (my arithmetic; no code check exists)

Stability limit: dt_max = dx²/(6κ_max). Overshoot: the D2c rule, q·dt/(ρc_p·dz) at the firing-time flux.

| deck | κ_max (m²/s) | dt_max (s) | dt / dt_max | flux estimate | overshoot (K/step), bar ≈ 6 |
|---|---|---|---|---|---|
| KANT_onset (constant: 1.5/(2750·790)) | 6.90e-7 | 0.241 | 0.05 | ≤ 0.69 MW/m² at t = 0 | ≤ 4.0 |
| KANT_dref (Grimsel at 293 K) | 1.71e-6 | 0.098 | 0.13 | 0.92 → 0.44 MW/m² | 5.7 → 1.8 |
| HOSER_0 (h0 7272) | 1.71e-6 | 0.098 | 0.08 | ≈ 7.3 MW/m² (h0 × (1866 − 862)) | **≈ 19–27** (Grimsel vs constant ρc_p): **outside the ladder** |
| HOSER_47 (h0 1190) / HOSER_0_S (727) | same | 0.098 | 0.08 | 1.2 / 0.73 MW/m² | 3.1 / 1.9 |
| MEIER (2 mm) | 1.71e-6 | 0.390 | 0.04 | ≈ 1.2 MW/m² (RM:387) | ≈ 3 |
| Kant at 0.5 mm (P2-2e M80_05) | 1.71e-6 | 0.024 | 0.51 | – | – |
| HU_onset | 1.14e-6 | 1.42 | **0.56** | prescribed T | – |
| MMW_ZHANG (implicit) | 1.32e-6 | (explicit 0.077) | dt = 26× the explicit limit | absorbed peak 1.0–4.1 MW/m² | n/a |

### 2.5 Validity envelope

| quantity | range where the evidence supports the model | evidence | outside that range |
|---|---|---|---|
| MMW thermal response (no removal) | Granite; Gaussian ω0 = 20 mm; incident peak 1.27–5.1 MW/m² (absorbed 1.0–4.1 with R 0.2); α = 300 m⁻¹ constant; 0–1750 s; reported ε·T 300–2870 K | Step 4b `A:231-250`; test bars `zhang_oglesby/test:8-15` | Untested. The three surface parameters are tuned to *this* case. Raw T exceeds the vapour window with no removal. |
| Beam deposition energy closure | constant α, any α·dz (analytic cell integral); exact for the uniform profile; Gaussian to O((dz/z_R)²) | `M:4415-4430`; unit/beam 0.00 %; beam_void_closure ≤ 1e-9 | α(T): not closed (F5). AMR: deposition lost outside the band (`M:1782`). |
| Absorption depth with spallation on | 1/α 0.5–10 mm, 1-D columns at 5 MW/m², **model-only** (no data) | P2-2e c4 (`CHECKS.txt:26-35`); S1 1/α 0.03–10 mm (old criterion) | Interior-peaked profiles make the reconstruction "flat" (F7). With d_ref = 2 mm, 1/α = 10 mm never fires (by design). |
| Melt | none in any removal deck (T_m = 1e30) | RAW_KEYS; P1 `A:2790-2795` (thermostat: melt unreachable at constant α) | Melt-aware and vapour abort under face_list. Melt evacuation is not modelled. |
| Vapour removal | synthetic unit logic only | unit/regime_low_high_power, unit/sp_v_n_regime | Untested against ablation data; energy consistency unverified (§1.3) |
| Mesh, firing temperature (surface heating) | 0.5–2 mm: ≤ 0.16 % (P0). With d_ref: 1–0.5 mm 0.38 %, 2 mm +25 K | §2.3 | > 2 mm untested; d_ref = 1 cell at 2 mm |
| Mesh, floor rate | 1–0.5 mm: 1.7 %; 2 mm: +9 % vs 1 mm (10.4 % overall) | P2-9 c4 | No cross-mesh claim tighter than about 10 % that is referenced to the absolute floor rate (RM:345-346) |
| Mesh, 3-D cavity (Höser) | 2/1/0.5 mm at 30 s: V within 5.8 %, Ø within 2.6 mm | P2-9 c3 | 60–120 s and the 47 MPa / G1 deck: no twin |
| Mesh, Meier ROP/shape | **none**; 2 mm only | D2g-1, D2i failed; SR 12 parks 1 mm; `A:19369` | Every Meier rate is provisional |
| Thermal skin vs cell | Meier and T00 floor ℓ = κ/v ≈ 1.6–1.9 mm, under-resolved at 2 mm (RM:345-346); ≈ 2 mm at 1.5 m/h (`A:21927`). Höser 1 mm: top-cell centre is 214–227 K below the reconstructed face (`A:22240-22243`); my exponential inversion gives ℓ ≈ 0.85–0.95 mm, about one cell. | as cited | MMW with 1/α < dz: the skin is sub-cell. The reconstruction clamps ℓ ≥ dz/4 (`M:8781`). |
| Timestep | ≤ ~6 K/step overshoot (16 ms at 2 mm on the D2c deck) | d2c ladder | HOSER_0 at ≈ 19–27 K/step is outside. Zhang dt 2 s is unstudied. |
| Property temperature range | Grimsel c_p measured 25–500 °C, extrapolation recommended only to 873 K; λ₁ to 500 °C; Table 3 to 650 °C, with E/E₀ and K_Ic flat ≥ 600 °C (873.15 K). Zhang κ and ρ capped at 1273 K. Enthalpy table 200–5000 K. | `validation/grimsel_properties/grimsel_granite_properties_reference.md:38-63, 86-118`; `Grimsel.H:78-85`; RM:296-308 (Meier m13: 22.6 % of floor firings by count, 23.7 % by volume, above the clamp; sensitivity ~1 K / 0.3 % price) | MMW surfaces at 1500–3500 K are far outside every Grimsel law. λ and c_p are extrapolated analytically with no clamp. |
| Jet heat correlation (flame only) | Goldstein: Re 61k–124k, H/D 6–12 | `A:22071-22075`, `A:24246-24248` | Höser `h_expr` applies `pow(0.0632/max(0.0632,s),1.3185)` with no upper clamp. At 120 s (axis ≈ 55 mm) s/D ≈ 15. Meier's map clamps s at 0.09 m (H/D = 12, `Claude_markdowns/2026-10-06.md:677`). |
| Domain, BCs | Neumann-0 everywhere; quarter symmetry on Höser and Meier (×4 volumes, SR 3) | RAW_KEYS | The ledger is invalid for non-Neumann-0 BCs. Floor must stay above the adiabatic base (RM:433-435). |
| Decomposition | byte-identical over np 1–4 and mgs 8–60 | §2.2 | Rebuilt binaries are not byte-reproducible (SR 9) |

---------------------------------------------------------------------------------------------------------------------
## (a) Items that would block a blind MMW prediction

1. **No energy-conserving α(T).** `BeamSource` uses the local α without the cumulative optical depth (`M:4423-4427`), so the defining MMW feedback (rising tan δ) cannot be run conservatively.
2. **No absorption, reflectivity or emissivity data on disk.** α = 300 m⁻¹ comes from main.tex's citation of Woskov 2009 (not on disk). R = 0.2, ε = 0.8 and h = 51 are tuned to one case whose reference says 0, 0.7/0.3 and 56.
3. **The only MMW validation is thermal.** It runs on a path (`material`) that cannot remove rock, and the removal path cannot use the Zhang high-T property law (`M:1082-1084`, Step 17b deferred).
4. **No high-temperature properties on the removal path.** Grimsel ρ, c_p and λ are measured to ~500 °C and extrapolated without a clamp. E, β and K_Ic are flat above 600 °C. Every deck defaults the melt and vapour scalars to off.
5. **The failure criterion is undefined for interior-peaked heating.** The P0 reconstruction returns "flat" when T_1 ≥ T_top (F7), and RM:360-363 says neither breaking rule handles volumetric heating.
6. **Today's removal rule (`face_list`) aborts with vapour, melt_aware or d_ref** and needs a flame patch (F6). An MMW deck would fall back to the older column path.
7. **The vapour front is untested against any data.** It ignores radiation losses and the uniform profile, and appears to spend beam power the field also receives (§1.3). No ablation-rate data are on disk (RM:328).
8. **Missing physics:** no melt evacuation (purge or shear), no vapour or plasma shielding, no bore waveguiding or wall absorption, and no softening s(T) (DIRECTION step B).
9. **The Zhang validation is stale and unguarded.** It was last run 2026-05-30 and is not in the 2639 tree gate, so nothing shows it still passes on today's binary.

## (b) Numerical settings without convergence evidence for the regime they are used in

1. MMW_ZHANG: dt 2 s (implicit, 26× the explicit limit, non-conservative latent linearisation) and 0.78 mm cells have no dt or mesh study.
2. HOSER_0 at measured heat: dt 8 ms at 1 mm gives an estimated 19–27 K/step against the ladder's ~6 K/step. There is no dt ladder at 7 MW/m².
3. HOSER_47 (G1, 120 s) and HOSER_0 (HG): no mesh twin. The only Höser mesh row is P2-9 c3, at 30 s on the 0 MPa M0 deck.
4. MEIER and MEIER_dref at 2 mm have no current mesh partner (D2g-1 and D2i failed on older decks; 1 mm is parked). ℓ ≈ 1.6–1.9 mm is under-resolved, and d_ref = 2 mm is one cell (+25 K).
5. The dt ladder (16 ms at 2 mm) was run on the D2d deck with the pinned closure and the pre-P0 criterion. "8 ms at 1 mm" is extrapolated from the rule, and no ladder was repeated on the P0b, P2-2 or face-list criterion.
6. KANT_onset and KANT_dref: dt 12.5 ms was never laddered. It was used unchanged at 2/1/0.5 mm.
7. HU_onset: 3.125 mm with two mesh points (40³ vs 32³); dt 0.8 s (0.56 of the stability limit) has no dt study.
8. Any beam run: the absorption depth vs cell size has no study under the current criterion. The only beam-on mesh study (S1) predates P0.

## (c) Contradictions with the context given

1. **"α(T) closure at `MMWSpalling.H` ~6166 (ROADMAP)"**: the substance holds (not done), but the line pointer is stale in both `RM:281` and `DIRECTION…:103`. The text now sits at `M:4423-4427`; line 6166 is in the cluster labeller.
2. **"κ/v ≈ 0.4–0.8 mm under one 1 mm cell (P2-7 note)"**: NOT FOUND in the repository as stated. The P2-7 note gives a top-cell centre 214–227 K below the face (`A:22240-22243`), which implies ℓ ≈ 0.85–0.95 mm by my inversion. Other repository values are 1.6–1.9 mm (RM:345) and ≈ 2 mm (`A:21927`).
3. **"P0's mesh-independent surface temperature ±0.5 %"**: ±0.5 % is the rev-5 **bar**. The measured spread was 0.161 / 0.116 %. It covers the **firing T** on the Kant cube only, not a rate: the floor rate later drifts 10.4 % (P2-9 c4).
4. **The 10.4 % floor drift** is confirmed, but the 0.5 mm leg also halved dt (8 → 4 ms, `unit/face_removal/test:55-59`). Mesh and dt are confounded at that point.
5. **"Meier at 2 mm (SR 12)"**: confirmed. SR 12 parks "AMR, the third mesh, the 1 mm run" (RM:89-90), and P2-2e says "Meier runs at 2 mm (SR 12)" (`A:21600-21601`).
6. **Zhang deck internals**: the deck calls 137 GHz the gyrotron frequency, but the reference names it as the **radiometer's** (summary :12). The deck comments state ε 0.7/0.3, h 56 and "forward Euler, 0.1 s" while the values are 0.8/0.8, 51 and implicit 2 s. The reference file contains no α value.
7. **Zhang "validated" curve**: the test plots ε·T with the tuned ε = 0.8 (the reference correction factor is 0.7), so the raw model peak is ≈ 3590 K, above T_vap_hi = 3503 K. Step 4b says it was tuned to Zhang's *simulation*, not the experiment (`A:246-249`).
8. **`microstructure.phase0.alpha_attenuation`** (2000 in Höser and Meier) looks like an MMW parameter but is never read (F4).
9. **"MMW has no validated case with the beam on"** (RM:309, :337) holds, but unvalidated beam-on spallation runs do exist (§1.5). P2-2e already made a first beam prediction (B100 never fires with d_ref) and is **not reviewed**.
