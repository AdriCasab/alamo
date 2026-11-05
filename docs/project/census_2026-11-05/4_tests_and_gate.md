# Census 2026-11-05, part 4 — tests and the tree gate

Read-only census of what every test, validation deck, gate deck and study exercises, so a deletion
packet knows which tests and gate inputs retire with each deletion. Nothing under `src/` or `tests/`
was edited. Binary audited: `bin/mmwspalling-3d-g++` sha256 `38c8623885c5b8b6…` (the same binary as
`docs/project/model_card/RUNNABLE.txt` and `studies/p2_16_span_rule/UNIT_p2_16/rc.txt`).

Paths below are relative to `tests/MMWSpalling/` unless they start with `src/` or `docs/`.

**Helper files** (all under `docs/project/census_2026-11-05/work/`):

| file | what it is |
|---|---|
| `classify.py` | static classifier: reads a deck (following `FILE =` includes, last value wins, argv last) and prints thermal path, removal path, lineages, closed switches still set, retired keys present |
| `decks_raw.txt` | `classify.py` run over every `input*` under `unit/`, `retired_units/`, `validation/`, `extensions/` |
| `sp_keys.txt` | every `surface_patch.*` key `SurfacePatch::Parse` reads (`src/Integrator/MMWSpalling.H:7897–8420`), used to spot strict-scope aborts |
| `runnable_all.py`, `RUNNABLE_ALL.txt` | **new run audit**: all 93 unit / validation / extension decks started on today's binary with `stop_time=1e-9`, `-np 4`, plot output to the scratchpad (same method as `docs/project/model_card/runnable.py`) |
| `gate_runs_p2_16.txt` | the 2,639 hashed paths of `witness/p2_16.md5` grouped by run (leg) directory |
| `studies_scan.py`, `STUDIES_SCAN.tsv` | static scan of every study and extension: recorded `argv*.sh` count, how many name a P2-10 or P2-16 retired key, unknown strict `surface_patch.*` keys, `face_list=1` count, key families named |

**Lineage tags used throughout**

| tag | meaning (and how it is detected) |
|---|---|
| LEFM | `spallation.model = sp_weibull` (Weibull flaws if `weibull.enabled = 1`, Grimsel law if `phase0.law = grimsel`) |
| DAMAGE | any of `damage.enabled = 1`, `damage.ic.*`, `cohesive.enabled = 1`, `spall.enabled = 1` with `spallation.model` absent or `damage_law` (**the default is `damage_law`**, `MMWSpalling.H:8693`, `:9140`), `spall.sample = column_max \| shallow_max`, `spall.{prescribed_ratio, alpha, t_spall, xi_s, C_h, damage_threshold}`, `conductivity.kappa_damage_alpha ≠ 0` |
| FLAME | `surface_patch.mode = convective_flame` (subtags: feet, jet = enthalpy, robin form) |
| PATCH-T | `surface_patch.mode = prescribed_T` (Hu, Rossi) |
| BEAM | `beam.P0 ≠ 0` (Beer–Lambert MMW source) |
| ZHANG | a non-microstructure thermal path: `material.type` present → `AdvanceMaterial` (explicit) or `AdvanceMaterialImplicit` (`solver = implicit`); neither `material.type` nor `microstructure.nphases` → `AdvanceConstAlpha` (dispatch `MMWSpalling.H:1457–1460`). Tagged ZHANG/MAT or ZHANG/CONST |
| MICRO | multi-phase, Voronoi (`microstructure.mode = voronoi`), `grain_expr`, alpha–beta, AMR regrid of phase fields |
| MECH | the mechanics operator: `el.type ≠ disable`, or `spall.stress_source` unset (default `mechanics`, `MMWSpalling.H:5416`) on a spall/onset deck |
| LTE | `spall.stress_source = local_thermoelastic` |
| DIAG / HARNESS | diagnostics only / source-lock or harness test |

**Removal path**: COLUMN = spall or vapor on with `surface_patch.face_list` unset/0 (`Removal.H`);
FACE = `surface_patch.face_list = 1` (`FaceRemoval.H`); none = no removal. The face list refuses, by
abort at parse (`MMWSpalling.H:972–1061`): vapor, bit, `foot_body_clearance`, `melt_aware`,
`surface_normal = 1`, a stress source other than `local_thermoelastic`, `amr.max_level > 0`, any
model other than `sp_weibull`, `spall.surface_T ≠ reconstructed`, `spall.removal_energy ≠ flake`,
`flake_length = flux` / `flake_fraction = weibull` / `size_volume = flux`, `flake_coherence_length > 0`,
`connected_cluster` with `A_crit > 0`; and it requires `mode = convective_flame` with `h_expr`,
`robin_form = surface`, the microstructure path and `surface.follow_mask = 1`. **So today no BEAM-only,
PATCH-T, pinned-robin or feet-with-clearance deck can be moved to the face list without the port.**

---

## 1. Kept unit tests

40 directories under `unit/`; verdicts on today's binary from `studies/p2_16_span_rule/UNIT_p2_16/rc.txt`
(5 pre-existing FAILs: `amr_microstructure_regrid`, `beam_void_closure` (by design, T3 is D2z-1's
documented debt, `witness/WHICH_TESTS.md`), `heterogeneous_kappa`, `sp_melt_skeleton`,
`surface_T_recon`). Every `input` in `unit/` parses and takes its first step today
(`work/RUNNABLE_ALL.txt`). **G** = in the 2,639-file tree gate (section 3).

| test | asserts (from the docstring) | runs (input + argv) | removal | lineage | disposition |
|---|---|---|---|---|---|
| `alpha_beta_transition` | quartz α–β eigenstrain appears in `thermal_free_strain` per phase | two-step: `input` (2-phase expression, `alpha_beta.enabled=1`, `el.type=static`) | none | MICRO, MECH | **KEEP** (deck names `conductivity.kappa_damage_alpha = 0.0`: edit if that key is put on an abort list) |
| `amr_microstructure_regrid` | AMR regrid keeps phase/grain/property fields sharp on every level; passive `D` strip sharp; `kappa_eff = kappa_phase·exp(−α_D·D)` | two-step: `input` (2-phase, `grain_expr`, `damage.ic` strip, `kappa_damage_alpha = ln 2`, `amr.max_level=1`) | none | MICRO + DAMAGE (passive `damage.ic`, κ-damage) | **RE-POINT** at DAMAGE: drop the `D`/`kappa_eff` clauses (`test:137–198`), keep phase/grain regrid repair. Already FAILs (rc 1) |
| `beam` | Beer–Lambert integral = `t·P0·(1−e^{−αL})` within 5 %; hot spot on axis | two-step: `input` (no material, no microstructure) | none | BEAM, ZHANG/CONST | **RE-POINT** at ZHANG merge: same check on a one-phase microstructure deck (`beam_void_closure` covers closure on that path but not the analytic integral) |
| `beam_void_closure` **G** | E1 beam/void gap closure: T1 counter, T2 plotfile audit, T3 φ-histogram (FAILs by design), T4 ≥ 3 layers removed | self-run `input` (one phase, `beam.P0=5816`, `spall.surface_normal=1`, `flake_coherence_length=0.004`) | COLUMN | BEAM, LEFM+Weibull, LTE; closed: `surface_normal`, `flake_coherence_length` | **RE-POINT** with the beam port (closure on the face list is still wanted); RETIRE-WITH column as written; also moves under `surface_normal` / `flake_coherence_length` deletion |
| `buckle_lock` | source lock B1–B7 on `spall.buckle_span` with mutation controls | source only (`MMWSpalling.H`, `FaceRemoval.H`, `Removal.H`) | FACE (source) | HARNESS, LEFM | **KEEP**; restate B2/B3 at column deletion (see 1.2) |
| `convective_patch` | Robin flame BC vs Carslaw & Jaeger half-space (transient, profile, bounds) | two-step: `input` (`material.type=constant`, explicit, `robin_form` default `cell`) | none | FLAME, ZHANG/MAT | **RE-POINT** at ZHANG merge (onto one phase; `robin_face` form legs already replay the micro-path Robin update); also moves if `robin_form = cell` is retired |
| `dp_yield` | Drucker–Prager yield, irreversible damage ODE, κ-damage | self-run `input` (4-slab, `damage.enabled=1`, `kappa_damage_alpha = ln 2`, mechanics) | none | DAMAGE, MECH | **RETIRE-WITH DAMAGE** |
| `equilibrium` | top-row T equals the radiation+convection balance within 5 % | two-step: `input` (const-α) | none | BEAM, ZHANG/CONST | **RE-POINT** at ZHANG merge (one-phase deck) |
| `face_buckle` | P2-15 buckling condition: F_0/F_S unmoved, N_S blocked, bk mesh-free, F_S_k0 aborts | self-run `unit/face_removal/input` + `face_list=1` (F_S_k0: `face_list=0` + span) | FACE | LEFM, LTE, FLAME | **KEEP**; leg F_S_k0 (abort when `face_list=0`) retires with the column path |
| `face_heat` | P2-11b `face_heat = path`: h ratio, path distance, flange blocking, march/ledger identities, refusals (vii) | self-run `validation/hoeser/input_lab` + G1's argv (`BASE`, `test:47–66`), seeded voids | FACE | FLAME, LEFM, LTE | **KEEP**. `BASE` names `spall.flake_coherence_length=0.0`, `spall.surface_normal=0`, `bit.enabled=0`: drop those tokens if the keys go on an abort list. Refusal legs R_pad (`pad_shield_open`), R_above (`jet_above_h`), R_robin (`robin_form=face`), R_feet/R_fs (feet) must be restated when those switches are deleted (then the abort is the strict-scope "unused key", not "P2-12's") or ported |
| `face_list` | P2-8 face census, per-face h, area rule; (iv) OFF0 == CAVU byte for byte | self-run `input` (from `make_input.py`) + `face_list=1` (CAV, STAIR), `=0` (OFF0), unset (CAVU) | FACE (OFF0/CAVU: COLUMN key-off) | FLAME, LEFM | **KEEP**; clause (iv) and legs OFF0/CAVU retire with the column path (the key itself) |
| `face_removal` | P2-9 one breaking rule per surface cell; `--offline` (a3 i) source: the column path and the line path call ONE `ReconFromX/ScanTAtZv/ScanKIcPa/ScanSigmaT/DepthScanH/RemapProfile` (`test:312–332`, reads `Removal.H`); floor T00 key 0 vs key 1 | self-run `input` (from `make_input.py`) + mesh/tilt/curved argv; `_k0` legs = column | FACE + COLUMN (k0 legs) | LEFM, LTE, FLAME | **KEEP** (core face-list test); at column deletion restate (a3 i) as "declared once, called by the line path" and drop the key-0 legs and the (a4) key-0-vs-key-1 comparison |
| `gb_cohesive` | bilinear cohesive envelope, irreversibility, ∫ = G_c; no-GB run all zero | self-run `input`, `input_no_gb` (`cohesive.enabled=1`, `grain_expr`, mechanics) | none | DAMAGE (cohesive), MICRO, MECH | **RETIRE-WITH DAMAGE** (`cohesive.*`) |
| `gradient_blind` | a uniform 820 K slab with no net flux MUST NOT fire (restraint reference = stored `T_init`) | self-run `input` (one leg, no argv) | COLUMN | LEFM, LTE, FLAME (`robin_form` default `cell`), closed: `flake_coherence_length=0.004` | **RE-POINT** onto the face list (the assertion is wanted). ⚠ Its `input` is also the runtime FIXTURE of `lateral_flux_lock` (`test:704`) and `retired_keys` (run half) |
| `grain_topology` | grain_id / phase / is_grain_boundary / is_phase_boundary from expressions | two-step: `input` (2-phase, `grain_expr`) | none | MICRO | **KEEP** |
| `heterogeneous_enthalpy` | `H = ρ_p·Cp_p·(T − T_ref)` per cell; thermal-mass ordering | two-step: `input` (2-phase laminate) | none | MICRO | **KEEP** |
| `heterogeneous_kappa` | series-resistance steady T; damaged half same profile; `kappa_eff = κ·e^{−α_D·D}`; `k_eff` | two-step: `input` (2-phase, `damage.ic` half, `kappa_damage_alpha = ln 2`) | none | MICRO + DAMAGE (passive κ-damage) | **RE-POINT** at DAMAGE: keep checks 1 and 4, drop 2 and 3. Already FAILs (rc 1) |
| `jet_d2c` | D2c march options (`jet_decay_diameter`, `jet_free_surface`, `jet_bin_update`, `jet_negative_flux`, `foot_body_clearance`), aborts A1–A19, `rec_fixed` legs KT/K15/A17–A19 | self-run `../robin_feet/input_feet`, `../jet_enthalpy/input_strip`, `input_quarter`, `../robin_jet/input_jet` | COLUMN | FLAME (feet, enthalpy), LEFM, LTE; closed: `melt_aware=1` (inputs), `jet_T_ent_mode = rec_fixed` (KT…A19) | **RE-POINT** with the feet port; KT/K15/A17–A19 **RETIRE-WITH `rec_fixed`**. ⚠ Run by the gate scripts but **not hashed** (section 3) |
| `jet_enthalpy` **G** | D2a2 enthalpy march: N1 == N0, discrete march to 1e-9, power identities, ROP mesh check | self-run `input_strip` (robin `cell`), `input_quarter` (pinned), `../robin_jet/input_jet` (N0/N1) | COLUMN | FLAME (enthalpy), LEFM, LTE; closed: `melt_aware=1`, `robin_form=cell` (strip) | **RE-POINT** on the face list (the march is kept flame capability) |
| `latent_heat` | the enthalpy table carries the latent heat under `law = constant` and `grimsel` | self-run `input` + `phase0.law=grimsel` | none | LEFM (Grimsel law), thermal | **KEEP** |
| `lateral_flux_lock` | source lock: stress-free reference, property route, size factor, conductivity, `LawSigma`, face-list one-builder/one-heat-rule/one-breaking-rule, P2-10 "nothing deleted survives"; runtime: retired keys unconsumed / group-A keys abort | source (`MMWSpalling.H`, `Removal.H`, `FaceRemoval.H`); runtime legs on `../gradient_blind/input` | COLUMN + FACE (source) | HARNESS | **KEEP**; heavy restatement at column deletion (see 1.2) |
| `pad_clear_ledger` | P2-2f pad-clearance removal books its own ledger term; N leg has no `foot_mech_E` column | self-run `../robin_feet/input_feet` + `foot_body_clearance=1/0`, jet_d2c's DEAD `h_expr` | COLUMN | FLAME (feet, clearance), LEFM | **RE-POINT** with the feet + clearance port (the face list refuses `foot_body_clearance`, `MMWSpalling.H:1033`) |
| `pin_clock` | the pin/idle-clock product: `rho_phase`, `Cp_phase` exact under `constant`, non-zero, differ under `grimsel` | self-run `input` (spall off) + `phase0.law=grimsel` | none | LEFM law, FLAME (pinned clock) | **KEEP** while `robin_form = pinned` exists; RETIRE-WITH pinned if the pinned form goes with the column path (the face list requires `surface`) |
| `regime_low_high_power` | Step 11 vapor vs spall regimes, closed-form `RoP_vap` | self-run `input_low_power`, `input_high_power` (damage-law spall, `spall.alpha/t_spall/xi_s/C_h/damage_threshold`, vapor, beam) | COLUMN | DAMAGE, VAPOR, BEAM, MECH | **RETIRE-WITH DAMAGE**; the vapor half is re-pointed by the vapour port (`sp_v_n_regime` is the LEFM vapor test) |
| `removal_increment` | D2z-1 whole-cell removal unconditional (W), prescribed-T keeps the scan (PRESC), SN carries `sn_recession_superseded` | self-run `input` (pinned flame), `input_presc` (PATCH-T), `input + spall.surface_normal=1` | COLUMN | LEFM, LTE, FLAME, PATCH-T; closed: `surface_normal` (SN) | **RETIRE-WITH column** (whole-cell column recession); SN leg **RETIRE-WITH `surface_normal`** |
| `retired_keys` | P2-16: constraint depth deleted (R1–R4 source, mutations, `--pre`); K1/K2 abort, K0 runs | source walk of `src/`; run half on `../gradient_blind/input`, one step | COLUMN (fixture) | HARNESS | **KEEP and EXTEND**: one R-clause + one K-leg per new retired list; re-point the fixture at column deletion (see 1.2) |
| `robin_face` **G** | (a) `robin_form = cell \| face` replayed in Python to 1e-10, Bi ratio, ledger; (b) face form + removal + `follow_mask` + ledger; `follow_mask=0` control | self-run `input_form` (4 legs: `robin_form=cell/face`, losses off/on), `input_ledger` (face) + `ledger_nofollow` (`surface.follow_mask=0`) | none (form) / COLUMN (ledger) | FLAME, LEFM; closed: `robin_form = cell, face` | form legs **KEEP** while `cell`/`face` exist; ledger legs **RE-POINT** at column deletion; whole test **RETIRE-WITH** `robin_form = cell/face` (re-point: a `surface`-form replica) |
| `robin_feet` **G** | D2b burner-on-feet descent, pad rule, liftoff cap, stall, rank identity, pinned closed form, exhaust recirculation | self-run `input_feet`, `input_column`, `../jet_enthalpy/input_quarter` | COLUMN | FLAME (feet, enthalpy, pinned), LEFM, LTE; closed: `melt_aware=1` | **RE-POINT** with the feet port |
| `robin_jet` **G** | D2a spatial `h_expr`/`T_flame_expr`, nozzle path, pinned counters, strip ROP, rank identity | self-run `input_jet`, `../robin_pinned/input_pinned` | COLUMN | FLAME (pinned), LEFM, LTE; closed: `melt_aware=1` | **RE-POINT** on the face list (expressions are the face list's `h_expr`) |
| `robin_pinned` **G** | C1 pinned closure: P0 == F, closed form, both branches, rank identity | self-run `input_pinned` + F (`robin_form=face`), P0/P/P4 (idle cycles), M | COLUMN | FLAME (pinned), LEFM, LTE; closed: `melt_aware=1`, `robin_form=face` (F) | **RETIRE-WITH** pinned/column (the face list requires `surface`); RE-POINT only if the pinned form is ported |
| `scalar_flaw` **G** | A1 scalar-flaw path; removal event log complete and exact; S == D; S4 rank identity | self-run `input` (beam, `weibull.enabled=0`) + D (degenerate Weibull), S4 | COLUMN | BEAM, LEFM, LTE; closed: `melt_aware=1` | **RE-POINT** with the beam port (event-log completeness and rank identity are wanted) |
| `seed_void` | D2q-2b interior void: beam invariant, path in solid cells, k_top scans, depth-scan clamp | self-run `input` (`spall.seed_void`, beam, `surface_normal=1`, `flake_coherence_length`) | COLUMN | BEAM, LEFM, LTE, DIAG | **RETIRE-WITH column** (its four call sites are column-top scans); re-point the beam-path claim (S2) onto the beam port. ⚠ `spall.seed_void` itself stays: `face_list` and `face_heat` use it as a fixture tool |
| `sp_lefm` | K_I Tada integral + Nasseri K_Ic closed form; Weibull flaw statistics; `flaw_a` bit-identical across regrids | self-run `input_kic_table`, `input_weibull`, `input_weibull_regrid` (prescribed σ; `el.type=disable`; `amr.max_level=1` on regrid) | none | LEFM, MICRO (Voronoi single phase) | **KEEP** (LEFM core) |
| `sp_melt_skeleton` | P1 melt-aware A/B on an imposed molten cap: ON removes 0, OFF drills | self-run `input` + `spall.melt_aware=1/0` | COLUMN | LEFM, MICRO (Voronoi n1); closed: `melt_aware` | **UNCLEAR**: RE-POINT if the melt port carries `melt_aware` onto the face list, otherwise **RETIRE-WITH `melt_aware`** (C.7.3 lists it closed). Already FAILs (rc 1) |
| `sp_v_n_regime` | Step 17 per-cell `regime_field`, local-Q `RoP_vap` radial variation | self-run `input` (Voronoi n1, beam, vapor) | COLUMN | LEFM, VAPOR, BEAM, MICRO | **RE-POINT** with the vapour port |
| `spall_event` **G** | Step 10 phi shift, detached-cell state, surface mask migration; no-event companion | self-run `input`, `input_no_event` (`spall.enabled=1`, no `spallation.model` ⇒ damage_law, `damage.ic`, `spall.alpha/t_spall/xi_s/C_h/damage_threshold`, `el.type=static`) | COLUMN | DAMAGE, MECH | **RETIRE-WITH DAMAGE** |
| `stefan` | enthalpy Stefan front within 5 % | two-step: `input` (const-α, `phase.*`) | none | ZHANG/CONST | **RE-POINT** at ZHANG merge (one-phase deck with `phase0.T_m/L_m`; `latent_heat` checks the table, not the front) |
| `surface_T_recon` | P0 reconstruction invertible; header gating; `diagnostic` logs without applying; partial recession under `cell` | self-run `input` A + B (`removal_energy=cell`), C (`surface_T=diagnostic`), D (`surface_T=cell_average`) | COLUMN | LEFM, LTE, FLAME (`robin_form` default `cell`); closed: `flake_coherence_length` | **RETIRE-WITH column** (`cell_average`, `diagnostic` and `removal_energy = cell` are refused by the face list); RE-POINT check (1) (invertibility) onto the face list. Already FAILs (rc 1) |
| `thermal_stress` | per-phase moduli, nodal free strain, stress magnitude ordering | two-step: `input` (2-phase, `el.type=static`) | none | MICRO, MECH | **KEEP** |
| `voronoi` | Voronoi phase fractions, `is_gb`, phase ids | two-step: `input` (3 phases, 200 grains) | none | MICRO | **KEEP** |

**Contrast — `retired_units/`** (P2-16 moved them out; never run by any harness):

| test | asserted | why retired |
|---|---|---|
| `retired_units/constraint_ref_lock` | P2-2e/P2-14a source lock: every stress-reference read uses the constraint reference under `spallation.sp.constraint_depth > 0` (Removal.H depth scan, two `Sp_field` branches, the face list's line scan) | its subject was deleted in P2-16; the key now aborts (`MMWSpalling.H:128–133`) |
| `retired_units/face_cref` | P2-14a constraint depth on the face list vs the column path (CROSS/TILT/NEG legs) | same; the run half would abort at parse |

P2-10's six deleted tests (`side_face_flux`, `side_face_area`, `side_face_step_height`, `side_staircase`,
`lateral_drive`, `wall_march`) are archived under `studies/p2_10_retire/retired_units/` (they show as
`D` in `git status`). This is the precedent for every retirement below: move, do not delete.

### 1.1 What `witness/WHICH_TESTS.md` says is the kept set

`studies/tree_sweep_0922/witness/WHICH_TESTS.md` is a site → test table for the per-build check, not
a full list. Rows: the enthalpy march (`jet_enthalpy`, `jet_d2c`), `BuildFlameColumns` (`jet_enthalpy`,
`robin_jet`), the feet rule (`robin_feet`), the top-face Robin (`robin_face`, `robin_pinned`), **removal
(`Removal.H`): `spall_event`, `removal_increment`, `scalar_flaw`**, the beam (`beam_void_closure`),
mesh/AMR (`s1_surface_resolution --smoke`), the retired-key lists (`retired_keys`, `lateral_flux_lock`),
the buckling condition (`buckle_lock`, `lateral_flux_lock`, `face_buckle`), the ledger (`robin_face`,
`scalar_flaw`). ⚠ Every row except the retired-key and buckling rows names a COLUMN-path test, and its
"removal" row names a DAMAGE test (`spall_event`) — the table must be rewritten at column deletion.
There is no row for the face list's breaking rule (`face_removal`) or builder (`face_list`).

### 1.2 The source-lock tests: what each deletion does to them

| lock | today | (a) DAMAGE | (b) COLUMN path | (c) ZHANG path | (d) closed switches |
|---|---|---|---|---|---|
| `lateral_flux_lock` | 104 `need()` sites; **48 read `Removal.H`** (`rem`, `both`, `allc`, `code_r`), and `open(REM)` at `test:89` is unconditional | no clause names damage; ADD a "no damage/cohesive name survives in `src/` outside the retired list" clause (P2-10's `DELETED` regex pattern, `test:931–946`) and unconsumed/abort legs for the damage keys | **raises at `test:89`** (FileNotFoundError). Delete clauses whose only site is in `Removal.H` (`:187` `(T_z − T0_z)`, `:189/:213` `T0_at_didx`, `:222`, `:257` `phase_beta_loc`, `:282–293` `Grimsel` reads, `:357` `ScanKIcPa` declared in `Removal.H`, `:398` thickness-scan flaw read, `:414` `if (sp_size_flux)`, `:435` `rc_ellq`, `:607`); restate the counts (`:311` "THREE stress expressions", `:316` "THREE E·β products", `:354` "both toughness sites", `:401` "all N flaw reads", `:857–904` ndecl/ncall over `code_r`); re-point the runtime fixture (`:704` `gradient_blind/input`) to a face-list deck; ADD "no column-path name survives" | ADD "`AdvanceMaterial`, `AdvanceMaterialImplicit`, `AdvanceConstAlpha`, `SolverKind::Implicit` gone; `material.type`, `solver` abort" | ADD one name clause per deleted switch (`spall_surface_normal`, `sn_recession_superseded`, `melt_aware`, `spall_flake_coh_len`, `pad_shield_open`, `jet_above_h`, `rec_fixed`, …) |
| `buckle_lock` | B1–B7 + mutations; `read()` opens `Removal.H` | none | B2's "`Removal.H` never reads `sp_buckle_span`" half and every mutation targeting `r` go; `read()` must stop opening `Removal.H`. B3 (abort under `face_list = 0`) retires if the key itself is retired (the face list is then the only path) | none | none |
| `retired_keys` | R1–R4 source (walks all of `src/`), `--pre` substitutes `Removal.H` (`SUBST`); run half K0/K1/K2 on `gradient_blind/input` | ADD R/K pair for the damage list (key strings appear only in the new retired list; one damage key aborts naming the step) | `--pre` `SUBST` entry for `Removal.H` must go; the K-legs' fixture must become a face-list deck (`gradient_blind/input` is column) | ADD R/K pair (`material.type`, `solver`, `implicit.*`, `phase.*` if the const-α keys go) | ADD R/K pair per retired switch |
| `face_removal --offline` (a3 i) | asserts the six scan/remap helpers are called by the column path AND the line path | none | restate: declared once, called by the line path | none | none |
| `face_heat` (vii) | refusals name `pad_shield_open`, `jet_free_surface`, feet, `jet_above_h`, `robin_form`, `jet_closure`, `barrier_expr`, `face_list = 0` | none | R_nolist (`face_list=0`) retires with the key | none | R_pad / R_above / R_robin become strict-scope "unused key" aborts when those keys are deleted: restate the expected message |

⚠ **The default-value trap** (standing rule 1: "only a run can retire a key: ALAMO does not abort on
an unused top-level key"). P2-10's retired list aborts on `pp.contains(k)` (`MMWSpalling.H:119–122`), i.e.
**even at the old default**. Kept decks name the to-be-deleted keys at their defaults far more often
than at live values (counted over the 95 kept deck files: unit, validation, shared includes, extensions,
`s1_surface_resolution/input`):

| key | decks naming it | values |
|---|---|---|
| `damage.enabled` (+ `damage.*`) | 59 | `=0` ×46, `=1` ×11 |
| `conductivity.kappa_damage_alpha` | 78 | `=0.0` ×75, `=ln 2` ×3 |
| `conductivity.h_gb0` | 78 | `=0.0` ×78 (grain-boundary conductance — MICRO, probably kept) |
| `solver` | 76 | `explicit` ×75, `implicit` ×1 |
| `spall.sample` | 43 | `top_cell` ×39, `shallow_max` ×4 |
| `spall.surface_normal` | 25 | `=0` ×21, `=1` ×4 |
| `spall.flake_coherence_length` | 24 | `0.004` ×12, `0.0` ×7, `0.002` ×5 |
| `spall.melt_aware` / `melt.lambda_crit` | 11 / 10 | `=1` ×9, `=0` ×2 / `0.4` ×10 |
| `material.type` | 11 | `constant` ×10, `zhang` ×1 |

So an abort-on-contains list for damage, `solver` or `kappa_damage_alpha` breaks most of the kept
decks until each is edited. P2-10 handled the same problem by editing the kept `make_input.py`
generators (`studies/p2_10_retire/TEST_EDITS.txt`); the alternative is to leave the parse of the
default value in place and abort only on a non-default value.

---

## 2. Validation decks

All decks run on today's binary except the two templates (`work/RUNNABLE_ALL.txt`; the five also in
`docs/project/model_card/RUNNABLE.txt` agree). None names a P2-10 or P2-16 retired key; the retired keys
live only in study argv (section 4).

| deck | lineage | removal | runs today | disposition |
|---|---|---|---|---|
| `hu/hu_conduction/input_{granite2,sandstone2}` | ZHANG/MAT (`material.type = constant`, explicit), PATCH-T | none | RUNS | **MERGE** with ZHANG: re-express as a one-phase microstructure deck (the micro path already carries PATCH-T: `hu_spall_onset`, `rossi`). Not damage |
| `hu/hu_thermoelastic/input_{granite2,sandstone2}` | ZHANG/MAT, PATCH-T, MECH (`el`) | none | RUNS | **MERGE** (as above); mechanics kept (decision 5) |
| `hu/hu_thermoelastic_amr/input_{granite2,sandstone2}` | ZHANG/MAT, PATCH-T, MECH, `amr.max_level = 1` | none | RUNS | **MERGE**; AMR on the micro path exists (`unit/amr_microstructure_regrid`) |
| `hu/hu_breakage_index/input_{granite2_37,sandstone2_37,sandstone2_89}` | ZHANG/MAT, PATCH-T, MECH, AMR; "intentionally leaves the Step 8 Drucker-Prager damage model disabled" (`test:6–7`) | none | RUNS | **MERGE** (a stress indicator, not damage) |
| `hu/hu_spall_onset/input_{granite2,sandstone2}` | DAMAGE (`damage.enabled = 1`, D-threshold onset proxy), PATCH-T, MECH, one phase | none | RUNS | **RETIRE-WITH DAMAGE** unless re-expressed as an LEFM onset (Sp_field under PATCH-T + mechanics or LTE, removal off — the Kant onset decks prove the mechanism) |
| `hu/hu_spall_onset/diag/input_{granite2_base,sandstone2_base,sandstone2_clamp,sandstone2_clampall}` | DAMAGE, PATCH-T, MECH, AMR | none | RUNS | **RETIRE-WITH DAMAGE** (diagnostic records) |
| `hu/hu_end_to_end/input_{granite2,sandstone2}` | DAMAGE (damage_law spall, `spall.sample = shallow_max`, `xi_s`, `C_h`, `damage_threshold`), PATCH-T, MECH | COLUMN | RUNS | **RETIRE-WITH DAMAGE** |
| `zhang/zhang_oglesby/input` | ZHANG/MAT-IMPL (`material.type = zhang`, `solver = implicit`), BEAM | none | RUNS | **MERGE**: the Zhang ρ·c_p(T), κ(T) as a phase law on the micro path; the test reads only `losses.*` and the plotfile, so it survives a deck rewrite |
| `kant/sp_kant_onset/input_{verification,verification_hfl,model}` | LEFM, LTE, FLAME (`robin_form` default `cell`), one phase | none (Sp_field onset only) | RUNS | **KEEP** (no removal path involved); moves only if `robin_form = cell` is retired |
| `kant/sp_kant_onset/input_weibull`, `input_weibull_cluster` | LEFM+Weibull, LTE, FLAME (`cell`), MICRO (Voronoi single phase); `_cluster`: `weibull.detach_mode = connected_cluster` | none | RUNS | **KEEP** (Voronoi kept as future study); `_cluster` retires with `connected_cluster` if that goes |
| `kant/sp_kant_pressure_sweep/input` | LEFM, LTE, FLAME (`cell`), `confining.p` by argv | none | RUNS | **KEEP** |
| `kant/sp_kant_onset_meiercfg/input_{model,weibull}` | LEFM, LTE, FLAME (`cell`; legs S150/S500 `surface`), closed `flake_coherence_length = 0.004` | COLUMN | RUNS | **RE-POINT** onto the face list (`kant/grimsel_calibration` is already its face-list successor); RETIRE-WITH column as written. Base of `studies/p2_2b_props`, `p2_2d_kic`, `p0_surface_T` |
| `kant/grimsel_calibration/input` | LEFM+Weibull + Grimsel law, LTE, FLAME (`h_expr`, `robin_form = surface` via `shared/model_form_flame.inputs`) | FACE | RUNS | **KEEP** (the production calibration deck) |
| `kant/grimsel_calibration/pre/input_lab.pre_p2_13` | LEFM, flake coherence | COLUMN | FAILS (template, no `amr.n_cell`) | record |
| `hoeser/input_lab` | LEFM+Weibull + Grimsel law, LTE, FLAME via the two shared includes | FACE | FAILS (template; runs with G1's argv, `RUNNABLE.txt` HOSER_47 G1 RUNS) | **KEEP** |
| `shared/rock_grimsel.inputs`, `shared/model_form_flame.inputs` | rock + model form (`spall.surface_T = reconstructed`, `removal_energy = flake`, `robin_form = surface`, `face_list = 1`) | FACE | (includes) | **KEEP** |
| `meier/sp_meier_pilot/input_feet`, `input_feet_d2c` | LEFM+Weibull, LTE, FLAME (feet, enthalpy, `robin_form = pinned`) | COLUMN | RUNS | **needs the feet port** (and `robin_form = surface` + `h_expr`); RETIRE-WITH column otherwise. `test_feet*` are check-only scorers of `studies/d2b_feet_rop` output |
| `meier/sp_meier_pilot/input_energy`, `input_drilling`, `input_drilling_25mm` | LEFM+Weibull, LTE, BEAM (a beam stand-in for the flame), `flake_coherence_length`; drilling: BIT | COLUMN | RUNS | needs the beam port; the drilling decks also need `bit.enabled` (regime 3, refused by the face list and **not on the port list**) ⇒ **UNCLEAR**. Base of most D2* and P1* studies |
| `meier/sp_meier_pilot/input_2d_dev` (gate deck `dev2d`), `input_2d_dev_melt` | LEFM+Weibull, LTE, BEAM; `surface_normal = 1`, `flake_coherence_length = 0.004`; `_melt`: `melt_aware = 1` | COLUMN | RUNS | needs the beam port; drops from the gate at column, `surface_normal` or `flake_coherence_length` deletion |
| `rossi/sp_rossi_damage_profile/input` | LEFM+Weibull (despite the name), LTE, PATCH-T, MICRO (Voronoi single phase); scores `h_col_events.csv` | COLUMN | RUNS | **UNCLEAR → RE-POINT or RETIRE-WITH column**: the face list refuses PATCH-T and does not write `h_col_events.csv`; keep only if PATCH-T is ported |
| `p0_steady_column/input` | LEFM, LTE, FLAME (`cell` base; legs `surface/cell/face/pinned`), `flake_coherence_length = 0.004`; its harness still passes `spall.driving_stress`/`driving_ref` (deleted P2-2a, silently unconsumed, `test:121`) | COLUMN | RUNS | **RE-POINT** to the face list's flat floor (`unit/face_removal` T00 legs already score the same rate) or RETIRE-WITH column |
| `spall_thickness/REFERENCES.md`, `grimsel_properties/grimsel_granite_properties_reference.md` | reference notes only | — | — | KEEP (no deck) |

---

## 3. The tree gate (`studies/tree_sweep_0922/witness/`)

### 3.1 What it hashes, exactly

- **Runner.** `guarded_post_all_par.sh <label>` fingerprints the binary (mtime, size, sha256), runs
  `post_all_par.sh` (three static lanes), re-fingerprints, exits 2 if the binary moved, and memoises a
  pass on the binary hash in `witness/.gate_cache/<sha256>`. `post_all.sh` + `wit.sh` are the serial
  reference ("if the parallel list ever disagrees with a serial one, the serial is right",
  `REFERENCE.md`).
- **What runs** (`post_all_par.sh`, `wit.sh`): `unit/jet_enthalpy`, `robin_pinned`, `robin_face`,
  `robin_feet`, `scalar_flaw`, `spall_event`, `jet_d2c`, `robin_jet`, `beam_void_closure` (each its own
  `test`, venv python), `mpirun -np 4 bin/mmwspalling-3d-g++ validation/meier/sp_meier_pilot/input_2d_dev`
  (leg `dev2d`), and `studies/s1_surface_resolution/run_sweep.py --smoke --force`. rc is reported,
  never fatal.
- **What is hashed** (`hash.sh`): under the directories `unit/{robin_face, robin_pinned, robin_jet,
  jet_enthalpy, beam_void_closure, scalar_flaw}/output`, `unit/spall_event/{output, output_no_event}`,
  `unit/robin_feet/output`, `validation/meier/sp_meier_pilot/output/dev2d`,
  `studies/s1_surface_resolution/output_smoke`, every file matching `*Level_0/Cell_D_*` (cell **and**
  node plotfiles), `*_removal_events.csv`, `*_h_col_events.csv`, `*_clusters.csv`, `*_jet_profile.csv`,
  `thermo.dat`; `*.old*` pruned; `md5 -r`, sorted unique by path. Not hashed: plotfile `Header`/`Cell_H`,
  logs, AMR levels > 0, **all of `unit/jet_d2c/output`** (it runs, but its directory is not in
  `hash.sh`'s list — its only contribution is its rc, which is not fatal), and **no face-list output at
  all**.
- **Composition of `p2_16.md5`** (2,639 paths): 2,466 `Cell_D`, 63 `thermo.dat`, 47
  `removal_events.csv`, 56 `h_col_events.csv`, 7 `jet_profile.csv`, 0 `clusters.csv`.
- **Reference history.** `witness/REFERENCE.md` still declares `d2z1.md5` the reference (stale, as
  `ARCHIVE_DONE.md:19260` notes). In practice each packet diffs against the previous packet's list. The
  list last moved at **P2-2a** (`p2_1.md5` → `p2_2a_r7.md5`: 26 `jet_enthalpy` strip `Cell_D` files,
  legs B05–B4, E, Eon, G, Gbug); every list from `p2_2a_r7` to `p2_16` is identical. The two "controls"
  every packet quotes are fixed old lists that prove the comparator is live: `d2aa.md5` differs in
  exactly those **26** files, `d2t.md5` in **863**. The only documented re-pin procedure is D2z-1's
  (`REFERENCE.md` "How to re-pin": new list → diff and classify every differing path → declare it in
  `REFERENCE.md` → keep the old list); its audit is `studies/d2z1_delete/tree_repin_audit.txt`.
- **How P2-10 handled it** (`ARCHIVE_DONE.md` heading "P2-10", `studies/p2_10_retire/`): P2-10 did
  **not** drop any gate deck. Every key it deleted was off in all 2,639 hashed runs, so the gate was
  run as-is on both deletion binaries — **0 of 2639** vs `p2_9_bB`, controls 26 / 863
  (`GATE_tree_c1.txt`, `GATE_tree_c2.txt`). Because "the 2639 sweep is blind to the deletions" (review
  DR-2, `ARCHIVE_DONE.md:22789`), it added a **separate witness set W**, hashed on the pre-deletion
  binary before the first source edit (`W_bB.md5`): `unit/face_removal` T00_1_k0, T00_1_k1, TL35_1,
  TL45_1, SPH_1, BORE_1; `unit/face_list` CAV, STAIR, OFF0, CAVU; and M0T (P2-9's M0 Hoeser leg minus
  its 13 retired tokens). It compared W column by column, excluding only the thermo columns it deleted,
  each declared in `DELETION_MAP.txt` and shown identically 0 in the reference. It also ran every kept
  unit test before and after (`UNIT_bB` vs `UNIT_c2`, rc and PASS/FAIL lines), archived the six retired
  units in `retired_units/`, itemised every lock edit in `TEST_EDITS.txt`, and ran one abort leg per
  retired key (`RETIRED.txt`). **That is the template, but it only works when the deleted code is off
  in every gate deck — which is false for every deletion below except ZHANG.**

### 3.2 The gate decks (66 legs in 10 output directories)

Files per leg are from `work/gate_runs_p2_16.txt`. "melt" = `spall.melt_aware = 1` in the base input;
"cell/face/pinned" = `surface_patch.robin_form`.

| test dir (files) | legs (files) | deck(s) | removal | lineage | closed switches set |
|---|---|---|---|---|---|
| `unit/robin_face/output` (798) | form_cell_off, form_cell_on, form_face_off, form_face_on (169 each = 676) | `input_form` + `robin_form` + losses argv | **none** | FLAME | `robin_form = cell / face` |
| | ledger (89), ledger_nofollow (33) | `input_ledger` (face form) | COLUMN | FLAME, LEFM+Weibull, LTE | `robin_form = face` |
| `unit/robin_pinned/output` (323) | F (45, `robin_form=face`), P0 (45), P (45), P4 (171), M (17) | `input_pinned` | COLUMN | FLAME (pinned), LEFM, LTE | melt; F: face |
| `unit/robin_jet/output` (448) | A (45), B (85), C (45, on `robin_pinned/input_pinned`), I (85), S (45), S3 (129), X (5), X2 (9) | `input_jet` | COLUMN | FLAME (pinned, `h_expr`), LEFM, LTE | melt |
| `unit/jet_enthalpy/output` (288) | strip B05, B1, B1e, B1n, B2, B4 (10 each), B3 (22), E, Eon, G, Gbug (9 each) = 118 | `input_strip` | COLUMN | FLAME (enthalpy), LEFM, LTE | melt; `cell` |
| | N0, N1 (25 each) | `../robin_jet/input_jet` (N1: `jet_closure = none`) | COLUMN | FLAME (pinned) | melt |
| | Q (17), Qh (17), Q3 (45), QN (9), QE (9), QF (23) = 120 | `input_quarter` | COLUMN | FLAME (enthalpy, pinned) | melt |
| `unit/robin_feet/output` (138) | D0, D1 (7 each, `jet_enthalpy/input_quarter`), F (13), F3 (33), FM (9), FA (9), FC (13), SM (13), SP (13), R1, R2 (5 each, runtime aborts), H (11, `input_column`) | `input_feet`, `input_column` | COLUMN | FLAME (feet, enthalpy, pinned), LEFM, LTE | melt |
| `unit/scalar_flaw/output` (261) | S (45), D (45), S4 (171) | `input` | COLUMN | BEAM, LEFM, LTE | melt |
| `unit/beam_void_closure/output` (169) | one run (21 cell + 21 node plotfiles × 4 + thermo) | `input` | COLUMN | BEAM, LEFM+Weibull, LTE | `surface_normal = 1` (adds thermo column `sn_recession_superseded`), `flake_coherence_length = 0.004` |
| `unit/spall_event/output`, `output_no_event` (4 + 4) | one run each | `input`, `input_no_event` | COLUMN | **DAMAGE**, MECH (`el.type = static`) | — |
| `validation/meier/sp_meier_pilot/output/dev2d` (64) | one run | `input_2d_dev` | COLUMN | BEAM, LEFM+Weibull, LTE | `surface_normal = 1`, `flake_coherence_length = 0.004` |
| `studies/s1_surface_resolution/output_smoke` (142) | F1, F1s, M-100, M-2000 (beam, 14 each); R-cell (`cell`), R-face, R-face-1400 (`face`) (14 each); rankcheck_np1 (10), rankcheck_np4 (34) (F1s) | `studies/s1_surface_resolution/input` + `run_sweep.py` overrides | COLUMN | BEAM or FLAME, LEFM+Weibull, LTE | melt (input sets `melt_aware = 1`); R-legs cell/face |

Every gate deck is on the microstructure path (none is ZHANG); none uses the face list; none uses the
current rate model (`spall.surface_T = reconstructed`, `removal_energy = flake`).

### 3.3 Which gate decks must leave the gate under each deletion

| deletion | gate legs whose output must change (drop or re-pin) | files | note |
|---|---|---|---|
| (a) DAMAGE lineage | `spall_event` output + output_no_event | **8** | ⚠⚠ plus a **whole-tree hazard**: `D`, `kappa_eff`, `k_eff` are registered **plotted** on every microstructure run (`MMWSpalling.H:502–504`) and `regime` on every removal run (`:534`, its comment calls it "damage_law back-compat"), so they are components of every hashed `Cell_D`. Deleting them from the plotfile changes all **2,466** `Cell_D` paths with no physics change (CODE_CENSUS B1: registration order is a hidden contract). Either keep them registered (constant) in the deletion step and remove them in a separate declared re-pin, or re-pin the whole tree with an audit showing only the component count moved. `kappa_damage_alpha = 0.0` in every other gate deck, so the κ-damage factor is `exp(0) = 1` exactly and numerically neutral. Also: `spallation.model` defaults to `damage_law`; the default must become `sp_weibull` or required — every gate deck except `spall_event` sets it explicitly |
| (b) COLUMN path (`Removal.H`) | every leg with removal on: everything except `robin_face` form_* | **1,963** of 2,639 (only the 676 `robin_face` form files survive, and those are `cell`/`face` controls) | **The 2,639 gate cannot pass and has no successor inside itself.** No gate deck runs the face list today, and none can until the port: the face list needs `robin_form = surface`, an `h_expr` flame, `surface_T = reconstructed`, `removal_energy = flake`, and refuses beam-only, pinned, melt, `surface_normal`, `flake_coherence_length > 0`, clearance. The h_col event CSVs (56 files) have no face-list writer |
| (c) ZHANG separate path | none | **0** | the only deletion the P2-10 template covers unchanged: run the gate as-is (0 of 2639) and run the re-pointed unit tests and Zhang/Hu decks as the witness |
| (d1) `spall.surface_normal` | `beam_void_closure`, `dev2d` | 233 | thermo header loses `sn_recession_superseded` (certain) and the criterion depth `/cos θ` (`MMWSpalling.H` criterion update) changes Cell_D |
| (d2) `spall.flake_coherence_length` (> 0) | `beam_void_closure`, `dev2d` | 233 (same legs) | face list already refuses it; naturally deleted with the column path |
| (d3) `spall.melt_aware` | `robin_pinned`, `robin_jet`, `jet_enthalpy`, `robin_feet`, `scalar_flaw`, `s1` smoke | 1,600 (if any leg melts) | registers no output; `f(Λ_L) = max(0, 1 − Λ_L/λ_c)` is exactly 1 where `Λ_L = 0`, so deletion is hash-neutral **iff no firing cell in any leg ever has `Λ_L > 0`** — must be measured (0 of 2639), not assumed; M-2000/M-100 at 15 MW/m² are the likely movers. ⚠ Contradiction to resolve first: C.7.3 lists it closed, but the plan ports "melt" to the face list |
| (d4) `robin_form = cell / face` (the C.7.3 "controls") | `robin_face` (all 798), `robin_pinned` F (45), `jet_enthalpy` strip legs (118), s1 R-cell / R-face / R-face-1400 (42) | **1,003** | ⚠ `cell` is the **parse default** (`MMWSpalling.H:7962`): retiring it needs a new default, which also moves every Kant onset deck, `p0_steady_column`, `gradient_blind`, `surface_T_recon`, `convective_patch` |
| (d5) `jet_T_ent_mode = rec_fixed` | none hashed (`jet_d2c` KT/K15/A17–A19 run but are not hashed) | 0 | retire the legs from `jet_d2c` |
| (d6) `pad_shield_open`, `jet_above_h` / `jet_T_above` | none | 0 | `face_heat` R_pad / R_above restated |
| (d7) `flake_length = flux`, `flake_fraction = weibull` | none | 0 | no unit test; `studies/p2_2d_kic` (80 legs) and `p2_2e_cref` (204 legs) used `flake_length = flux` — C.7.3 calls it closed at P1g-3 but two later studies ran it |
| (d8) `weibull.detach_mode = connected_cluster` (C.7.2) | none | 0 | `kant/sp_kant_onset/input_weibull_cluster`, extension `input_pulsed_brittle_400w_cc` |
| P2-10 keys, P2-16 `constraint_depth` | none (already deleted, 0/2639 held) | 0 | — |

**Consequence for the order of work.** Because (b) empties the gate, the packet sequence that keeps an
oracle at every step is: (1) port beam / vapour / melt / feet onto the face list behind keys that are
off in every gate deck (0 of 2639 holds, P2-10 style); (2) build a **new face-list gate** — P2-10's W
set without its key-0 members (`face_removal` T00_1_k1, TL35_1, TL45_1, SPH_1, BORE_1; `face_list` CAV,
STAIR; M0T) plus face-list versions of the ten gate tests (`scalar_flaw` and `beam_void_closure` on the
beam port, `robin_feet`/`jet_enthalpy`/`robin_jet` on the feet port, `face_buckle` legs F_0/F_S) — hash
it on the pre-deletion binary and declare it in `REFERENCE.md` as the new reference, keeping the 2,639
list as the record of the column era; (3) delete the column path and require 0 of N on the new list.
Do DAMAGE before COLUMN (it removes only 8 gate files if the plot fields are kept), and ZHANG at any
time (gate-neutral).

---

## 4. Extensions and studies (records — never edited)

### 4.1 Extensions

| extension | decks | lineage | removal | runs today | after the deletions |
|---|---|---|---|---|---|
| `extensions/hu_spall_onset_mmwbeam` | `input_granite2` | DAMAGE (`damage.enabled = 1`, damage_law, `shallow_max`, `xi_s`, `C_h`, `damage_threshold`), BEAM, VAPOR, MICRO (3-phase Voronoi), MECH | COLUMN | RUNS | stops at DAMAGE |
| `extensions/hu_spall_onset_mmwbeam_lefm` | `input_granite2`, `input_pulsed_brittle{,_400w,_400w_cc,_400w_hires,_hires,_mid,_ultra}`, `input_wellbore_3d{,_fine,_smallgrain}` (12) | LEFM+Weibull, BEAM, VAPOR, MICRO (3-phase Voronoi), LTE; `flake_coherence_length` on 5, `connected_cluster` on `_cc`; all name `damage.enabled = 0` | COLUMN | all 12 RUN | stop at COLUMN unless re-expressed after the beam/vapour port; abort earlier if `damage.enabled` goes on an abort-on-contains list |
| `extensions/hu_spall_onset_voronoi` | `input_granite2` | DAMAGE, PATCH-T, MICRO (Voronoi), MECH | COLUMN | RUNS | stops at DAMAGE |

### 4.2 Studies

Static scan (`work/STUDIES_SCAN.tsv`): a recorded `argv*.sh` naming a P2-10 or P2-16 retired key, or a
`surface_patch.*` key the strict parser no longer reads, **aborts today** at parse. Older studies (D2*,
P0*, P1*) recorded no `argv.sh`; their legs are built inside `run.py` from base decks, so their
"driver names a retired key" flag means the legs that pass it abort.

**(i) Cannot run today, wholly or partly (34):**

| study | recorded argv | abort today | other to-be-deleted families named |
|---|---|---|---|
| `p2_2e_cref` | 221 | 212 name `spallation.sp.constraint_depth` (P2-16), 12 name P2-10 keys, `side_face_h_scale` (strict) | `flake_length = flux` (204) |
| `p2_14a_cref_face` | 24 | 11 (`constraint_depth`) | — |
| `p2_7_hoeser_heat` | 15 | 15 (P2-10), 3 also P2-16 | — |
| `p2_5_hoeser`, `whatif_uptake` | 11, 6 | all (P2-10) | — |
| `p2_8_face_list`, `p2_9_face_removal` | 18, 16 | 18, 14 (P2-10: their legs set `spall.lateral_*` / `side_face_*`) | — |
| `p2_3_wall`, `p2_4_body` | 6, 8 | all (P2-10; `side_face_h_scale` strict) | — |
| `p2_2c_meier` | 8 | 6 (P2-10) | — |
| `p1g3_thickness` | 4 | 4 (P2-10) | `flake_length = flux`, `flake_fraction = weibull` |
| `p2_10_retire` | 33 | 6 — **deliberate** retired-key legs (`retired/R_*`) | surface_normal, damage (driver) |
| `p2_12_retire_path` | 20 | 0 retired; one deliberate `bogus_key_xyz` strict-abort leg | — |
| driver-only (no argv recorded; `run.py` names P2-10 keys): `d2aa_heatdelivery` (+ `pad_shield_open`), `d2ab_gaph`, `d2ac_feed`, `d2k_sideface`, `d2l_scan`, `d2q2c_arbiter`, `d2r1_aboveplane_sides`, `d2t_radial_taper`, `d2w_attribution`, `d2x_areasplit`, `d2z0_increment`, `d2z2_multicell`, `d2z_stepheight`, `p1b_support_rule`, `p1c_isotropy`, `p1e_staircase`, `p1f_lateral`, `p1g2b_subcell`, `p1g2c_shiftfix`, `p2_0_annulus`, `p2_1_wallh` | — | legs passing the P2-10 keys | `p1b`/`p1c`/`p1e` also pass `spall.driving_*`/`support_*` (silently unconsumed since P2-2a) |

`docs/project/model_card/RUNNABLE.txt` confirms by running: P2-2e G2.0_130 / V80 / SD, P2-7 HG,
whatif s100, P2-2f ctrl8/m13 all **ABORT (rc 6)**.

**(ii) No retired key found; column path (stop running at COLUMN deletion) (35):** `amr_band`,
`c1_pinned_closure`, `d2a2_walljet`, `d2a_jet_face`, `d2ad_fluxcal`, `d2b_feet_rop` (+ `surface_normal`),
`d2c_steady`, `d2e_mesh` (+ `rec_fixed`), `d2f_support`, `d2g_mesh`, `d2h_openloop` (+ `rec_fixed`),
`d2i_gate`, `d2j0_hotdisk`, `d2j0b_plane` (+ `surface_normal`), `d2m_feet`, `d2o0_wallsplit`,
`d2o0b_wallsplit`, `d2q1_percellvoid`, `d2q2b_percell`, `d2r0_annulus`, `d2r0c_aboveplane`, `d2v_hmap`,
`d2z1_delete`, `heat_budget`, `law_form`, `p0_surface_T`, `p0b_one_surface_T`, `p1_gradient_criterion`,
`p2_11_research`, `p2_2_props` (driver passes `spall.driving_*`), `p2_2b_props` (24 legs on
`kant/sp_kant_onset_meiercfg/input_model`), `p2_2d_kic` (83 legs, same base, 80 with `flake_length =
flux`), `s1_surface_resolution` (its smoke legs are gate legs; `melt_aware = 1`), `s1b_weibull_firing`
(+ `surface_normal`), `tree_sweep_0922` (the gate itself). Their bases are the Meier decks
(`input_feet`, `input_feet_d2c`, `input_drilling`, `input_2d_dev`), `hu/hu_end_to_end/input_granite2`
(⇒ DAMAGE) for several D2*/P1* drivers, and `kant/sp_kant_onset_meiercfg/input_model`.

**(iii) Face list (survive COLUMN deletion if no other deleted key is set) (5):** `p2_11_heat` (3/3
legs `face_list=1`), `p2_11b_path` (13/14), `p2_15_buckle` (26 legs; 21 on
`kant/grimsel_calibration/input`, which sets `face_list = 1` via the shared include), `p2_16_span_rule`
(7), `p2_17_buckle_on` (34; `face_heat_output_pre/R_pad`, `R_above`, `R_robin` are face_heat refusal
legs naming `pad_shield_open`, `jet_above_h`, `robin_form=face`). ⚠ Many of their argv name
`spall.flake_coherence_length=0.0` and `spall.surface_normal=0` (14 of P2-17's): an abort-on-contains
list for those two keys stops them too.

---

## Gate and test retirement plan

| deletion | unit tests to retire (→ `retired_units/`) | unit tests to re-point | validation decks | gate decks to drop from `hash.sh` / lanes | new lock clause replacing the old assertion |
|---|---|---|---|---|---|
| **(a) DAMAGE lineage** (`damage.*`, `cohesive.*`, damage_law, `column_max`/`shallow_max`, `prescribed_ratio`, `spall.alpha`, `t_spall`, `xi_s`, `C_h`, `damage_threshold`, κ-damage) | `dp_yield`, `gb_cohesive`, `spall_event`, `regime_low_high_power` | `amr_microstructure_regrid` (drop D/`kappa_eff` clauses), `heterogeneous_kappa` (drop damaged half) | retire `hu_spall_onset` (+ `diag/`), `hu_end_to_end`; extensions `hu_spall_onset_mmwbeam`, `hu_spall_onset_voronoi` become dead records | `spall_event` output + output_no_event (8 files); ⚠ keep `D`, `kappa_eff`, `k_eff`, `regime` registered or re-pin all 2,466 `Cell_D` | `retired_keys`: R-clause "no `DamageGrowthRHS`, `IntegrateDamageRK4`, `UpdateDamageAfterMechanics`, `UpdateCohesiveAfterMechanics`, `PrincipalStressRatio`, `damage_enabled`, `cohesive_enabled` in `src/` outside the list" + K-leg "`damage.enabled=1` aborts naming the step"; `spallation.model=damage_law` and `spall.sample=column_max` refused by value; `spallation.model` default/required asserted |
| **(b) COLUMN path** (`Removal.H`), after the beam/vapour/melt/feet port | `removal_increment`, `seed_void` (k_top scans), `surface_T_recon` (cell_average/diagnostic/cell modes), `robin_pinned` (unless pinned is ported), `p0`-style controls | `beam_void_closure`, `scalar_flaw` (beam port); `jet_enthalpy`, `robin_jet`, `robin_feet`, `jet_d2c`, `pad_clear_ledger` (feet/flame port); `sp_v_n_regime` (vapour port); `sp_melt_skeleton` (melt port); `gradient_blind` (+ its role as lock fixture); `robin_face` ledger legs; `face_removal` (drop key-0 legs and (a3 i) column half); `face_list` (drop OFF0/CAVU (iv)); `face_buckle` (drop F_S_k0) | re-point `kant/sp_kant_onset_meiercfg` (successor: `grimsel_calibration`), `p0_steady_column`, Meier `input_feet*` (feet port), `input_2d_dev*`/`input_energy` (beam port); `input_drilling*` UNCLEAR (bit not on the port list); `rossi` UNCLEAR (PATCH-T not on the face list) | **all 1,963 removal files** — retire the 2,639 list as a record and switch to a new face-list reference built **before** the deletion (W without k0 + the ported gate tests); `WHICH_TESTS.md` rewritten | `lateral_flux_lock`: stop opening `Removal.H`, delete/restate the 48 `Removal.H`-reading sites, re-point its fixture; `buckle_lock`: drop B2's `Removal.H` half and B3 if the key goes; `retired_keys`: drop the `Removal.H` `SUBST`, re-point fixture; NEW clause "no `UpdateRemovalAfterCohesive`, `Removal.H`, column-only names in `src/`"; `face_removal` (a3 i) restated as "declared once, called by the line path" |
| **(c) ZHANG separate path** (`AdvanceMaterial`, `AdvanceMaterialImplicit`, `AdvanceConstAlpha`, `material.type`, `solver`) | none | `beam`, `equilibrium`, `stefan`, `convective_patch` (onto one-phase microstructure) | MERGE `zhang_oglesby`, `hu_conduction`, `hu_thermoelastic`, `hu_thermoelastic_amr`, `hu_breakage_index` onto the microstructure path | **none** (0 files; gate as-is is the witness) | `lateral_flux_lock` or `retired_keys`: the three Advance functions and `SolverKind::Implicit` are gone; `material.type` aborts; ⚠ 75 decks name `solver = explicit` — abort only on `implicit`, or edit them |
| **(d1) `spall.surface_normal`** | — (SN leg of `removal_increment`) | `beam_void_closure` (re-pin), `face_heat` BASE token | `meier/input_2d_dev*` re-pin | `beam_void_closure`, `dev2d` (233) — re-pin, or fold into (b) | name clause: `spall_surface_normal`, `sn_recession_superseded`, `icth` absent; K-leg aborts; ⚠ 21 decks name it `= 0` |
| **(d2) `flake_coherence_length`** | — | as (d1) | Meier beam decks, `kant/sp_kant_onset_meiercfg`, `p0_steady_column` | same 233 — fold into (b) (column-only already) | name clause `spall_flake_coh_len`; ⚠ 7 decks + face-list argv name it `= 0.0` |
| **(d3) `melt_aware`** | `sp_melt_skeleton` (unless ported) | — | `meier/input_2d_dev_melt` | measure: 0 of 2639 iff `Λ_L = 0` at every firing; else re-pin the 1,600 melt-deck files | name clause `melt_aware`, `melt_lambda_crit`; resolve the "closed vs ported" contradiction first |
| **(d4) `robin_form = cell / face`** | `robin_face` (or re-point to `surface`) | `robin_pinned` F leg, `jet_enthalpy` strip legs, `convective_patch`, Kant onset decks (default) | Kant onset ×7, `p0_steady_column` | 1,003 files (`robin_face` 798, `robin_pinned` F 45, strip 118, s1 R-legs 42) — re-pin after choosing a new default | value-refusal clause ("must be pinned or surface"), new default asserted |
| **(d5) `rec_fixed`**, **(d6) `pad_shield_open`, `jet_above_h`/`jet_T_above`**, **(d7) `flake_length=flux`/`flake_fraction=weibull`**, **(d8) `connected_cluster`** | `jet_d2c` legs KT, K15, A17–A19 (d5) | `face_heat` R_pad, R_above (d6) | `kant/sp_kant_onset/input_weibull_cluster` (d8) | none (0 files) | one name/value clause each in `retired_keys`; gate as-is (0 of 2639) |

---

When a deletion empties a gate deck, the packet has to say so before it runs, list every affected path,
and re-pin deliberately (D2z-1's procedure). A deletion can only claim "0 of 2639" when its code is off
in every gate deck, and among the deletions above only ZHANG and (d5)–(d8) meet that.
