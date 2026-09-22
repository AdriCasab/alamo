# ROADMAP.md

This is the short project map. Claude should read this at startup, then read
`ACTIVE_STEP.md` for the detailed coding packet. The full historical plan remains
in `in-main-tex-you-will-quizzical-treasure.md`.

## Current Status

Completed and passing:

- Step 1: MMWSpalling skeleton integrator and executable.
- Step 2: Enthalpy formulation, phase fractions, and Stefan validation.
- Step 3: Gaussian/Beer-Lambert MMW beam source and beam-energy validation.
- Step 4: Radiation/convection surface losses and equilibrium validation.
- Step 4b: Zhang/Oglesby granite-heating thermal validation with material
  abstraction, time-varying beam power, T-dependent emissivity, tabulated H/T
  inversion, and explicit/implicit thermal switch.
- Step 5: Voronoi mineral microstructure, phase-property fields, and
  grain-boundary flags.
- Step 6: Heterogeneous and damage-modified conductivity, passive damage field,
  GB conductance hook, and effective-conductivity diagnostics.
- Step 6b: Microstructure H/T consistency with reusable material enthalpy
  tables, conservative heterogeneous enthalpy update, and H-based phase
  fractions. Fast regressions passed; `zhang_oglesby` deferred.
- Step 6c: Grain topology correction with explicit `grain_id`,
  `is_grain_boundary`, `is_phase_boundary`, `is_gb` kept as a phase-boundary
  alias, and `h_gb0` keyed on true grain boundaries. Fast regressions passed;
  `zhang_oglesby` deferred.
- Step 6d: Hu prescribed-surface-temperature conduction validation with
  `surface_patch.*`, `material.type = constant`, Granite 2 / Sandstone 2
  runs, and ALAMO-only Hole 1 comparison plot. Fast regressions passed;
  `skeleton` and `zhang_oglesby` deferred.
- Step 7: Heterogeneous thermoelastic mechanics with microstructure-driven
  `mu_phase`, `T_ref_phase`, nodal stiffness/eigenstrain refresh, and tightened
  `thermal_stress` regression. Hu conduction and fast regressions passed;
  `zhang_oglesby` deferred.
- Step 7b: Hu Granite 2 vs Sandstone 2 thermoelastic stress validation,
  single-level. Homogeneous thermoelastic inputs and regression pass with
  documented bottom-clamp BC and L3 reporting-only caveats.
- Step 7c: Early AMR correctness for microstructure regrid, including regrid
  and post-average-down repair of discrete topology, phase-property, passive
  `D`, `kappa_eff`, and `k_eff` fields. AMR microstructure regression passes.
- Step 7d: Hu thermoelastic stress validation with AMR, using homogeneous
  Granite 2 / Sandstone 2 `20^3 + max_level=1` reruns, raw-cell checks against
  uniform `40^3`, and preserved bottom-clamp/L3 caveats.
- Step 8: Drucker-Prager failure criterion and continuous damage evolution with
  default-off parser/settings, irreversible RK4-style `D` evolution after
  mechanics, damage-scaled microstructure stiffness, evolved-D AMR repair, and
  the `dp_yield` regression. Hu and nearby regressions passed;
  `zhang_oglesby` deferred.
- Step 8b: Hu breakage-probability indicator without damage evolution. Granite
  and Sandstone surface `f_b` targets pass after AMR refinement, surface
  extrapolation, and a top-node `surface_patch` eigenstrain fix; L6 deep check
  remains relaxed under the bottom-clamp BC. Also added default-off quartz
  alpha-beta transformation strain with `alpha_beta_transition` regression.
- Step 9: Grain-boundary cohesive-zone unit test. Added default-off bilinear
  CZM law/history/diagnostics on true `is_grain_boundary` cells with the
  `gb_cohesive` regression; traction feedback into the elastic operator is
  explicitly deferred.
- Step 10: Surface advancement and spall detachment. Added default-off
  `spall.*` removal, level-set `phi`, `removed`, event/RoP diagnostics,
  void-aware microstructure heat advance, and `spall_event` regression.
  Verification is single-level with per-column connectivity.
- Step 11: Vaporisation removal and unified rate of penetration. Added
  default-off `vapor.*` removal, shared `phi`/`removed`/`regime`/`RoP`
  diagnostics, per-mode spall/vapor diagnostics, and the 4-rank
  `regime_low_high_power` regression. Per-column removal now uses domain-wide
  MPI reductions, so z-split columns are handled correctly.
- Step 11b: Hu LRST/onset validation with ALAMO DP damage. Added
  `hu_spall_onset`, one-phase homogeneous Granite 2 / Sandstone 2 inputs,
  microstructure `surface_patch` heat support, and a 4-rank validation script.
  LRST/onset/order pass in v1; Sandstone damage depth is warning-only and AMR
  onset/removal remains deferred. (Test inputs were retargeted to
  `amr.n_cell = 32 32 32` and the sandstone tolerance window to the
  Hu-primary case during Step 12; see Step 11b notes below.)
- Step 12: Spallation-model toggle (`damage_law` default / `sp_weibull`),
  Sp = K_I/K_Ic(T) onset criterion (Tada weight + Nasseri table), and a
  Kant 2017 closed-form verification at machine precision. EOS dropped
  from the step indefinitely. Working-tree mechanics-stack cleanup
  (Operator/, Solver/Nonlocal/, partial Mechanics.H + Newton.H) was
  performed hunk-by-hunk under user authorisation to restore Step 11b's
  calibration baseline. Full green-light regression sweep PASS 9/9.
- Step 14: Hu flame-jet end-to-end validation under `damage_law`
  (combines 6d, 7b/7d, 8b, 11b on Granite 2 / Sandstone 2 plus Step 10
  spall removal on a single integrated run). 32³ single-level.
  Granite onset 40.0 s, sandstone onset 90.5 s — both inside Hu-primary
  windows; ordering correct. Test-only step (zero `src/` changes). 10/10
  green-light regressions PASS. `spall_event > 0` firing initially
  relaxed to warning-only; Step 14b resolved it.
- Step 14b: Damage-law spall sampling knob + `PrincipalStressRatio`
  patch. Added `spall.sample = top_cell | column_max | shallow_max`
  (default `top_cell`, byte-identical to Step 10/11) and
  `spall.shallow_depth`; `UpdateRemovalAfterCohesive` reduces `D` over
  the configured per-column scan window. Post-completion patch:
  `PrincipalStressRatio` returns a `1e30` sentinel on mixed/tensile
  principal-stress states instead of 0. `hu_end_to_end` now fires
  spall for both rocks (hard pass); 10/10 green-light regressions PASS.
- Step 15a: convective Robin surface-patch BC
  (`surface_patch.mode = prescribed_T | convective_flame`,
  `surface_patch.T_flame`) + the `convective_patch` Carslaw-Jaeger
  verification. Robin flux `h_fl*(T_flame-Tc)*inv_dz` on the in-patch
  top row in both explicit heat paths; Step 8b top-node eigenstrain
  override gated to `prescribed_T`. Verification PASS (<=2.6% error);
  prescribed_T regressions bit-identical. 11/11 regressions green + the
  new test.
- Step 15b (Sp mechanics): the machinery half of long-plan §15b. Added
  `weibull.{enabled,a0_gb,a0_ig,m,seed}` (default-off), the cell-centred
  `flaw_a_mf` field sampled at Initialize via a global-cell-index
  splitmix64 hash, the `-sigma_xx` integrand sign fix, depth-resolved
  physical-path `sigma_xx` sampling, and per-cell `a_f` in
  `UpdateSpAfterMechanics`. Verified by the new `sp_weibull_unit`
  synthetic test; `sp_onset_kant_closed_form` (§12c) updated for the
  sign convention (the only existing test changed). `damage_law` sweep
  bit-identical (`hu_spall_onset` 37.6/93.0 s, `hu_end_to_end`
  40.0/90.5 s). 13/13 green.
- Step 15b (Kant onset validation, test-only — zero `src/` changes):
  added `tests/MMWSpalling/sp_kant_onset/` with 4 inputs
  (`input_verification` frozen `K_Ic⁰=1.5`/`h_fl=150`,
  `input_verification_hfl` `h_fl=250`, `input_model` `K_Ic(T)` Nasseri
  `×1.05`, `input_weibull` single-phase Voronoi + `weibull.enabled`)
  and a self-running test. Authorized packet amendment with the user:
  the onset closed form is **Kant Eq. 15** `Θ_onset = K_Ic·(1−ν)/(2·
  1.763·√(a/π)·E·α)` (`h_fl`- and `λ`-independent), NOT Eq. 14
  `Sp_red`; Kant's 390–560 °C Table-2 band was computed with frozen
  room-T `K_Ic⁰`, so target 1 is scored against frozen `K_Ic⁰` (the
  `K_Ic(T)` run is the model's honest prediction, reported but not
  banded). Setup: lateral roller box (`el.bc.type = constant`, all 26
  regions) + full-face `convective_flame` enforces `ε_xx=ε_yy=0`
  exactly. 6/6 targets PASS — onset ΔT 461.3 °C ∈ [390,560];
  FEM-vs-Eq.15 3.5%/3.2% (tol 8%); a-sensitivity √2;
  `h_fl`-independence 2.2%; weibull patchy first firing 14/576 cells
  with GB fraction 1.00. `hu_spall_onset` bit-identical.
- Step 15c (surface labeller + detach_mode A/B + regrid-repair + v2
  dual-scoring sweep): two implement passes. v1 shipped the
  `DetachMode { PerFace, ConnectedCluster }` enum,
  `weibull.detach_mode`/`weibull.A_crit` parser, `Sp_cluster_id_mf`
  top-z plane plotfile field, `UpdateSpClusters` 4-connected
  union-find labeller, per-step `<plot_file>_clusters.csv` diagnostic,
  AMR regrid-repair for `flaw_a_mf` (splitmix64 hash of global cell
  index → bit-identical), new `sp_kant_onset/input_weibull_cluster`
  A/B partner and `sp_weibull_unit/input_regrid` AMR check, plus the
  25-case `(a0_gb, m)` sensitivity sweep machinery. v1 sweep PASS gate
  was rejected by review on a physics mismatch (Weibull first-fire is
  governed by `a_max ≈ a0·(ln N)^(1/m)`, not `a0_gb`); v2 (test-only,
  zero `src/`) rewired `run_sensitivity_sweep` to score
  `|dT_FEM − Eq.15(a_max_realized, K_Ic(T))| / Eq.15(a_max_realized)
  ≤ 8%` on the 9 central-third cases (`m ∈ {12,15,20}` ×
  `a0 ∈ {0.75×, 1.0×, 1.5×}·20 µm`). 9/9 central-third PASS at
  2.8–4.5% rel_err; realized `a_max/a0` matches Galambos asymptotic
  to ~1%; 6/6 §15b targets unchanged; A/B onset timing identical,
  cluster_count 2→0, total_firing_area 2.0→0.0 mm² at A_crit; AMR
  regrid-repair bit-identical. `damage_law` sweep bit-identical.

- Step 16 (Sp-vs-depth `h_spall` mechanics, revised-approach Stage 3,
  **mechanics half only**). Wires `sp_weibull` material removal end-
  to-end: per-column K_I-vs-K_Ic depth scan, Sp + zero-crossing path
  replaces the `damage_law` `D`-threshold + `C_h·√(α·t)` scaling under
  `spallation.model = sp_weibull`, new `h_spall_field_mf` running-max
  diagnostic, gate consumes the §15c `Sp_cluster_id_mf > 0` labeller
  output. `damage_law` byte-identical. Geometry deviated from packet
  (MLMG conditioning + free-lateral physics): working config is
  `32^3 cubic cells in 50×50×50 mm + Kant-style lateral roller box +
  full-face heating`. P1-P4 mechanics gates PASS at end of 130 s
  ramp; Rossi quantitative bands not gated (running-max diagnostic
  ≠ Rossi's spatial distribution).
- Step 16b (Rossi crack-event count + sub-cell binary search). All
  5 deliverables shipped per design note + user-approved decisions
  (edge-triggered first-crossing per cell, option (ii) freeze
  re-firing without σ touch, 8-iter bisection bracketing
  `s_last_pass` and `s_last_pass+1` refining `h_col` to dz/256).
  Test FAILs R1 (peak at 47 µm not 100–200 µm) and R3 (no events in
  100–520 µm zone): 528 firing columns all land in the top z-cell
  because at dz = 1.5625 mm, K_I drops sharply across one cell — sub-
  cell crack tips are recorded in `h_spall_field_mf` (bisection-
  refined `h_col ∈ [1.538, 1.556] mm`) but NOT in the cell-centred
  `cec`. Two new findings: (1) σ_xx local-per-z under 1-D
  confinement makes option (i) ineffective on this geometry; (2)
  24:1 GB/IG (bisection lifted Step 16's 100 %) vs Rossi 5:1 is a
  separate gap that depth-binning won't fix. §15c + `damage_law`
  regressions PASS unchanged.
- Step 16c (per-step h_col CSV log + linear cell-to-cell σ_xx,
  option (L)). All 5 deliverables shipped per spec; **R1 BINDING
  gate FAILs on sub-criterion (1)** — a new physics-mechanism
  finding, not a delivery gap. 528 first-firing events captured
  in `<plot_file>_h_col_events.csv`; h_col distribution narrow at
  757-787 µm (mean 769, std 5.8) — ALL events in bin centred at
  799 µm. The flat-extrapolation of σ_xx in `z ∈ [0, dz/2]` (no
  adjacent cell ABOVE the top cell) pins `h_col ≈ dz/2 − a_f`
  for every column (verified arithmetically to within bisection
  tolerance). Sub-criteria R1 (2) + (3) PASS; R2 WARN. Truncation
  analysis on the CSV verified peak invariant under stop_time
  → packet protocol step 2 (flag for reviewer) the active path.
  Path forward = Step 16d (option (a), T-derived σ_xx);
  `rossi-validation-diagnostic-design.md` §7-§10 captures the
  (L) artefact, Step 16d scope, mesh-resolution limit reasoning,
  and anti-patterns. §15c + `damage_law` regressions PASS unchanged.

- Step 16d (T-derived σ_xx with prescribed-T face BC, option (a),
  revised-approach §4). One-lambda swap in `k_i_at_zcrack` inside
  the SpWeibull spall arm of `UpdateRemovalAfterCohesive`: T(z)
  linearly interpolated between `surface_patch.T_f(time)` at z=0
  and cell-centred `temp_mf` at z = (n+0.5)·dz; σ_xx from 1-D-
  confinement thermoelastic relation `+E·β·(T−T_ref)/(1−ν)`
  (compression-positive double-negation matching
  `UpdateSpAfterMechanics`). Material-param capture from phases[0]
  with four fail-fast aborts (empty phases / missing E or μ /
  unphysical ν / mode != prescribed_T). All 4 deliverables shipped
  per spec; **R1 FAILs on sub-criterion (1) with peak at 705 µm**
  — the §8 decision tree's "peak >> 200 µm → accept mesh-
  resolution conclusion" branch is the active path. Sign-
  convention bug caught during implementation (returned negative,
  should be positive to match K_I integrator's compression-positive
  convention) — fixed before release. R1 (2) + (3) PASS; R2 WARN;
  R3 informational (24:1 unchanged). 528 first-firing events,
  h_col distribution broad [573, 1556] µm (vs 16c's narrow pile
  at 757-787 µm — option (L) artefact is gone). §15c + damage_law
  regressions PASS unchanged. **Three-packet Rossi chain (16b/
  16c/16d) closed; R1 deferred per §10 anti-patterns + §11 final
  status in `rossi-validation-diagnostic-design.md`.**

- Step 17 (per-cell `v_n` with simultaneous spall + vap, revised-
  approach Stage 4 — mechanics half). All 5 deliverables shipped
  per spec, accepted by review. Under `sp_weibull`: replaced
  winner-takes-all-on-h with spall-priority precedence in PASS 2's
  column-winner block; replaced bulk-averaged `RoP_vap = oneR·P/
  (rho·L_tot·A_beam)` with per-cell local-Gaussian `Q(x,y,t) =
  (2P/(π·w0²))·exp(-2·r²/w0²)` in the vapor-RoP path
  (three-branch ladder: `vapor_prescribed` > `sp_weibull` >
  `damage_law` bulk). New `regime_field_mf` cell-resolved
  diagnostic ({0/1/2}) registered alongside the existing scalar
  `regime_mf`. Damage_law else-branches preserved bit-for-bit
  (verified by `hu_end_to_end` PASS at the exact Step 16d
  baseline — granite 40.0 s, sandstone 90.5 s, any_spall_event
  for both). Synthetic test `tests/MMWSpalling/sp_v_n_regime/`
  PASS: 524 firing-cell consistency observations, Pearson r =
  +0.988 between -r² and RoP_vap (perfect Gaussian profile
  confirming local-Q wired); the "simultaneous regimes" sub-check
  downgraded to informational because Gaussian-beam-in-lateral-
  roller geometry can't produce vap-without-Sp (Step 16 patch-vs-
  1D σ_xx issue). §15c + damage_law regressions PASS unchanged.
  **Step 17b (full Zhang/Oglesby reproduction under sp_weibull)
  deferred** — two specific obstacles in the next planner's
  scope: Zhang reports T_surface not ROP; `material.type = zhang`
  doesn't map to `microstructure.phase0` scalars.

- Step 18 (Kant 2017 confining-pressure sweep validation,
  revised-approach Stage 5). New parser key `confining.p`
  (default 0, byte-identical to all pre-Step-18 inputs);
  `+p·ν/(1−ν)` offset added to compression-positive σ_xx in
  both `UpdateSpAfterMechanics` and Step 16d's `k_i_at_zcrack`
  lambdas. Single base input + command-line override drives
  three sims (p ∈ {0, 27, 48} MPa). **All gates PASS** — R1
  binding 3/3 within 8% (3.5/3.8/4.1% rel err vs Eq.15), R2
  binding monotone (461.3 > 436.6 > 417.4 °C), R3 info p=0 in
  Kant [390, 560] band, R4 info FEM slope −0.913 vs Eq.15
  −0.929 K/MPa (1.7%). Mechanics solve unchanged — confining
  pressure is a post-solve offset, not a new BC. Caveat:
  `stress_xx` plotfile still reflects thermal-only σ_xx; the
  confining contribution lives only in the σ_xx-lambda returns
  feeding the K_I integrator. **Kant validation chain (12c →
  15b → 15c → 18) closed; revised-approach Stage 5 complete.**
- Step 19: sp_weibull spall under MMW beam (relax Step 16d
  surface_patch gate). Two-hunk src change in
  `UpdateRemovalAfterCohesive`: removed the unconditional
  prescribed_T abort, added `sp_weibull_use_T_face` flag
  (requires `surface_patch.enabled` — fixes a hidden T_face=∞ bug
  from calling an uncompiled parser) + fail-fast on missing
  stress, and forked `k_i_at_zcrack`'s σ_xx between the Step 16d
  T-derived form and `stress_mf` NodeToCellAverage. New
  model-extension test `hu_spall_onset_mmwbeam_lefm` (NOT a
  validation). G1+G2+G3+G4 PASS (first spall at t=2.0 s, Sp_max
  2.91, 207 cells removed); all prescribed_T/damage_law
  regressions byte-identical. Accepted by review.

- **Refactor R1 — extract `UpdateRemovalAfterCohesive` into a partial
  header (COMPLETED, archived).** Moved the ~957-line method *verbatim*
  into `src/Integrator/MMWSpalling/Removal.H` (`#include`d inside the
  class body; stays an inline member). Byte-exactness proven by
  range-diffs vs HEAD; all removal witnesses byte-identical
  (`sp_rossi_damage_profile`, `hu_end_to_end`, `spall_event`,
  `regime_low_high_power`, `sp_v_n_regime`, `sp_kant_onset`). Data members
  + `PrincipalStressRatio` deferred to R2. See `docs/project/ARCHIVE_DONE.md`.

- **P1 — Melt-aware spallation (COMPLETED, archived).** Default-off
  `spall.melt_aware` + `melt.lambda_crit` scale the sp_weibull driving stress
  σ_xx by a load-bearing skeleton `f(Λ_L)=max(0,1−Λ_L/Λ_crit)` (Λ_crit=0.4) in
  `Removal.H`'s `k_i_at_zcrack` (both branches) so the criterion stops spalling
  liquid; the scan self-anchors at the solid–liquid interface (explicit anchor
  dropped — would evacuate the melt cap = P2). Byte-identical when off (6
  witnesses PASS). New clean A/B unit test `unit/sp_melt_skeleton/`. KEY
  finding: the spall thermostat pins the surface at ~850 K, so mid-block melt is
  unreachable without α(T) (P3) — the Meier melt demo is bottom-confounded.
  See `docs/project/ARCHIVE_DONE.md`.

- **E1 — Beam/void gap energy leak fix + deposition closure counter (COMPLETED,
  review-accepted, archived).** Beam path = `(n_above+½)·dz` with n_above
  counted from the `removed` mask (`ColumnTopSolid` + `SolidFacePath`). The
  `floor(φ/dz)` form failed at round-off ties (φ = dz exactly, which is
  common). Per-step `BeamEnergyBalance` thermo counter (closure ~1e-15), φ/mask
  invariant check, and new `unit/beam_void_closure`. Meier 2D: closure 0.995,
  ROP flat 6.3 m/h (was 3.4 m/h and declining); ROP = f_flake·q/(ρ·Cp·ΔT_rem)
  holds to ratio 1.00 at ΔT_rem = 519 K. Kant/Rossi/Hu byte-identical. See
  `docs/project/ARCHIVE_DONE.md`.

- **S1 — 1-D surface-resolution study (COMPLETED, review-accepted, archived).**
  Default-off `surface.follow_mask`, `surface_patch.robin_form = cell|face`,
  `energy_ledger.enabled` (exact H ledger), new `unit/robin_face`, study
  harness `studies/s1_surface_resolution/` (35 runs). Findings:
  - the criterion is a cell-average T threshold (656 K at a = 20 µm);
  - h_col is a removal increment ∝ dz, not a flake;
  - fixed-flux ROP is converged at 2 mm;
  - Robin is not converged at 2 mm (cell +42 %, face −82 %; finest-dz bracket
    14.0–18.6 m/h), and the converged form is the pinned surface
    q = h(T_gas − T_fire);
  - no melt anywhere, so P1 stands.
  See `docs/project/ARCHIVE_DONE.md` and `Claude_markdowns/2026-09-15b.md`.

- **A1 — Removal fixes + S1b firing temperature (COMPLETED, review-accepted,
  archived).**
  - `weibull.enabled = 0` + removal works (scalar `a_f`, `Sp_top >= 1` gate;
    bit-identical to the degenerate draw).
  - New `spall.removal_events_csv` logs every removal (rank-count identical;
    key on `regime`).
  - Face form aborts if ε_high < ε.
  - S1b: the Meier Weibull block fires at 810.8 K at 2 mm (Meier 2D 812 K).
    V0 = 1e-9 drifts +9–11 K per halving; V0 = dz³ is flat at ≈ 821.6 K.
  - **User decision: V0 = V_cell for all C/D inputs** (defaults unchanged).
  - See `docs/project/ARCHIVE_DONE.md`.

- **C1 — Pinned-surface Robin closure (COMPLETED, review-accepted,
  archived).**
  - `surface_patch.robin_form = pinned`: T_s = min(T_face, T_pin(col)),
    face fallback after `pinned_idle_cycles` (default 2).
  - 1-D 2 mm: 15.67 m/h, inside S1's 14.0–18.6 bracket and 0.999× closed
    form, < 0.6 % per halving.
  - 2-D Meier harness (V0 = V_cell): T_fire 821–823 K, centre ROP = closed
    form, 2 → 1 mm −0.7 %.
  - Caveats: the closed-form match is mostly a consistency check; the min rule
    was barely exercised; a sharp patch edge leaves a never-firing rim pair.
  - See `docs/project/ARCHIVE_DONE.md`.

- **D2a — Impinging-jet face source, prescribed nozzle path (COMPLETED,
  review-accepted, archived).**
  - `surface_patch.h_expr` / `T_flame_expr` in (x, y, r, s, t);
    `nozzle_z0`, `nozzle_feed`, `nozzle_collision_radius`;
    `BuildFlameColumns` evaluates once per column per step and the kernel,
    pin rule, ledger and diagnostics all read those vectors.
  - New thermo `patch_min_standoff`, `patch_P_robin`,
    `pinned_cols_face_unset/idle/qneg/minrule`. New `unit/robin_jet`.
    25,887 witness files hash-identical.
  - Study `studies/d2a_jet_face/`: Martin (1977) h(r, s), anchors J-M / J-5 /
    J-10, 7 × 2-D + 5 × 3-D quarter-domain runs.
  - **Findings:** under a prescribed feed the centre ROP is an identity and
    carries no information about h or T_gas; steady bowls form at 1.5 and
    3 m/h after ~150 s; C1's min-rule and rim nits closed; the AMReX
    `a/max(x,c)` parser bug found.
  - **Two readings superseded by the 09-16 review** — see the wall-jet note
    in Known Notes.
  - See `docs/project/ARCHIVE_DONE.md`.

- **D2a2 — Energy-conserving wall-jet T_gas(r, s) (COMPLETED,
  review-accepted, archived).**
  - `surface_patch.jet_closure = enthalpy` (default `none`): `JetEnthalpyMarch`
    marches `m·cp·dT_gas/dr = −q·2πr` from `T_stag = T_ent + (T_noz −
    T_ent)·min(1, 5D/s_c)`; `jet_stagnation = decay | nozzle`,
    `jet_entrained_mass`; s ≤ 0 → no flux; in-code energy invariant; shared
    `PinRule`; 8 `jet_*` thermo columns. New `unit/jet_enthalpy`. 26,327
    witness files hash-identical.
  - **Findings:** the stagnation decay is essential (without it the centre
    runs away); face power 5.5–8.3 kW vs D2a's 21–29; in 3-D the centre leads
    the nozzle (s_c 60–90 mm), so hand tables at s_c = SOD overstate the ring
    flux 2–3×; ring ROP at r = 40 mm 0.5–1.3 m/h (early transient window,
    below Meier's band); 1900 K holes Ø 84–103 mm (≥ burner), 1436 K Ø 71–79.
  - **Recommended for D2b, before D2b:** `decay` + entrained mass.
  - See `docs/project/ARCHIVE_DONE.md`.

- **D2b — Burner-on-feet descent and the scored Meier test (COMPLETED,
  review-accepted, archived 09-17).**
  - `nozzle_descent = feet` (tripod pads, annulus [0.028, 0.040], 0.9
    quantile), `jet_T_ent_mode = exhaust`, `foot_stall_time`; new
    `unit/robin_feet`; `sp_meier_pilot/input_feet` + check-only `test_feet`.
  - **J-M untuned: ROP 1.36 m/h in band, ΔT_fire 528 K pass; but no steady
    window** (the 12 D h clip in `walljet.h_expr` lets the centre pit run to
    the domain bottom), **mouth funnel** (whole-excavation volume 2.3×
    measured), min Ø over the drilled depth 78 mm (no burner body). Fixed
    293 K stalls; 1436 K 0.67 m/h; J-5/J-10 above band (retired).
  - Mesh on a trimmed domain (−4.8 % ring) is not like for like under exhaust.
  - See `docs/project/ARCHIVE_DONE.md`.

- **D2c — Steady-state and no-funnel closures + pre-registration (COMPLETED,
  review-accepted, archived 09-17).**
  - Default-off keys: `jet_decay_diameter = momentum` + `jet_De_ref` (full
    jet), `jet_free_surface` (+ aspect, exponent; T_rec from the mouth),
    `jet_bin_update = exponential`, `jet_negative_flux = count`,
    `foot_body_clearance` (regime 4); helper far law `power`. New
    `unit/jet_d2c` 26/26; 2639 witnesses identical; `input_feet_d2c`.
  - Hashed `PREDICTIONS.md`: steady ring ROP 1.93 m/h (just above band),
    steady window from 750 s needs a ≥ 0.65 m domain; ring insensitive to n;
    T sweep 1.34/1.64/1.93 m/h at 1600/1750/1900 K; mouth predicted to fail
    (134/118 mm vs 92/88). Probes: 6.7 (2 mm) / 55.6 (1 mm) wall-s per s.
  - See `docs/project/ARCHIVE_DONE.md`.

- **D2d — Meier campaign (COMPLETED, review-accepted, archived 09-18).**
  - R1 (2 mm, 0.70 m): steady window 470–694 s (marginal), ROP 1.65 m/h and
    whole volume 2.59 cm³/s PASS; mouth 130/114 vs 92/88 mm and depth-mean
    Ø 101 FAIL; min Ø to feet 77.6 (unresolved).
  - **Mesh FAIL: 1 mm +42.8 % ROP (150–300 s, 0.40 m full domain).** All 2 mm
    Meier numbers D2b–D2d are conditional on it; the T_nozzle inference is not
    to be quoted. P6: 11 held, 1 refuted (P3 1600 K). Free-surface dilution
    closes ~half the mouth excess; the far law is a centre closure.
  - Diagnosis (09-18b): seed = whole-cell clearance clip of the 0.9-quantile
    rim band; amplifier = flank exclusion → hotter recirculation.
  - See `docs/project/ARCHIVE_DONE.md`.

- **D2e — Speed-up acceptance + mesh-seed diagnosis (COMPLETED, archived
  09-18; NOT /verify-ed).** Perf edits accepted (2639/2639), mgs 60, dt 16 ms
  (2 mm) / 8 ms (1 mm) gated: 1 mm cost 54 → 19.5 wall-s/s. New default-off
  `jet_T_ent_mode = rec_fixed`. **Seed = the whole-cell clearance clip**: with
  it off, 2 mm vs 1 mm +0.1 % (B0 +42.7 %). P-b/P-c at 1 mm stopped by the user
  to move to a physical support rule; amplifier untested. See ARCHIVE_DONE.

- **D2f — Contact support rule (COMPLETED, review-accepted, archived 09-19).**
  Outcome (b): D2g scores pads, q = 0.9, clip off. q barely matters (2 mm ROP
  1.535 / 1.491 / 1.448 m/h at q 0.8 / 0.9 / 1.0, ±2.9 %; carry ±3 %). q = 1.0
  failed one late block (−9.7 %). Q09 converged to 300 s but trending
  (−3.0, −4.3 % in the last blocks): late 1 mm drift, cause unknown. Funnel
  unchanged (gas-side). D2e reviewed retroactively 09-19 (accepted); commit deferred.

- **D2g-1 — Long mesh check of q = 0.9 (COMPLETED, review-accepted, archived
  09-20; GATE FAILED).** Matched window 150–573 s −28.5 %; blocks −1.8, −12.5,
  −40.2, −64.8 %. The centre keeps pace while the ring stalls at 1 mm (feet 229
  vs 183 mm at 550 s). **D2b–D2d Meier numbers were read where the pair is
  40–65 % apart; none may be quoted.** The 2 mm steadiness is mesh-dependent.
- **D2h — key-only diagnostic pairs (2026-09-20, reviewer; no packet, archived).**
  **The mesh-dependent element is the idle clock** (`pinned_idle_cycles`,
  timeout ∝ dz): unpinning drops a column to the self-limiting face form, lateral
  loss (∝ 1/dx) parks it sub-critical, and a silent front walks inward. Idle off:
  ROP gap −4.0 % (ref −12.5 %), silent pad columns 0 %. Fixed recirculation:
  −15.0 % — **the exhaust ratchet is not the amplifier (closes D2e's P-b)**. The
  open-loop pair aborted on a collision-guard default; not rerun.

- **D2i — Support-rule gate with the idle clock off (COMPLETED, review-accepted,
  archived 09-21; GATE FAILED).** Matched window 150–572 s **−13.9 %**; 100 s
  blocks −2.1, −4.0, −16.3, −36.1 %. Identity to D2h PASS (byte-identical to
  349.9 s). **The idle clock was a real mechanism but not the cause:** with it
  off the divergence is delayed ≈ 150 s and its slope halves, but it still runs
  away, and **idle columns are 0 at both meshes in every block** — a third
  mechanism with the same signature is underneath (ring stalls, centre keeps
  pace; feet 23 mm apart at 550 s vs a 9 mm centre floor). `pinned_idle_cycles`
  is a **threshold, not a dial**: `= 20` is byte-identical to `= 1e9` at 1 mm.
  **Stage 2 must not be planned on this rule.** See ARCHIVE_DONE.

- **D2j-0 — Hot-disk edge test, event model (COMPLETED 09-21, archived; its
  verdict WITHDRAWN the same evening).** A fixed Robin disk with no jet, feet or
  descent, at 4/2/1 mm, sharp and tapered h. As scored, P2 failed (`f_stall`
  mesh-independent, `w_edge` converging), but the 150 s window hid the growth and
  D2j-0b withdrew the "bounded edge" reading. **Durable:** sharp ≡ tapered, so
  the recession edge is not set by the sharpness of h and the gas closure
  (`jet_r_reach`) is ruled out; and the **disk-mean rate**, not `w_edge`, is the
  quantity a mesh statement must pre-register.
- **D2j-0b — Hot disk under a descending exclusion plane (reviewer diagnostic
  09-21, no packet, archived). THIS FOUND THE CAUSE.** The D2i cascade reproduces
  with **no feet, no jet, no gas closure, and the exclusion is not even
  necessary** (C0 cascades three times slower). **Mechanism: a firing top cell
  beside a taller, colder column loses k(T_top − T_wall)/dx through a step face
  the model never heats** (≈ 240 kW/m² at 2 mm, ≈ 480 at 1 mm, against q_pin
  ≈ 520); it stalls, then becomes the next wall. The threshold scales with 1/dx,
  so refinement makes it worse, at any lateral hot/cold edge. **Withdrawn by it:**
  wall heating above the plane (C1 ≡ C2 ≡ C3), `spall.surface_normal` as a remedy
  (inert at both meshes), the D2h open-loop pair, and D2j-0's "bounded edge".
  The half-cell increment is **not** a second mechanism (flat interior agrees to
  0.5 %; the whole offset lives where slope × dx ≳ the thermal length).

- **D2k — Side-face flux on exposed vertical faces (COMPLETED 09-22,
  review-accepted; arbiter PASS, Meier gate FAIL and inconclusive).** Default-off
  keys; **2639/2639 identity PASS**; unit 12/12. **The arbiter passes**: the
  D2j-0b cascade goes from −13.5 / −28.7 / −45.5 % to **+3.4 / +1.5 / −3.2 %**
  with no stalls — the step-face mechanism is real and this removes it. **The
  Meier gate fails**: rate **up 2.7×** (4.07 vs 1.49 m/h) against a pre-registered
  fall, hole only +1.4…+4.0 mm wider at matched feet height, and both legs drilled
  out of the 0.40 m box so the failing block is a bottom stop. **Not a budget
  leak** (`jet_P_face` = `ledger_P_robin` to 0.1 W, cap ratio 0.41): the jet is
  **area-limited, not enthalpy-limited**. Power 1639 → 3204 W (2.0×) while the
  rate went 2.7×, because warming the wall also removes the lateral conduction
  sink — **that second factor is correct physics and must survive any closure**.
  See ARCHIVE_DONE.

- **D2l — Gate 0 scan: why the model drills 2.7× too fast (COMPLETED 09-22;
  window EMPTY ⇒ pre-registered STOP, no code written).** `f_min` = 0.5 (the
  side-face dial must stay at half or above to keep the D2j-0b cascade
  suppressed) but `f_rate` = 0 (the rate returns to 1.3–1.6 m/h only with the
  term off). No overlap — **the wall-`h` closure is dead on its own terms.**
  **The finding that outranks the verdict:** with the side term entirely off,
  the model already removes **4.58 cm³/s against Meier's 2.47 (≈ 1.9×)** while
  the burner rate looks right — so **the D2b–D2d rate agreement was partly a
  coincidence of two errors**, which is why "ROP in band" kept coexisting with
  "Ø and volume fail". Cut Ø **78 mm** (Meier 85–93) flaring to 145 mm at the
  mouth. Method corrections: **the cascade is invisible at a single mesh at any
  f — every future arbiter run must be a pair**; the 4–6 J/mm³ "to-rock band" is
  an assumed 30 % delivery, not a measurement, and is not a target; pre-register
  **Ø(z) plus the minimum**, never one depth. See ARCHIVE_DONE.

Active step:

- **D2m — Is the cut diameter set by the feet or by the flame?** Key-only, no
  build, ~1 h. D2l's cut diameter was **flat to 0.9 mm (77.1–78.0) while
  absorbed power doubled** — a thermal diameter would move with power, so
  something geometric is pinning it, and the feet span Ø 80. Five 2 mm Meier
  legs to 250 s varying only `foot_r_outer` (32/40/44/48 mm) plus one case
  shifting the ring outward at constant width. **H-feet:** min Ø ≈ 2 ×
  `foot_r_outer` within 5 mm. **H-thermal:** min Ø stays 78 ± 3 mm. Hashed
  before launch, with the ROP direction pre-registered as a sign test on the
  reasoning itself. Meier's hole was 85–93 mm with the *same* Ø 80 feet, so if
  ours tracks the stance, the mechanism that removes rock **beyond** the feet is
  missing and no wall-`h` work will supply it. See `ACTIVE_STEP.md`.

Withdrawn (not implemented):

- **P4 — Cuttings/debris shielding.** Premise ("model ~1.9× too fast") invalid:
  the leak biased ROP low. Also, shielding η is exactly degenerate with delivery
  efficiency under a fixed-flux BC, so Meier cannot calibrate Σ₀ or show that
  shielding exists. Design kept in Known Notes for a future case with an
  independent anchor.

Previous active step (archived to ARCHIVE_DONE.md):

- Step 19 (sp_weibull spall under MMW beam). Removed the
  prescribed_T hard-abort in `UpdateRemovalAfterCohesive` and
  forked `k_i_at_zcrack`'s σ_xx reconstruction so `sp_weibull` +
  `spall.enabled` works under arbitrary heat sources (MMW beam,
  convective_flame, surface_patch off). `sp_weibull_use_T_face`
  flag gates the byte-identical Step 16d T-derived branch; the new
  branch samples `stress_mf` via NodeToCellAverage (the
  `UpdateSpAfterMechanics` pattern). New model-extension test
  `tests/MMWSpalling/hu_spall_onset_mmwbeam_lefm/`. No new src/
  files, parser keys, or plotfile fields. See ARCHIVE_DONE.md
  Step 19 for full results + the 7 takeaways (esp. #1, the
  T_face=∞ / uncompiled-parser bug).

Next milestones after the active step:

Programme source: `Claude_markdowns/2026-09-15b.md`, with the planner
amendments noted here (user decisions 2026-09-15: flake size retired, Meier ROP
next, MMW after):

- **D2n — the altered-band anchor: ATTEMPTED AND CLOSED (2026-09-22).** Two
  reasons, either sufficient:
  1. **It is not measurable.** Fig. 8.8 embeds the same greyscale JPEG in every
     available copy, and the alteration is a reddish oxidation tint with little
     luminance contrast, so greyscale loses it. A first pass giving 1.2 cm /
     1.5–2.7 kW was **withdrawn**: against a 600-position control the shadow-side
     fringe is p = 0.04 (0.11 in the lower half, control sd 3.0) and the lit-side
     fringe follows the rim lip exactly — it is lighting, not rock. Full record in
     `Claude_markdowns/2026-09-22a.md`.
  2. **It should never have been a target.** The altered band is an **output** a
     correct model predicts, not an input to calibrate against. Treating it as an
     anchor was planner error (mine): "we need a number to split delivery from
     cost" turned a validation check into a calibration knob. **The closure must
     instead be built with no free wall parameter** — gas temperature from
     conservation along the flow path, heat transfer from a standard annular-duct
     correlation — so the wall heat is a **prediction**, and the band (if a colour
     original ever appears) is a check on it.
  - **What survived is worth more than what failed: a measured Ø(z) for Meier's
    hole**, from rim texture and shadow rather than tint — ≈ 107 mm at 3–12 cm
    (his own deliberate underream), **95 ± 5 mm constant from 14 to 38 cm**, and a
    taper over the last 8 cm. **Meier's hole has a constant-diameter shaft; ours
    has no constant section anywhere** (145 mm at the top narrowing to 78 mm at
    the face). That is the scoring target, and it needs no energy accounting.
  - **Unresolved conflict, worth 20 minutes before it enters a criterion:** the
    bottom taper reads **60–70 mm** here and **78–80 mm** in the 09-18 pass of the
    same figure. That is the *freshly cut* diameter, so it decides whether our
    78 mm cut is too wide or too narrow — and it softens or strengthens D2m's
    "our hole can never exceed the tool's stance". The shaft also reads 95 mm in
    the photo against Table 8.2's 85 mm (known lighting bias on one rim), so the
    target is a **band, 85–95 mm**, not a number.
  - **Standing consequence:** Meier's data fix only delivery × (volume per joule);
    **nothing in Ch. 8 breaks that tie.** Score the closure on **shape**, report
    wall heat as a prediction, and never tune it to a target.

- **D2o — The wall closure, one packet with a shape score.** Only after D2m and
  D2n. **The physics to encode:** a surface element gets impinging-jet heating
  where the flame strikes it, and grazing-flow heating from **already-cooled gas**
  where it is a wall the exhaust slides past; the **enthalpy march continues along
  the wall with the gas it has left after the floor**, so the wall's gas
  temperature comes from conservation rather than a new free parameter.
  - **Orientation must come from a smoothed profile, not single-cell steps.** A
    cell-level discrete gradient at a surface edge is exactly where the last three
    mesh failures lived.
  - **The score is the diameter profile against depth**, specifically the **ratio
    of wall recession to floor recession, target ≤ 0.1** (Meier's hole is
    cylindrical and calibrated). The model is at roughly **0.3** today.
  - **The altered band from D2n is the check that the wall term was not simply
    zeroed** — passing the shape score with no wall heat at all is a failure, not
    a pass.
  - **Two meshes from the first day** (D2l: the cascade is invisible at a single
    mesh). Cost: a few days.

- **Deliberately NOT now (2026-09-22 decision):** gas temperature, cuttings
  absorption and melting. All three scale every flux roughly equally, so they move
  how fast the model drills but not the shape of the hole — and **shape is the
  only thing we can currently score against Meier without circularity**. They
  return once shape is right and the model can be scored on rate again. *(Partial
  caveat: debris accumulates preferentially in the annulus, so cuttings shielding
  is not perfectly uniform and could touch shape — second-order, and unmeasurable
  until the wall term exists.)*
- **Standing until D2o exists: no diameter, rate or energy number from the model
  may be quoted as a Meier comparison.**

- **AMR surface-band refinement — its own packet, reviewed before merge**
  (working tree, uncommitted, `bin/mmwspalling-3d-g++-amr`, notes in
  `Claude_markdowns/2026-09-22.md`; study `studies/amr_band/`). The finest level
  owns the surface (`SurfaceOwner`), a band follows each column's top cell,
  `RepairRemovalStateAfterRegrid` rebuilds the column-linear level set, and the
  single-level aborts were relaxed. Already gated on the hot-disk exclusion case:
  key-off byte identity on one run, whole-domain level-1 ≡ uniform 1 mm at 60 s
  (removal log byte-identical), and a 20 mm band within 1.4 % at 250 s. **Four
  gates still owed, in this order:**
  1. **the 2639-file reference sweep on the working tree** (D2l Gate 1 — on the
     critical path for *any* rebuild, AMR or not);
  2. **LI0 identity under whole-domain refinement** with side-face flux on,
     60 s, removal log byte-identical — this is a *wiring* test (which level owns
     each per-column gather), **not** a drift test: the feet/jet loop has diverged
     at 300–450 s every time, so a 60 s identity must never be quoted as "safe
     for Meier";
  3. **an LI0 band gate over a long window.** Note the existing D2k 1 mm leg is
     only usable to ≈ 350 s (it drilled out of the box by ~400 s), and it will be
     superseded by the D2l closure run — so this gate is best run *after* the
     closure lands, against a reference that will still be current;
  4. **a two-level energy ledger.** It aborts under AMR today, and it is the check
     that caught the E1 beam leak and proved D2k's budget exact. **Not a plain sum
     over levels:** with average-down the covered coarse cells duplicate the fine
     ones, so the sum must exclude covered regions via the level mask — a naive
     sum closes spuriously, which is worse than no ledger. No scored AMR run
     without it.
  - **Speed target before it goes near a scored run.** Measured wall ratio is
    **1.00** with 3.3× fewer cells (per-column gathers, two-level FillPatch, 339
    small level-1 boxes), so AMR currently buys correctness and no speed. State
    the target as wall time for a named job — the 0.5 mm arbiter leg (5–7 h
    uniform) or a taller Meier domain — not in the abstract. Known fixes: level-1
    `max_grid_size`/`blocking_factor`, `grid_eff`, coverage check at regrid only.
  - **Untested:** `max_level ≥ 2`, mechanics-on paths (dormant only because the
    D2 line runs `el.type = disable`), and **beam + AMR**, which has a silent
    trap: `P_now = owner ? … : 0`, so deposition below the band vanishes if the
    Beer-Lambert depth exceeds `surface.refine_depth`, with nothing tying the two
    together. Worth an abort while the code is fresh — it will bite the MMW anchor.
  - **AMR does not fix the mesh problem.** By the step-face mechanism refinement
    makes the cascade *worse*; AMR is a cost tool, not a remedy.

- **Overnight confirmation of D2k, once it passes:** `d2j0b_plane` C1 at
  **0.5 mm** with the key on, in the 0.15 m box (7.7 M cells, ≈ 5–7 h on four
  ranks). A third mesh point on the *fixed* configuration, not a gate on the
  2 mm / 1 mm result — the week's recurring failure was "converged at two meshes
  over a short window".
- **FALLBACK, only if D2k's arbiter fails — D2j-a/b/c, the continuous ablation
  front** (`Claude_markdowns/2026-09-21b.md`): `spall.removal = events |
  isotherm`, default-off and byte-identical across the 2639-file set; recession
  to the `T_abl` isotherm in Removal.H PASS 2; the 1-D ablation-column test; then
  the hot disk and the Meier pair. **No slope factor** (see Known Notes). It is
  the fallback **plus** side-face flux, never either alone: the front alone still
  leaves a staircase whose step faces are unheated.
- **WITHDRAWN by D2j-0b:** the D2h open-loop pair (it would only re-measure the
  same cascade without the feet), wall heating above the nozzle plane, and
  `spall.surface_normal` as a remedy.
- **D2g Stage 2 — Meier rescored under the D2f support rule** (blocked until a
  mesh gate passes; fresh pre-registration;
  0.70 m; T_nozzle sweep on a single basis; full-domain mesh pair). Every D2b–D2d Meier number
  is superseded by it. Tube-wall confinement of the gas (slots) stays a
  candidate if the mouth still fails.
- **Event location for the firing** (`2026-09-18a.md`): interpolate to the
  Sp = 1 crossing; removes the flake dt drift (32 ms+ fails today); re-run the
  dt ladder after. Default-off.
- **Outside the code (user):** ask Meier for the 160 L/h cooling-water ΔT (or
  the chamber design heat load); recompute adiabatic T for CH₄/air λ = 1.2
  with dissociation. **Commit the tree** (E1 … D2c) with the full sweep
  before the D2d campaign — pending user go-ahead.
- **D3 — Report rewrite** (after D2g): drop the 1.9×/P4 and flake/PSD drilling claims,
  state the 09-15b §1 energy accounting, and reposition as a validated
  moving-boundary thermoelastic ablation model with Weibull onset statistics.
  Must also state: the D2a/D2a2 supersessions and the 09-16 corrections to
  09-15b §1; `V0 = V_cell` switches the size effect off; the slope area factor
  and side-wall/exhaust heating omitted; face-form mesh sensitivity of the
  jet march; `T_ent` and entrainment as first-order unknowns; the D2b verdict
  as it came out.
- **Then MMW**: known-power anchor (Woskov/Oglesby), E2 α(T), R(T). MMW melt
  needs missing physics (S1 Q3), not resolution.
- **E2 — α(T) optical-depth closure**: per-cell `a_k` breaks column telescoping
  (closure not exact, can create energy). Accumulate τ = Σ α_k·dz down the
  column. Required before any real-α(T) MMW work.
- **Sub-cell *criterion* closure (flake size): RETIRED** with the flake claim.
  AMR + removal is not needed for ROP (S1 decision (i)).
- `input_2d_dev` drills ~100 of 120 mm in 56 s, and its header (~2.9 m/h,
  depth margin) is stale. Fix it before any longer Meier 2D runs.
- **Meier physics program (remaining)** — `tests/MMWSpalling/validation/meier/TODO.md`.
  P4 withdrawn (above). **P2 melt evacuation and P3 α(T) remain deprioritized**
  for Meier (both gated on the melt regime being reachable). P3 is now
  motivated from the MMW side (see E2). P5 subcritical crack growth; P6 wall
  stress anisotropy; P7 elastic-solve A/B; P8 melt-regime radiation.
- **Scope decision pending (user)**: commit to grain-scale (revive elastic
  solve P7, close CZM loop, sub-mm AMR) vs reposition as a "validated
  moving-boundary LEFM ablation model". See review Finding 3.
- **Step 19b (potential follow-on, only if Step 19 G4 FAILs)**: a
  lateral-roller-box + MMW beam variant that recovers 1-D-confinement
  σ_xx (Step 16 working geometry but with MMW as the heat source
  instead of `convective_flame`). Only spin up if Step 19's
  free-lateral `zlo_roller_321` test informally records "Sp fires but
  σ_xx too low for spall removal in the ramp" — that finding is the
  motivation. If Step 19 G4 informally PASSes, no Step 19b needed.
- Step 13: production AMR refinement criteria + adaptive timestep
  (eq. 40). **Deferred indefinitely** — the AMR-mechanics smoother
  limitation in `hu_spall_onset` is documented as permanent (five
  hypotheses refuted, see `AMR_FAILURE_ANALYSIS.md` §9–§10). Re-promote
  only if/when that limitation is resolved upstream.
- **Step 17b (Zhang/Oglesby full reproduction under sp_weibull):
  DEFERRED to Step 19+ scope.** Two specific obstacles for the
  packet that takes this on: (i) the Zhang/Oglesby reference
  reports T_surface_center(t), NOT ROP — long-plan §17a's "ALAMO
  ROP within 15%" mis-cites the source; (ii) `material.type =
  zhang` (custom plugin with T-dependent ρ, κ, Cp) doesn't map
  cleanly to `microstructure.phase0` (which stores scalars).
  Step 17b's planner needs a scope decision: re-validate
  T_surface under sp_weibull, find a different ROP-reporting
  reference, or add T-dependent property support to the
  microstructure phase struct.
- **Rossi 2018 R1 quantitative validation: DEFERRED to Step 19+ scope.**
  Three packets (16b, 16c, 16d) converged on the mesh-resolution
  limit predicted in `rossi-validation-diagnostic-design.md` §9
  — at the working 32³ in 50×50×50 mm geometry (dz = 1.5625 mm),
  no two-node sub-cell σ_xx reconstruction can produce a
  K_I = K_Ic crossing in Rossi's [100, 200] µm peak band. Revisit
  conditions in design note §11: (a) anisotropic z-refinement
  (`dz_top ≈ 20 µm` — needs MLMG conditioning improvement first),
  (b) a different reference experiment with cell-resolution-
  appropriate observables (Friedrich-Wong slow-rate damage is one
  candidate). The Rossi mechanics test
  (`tests/MMWSpalling/sp_rossi_damage_profile/`) stays as an
  ongoing P1-P4 mechanics regression; R1/R2/R3 are reported but
  not project-gated until one of the revisit paths is taken.

## Validation Story

- Analytical/unit tests validate individual numerical components.
- Zhang/Oglesby validates the MMW source and thermal model.
- Hu validates the prescribed-temperature thermal-spallation chain on
  the `damage_law` branch: 6d conduction -> 7b/7d thermoelastic stress
  -> 8b breakage indicator -> 11b/14 onset/removal.
- Kant 2017 Central Aare granite validates the LEFM Sp onset criterion
  on the `sp_weibull` branch: 12c integrator unit -> 15b `p = 0` onset
  + Eq. 15 cross-check + Weibull patchiness (done) -> 15c labeller +
  sensitivity sweep (active) -> 18 confining-pressure sweep.
- Rossi 2018: a prescribed-T mechanics regression of the Sp depth scan
  (isotherm depth; Step 16, `sp_weibull`). It is **not** a flake-size
  validation for self-heated drilling (retired 2026-09-15, S1).
- Meier 2017 pilot: the ROP and hole-diameter validation, once the D2
  face-source + burner-descent model exists. Until then only volume rate and
  ΔT_fire are meaningful.
- The final depth scan is exploratory; it has no direct experimental
  validation.

## Context Files

- `CLAUDE.md`: startup router, project rules, commands, and ALAMO gotchas.
- `ACTIVE_STEP.md`: the only detailed task packet Claude should read by default.
- `ARCHIVE_DONE.md`: completed-step history. Read only for debugging regressions
  or when Codex is updating project context.
- `in-main-tex-you-will-quizzical-treasure.md`: full plan and historical
  checklist. Do not read by default.
- `main.tex`: full physics proposal and equations. Read only relevant sections.
- Validation `.txt` files: source data and pass criteria for named validation
  steps only.

## Known Stale/Important Notes

- **Idle clock (`pinned_idle_cycles`) is mesh-dependent (D2h, 2026-09-20).** Its
  timeout is measured in cell-heating times (t_cell ∝ dz), and the unpinned face
  form is self-limiting, so at finer meshes edge columns park sub-critical and a
  silent front walks inward. Set the key per input (1e9 = off); **changing the
  source default rebases every pinned-closure result since C1** and needs its own
  packet with the 2639-file reference set re-checked.

- **D2e retroactive review (2026-09-19, accepted) — carry to the commit.** The
  2639-file witness set is single-level MMWSpalling only; the Perf edits touch
  shared code (`BC::Constant`, `Numeric/Material/Table.H`) and an AMR regrid
  path, so the deferred full sweep must include the AMR microstructure test and
  other `BC::Constant` users. Nits: the `Table.H` seed can in principle pick a
  bracket one to the right at an ulp tie (step `lo` back while `m_H[lo] > H`);
  `Constant.cpp` `Domain().contains(box)` mixes index types (debug-build
  `sameType` assert on node-centred FABs).

- **How Meier ran the pilot test vs the model (2026-09-18, full Ch. 8
  reread; `Claude_markdowns/2026-09-18c.md`).** Decisions for D2g, before any
  scoring:
  1. **The mouth was inside a sealed, water-cooled wellhead** (Fig. 8.1,
     p. 213), glued to the rock, with exhaust to a cyclone and a bore of
     roughly twice the burner Ø (visual estimate; no dimension given). In
     phase 1 the gas over the surface was confined, cooled exhaust, not
     293 K room air. D2c's free-surface ambient dilution, which D2d credits
     with about half the mouth fix, rests on an open surface that did not
     exist. The wellhead also caps any funnel. **Confirmed by the user against
     independent thesis notes, 2026-09-20: this is a gas-side finding, not
     only a mouth-scoring one — it is the same physics as the tube-wall
     confinement, so Stage 2 must carry at minimum a free-surface-off variant
     as a sensitivity rather than treat dilution as settled.**
  2. **Two drilling phases** (647 s + 736 s, p. 219). Phase 2 restarted in a
     cooled, part-drilled hole (about 0.25 m, estimate) with the hole in
     rough vacuum and room air drawn in through the burner–wellhead gap. The
     model is one continuous run at 1 atm.
  3. **The thesis's ROP figures disagree with each other:** 0.5 m / 1383 s =
     **1.30 m/h**, against the text's 1.5. TSE 15.4 J/mm³ confirms 1383 s.
     Quote 1.30 as the time average, and state the discrepancy.
  4. **The measured ROP is an operator choice at or below the contact
     limit** (crane advance, hook load, clearance set by advance rate,
     pp. 214–223). The model's contact rule is an upper bound. Read the ROP
     one-sided: a model ROP below the measurement is a miss, while a
     somewhat higher ROP with a somewhat narrower hole is consistent.
  5. **The cooling-water ΔT was instrumented** (thermocouples up- and
     downstream, p. 208) but not reported, so the request to Meier is for
     recorded data. Also ask for the wellhead bore and height, the depth at
     the end of phase 1, and what the 1.5 m/h is.
  - **D2f wording corrected 2026-09-18** (the gate is unaffected):
    - "jolts clear isolated high rock" and "gas reached under the rim" are
      our inferences;
    - "wider than the tube at every depth" is overstated (Tab. 8.2: 85 mm;
      Fig. 8.8 visible width 80 → 78 mm near the exit, a lower bound).
- **D2f planning decisions (2026-09-18, revised same day).** A shielded
  "ledge" rule (rock under the rim gets no gas; broken off when undercut) was
  planned and **withdrawn before implementation**: Meier's hole is wider than
  the tube at every depth and he calls the tool contactless, so the rock under
  the feet was removed thermally. The support is the contact rule with the rock
  heated and no clip; the per-pad quantile is read as the operator-jolt
  allowance and must be shown not to control the rate. The model never lifts
  the burner (Meier's holds and jolts shape the hole locally).
- **D2e planning decisions (2026-09-18).** 09-18b's P-b as written
  (`jet_T_ent = 1550` under `fixed`) would also move the march floor and the
  free-surface ambient, so D2e adds `rec_fixed` (T_mix only). P-c
  (`foot_rule = mean`) must run with clearance off (parse abort), so it is read
  against P-a. 8 ms at 1 mm is extrapolated (K/step ∝ dt/dz); B0 checks it
  against D2d's 4 ms R7. The Perf edits are unreviewed until D2e Goal 0.
  Regression sweep 09-17: `heterogeneous_kappa`, `amr_microstructure_regrid`
  fail under `mpirun -np 4` (k_eff diagnostic not reduced per rank; one-line
  `ReduceLongSum` fix when that file is next open). `pkill -f <run name>` kills
  shell waiters too. Advisor's Hutchinson–Suo buckling release (09-18) is the
  one import worth considering later, not now.

- **D2d planning decisions (2026-09-17, planner, before any run).** The
  steady-window mesh criterion is replaced by a matched transient window on a
  0.40 m full-domain pair (a steady 1 mm run needs ≈ 24 h). Matrix domain
  0.70 m (hand needs 0.65 m, runs 13 % fast). Min Ø scored to the feet depth
  (D2c's nozzle-plane range missed the 50 mm feet zone). T sweep lower end
  (1.3 m/h) is below 1600 K and reported open-ended. Run 3 (n = 0.5) is hand
  only. "D2b closures + n = 1" keeps exponential/count/clearance.

- **D2c planning decisions (2026-09-17, user-reviewed).**
  - Feet-depth blending of the stagnation gas is **rejected**: the jet runs
    inside the Ø 56 collar from t = 0, and a blend starting at ambient
    reproduces the D2b 293 K stall.
  - The 12 D clamp lives in the Python `h_expr` string, not in C++.
  - Nozzle (γ 1.3, ṁ measured, 1900 K): p0 5.95 bar, M_e 1.86, p_e 0.97 bar
    (perfectly expanded), u_e ≈ 1277 m/s, J ≈ 19.3 N; D_e 3.5–3.6 mm in 293 K
    air, 8.0 mm in 1500 K exhaust. 8·D_e = 29 mm ambient (not 37), 64 mm in
    exhaust. The quarter domain must use the full-jet D_e (mdot/4 would halve
    it).
  - Fig. 8.8 is digitised (`validation/meier/meier_fig8_8_hole_profile.csv`,
    visible width = lower bound, ±5 mm): 96 mm at 10 mm, 88 at 50, 87 to
    125, 84 at 300, 78 at the exit. A mild collar (~10 mm) is real.
  - Volume is scored on the whole excavation (water filling); the 09-16
    inside-Ø 93 addendum is withdrawn. "1436 K" is not a gas temperature.

- **P1 melt-aware spallation is DONE (default-off `spall.melt_aware`).** The lever
  is the DRIVING STRESS, not K_Ic: σ_xx scaled by `f(Λ_L)=max(0,1−Λ_L/Λ_crit)`
  (Λ_crit=0.4) in `Removal.H`'s `k_i_at_zcrack` (≡ `E_eff=f·E`). Do NOT reduce
  K_Ic in melt (backwards). KEY finding for all melt work: **the spall thermostat
  pins the surface at the ~850 K spall onset, so melt is unreachable in real
  drilling without α(T)** — melt only appears where removal is throttled (domain
  bottom, confounding the Meier melt demo) or via an imposed IC (the clean
  `unit/sp_melt_skeleton/` A/B test). `prescribed_T` applies at the FACE and does
  NOT follow surface recession. Test the removal split (`removed`), not the
  melt-blind `Sp_field`. See `mmwspalling-melt-aware-p1` memory + ARCHIVE_DONE.
- **Beam/void gap energy leak — FIXED by E1.** Every beam+removal result from
  before E1 deposited only 56–70% of beam energy: old Meier ROPs, the
  "1.9× too fast" framing, the N1 numbers, and memory's "implied efficiency
  ~17%" are all biased LOW or confounded. Invariant: solid ⇔ φ ≥ 0, so
  φ_top ∈ [0, dz]. **Round-off ties φ_top = dz are common**, so any new
  "top solid cell" logic must use the `removed` mask (`ColumnTopSolid`), never
  `floor(φ/dz)` or a φ-window. **The criterion and the beam now share one depth
  origin**: the K_I scan's `z_crack = 0` is the top face of `k_top` (E1
  takeaway 3), so the older "criterion measures from the level set" note was
  wrong. Deposition closure is exact only for constant α + uniform profile
  (see E2). `Claude_markdowns/scripts/phi_gap.py` still models the pre-E1
  formula; `energy_audit.py` is still valid. `unit/beam_void_closure` T3 (lower
  half-cell hit) is fragile: re-check it if that input changes. The
  `AMREX_DEBUG` Abort branch has never been compiled.
- **`surface_mf` is φ-windowed** (`UpdateSurfaceMaskFromPhi`, φ ∈ [0, dz),
  refreshed each step at Removal.H ~1295). At a φ_top = dz tie a column has NO
  surface cell, so radiation/convection losses drop out for that step.
  **Measured by S1 as negligible**: on Meier 2D, 1 of 280 columns for ~1.1 s,
  7e-5 of column-steps, ~0.03 J of 77 kJ. Use `surface.follow_mask = 1` whenever
  a source must track the receding surface.
- **The S1 criterion finding (2026-09-15).**
  - `sp_weibull` + `spall.sample = top_cell` + `local_thermoelastic` fires at
    a fixed **cell-average** T: 656 K at a = 20 µm, p = 1 MPa; ≈ 812 K for the
    Meier Weibull draw (to be confirmed by A1).
  - **h_col ≈ dz/2 − a_f is a removal increment, not a flake.** Two firings
    empty one cell. Flake/PSD claims are retired for self-heated drilling.
    Rossi (prescribed_T) keeps its isotherm-depth reading.
  - At fixed flux, ROP = q/(ρ·Cp·ΔT_fire) at any dz.
  - Robin at 2 mm is not converged (cell +42 %, face −82 %). The converged
    1-D answer is the pinned surface q = h(T_gas − T_fire) − εσ(T_fire⁴ −
    T_a⁴).
  - **Quote Robin as a bracket.** The finest-dz data bracket is 14.0–18.6 m/h
    at (1e4, 1000 K); "15–16 m/h" is an extrapolation.
  - Do not reconstruct a sub-cell T from the 2 mm cell average for the Robin
    flux (it gives ~3000 K).
- **Overshoot time-step rule for removal runs** (binding over the conduction
  limit): ΔT_step = q·dt/(ρ·Cp·dz) ≤ 5 K. That means dt ≤ 11 ms at 2 mm and
  2 MW/m²; ≤ 1.4 ms at 2 mm and 15 MW/m²; ≤ 0.7 ms at 0.125 mm and 2 MW/m².
- **`spall.h_col_events_csv` is edge-triggered.** It logs only the first firing
  per new top cell, so never sum it for recession (≈ 2× undercount). Rossi's
  test depends on these semantics. Use `spall.removal_events_csv` (A1) for
  recession and T_fire: one row per (step, column) removal. Key on `regime`,
  because `h_scan` can be nonzero on regime 2/3 rows.
- **`weibull.enabled = 0` + `spall.enabled` segfault — FIXED by A1.** The
  scalar path's per_face gate is `Sp_field(k_top) >= 1`. `connected_cluster`
  still needs Weibull. The S1 degenerate-draw workaround is now optional.
- **Weibull `V0` (user decision 2026-09-15): use `weibull.V0 = V_cell` in all
  C/D inputs.**
  - With the default 1e-9, T_fire drifts ~+10 K per halving and never
    converges (S1b).
  - The default is left unchanged, so Kant/Rossi/Meier existing inputs are
    unchanged.
  - The IG vol_factor applies only to IG cells (GB = per-facet).
  - The Meier 2D 812 K / 6.3 m/h numbers are V0 = 1e-9; with V_cell, 2 mm
    fires at ≈ 820 K.
- **Pinned Robin closure (C1) caveats.**
  - The closed-form match is mostly a consistency check; the independent
    evidence is that the result lies inside S1's fine-mesh bracket.
  - A sharp uniform patch edge leaves a never-firing rim column pair, so score
    non-rim columns or use a smooth profile.
  - Time to first firing is face-form and mesh-dependent: score by windows and
    `removal_events.csv` `h_applied` fits, never by end depth or plotfile
    staircases.
  - Pin state is not checkpointed.
  - `pinned_idle_cycles` (default 2) is not a knob: 1-D refire gap
    1.004·t_cell, 2-D up to 1.32·t_cell.
- **Burner drawings (planner, 2026-09-16) — the D2b blocker is RESOLVED.**
  Source: `thierrymeier.ch/documents/drillBabydrill_VT5_R4.pdf`
  ("Zusammenstellung Drill baby drill VT5", CATIA V5, 19.11.2015, 17 A2
  sheets, 797 mm long, ca. 14 kg). **Not in the repo** — re-download to
  re-check.
  - Sheet 15/17 **"Abstand Halter F VT5"** = the feet: **OD Ø 80 mm**, bore
    Ø 45 (+0.2/+0.05), plus Ø 56 and a Ø 68 bolt circle with Ø 6.2 holes;
    overall length **58 mm** with an **8 mm** collar, so it protrudes
    **50 mm**; the isometric shows **3 slots**, so the bearing surface is an
    interrupted annulus (a tripod), not a continuous ring.
  - Sheet 14/17 **"Flame jet nozzle D VT5"**: flange **OD Ø 80 mm**, Ø 68
    bolt circle, Ø 6.2 holes at 120° (it bolts to the holder); central bore
    **Ø 7.5 mm** at the outlet; converging contour R12/R20/R10 from a Ø 46
    bore, Ø 56 f7 spigot.
  - Sheet 1/17 assembly, section A-A: confirms the Ø 7.5 bore and a **50 mm**
    nozzle-exit-to-foot-tip dimension. The mantle above the nozzle is a
    67 × 2.5 tube (Ø 67), narrower than the Ø 80 drilling end.
  - **Derived:** burner Ø at the drilling end **80 mm** (r = 0.040 m); foot
    annulus **[0.028, 0.040] m** — section N-N, re-read at higher resolution,
    shows the feet are a **tube OD Ø 80 / ID Ø 56** running 50 mm below an
    8 mm flange (bore Ø 45), with 3 slots at the bolt positions (the earlier
    [0.034, 0.040] "nominal" misread the Ø 68 bolt circle). The jet expands
    inside this collar for its full 50 mm and exhausts through the slots;
    `foot_standoff` **0.050 m**; `D = 7.5e-3 m`, so `SOD/D = 6.67` (not 7).
    Meier's Ø 85–93 mm hole is 5–13 mm larger than the burner, consistent
    with thesis p. 223.
  - **Two discrepancies, do not silently resolve.** (1) The drawing's
    **D = 7.5 mm** vs the reference `.md`'s **7.1 mm** (thesis text): the
    drawing wins; consequences are Re ≈ 5.66e4, h_Martin ~8 % lower
    (h ∝ D^−1.5), potential core 5D = 37.5 mm. D2a's `jet.py` used 7.1 mm.
    (2) The drawing's flows (50 + 5 m³ₙ/h ≈ 62–64 + 3.4–3.6 kg/h) are
    **design** values, above the thesis's **measured** 52 + 2.47 kg/h. Use
    the measured `ṁ = 0.01513 kg/s`.
- **D2a's wall jet violated the jet energy budget — D2a2 fixes it**
  (reviewer, `Claude_markdowns/2026-09-16.md`; planner-verified arithmetic).
  - D2a's study `T_gas` depended on **stand-off only**, with `s` clamped to
    [2D, 12D]. So rock far off-axis but close to the nozzle plane received
    potential-core gas. The reported J-5/J-10 `R_h` (149/170 mm) is
    reproduced exactly by the closed form **at the clamp** — the `RESULTS.md`
    explanation ("h at r = 100 mm is still ≈ 2 kW/m²K") is not the mechanism.
  - With `h ∝ 1/r` and a radially uniform `(T_gas − T_s)`, `∫q·2πr dr` grows
    **linearly** with radius, so face power is unbounded. D2a's 20.9–29.0 kW
    exceeds the cap `ṁ·cp·(T_nozzle − T_s)` ≈ **20.4 kW** at 1900 K and
    **11.6 kW** at the 1436 K chamber TC. Meier's rock-side removal power is
    **2.83 kW**.
  - **Both readings are superseded.** D2a's machinery, the "centre ROP is an
    identity under a prescribed feed" result, the steady-bowl timing, the
    closed C1 nits and the parser-bug find all stand.
  - Also corrects 09-15b §1: "low-bracket stagnation flux fits the volume
    rate within 35 % untuned" was a **stagnation-point** statement; integrated
    over the Ø 85 face the bracket anchors give 2.4–2.5× the data and Martin
    1.4×.
  - **Report face power against the jet cap in kW, never as a % of 38 kW.**
- **D2a2 lessons that constrain D2b and later jet work (2026-09-16).**
  - **Keep the stagnation decay.** Without it the centre runs away from any
    feed. Baseline treatment for D2b: `jet_stagnation = decay` +
    `jet_entrained_mass = 1` (energy-conserving; nozzle mass discards 12–17 kW
    of the 30.4 kW nozzle budget as `jet_P_decay`).
  - **The centre leads the nozzle** (s_c 60–90 mm, T_stag 1000–1270 K at
    1900 K). Hand predictions must use the centre's equilibrium stand-off
    (`walljet.py` `s_equilibrium`), never s_c = SOD.
  - **Planner hand-table errors, recorded.** D2a2's first table used a
    non-computable criterion (T_gas decays exponentially toward T_s and never
    reaches T_fire) and a linear in-core depletion; the revised table's decay
    rows were computed at **D = 7.1 mm** (`jet.py`), not 7.5. Hand numbers go
    through `walljet.py` (D = 7.5 mm) and are written by the implementer
    before runs.
  - **Cold entrainment lowers flux at every radius** (reviewer, integrated):
    neglecting it is an upper bound on flux and on removal reach.
  - `jet_mdot` is the modelled sector's share (quarter domain = mdot/4).
  - **Wall treatment is open:** `input_drilling`'s
    `spall.flake_coherence_length = 0.004` caps walls at 63–65° and capped
    cells overheat under the pinned flux; `spall.surface_normal = 1` recedes by
    h/cos θ with no matching flux on horizontal area. D2b decides by a
    pre-stated energy-consistency rule.
  - **Most of the march is face-form uptake by unfired rock** (pinned share
    0.01–0.12), so T_gas(r) and the ring flux carry the 2 mm face-form mesh
    sensitivity.
  - Under a prescribed feed the flank freezes once the nozzle plane passes
    it; read prescribed-feed runs before that time. The feet rule removes it.
  - `jet_dr = dx` is first-order in the march (16 % excess-T error at 2 mm in
    a strip) but moves the quarter-disk ROP only 0.4–1.4 % per halving.
  - `patch_P_robin − jet_P_face` is not a clean overdraw measure; a direct
    counter is deferred.
- **Process (planner decision 2026-09-16, from D2a2 takeaway 12).** Per step:
  targeted bit-identity witnesses for the touched paths, with the full
  26k-file sweep once at commit; no 2-D slab parts for jet closures (the slab
  needs an invented jet share); a few 3-D runs rather than full grids; sparse
  `amr.plot_int` (D2a2 wrote 11 GB). Implementers re-read `ACTIVE_STEP.md`
  before long runs (D2a2 was revised mid-phase).
- **AMReX 25.12 parser bug (in `ext/`, not touched; found by D2a).**
  `parser_ast_optimize` rewrites `f / F2(x, number)` as `f * F2(x, −number)`
  for **every** two-argument function, not just `pow`. So `1527/max(r, 0.5)`
  evaluates to 0, silently. **Numeric arguments must come first.** D2a's first
  probe and Part 1 set ran with the flame effectively off and were discarded.
  D2a2 adds an in-ALAMO regression guard in `unit/jet_enthalpy`.
- **`weibull.V0 = V_cell` switches the size effect off, not just its mesh
  dependence** (reviewer, 2026-09-16). `vol_factor ≡ 1`, so `a0_ig` becomes
  "the characteristic flaw of a cell". The right choice, but **D3 must not
  describe the model as size-scaling.**
- **Jet h magnitude conflict (planner, 2026-09-15; D2a2 narrowed it, D2b settles).** Under the D2a2 closure every anchor drilled the r = 40 mm ring at only 0.5–1.3 m/h, so the Martin-vs-bracket difference is now a factor ~2 in ring ROP, not in face power. Martin
  (1977) with Meier's ṁ = 0.0151 kg/s, D = 7.1 mm, Re ≈ 5×10⁴ gives
  h ≈ 1.5 kW/m²K at 2.5 D, SOD 7 D. The reference bracket is 5–30 kW/m²K, and
  Ch. 7's 10 kW/m²K comes from Potter Drilling, not measured. Correlation
  constants were quoted from memory; D2a must verify them against the source.
- **AMReX rank-count identity:** Cell_D is written one file per rank. Compare
  the rank-ordered concatenation or the field arrays, not a file-by-file `cmp`.
- **Energy ledger (`energy_ledger.enabled`)** is exact only with Neumann-0
  outer BCs, and its state is not checkpointed (it re-baselines on restart).
- **Kant p = 0 onset now prints 446.1 °C** (Eq. 15 445.6 °C) in
  `sp_kant_pressure_sweep` (reported by E1), both pre and post E1 on the
  current tree, vs the 461.3 °C quoted in older notes. The cause has not been
  identified (it predates E1). Treat the live test output as authoritative
  until reconciled.
- **P4 cuttings/debris shielding — WITHDRAWN 2026-09-15 (design kept for
  reference; premise invalid, Σ₀ uncalibratable on Meier).** Default-off `shield.*`.
  Per-column debris areal mass Σ(i,j) [kg/m²]: `Σ ← Σ·exp(−dt/τ) + ρ·h_spall`
  each step (SPALL only, not vapor — vaporized rock leaves as gas); the beam is
  attenuated `Pi → Pi·exp(−Σ/Σ₀)` in `AdvanceMicrostructure` only (the sole
  surface-following/removal path). Removal runs BEFORE the beam in `Advance()`
  (747 vs 1229), so the beam reads fresh debris (stable — `h_spall` derives from
  the prior-step temperature, no implicit loop). Σ stored broadcast down each
  column so the per-cell beam read is z-decomposition-safe; `h_s_cand_col` is
  already MPI-gathered. Single-level only (max_level=0 dev harness); AMR debris
  repair deferred. Steady η is set by τ/Σ₀ (degenerate at steady state; τ also
  sets the smoothing timescale — pick τ physical ~1 s, calibrate Σ₀).
- **Refactor track (R-series) is decoupling `src/Integrator/MMWSpalling.H`
  (~4377 lines) into per-concern pieces WITHOUT changing its virtual
  multiple-inheritance graph** (CLAUDE.md Lessons #1 — base-class edits caused
  AMR vtable crashes before). Planned decomposition by field-ownership cluster:
  Thermal/enthalpy, Microstructure, ContinuousDamage (`damage_law`), Cohesive,
  LefmSp (`sp_weibull`), Removal. **R1 (active)** extracts the ~957-line
  `UpdateRemovalAfterCohesive` into `src/Integrator/MMWSpalling/Removal.H`,
  `#include`d *inside the class body* — the method stays an inline member, pure
  code-motion, byte-identical regressions. The partial-header trick is safe
  because the build (Makefile) compiles objects only from `*.cpp` (mindepth ≥2)
  and top-level `*.cc`; `*.H` files are include-only and never standalone-
  compiled, and incremental rebuilds still trigger via the generated `.d` deps.
  Such a fragment must contain NO `#include`/`namespace`/`class` (it is pasted
  into `class MMWSpalling { ... }`). Data-member / helper relocation is
  deferred to R2+. New pieces go under `src/Integrator/MMWSpalling/` (mirrors the
  existing `src/Integrator/Base/` subdir convention).
- The long plan still contains an old Zhang/Oglesby note with `omega_0 = 0.20 m`.
  The corrected validation reference says `omega_0 = 0.02 m`.
- Step numbering contains inserted/moved validation steps such as 4b, 6d, 7b,
  7c, 7d, 8b, and 11b. Keep those labels because tests and notes refer to them.
- Homogeneous Hu reproductions should stay homogeneous. Voronoi heterogeneity is
  a model-extension study for Hu, not the direct reproduction.
- Step 8 changed AMR microstructure repair: passive `D` is still re-evaluated
  from `damage.ic.*` when `damage.enabled = 0`, but evolved `D` is no longer
  reset from ICs when `damage.enabled = 1`; repair refreshes derived fields
  such as `kappa_eff`.
- Step 8b is complete and must remain an indicator validation with
  `damage.enabled = 0`: it validates `f_b = sigma_v / sigma_s` from
  homogeneous thermoelastic stress, not the DP damage law.
- Step 8b fixed a mechanics coupling issue: `surface_patch` prescribed
  temperature must feed top-patch nodal eigenstrain, not only the cell-centred
  heat flux. Preserve that coupling for future surface-T / moving-surface work.
- Step 8b passes Hu surface `f_b` targets but relaxes L6 deep `f_b` to `<0.25`
  because the bottom clamp inflates deeper stresses; strict `<0.1` remains
  blocked on a stable pure-roller elastic solve.
- Quartz alpha-beta transformation strain is implemented only on the
  microstructure mechanics path, default-off per phase via
  `microstructure.phaseN.alpha_beta.*`, and verified by
  `tests/MMWSpalling/alpha_beta_transition/`.
- Step 9 CZM is default-off under `cohesive.enabled`. Its fields
  `gb_delta`, `gb_delta_max`, `gb_traction`, and `gb_cohesive_damage` are
  registered only when enabled and key on `is_grain_boundary`, not `is_gb`.
  Mechanics traction feedback is not active yet.
- Step 11 removal is default-off under `spall.enabled || vapor.enabled`.
  Shared fields `phi`, `removed`, `regime`, and `RoP` register whenever either
  mode is enabled; per-mode diagnostics register only for their active mode.
  Per-column removal is MPI-safe for arbitrary BoxArray partitioning, including
  z-decomposed columns. AMR moving-surface validation is still deferred.
- Step 11 void handling exists only on the microstructure heat path; the
  constant-alpha, global-material explicit, and global-material implicit paths
  do not yet honor `removed`.
- For Hu onset/removal work that needs damage or spall fields, use homogeneous
  one-phase `microstructure.mode = expression` inputs rather than Voronoi. The
  direct Hu reproduction stays homogeneous; heterogeneity is a later extension.
- Step 11b is a calibrated v1 onset validation, not a fully predictive damage
  closure. It uses single-level, spall removal off, `damage.stiffness_floor
  = 1.0`, and material-specific `A_D` values to match Hu onset/LRST. Sandstone
  still develops a deep DP damage zone; the test reports this as a warning.
  Re-enable AMR and strict shallow-depth failure only after the mechanics/regrid
  and DP localization issues are resolved. As of Step 12, `hu_spall_onset` runs
  at `amr.n_cell = 32 32 32` (was `40 40 40`) — the 32³ partition is a clean
  16+16 box split under `max_grid_size = 16` and dodges the 40³ intermittent
  SIGBUS that BUS_ERROR_ANALYSIS.md documents. The sandstone tolerance window
  was retargeted from the Hu-average onset (76.7 s) to the Hu-primary onset
  (89.0 s), keeping ±20%. The live 32³ result is granite 37.6 s / sandstone
  93.0 s, both inside the tolerance windows; ARCHIVE_DONE.md keeps the
  historical 40³ numbers (37.0 s / 87.0 s) with an explicit post-Step-12
  update note.
- AMR multi-level mechanics for `hu_spall_onset` is a known limitation, not a
  bug to chase. Five hypotheses (Jacobi smoother bandwidth, bottom-solver type,
  CF interface placement via `n_error_buf`, RBM null-space deflation, F0
  cross-level inconsistency) were refuted with quantitative numbers. See
  `tests/MMWSpalling/hu_spall_onset/AMR_FAILURE_ANALYSIS.md` §9–§10 and
  `AMR_SMOOTHER_LIMITATION.md` before re-opening; surviving candidates are
  per-level coefficient consistency and vector-elastic multigrid machinery vs
  `MLNodeLaplacian`. The in-tree diagnostic toolkit (`DiagnoseAMRResidual`,
  `ApplyAMRF0Consistency`, `CheckAMRMechanicsFieldsFinite` and the
  `Check*FieldFinite` family, plus the `mmw.debug_amr_mechanics` /
  `ALAMO_AMR_DIAG` / `ALAMO_LEV1_F0_FROM_COARSE` gates) was removed during
  post-Step-18 cleanup — see [git history of `src/Integrator/MMWSpalling.H`]
  for the prior implementation if reopening this work.
- `revised-model-approach.md` (Stages 1–5 of the LEFM Sp + Weibull
  re-architecture) is in the repository root and is the authoritative spec
  for the `sp_weibull` branch. The long plan §12–§18 is the implementation
  staging; cross-check both when planning sp_weibull work.
- Step 12 ships the LEFM Sp + K_Ic infrastructure (toggle, Tada-weight
  K_I integral via trig substitution, Nasseri K_Ic table, optional
  prescribed σ_xx / T_surface parser overrides for verification) under
  `spallation.model = sp_weibull`. Default `damage_law` preserves
  Step 8–11b behavior byte-for-byte. EOS work is dropped indefinitely;
  no `eos.type`, no `src/Numeric/EOS/*`, no `tests/MMWSpalling/eos/`.
  DP "logged-only under sp_weibull" (revised-model-approach §10 #4) is
  NOT implemented yet — `damage.enabled = 1` is a hard parser error
  alongside `sp_weibull`.
- The K_I sign convention is resolved as of Step 15b (Sp mechanics): the
  long-plan §12b `-sigma_xx` (compression-positive) integrand is now
  applied in `UpdateSpAfterMechanics` on both the prescribed and physical
  paths, the physical path samples `stress_mf` depth-resolved
  (`k_d = top_k - floor(depth/dz)`, clamped), and per-cell `a_f` comes
  from `flaw_a_mf` when `weibull.enabled`. `sp_onset_kant_closed_form`
  (§12c) was the only existing test changed (input now prescribes a
  compressive `sigma_xx`, Python integrates `-sigma`, explicit `Sp > 0`
  assertion). All other regressions stayed byte-identical because
  `damage_law` early-returns from `UpdateSpAfterMechanics`. The long-plan
  §15b was split into a mechanics half and a Kant-onset-validation half
  per user direction; both halves are now done.
- **Kant onset closed form is Kant Eq. 15, NOT Eq. 14 `Sp_red`.** Eq. 14
  `Sp_red = Sp_rock·Sp_fluid + Sp_conf` is an `h_fl`-dependent
  rock-classification number, not an onset temperature. The onset
  equation is Kant Eq. 15
  `Θ_onset = K_Ic·(1−ν)/(2·1.763·√(a/π)·E·α)` at `p=0` —
  `h_fl`-independent and `λ`-independent, a material threshold. `h_fl`
  sets only the onset *time*. Step 18 (confining-pressure sweep) must
  use Eq. 15 with the `−p·ν/(1−ν)` confining term, NOT Eq. 14.
  `kant-2017-central-aare-validation-reference.md` §5 was corrected
  during Step 15b (Kant onset) accordingly.
- Kant's Table-2 390–560 °C band was computed with **frozen room-T
  `K_Ic⁰ = 1.5`** (Kant Table 1 / §4.3), not `K_Ic(T)`. The
  `sp_kant_onset/input_verification` runs use frozen `K_Ic⁰` and are
  the ones scored against the Kant band; `input_model` uses the
  Nasseri-`×1.05` `K_Ic(T)` table and is reported as the production
  model's honest prediction, not banded. Honest agreement at Kant
  onset is ~3.2–3.5% Eq.15-vs-FEM; the FEM `σ_xx` residual is
  ~0.96–0.97× the thin-film ideal (O(dz/thermal-layer) discretization,
  approaches 1 as the thermal layer thickens with slower `h_fl`) —
  hence the 8% target-3 tolerance.
- The **lateral roller box** (`el.bc.type = constant`, all 26 regions:
  `u_x = 0` x-faces, `u_y = 0` y-faces, `u_z = 0` zlo, +z free; all 12
  edges + 8 corners specified consistently) + **full-face**
  `convective_flame` heating is the canonical Kant-style 1-D-confinement
  setup — `ε_xx = ε_yy = 0` exactly, `σ_xx = -EαΔT/(1-ν)` per cell,
  MLMG well-conditioned, coarse cubic 1 mm mesh sufficient. Reuse this
  layout for Step 18 (confining-pressure sweep). The patch + cold-rim
  layout was abandoned during Step 15b Kant validation because the
  sub-cell thermal layer under-resolves `σ_xx` (FEM only ~0.45× ideal)
  and fixing it requires anisotropic fine-`dz` cells that the elastic
  MLMG solver does not converge on (see next note).
- **The elastic MLMG solver fails on anisotropic cells** — including a
  modest 2:1 aspect ratio. Keep mechanics-active meshes cubic. Sub-key
  trap with `el.bc.type = constant`: BC-type goes in the prefix, so the
  sub-keys are `el.bc.constant.type.<region>` (one per all 26 regions).
- `weibull.{enabled,a0_gb,a0_ig,m,seed}` (Step 15b mechanics) are
  default-off, parsed only under `spallation.model = sp_weibull`.
  `weibull.a0_gb`/`a0_ig` are the Weibull *scale* parameter, not the
  literal mean — sample mean is `a0·Γ(1+1/m)` (≈ `0.965·a0` for `m=15`).
  Kant validation used `a0_gb = 20 µm` as the *scale* (the ≈3.5% scale
  vs mean delta is well inside the FEM residual). `flaw_a_mf` is
  cell-centred (the long-plan §15b "face-centred" wording was deviated
  from, with rationale: ALAMO microstructure topology is cell-centred
  and Sp is per surface cell). The per-cell RNG is a splitmix64 hash of
  the GLOBAL cell index — bit-identical across MPI rank counts.
  `InitializeFlaws` runs at the tail of `InitializeMicrostructure` so it
  also fires on AMR regrid; because it is keyed on the global cell
  index it resamples deterministically. **AMR regrid-repair of
  `flaw_a_mf` was landed in §15c** (`sp_weibull_unit/input_regrid`
  asserts bit-identity across regrid events; max level observed = 1).
- `sp_weibull` material removal is **wired by Step 16 (mechanics
  done)**. `UpdateRemovalAfterCohesive` PASS 2 spall arm is branched
  on `spallation_model`: under `SpWeibull`, the gate is
  `Sp_cluster_id_mf(i, j, k_top) > 0` (consumes §15c labeller) and
  `h_col` is the deepest `z_crack` from a per-column K_I-vs-K_Ic scan
  with `a_f` fixed at the surface cell's `flaw_a`. New diagnostic
  `h_spall_field_mf` (cell-centred, registered only under
  `sp_weibull && spall.enabled`) records the **running max** of
  `h_col` per column at the **immovable initial top-z slice**
  (`top_k_dom`), NOT reset per step — so a single last-plotfile read
  recovers the deepest h_col every column has produced across the
  run. New parser keys `spall.h_spall_max` (default 10 mm safety
  cap) and `spall.h_spall_n_segments` (default 0 = fall back to
  `sp_n_segments`). `damage_law` else-branch is bit-for-bit
  unchanged.
- **Weibull first-fire is governed by order statistics, not the scale
  parameter.** Under `weibull.enabled = 1` with `N` surface cells, the
  cell that first satisfies `Sp ≥ 1` holds the realized maximum flaw
  length `a_max ≈ a0·(ln N)^(1/m)` (Galambos extreme-value asymptotic),
  NOT `a0` itself. Sp scales as `√a`, so first-fire onset ΔT is
  systematically *below* the deterministic Kant Eq. 15 prediction at
  `a = a0` (~7% below for `N = 576, m = 15`; ~80–110 °C below the Kant
  `[390, 560] °C` band that was derived from a single deterministic
  crack). This is why the Step 15c v1 sweep PASS gate was rejected —
  it scored an order-statistic quantity against a deterministic band.
  Resolved in §15c v2: score Weibull-on FEM first-fire against
  `Eq.15(a_max_realized)` computed from `max(flaw_a_mf)` on the t=0
  plotfile, with the deterministic Kant band kept as a heatmap
  colorbar reference for the Weibull-off verification run
  (`sp_kant_onset/input_verification`, 461.3 °C ∈ [390, 560],
  unchanged). **Always score Weibull-on first-fire against
  `Eq.15(a_max_realized)`, never `Eq.15(a0)` or the Kant band.**
- The Sp criterion stays sign-agnostic in `src/Numeric/SpCriterion.H`:
  `K_I` consumes whatever the `sigma_at_depth` callable returns; the
  `-sigma_xx` convention lives in the `UpdateSpAfterMechanics` caller, not
  the header. Do not move it into the header. Step 16's column h_spall
  scan likewise builds its `sigma_at_depth` callable in the integrator
  with the compression-positive convention; the header stays generic.
- **`h_spall_field_mf` and Rossi's measured observable are different
  physical quantities.** Step 16's diagnostic is a per-column running
  max of the deepest K_I-vs-K_Ic crossing across all steps — under a
  monotonic ramp this captures the LATE-ramp value, not the
  integrated event population. Rossi 2018 measures the spatial
  distribution of distinct cracks across all events during heating.
  Step 16's Rossi test reports the bands as informational only; the
  quantitative validation is Step 16b's scope using a new
  time-integrated 3-D `crack_event_count_mf` field. Edge-triggered
  counting (first crossing per cell only) and stress-relief
  option (ii) (freeze re-firing, no σ touch) are user-approved and
  documented in `rossi-validation-diagnostic-design.md` (repo root).
- **Step 16 geometry was forced off-packet by MLMG conditioning + free-
  lateral physics.** Working config for sp_weibull + spall.enabled +
  prescribed-T heating: `32^3 cubic cells in 50×50×50 mm + Kant-style
  lateral roller box + full-face heating` (heat patch radius >=
  face half-diagonal). Avoid: thin slabs (e.g. 30 mm deep at 1 mm
  cells — intermittent MLMG failures), non-power-of-2 meshes (can't
  coarsen multigrid past 1-2 levels — MLMG stalls at residual ~1e-5),
  free-lateral BCs with finite heated patch (central-column σ_xx
  drops to ~0.3× of 1-D ideal — Sp never fires within the ramp). The
  lateral-roller-box config is a "central column of an infinite
  heated zone" reading of Rossi, not a literal patch-plus-rim
  reproduction.
- **Sub-cell h_col events are scan-quantization-suppressed.** Step 16
  walks the depth scan in integer-cell strides. A column where K_I ≥
  K_Ic only at `s = 0` (h_col < dz) registers as h_col = 0 and is
  NOT counted as firing. At 1.5625 mm cells, this biases the GB/IG
  firing ratio to 100% GB (IG cells with `a0_ig = 4 µm` have Sp ≥ 1
  at the surface but fail the scan at `s = 1`). Step 16b added an
  8-iter sub-cell binary search bracketing `s_last_pass` and
  `s_last_pass + 1`; refines `h_col` to dz/256 (~6 µm). This partly
  lifts the quantization — Step 16b's GB fraction dropped from 100 %
  (Step 16) to 96 % (21 IG firings of 528). The remaining 96:4 ≈
  24:1 GB/IG ratio is **not** Rossi's 5:1 and **will not be fixed by
  Step 16c's depth-binning**; the path to 5:1 is higher peak σ_xx
  (longer/hotter ramp), a re-examined `a0_ig`, or eventually finer
  mesh (option (a)). Step 16c calls this out as a scoped-out
  limitation.
- **σ_xx(z) is local-per-z under 1-D confinement** (Kant-style
  lateral roller box + full-face heating, as used in Step 16's
  working geometry). `σ_xx = -E·α·(T(z) - T_ref)/(1-ν)` depends only
  on the local cell's temperature — zeroing σ at fired cells (design
  note §3 option (i)) does NOT change σ at deeper cells on this
  geometry, so option (i) is ineffective as an R1 fix here. Option
  (i) might help on patch+free-lateral geometries (stress
  redistribution through the cold rim) but that reopens Step 16's
  MLMG conditioning failures.
- **The `crack_event_count_mf` field is cell-centred and is
  insufficient for the Rossi depth profile at dz = 1.5625 mm.** Step
  16b's `cec` flips 0→1 on the first K_I ≥ K_Ic crossing per cell
  (edge-triggered, option (ii) freeze, no σ touch). But at this
  mesh, K_I crashes from above-K_Ic at z = 0 to below-K_Ic at z = dz
  in one cell, so only the top cell ever flags. The Rossi-comparable
  observable is Step 16c's per-step **h_col CSV log + post-process
  histogram** instead of a cell-centred event field. Sub-cell
  crack-tip depths from the bisection (`h_col ∈ (0, dz)`) are
  recorded in `h_spall_field_mf` as a running max per column, and
  Step 16c logs them per step for the histogram. Sub-cell σ_xx
  reconstruction (linear cell-to-cell vs T-derived vs node-based
  interp) is the key modelling knob Step 16c locks down explicitly.
- **Step 16c's option (L) σ_xx reconstruction has a sampling
  artefact at `z = dz/2 − a_f`.** Under linear cell-to-cell σ_xx,
  the top cell has no adjacent cell ABOVE the surface to
  interpolate from, so σ_xx is flat in `z ∈ [0, dz/2]`. K_I
  integration over `(z_crack, z_crack + a_f)` is then constant
  for `z_crack ∈ [0, dz/2 − a_f]` and starts dropping only past
  that. The bisection brackets `(0, dz)` and converges to the
  START of the K_I drop — h_col ≈ dz/2 − a_f for every firing
  column, regardless of a_f scale or firing time. At dz = 1.5625
  mm and a_f ∈ [4, 22] µm: h_col ∈ [759, 779] µm — verified
  arithmetically against 16c's observed [757, 787] µm range
  (matches to within bisection tolerance dz/256 ≈ 6 µm).
  Physically, option (L) discards the 20 % surface-T jump: at
  first-firing time t ≈ 80 s, `T_face = 693 K` (from
  prescribed-T `T_expr`) but `T_top_cell_centre ≈ 625 K`
  (thermal-diffusion gap across dz/2). **Step 16d picks option
  (a) — T-derived σ_xx with the prescribed-T face BC** — to
  address the root cause.
- **`rossi-validation-diagnostic-design.md` is the canonical
  reference for the Rossi validation chain.** §1-§5 cover the
  16b/16c design decisions (edge-triggered counting, stress-
  relief options i/ii/iii, sub-cell bisection, open questions).
  §6 is the Step 16 mechanics status. §7-§10 (added during 16c)
  cover the (L) artefact mechanism, Step 16d scope (option (a)
  only, with the K_I-not-K_Ic-driven argument against (c)), the
  mesh-resolution limit reasoning (two reconstruction points
  between face and top-cell centre means NO sub-cell strategy
  can manufacture features at <100 µm depths — only anisotropic
  z-refinement can), and anti-patterns (no sub-cell tweaks
  beyond (a), no K_Ic-only smoothing, no stop_time tuning past
  Rossi conditions). **§11 (added post-16d) marks the chain as
  CLOSED with R1 DEFERRED**: empirical convergence of 16b/16c/
  16d on three different reconstructions, all peaks outside
  [100, 200] µm; revisit conditions enumerated (anisotropic
  z-refinement, different reference experiment, or accepting
  the shape-without-band-location reading). The Rossi mechanics
  test stays as an ongoing P1-P4 regression; R1/R2/R3 reported
  but not project-gated. Step 17 (per-cell v_n + Zhang/Oglesby)
  and Step 18 (Kant confining-pressure) advance without the
  Rossi R1 gate.
- **Rossi R1 deferral checklist (do not re-litigate without
  reopening one of the §11 revisit paths)**: do NOT chase the
  Rossi peak on this mesh with further σ_xx reconstruction
  tweaks; do NOT smooth K_Ic-vs-z to chase the peak (K_I is
  σ_xx-driven, K_Ic continuity is second-order — verified by
  the 16c/16d empirical chain); do NOT extend stop_time past
  130 s; do NOT modify Weibull parameters (`a0_gb`, `a0_ig`, m)
  to shift R3 toward 5:1 — those are calibrated against Kant.
  The 24:1 R3 ratio reflects late-ramp σ_xx skew on this
  geometry, not a model bug.
- **Step 16d's σ_xx reconstruction is gated to `surface_patch.mode =
  prescribed_T`** because the T(z) lerp uses `surface_patch.T_f(time)`
  to pin T at z=0. Step 19 (active) forks this lambda: under
  prescribed_T it keeps the byte-identical T-derived form; under
  any other heat source (MMW beam, convective_flame, surface_patch
  disabled) it reads σ_xx from `stress_mf` via NodeToCellAverage —
  same pattern as the diagnostic Sp_field path in
  UpdateSpAfterMechanics. The mechanics solve must be active
  (`el.time_evolving = 1`) for the non-prescribed_T branch to work;
  Step 19 adds a fail-fast abort if it isn't. The prescribed_T
  branch byte-identity is the canonical witness for Step 19
  correctness (Step 18 `sp_kant_pressure_sweep` p=0 = 461.3 °C
  must stay unchanged).
- **Step 17 per-cell v_n + regime_field + local-Q infrastructure.**
  Under `spallation.model = sp_weibull && (spall.enabled ||
  vapor.enabled)`, the vapor-RoP path is a three-branch ladder:
  (1) `vapor_prescribed` (synthetic-test priority, unchanged),
  (2) `sp_weibull` local-Gaussian `Q(x,y,t) = (2P/(π·w0²))·
  exp(-2·r²/w0²)`, `rop_v = oneR·Q/(rho·L_tot)` (Step 17 NEW),
  (3) damage_law bulk-averaged `oneR·P/(rho·L_tot·A_beam)` (Step
  11 byte-identical). Column-winner precedence under sp_weibull:
  spall-priority (if Sp fires, regime=1; else if H≥H_vap,
  regime=2); damage_law's winner-on-h preserved. New
  `regime_field_mf` cell-centred plotfile diagnostic carries
  the canonical per-cell regime under sp_weibull; the scalar
  `regime_mf` stays for damage_law back-compat (same per-cell
  value, different field name). Beam class NOT refactored —
  Gaussian formula inlined at the vapor-RoP call site. No new
  parser keys. Local-Q at r=0 exactly equals the bulk-RoP value
  by construction (`A_beam = 0.5·π·w0²`); the local-Q signature
  is the radial decay (Pearson r = +0.988 verified in
  sp_v_n_regime). `hu_end_to_end` is the canonical damage_law
  simultaneous-spall+vap test — must stay byte-identical to
  Step 16d baseline (granite 40.0 s, sandstone 90.5 s) for any
  future change in this area.
- §15c shipped `DetachMode { PerFace, ConnectedCluster }` +
  `weibull.detach_mode` / `weibull.A_crit` parser, `Sp_cluster_id_mf`
  (top-z plane only, 0 = not firing or below-A_crit cluster, positive
  int = labelled), `UpdateSpClusters` 4-connected union-find with
  ReduceIntMax across MPI ranks, and per-step
  `<plot_file>_clusters.csv` (time, cluster_count, max_cluster_area,
  total_firing_area). `Sp_field_mf` is never modified by the labeller;
  cluster gating moves *which* cells detach but not *when* first-fire
  happens. Cross-level AMR connectivity is deferred (current `sp_*`
  runs are `max_level = 0`).
- The convective Robin BC on the **microstructure** heat path
  (`AdvanceMicrostructure`) and the depth-resolved **physical** stress
  path in `UpdateSpAfterMechanics` were both first exercised end-to-end
  by `sp_kant_onset/input_weibull` (single-phase Voronoi + convective
  Robin BC + `weibull.enabled`) and passed without code changes.
  `AdvanceMaterialImplicit` still does not consume `surface_patch` at
  all — implicit-path convective BC is out of scope.
- Step 14b patched `PrincipalStressRatio` (Step 10 production spall
  path) to return a `1e30` sentinel on mixed/tensile principal-stress
  states instead of 0 — physically, a free surface with tensile sigma_zz
  is more spall-prone, so the Hoek-Brown anisotropy gate should pass.
  Production Hu / Kant / Rossi validations no longer need the
  `spall.prescribed_ratio` workaround; that hook stays only for
  synthetic unit-style spall tests (`spall_event`,
  `regime_low_high_power`).
- Step 14 discovered a geometry mismatch in the production
  `damage_law` spall pipeline: `UpdateRemovalAfterCohesive` reads `D`
  at `k_top` only, but Hu's prescribed-T patch damages cells ~5 mm
  subsurface. Step 14b adds `spall.sample = top_cell | column_max |
  shallow_max` (default `top_cell` preserves Step 10/11/14 behavior;
  every prior input keeps its current behavior).
- The Hu validation chain (6d → 7b/7d → 8b → 11b → 14) is closed on
  the `damage_law` branch. Future Hu-style validations should use
  Step 14's `tests/MMWSpalling/hu_end_to_end/` as the canonical setup
  rather than spinning up new test directories.
- ARCHIVE_DONE.md and the Step 11b prose still mention
  `A_D = 0.000667` for sandstone, but the live
  `hu_spall_onset/input_sandstone2` and `hu_end_to_end/input_sandstone2`
  ship `0.005`. Treat the live test inputs as authoritative; the
  archive prose is historical.
- `src/BC/Operator/Elastic/ZloRoller321.H` is untracked in git but is
  load-bearing for `hu_spall_onset` (`el.bc.type = zlo_roller_321`)
  shipped by Step 11b. Treat as "intentionally untracked, awaiting a
  proper checkin"; do not delete from the working tree.
- The kept working-tree changes to `src/Solver/Nonlocal/Newton.H` and
  `src/Integrator/Base/Mechanics.H` (after Step 12's hunk-by-hunk
  cleanup) are: the `ZloRoller321` BC registration; the simplified
  `solver.solve(disp, rhs, model)` API with default-tol fallback in
  Newton; the post-solve `Util::RealFillBoundary(*disp_mf[lev],
  geom[lev])`. These are net improvements and stay. The reverted
  pieces were the `bx → domain` stencil swap, the `domain` arg to
  `GetBC()`, the rhs_mf defensive zeroing, and the AMR-investigation
  escape hatches in `Operator/`.
- Run test and validation simulations with 4 MPI ranks by default on this
  machine, e.g. `mpirun --oversubscribe --bind-to none -np 4 ...`, unless the
  test metadata or user request explicitly requires a different rank count.
- Step 6b fixed the old microstructure path that evolved `Temp` directly and
  skipped `InvertHtoT`; the microstructure path now advances conserved `H`.
- Step 6c split mineral topology from grain topology: `phase` is a mineral ID,
  `grain_id` is a grain label, `is_grain_boundary` is the future CZM/damage
  flag, and `is_gb` remains only a phase-boundary compatibility alias.
- `microstructure.grain_expr` is optional in expression mode; avoid
  comma-bearing parser expressions in inputs because ParmParse splits on
  commas.
- Old Step 4c is now Step 6d so Hu conduction validates the thermal driver
  before mechanics.
- The current Hu validation text has constants, BC trendlines, and qualitative
  field-comparison numbers, but not full digitized hole-1 time-series arrays.
  Do not invent Hu reference curves; add real digitized CSVs if curve checks
  are needed.
- Step 6d added `surface_patch.*` on the explicit global-material path;
  Step 15a extended the convective-flame branch to `AdvanceMicrostructure`
  as well. `AdvanceMaterialImplicit` still does not consume `surface_patch`
  at all (out of scope, ROADMAP Known Note).
- `surface_patch.T_expr` must be one ParmParse token with no whitespace, e.g.
  `638.22+3.74*t`.
- Step 15a wired `surface_patch.mode = prescribed_T | convective_flame`
  (default `prescribed_T`). Under `convective_flame` the Robin patch flux
  `h_fl*(T_flame-Tc)*inv_dz` is applied on the in-patch top row;
  `surface_patch.h_conv` is reused as `h_fl` and `surface_patch.T_flame` is
  a required new key. `surface_patch.T_ext` remains parsed-but-unused
  (reserved for a future non-patch convective exterior). The Step 8b
  top-node eigenstrain override is gated to `prescribed_T`; under
  `convective_flame` the top node reads solved `temp_mf`.
- Step 7b is a homogeneous Hu reproduction and now passes single-level. It ships
  with a bottom clamp instead of Hu's ideal bottom roller because the pure
  roller elastic solve lost convergence by `t = 6.2 s`.
- Do not tighten the Step 7b L3 stress ordering without revisiting the
  mechanical BC; the bottom-clamp case reverses Hu's L3 ordering. Keep stress
  sampling on node plotfiles for L3/L4.
- Step 7d is complete with shipped `20^3 + max_level=1` AMR Hu reruns. Keep the
  bottom-clamp and L3 reporting-only caveats unless the elastic
  BC/preconditioner changes.
- Exploratory `16^3 + max_level=2` Hu AMR runs showed static mechanics
  displacement reuse/prolongation can blow up after regrid; setting
  `el.zero_out_displacement = 1` stabilized Granite. Sandstone still needed
  `amr.abort_on_nan = 0` because dirty `Temp` ghost data triggered plotfile
  warnings, so do not promote that two-level setup until the ghost-fill warning
  is understood or fixed.
- For yt/numpy regressions on this machine, use
  `/Users/tzetze20/Desktop/code/.venv/bin/python`; `/usr/bin/python3` lacks
  `numpy`.

- **The continuous ablation front — design record (2026-09-21).** Source:
  `Claude_markdowns/2026-09-21b.md` (scope = recession law only; onset stays the
  K_I/Weibull test as a one-time per-column switch), reviewed by the user the
  same day. The central energy argument holds: with the flux BC computed from a
  fixed `T_abl` and the removed slab carrying its enthalpy, the steady rate is
  `q/(ρc(T_abl − T_∞))` with no mesh-dependent term, so the residual sensitivity
  really is lateral conduction only. Decisions taken, binding on D2j-a:
  1. **The slope-aware variant is DROPPED** (`spall.isotherm_normal` is not to
     be added). Three reasons: (i) the doc's gradient relation is inverted —
     `∂T/∂z = cos θ ∂T/∂n`, so `∂T/∂n = (∂T/∂z)/cos θ`, and as written it would
     over-remove by 1/cos²θ (2× at 45°, 4× at 60°), a silent slope-dependent
     error in exactly the 80° cone geometry under diagnosis; (ii) done correctly
     it is a **no-op** — moving each column's face down its own vertical line to
     `T = T_abl` puts the face set on the isosurface whatever the slope, so the
     cos θ factors cancel identically; (iii) D2b already measured the Step 20
     normal kernel **inert** under the pinned closure, and the needle artefact it
     was built for belongs to the per-column event increment this packet retires.
     **The slope correction that does matter is the one in the flux** (per
     horizontal vs true area) and it stays deferred to the cone packet.
  2. **Linear two-cell reconstruction is settled** (the quadratic is optional):
     its error scales with the *per-step* increment, not with dz. At Meier fluxes
     the overshoot rule gives ≈ 5 K/step ≈ 0.03 mm against a thermal length of
     3.2 mm, so the secant-vs-tangent bias costs ≈ 0.1 % in rate. A **dt**
     question, covered by the dt ladder.
  3. **The `δh ≤ dz/2` cap is a gate quantity, not a log line.** Pre-register
     max(δh/dz) and cap-hit counts per block: if the cap binds at 1 mm and not at
     2 mm, the run is void. The dangerous case is `T_0` just above `T_abl` with a
     near-flat vertical gradient, where `z_iso` divides by something small.
  4. **Registration must be conditional** (the `jet_P_face` pattern at
     `src/Integrator/MMWSpalling.H:600`). An unconditional
     `RegisterIntegratedVariable` or new plotfile field changes every
     `thermo.dat`/plotfile header even with the key off, and the 2639-file
     identity would fail on headers alone — discovered late, in the full sweep.
  5. **Abort on checkpoint restart under the front.** `abl_on`/`T_abl` are host
     state and are not saved (same class as `pin_T`), so a restart silently
     resets the surface to pre-onset. With 3 h+ runs, document-and-hope is not
     enough.
  6. **Diagnostics parity before the Meier pair.** `rim.py` and the d2e→d2i
     analysis chain read `removal_events.csv`, which the front would make
     unusable (a row per column per step). Ship the cumulative-recession field
     plus a converter, and show it reproduces D2i's rim numbers in `events`
     mode, or the "compare to D2i within 5 %" test is not like-for-like.
  7. **The pre-registered expectation straddles its own gate** (the doc expects a
     5–8 % residual against a 5 % limit). Do **not** widen the gate after the
     fact. Pre-register the *interpretation* instead: a failure at a gap
     consistent with D2j-b's measured edge-width sensitivity is diagnosed as
     lateral resolution, and the response is surface-refined AMR or the
     true-area flux — never a re-tune.
  8. **Resolution is marginal everywhere, not only at the edge.** κ/v ≈ 3.2 mm,
     so dz = 2 mm spans 1.6 cells of the ablation layer and 1 mm spans 3.2 —
     consistent with D2i's centre also drifting 9 mm by 550 s. A uniform 0.5 mm
     Meier run is a day or more, so **surface-refined AMR is a named decision
     point**, not a parenthesis, and is probably what makes this line finishable.
  9. **Cost note:** after onset the K_I scan can be skipped per column, so the
     front should be somewhat *cheaper* than the event model, not dearer.
  10. **Packet split** (the doc's single D2c-sized packet is too big):
      D2j-0 (the hot-disk test, no code) → D2j-a (vertical front + 1-D test +
      identity + regressions) → D2j-b (hot disk with the front) → D2j-c (the
      Meier pair).

- **THE STEP-FACE MECHANISM — the cause of the D2d–D2i mesh failures
  (D2j-0b, 2026-09-21; supersedes the D2d and D2h diagnoses).** Flux and losses
  act on a column's **top face only**; an exposed vertical step face is
  adiabatic. A firing top cell beside a taller, colder column therefore loses
  `k(T_top − T_wall)/dx` into that column's *interior at the step depth*
  (≈ 240 kW/m² at 2 mm, ≈ 480 at 1 mm, against q_pin ≈ 520). It stalls, goes
  cold, and becomes the next wall — an inward cascade. **Reproduced with no
  feet, no jet, no gas closure and without the `s ≤ 0` exclusion** (C0 cascades
  too, three times slower).
  - **Why the mesh dependence.** A **one-cell** step has ΔT ∝ dz, the two cancel,
    and the loss converges — which is why the static D2j-0 disk looked bounded.
    Once the step exceeds the thermal length (κ/v ≈ 3.2 mm) the neighbour sits at
    bulk temperature at any mesh, ΔT is pinned, and the flux density goes as
    1/dx. This is the continuum **corner singularity**: a pointwise firing
    criterion evaluated at a geometric corner is mesh-divergent by construction,
    and **refinement makes it worse**.
  - **Measured, not argued:** per-bin `v_1mm/v_2mm` is 0.996–1.007 wherever the
    local slope < 0.7 and falls to 0.95 … 0.53 in the six bands where the slope
    runs 0.7 → 8.2. The whole ≈ 10 % disk-mean offset lives where
    slope × dx ≳ 1.5 mm. **The half-cell removal increment is not a second
    mechanism** (the flat interior is where it would show, and S1 already had the
    flat 1-D case converging).
  - **Dead remedies, each on measurement, not argument:** wall heating above the
    plane (C1 ≡ C2 ≡ C3 band for band — the sink is 10–20 mm below the exposed
    surface); `spall.surface_normal = 1` (inert at 2 and 1 mm, the third
    independent measurement of its inertness); true-area flux × 1/cos θ on the
    top face (deposits the energy in the column's own top cell, never reaching
    the sink); a longer idle clock; more mesh pairs on the same configuration.
  - **The fix must act on the step**, i.e. heat the sink cell directly. Hence
    D2k. Design constraints found in review and binding on it: the side power
    must come **out of the enthalpy march's budget** (the D2a failure D2a2 exists
    to prevent, invisible in a fixed-h hot disk); a side-heated cell must obey
    `T_s = min(T_f, T_pin)` or an unremovable mid-wall cell becomes a heat
    reservoir that flatters the fix; the `s ≤ 0` exclusion must apply identically
    to side faces or the C1/C2/C3 result is silently invalidated; and `f` is
    pre-registered, never fitted (a staircase of one-cell steps has area
    `dx + Δz` against a true `√(dx² + Δz²)`, over-heating by ≤ 41 % at 45° and
    ≈ 16 % at the observed 80° cone).
  - **The arbiter already exists and is cheap:** `d2j0b_plane` C1 at 2 and 1 mm,
    3 min + 26 min. Any candidate remedy is tested there first. An acceptance of
    "the gap stands" is not enough — it must be **within 5 % in every block**, on
    the same estimator as the Meier gate, or a pass will not transfer.

- **D2k's outcome and the two live explanations (2026-09-22).** Side-face flux
  fixes the step-face cascade (arbiter −45 % → −3 %) but drives the Meier rate to
  4.07 m/h. **The decomposition is the thing to carry:** absorbed power 1639 →
  3204 W (2.0×) while the rate went 2.7×, because warming the wall also removes
  the lateral conduction sink (240–480 of ~520 kW/m²), so far more of the absorbed
  power converts to recession.
  - **The second factor is correct physics and must survive any closure.** Only
    the first — a stagnation-impingement `h` on a near-vertical face — is suspect.
  - **Explanation 2 is live and has not been tested:** over a *smooth* steep cone
    the wall jet is attached and the area really is ~6× the projection, so large
    absorbed power is not obviously wrong. What may be wrong is that `h(r, s)` is
    too peaked, so the model digs a pit with an 80° skirt instead of Meier's "half
    elongated spheroid" at Ø 85–93 mm. **Its discriminating signature is rate down
    *and* Ø up together**, which no cut to side flux can produce.
  - **The coupling that may kill explanation 1:** the wall stays warm because it
    is heated, so cascade suppression and the power excess ride the *same*
    parameter in opposite directions. If the f that restores 1.3–1.6 m/h is below
    the f at which the cascade returns, the wall closure cannot satisfy both gates
    and no closure should be written. D2l Gate 0 measures both bounds in ~1 h,
    key-only.
  - **`f` is an area correction and a diagnostic dial only.** It scales `q_side`
    linearly where an `h` scale would, which is why it maps the window — but no
    value found by scanning it may be kept as a setting. A closure value must be
    derived from a named heat-transfer source and pre-registered before the Meier
    run, then checked against Gate 0's window.
  - **Score three axes, never rate alone:** rate one-sided against 1.30 m/h, **Ø
    at matched feet height** against 85–93 mm, and **J/mm³ against Meier's 15.4**.
    Cutting side flux to recover ~1.4 m/h would give the right rate for the wrong
    reason — a narrower hole drilled more slowly, which is D2b's failure again.
  - **Hole shape is compared at matched feet height, never at matched time.**
    D2k's hole looks 30 mm wider at matched time purely from drilling 3× further.

- **Shape, not rate, is the honest score (D2l, 2026-09-22).** With the side-face
  term off — the D2i/D2d configuration — the model removes **1.9× too much rock**
  (4.58 vs 2.47 cm³/s) while the burner rate reads 1.49 m/h, in band. **The rate
  agreement was partly a coincidence of two errors.** Cut Ø is 78 mm against
  Meier's 85–93; the mouth flares to 145 mm. Freshly cut rock is ~78 mm at every
  setting — everything wider is removed *afterwards*, and scales with dwell.
  - **The excess volume and the too-cheap removal are one thing seen twice.** The
    excess 2.11 cm³/s costs ≈ 3.0 kW, ~46 % of absorbed power; remove it and the
    cost lands near 2.6 J/mm³, about what Meier's numbers imply at a comparable
    delivery. In his hole that heat soaked into the wall (the **thermally altered
    band** in the sawn section); in ours it removes rock.
  - **Scoring rules now standing:** score Ø(z) and volume, never rate alone;
    compare at matched feet height, never matched time; quote the volume (water
    filled) and the cut Ø, and keep the **mouth reported rather than scored** (it
    formed under the wellhead across two phases).
  - **Anchors named for the delivery-vs-cost question** (one equation, two
    unknowns, not separable by rate data): the altered-band width, and the
    mouth-flank behaviour with and without dilution.
  - **Reopened:** cuttings absorbing heat was withdrawn as degenerate with
    delivery efficiency *under a fixed-flux BC*; the enthalpy march made the
    budget finite, so that premise has expired.
  - **Largest unverified input:** the ~1900 K driving gas temperature. Meier's
    burner cooling water (160 L/h, thermocouples both ends, ΔT never reported) is
    186 W per kelvin — a 10–30 K rise is 1.9–5.6 kW, i.e. **100–300 K off the
    nozzle gas temperature**, and heat flow scales with (T_gas − 821 K). Cheapest
    question to ask Meier, largest leverage.

- **How the feet reach the heat distribution (code path, for D2m).**
  `FeetDescent` reads only live columns with r in [`foot_r_inner`,
  `foot_r_outer`] (28–40 mm), takes the 0.9 quantile per pad, and sets
  `z_target = z_foot + foot_standoff` with monotone descent.
  `BuildFlameColumns` then derives every column's `s = z_nozzle − z_face` and
  `h = h_expr(r, s)`, with `h = 0` where `s ≤ 0`. **The feet touch heating through
  exactly one scalar — the nozzle height.** Consequence: the 28–40 mm ring is the
  only rock held at the design stand-off for the whole run, so pad recession and
  ROP are locked together; the centre runs away (s_c 157–178 mm against a design
  50 mm) and rock outside 40 mm is left behind with a growing stand-off and
  falling heat. That predicts an erosion-rate transition at r ≈ `foot_r_outer`,
  which is what D2m tests.
