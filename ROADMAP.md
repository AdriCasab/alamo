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

Active step:

- **None** — project at wrap-up. The headline `sp_weibull`
  validation chain is closed and the headline `damage_law`
  chain (Hu) was closed at Step 14b. No new code packet is
  active; the next deliverable is a project-wide summary
  document for field engineers / geoscientists
  (`output/PROJECT_SUMMARY.html`). See the prior
  `ACTIVE_STEP.md` (Step 18) for the most recent completed
  packet's takeaways and review findings.

Previous active step (archived to ARCHIVE_DONE.md):

- Step 18 (Kant 2017 confining-pressure sweep validation,
  revised-approach Stage 5). One-parser-key src change:
  `confining.p` (default 0 — preserves all existing tests byte-
  identical) added to the σ_xx-reconstruction lambdas in both
  `UpdateSpAfterMechanics` (Sp_field path) and Step 16d's
  `k_i_at_zcrack` (depth-scan path) — adds `+p·ν/(1−ν)` to the
  compression-positive return so the K_I integrand sees the
  Kant Eq. 15 confining contribution. New test
  `tests/MMWSpalling/sp_kant_pressure_sweep/` with three input
  files (`p ∈ {0, 27, 48} MPa`) mirroring sp_kant_onset's lateral
  roller box + frozen K_Ic⁰ verification setup. PASS criteria:
  (R1, binding) at each `p`, FEM onset ΔT matches Kant Eq.15 with
  confining term within 8% (same tolerance as §15b target 3);
  (R2, binding) monotone decreasing across p; (R3, informational)
  at p=0, ΔT in Kant's 390–560 °C band; (R4, informational)
  closed-form slope `dΔT/dp = ν/(E·α) ≈ -0.93 K/MPa` matched
  within 20%. See `ACTIVE_STEP.md`.

Next milestones after the active step:

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
- Rossi 2018 validates Sp spall-depth profile (Step 16, `sp_weibull`).
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
  `MLNodeLaplacian`. Diagnostic helpers `DiagnoseAMRResidual` and
  `ApplyAMRF0Consistency` live in `src/Integrator/MMWSpalling.H` env-var-gated
  (`ALAMO_AMR_DIAG=1`, `ALAMO_LEV1_F0_FROM_COARSE`); default behavior unchanged.
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
