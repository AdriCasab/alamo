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

Active step:

- Step 10: Level-set surface advancement and spall detachment.
  See `ACTIVE_STEP.md`.

Next milestones after the active step:

- Step 11: Vaporisation removal and unified rate of penetration.
- Step 11b / 14 / 15: Hu onset/end-to-end and Kant/von Rohr threshold
  validations.
- Step 12 / 13 / 16: EOS, advanced AMR/adaptive timestep, and final
  depth/fracture-toughness parametric study.

## Validation Story

- Analytical/unit tests validate individual numerical components.
- Zhang/Oglesby validates the MMW source and thermal model.
- Hu validates the prescribed-temperature thermal-spallation chain:
  6d conduction -> 7b/7d thermoelastic stress -> 8b breakage indicator ->
  11b/14 onset/removal.
- Kant/von Rohr validates convective flame-threshold behavior.
- The final depth scan is exploratory; it has no direct experimental validation.

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
- Step 6d added `surface_patch.*` only on the explicit global-material path.
  `AdvanceMaterialImplicit` and `AdvanceMicrostructure` do not consume it yet.
- `surface_patch.T_expr` must be one ParmParse token with no whitespace, e.g.
  `638.22+3.74*t`.
- `surface_patch.T_ext` and `surface_patch.h_conv` are parsed but not used;
  non-patch exterior convection is deferred.
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
