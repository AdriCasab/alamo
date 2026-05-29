# MMWSpalling test suite

Tests for the millimetre-wave (MMW) thermal-spallation integrator
(`src/Integrator/MMWSpalling.H`, binary `bin/mmwspalling-3d-g++`). See the
top-level `README.rst` for what the physics model itself does.

The suite is split into three tiers by **purpose**, so it's clear at a glance
what a directory is for:

```
tests/MMWSpalling/
  unit/          fast synthetic / analytic checks of one component each
  validation/    quantitative comparison against published experiments/theory
    zhang/         MMW source + thermal model        (Zhang & Oglesby)
    hu/            prescribed-T thermal-spallation    (Hu, damage_law branch)
    kant/          LEFM Sp onset + confining pressure (Kant 2017, sp_weibull branch)
    rossi/         spall depth distribution           (Rossi 2018, sp_weibull branch)
  extensions/    model-extension demos — NOT validations (no reference data)
```

If you are new here: start with the **unit** tier (each test isolates one
piece and runs in seconds), then read the **validation** tables below to see
how the assembled model is checked against the literature. The **extensions**
are exploratory runs, not pass/fail science.

---

## Running tests

Build first (from the repo root):

```bash
EIGEN=$PWD/ext \
  CPLUS_INCLUDE_PATH=/opt/homebrew/include \
  LIBRARY_PATH=/opt/homebrew/lib:/opt/homebrew/Cellar/gcc/15.2.0_1/lib/gcc/current \
  make -j8
```

The Python checkers use `yt`/`numpy`; on this machine use the venv interpreter:

```bash
VENV=/Users/tzetze20/Desktop/code/.venv/bin/python
```

There are two invocation styles (the test's header comment says which):

- **Self-running** (most `sp_*`, `spall_event`, `regime_low_high_power`,
  `sp_v_n_regime`, and the `hu_*` validations): the checker launches `mpirun`
  itself, then scores. Just run the checker:

  ```bash
  $VENV tests/MMWSpalling/unit/spall_event/test
  $VENV tests/MMWSpalling/validation/hu/hu_end_to_end/test
  ```

- **Two-step** (the older thermal/microstructure unit tests, e.g. `stefan`,
  `beam`, `equilibrium`, `voronoi`): run the sim, then score the output dir:

  ```bash
  mpirun --oversubscribe --bind-to none -np 4 \
    bin/mmwspalling-3d-g++ tests/MMWSpalling/unit/stefan/input
  $VENV tests/MMWSpalling/unit/stefan/test          # defaults to ./output
  ```

Use 4 MPI ranks by default. Each test self-locates the repo root and binary
from its own path, so the tier nesting does not need any extra configuration.

---

## Unit tier — component checks

Fast, synthetic or analytic. Each isolates one numerical/physical building
block. These are the first line of defence and should always pass.

| Test | What it checks |
|---|---|
| `skeleton` | Executable + integrator smoke test (builds, runs one step). |
| `stefan` | Enthalpy + phase fractions vs the analytic Stefan moving-front solution. |
| `beam` | Gaussian / Beer-Lambert MMW source energy deposition. |
| `equilibrium` | Radiation + convection surface-loss steady state. |
| `convective_patch` | Convective flame-jet Robin BC vs the Carslaw & Jaeger half-space solution. |
| `voronoi` | Voronoi mineral microstructure generation. |
| `grain_topology` | `grain_id` vs `is_grain_boundary` vs `is_phase_boundary` flags. |
| `heterogeneous_kappa` | Per-phase (heterogeneous) conductivity. |
| `heterogeneous_enthalpy` | Per-phase enthalpy / `H↔T` table consistency. |
| `thermal_stress` | Thermoelastic stress regression. |
| `dp_yield` | Drucker-Prager yield + irreversible damage evolution. |
| `gb_cohesive` | Bilinear grain-boundary cohesive-zone law + irreversibility. |
| `alpha_beta_transition` | Quartz α–β transformation strain. |
| `amr_microstructure_regrid` | AMR repair of discrete microstructure/derived fields on regrid. |
| `spall_event` | Deterministic spall detachment + surface advance (closed-form `h_spall`). |
| `regime_low_high_power` | Spall-vs-vaporisation regime selection + unified rate-of-penetration. |
| `sp_weibull_unit` | LEFM `K_I` integrator + Weibull per-cell flaw-distribution statistics. |
| `sp_onset_kant_closed_form` | Closed-form `Sp = K_I/K_Ic(T)` (prescribed stress/T). |
| `sp_v_n_regime` | Per-cell normal velocity + `regime_field` + local-Q (Gaussian) reconstruction. |

---

## Validation tier — literature comparison

Quantitative checks against published data/theory. **Status** reflects the
most recently recorded run; `*` = passes with documented caveats (see ROADMAP
`Known Stale/Important Notes`). Re-run to confirm before relying on a status.

### `validation/zhang/` — Zhang & Oglesby (MMW thermal)

| Test | Validates | Status |
|---|---|---|
| `zhang_oglesby` | MMW granite-heating: T-dependent `ρ,κ,c_p`, beam source, surface losses, `H↔T` inversion vs Zhang/Oglesby. | DEFERRED — not in the routine rotation (long-running; repeatedly skipped during fast iteration). Run to obtain a current verdict. |

### `validation/hu/` — Hu (prescribed-T thermal spallation, `damage_law` branch)

| Test | Validates | Status |
|---|---|---|
| `hu_conduction` | Prescribed-temperature surface-patch conduction (Granite 2 / Sandstone 2). | PASS |
| `hu_thermoelastic` | Granite 2 vs Sandstone 2 thermoelastic surface stress, single level. | PASS\* (bottom-clamp BC; L3 stress-ordering reporting-only) |
| `hu_thermoelastic_amr` | AMR variant of `hu_thermoelastic`. | PASS\* (same caveats) |
| `hu_breakage_index` | Breakage indicator `f_b = σ_v/σ_s` from thermoelastic stress. | PASS\* (surface targets pass; deep L6 `f_b` relaxed under bottom clamp) |
| `hu_spall_onset` | LRST/onset via DP-damage threshold proxy (homogeneous, removal off). | PASS\* (v1; sandstone damage-depth warning-only; AMR onset deferred) |
| `hu_end_to_end` | **Canonical** Hu chain end-to-end (conduction → thermoelastic → DP damage → spall removal). | PASS (granite onset 40.0 s, sandstone 90.5 s) |

### `validation/kant/` — Kant 2017 Central Aare granite (LEFM, `sp_weibull` branch)

| Test | Validates | Status |
|---|---|---|
| `sp_kant_onset` | Spall onset under a convective flame-jet BC at confining pressure p = 0; FEM vs Kant Eq. 15 + Weibull patchiness. | PASS (6/6 targets; verification onset 461.3 °C) |
| `sp_kant_pressure_sweep` | Confining-pressure dependence p ∈ {0, 27, 48} MPa. | PASS (R1+R2 binding; monotone, within 8% of Eq. 15) |

### `validation/rossi/` — Rossi 2018 (spall depth distribution, `sp_weibull` branch)

| Test | Validates | Status |
|---|---|---|
| `sp_rossi_damage_profile` | Spall crack-depth distribution vs Rossi's 100–200 µm peak. | **FAIL / DEFERRED** — mechanics gates P1–P4 PASS, but the depth-band (R1) peak lands at ~705 µm, not [100, 200] µm. This is a mesh-resolution limit, not a code bug; deferred per `rossi-validation-diagnostic-design.md` §11. The test stays as an ongoing P1–P4 mechanics regression; R1/R2/R3 are reported but not gated. |

---

## Extensions tier — model-extension demos (NOT validations)

No published reference exists for these configurations, so they are **sanity /
diagnostic** runs, not quantitative validations. They pass on "the run is
physically sane and damage evolved," and print drilling-progress diagnostics.
Do not treat a pass here as scientific validation.

| Test | Model configuration exercised | Pass criterion (sanity only) |
|---|---|---|
| `hu_spall_onset_voronoi` | Voronoi-microstructure granite spall-onset + drilling on the DP `damage_law` branch (heterogeneous, vs the homogeneous `hu_spall_onset` validation). | Plotfiles exist; `Temp`,`D` finite; `D ∈ [0,1]`; max `D` grows over the run. |
| `hu_spall_onset_mmwbeam` | Voronoi granite drilled by a 50 kW MMW beam (Beer-Lambert source, no prescribed-T patch), `damage_law`. | Same sanity gates; reports max T, deepest `D≥0.5`, removed cells, drill depth. |
| `hu_spall_onset_mmwbeam_lefm` | `sp_weibull` LEFM spall under an MMW beam (Step 19 — relaxes the prescribed-T gate so σ_xx comes from `stress_mf`). | Binding code-correctness gates: G1 no abort, G2 `max(Sp_field)>0`, G3 fields registered. G4 (`spall_event>0`) and drilling metrics are informational. |

---

## Where to start / what matters most

- **Quick smoke (seconds each):** `unit/stefan`, `unit/beam`, `unit/dp_yield`,
  `unit/gb_cohesive`, `unit/spall_event`, `unit/regime_low_high_power`,
  `unit/sp_weibull_unit`, `unit/sp_onset_kant_closed_form`.
- **Canonical "is the science right" set:** `validation/hu/hu_end_to_end`
  (damage_law), `validation/kant/sp_kant_pressure_sweep` (LEFM onset),
  `validation/rossi/sp_rossi_damage_profile` (Sp depth mechanics — P1–P4),
  `validation/zhang/zhang_oglesby` (MMW thermal).

Everything else is either a narrower slice of one of these or an exploratory
demo. The `damage_law` removal path is anchored by `hu_end_to_end`; the
`sp_weibull` removal path is anchored by `sp_rossi_damage_profile`.
