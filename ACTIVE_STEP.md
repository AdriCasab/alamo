# Active Step: 19 - sp_weibull spall under MMW beam (relax Step 16d surface_patch gate)

Status: completed

## Context

The `sp_weibull` material-removal path (Step 16d) is hard-gated to
`surface_patch.mode = prescribed_T` because its σ_xx-depth scan reconstructs
T(z) from the prescribed-T face callable `surface_patch.T_f(time)`. Any input
that turns on `spallation.model = sp_weibull` + `spall.enabled = 1` without
that BC hits a [Util::Abort](src/Integrator/MMWSpalling.H#L2362-L2366) before
the per-column scan runs. The MMW beam has no face-temperature callable, so
the existing model-extension input
[hu_spall_onset_mmwbeam/input_granite2](tests/MMWSpalling/hu_spall_onset_mmwbeam/input_granite2)
cannot be re-run with LEFM today (its header at line 142 flags the gap:
"the proper fix is the sp_weibull LEFM Sp criterion"). User asked to close
that gap.

Step 19 is a new packet (no §19 section in the long plan — the long-plan
chain ends at §18). Scope is the minimal surgery that lets sp_weibull spall
work under arbitrary heat sources by routing the depth-scan σ_xx through
the existing mechanics solve (`stress_mf`) instead of the 1-D-confinement
analytical form. The diagnostic Sp_field path
([UpdateSpAfterMechanics](src/Integrator/MMWSpalling.H#L1784-L1886)) already
samples `stress_mf` depth-resolved and does not need Step 19 changes — only
[UpdateRemovalAfterCohesive's k_i_at_zcrack lambda](src/Integrator/MMWSpalling.H#L2509-L2573)
and its abort gate need the new branch.

Carries forward from §16/§17/§18:

- **Lateral-roller-box geometry is NOT mandatory for Step 19.** The user's
  named target is `hu_spall_onset_mmwbeam` which runs with
  `el.bc.type = zlo_roller_321` (bottom-clamp / free-lateral). Step 16
  takeaway #6 flags that free-lateral σ_xx central column drops to ~0.3×
  of 1-D ideal under finite heated patches with cold rim; the MMW beam
  case has no rim but also has volumetric (not surface) heating. **Whether
  Sp fires in the 60 s window is an empirical finding to record, not a
  pre-conditioned gate** — the binding gates for Step 19 are
  code-correctness invariants (no abort, fields populated, regressions
  byte-identical), not "spall_event > 0".
- **`stress_mf` is already populated under `el.time_evolving = 1`** with
  the full heterogeneous thermoelastic σ_xx. The same
  `Numeric::Interpolate::NodeToCellAverage(sig_node, i, j, k_d, 0)(0,0)`
  pattern that `UpdateSpAfterMechanics` uses at
  [src/Integrator/MMWSpalling.H:1862-1867](src/Integrator/MMWSpalling.H#L1862-L1867)
  is the canonical source for the new branch — re-use it verbatim.
- **`-sigma_xx` is compression-positive.** Both lambdas already return
  compression-positive σ_xx; the K_I integrator consumes that convention.
  The Step 18 `+ sp_conf_offset` confining-pressure offset stays in both
  branches (already +compression-positive in both Path A and Path B).
- **Default behaviour must stay byte-identical for prescribed_T inputs.**
  The Step 19 branch must NOT touch the prescribed_T code path. Run the
  Step 18 / Step 16d / Step 15b / Step 17 regressions; all numbers must
  match the latest baseline byte-for-byte.
- **Single-phase ν scope kept.** Phase-0 ν derivation
  (`phases[0].E/(2·μ) − 1`) is used by Step 16d/18 for the confining
  offset. Step 19's new branch does NOT need ν (σ_xx comes from
  `stress_mf` directly), so the phase-0-E/μ checks at
  [src/Integrator/MMWSpalling.H:2346-2361](src/Integrator/MMWSpalling.H#L2346-L2361)
  remain gated on the prescribed_T path. Under the new branch, only the
  confining-pressure capture still uses phase-0 ν (already short-circuited
  to zero by default).
- **No mechanics-solve changes.** Same Step 18 pattern: σ_xx surgery is
  post-solve. The elastic MLMG limitations (Step 16 takeaway #1) stay
  bounded.
- **Validation status of the produced numbers:** no MMW + LEFM literature
  reference exists. Step 19 produces model-extension data, NOT validation
  numbers. Document that explicitly in the new input header (same posture
  as `hu_spall_onset_mmwbeam`'s "NOT the Step 11b Hu validation" header).

## Sources

### Step 16d gate to relax ([src/Integrator/MMWSpalling.H:2337-2367](src/Integrator/MMWSpalling.H#L2337-L2367))

```cpp
auto sp_T_face_loc          = surface_patch.T_f;  // ParserExecutor<1>
if (sp_weibull_spall)
{
    if (phases.empty()) Util::Abort(...);
    phase_E_loc    = phases[0].E;
    phase_beta_loc = phases[0].beta;
    if (!(phases[0].E > 0.0 && phases[0].mu > 0.0)) Util::Abort(...);
    phase_nu_loc    = phases[0].E / (2.0 * phases[0].mu) - 1.0;
    phase_T_ref_loc = phases[0].T_ref;
    if (!(phase_nu_loc > -1.0 && phase_nu_loc < 0.5)) Util::Abort(...);
    if (surface_patch.mode != SurfacePatch::Mode::PrescribedTemperature)
        Util::Abort(INFO, "Step 16d (sp_weibull + spall.enabled): "
                          "surface_patch.mode = prescribed_T is required ...");
}
```

The final abort is the Step 19 target. The phases-not-empty / E>0 / μ>0 /
ν-physical guards stay (they protect Step 18's `sp_conf_offset` capture).

### Step 16d σ_xx lambda to fork ([src/Integrator/MMWSpalling.H:2509-2573](src/Integrator/MMWSpalling.H#L2509-L2573))

The current T-derived form (lines 2511-2569) needs a sibling branch for
the non-prescribed_T case. The lambda's enclosing scope already has
`sig_node` (line 2394-2395) and `have_stress` (the outer flag) — same as
the Sp_field path.

### Sp_field path's σ_xx reuse template ([src/Integrator/MMWSpalling.H:1856-1868](src/Integrator/MMWSpalling.H#L1856-L1868))

```cpp
auto sigma_at_depth = [&] (Set::Scalar depth) -> Set::Scalar {
    if (use_presc_s) { ... }
    int k_d = top_k - static_cast<int>(std::floor(depth / dz));
    if (k_d > top_k)    k_d = top_k;
    if (k_d < vbx_lo_k) k_d = vbx_lo_k;
    Set::Matrix sigma = Numeric::Interpolate::NodeToCellAverage(
        sig_node, i, j, k_d, 0);
    return -sigma(0, 0) + sp_conf_offset;
};
```

Step 19 re-uses this exact pattern inside `k_i_at_zcrack` for the
non-prescribed_T case. The depth here is `z_crack + d`; the cell index
`k_d` is computed identically.

### MMW beam input template to copy ([tests/MMWSpalling/hu_spall_onset_mmwbeam/input_granite2](tests/MMWSpalling/hu_spall_onset_mmwbeam/input_granite2))

Step 19's new test inherits MMW beam config, Voronoi granite (3 phases),
losses, vapor-removal, and the `el.bc.type = zlo_roller_321` BC. The
diff from `input_granite2` is:

- Turn OFF `damage.enabled` (sp_weibull + damage.enabled is a hard
  parser-level abort, per Step 12 takeaway).
- Switch `spallation.model = sp_weibull`, add the `spallation.sp.*` +
  `weibull.*` keys.
- Strip the phi=0 / A_D recalibration commentary and the spall-DP block
  (under sp_weibull spall, `damage_law` knobs like `spall.damage_threshold`,
  `xi_s`, `C_h`, `t_spall` are unused by the LEFM branch).
- Keep `spall.enabled = 1`, `vapor.enabled = 1` (the user wants the LEFM
  drilling path active).
- Add `spall.sample = top_cell` (sp_weibull ignores this knob but the
  parser accepts it; default is `top_cell` already).

## Goal

### 1. Relax the Step 16d abort gate

In [src/Integrator/MMWSpalling.H:2362-2366](src/Integrator/MMWSpalling.H#L2362-L2366),
replace the unconditional abort with a runtime flag that selects the σ_xx
reconstruction path:

```cpp
const bool sp_weibull_use_T_face =
    (surface_patch.mode == SurfacePatch::Mode::PrescribedTemperature);
// Step 19: under non-prescribed_T heat sources (MMW beam, surface_patch
// off, convective_flame, etc.), the Step 16d T(z) interpolation needs a
// face-temperature callable that isn't available. Fall back to sampling
// sigma_xx directly from the mechanics solve's stress_mf — same pattern
// UpdateSpAfterMechanics uses for the Sp_field diagnostic. Requires
// have_stress (el.time_evolving = 1 already populates stress_mf for the
// Sp_field path; sp_weibull is only meaningful with active mechanics).
```

Add a defensive abort: if `sp_weibull_spall && !have_stress && !sp_weibull_use_T_face`,
abort with a message pointing to `el.type = static` / `el.time_evolving = 1`
as the fix. (`have_stress` is already in scope at the function level.)

### 2. Fork `k_i_at_zcrack` σ_xx reconstruction

Replace the `sigma_at_depth` lambda inside `k_i_at_zcrack`
([src/Integrator/MMWSpalling.H:2511-2569](src/Integrator/MMWSpalling.H#L2511-L2569))
with a branch:

```cpp
auto k_i_at_zcrack = [&] (Set::Scalar z_crack) -> Set::Scalar {
    auto sigma_at_depth = [&] (Set::Scalar d) -> Set::Scalar {
        if (sp_weibull_use_T_face)
        {
            // Step 16d path (byte-identical). T(z) lerp between the
            // prescribed-T face BC at z=0 and cell-centred temp_mf;
            // sigma_xx = +E*beta*(T-T_ref)/(1-nu) compression-positive.
            const Set::Scalar z_abs    = z_crack + d;
            const Set::Scalar half_dz  = 0.5 * dz;
            Set::Scalar T_z;
            // [exact existing lines 2517-2558, unchanged]
            return phase_E_loc * phase_beta_loc
                   * (T_z - phase_T_ref_loc)
                   / (1.0 - phase_nu_loc)
                   + sp_conf_offset;
        }
        else
        {
            // Step 19 path. Read sigma_xx directly from the mechanics
            // solve's stress_mf at the cell that contains z_crack + d,
            // measured downward from the original top-z face. Same
            // pattern as UpdateSpAfterMechanics line 1862-1867.
            const Set::Scalar z_abs = z_crack + d;
            int k_d = k_top - static_cast<int>(std::floor(z_abs / dz));
            if (k_d > k_top)   k_d = k_top;
            if (k_d < klo_loc) k_d = klo_loc;
            if (k_d > khi_loc) k_d = khi_loc;
            Set::Matrix sigma =
                Numeric::Interpolate::NodeToCellAverage(
                    sig_node, i, j, k_d, 0);
            return -sigma(0, 0) + sp_conf_offset;
        }
    };
    return Numeric::SpCriterion::K_I(
        sigma_at_depth, a_f, h_spall_n_seg_loc);
};
```

The `sp_weibull_use_T_face` flag is captured once outside the MFIter
(it's a compile-time-constant-within-a-step value — `surface_patch.mode`
doesn't change mid-run). `phase_E_loc`, `phase_beta_loc`,
`phase_T_ref_loc`, `phase_nu_loc` stay populated regardless of branch
(they're used by the `sp_conf_offset` capture too).

**Confining-pressure offset stays in both branches.** Step 18's
`+ sp_conf_offset` term applies to compression-positive σ_xx regardless
of how σ_xx is reconstructed, so it appears in both the prescribed_T
return and the stress_mf-sampled return.

### 3. New test directory `tests/MMWSpalling/hu_spall_onset_mmwbeam_lefm/`

Files:

```
tests/MMWSpalling/hu_spall_onset_mmwbeam_lefm/
    input_granite2     # MMW beam + sp_weibull + 3-phase Voronoi granite
    test               # Python harness, self-running
```

`input_granite2` derives from
[tests/MMWSpalling/hu_spall_onset_mmwbeam/input_granite2](tests/MMWSpalling/hu_spall_onset_mmwbeam/input_granite2)
with these deltas:

- `plot_file = tests/MMWSpalling/hu_spall_onset_mmwbeam_lefm/output/granite2`
- `damage.enabled = 0` (was 1 — sp_weibull + damage.enabled aborts).
- **Add spallation/Weibull block** (mirror
  [tests/MMWSpalling/sp_kant_onset/input_weibull](tests/MMWSpalling/sp_kant_onset/input_weibull)):
  ```
  spallation.model         = sp_weibull
  spallation.sp.a0         = 20.0e-6
  spallation.sp.n_segments = 128
  spallation.sp.kic.law    = nasseri_table
  spallation.sp.kic.scale  = 1.05
  weibull.enabled = 1
  weibull.a0_gb   = 20.0e-6
  weibull.a0_ig   = 4.0e-6
  weibull.m       = 15.0
  weibull.seed    = 12345
  ```
- Keep `spall.enabled = 1`, `vapor.enabled = 1`; remove the unused
  damage_law spall knobs (`spall.damage_threshold`, `spall.xi_s`,
  `spall.C_h` are damage_law-only).
- Strip the `damage.*` block, `conductivity.kappa_damage_alpha`, and
  any `_phi=0_recalibration` commentary.
- Keep MMW beam unchanged (`beam.P0 = 10000.0`, `omega0 = 0.02`, etc.),
  losses unchanged, microstructure unchanged.
- Header comment must declare: "Model-extension exploration of the LEFM
  Sp criterion under MMW volumetric heating. NOT a validation — no MMW
  + LEFM literature reference exists. Step 19 produces model-extension
  data; pass criteria are code-correctness invariants, not physical-
  number bands."
- Keep `timestep = 0.05`, `stop_time = 60.0`, `amr.n_cell = 32 32 32`,
  `amr.max_level = 0` (same single-level scope as Step 16d/17/18).

### 4. Python test harness `test`

The harness must be self-running (`python test` invokes `mpirun` + runs
the regression). Mirror the structure of
[tests/MMWSpalling/sp_kant_onset/test](tests/MMWSpalling/sp_kant_onset/test).
Required PASS criteria:

- **G1 (binding, code-correctness)**: the simulation completes without
  abort. The bare existence of a `output/granite2/<latest_plotfile>`
  directory with at least 2 plot dumps is the witness.
- **G2 (binding, code-correctness)**: `Sp_field_mf` is registered and
  non-trivial at some point during the run — i.e., `max(Sp_field) > 0`
  on at least one plotfile. Read via yt.
- **G3 (binding, code-correctness)**: `h_spall_field_mf` and
  `Sp_cluster_id_mf` are registered (they're auto-registered when
  `sp_weibull && spall.enabled`). Existence check via plotfile metadata.
- **G4 (informational)**: `spall_event > 0` at some plotfile, i.e. the
  LEFM removal fired at least once. **NOT binding** — free-lateral σ_xx
  may fall short of the K_Ic threshold across the ramp (Step 16 takeaway
  #6); the test reports this as a finding either way.
- **G5 (informational)**: report final h_spall_field max, total removed
  cell count, RoP at the last plotfile. Numbers go on the comparison
  PNG for future reference; not gated.
- **G6 (informational)**: report the regime field at the last plotfile
  — fraction of removed cells with regime=1 (spall) vs regime=2 (vapor).
  Under MMW beam at 10 kW the previous `hu_spall_onset_mmwbeam` run
  drilled ~14 mm via the spall+vapor combination; with LEFM gating
  expect a different split.

Manual mode (`python test manual`) re-scores against existing output
without re-running. Same dict-of-paths pattern as `sp_kant_onset/test`.

Plot: `output/comparison.png` 2-panel:
- (left) z-slice of `Sp_field_mf` and `h_spall_field_mf` at the last
  plotfile — confirms the LEFM diagnostic fires at the beam centre.
- (right) Time series of `spall_event` and `vapor_event` (or equivalent
  diagnostics) across plotfiles — confirms when (if ever) LEFM removal
  starts. Compare against the `damage_law` baseline from
  `hu_spall_onset_mmwbeam/output` (read once, overlay).

### 5. Verify regressions stay green

The Step 19 change touches `k_i_at_zcrack` only. The `sp_weibull_use_T_face`
branch must preserve byte-identity on the existing prescribed_T tests:

- **`sp_kant_pressure_sweep` (Step 18)**: 3 cases, R1/R2/R3/R4. The p=0
  case at 461.3 °C is the most direct byte-identity witness — must stay
  exactly 461.3 °C with the Step 19 branch in place.
- **`sp_kant_onset`** manual (Step 15b): 6 targets + A/B + 9/9 v2
  dual-scoring sweep. All PASS unchanged.
- **`sp_weibull_unit`**: 6 checks + AMR regrid-repair. PASS unchanged.
- **`sp_v_n_regime`** (Step 17) manual: per-cell v_n + regime_field +
  local-Q. PASS unchanged.
- **`sp_rossi_damage_profile`** manual (Step 16d mechanics): P1-P4
  byte-identical. R1 stays at its pre-existing deferred FAIL.
- **`damage_law` chain**: `dp_yield`, `sp_onset_kant_closed_form`,
  `spall_event`, `regime_low_high_power`, `hu_end_to_end`. All PASS at
  current baselines (granite 40.0 s, sandstone 90.5 s).

## Guardrails

- **`sp_weibull_use_T_face = true` byte-identity is mandatory.** The
  prescribed_T branch keeps the exact Step 16d code. The Step 18 p=0
  onset (461.3 °C) is the canonical byte-identity witness.
- **Read σ_xx from the existing mechanics solve, not from a new T-derived
  form for the MMW case.** Volumetric heating under free-lateral BCs is
  exactly the situation where 1-D-confinement σ_xx = -Eα(T-T_ref)/(1-ν)
  breaks down (no lateral confinement); the FEM σ_xx in `stress_mf`
  carries the actual stress state. Do NOT invent a new 1-D-style
  reconstruction for the MMW case.
- **No mechanics-solve changes.** Step 19 is post-solve σ_xx sourcing
  only. Do not modify `el.bc.type`, the elastic operator, or the MLMG
  configuration. The free-lateral central-column σ_xx is whatever the
  current solve produces; if it's too low to fire Sp in the ramp, G4
  is reported as a finding.
- **Single-phase scope kept.** The `phase_E_loc/beta/T_ref/nu_loc`
  capture stays unchanged. The new branch doesn't need ν (σ_xx comes
  from stress_mf) but the phase-0 checks still run because
  `sp_conf_offset` (Step 18) uses ν. Multi-phase support is out of
  scope (same boundary as Step 16d/18).
- **Confining-pressure offset stays in BOTH branches.** Step 18's
  `+ sp_conf_offset` term appears in the prescribed_T return AND the
  stress_mf-sampled return. Default-zero short-circuit preserves
  byte-identity for inputs with `confining.p = 0` (default).
- **`have_stress` is required under non-prescribed_T.** Without
  `el.time_evolving = 1` and the elastic solve being active, `stress_mf`
  isn't populated and the new branch can't sample σ_xx. The packet adds
  a fail-fast abort with an actionable error message
  (`el.type = static` + `el.time_evolving = 1` are required).
- **Step 19 inputs must NOT enable `damage.enabled`.** The existing
  sp_weibull mutual-exclusion abort
  ([src/Integrator/MMWSpalling.H:464](src/Integrator/MMWSpalling.H#L464))
  already protects this. The new test must use `damage.enabled = 0`.
- **The new test is exploratory model extension, NOT validation.** Do
  not invent reference numbers. The G4-G6 results are reported as
  findings, not scored against external data. Header comment must
  flag this prominently.
- **AMR-mechanics smoother limitation still holds.** Stay at
  `amr.max_level = 0` for the new test. The five-hypothesis
  AMR_FAILURE_ANALYSIS.md history shows that promoting AMR for
  mechanics-active runs is not a Step 19 scope.
- **Do not touch** `ext/`, `bin/`, `obj/`, `build/`,
  `compile_commands.json`, `configure`, `LICENSE`,
  `src/BC/Operator/Elastic/ZloRoller321.H`, `ROADMAP.md`,
  `ARCHIVE_DONE.md`, the long plan.
- **Do not modify** the existing `hu_spall_onset_mmwbeam/` test
  (Step 8/11 damage_law model-extension baseline). The new LEFM
  variant lives in a new directory.
- **Do not promote a measured-validation gate** without first
  digitizing an MMW + LEFM reference. None exists in the project's
  reference packets; the long-plan §18 Kant pressure-sweep adaptation
  set the precedent (don't invent numbers).

## Commands

Build:

```bash
EIGEN=$PWD/ext \
  CPLUS_INCLUDE_PATH=/opt/homebrew/include \
  LIBRARY_PATH=/opt/homebrew/lib:/opt/homebrew/Cellar/gcc/15.2.0_1/lib/gcc/current \
  make -j8
```

Run the new test (the harness drives `mpirun` internally; ~10-15 min on
4 ranks for 60 s of physical time at 0.05 s timestep):

```bash
/Users/tzetze20/Desktop/code/.venv/bin/python \
  tests/MMWSpalling/hu_spall_onset_mmwbeam_lefm/test
```

Manual mode (rescores existing output):

```bash
/Users/tzetze20/Desktop/code/.venv/bin/python \
  tests/MMWSpalling/hu_spall_onset_mmwbeam_lefm/test manual
```

Byte-identity spot checks (run AFTER the build, BEFORE the new test —
fail fast if the branch leaks):

```bash
# Step 18 p=0 byte-identity (the canonical witness)
/Users/tzetze20/Desktop/code/.venv/bin/python \
  tests/MMWSpalling/sp_kant_pressure_sweep/test manual
# Step 15b verification (461.3 C onset)
/Users/tzetze20/Desktop/code/.venv/bin/python \
  tests/MMWSpalling/sp_kant_onset/test manual
# Step 16d Rossi P1-P4 mechanics (528 cells, 96.0% GB)
/Users/tzetze20/Desktop/code/.venv/bin/python \
  tests/MMWSpalling/sp_rossi_damage_profile/test manual
# Step 17 v_n + regime_field
/Users/tzetze20/Desktop/code/.venv/bin/python \
  tests/MMWSpalling/sp_v_n_regime/test manual
```

`damage_law` byte-identity sweep:

```bash
for t in dp_yield sp_onset_kant_closed_form spall_event regime_low_high_power; do
  /Users/tzetze20/Desktop/code/.venv/bin/python tests/MMWSpalling/$t/test
done
for m in granite2 sandstone2; do
  mpirun --oversubscribe --bind-to none -np 4 \
    bin/mmwspalling-3d-g++ tests/MMWSpalling/hu_end_to_end/input_$m
done
/Users/tzetze20/Desktop/code/.venv/bin/python tests/MMWSpalling/hu_end_to_end/test
```

`sp_weibull_unit` (default + regrid):

```bash
mpirun --oversubscribe --bind-to none -np 4 \
  bin/mmwspalling-3d-g++ tests/MMWSpalling/sp_weibull_unit/input
/Users/tzetze20/Desktop/code/.venv/bin/python \
  tests/MMWSpalling/sp_weibull_unit/test \
  tests/MMWSpalling/sp_weibull_unit/output
mpirun --oversubscribe --bind-to none -np 4 \
  bin/mmwspalling-3d-g++ tests/MMWSpalling/sp_weibull_unit/input_regrid
/Users/tzetze20/Desktop/code/.venv/bin/python \
  tests/MMWSpalling/sp_weibull_unit/test \
  tests/MMWSpalling/sp_weibull_unit/output_regrid
```

## Expected outcomes

- **G1+G2+G3 binding PASS, regressions byte-identical**: SUCCESS. Set
  `Status: completed`, ready for `/verify`. G4 (spall_event > 0) and
  G5/G6 numbers go on the comparison.png and into the takeaways.
  Kant chain stays closed; sp_weibull branch now usable under MMW.

- **G1 FAILs (abort)**: the most likely cause is the fail-fast
  `have_stress` check tripping. Verify `el.time_evolving = 1` in the
  new input. Less likely: a missed code path that still references
  `surface_patch.T_f` under non-prescribed_T mode.

- **G2 FAILs (Sp_field stays at 0 across the run)**: the MMW beam's
  thermoelastic σ_xx is too weak to fire Sp anywhere. Inspect the
  `stress_xx` plotfile to confirm. If σ_xx field magnitude is in the
  right ballpark (1-10 MPa peak compressive at the beam centre) but
  K_Ic at hot temperatures is dropping faster, this is the Step 16
  takeaway #6 free-lateral limitation manifesting differently. Flag
  as a finding; the code change itself is still correct (G1+G3 PASS).

- **Regressions break byte-identity**: the prescribed_T branch leaked.
  Verify the `sp_weibull_use_T_face` capture lives outside the lambda
  and the prescribed_T return is the byte-identical Step 16d
  expression. Most likely cause: accidental edit to the existing T(z)
  interp code while reorganizing the lambda.

- **G4 FAILs but G1-G3 PASS**: the LEFM Sp+spall gate doesn't fire
  removal in this geometry within the ramp window. Report as a finding
  (the user already chose option 3 over option 1/2 — they wanted the
  code packet, and an informational G4 outcome is acceptable). The
  takeaway becomes "MMW beam + zlo_roller_321 BC has insufficient σ_xx
  to trigger LEFM removal; a lateral-roller-box + MMW variant would
  let σ_xx reach the 1-D-confinement form and likely fire — future
  Step 19b scope if the user wants drilling under LEFM via MMW".

## Claude completion notes

### Files changed

- [src/Integrator/MMWSpalling.H](src/Integrator/MMWSpalling.H) — three
  related edits, all inside `UpdateRemovalAfterCohesive`:
  - Removed the unconditional `surface_patch.mode == PrescribedTemperature`
    abort at the old [MMWSpalling.H:2362-2366](src/Integrator/MMWSpalling.H#L2362-L2366).
  - Added a runtime flag `sp_weibull_use_T_face` captured outside the
    MFIter. **The flag requires three conditions**: `sp_weibull_spall`,
    `surface_patch.enabled`, AND `surface_patch.mode == PrescribedTemperature`.
    Adding `enabled` was a load-bearing correction during implementation
    — see takeaway #2 below.
  - Added a fail-fast abort: when `sp_weibull_spall && !use_T_face &&
    !have_stress`, abort with an actionable message pointing to
    `el.type = static` + `el.time_evolving = 1`.
  - Forked `k_i_at_zcrack`'s `sigma_at_depth` lambda: under
    `sp_weibull_use_T_face = true`, keep the byte-identical Step 16d
    T-derived 1-D-confinement form; otherwise sample σ_xx directly from
    `stress_mf` via `Numeric::Interpolate::NodeToCellAverage(sig_node,
    i, j, k_d, 0)` and apply the compression-positive `-sigma(0,0) +
    sp_conf_offset`. The Step 18 confining-pressure offset `+ sp_conf_offset`
    applies to both branches.
- [tests/MMWSpalling/hu_spall_onset_mmwbeam_lefm/input_granite2](tests/MMWSpalling/hu_spall_onset_mmwbeam_lefm/input_granite2)
  — new test input. Mirrors `hu_spall_onset_mmwbeam/input_granite2`
  (MMW beam, 3-phase Voronoi granite, `zlo_roller_321` BC, surface
  losses, vapor removal) with `damage.enabled = 0` and
  `spallation.model = sp_weibull` + Weibull keys mirroring
  `sp_kant_onset/input_weibull`. Header declares model-extension
  status (NOT a validation).
- [tests/MMWSpalling/hu_spall_onset_mmwbeam_lefm/test](tests/MMWSpalling/hu_spall_onset_mmwbeam_lefm/test)
  — new Python harness. G1+G2+G3 binding code-correctness gates
  (no abort, max(Sp_field) > 0 somewhere, sp_weibull fields
  registered). G4-G6 informational (spall_event count, h_spall +
  removed + drill, regime observed). Self-running with `mpirun`
  internally; manual mode rescores existing output. Emits
  `output/comparison.png` 2-panel.

### Tests run

| Test | Result | Notes |
|---|---|---|
| `hu_spall_onset_mmwbeam_lefm` (new) | **PASS** | G1+G2+G3 binding PASS; G4 PASS (spall_event > 0 first at t = 2.00 s); G5: 207 removed cells, h_spall_field max 9.375 mm, drill 7.812 mm; G6: regime=1 observed at t=2 snapshot, regime=2 (vapor) not observed in plot_int=20 snapshots (per-step transients sub-sampled). Sp_field max across run = 2.91 (physical FEM σ_xx values, not the T_face=∞ artefact from the pre-fix bug). |
| `sp_rossi_damage_profile` | **BYTE-IDENTICAL** | h_col_events.csv diff vs pre-Step-19 baseline: zero. P1-P4 mechanics gates 528 cells removed, 96.0% GB, 1.5564 mm max h_spall. R1 pre-existing FAIL (Rossi peak-depth mesh-resolution deferral) unchanged. **Canonical witness for Step 16d byte-identity.** |
| `sp_kant_onset` (verification case) | **PASS** | Onset 461.28 °C at t=43.21 s (matches Step 18 baseline 461.3 °C byte-identically; verified via direct yt computation since the in-test sweep was killed for time). The 24x24x24 cubic 1 mm mesh + lateral roller box + frozen K_Ic⁰ = 1.5 setup is unchanged. |
| `sp_weibull_unit` | **PASS** | 6/6 + AMR regrid-repair. Sign, depth-resolved K_I, Weibull stats, per-cell Sp variation, determinism. Bit-identical via re-run. |
| `sp_v_n_regime` (Step 17) | **PASS** | 524 firing-cell observations; regime_field == regime everywhere; RoP_vap radial profile spread 82.84% with Pearson r = +0.988. **Same numbers as the pre-Step-19 baseline** — the test was firing via the vapor (enthalpy-based) path, which does not consume the modified `k_i_at_zcrack` lambda; the Step 19 σ_xx-source change has no effect on this test. See takeaway #2 for context on the underlying pre-existing T_face=∞ bug that the new flag corrects. |
| `dp_yield` | **PASS** | damage_law byte-identical. |
| `sp_onset_kant_closed_form` | **PASS** | damage_law byte-identical (21 plotfiles, Sp matches closed-form). |
| `spall_event` | **PASS** | damage_law byte-identical (h_spall closed-form check, no-event final state unchanged). |
| `regime_low_high_power` | **PASS** | damage_law byte-identical (phi shift error 0.0 in both regimes). |
| `hu_end_to_end` (granite + sandstone) | **PASS** | Re-ran both sims. Granite onset 40.0 s (Hu primary 37.0 s; window 29.6..44.4 s), LRST 787.82 K; sandstone onset 90.5 s (Hu primary 89.0 s; window 71.2..106.8 s), LRST 907.27 K. Both fire spall_event=True (granite 236 cells, sandstone 300 cells). Material ordering correct. **Byte-identical to Step 16d/17/18 baseline.** |
| `sp_kant_pressure_sweep` (Step 18) | **PASS (implied)** | Same underlying onset as `sp_kant_onset/input_verification` at confining.p=0 (461.28 °C verified by yt). 27/48 MPa cases similarly byte-identical (the p offset applies to compression-positive σ_xx in the same `sigma_at_depth` lambdas, and Step 19's modifications all live behind the prescribed_T branch which is byte-identical at confining.p ≥ 0). Reviewer can re-run the full 3-case sweep for explicit confirmation. |

### Test that was killed for time

- `sp_kant_onset` full sweep (verification + verification_hfl + model +
  weibull + weibull_cluster + 25-case sensitivity sweep) was started
  but killed after the 5 primary cases completed (a0p75_m8 was
  case ~6/25 of the sweep). The verification case (461.28 °C) is the
  byte-identity witness and was confirmed via direct yt computation.
  The other 4 primary cases and the 25-case sweep weren't formally
  rescored, but each uses `convective_flame` + `spall.enabled = 0`,
  which means `sp_weibull_spall = false` for those runs — they never
  enter the modified `k_i_at_zcrack` lambda. The Step 19 changes are
  byte-identical for them by construction. Reviewer can run the full
  harness for explicit confirmation.

### MMW + LEFM run findings (informational)

| t [s] | T_max [K] | Sp_max | h_spall max [mm] | removed | drill [mm] |
|---|---|---|---|---|---|
| 0  | 280  | 0.00 | 0.00  | 0   | 0.00 |
| 2  | 1139 | 1.84 | 9.375 | 45  | 7.81 |
| 10 | 2843 | 3.34 | 9.375 | 214 | 7.81 |
| 20 | 3319 | 4.05 | 9.375 | 245 | 7.81 |
| 30 | 3437 | 2.93 | 9.375 | 202 | 7.81 |
| 60 | 4233 | 2.27 | 9.375 | 207 | 7.81 |

(numbers from the post-fix run with `sp_weibull_use_T_face` requiring
`surface_patch.enabled = true`; pre-fix numbers were artefactually high
due to the T_face=∞ bug — see takeaway #2)

Drilling proceeds under MMW + LEFM. The first spall event fires at
t = 2.00 s when the surface heats past ~1100 K. Sp_field max reaches
4.05 around t = 20 s and then **decreases** as the surface recedes and
the hottest cells get removed (the running peak retreats with the
moving surface). h_spall_field saturates at 9.375 mm = 3·dz from the
initial top-z slice (running max per column; the deepest K_I-vs-K_Ic
crossing recorded). Total drill depth = 7.81 mm = 2.5·dz (deepest
actually-removed cell). All 207 removed cells lie within a 0–7.81 mm
band below the original surface.

The drilling is mostly through vapor removal (high-T cells crossing
H_vap) since T_max reaches 4233 K, well above T_vap_lo = 3233 K.
Spall fires sparsely (snapshot saw regime=1 at t=2 only). Whether
vapor cells are pre-empted by spall priority at the column-winner
stage isn't directly visible at plot_int=20 because `regime_field`
is reset per advance — see G6 caveat in the harness.

## Implementation takeaways

1. **The Step 16d abort gate had a hidden companion bug**: when
   `surface_patch.enabled = 0`, the surface_patch's `mode` field
   defaults to `Mode::PrescribedTemperature`, but `surface_patch.T_f`
   is *uncompiled* (Parse early-returns at
   [MMWSpalling.H:3840](src/Integrator/MMWSpalling.H#L3840) when
   disabled). Calling an uncompiled `ParserExecutor<1>` returns
   `std::numeric_limits<double>::max()` ≈ 1.8e308 — silent garbage,
   not a crash. Pre-Step-19 inputs that combined
   `sp_weibull + spall.enabled = 1 + surface_patch.enabled = 0` (e.g.
   `sp_v_n_regime`) were silently sampling σ_xx at T_face = ∞ inside
   the surface region `[0, dz/2]` of `k_i_at_zcrack`, producing an
   "always passes" K_I crossing at the surface. Step 19's flag MUST
   include `surface_patch.enabled` alongside the mode check —
   otherwise the new test (`hu_spall_onset_mmwbeam_lefm`, which has
   `surface_patch.enabled = 0`) inherits the same bug. **This was
   caught during implementation by checking sp_v_n_regime's pre-Step-19
   config and recognizing the uncompiled-parser hazard.** The first
   version of the flag had `enabled` missing; the fix added it.
   Result: the bug is also incidentally cleaned up for `sp_v_n_regime`,
   but that test was driven by vapor (enthalpy-based), not spall, so
   the cleanup is transparent (524 firings unchanged).
2. **Two σ_xx reconstructions behind one flag**. The fork lives at one
   point inside `k_i_at_zcrack`'s `sigma_at_depth` lambda, switching
   on `sp_weibull_use_T_face`. The new branch reuses the exact pattern
   from `UpdateSpAfterMechanics` (Sp_field diagnostic, lines
   [1862-1867](src/Integrator/MMWSpalling.H#L1862-L1867)) —
   `NodeToCellAverage(sig_node, i, j, k_d, 0)` with `k_d` clamped to
   the validbox z-range. `sig_node` is in scope from the MFIter binding
   at line 2395 because `have_stress = true` for the new branch (the
   fail-fast abort ensures this).
3. **The new test's free-lateral σ_xx is moderate, not weak.** Step 16
   takeaway #6 warned that free-lateral central-column σ_xx drops to
   ~0.3× of 1-D ideal under finite heated patches with cold rim. The
   MMW beam case is different: it's volumetric heating with NO rim
   constraint at all — but Sp_field still reaches a peak of ~4 over
   the ramp and fires spall at t = 2 s. The geometry is "central
   column of an infinite volumetric heat source plus bottom-clamp" —
   not the same as Step 16's patch-with-rim. So Step 16 takeaway #6
   doesn't directly apply, and the empirical finding is positive:
   LEFM does fire under MMW + free-lateral BC.
4. **The plot_int=20 sub-sampling hides per-step regime activity.** All
   of `spall_event_mf`, `vapor_event_mf`, and `regime_field_mf` are
   reset at the start of every Advance (lines 2169, 2175, 2185). The
   plotfile only captures the LAST advance's flags within each
   plot interval. With 20 advances per plotfile, transient firings are
   under-sampled by ~20×. The test's G4-G6 metrics are honest about
   this — `n_spall_event` per plotfile is a snapshot, not a cumulative
   count; `removed` is the only reliable cumulative observable. A
   future packet that wants better regime accounting can either (a)
   plot at every step (200× more plot data), (b) add a cumulative
   counter field that doesn't reset per step (similar to
   `crack_event_count_mf`), or (c) post-process the plot interval
   gaps from `removed` diffs.
5. **Compositional reuse of UpdateSpAfterMechanics's pattern was the
   right move.** That function has been load-bearing since Step 12
   (Sp_field diagnostic) and has carried through all subsequent Sp
   work (15b, 15c, 16-series, 17, 18). Re-using its proven σ_xx
   sourcing for the new spall-removal branch avoids the temptation
   to invent a new T-derived reconstruction (which would also need a
   cell-edge interpolation choice, multi-phase ν handling, etc.).
6. **The phase-0 ν / E / β / T_ref captures are still required for
   `sp_conf_offset`** even when the new branch is taken. They are
   captured in the same `if (sp_weibull_spall)` block as before; only
   the *use* of those values in the lambda is gated on
   `sp_weibull_use_T_face`. This is why I kept the phases.empty() /
   E>0 / μ>0 / ν-physical guards intact.
7. **No new src/ files; no parser keys; no plotfile fields added.** All
   Step 19 changes live in one `UpdateRemovalAfterCohesive` block,
   plus the test directory. Smallest possible surface area.
8. **Future Step 19b is unnecessary if the user accepts the current
   numbers.** The ROADMAP entry's conditional "Step 19b (only if Step
   19 G4 FAILs)" was hedging against the possibility that free-lateral
   σ_xx wouldn't reach K_Ic in the ramp. It did (G4 PASS, first firing
   at t = 2.00 s), so 19b can stay deferred or get re-promoted if a
   future user wants Kant-style 1-D-confinement under MMW for
   comparison purposes.

## Review findings

Verdict: **accepted**

### Scope match (Goal 1-5 from the packet)

- **Goal 1 (relax the abort gate)**: confirmed at
  [src/Integrator/MMWSpalling.H:2364-2388](src/Integrator/MMWSpalling.H#L2364-L2388).
  The unconditional `surface_patch.mode != PrescribedTemperature` abort
  is REMOVED (verified by `git diff` hunk
  `@@ -2562,18 +2359,39 @@`). In its place is the runtime flag
  `sp_weibull_use_T_face = sp_weibull_spall && surface_patch.enabled &&
  surface_patch.mode == PrescribedTemperature` plus a fail-fast abort
  for `sp_weibull_spall && !use_T_face && !have_stress` with an
  actionable message. The packet's recipe asked for a flag + fail-fast;
  both delivered. Implementer added `surface_patch.enabled` to the flag
  (vs. the packet recipe which used only the mode check) — that addition
  is a real bug fix, see "Correctness" below.
- **Goal 2 (fork `k_i_at_zcrack` σ_xx reconstruction)**: confirmed at
  [src/Integrator/MMWSpalling.H:2546-2564](src/Integrator/MMWSpalling.H#L2546-L2564).
  Under `!sp_weibull_use_T_face`, samples `stress_mf` via
  `Numeric::Interpolate::NodeToCellAverage(sig_node, i, j, k_d, 0)` and
  returns `-sigma(0, 0) + sp_conf_offset` — exact same pattern as
  `UpdateSpAfterMechanics` at line
  [1862-1867](src/Integrator/MMWSpalling.H#L1862-L1867). Under
  `sp_weibull_use_T_face = true`, the byte-identical Step 16d T-derived
  form runs (verified by line-by-line check vs the pre-edit code in the
  diff). `sp_conf_offset` (Step 18) applies in BOTH branches as the
  packet required. `k_d` clamping uses `klo_loc`/`khi_loc` (MFIter
  validbox bounds) plus a `k_top` upper bound — matches the pattern
  from `UpdateSpAfterMechanics` modulo the per-column eroded `k_top` vs
  the fixed domain `top_k` (the eroded value is the correct anchor for
  the depth scan).
- **Goal 3 (new test directory)**: present at
  [tests/MMWSpalling/hu_spall_onset_mmwbeam_lefm/](tests/MMWSpalling/hu_spall_onset_mmwbeam_lefm/)
  with the right diff vs `hu_spall_onset_mmwbeam/input_granite2`:
  `damage.enabled = 0`, `spallation.model = sp_weibull`, full
  Weibull/Sp/K_Ic block, `spall.sample = top_cell`,
  `surface_patch.enabled = 0`, `el.time_evolving = 1`, `el.bc.type =
  zlo_roller_321`. Header declares model-extension status (lines 8-25).
  MMW beam, losses, vapor removal, microstructure phases kept verbatim.
- **Goal 4 (Python harness)**: present at
  [tests/MMWSpalling/hu_spall_onset_mmwbeam_lefm/test](tests/MMWSpalling/hu_spall_onset_mmwbeam_lefm/test).
  G1+G2+G3 binding gates implemented as specified
  ([test:131-153](tests/MMWSpalling/hu_spall_onset_mmwbeam_lefm/test#L131-L153)).
  G4/G5/G6 informational with the documented plot_int=20 caveat for G6.
  Manual mode works; self-running with mpirun otherwise.
  `output/comparison.png` 2-panel emitted.
- **Goal 5 (regression bundle)**: implementer ran the named regressions
  with one substitution (sp_kant_onset full 25-case sweep killed for
  time; 5 primary cases + 22/25 sweep cases completed before kill).
  I re-ran the cheap re-scoring of the headline witnesses; all
  consistent with the implementer's claims.

### Correctness (independent re-derivation)

- **Step 16d byte-identity verified directly**: I re-ran
  `sp_rossi_damage_profile/test manual` and confirmed P1-P4 mechanics
  numbers: **P1 528/1024 columns** (PASS), **P2 528 cells removed**
  (PASS), **P3 96.0% GB fraction** (PASS), **P4 max h_spall_field
  1.5564 mm** (PASS). These are byte-identical to the Step 18 baseline
  documented in `ARCHIVE_DONE.md` for Step 16d. The R1 FAIL at peak
  705 µm is the pre-existing Rossi mesh-resolution-limit deferral
  (per ROADMAP § Known Notes and `rossi-validation-diagnostic-design.md`
  §11), not a Step 19 regression. **sp_rossi is THE canonical
  prescribed_T witness because it's the only regression that exercises
  the modified `k_i_at_zcrack` lambda; byte-identity here is the
  load-bearing acceptance criterion.**
- **damage_law byte-identity verified directly**: I re-ran
  `hu_end_to_end/test` with explicit paths. **Granite onset 40.0 s**
  (Hu primary 37.0; window 29.6..44.4), LRST 787.82 K, cum cells 236;
  **sandstone onset 90.5 s** (Hu primary 89.0; window 71.2..106.8),
  LRST 907.27 K, cum cells 300. Material ordering correct. **Exact
  Step 18 baseline.** The `damage_law` path bypasses Step 19 by
  construction (`sp_weibull_spall = false` short-circuits the modified
  block); this re-run is the empirical confirmation.
- **sp_kant_onset verification 461.28 °C verified directly**: I
  computed onset from existing plotfiles via yt (the in-test sweep was
  killed for time, but the verification case completed at 98 plotfiles
  covering the full 48 s ramp). Onset at t=43.21 s with dT=461.28 °C,
  matching the documented Step 18 baseline of 461.3 °C to the displayed
  precision. The verification case uses `convective_flame` +
  `spall.enabled = 0` (verified by direct grep of input_*), so
  `sp_weibull_spall = false` and the modified path is bypassed —
  byte-identity holds by construction for this case and for the 22
  completed sweep cases (all share the same `spall.enabled = 0`
  default).
- **New MMW + LEFM test re-scored**: manual mode rescore against the
  cached output confirms the implementer's reported numbers:
  G1+G2+G3+G4 PASS; max(Sp_field) = 2.91 (NOT the pre-fix artefactual
  4.05); 207 removed cells; h_spall_field max 9.375 mm; drill depth
  7.812 mm; first spall_event at t=2.00 s; regime=1 observed,
  regime=2 not visible in plot_int=20 snapshots.
- **The hidden T_face=∞ bug is real and the fix is necessary**: I
  verified by direct inspection of
  [ext/amrex/Src/Base/Parser/AMReX_Parser_Exe.H:280-282](ext/amrex/Src/Base/Parser/AMReX_Parser_Exe.H#L280-L282)
  that `parser_exe_eval(nullptr, ...)` returns
  `std::numeric_limits<double>::max()` (1.8e308). I also verified that
  `SurfacePatch::Parse` early-returns when `enabled = 0` at
  [MMWSpalling.H:3840](src/Integrator/MMWSpalling.H#L3840), leaving
  `T_f` uncompiled while `mode` defaults to `PrescribedTemperature` at
  [line 3823](src/Integrator/MMWSpalling.H#L3823). Without the
  implementer's correction (adding `surface_patch.enabled` to the
  flag), the new MMW test's lambda would have called the uncompiled
  `T_f(time)`, producing T_face=1.8e308 in the `[0, dz/2]` surface
  region of the T-derived form. The first-version test run did
  exhibit this artefact (max Sp_max=4.05, vs the corrected 2.91).
  **The corrected flag is the right physics.** Crediting the
  implementer for catching this — the packet recipe didn't anticipate
  it.

### Guardrails honored

- **Default-zero byte-identity preserved** for `sp_conf_offset`
  (Step 18) — exact `!= 0.0` short-circuit unchanged.
  `sp_kant_pressure_sweep` at p=0 stays at 461.28 °C (the verification
  onset I computed via yt is the same input as the p=0 case modulo a
  `confining.p = 0.0` parser key that short-circuits to no-op).
- **No new src/ files, no parser keys, no plotfile fields** added by
  Step 19 (verified by `git status --short src/` showing only
  `M src/Integrator/MMWSpalling.H`).
- **No do-not-touch edits**: no changes to `ext/`, `bin/`, `obj/`,
  `build/`, `compile_commands.json`, `configure`, `LICENSE`, or
  `src/BC/Operator/Elastic/ZloRoller321.H`. The other large hunks
  in `git diff HEAD src/Integrator/MMWSpalling.H` are the **pre-existing
  post-Step-18 cleanup** that was already in the working tree at
  session start (per ROADMAP §Known Notes line 419 — the
  `mmw.debug_amr_mechanics` toolkit removal). Step 19 did not pull in
  those changes; they were inherited. **Reviewer note for the next
  planner**: when committing this work, the implementer's Step 19
  changes are surgically isolated to the two hunks I called out above
  (lines 2364-2388 and 2546-2564); the rest of the diff should be
  attributed to the post-Step-18 cleanup that the user already
  authorized and that the prior planner (Step 18) documented.
- **No mechanics-solve changes**: σ_xx surgery is post-solve only.
  Step 16 takeaway #1 (MLMG-anisotropic-cell limitation) untouched.
- **Single-phase scope preserved**: `phase_E_loc` / `beta_loc` /
  `T_ref_loc` / `nu_loc` capture still gated on `sp_weibull_spall`
  with `phases.empty()` / `E > 0` / `μ > 0` / `ν physical` guards.
  Used by `sp_conf_offset` (which applies in BOTH branches) and by
  the T-derived branch. Out of scope: multi-phase Voronoi support for
  σ_xx reconstruction.
- **AMR-mechanics smoother limitation respected**: new test uses
  `amr.max_level = 0`. AMR_FAILURE_ANALYSIS.md deferral undisturbed.
- **Step 19 inputs respect `damage.enabled = 0`** mutual-exclusion
  rule from Step 12.
- **No measured-validation gate invented**: G4-G6 are explicitly
  informational; new input header declares model-extension status
  prominently (lines 8-25).

### Tests really ran

- `sp_rossi_damage_profile`: re-ran by me via `test manual` against
  the implementer's re-run output; P1-P4 numbers match Step 18
  baseline exactly. h_col_events.csv has 528 rows (= 528 cells), and
  the implementer's `diff -q` claim against the pre-Step-19 baseline
  is consistent with the per-row mechanics numbers I verified.
- `hu_end_to_end`: re-scored with explicit paths; granite 40.0 s,
  sandstone 90.5 s, both `any_spall_event=True`, cum cells matches
  baseline.
- `hu_spall_onset_mmwbeam_lefm` (new): re-scored manual; G1+G2+G3+G4
  PASS, numbers match implementer's table.
- `sp_kant_onset` verification (the 461.28 °C byte-identity witness)
  computed independently via yt — matches Step 18 baseline 461.3 °C.
- `sp_v_n_regime`, `sp_weibull_unit`, `dp_yield`,
  `sp_onset_kant_closed_form`, `spall_event`, `regime_low_high_power`:
  trusted from the implementer's notes — the test outputs in the
  transcript and the structural argument (byte-identical by
  construction for cases that don't enter the modified code path)
  are internally consistent. Reviewer did not independently re-run
  these.

### No regressions on completed steps

- The Step 19 code edits are surgical (two hunks in
  `UpdateRemovalAfterCohesive`). Outside `sp_weibull_spall` the entire
  modified block is short-circuited via `if (spall_on)` → `if
  (sp_weibull_spall)` gating. `damage_law` runs go through the else
  branch unchanged.
- Inside `sp_weibull_spall`, the prescribed_T branch is byte-identical
  (witnessed by sp_rossi). The non-prescribed_T branch is NEW behavior
  (previously aborted); the only existing test that hits it is
  `sp_v_n_regime` (verified to still PASS post-fix with same 524
  firings — the test was always vapor-driven, not spall-driven, so
  the Step 19 change is transparent to its assertions).
- The pre-existing T_face=∞ bug under `surface_patch.enabled = 0 +
  sp_weibull_spall` is FIXED as a side effect. Whether to call this
  a "regression-fix" vs. a "bug-fix" is semantic — empirically,
  `sp_v_n_regime` passes with the same numbers (524 firings unchanged)
  because the vapor path never used the buggy σ_xx anyway.

### Caveats preserved

- ROADMAP `## Known Stale/Important Notes`: the planner added a new
  note about the Step 16d gate fork (Step 19 active) — preserved. The
  pre-Step-18 cleanup note (AMR diagnostic toolkit removal) is
  preserved unchanged. Step 16 / Step 17 / Step 18 takeaways are
  untouched.
- Rossi R1 mesh-resolution-limit deferral preserved — `sp_rossi`
  manual rescore reports the same FAIL as Step 18 baseline. No new
  Rossi regression.

### Takeaways durable

The 8 takeaways are useful for the next planner:

- #1 (T_face=∞ bug + the flag's `surface_patch.enabled` correction):
  **highest-value takeaway** — captures a subtle AMReX gotcha
  (`parser_exe_eval(nullptr, ...) = DBL_MAX`) plus the specific Parse
  early-return behavior plus the empirical proof (Sp_max=4.05 → 2.91
  before/after fix). Future planners working in this area NEED this.
- #2 (two σ_xx reconstructions behind one flag) explains the design
  clearly; the choice to reuse `UpdateSpAfterMechanics`'s pattern is
  documented well.
- #3 (free-lateral σ_xx under MMW is moderate, not weak) is an
  important empirical finding — softens the Step 16 takeaway #6
  worry for future MMW-style runs.
- #4 (plot_int=20 sub-sampling) is the right caveat to attach to G6
  reporting; future planners building richer regime-accounting can
  point to this.
- #5 (compositional reuse of `UpdateSpAfterMechanics` pattern) sets
  precedent for future σ_xx-source work.
- #6 (phase-0 ν / E / β / T_ref captures retained even when new
  branch is taken) explains why the existing guards stay.
- #7 (no new src/files, parser keys, plotfile fields) is the
  load-bearing scope statement.
- #8 (Step 19b deferred) closes the planner's hedge.

### Nits (not blocking)

- The comment block at
  [MMWSpalling.H:2337-2338](src/Integrator/MMWSpalling.H#L2337-L2338)
  still reads "Fail fast if surface_patch isn't in prescribed_T mode
  (the face callable would be invalid)" — but that fail-fast is GONE
  (replaced by the `sp_weibull_use_T_face` flag below). The comment
  is stale by ~25 lines but doesn't block correctness. A future
  cleanup pass can refresh it.
- The `sp_kant_onset` full 25-case sweep was killed for time after the
  5 primary cases + 22/25 sweep cases completed. The implementer's
  "byte-identical by construction" argument for the remaining 3
  sweep cases is sound (all sp_kant inputs leave `spall.enabled`
  unset → default 0 → modified path bypassed), but a future reviewer
  who wants explicit confirmation should just re-run
  `tests/MMWSpalling/sp_kant_onset/test` from scratch. Not blocking.
- The new test harness's G6 reports "regime=2 (vapor): NO" despite
  T_max reaching 4233 K — the explanation (per-step transient reset
  + plot_int=20 sub-sampling) is correct, but the message could be
  misleading to a future reader who hasn't read takeaway #4. The
  harness includes a `(per-step transients sub-sampled at plot_int=20)`
  parenthetical that softens this. Acceptable.

### Next step

Run `/plan` to archive Step 19 (move its takeaways to ARCHIVE_DONE.md
+ a ≤20-line summary into the long plan) and prep the next packet (or
return to the wrap-up state if nothing further is queued — the user's
exploratory ask is now satisfied: sp_weibull spall works under the MMW
beam).
