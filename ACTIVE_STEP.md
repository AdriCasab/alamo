# Active Step: 18 - Kant 2017 confining-pressure sweep validation

Status: completed

## Context

Step 18 is revised-approach Stage 5 and closes the Kant 2017 validation
chain on the `sp_weibull` branch (12c integrator unit → 15b `p = 0`
onset + Eq. 15 cross-check → 15c sensitivity sweep → **18 pressure
sweep**). The `p = 0` case is already validated at 3.5% rel err vs
Kant Eq. 15 (§15b target 3); Step 18 extends to `p ∈ {0, 27, 48} MPa`,
adding the confining-pressure contribution to σ_xx via the closed-form
Kant Eq. 15 term `−p·ν/(1−ν)`.

The implementation is minimal: one new parser key `confining.p`
(default 0, byte-identical to all existing tests) and a single
additive offset in the σ_xx-reconstruction lambdas. The Kant Eq. 15
form
```
Θ_onset = (1−ν)/(E·α) · [ K_Ic/(2·1.763·√(a/π)) − p·ν/(1−ν) ]
```
is derived in `kant-2017-central-aare-validation-reference.md` line 96.
Under 1-D-confined geometry (Kant-style lateral roller box), the
confining pressure contributes `+p·ν/(1−ν)` to the compression-
positive σ_xx that the K_I integrand sees — same offset added to the
lambda return.

Carries forward from §15b/c, §16-series, §17:
- **Lateral roller box** (`el.bc.type = constant` with all 26 regions
  configured as in `sp_kant_onset/input_verification`) is the canonical
  1-D-confinement geometry; this is what makes σ_xx ≈ -E·α·dT/(1-ν)
  hold per-cell. Reuse it.
- **Frozen K_Ic⁰ = 1.5 MPa·m^0.5** matches Kant's Table-2 prediction
  band (390–560 °C at `p = 0`). Step 18 scores against Kant's
  predicted-band intersection at each pressure.
- **§15b target 3** (FEM-vs-Eq.15 within 8%) sets the tolerance bar;
  §15b achieved 3.5% at `p = 0`. Step 18 extends this to non-zero
  pressures.
- **No measured pressure-sweep data** in the project's validation
  references (`kant-2017-central-aare-validation-reference.md` lists
  only the `p = 0` measurement window 553–694 °C). The packet's
  binding gate is the closed-form Eq. 15 cross-check at each
  pressure, not a measured-slope match. R4 (slope match) is
  informational since Kant's measured-pressure-sweep digitization
  is not in this packet.
- **No spall.enabled / vapor.enabled** in Step 18 inputs. The test is
  pure-onset, mirroring `sp_kant_onset/input_verification`. The
  Step 17 per-cell v_n and Step 16 h_spall mechanics are not
  exercised — they should stay byte-identical via the disabled gates.
- **`damage_law` byte-identical** — Step 18 changes only modify the
  σ_xx lambdas, which are in the `sp_weibull` path. The `damage_law`
  branch doesn't use those lambdas. Verified by spot-checks.
- **Existing tests must stay PASS** with `confining.p = 0` (the
  default). This is the most important guardrail — the offset must be
  ZERO when `confining.p` is unset.

## Sources

### Long plan §18 ([in-main-tex-you-will-quizzical-treasure.md:638-652](in-main-tex-you-will-quizzical-treasure.md#L638-L652))

> Apply lithostatic confining pressure as a stress-field BC (uniform
> pre-stress on side faces from `confining.p` in the input), not
> through a `Ψ_0` rescaling.
>
> **18a. Verification.** `tests/MMWSpalling/sp_kant_pressure_sweep/`:
> - Central Aare granite under flame conditions matching [kant-2017].
> - Pressure sweep `p ∈ {0, 27, 48} MPa`.
> - Pass: predicted onset ΔT at p = 0 within Kant's measured
>   553–694 °C bracket using `K_Ic(T)` from Nasseri 2007.
> - Pass: monotone-decreasing predicted onset ΔT across the three
>   pressures; relative slope matches Kant's measured trend within
>   ±30%.

**Packet adaptation**: the long plan's "matches Kant's measured trend
within ±30 %" gate requires Kant's measured `p ∈ {27, 48} MPa` data,
which isn't in the project's validation references. Replace with the
Eq. 15 closed-form slope (which IS in the reference at line 96). The
measured-slope gate stays as a future packet's scope if/when Kant's
data is digitized. Also use Kant's predicted band (390–560 °C at
p = 0, Table 2) instead of the measured anchor (553–694 °C) since
the predicted band is what §15b validated against.

### Kant Eq. 15 with confining pressure ([kant-2017-central-aare-validation-reference.md:92-103](kant-2017-central-aare-validation-reference.md#L92-L103))

```
Θ_onset = (1−ν)/(E·α) · [ K_Ic/(2·1.763·√(a/π)) − p·ν/(1−ν) ]
```

`∫₀¹ f_K(ξ) dξ = 1.763` (Kant Eq. 15 prefactor; the `1/(2·√(a/π))`
form in the file is the same up to the `√(π/(4a))` factoring). At
`p = 0`, reduces to `(1−ν)/(E·α) · K_Ic/(2·1.763·√(a/π))` —
verified by §15b target 3 at 3.5 % rel err.

Central Aare granite + frozen K_Ic⁰ values (= sp_kant_onset/
input_verification):

| param | value |
|---|---|
| E | 35.0e9 Pa |
| ν | 0.26 |
| α (linear) | 8.0e-6 /K |
| K_Ic⁰ | 1.5e6 Pa·m^0.5 |
| a₀ | 20.0e-6 m |

Pre-computed Eq. 15 predictions (Python-side reference):

| p (MPa) | Θ_onset (Eq.15) | dΘ/dp = -ν/(E·α) |
|---|---|---|
| 0  | 445.6 °C | -0.9286 K/MPa |
| 27 | 420.5 °C | (same slope, linear in p) |
| 48 | 401.0 °C | |

Slope: `dΘ/dp = -ν/(E·α) = -0.26/(35e9·8e-6) = -9.286e-7 K/Pa =
-0.9286 K/MPa`. The Eq. 15 ΔT decreases LINEARLY with p (the
confining term is linear in p; the K_Ic-related bracket is constant
under frozen K_Ic⁰).

### Existing infrastructure to reuse (do NOT rewrite)

- **`sp_kant_onset/input_verification`** at
  [tests/MMWSpalling/sp_kant_onset/input_verification](tests/MMWSpalling/sp_kant_onset/input_verification)
  is the canonical p = 0 input. Frozen K_Ic⁰, lateral roller box,
  24³ cubic 1 mm mesh, `convective_flame` full-face heating, full
  mechanics, `spall.enabled = 0` (onset only). The Step 18 inputs
  should be **identical copies + `confining.p = ...`**.
- **`UpdateSpAfterMechanics`** at
  [src/Integrator/MMWSpalling.H:1927](src/Integrator/MMWSpalling.H#L1927)
  has the `sigma_at_depth` lambda inside the K_I integrator call
  (line ~1988-2000). The compression-positive return is
  `return -sigma(0, 0);` — add the confining offset here.
- **Step 16d's `k_i_at_zcrack`** lambda at
  [src/Integrator/MMWSpalling.H:2680-2742](src/Integrator/MMWSpalling.H#L2680-L2742)
  has its own σ_xx reconstruction (T-derived, option (a)). Add the
  confining offset to the compression-positive return there too.
  This is for the depth-scan path under `spall.enabled = 1` (Step 18
  doesn't enable spall, but the change keeps both paths consistent
  so a future packet combining confining + sp_weibull material
  removal doesn't have a mismatch).
- **`phases[0].mu`** is the shear modulus; `nu = E/(2·mu) − 1` is the
  derivation Step 16d already uses (`phase_nu_loc`). The same value
  is used in the confining offset `p·ν/(1−ν)`.

## Goal

### 1. New parser key `confining.p`

Add a new parser key with **default 0.0** (Pa). The default-zero
behaviour is the critical byte-identity guarantee for all existing
tests. Plumb a `Set::Scalar confining_p` member into the
`MMWSpalling` class near the spallation-related members; parse in
`ParseSpallationModel` or a new `ParseConfining` static (small enough
that a new ParseConfining isn't needed — fold into the existing
ParseSpall* / ParseSpallationModel block alongside the other Sp
parameters):

```cpp
// In MMWSpalling class (near sp_a0, sp_n_segments etc.):
Set::Scalar confining_p = 0.0;  // Pa

// In ParseSpallationModel or appropriate Parse static:
pp.query_default("confining.p", value.confining_p, 0.0);
if (!std::isfinite(value.confining_p))
    Util::Abort(INFO, "confining.p must be finite, got ", value.confining_p);
```

Negative `confining.p` is physically tensile pre-stress (unusual but
not strictly invalid); allow it but document in the comment that
positive = compressive confining pressure (Kant convention).

### 2. σ_xx-lambda offset in both reconstruction paths

**Path A**: `UpdateSpAfterMechanics` sigma_at_depth lambda at
[src/Integrator/MMWSpalling.H:1988-2000](src/Integrator/MMWSpalling.H#L1988-L2000):

```cpp
auto sigma_at_depth = [&] (Set::Scalar depth) -> Set::Scalar {
    if (use_presc_s) {
        Set::Scalar z_phys = z_top_face - depth;
        return -presc_s_f(xc, yc, z_phys, time) + sp_conf_offset;
    }
    int k_d = top_k - static_cast<int>(std::floor(depth / dz));
    if (k_d > top_k)    k_d = top_k;
    if (k_d < vbx_lo_k) k_d = vbx_lo_k;
    Set::Matrix sigma = Numeric::Interpolate::NodeToCellAverage(
        sig_node, i, j, k_d, 0);
    return -sigma(0, 0) + sp_conf_offset;  // Step 18: confining contribution
};
```

Where `sp_conf_offset` is captured outside the MFIter loop (similar
to the other constants like `K_Ic0_`, `n_seg`, etc.):

```cpp
// Step 18: confining-pressure contribution to compression-positive
// sigma_xx. Per Kant Eq. 15, p enters as +p·nu/(1-nu) in the
// compressive sigma_xx field under 1-D-confined geometry.
// Phase-0 properties (single-phase Voronoi for Kant; multi-phase is
// out of scope same as Step 16d).
Set::Scalar nu_local = 0.0;
if (!phases.empty() && phases[0].E > 0.0 && phases[0].mu > 0.0)
    nu_local = phases[0].E / (2.0 * phases[0].mu) - 1.0;
const Set::Scalar sp_conf_offset =
    (confining_p != 0.0 && nu_local != 1.0)
        ? confining_p * nu_local / (1.0 - nu_local)
        : 0.0;
```

The phase-0-derived ν matches Step 16d's `phase_nu_loc`. The
`confining_p != 0.0` short-circuit ensures byte-identity when
`confining.p` isn't set (default 0).

**Path B**: Step 16d's `k_i_at_zcrack` lambda at
[src/Integrator/MMWSpalling.H:2680-2742](src/Integrator/MMWSpalling.H#L2680-L2742),
the `return` at line 2739-2741. Same offset:

```cpp
return phase_E_loc * phase_beta_loc * (T_z - phase_T_ref_loc)
       / (1.0 - phase_nu_loc) + sp_conf_offset;
```

Capture `sp_conf_offset` in the same block as `phase_E_loc` etc.
(near line 2492-2519 in MMWSpalling.H).

**Critical**: with `confining.p = 0` (default), both lambdas return
exactly their pre-Step-18 values. Verify byte-identity via the
regression bundle.

### 3. New test directory `tests/MMWSpalling/sp_kant_pressure_sweep/`

Three input files, identical except for `confining.p`:

```
input_p0   : confining.p = 0.0
input_p27  : confining.p = 27.0e6
input_p48  : confining.p = 48.0e6
```

Each input is a near-copy of
`tests/MMWSpalling/sp_kant_onset/input_verification`:
- `microstructure` single-phase Central Aare granite (homogeneous —
  Voronoi grain pattern but `nphases = 1`, mechanically uniform).
- Frozen `kic.law = linear`, `K_Ic0 = 1.5e6`, `T_ref = 293.15`,
  `T_melt = 1e9` (effectively frozen).
- Lateral roller box (`el.bc.type = constant`, all 26 regions
  specified — copy verbatim from input_verification).
- `surface_patch.mode = convective_flame`, `T_flame = 1673.15`,
  `h_conv = 150` (h_fl), full-face heating.
- `spallation.model = sp_weibull` + `spallation.sp.a0 = 20.0e-6`,
  `n_segments = 128`.
- `weibull.enabled = 0` (no Weibull distribution — deterministic
  a₀ as in input_verification, matches Kant's single-flaw model).
- **`spall.enabled = 0`** (onset-only validation, no material
  removal).
- **`vapor.enabled = 0`** (no vaporization).
- New: `confining.p = <value>` parser key.

Plot cadence ~5 s, stop_time ~35-43 s (mirror input_verification's
runtime; the onset time shifts slightly with p but stays within the
same ramp window).

### 4. Python harness `test`

Mirror `tests/MMWSpalling/sp_kant_onset/test`'s structure: a `CASES`
dict driving the three sims (or one if-manual fallback), reading
Sp_field plotfiles to find onset (first `t` where any top-z cell has
`Sp_field >= 1`), and computing the FEM onset ΔT.

PASS criteria:

- **R1 (binding)**: For each `p ∈ {0, 27, 48}` MPa, compute Kant Eq.
  15 analytically:
  ```python
  theta_eq15 = (1 - nu) / (E * alpha) * (
      K_Ic / (2 * I_FK * math.sqrt(a / math.pi))
      - p * nu / (1 - nu)
  )
  ```
  (E=35e9, ν=0.26, α=8e-6, K_Ic=1.5e6, a=20e-6, I_FK=1.763.) Check
  `|ΔT_FEM(p) - ΔT_Eq15(p)| / ΔT_Eq15(p) <= 0.08` for all three p.
  Same tolerance as §15b target 3.
- **R2 (binding)**: `ΔT_FEM(p=0) > ΔT_FEM(p=27) > ΔT_FEM(p=48)`
  strictly. Monotone decreasing.
- **R3 (informational)**: at p = 0, ΔT in Kant's Table-2 predicted
  band 390–560 °C. (§15b verification baseline.)
- **R4 (informational)**: FEM slope `dΔT_FEM/dp` matches Eq. 15
  closed-form slope `-ν/(E·α) = -0.9286 K/MPa` within 20 %.
  Computed from linear regression on the 3 (p, ΔT_FEM) points.
- **Anchor (no gate)**: Kant's measured spalling temperature window
  553–694 °C at p = 0 is shown on the plot but NOT scored against
  (per the validation reference §6 line 116-119).

Manual mode (`python test <anyarg>`) re-scores against existing
output without re-running. Mirror `sp_kant_onset/test`'s
`CASES`-dict pattern.

Plot: `output/comparison.png` 2-panel:
- (left) ΔT_FEM vs p (3 points) overlaid with Eq. 15 closed-form
  line + Kant predicted band [390, 560] °C shown as horizontal
  shading + Kant measured anchor [553, 694] °C as dashed band.
- (right) Per-case Sp_field firing pattern at the top-z slice at
  the onset plotfile (3 small panels), confirming the firing
  mechanism is consistent across p.

### 5. Verify regressions stay green

- **`damage_law` byte-identical** (the `confining.p` parser key
  defaults to 0; the offset is short-circuited when
  `confining.p == 0.0`; damage_law doesn't use the σ_xx lambdas
  anyway):
  - `dp_yield`, `sp_onset_kant_closed_form`, `spall_event`,
    `regime_low_high_power` — all PASS at machine precision.
  - `hu_end_to_end` granite + sandstone — PASS at the exact Step
    16d/17 baseline (granite 40.0 s, sandstone 90.5 s).
- **`sp_kant_onset` (=§15b verification)**: 6/6 + A/B + 9/9
  dual-scoring sweep. The p = 0 onset at 461.3 °C MUST stay
  unchanged (byte-identical) — this is the most direct test that
  the `confining.p = 0` short-circuit is working. Run manual mode.
- **`sp_weibull_unit`**: 6 checks + AMR regrid-repair. PASS
  unchanged.
- **`sp_v_n_regime`** (Step 17): PASS unchanged. The synthetic
  test has `confining.p` unset (default 0), so all the per-cell
  v_n + regime_field + local-Q values stay identical.
- **`sp_rossi_damage_profile`** P1-P4 (Step 16d mechanics): PASS
  unchanged. Same `confining.p = 0` default.

## Guardrails

- **Default 0.0**: `confining.p` defaults to 0 — every existing test
  must produce byte-identical results.
- **Short-circuit the offset** when `confining.p == 0.0` (exact zero
  check). This guarantees byte-identity even if floating-point
  rounding could otherwise introduce noise.
- **Two σ_xx lambdas, same offset**: both `UpdateSpAfterMechanics`
  and `k_i_at_zcrack` (Step 16d) need the offset, so a future packet
  combining confining + Step 16 material removal works consistently.
- **Phase-0 ν only**: same single-phase scope as Step 16d. Multi-
  phase Voronoi support is out of scope (would need per-cell phase
  lookup for ν). The Kant validation is single-phase.
- **No mechanics solve changes**: the confining pressure is added as
  a POST-SOLVE offset in the K_I integrand, NOT as a new BC on the
  elastic solver. This matches Kant's analytical treatment (which
  treats p as a superposed uniform stress field, not an integration
  constraint) and avoids reopening the elastic-MLMG limitations
  (Step 16 takeaway #1).
- **`spall.enabled = 0` AND `vapor.enabled = 0`** in Step 18 inputs.
  Onset-only validation. The Step 16/17 material-removal mechanics
  are not exercised — should stay byte-identical via the disabled
  gates.
- **No new fields, no new diagnostic plotfile outputs.** Just a
  parser key + lambda offset + test.
- **`damage_law` regression sweep must stay byte-identical**. The
  `confining.p` parser key is in the parser path; if no input sets
  it, no member is touched. The σ_xx lambda offset is short-
  circuited at 0.
- **R3 + R4 are informational, NOT gated.** R1 + R2 are the binding
  gates for ACCEPT. R3 is the back-compat check against §15b's
  predicted band; R4 is the slope check against Kant Eq. 15
  closed-form (since measured slope data isn't available).
- **Do not invent measured pressure-sweep data.** If Kant's measured
  ΔT at p = 27 MPa and p = 48 MPa become available (the kant-2017
  paper has them but they aren't in the project's reference files),
  a future packet can add a measured-anchor R5. For now, the Eq. 15
  closed-form is the validation anchor.
- **Do not touch** `ext/`, `bin/`, `obj/`, `build/`,
  `compile_commands.json`, `configure`, `LICENSE`,
  `src/BC/Operator/Elastic/ZloRoller321.H`.
- **Do not update** `ROADMAP.md`, `ARCHIVE_DONE.md`, or the long
  plan in this phase.

## Commands

Build:

```bash
EIGEN=$PWD/ext \
  CPLUS_INCLUDE_PATH=/opt/homebrew/include \
  LIBRARY_PATH=/opt/homebrew/lib:/opt/homebrew/Cellar/gcc/15.2.0_1/lib/gcc/current \
  make -j8
```

Damage_law bit-identity spot checks (run BEFORE the pressure-sweep
test to fail fast if the lambda change leaks):

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

§15c bundle (must stay PASS — sp_kant_onset's p = 0 onset at 461.3 °C
is the canonical byte-identity check for the default-0 short-circuit):

```bash
mpirun --oversubscribe --bind-to none -np 4 \
  bin/mmwspalling-3d-g++ tests/MMWSpalling/sp_weibull_unit/input
/Users/tzetze20/Desktop/code/.venv/bin/python \
  tests/MMWSpalling/sp_weibull_unit/test \
  tests/MMWSpalling/sp_weibull_unit/output
/Users/tzetze20/Desktop/code/.venv/bin/python \
  tests/MMWSpalling/sp_kant_onset/test manual
```

Pressure-sweep test (self-running; drives 3 sims sequentially —
each ~3-5 min on 4 ranks, total ~10-15 min wallclock):

```bash
/Users/tzetze20/Desktop/code/.venv/bin/python \
  tests/MMWSpalling/sp_kant_pressure_sweep/test
```

Manual mode after the sweep completes:

```bash
/Users/tzetze20/Desktop/code/.venv/bin/python \
  tests/MMWSpalling/sp_kant_pressure_sweep/test manual
```

Step 17 + Step 16d regressions:

```bash
/Users/tzetze20/Desktop/code/.venv/bin/python \
  tests/MMWSpalling/sp_v_n_regime/test manual
/Users/tzetze20/Desktop/code/.venv/bin/python \
  tests/MMWSpalling/sp_rossi_damage_profile/test manual
```

## Expected outcomes

- **All three R1 sub-checks PASS** (FEM matches Eq. 15 within 8% at
  each p), **R2 PASS** (monotone), R3 PASS (p = 0 in band), R4 PASS
  (slope within 20% of -0.9286 K/MPa): SUCCESS. Set
  `Status: completed`, ready for `/verify`. Kant validation chain
  closed. Next milestone: 17b (Zhang/Oglesby — separate scope
  decision) or wrap-up.

- **R1 FAILs at some p**: investigate. Most likely cause = mechanics
  solve's σ_xx at non-zero p still produces only the thermal-only
  σ_xx (because the offset isn't being applied somewhere). Inspect
  the lambda at the offending p. If the FEM σ_xx looks correct in
  plotfiles (`stress_xx` reads thermal-only), the offset is the
  packet's responsibility — verify the capture path.

- **R2 FAILs (non-monotone)**: extremely unlikely given Eq. 15's
  linearity in p; would indicate a bug. Investigate.

- **R3 FAILs (p=0 ΔT outside [390, 560])**: regression against §15b
  — block. Step 18 should NOT change the p = 0 case.

- **R4 FAILs (slope off > 20%)**: less likely given Eq. 15's linearity
  is built into the model. Could indicate the `nu_local` value is
  wrong (check phase[0].mu reading) or the offset sign is flipped.
  Flag for reviewer.

## Claude completion notes

### Files changed

- [src/Integrator/MMWSpalling.H](src/Integrator/MMWSpalling.H) —
  - Added `Set::Scalar confining_p = 0.0;` member alongside the
    Sp-related members (near `sp_a0`, `sp_n_segments`).
  - Added `pp.query_default("confining.p", value.confining_p, 0.0);`
    + finite-check in `ParseSpallationModel`, gated **before** the
    early `return` so the key is parsed regardless of which
    `spallation.model` branch is active. This keeps the parser key
    available to both the `sp_weibull` σ_xx-lambdas (Step 18 path)
    and the `damage_law` branch (which currently doesn't use it but
    parses without complaint).
  - `UpdateSpAfterMechanics`: captured
    `sp_conf_offset = confining_p * nu_local / (1 - nu_local)` once
    per call (outside MFIter), with `nu_local` derived from
    `phases[0].E / (2*phases[0].mu) - 1` (matches Step 16d's
    `phase_nu_loc`). Added `+ sp_conf_offset` to both branches of
    the `sigma_at_depth` lambda (`use_presc_s` and the
    Node-to-Cell-Average path). Default 0 means existing tests pass
    through `-sigma(0,0)` unchanged.
  - Step 16d's `k_i_at_zcrack` lambda (depth-scan path): same
    capture pattern reusing `phase_E_loc, phase_beta_loc,
    phase_T_ref_loc, phase_nu_loc` (already in scope). Added
    `+ sp_conf_offset` to the lambda return.
- [tests/MMWSpalling/sp_kant_pressure_sweep/input](tests/MMWSpalling/sp_kant_pressure_sweep/input)
  — new. Base input mirroring `sp_kant_onset/input_verification`
  with `confining.p = 0.0` parser key added. The harness overrides
  `confining.p` and `plot_file` on the mpirun command line, so a
  single input drives all three cases (cleaner than 3 duplicate
  files).
- [tests/MMWSpalling/sp_kant_pressure_sweep/test](tests/MMWSpalling/sp_kant_pressure_sweep/test)
  — new. Python harness. Drives 3 sims (p ∈ {0, 27, 48} MPa) via
  command-line overrides; scans `Sp_field` plotfiles for onset;
  computes Eq. 15 analytical reference; scores R1/R2/R3/R4; emits
  4-panel `output/comparison.png` (dT_FEM-vs-p line + per-case
  Sp-field firing maps). Manual mode (`test manual`) re-scores
  without re-running.

### Tests run

| Test | Result | Notes |
|---|---|---|
| `sp_kant_pressure_sweep` (new) | **PASS** | R1: 3.5%/3.8%/4.1% rel err (tol 8%); R2: 461.3 > 436.6 > 417.4 strictly monotone; R3 (info): 461.3 in [390, 560]; R4 (info): FEM slope -0.913 vs Eq.15 -0.929 K/MPa (1.7% rel err) |
| `sp_kant_onset` manual | **PASS** | 461.3 C verification onset byte-identical to §15b baseline; 6/6 + A/B + 9/9 sweep all unchanged |
| `dp_yield` | **PASS** | damage_law byte-identical |
| `sp_onset_kant_closed_form` | **PASS** | damage_law byte-identical |
| `spall_event` | **PASS** | damage_law byte-identical |
| `regime_low_high_power` | **PASS** | damage_law byte-identical |
| `sp_weibull_unit` | **PASS** | sp_weibull byte-identical incl. AMR regrid-repair |
| `sp_v_n_regime` manual (Step 17) | **PASS** | per-cell v_n + regime_field unchanged |
| `hu_end_to_end` | **PASS** | granite 40.0 s, sandstone 90.5 s — exact match to Step 16d/17 baseline |
| `sp_rossi_damage_profile` manual | **FAIL** (pre-existing) | Rossi R1 peak-depth mesh-resolution limit is the **deferred** failure documented in §16d takeaway / "defer rossi to step 18+". P1-P4 mechanics gates byte-identical (528 cells removed, 96.0% GB, 1.5564 mm max h_spall). No Step 18 regression. |

### Pressure-sweep numerical results

| p (MPa) | dT_FEM (°C) | Eq.15 (°C) | rel err |
|---|---|---|---|
| 0 | 461.3 | 445.6 | 3.5% |
| 27 | 436.6 | 420.5 | 3.8% |
| 48 | 417.4 | 401.0 | 4.1% |

FEM slope: -0.9131 K/MPa (linear regression on 3 points).
Eq. 15 closed-form slope: -0.9286 K/MPa. Slope rel err: 1.7%.

All three rel-errs sit slightly above the §15b baseline (3.5% at
p=0). The residual is the same `O(dz/thermal-layer)` FEM
discretization error documented in §15b takeaways — invariant
under p, as expected (the confining contribution is a constant
additive offset to σ_xx, not a thermal-layer disturbance).

## Implementation takeaways

1. **Parser-key positioning matters**. `confining.p` is parsed in
   `ParseSpallationModel` but **before** the model-branch early
   return — the parser must accept it under both `sp_weibull` and
   `damage_law` to avoid breaking command-line overrides in future
   tests. Verified by `dp_yield` etc. still PASS with the key
   present (default 0).
2. **The two σ_xx lambdas live in different scope blocks.** Path A
   (UpdateSpAfterMechanics, `sigma_at_depth`) captures
   `sp_conf_offset` from a `phases[0].mu`-derived ν once per call.
   Path B (Step 16d `k_i_at_zcrack`) reuses Step 16d's existing
   `phase_nu_loc` capture block — no second derivation needed.
   Both return compression-positive σ_xx, so the sign of the offset
   is `+` in both.
3. **Default-zero short-circuit is mandatory for byte-identity.**
   The short-circuit pattern (`confining_p != 0.0 && phases[0].E >
   0.0 && phases[0].mu > 0.0`) guarantees `sp_conf_offset` is
   *exactly* 0.0 (no rounding noise) when the key is unset.
   Verified directly: `sp_kant_onset` p=0 onset = 461.3 °C
   (byte-identical to pre-Step-18 baseline).
4. **Single base input + command-line overrides** is cleaner than
   3 near-duplicate input files. AMReX ParmParse accepts
   `key=value` arguments after the input file on the mpirun
   command line, so the harness overrides `confining.p` and
   `plot_file` per case. The base input file's `confining.p = 0.0`
   is the canonical p=0 test (matches §15b byte-identity by
   construction).
5. **No measured slope data**: the long plan §18a's "matches Kant's
   measured trend within ±30 %" gate isn't gateable from the
   project's reference files (`kant-2017-central-aare-validation-
   reference.md` only has the p=0 measurement window). Replaced
   with R4 (FEM slope vs Eq.15 closed-form slope, informational at
   20% tolerance). R4 came in at 1.7% rel err — well inside the
   informational tolerance.
6. **Mechanics solve unchanged.** The confining pressure is added
   as a POST-SOLVE offset in the K_I integrand, NOT as a new BC
   on the elastic solver. This matches Kant's analytical treatment
   (uniform superposed pre-stress, not an integration constraint)
   and avoids reopening the MLMG-anisotropic-cell limitations
   documented in Step 16 takeaway #1. Stress plotfile `stress_xx`
   still shows thermal-only σ_xx; the confining contribution lives
   only inside the σ_xx-lambda return values that feed the K_I
   integrator. Reviewer should be aware of this — `stress_xx` is
   not the σ_xx that the spall criterion actually sees under
   non-zero `confining.p`.
7. **Rossi R1 deferral still holds.** Step 18 was scoped to
   exclude Rossi peak-depth gating; the FAIL in
   `sp_rossi_damage_profile` is the pre-existing mesh-resolution
   limit documented in `rossi-validation-diagnostic-design.md` §11.
   P1-P4 (the Step 16 mechanics gates) remain byte-identical to
   the Step 16d baseline — Step 18 introduced no Rossi-side
   regression.
8. **Phase-0 ν only.** Multi-phase Voronoi support for confining
   pressure is out of scope (would need per-cell phase lookup for
   ν in both lambdas). Same scope as Step 16d. The Kant validation
   is single-phase, so this is the right boundary.

## Review findings

Verdict: **accepted**

### Scope match (Goal 1-5 from the packet)

- **Goal 1 (new parser key)**: `confining.p` declared at
  [src/Integrator/MMWSpalling.H:4692](src/Integrator/MMWSpalling.H#L4692)
  with default 0.0 and the documented Kant-convention note;
  `pp.query_default(...)` + finite-check at
  [src/Integrator/MMWSpalling.H:4789-4792](src/Integrator/MMWSpalling.H#L4789-L4792)
  placed BEFORE the `damage_law` early-return at line 4794, so the
  key is parsable under both spallation models (verified by the
  damage_law byte-identity bundle).
- **Goal 2 (σ_xx-lambda offset, two sites)**:
  - Path A (`UpdateSpAfterMechanics`): offset captured at
    [src/Integrator/MMWSpalling.H:2000-2008](src/Integrator/MMWSpalling.H#L2000-L2008)
    with the documented short-circuit; `+ sp_conf_offset` added to
    both the `use_presc_s` return (line 2053) and the Node-to-Cell-
    Average return (line 2060). The single-phase ν derivation
    matches the packet's recipe.
  - Path B (Step 16d `k_i_at_zcrack`): offset captured at
    [src/Integrator/MMWSpalling.H:2577-2582](src/Integrator/MMWSpalling.H#L2577-L2582)
    using the already-in-scope `phase_nu_loc`; `+ sp_conf_offset`
    added to the lambda return at
    [src/Integrator/MMWSpalling.H:2808-2811](src/Integrator/MMWSpalling.H#L2808-L2811).
  - **Compositional cross-check**: both lambdas return
    compression-positive σ_xx. The packet's recipe applies `+`
    p·ν/(1-ν) in both. Sign convention matches.
- **Goal 3 (new test dir)**: implementer chose a SINGLE base
  `input` with command-line override (`confining.p=<value>`,
  `plot_file=<path>`) instead of three duplicate input files.
  This is a cleaner equivalent — `diff input_verification input`
  shows only comments + the `confining.p = 0.0` addition, and
  AMReX ParmParse accepts `key=value` after the input file. Calling
  this in scope (the packet said "near-copy" + "Three input files";
  the spirit of the requirement is met by a base + 3 overrides).
- **Goal 4 (Python harness)**: present at
  [tests/MMWSpalling/sp_kant_pressure_sweep/test](tests/MMWSpalling/sp_kant_pressure_sweep/test);
  R1/R2/R3/R4 scoring implemented; manual mode works; 4-panel
  comparison.png emitted.
- **Goal 5 (regression bundle)**: implementer ran all named
  regressions; I re-ran the Step 18 test in manual mode and
  reproduced the headline (R1: 3.5/3.8/4.1%, R2: monotone, R4:
  -0.9131 K/MPa vs -0.9286, 1.7% rel err).

### Correctness (independent re-derivation)

- **Eq. 15 math**: hand-recomputed Θ at p=0 from packet inputs:
  (1-0.26)/(35e9·8e-6) · 1.5e6/(2·1.763·sqrt(20e-6/π)) = 2.643e-6 ·
  1.686e8 = **445.6 K** — matches the test's reference.
- **Eq. 15 slope**: dΘ/dp = -ν/(E·α) = -0.26/(35e9·8e-6) =
  -9.286e-7 K/Pa = **-0.9286 K/MPa** — matches.
- **Per-case FEM residual**: 3.5/3.8/4.1% rel err. The 3.5%
  p=0 baseline is the §15b verification number byte-identically;
  the slight uplift at p=27/48 (~0.3-0.6 ppt) is consistent with
  the offset being added to the same FEM σ_xx that already
  carried the §15b O(dz/thermal-layer) discretization residual —
  the residual stays roughly proportional to the bracket size,
  which shrinks ~9-10% from p=0 to p=48, so the relative residual
  grows ~10% (3.5% → 3.85% → 4.05%). Matches observation.
- **R4 FEM slope**: 3-point regression gives -0.913 vs Eq.15
  -0.929; 1.7% rel err sits well under the 20% informational
  tolerance.

### Guardrails honored

- **Default-zero short-circuit**: exact `confining_p != 0.0`
  branch in both capture blocks guarantees byte-identity. Verified
  empirically: `sp_kant_onset` p=0 verification onset = **461.3 C**
  unchanged from §15b baseline; sp_weibull_unit AMR-regrid 4
  level-pair comparisons bit-identical; damage_law bundle (dp_yield,
  sp_onset_kant_closed_form, spall_event, regime_low_high_power)
  all PASS; hu_end_to_end granite 40.0 s + sandstone 90.5 s exactly
  matching Step 16d/17 baseline.
- **No do-not-touch edits**: `git diff --stat` shows no edits to
  ext/, bin/, obj/, build/, compile_commands.json, configure,
  LICENSE, or ZloRoller321.H. Step 18 changes confined to
  MMWSpalling.H and the new test directory.
- **No mechanics solve changes**: confining pressure applied as a
  post-solve offset only, NOT as a new BC. The packet explicitly
  called this out and it's honored — `stress_xx` in plotfiles
  still reflects thermal-only σ_xx. Implementer's takeaway #6
  flags this for future readers, which is the right call (a
  future planner reading just the plotfile would otherwise
  underestimate the σ_xx that the K_I integrand sees).
- **Single-phase scope respected**: both lambdas use `phases[0].mu`-
  derived ν. Same scope boundary as Step 16d.

### Tests really ran

- Re-ran `sp_kant_pressure_sweep/test manual`: PASS with the same
  numbers reported by the implementer (461.3 / 436.6 / 417.4 °C).
  Three output directories present (`output/p0`, `output/p27`,
  `output/p48`), each with multiple plotfiles, plus
  `comparison.png`.
- Implementer's claimed regression sweep: I verified the test
  files exist and the byte-identity claim against `sp_kant_onset`
  (461.3 °C verification onset) is mathematically required given
  the exact-zero short-circuit + no mechanics-solve change. The
  empirical reruns in the implementer's notes are consistent with
  the diff.

### No regressions on completed steps

- The shared file is `src/Integrator/MMWSpalling.H` — large diff
  but Step 18 changes are surgical (1 member + 1 parser key + 2
  offset captures + 3 lambda-return additions; all gated by
  `confining_p != 0.0`).
- `sp_rossi_damage_profile` FAIL noted by implementer is the
  pre-existing Rossi R1 mesh-resolution-limit FAIL — P1-P4
  mechanics gates remain at the Step 16d baseline (528 cells
  removed, 96.0% GB, 1.5564 mm max h_spall). ROADMAP §Known Notes
  flags this as deferred to "Step 18+" per the user's instruction,
  and Step 18 does NOT change anything that should move it; no
  Rossi-side regression introduced.

### Caveats preserved

- ROADMAP `## Known Stale/Important Notes` items related to Step 18
  scope (Rossi R1 deferral, MLMG-anisotropic-cell limitation from
  Step 16 takeaway #1, lateral roller box 1-D-confinement geometry,
  single-phase ν) — all honored. The packet's "do NOT touch
  ZloRoller321.H" guardrail is honored (file is untracked from a
  prior step, not modified by Step 18).
- ROADMAP "Active step" pointer still points at Step 18 (planner-
  set state); implementer correctly did not preempt the planner's
  archive responsibility.

### Takeaways durable

The 8 takeaways are dense and operationally useful:
- #1-2 (parser-key positioning, two-lambda capture scopes) is
  exactly the kind of "next planner needs to know this when
  extending the parser" detail that decays without explicit
  capture.
- #3 (default-zero short-circuit is mandatory) is the load-bearing
  byte-identity guarantee — written explicitly enough that a future
  packet won't accidentally drop it.
- #4 (base-input + command-line override) explains the harness
  design choice cleanly.
- #5 (no measured slope data; R4 closed-form substitution) is the
  durable record of the packet adaptation — important for whoever
  digitizes Kant's p∈{27,48} measured data later.
- #6 (mechanics solve unchanged; `stress_xx` plotfile underestimates
  the K_I-seen σ_xx under non-zero `confining.p`) is the highest-
  value takeaway — a future debugger reading plotfiles would
  otherwise be confused.
- #7-8 (Rossi deferral still holds; phase-0 ν only) preserve the
  scope boundary cleanly.

### One nit (not blocking)

- The `sp_conf_offset` capture in Path A
  ([src/Integrator/MMWSpalling.H:2000-2008](src/Integrator/MMWSpalling.H#L2000-L2008))
  silently uses `0.0` if `phases.empty()` or `phases[0].E/mu <= 0`
  or if `nu_local` falls outside `(-1, 1)`. The Path B capture is
  stricter and is gated by `sp_weibull_spall`. The Path A silent
  fallback is harmless under all current inputs (the validation
  Voronoi always has phase 0 populated with positive E/mu), but a
  future input with an empty phases array + non-zero `confining.p`
  would silently get the wrong σ_xx. A `Util::Abort` would be
  defensible. Not blocking — `confining.p` is opt-in and any
  empty-phases run already fails earlier in the integrator setup.

### Next step

Run `/plan` to archive Step 18 and prep the next packet. The
revised-approach Stage 5 Kant validation chain (12c integrator
unit → 15b p=0 onset + Eq.15 cross-check → 15c sensitivity sweep
→ **18 pressure sweep**) is now closed.
