---
title: Revised Modeling Approach — Sp + Weibull on Voronoi (replacing the eq. 27 damage chain)
type: analysis
tags: [model-design, spallability, weibull, voronoi, kant-2017, walsh-2013, main-tex-revision]
created: 2026-05-07
updated: 2026-05-07
sources: [[[sources/kant-2017]], [[sources/walsh-lomov-2013]], [[sources/rossi-2018]], [[sources/liu-2024-spallation]], [[sources/zhang-2023-mmw-vaporization]], [[sources/nasseri-2007]], [[sources/nasseri-2009]], [[sources/friedrich-wong-1986]]]
---

# Revised modeling approach for [[analysis/main.tex]]

This document captures the planned revision to the damage / spallation block of the unified MMW drilling model, before edits to main.tex. The motivation is the user's flag that eq. 27's continuous-damage law has too many uncertain free parameters (`A_D, n, m, Ψ_0`). The revision replaces that block with an LEFM-based criterion ([[sources/kant-2017]]) coupled to a Weibull flaw-length distribution on the Voronoi grid (after [[sources/walsh-lomov-2013]]).

## Scope: what stays vs what changes

**Kept (no change):**
- 3D simulation domain (heterogeneity is a 3D effect; 2D axisymmetric suppresses grain-scale stress concentrations).
- AMReX / ALAMO framework.
- Voronoi microstructure with mineral phases assigned per cell.
- Static mechanical solve (operator-split from thermal).
- Enthalpy-based heat conservation with Beer–Lambert volumetric source.
- Beer–Lambert MMW absorption with `α(T)` and Gaussian beam.
- Drucker–Prager *as a stress criterion for plastic / shear yield* (separate from the spall criterion — the two cover different physics).

**Replaced:**
- Eq. 27 continuous-damage rate law `∂D/∂t = A_D (Ψ/(1+Ψ_0))ⁿ (1−D)ᵐ` → **Sp-based onset criterion + Weibull flaw distribution on grain boundaries** (this document, §1–§3).
- Eq. 36 spall-thickness scaling `h_spall = C_h √(α·t_spall)` → **Sp-vs-depth zero-crossing** (§4).
- Eq. 39 `RoP = max(RoP_spall, RoP_vap)` → **per-cell surface velocity**: simultaneous local spall + vap (§5).
- Mie–Grüneisen EOS for `E(P)` → **kept as a note**, replaced in numerics by linear `E(P) = E_0 + (∂E/∂P)·P` (§7). The full EOS adds little at lithostatic pressures relevant for ≤10 km drilling but the formulation in main.tex remains as a reference.

**Dropped:**
- Cohesive-zone interface elements at grain boundaries. Grain-boundary mechanics enters statistically, through a Weibull flaw-length distribution biased to grain-boundary cells (§3), justified empirically by [[validation/rossi-2018__central-aare-granite-strength-vs-T-rate]] (ρ_B/ρ_I = 5.8 under flame).

## 1. New onset criterion: Kant 2017 Sp number, generalized to volumetric heating

[[sources/kant-2017]] derives the spalling criterion `Sp = K_I / K_Ic ≥ 1` for a plane-strain edge crack of length `a` in a half-space. The closed-form decomposition `Sp_red = Sp_rock · Sp_fluid + Sp_conf` is specific to **convective surface heating** (Carslaw–Jaeger temperature distribution). The user's model uses **volumetric Beer–Lambert deposition**, so the closed-form simplification does not transfer.

What does transfer is the underlying LEFM relation:

```
K_I(x,y,t) = 2 √(a/π) · ∫₀¹ σ_xx(x, y, ξa, t) · f_K(ξ) dξ
```

with `f_K(ξ) = (1.3 − 0.3 ξ^(5/4)) / √(1 − ξ²)` (Tada-et-al. weight function for an edge crack in a half-space; Kant Eq. 6).

Onset condition: `K_I(x,y,t) ≥ K_Ic(T(x,y,0,t))` evaluated locally over the surface.

**Key adaptation steps:**

1. The stress field `σ_xx` is computed by the model's static mechanical solve (kept) from the actual MMW-driven temperature field — *not* approximated by Kant's `σ_xx = −EαΘ/(1−ν) − pν/(1−ν)` analytical form. The integral above uses whatever σ_xx the FEM gives.
2. `K_Ic(T)` enters at the local surface temperature. Direct measurement values from [[validation/nasseri-2007__westerly-kic-vs-T]] for Westerly granite: 1.43 (RT) → 1.35 (250 °C) → 0.98 (450 °C) → 0.43 (650 °C) → 0.22 (850 °C) MPa·m^0.5 — concave curve with the steepest drop bracketing the quartz α–β transition (573 °C). The linear `K_Ic(T) = K_Ic⁰·(1 − (T−T_ref)/T_melt)` form in main.tex eq. 30 does not fit this shape; a piecewise interpolation through the table or a sigmoid centered near 573 °C is more honest. For onset prediction (rock first reaches spalling temperature ~400–600 °C), `K_Ic ≈ 0.6–1.2 MPa·m^0.5` brackets the value to use. **Open assumption**: Nasseri's measurement is at slow heating rate; whether the same `K_Ic(T)` applies under MMW gradient-concentration heating is unverified.
3. Confining pressure enters the *stress field* directly via the in-situ stress component on the model's mechanical BCs, *not* through a separate `Sp_conf` term. This is cleaner — one stress field, one integral.

**Why this is parameter-light.** Inputs: `E, α, ν, K_Ic, λ, κ` (all measurable); plus `a` (handled stochastically; §3); plus the MMW source parameters which are independently set. **Zero fitted rate constants.** `A_D, n, m, Ψ_0`, `C_h` all gone.

## 2. Drucker–Prager stays — but for a different role

DP yield (eq. 11–12 in main.tex) **is not** the spall criterion under this revision. It is the criterion for **plastic / shear yield** at depth, where confining pressure suppresses brittle spall and the rock yields plastically (Byerlee 1968 brittle–ductile transition). The hierarchy in §3.4.4 of main.tex collapses to:

- If `K_I(x,y,t) ≥ K_Ic(T)` over a connected surface zone (§4) → **spall detach**.
- Else if `F_DP(σ) ≥ 0` AND surface T > brittle-ductile transition T → **plastic yield / no detachment**, allow plastic strain to develop.
- Else if `H ≥ H^V` at surface → **vaporization** (§5).

The Hoek–Brown spalling-ratio cutoff is dropped; it was acting as a sanity filter on the damage path that no longer exists.

## 3. Weibull flaw-length distribution on grain boundaries

The Sp criterion's strongest sensitivity is to flaw length `a` (a^(3/2) in Sp_rock; Kant 2017 Fig. 4). A single deterministic `a` gives a uniform onset prediction across the surface and misses the experimentally observed **patchy, grain-scale** spall pattern. Following [[sources/walsh-lomov-2013]]'s use of statistical critical flaws and the [[sources/rossi-2018]] empirical result that flame-rate heating concentrates damage 5.8× at grain boundaries:

**Implementation (decided):**
- For each Voronoi face that lies on a grain boundary, assign a flaw length `a_f ~ Weibull(a_0, m_W)`.
- Intragranular faces: assign `a_f` from the same Weibull but with characteristic length reduced by a factor (~5×, mirroring the ρ_I/ρ_B partition under fast heating).
- Sp evaluated per-face using the local `a_f`. Surface "spall map" = locus where `K_I(a_f, σ_xx, T) ≥ K_Ic`.

**Parameter values for Stage-1 (Westerly granite, after [[sources/nasseri-2007]] Table 4):**
- `a_0` (characteristic boundary flaw length, RT): **0.26 mm** (mean GB crack length in untreated Westerly granite). For intragranular faces, characteristic length **0.13 mm** (mean IG crack length).
- These are an order of magnitude larger than Kant 2017's 20 µm for Central Aare granite — different rocks have different flaw populations. Per-rock parameterization is necessary.
- `m_W` (Weibull modulus): granite tensile-strength literature gives m = 8–25. [[sources/lyu-2018-mineralogy]] used m = 25 for the Rauenzahn–Tester ROP formula on granite. **Default `m_W = 15`** as a midpoint, with sensitivity sweep over [8, 25] as part of the parametric study.

**Bracket from [[sources/friedrich-wong-1986]]** for sanity-check: virgin Westerly granite has a maximum crack length ~600 µm with arithmetic mean ~18 µm (Hadley 1976, cited in Friedrich-Wong). The Nasseri 2007 mean values (0.13 and 0.26 mm) sit well above the 18 µm Hadley figure — reflecting that Nasseri's "average length" is over the full population while Hadley's "mean" refers to actively propagating cracks. The Weibull fit should reproduce both.

**This is the one place the model retains parameters that are not directly measured at the application scale.** They are, however, well-defined experimentally — a uniaxial tensile-strength Weibull fit on the rock specimen gives `(a_0, m_W)` directly. Two parameters, both standard rock-mechanics measurements; vs four phenomenological constants in eq. 27.

## 4. Spall thickness from Sp-vs-depth

For a surface point where the Sp onset condition is satisfied, the spall thickness is set by the depth at which the LEFM driving force falls below the toughness:

```
h_spall(x,y) = max{ z : K_I(x, y, z, t; a) ≥ K_Ic(T(x,y,z,t)) }
```

i.e., evaluate the Tada-weight integral for hypothetical edge cracks at increasing depth and find where it falls below threshold. This is **internally consistent** with the onset criterion (one criterion, not two with different scalings), and **parameter-free** beyond what's already in §1 and §3.

Implementation note: this only needs to be evaluated once per spall event, not every step. The thermal-gradient profile gives the dominant z-dependence; computing K_I(z) over a vertical line of depths at that (x,y) is cheap.

For a sanity check against existing scaling: at the diffusion-dominated limit, this should recover `h_spall ~ √(κ·t_spall)` because that's the only length scale in the temperature profile — but with a coefficient that's now derived, not fitted to `C_h`.

**Additional sanity check from [[sources/friedrich-wong-1986]] threshold-T(a) curves:** in the limit of zero thermal gradient (uniform-T expansion-stress regime), the Sp model with Beer-Lambert source set to zero attenuation should produce onset T as a function of flaw length matching Friedrich-Wong's Fig. 15a for Westerly granite (e.g. a = 10 µm → T_onset ≈ 120 °C GB / 180 °C IG). This is a self-consistency test for the Sp implementation independent of any MMW-specific feature.

**Additional empirical anchor from [[sources/rossi-2018]]:** at flame conditions (650 °C / 5 °C/s) on Central Aare granite, the measured damage-zone depth is 100–500 µm (peak crack density at 100–200 µm). The Stage-3 implementation should reproduce this depth bracket within ~30%.

## 5. Spall + vaporization happen simultaneously on the surface

User confirmed eq. 39 `RoP = max(RoP_spall, RoP_vap)` is wrong — both modes operate at different points on the surface at the same time. Replace with a per-cell surface velocity:

```
v_n(x,y,t) = v_spall(x,y,t)  if Sp criterion fires this step
           = v_vap(x,y,t)    if H(x,y,0,t) ≥ H^V
           = 0               otherwise
```

with:
- `v_spall · Δt = h_spall(x,y)` from §4 (instantaneous detachment of an h_spall-thick patch above the firing surface point in the timestep where Sp first crosses 1).
- `v_vap` from the existing Hertz-Langmuir / enthalpy-based recession in main.tex (eq. 47 with `RoP_vap` interpreted locally rather than averaged).

Bulk RoP: `RoP_bulk = (∫_surface v_n dA) / A_beam`.

This naturally handles the regime continuum: at low power, most of the beam footprint is in the spall mode; at very high power, the central hot spot is vap-mode while the cooler rim is spall-mode; in between, both operate on different patches simultaneously. The ratio of vap-area to spall-area becomes a useful diagnostic for the regime map (planned result in main.tex §5).

## 6. Two cracking regimes (heating-rate awareness)

[[sources/rossi-2018]] establishes that **MMW heating rates sit firmly in the "gradient-concentration" regime** (vs the "differential-expansion" regime of slow oven heating). [[sources/friedrich-wong-1986]] is the canonical reference for the slow-rate / expansion-stress regime; its data sit on the opposite end of the rate spectrum from Rossi's flame data. Implications for the model:

- The σ_xx field driving K_I in §1 should be computed from the **transient gradient-driven thermal stress**, not from a uniform-heating differential-expansion calculation. This is what main.tex already does, so no change — but worth noting why.
- The grain-boundary bias of the Weibull distribution (§3, factor ~5× shorter `a` on intragranular faces) is empirically grounded by Rossi's ρ_B/ρ_I = 5.8 — not just a modeling convenience.
- **Friedrich-Wong's `S_v − S_v⁰ = k·(T−T₀)²` quadratic damage law** with `k = E·(Δα)²/[8·(1−ν)·G_Ic]` is a **back-of-envelope sanity check** on the model's predicted total damage volume in the slow-rate limit. All four parameters in `k` are measurable; `k ≈ 0.0025 mm⁻¹/°C²` for Westerly granite. If the user's full FEM in slow-rate mode predicts total `S_v(T)` deviating from this quadratic by >2× across the 100–600 °C range, something in the implementation is wrong. The quadratic does not apply directly under MMW (different regime), but the sub-model is useful for code verification.
- Both [[sources/nasseri-2007]] (slow-rate post-treatment) and [[sources/friedrich-wong-1986]] (slow-rate post-treatment) measure damage on cooled samples. The model's own predictions of *post-cooling* damage are the right comparison point for those datasets — not the in-situ damage at peak T.

## 7. EOS — kept as a note, simplified in numerics

The Mie–Grüneisen formulation (main.tex §3.3, eq. 13–15) is retained as a reference but not implemented as such. At lithostatic pressures relevant for ≤10 km drilling (≤300 MPa for granite), `E(P)` changes by a few percent in solid rock, well within the property uncertainty. A linear

```
E(P) = E_0 + (∂E/∂P) · P
```

with `(∂E/∂P) ≈ 5–10` from standard granite measurements gives the same accuracy with one parameter. The EOS section in main.tex stays as written for completeness; the numerical implementation uses the linear form.

Confining pressure enters the *stress field* through the mechanical BCs (already in main.tex), not through a modulus correction. Kant 2017's `Sp_conf` term is implicitly captured this way — `p` in the static stress field shows up in the σ_xx integrand naturally.

## 8. Summary of parameter accounting

| Parameter | Old (eq. 27, 30, 36 chain) | New (Sp + Weibull) | Source / value |
|-----------|-----------------------------|---------------------|------------------|
| `A_D` (damage rate) | required, uncalibrated | — gone — | — |
| `n` (damage exponent) | required, uncalibrated | — gone — | — |
| `m` (saturation exponent) | required, uncalibrated | — gone — | — |
| `Ψ_0` (initiation threshold) | required, calibrated against Hu f_b | — gone — | — |
| `C_h` (spall thickness coefficient) | required, fitted | — gone — | — |
| `K_Ic(T)` | implicit in Ψ_0 | **measured table** | [[sources/nasseri-2007]] — 1.43 → 0.22 MPa·m^0.5 over RT–850 °C, concave with quartz-transition knee |
| `a_0` (Weibull characteristic flaw) | — | **bracketed** | [[sources/nasseri-2007]] — 0.26 mm GB / 0.13 mm IG for Westerly granite at RT |
| `m_W` (Weibull modulus) | — | **default 15, sweep [8, 25]** | granite-tensile literature; [[sources/lyu-2018-mineralogy]] used 25 |
| `(∂E/∂P)` | full EOS (Hugoniot c_0, s, Γ) | one number, ~5–10 | granite handbook values |

Net: **5 fitted parameters out → 3 measurable parameters in (K_Ic, a_0, m_W)** plus the linear EOS slope. All three new parameters are now numerically bracketed for Westerly granite from this batch ingest; per-rock values needed for Central Aare granite.

## 9. Validation strategy under the revision

### Priority validation tests (recommended core set)

Three tests give the highest information per validation run, covering orthogonal axes of the model on three different rocks under three different heating regimes. Each is reproducible from a published, fully-specified experimental setup.

**Priority Test A — [[validation/kant-2017__central-aare-granite-spalling-temperature]]: onset criterion + confining pressure.**
- *What it tests*: the core Sp ≥ 1 onset criterion (§1); the confining-pressure entry through the stress field (§2 hierarchy + §7); the LEFM integral itself.
- *Why it's the right priority test for these axes*: Kant 2017 is the criterion's home turf, with a closed-form analytical prediction available as a cross-check; three discrete pressures (p ∈ {0, 27, 48} MPa) give a *trend*, not just a single number — important because confining-pressure entry is one of the more subtle changes vs main.tex's `Sp_conf` term; the property table is fully published.
- *Pass criteria*:
  - Predicted onset ΔT at p = 0 within the measured 553–694 °C bracket using T-dependent K_Ic from [[sources/nasseri-2007]].
  - Monotone-decreasing predicted onset ΔT across p ∈ {0, 27, 48} MPa, with relative slope matching the measured trend within ±30%.
  - As a code-verification cross-check: in the Bi < 0.5 limit with Beer-Lambert α set to mimic a surface flux, the FEM integral should reproduce Kant's closed-form `Sp_red = Sp_rock·Sp_fluid + Sp_conf` value within ~5%.
- *Gotcha*: Kant's measurement is post-first-spall (pyrometer can't catch the first one), so measured > predicted is expected; the comparison metric is the *trend with p*, not the absolute number.

**Priority Test B — [[validation/hu-2018__lrst-granite-sandstone]]: onset timing + Weibull statistics + cross-rock portability.**
- *What it tests*: onset *timing* (not just threshold T) under continuous flame heating; the Weibull `a`-distribution's role in driving spatially patchy first-spall events (§3); cross-rock parameter portability (granite + sandstone with the same modeling framework, different `a_0, m_W, K_Ic`).
- *Why it's the right priority test for these axes*: Hu 2018 is the **only published direct LRST measurement with both granite AND sandstone**, giving a two-rock check on the same physics; the onset *time* (31.5 s granite, 76.7 s sandstone at known heat flux ~7×10⁴ W/m²) tests the dynamics, not just the threshold — meaning it tests both the temperature evolution and the K_I-vs-time integral together; the spatial pattern of first-spall events should be patchy (grain-scale clusters), not an axisymmetric ring.
- *Pass criteria*:
  - Predicted onset time within ±20% of 31.5 s (granite) and 76.7 s (sandstone) at the published heat flux.
  - Surface T at onset matches LRST values 791 K (granite) and 842 K (sandstone) within ±10%.
  - First spall is *patchy* on the grain scale — not an annular or axisymmetric ring. Statistical signature: at least 5 distinct firing locations within the beam footprint at first onset, distributed irregularly.
  - Granite-vs-sandstone ordering correctly reproduced (granite spalls first, sandstone later) — cross-rock portability test.
- *Gotcha*: Hu's sandstone is high-porosity and high-quartz; the porosity-mediated stress relaxation noted in [[concepts/thermal-spallation]] may matter — flag if the model overpredicts sandstone spallability.

**Priority Test C — [[validation/rossi-2018__central-aare-granite-strength-vs-T-rate]]: spall thickness + heating-rate regime.**
- *What it tests*: the §4 Sp-vs-depth spall-thickness criterion; the gradient-concentration cracking regime under fast heating (§6); grain-boundary-vs-intragranular damage partitioning.
- *Why it's the right priority test for these axes*: Rossi 2018 gives a **spatial damage profile** (crack density vs depth, binned at 94 µm) under documented flame conditions on Central Aare granite — not just an integrated number. This is the only validation case in the wiki that directly probes spall thickness without going through a coupled ROP measurement. The same rock as Kant 2017 (Test A), so material parameters port directly between A and C.
- *Pass criteria*:
  - Predicted damage-zone depth peaks at 100–200 µm with falloff to baseline by ~520 µm (Stage-3 §4 prediction).
  - Crack-density partitioning ρ_B/ρ_I ≈ 5–6 in the damage zone (justifies the §3 grain-boundary bias factor).
  - Total damage zone ≈ 4.4–7.6 mm/mm² peak crack density.
- *Gotcha*: Rossi reports post-cooling damage. The model's predicted damage at peak T must be projected forward to post-cooling state for a fair comparison (some cracks close on cool-down, some new ones open from cooling-rate stresses). For a first-pass validation, the *depth profile shape* is more diagnostic than the absolute crack density.

### Coverage matrix for the priority set

| Model component | Test A (Kant) | Test B (Hu) | Test C (Rossi) |
|------------------|----------------|--------------|------------------|
| Sp ≥ 1 onset criterion (§1) | ✓ direct | ✓ direct | indirect |
| Weibull statistical `a` (§3) | indirect (single largest a) | ✓ direct (patchiness, timing dispersion) | indirect (damage-zone heterogeneity) |
| Spall thickness (§4) | — | — | ✓ direct |
| Per-cell spall + vap (§5) | — | — | — (covered by Test E supporting) |
| Heating-rate regime (§6) | flame (medium) | flame (medium) | flame (medium-fast) — within regime |
| Confining pressure (§2 + §7) | ✓ direct (3 pressures) | — | — |
| Cross-rock portability | — | ✓ direct (gr + ss) | — |
| Onset timing dynamics | indirect | ✓ direct | — |

Tests A + B + C together exercise every section of §1–§4 directly and §6 indirectly. The single coverage gap is §5 (per-cell simultaneous spall + vap) — handled by supporting Test E below.

### Supporting tests (run after the priority set)

- **Supporting Test D — [[validation/hu-2018__lrst-granite-sandstone]] FEM cross-check** (numerical-vs-numerical): also test against [[validation/hu-2018__numerical-model-validation]] subsurface-T data (8% / 15% error against thermocouples). Anchors the model's thermal field independently of the damage block.
- **Supporting Test E — [[validation/zhang-2023__oglesby-woskov-validation]]: full-system MMW under high power.** Validates the per-cell spall + vap implementation (§5) at conditions where both modes are active simultaneously. Pass: model agrees with Zhang in the high-power vap-dominated regime AND produces non-zero spall contribution at the lower-power steps. The spall contribution is the model's prediction (no published spall data at OW conditions), so this is one-sided validation.
- **Supporting Test F — [[validation/lyu-2018__spall-morphology]] + [[validation/walsh-2013__pore-fluid-spallation-mechanism]]: spall-size distribution.** Drives the §10 open question on connected-patch detachment area `A_crit`. Compare predicted spall-size histogram against granite's 5×2 mm disks + 5 µm micro-particles + log-normal Walsh distribution. The `A_crit` value that reproduces both modal size and distribution shape is the answer.
- **Supporting Test G — [[validation/friedrich-wong-1986__westerly-crack-density-and-thresholds]]: zero-gradient self-consistency.** Run the Sp model with Beer-Lambert α set to zero (uniform deposition); predicted onset T(a) should match Friedrich-Wong Fig. 15a within ~50 °C. Pure code-verification test, independent of any MMW-specific feature.

## 10. Decisions and remaining open questions

**Decided (2026-05-07):**

1. ✅ **Weibull on `a`** (not on σ_t). Sp's a^(3/2) sensitivity makes `a` the dominant stochastic axis. Mathematically equivalent to Weibull-on-σ_t when `K_Ic` is held fixed at a single value, but `a` is the better-tabulated axis in the literature (e.g. [[sources/nasseri-2007]] Table 4 gives mean GB and IG flaw lengths; uniaxial-tensile-strength Weibull fits also report this).
2. ✅ **Spall thickness from Sp-vs-depth zero-crossing** (§4). Internally consistent with the onset criterion (one criterion, no separate scaling). Alternatives (stress-zero-crossing, energy integral) are equivalent in the simplest cases but introduce a second tunable in pathological geometries.

**Still open — let validation data decide:**

3. **Connected-patch detachment area `A_crit`.** Sp ≥ 1 must hold over a connected surface zone of area > A_crit before a spall detaches. Two extremes:
   - **Per-face detachment**: any single face with Sp ≥ 1 produces a spall. Lets sub-grain flakes pop off; produces a spall-size distribution biased toward fines.
   - **Connected-cluster detachment**: A_crit ≈ (mean grain diameter)² ~ 0.5 mm² for Westerly. Suppresses sub-grain spalls; predicts spall sizes ≥ grain scale.
   - **Plan**: implement both as runtime options. Compare against [[validation/lyu-2018__spall-morphology]] (granite spalls 5 mm × 2 mm disks; 5 µm micro-particles; <1 µm cracks — i.e., a multi-scale distribution that suggests neither extreme alone is right) and [[validation/walsh-2013__pore-fluid-spallation-mechanism]] (log-normal spall-size distribution). The validation that reproduces both the modal spall size and the size-distribution shape sets `A_crit`.
4. **Whether to keep a parallel von Mises / DP "sanity check"** — i.e., flag in post-processing whenever `F_DP ≥ 0` somewhere on the surface, to catch cases where the rock is yielding plastically before Sp fires. Current plan: yes, log it; do not let it suppress Sp firing.

## 11. Outstanding lit-review gaps for this approach

**Status update (2026-05-07).** The two highest-priority gaps from the original draft of this document have been closed by the [[sources/nasseri-2007]], [[sources/nasseri-2009]], [[sources/friedrich-wong-1986]] ingest:

- ✅ **`K_Ic(T)` for granite** — closed by [[sources/nasseri-2007]]. Lookup values: 1.43 (RT) → 1.35 (250) → 0.98 (450) → 0.43 (650) → 0.22 (850) MPa·m^0.5 for Westerly granite. **Note: linear form in main.tex eq. 30 is the wrong functional shape** — actual curve is concave with a knee at the quartz transition. Piecewise or sigmoid recommended. For onset prediction at typical spalling T (~400–600 °C), use `K_Ic ≈ 0.6–1.2 MPa·m^0.5` from the table; for first-spall on undamaged rock, room-T `1.43` is appropriate.
- ✅ **Confining-pressure correction** — closed by [[sources/nasseri-2009]]. At drilling-relevant P_eff (≥15 MPa), pre-existing cracks close substantially. **Net implication for the Sp model: no `K_Ic(P_eff)` correction needed for onset on initially undamaged rock at the drilling face**. The Kant 2017 `Sp_conf` term applied directly to the stress field is sufficient.
- ✅ **Mineralogical mismatch parameters** — closed by [[sources/friedrich-wong-1986]] Appendix 2. Direct values for E, ν, α(matrix), α(plag), α(augite), Δα, Δα', L for Westerly granite + Frederick diabase + Oak Hall limestone. Can seed any Voronoi-mismatch extension.
- ✅ **Sanity-check sub-model** — closed by [[sources/friedrich-wong-1986]] `S_v ∝ (T−T₀)²` quadratic damage law with all-measurable parameters. Use as fast back-of-envelope check on full-FEM damage-volume output.

**Still open:**

- **Direct MMW-spallation onset measurement** is genuinely missing in the published literature ([[sources/oglesby-woskov-2014]] is too high-flux; [[sources/zhang-2023-mmw-vaporization]] / [[sources/gokhale-2026-multiphysics-mmw]] don't address spallation). This is a *validation gap*, not a *modeling gap*. The user's model output for MMW spallation onset is therefore a prediction, not a back-fit.
- **Granite Weibull `(a_0, m_W)` flaw-distribution parameters** — partially constrained by Nasseri 2007 crack-length data (avg IG length 0.13 mm, GB length 0.26 mm at RT for Westerly granite; 1 order of magnitude larger than Kant 2017's 20 µm for Central Aare granite). The full distribution is still unmeasured; a uniaxial-tensile Weibull fit on the specimen is the cleanest source. **Not blocking model implementation** — defensible default values are now bracketed.
- **Heating-rate-dependent K_Ic**: all ingested K_Ic data are post-heating, slow-rate. MMW (gradient-concentration regime per [[sources/rossi-2018]]) may show a different K_Ic(T) curve due to dominantly grain-boundary cracking. **An open assumption** that the slow-rate K_Ic(T) applies under MMW; the model output will inherit this assumption.

### Stage-1 parameter values

For **Westerly granite** (Stage-1 default; full Nasseri 2007 + Friedrich-Wong validation chain available):

| Parameter | Value | Source | Notes |
|-----------|-------|--------|-------|
| `K_Ic(T)` | tabulated 1.43 → 0.22 MPa·m^0.5 over RT–850 °C | [[sources/nasseri-2007]] Table 2 | piecewise interpolation; concave with quartz-transition knee |
| `E` | 100 GPa (room T) | [[sources/friedrich-wong-1986]] App. 2 | |
| `ν` | 0.23 | [[sources/friedrich-wong-1986]] App. 2 | |
| `Δα` (mismatch) | 1.9 × 10⁻⁶/°C | [[sources/friedrich-wong-1986]] App. 2 | plag inclusion in qtz/microcline matrix |
| `Δα'` (anisotropy) | 1.4 × 10⁻⁶/°C | [[sources/friedrich-wong-1986]] App. 2 | |
| `a_0` (Weibull, GB) | 0.26 mm | [[sources/nasseri-2007]] Table 4 | mean GB crack length, RT |
| `a_0` (Weibull, IG) | 0.13 mm | [[sources/nasseri-2007]] Table 4 | mean IG crack length, RT |
| `m_W` (Weibull) | 15 default; sweep [8, 25] | [[sources/lyu-2018-mineralogy]] + literature | parametric sensitivity |
| `(∂E/∂P)` | 5–10 | granite handbook | linear EOS simplification |

For **Central Aare granite** (Kant 2017 validation case):
- `K_Ic` = 1.5 MPa·m^0.5 (room T; Kazerani & Zhao 2014). T-dependent values not measured for this specific rock; use the Nasseri 2007 Westerly curve scaled to the room-T value as a first approximation (factor 1.5/1.43 = 1.05).
- `a` = 20 µm (Moeri 2003 thin-section; ~order of magnitude smaller than Westerly — different rock).
- `E` = 30–40 GPa, `ν` = 0.19–0.33, `α` = 7.8–8.2 × 10⁻⁶/°C from Kant 2017 Table 1.
- Mineralogy 46% qtz / 31% feldspar / 20% plag (Rossi 2018) — different mismatch values would need to be measured separately for an exact match; Friedrich-Wong's Westerly Δα is a reasonable surrogate for spallability-order purposes.

## 12. Implementation sequence

A staged build that lets each piece be validated independently:

1. **Stage 1 — Sp criterion with deterministic `a`** on the existing Voronoi grid; reproduce Kant 2017's Central Aare validation in 3D MMW geometry. Sanity check: predicted onset ΔT in the Kant range at flame conditions; well-behaved at MMW conditions.
2. **Stage 2 — Add Weibull `a` on grain-boundary faces.** Validate against the [[validation/hu-2018__lrst-granite-sandstone]] onset time. Tune `a_0, m_W` if needed within measured bracket.
3. **Stage 3 — Add depth-resolved h_spall (§4).** Validate against [[validation/rossi-2018__central-aare-granite-strength-vs-T-rate]] damage profile.
4. **Stage 4 — Per-cell v_n with simultaneous spall + vap (§5).** Validate against [[validation/zhang-2023__oglesby-woskov-validation]] ROP at varying power, expecting the model to agree with Zhang in the high-power vap-dominated regime AND produce non-zero spall contribution at lower power.
5. **Stage 5 — Confining-pressure sweep.** Reproduce [[validation/kant-2017__central-aare-granite-spalling-temperature]] decreasing ΔT trend with p.

After Stage 5 the model is the unified two-mode framework promised in main.tex §1 — but with a parameter-light damage block instead of the eq. 27 chain.
