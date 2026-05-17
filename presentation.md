---
marp: true
theme: default
paginate: true
title: Thermomechanical MMW Spallation Model in ALAMO
---

# Thermomechanical MMW Spallation Model in ALAMO

A unified continuum model for millimetre-wave rock drilling

Implementation plan summary — 2026-05-13

---

## What we are modelling

**Millimetre-wave drilling** — a high-power microwave beam (~30–300 GHz, kW–MW)
heats rock surfaces fast enough that the rock fails *before* it can melt or
vaporise efficiently. Two competing removal mechanisms:

- **Spallation** — thermal stresses crack and detach mm-thick chips.
- **Vaporisation** — surface boils away once enthalpy exceeds latent heat.

The right model has to handle **both regimes** plus the transition between them.

---

## Six coupled physics ingredients

1. **Enthalpy heat conduction** with solid / liquid / vapour phase tracking.
2. **MMW beam source** — Gaussian profile, Beer–Lambert absorption.
3. **Surface losses** — radiation + convection.
4. **Voronoi mineral microstructure** — per-grain properties, grain boundaries.
5. **Thermoelastic mechanics** with damage-degraded stiffness.
6. **Material removal** — spall detachment OR vaporisation recession.

All on an **adaptive mesh refinement (AMR)** grid using the ALAMO framework.

---

## Why ALAMO

ALAMO is an AMReX-based multi-physics framework with:

- A composable **integrator** pattern (multiple inheritance for coupled physics).
- Built-in nodal **elasticity operator** with multigrid solver.
- Voronoi and parser-based **initial conditions**.
- **AMR** with regrid/average-down already wired up.

We inherit from `HeatConduction` + `Mechanics<…>` and add the MMW-specific
coupling. **No new executable framework** — just one new integrator and one
new top-level binary, `mmwspalling-3d-g++`.

---

## The plan: 18 incremental steps

Each step ends with a **verification or validation test** with a pass criterion.

- **Steps 1–4** — Pure thermal physics (no mechanics).
- **Step 4b** — *Validation milestone* against Zhang/Oglesby MMW heating data.
- **Steps 5–6d** — Heterogeneous microstructure + Hu thermal validation.
- **Steps 7–7d** — Thermoelasticity, AMR, Hu thermoelastic validation.
- **Steps 8–11b** — Damage, cohesive zones, removal, Hu LRST/onset validation.
- **Steps 12–18** — Revised LEFM **Sp + Weibull** branch + Kant/Rossi data.

---

## Validation campaigns

Four independent experimental reference datasets, each anchoring a different
piece of the model:

| Source | Validates | Status |
|---|---|---|
| Zhang & Oglesby (granite MMW heating) | Heat + beam + losses | ✅ done (Step 4b) |
| Hu (prescribed-T flame heating, granite/sandstone) | Thermal, stress, damage chain | ✅ 6d, 7b/7d, 8b, 11b |
| Kant 2017 (flame jet, confining pressure) | LEFM onset, pressure dep. | ⏳ Step 12c, 18 |
| Rossi 2018 (post-flame thin-section crack density) | Spall depth profile | ⏳ Step 16 |

---

## Where we are today

✅ **Completed (Steps 1–11b):** full thermo-mechanical-damage pipeline,
including Voronoi microstructure, AMR for thermal-stress, Drucker–Prager
damage, cohesive zones (unit test), and a working spall + vaporisation
removal loop. Three of four validation anchors satisfied.

🔜 **Next (Steps 12–18):** a runtime toggle adds an alternative
LEFM-based spallation block based on a stress-intensity-factor `Sp = K_I/K_Ic`
criterion + per-face Weibull flaw statistics, then closes Kant and Rossi
validations.

---

## Steps 1–4: Thermal foundation

| Step | Adds | Verification |
|---|---|---|
| 1 | Skeleton integrator + build | Runs to t=1 |
| 2 | Enthalpy form (solid/liquid/vapour) | Stefan analytic, **4.23%** front-position error |
| 3 | Gaussian + Beer–Lambert MMW source | Domain-integrated power, **1.61%** error |
| 4 | Surface radiation + convection | Equilibrium temperature, **0.00%** error |

Forward-Euler explicit; latent heat handled by **closed-form T↔H inversion**
(no Newton iteration). The Voller–Prakash enthalpy formulation makes phase
plateaux exact.

---

## Step 4b: Zhang/Oglesby validation (mid-plan milestone)

**Setup:** 10 × 10 × 2.5 cm granite, 28 GHz beam (ω₀ = 2 cm), piecewise-constant
power schedule 800 W → 3200 W over 1750 s.

**Pass:** ALAMO peak surface temperature **2871 K** vs. Zhang's simulation
**2850 K** at t = 1250 s — **0.7%** error (criterion < 15%).

**What this exercises end-to-end:** enthalpy w/ latent heat, T-dependent
material properties (ρ(T), κ(T), c_p(T) polynomial fits), Beer–Lambert
absorption with reflectivity, radiation + convection losses, **and** the
implicit backward-Euler thermal solver added during tuning.

---

## Steps 5–6c: Heterogeneous microstructure

- **Step 5** — Reproducible seeded Voronoi grains; per-phase
  `{κ, β, E, ρ, c_p, A_D}`; `is_gb` flag at grain boundaries.
- **Step 6** — Damage-degraded conductivity `κ_eff = κ_phase · exp(−α_D D)`;
  optional expression-mode microstructure for deterministic laminate tests.
- **Step 6b** — H/T consistency: every phase has its own enthalpy table; the
  heterogeneous path now advances conserved enthalpy, not temperature.
- **Step 6c** — **Grain topology ≠ phase topology**. Same-mineral neighbouring
  grains can now have grain-boundary effects (cohesive zones, elevated `A_D`).

---

## Step 6d: Hu prescribed-T thermal validation

First Hu reproduction: **flame-jet patch BC** on a granite/sandstone cube,
mechanics off, heat conduction only.

- Hole-1 ΔT range: granite 0 → **35.3 K**, sandstone 0 → **73.2 K**.
- At t = 30 s, L1/L2 line-temperature ratios match Hu within ~3 percentage
  points (sandstone L1 +15.84% vs Hu +12.58%; L2 +17.38% vs +18.90%).

Establishes the **prescribed-T patch BC** that all subsequent Hu validations
use, ahead of bringing mechanics online.

---

## Steps 7–7d: Thermoelasticity + AMR

- **Step 7** — Per-phase `(E, μ, β, T_ref)` mechanical fields, nodal averaging,
  full microstructure stiffness rebuild every step.
- **Step 7b** — Hu Granite-2 / Sandstone-2 thermoelastic to 30 s, single level.
  Peak top von-Mises **granite > sandstone**, as Hu reports.
- **Step 7c** — AMR regrid repair for discrete fields (`phase`, `grain_id`,
  boundary flags, properties). Linear interpolation of integer phase IDs would
  smear them; we **re-evaluate from physical coordinates** instead.
- **Step 7d** — Hu thermoelastic rerun with AMR: ~28× fewer cells than uniform
  refinement, agreement to 7b within tolerance.

---

## Steps 8 + 8b: Damage and breakage indicator

**Step 8** — Drucker–Prager yield `F_DP = √J₂ + α_DP I₁ − k_DP`. Continuous
damage evolution (RK4, irreversible):
$$ \dot D = A_D\,\bigl(\Psi/(1+\Psi_0)\bigr)^n (1-D)^m \cdot \mathbf 1[\Psi \ge \Psi_0] $$
Stiffness degradation `E(D) = E (1 − D)`.

**Step 8b — Hu breakage index** `f_b = σ_v / σ_s`, *damage evolution off*:

- Granite f_b = **0.86** (Hu: 0.91), Sandstone **0.84** (Hu: 0.82).
- Granite/sandstone margin **+16.8%** vs Hu +19.7% — correct ordering.

---

## Step 9: Grain-boundary cohesive zones

Bilinear traction–separation law on **true grain-boundary faces** (not just
phase boundaries):

- Linear loading to peak traction f_t at δ_c, linear softening to zero at δ_max.
- Calibrated so integrated work = input fracture energy `G_c = ½ f_t δ_max`.
- GB tensile strength set to **40%** of intra-grain strength.

**Status:** unit-test verified (traction–separation curve overlays analytic
bilinear law to plotting precision). **Traction feedback into the elastic
operator deferred** — currently a passive diagnostic.

---

## Steps 10–11: Material removal

**Step 10 — Spall detachment.** Level-set field `φ` tracks the surface; columns
where `D ≥ threshold` *and* `σ₁/σ₃ ≥ ξ_s ≈ 8` (Hoek–Brown) detach with
thickness `h_spall = C_h √(α t_spall)`. Detached cells reset to ambient.

**Step 11 — Vaporisation.** Once `H ≥ H_vap` at a surface cell, advance with
`RoP_vap = (1−R) P₀ / (ρ L_tot A_beam)`. Global mode chosen by
`max(RoP_spall, RoP_vap)`.

Both modes share fields; per-column decisions use **MPI reductions** so
z-decomposed columns stay coherent across ranks.

---

## Step 11b: Hu LRST + onset validation

First validation that exercises the **full damage pipeline**.

| | Granite | Sandstone |
|---|---|---|
| ALAMO onset | 37.0 s | 87.0 s |
| Hu target onset | 31.5 s (±20%) | 76.7 s (±20%) |
| ALAMO LRST | 776.6 K | 898.7 K |
| Hu target LRST | 791 K (±10%) | 842 K (±10%) |

✅ All four metrics pass. **Granite spalls before sandstone**, as observed.

**Caveat:** sandstone damage zone is deeper than Hu reports (warning-only in v1)
and per-rock `A_D` values are explicit calibration knobs.

---

## Step 12: The branching point — spallation-model toggle

```
spallation.model = damage_law   (default — Steps 8–11b preserved)
spallation.model = sp_weibull   (revised LEFM approach — Steps 12–18)
```

Toggle is a **full block swap**, not per-component:

| Component | `damage_law` | `sp_weibull` |
|---|---|---|
| Onset criterion | Drucker–Prager | `Sp = K_I / K_Ic(T)` |
| Damage evolution | eq.27 RK4 | none (Weibull statistics) |
| Spall thickness | `C_h √(α t_spall)` | Sp-vs-depth zero crossing |
| GB statistics | CZM (Step 9) | Per-face Weibull flaws |
| Removal granularity | per-column | per-cell |

**EOS work (Mie–Grüneisen, linear E(P)) dropped 2026-05-13.**

---

## Why the revised approach?

The `damage_law` branch works but has **two per-rock free parameters** (`A_D`,
`Ψ₀`) tuned to fit Hu. The revised approach replaces both with quantities that
can be **read from independent literature**:

- `K_Ic(T)` from Nasseri 2007 piecewise-linear table (5 measured points).
- Weibull `(a₀, m)` from Westerly granite microcrack statistics.

The goal is **predictive, not fitted**: same parameters work across Hu, Kant,
and Rossi without per-test recalibration.

---

## Steps 12–13: Sp criterion + AMR refinement

**Step 12 — Sp onset.** Tada weight function
$$ K_I(x,y,t;a) = 2\sqrt{a/\pi}\int_0^1 \sigma_{xx}(x,y,\xi a,t)\, f_K(\xi)\,d\xi $$
with `f_K(ξ) = (1.3 − 0.3 ξ^{5/4})/√(1−ξ²)` and Simpson quadrature with
logarithmic refinement near the crack tip.

**Verification:** Kant 2017 closed-form Bi < 0.5 regime, ALAMO `Sp(t)` matches
analytic `Sp_red` within **5%** over the first 10 s.

**Step 13** — Production AMR refinement criteria (`|∇H|`, `|∇D|` or `|∇Sp|`,
beam footprint) + adaptive Δt per CFL bound.

---

## Steps 15–17: Weibull, depth, per-cell removal

**Step 15** — Per-face flaw lengths sampled at init: GB faces
`a_f ~ Weibull(a₀_gb, m)`, intragranular `a_f ~ Weibull(a₀_ig, m)` with
`a₀_ig ≈ a₀_gb / 5` (Rossi 2018). Re-runs Hu LRST/onset without per-rock
fitting.

**Step 16** — `h_spall(x,y) = max{ z : K_I(x,y,z;a) ≥ K_Ic(T(x,y,z)) }`
evaluated by re-running the Tada integral down a column. Validated against
Rossi 2018 crack density vs depth (peak at 100–200 µm).

**Step 17** — Per-cell normal velocity `v_n(x,y,t)`: spall and vaporisation can
fire **at different (x,y) within the same beam footprint in the same step**.
Re-runs Zhang/Oglesby.

---

## Step 18: Confining-pressure validation (Kant 2017)

Lithostatic pressure applied as a **stress-field BC** on side faces — not
through a `Ψ₀(P_c)` empirical correction. Pressure enters the LEFM criterion
implicitly through `σ_xx`.

**Pass criteria:**
- Onset ΔT at p = 0 inside Kant's measured 553–694 °C bracket.
- Monotone-decreasing onset ΔT across p ∈ {0, 27, 48} MPa.
- Slope vs pressure within ±30% of Kant's measured trend.

---

## Step 14: End-to-end Hu reproduction (both branches)

Runs **twice**, once per branch, on the same Hu setup:

- Under `damage_law`: consolidates 6d → 7b/7d → 8b → 11b.
- Under `sp_weibull`: consolidates 6d → 7b/7d → 12c → 15a — must reproduce
  Hu LRST/onset **without** the per-rock fitting that 11b relies on.

This is the final non-MMW spallation validation. Two passing branches give
us an apples-to-apples comparison of the empirical and LEFM approaches on
the same data.

---

## Open caveats and deferred work

- **Cohesive-zone traction feedback** into the elastic operator (Step 9).
- **AMR moving-surface validation** with spall/vap removal (Steps 10–11).
- **Sandstone damage depth** in Step 11b is larger than Hu — warning only.
- **Implicit thermal solver** is single-level + Neumann BCs; not yet
  compatible with fracture mechanics.
- **Mie–Grüneisen / linear E(P)** dropped 2026-05-13 — no pressure-dependent
  stiffness in either branch.
- Step 11b damage parameters (`A_D`) are **fitted**, not predictive — Step 15
  is the answer.

---

## Summary

- **11 of 18 steps complete**, including three of four validation anchors
  (Zhang/Oglesby, Hu thermal/thermoelastic, Hu LRST/onset).
- The model is a faithful re-implementation of the unified MMW drilling theory
  from main.tex, on top of ALAMO's AMR multi-physics infrastructure.
- The remaining **Sp + Weibull branch (Steps 12–18)** replaces fitted damage
  parameters with literature-anchored LEFM + flaw statistics, validated
  against Kant 2017 and Rossi 2018 alongside the existing Hu data.
- The **two-branch toggle** lets us compare the empirical and LEFM
  formulations head-to-head on identical setups.
