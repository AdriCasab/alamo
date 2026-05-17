---
title: Kant 2017 Central Aare granite — validation reference packet (sp_kant_onset)
type: analysis
tags: [validation, kant-2017, central-aare-granite, spallability-number, lefm, sp-weibull, reference-data]
created: 2026-05-14
updated: 2026-05-14
sources: [[[sources/kant-2017]], [[validation/kant-2017__central-aare-granite-spalling-temperature]]]
---

# Kant 2017 Central Aare granite — validation reference packet

Self-contained data for the `tests/MMWSpalling/sp_kant_onset/` validation case
(Step 15b of [in-main-tex-you-will-quizzical-treasure.md](in-main-tex-you-will-quizzical-treasure.md),
`sp_weibull` branch). Distilled from [[sources/kant-2017]] and
[[validation/kant-2017__central-aare-granite-spalling-temperature]] — read those for
full derivation and context. Everything an implementer needs to set up and score the
run is here.

## 1. What this validates

The LEFM Spallability-number onset criterion `Sp = K_I / K_Ic ≥ 1` on Central Aare
granite under a flame-jet (convective) surface heat load, at confining pressure
`p = 0`. The confining-pressure sweep (`p ∈ {0, 27, 48} MPa`) is a separate case
(Step 18). This case exercises: the §15a convective Robin BC, the §12b `K_I`
integrator with the `−σ_xx` sign convention, the `K_Ic(T)` table, and the §15b
Weibull flaw distribution.

## 2. Geometry

- Experimental sample: Central Aare granite cylinder, 85 mm diameter × 150 mm length.
- ALAMO approximation: a Cartesian block is fine (Kant's own model is a 1-D
  half-space; exact lateral shape is not critical). Suggested ≈ 85 × 85 × 150 mm,
  or any block whose lateral extent is ≥ ~3× the heated-patch diameter so the
  surrounding cold rock supplies the `ε_xx = ε_yy = 0` confinement on the heated zone.
- Heated patch: circular, diameter ≈ 46 mm (flange aperture, Fig. 5), centred on one
  85 × 85 face. The flame-jet flux distribution within the patch is **not specified**
  — use a uniform `h` over the patch; the revised-approach plan suggests a
  20 / 30 / 46 mm patch-diameter sensitivity check on at least one run.

## 3. Material — Central Aare granite (room-temperature properties, Kant Table 1)

| Property | Symbol | Value (bracket) | Use |
|---|---|---|---|
| Young's modulus | `E` | 30–40 GPa (mean 35) | bracket drives the predicted-ΔT range |
| Poisson's ratio | `ν` | 0.19–0.33 (mean 0.26) | " |
| Mode-I fracture toughness (room T) | `K_Ic⁰` | 1.5 MPa·m^0.5 | Kazerani & Zhao 2014 |
| Largest pre-existing crack length | `a` | 20 µm | Moeri 2003; dominant lever (`a^{3/2}`) |
| Linear thermal expansion | `α` | (7.8–8.2)×10⁻⁶ K⁻¹ (mean 8.0e-6) | |
| Initial temperature | `T_i` | 20 °C (293.15 K) | |

**Not in Kant Table 1 — use standard granite values** (needed for the transient heat solve):
- Thermal conductivity `λ ≈ 2.5 W/(m·K)`
- Thermal diffusivity `κ ≈ 1.1×10⁻⁶ m²/s`
- Density `ρ ≈ 2630 kg/m³`, specific heat `c_p ≈ 800 J/(kg·K)`

**Temperature-dependent `K_Ic(T)`** — Nasseri 2007 Westerly-granite table scaled to
Central Aare's room-T value (×1.05 = 1.5/1.43):

| T (°C) | 20 | 250 | 450 | 650 | 850 |
|---|---|---|---|---|---|
| `K_Ic` (MPa·m^0.5) | 1.50 | 1.42 | 1.03 | 0.45 | 0.23 |

Piecewise-linear interpolation; concave with a knee near the 573 °C quartz transition.

## 4. Boundary conditions

- **Heated patch:** convective Robin BC, `q = h_fl·(T_flame − T_surface)`.
  - `T_flame ≈ 1400 °C` (1673 K) — CH₄/air premixed flame jet.
  - `h_fl` is **not reported numerically**. Bracket `h_fl ~ 10³–10⁴ W/(m²·K)` and
    confirm the Biot regime `Bi = h_fl·a/λ < 0.5` (Kant's simplification assumed
    `Bi ~ 0.08` with `a = 20 µm`). `h_fl` is the free knob for this case.
- **Confining pressure:** `p = 0` for this case (sweep is Step 18). Applied
  uniaxially orthogonal to the heated face when non-zero.
- **Top face outside the patch / sides / bottom (thermal):** Kant does not specify;
  adiabatic or weak convective to `T_ext = 20 °C` is acceptable — the criterion is
  near-surface and early-time, insensitive to far BCs.
- **Mechanical:** minimal bottom support to remove rigid-body modes; other faces
  traction-free. Initial displacement zero. Heating continues until first spall.

## 5. The criterion (for scoring and the closed-form cross-check)

Stress-intensity factor (Kant Eq. 6), with the **§12b sign convention**
(`σ_xx` from ALAMO is tension-positive; the integrand consumes `−σ_xx`,
compression-positive, because Kant's criterion is compression-driven):

```
K_I = 2√(a/π) · ∫₀¹ (−σ_xx(ξa)) · f_K(ξ) dξ
f_K(ξ) = (1.3 − 0.3 ξ^(5/4)) / √(1 − ξ²)
Sp = K_I / K_Ic(T_surface) ≥ 1   at onset
```

**Onset-temperature closed form (Kant Eq. 15) — this is the equation to validate against.**
Setting `K_I = K_Ic` with the temperature assumed constant over the (small) crack length:

```
Θ_onset = (T − T_i) = (1−ν)/(Eα) · [ K_Ic / (2√(a/π) · ∫₀¹ f_K dξ) − p·ν/(1−ν) ]
∫₀¹ f_K dξ = 1.763        (Kant writes the prefactor as √(π/(4a)) = 1/(2√(a/π)))
```

At `p = 0`: `Θ_onset = K_Ic·(1−ν) / (2 · 1.763 · √(a/π) · Eα)`. It is **`h_fl`-independent
and `λ`-independent** — onset is a material threshold; `h_fl` sets only the onset *time*.
Scales as `a^(−1/2)` (halving `a` → onset ΔT × `√2 ≈ 1.41`).

> **Do not use Kant's Eq. 14 `Sp_red` as the onset criterion.** `Sp_red = Sp_rock·Sp_fluid + Sp_conf`
> is a *time-supremum rock-classification number* — it answers "can this rock/fluid/pressure
> configuration spall at all," and its `h_fl`/`λ` dependence enters only through the `Bi`
> substitution in Kant's Eq. 12. It is not parameterized by surface temperature and has no
> `Θ_onset`. An earlier draft of this packet conflated the two; corrected 2026-05-14 after
> ingesting the full Kant 2017 paper. See [[sources/kant-2017]].

## 6. Validation targets

1. **Binding target — predicted onset ΔT at `p = 0` in the range 390–560 °C**
   (mean ≈ 475 °C). This is Kant's *closed-form predicted* range (Table 2); ALAMO
   runs the equivalent onset calculation, so this is what the run must reproduce.
2. **Context anchor — measured spalling temperature 553–694 °C at `p = 0`.**
   Measured > predicted is *expected*: the 1.5 mm pyrometer spot cannot catch the
   first (smallest) spall, so the reported temperature is shortly after first spall.
   Do **not** score against this; it is the experimental upper anchor only.
3. **Closed-form cross-check** — the FEM onset ΔT reproduces Kant Eq. 15
   (`Θ_onset = K_Ic·(1−ν)/(2·1.763·√(a/π)·Eα)` at `p = 0`) within a few %. The FEM,
   integrating the actual `σ_xx` over a sub-grid crack, is structurally in Eq. 15's
   constant-Θ-over-crack regime, so this is a well-posed check.
4. **`a`-sensitivity** — halving `a` raises predicted onset ΔT by `√2 ≈ 1.41×`
   (Eq. 15 scales as `a^(−1/2)`).
5. **`h_fl`-independence** — onset ΔT is unchanged across `h_fl` values (Eq. 15 has no
   `h_fl`). Run two `h_fl` values and confirm the same onset ΔT; `h_fl` shifts only the
   onset *time*. Pick `h_fl` for transient resolvability, not to hit a target.
6. **Weibull patchiness (§15b)** — with the flaw distribution active, first firing is
   grain-scale patchy: ≥ 5 distinct firing locations within the heated footprint at
   onset, not an axisymmetric ring.

## 7. Implementer assumptions / under-specified items

- `h_fl` not given — bracket `10³–10⁴ W/(m²·K)`, confirm `Bi < 0.5`.
- Flame stand-off distance not in this paper (see Kant & Rudolf von Rohr 2016).
- T-dependence of `E, α, λ` not characterised for Central Aare granite; room-T
  values + the scaled `K_Ic(T)` table (§3) are the agreed inputs.
- `a = 20 µm` is from thin-section analysis, not measured on the spalled sample; it
  is the *largest* flaw, so the deterministic prediction is a worst-case lower bound
  on ΔT. The §15b Weibull distribution replaces the single-`a` assumption.
- Kant's model is a single plane-strain edge crack normal to the surface; the real
  rock has a distributed grain-boundary flaw population — that is exactly what the
  `sp_weibull` extension adds.
- Cylinder → Cartesian-box approximation is acceptable (see §2).

## 8. Source

[[sources/kant-2017]] — Kant et al. (2017), *JGR: Solid Earth* 122, 1805–1815,
doi:10.1002/2016JB013800. Section 4 (validation), Table 1 (properties), Table 2
(predicted ΔT), Figs. 5–7. See also [[validation/kant-2017__central-aare-granite-spalling-temperature]].
