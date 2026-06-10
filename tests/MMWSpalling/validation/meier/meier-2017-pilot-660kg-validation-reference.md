---
title: Meier 2017 Grimsel pilot demo (660 kg block) — validation reference packet
type: analysis
tags: [validation, meier-2017, grimsel-granite, flame-jet, pilot-scale, rop, tse, psd, reference-data]
created: 2026-06-02
updated: 2026-06-02
sources: [[[sources/meier-2017]], [[validation/meier-2017__pilot-scale-grimsel-660kg]], [[validation/meier-2017__grimsel-granite-properties]]]
---

# Meier 2017 Grimsel pilot demo (660 kg block) — validation reference packet

Self-contained data for reproducing the Ch. 8 pilot-scale flame-jet demonstration on a 660 kg Grimsel granite block. Distilled from [[sources/meier-2017]] (Chapter 8) and the two related validation cases — read those for narrative context and conflict tracking. Everything an implementer needs to set up and score the run is here.

## 1. What this validates

End-to-end thermal-spallation drilling under a methane/air diffusion flame on hard crystalline granite at near-ambient horizontal stress. The single tightest "what is the energy efficiency of granite spallation, and what does it produce?" datum at the meter scale in the wiki. Exercises: the convective surface heat-flux BC, the rock thermomechanical kernel (whichever you use), the failure/removal criterion, and — critically — the **moving boundary** (the hole advances by 50 cm during the run, ~83 % of the original 60 cm bit reach).

This is not an onset-only test (cf. [[analysis/kant-2017-central-aare-validation-reference]]). It runs to **steady-state advancing drilling** for 23 min and tests whether the model reproduces sustained ROP without decay.

## 2. Geometry

### Rock domain

- **Sample:** monolithic Grimsel granite block, **0.7 m × 0.7 m × 0.5 m**, mass 660 kg.
- **ALAMO / FEM approximation:** Cartesian box, dimensions as above. Lateral extent (0.7 m) is ≥ ~10× the heated-patch diameter so far-field temperature is unaffected; a smaller domain (e.g. 0.2 × 0.2 × 0.5 m) is acceptable provided side BCs are adiabatic or weakly convective.
- **Heated patch:** circular footprint at the bit tip, **diameter ≈ 50 mm** (full bit / wellhead clearance — the hole ends up cylindrical at ~85 mm Ø, but the *direct* impingement is over the nozzle's near field, which at SOD = 7 D ≈ 50 mm has spread to roughly the nozzle-bore + entrainment cone, ≈ 30–50 mm. Use 30 mm if you want a tighter Gaussian or 50 mm uniform.)
- **Hole geometry as it evolves:** cylindrical, ~85 mm Ø, slightly underreamed (the bit advances under gravity through the hole it carves; clearance maintained by drawworks rate, not bit geometry).

### Burner / heat source

- **Nozzle:** converging–diverging Laval (Zimmermann 2016 design), **outlet Ø = 7.1 mm** ("D_nozzle"). An alternative throat diameter of 6 mm is referenced in the cuttings discussion; use 7.1 mm unless explicitly varying.
- **Stand-off distance: SOD = 7 × D_nozzle ≈ 50 mm** (set by burner feet — Meier picked this from Zuckermann's optimal-impingement recommendations).
- **Combustion chamber:** internal Ø 56 mm × 545 mm long (1.34 L), 100 bar rated. Axial swirler 62° trailing-edge angle (radial/axial velocity ratio 1.88). Fuel injection: 12 holes Ø 1 mm at 30° from chamber axis. Pure turbulent diffusion flame (no premixing spacers in this test). **Combustion chamber is internal to the burner; the model needs only the jet exiting the nozzle, not the chamber.**

### Confinement

- **Lateral horizontal stress:** SHmin = SHmax ≈ **1 MPa** (symmetric pads, anti-splitting purpose only).
- **Vertical stress:** atmospheric + block weight ≈ 0.01 MPa. **Effectively unconfined** by downhole standards.
- This is the major caveat for transferring TSE to downhole conditions — Ch. 7 of Meier predicts TSE drops with overburden (depth amplifies tension).

## 3. Material — Grimsel granite (T-independent, Meier Table 7.1)

| Property | Symbol | Value | Use |
|---|---|---|---|
| Thermal conductivity | `k` | 1.5 W/(m·K) | heat solve |
| Density | `ρ` | 2750 kg/m³ | heat + momentum solves |
| Specific heat | `Cp` | 790 J/(kg·K) | heat solve |
| Young's modulus | `E` | 30 GPa | stress solve |
| Poisson's ratio | `ν` | 0.3 | stress solve |
| Linear thermal expansion | `α` | 8 × 10⁻⁶ K⁻¹ | thermoelastic coupling |
| Initial temperature | `T_i` | ~20 °C (293.15 K) | stress-free reference |

**Lithology (vol %):** quartz 28 / plagioclase 29 / K-feldspar 24 / sheet-silicates 18. ~2 mm grain size. (Per [[validation/meier-2017__grimsel-granite-properties]].)

### Not in Meier Table 7.1 — supplement as needed

**Tensile strength** (needed for any tensile-failure criterion):
- Meier states generically "≪ 10 MPa." Use **σ_t ≈ 5–10 MPa**.

**Mode-I fracture toughness K_Ic(T)** (needed for any LEFM criterion):
- Meier does not report. Use [[validation/nasseri-2007__westerly-kic-vs-T]] Westerly granite values scaled if necessary:

  | T (°C) | 20 | 250 | 450 | 650 | 850 |
  |---|---|---|---|---|---|
  | `K_Ic` (MPa·m^0.5) | 1.43 | 1.35 | 0.98 | 0.43 | 0.22 |

**Largest pre-existing crack length** (needed for Kant-style Sp):
- Meier does not report. Use **a ≈ 20 µm** (Kant's Central Aare granite value, [[sources/kant-2017]]) or the grain-boundary length (~2 mm).

**Weibull heterogeneity parameter** (needed for Rauenzahn-Tester analytical):
- Meier does not report. Use **m = 20–25** (Augustine 2009 / Rauenzahn 1986 calibration band; m = 25 is the standard, see [[validation/lyu-2018-se__heterogeneity-m-sensitivity]]).

**Per-mineral properties** (needed for Voronoi/grain-scale models):
- Meier does not provide. Substitute [[validation/walsh-2013__three-mineral-property-table]] (qtz + plag + K-spar) or [[validation/vogler-2020__lac-du-bonnet-microstructure]] (4-mineral with biotite — better match for sheet-silicate-rich Grimsel).

**Granite–stainless steel friction** (only if modeling Meier's specific confined-core setups; irrelevant for this pilot):
- µ_S = 0.4 ± 0.1 (static), µ_D = 0.32 ± 0.07 (dynamic). Measured.

**Derived constants (Meier Eq. 7.4–7.5):**
- Dilatational wave speed `vₑ = 3832 m/s`
- Coupling parameter `δ = 0.0012` ≪ 1 → **uncouple momentum from energy** (Meier does; recommended).

## 4. Boundary conditions — heat source

### Power input

- **Fuel:** **methane (CH₄)** at 2.47 kg/h ≡ 6.86 × 10⁻⁴ kg/s.
- **Oxidizer:** **air** at 52 kg/h, **air/fuel ratio λ = 1.2** (fuel-lean).
- **Combustion power: 38 kW** (HHV basis; Meier Table 8.2). Lower heating value of CH₄ = 50.0 MJ/kg → 34.3 kW LHV; HHV = 55.5 MJ/kg → 38.1 kW HHV.

### Heat flux on the rock

This is the most under-specified part of the test. Two options:

**Option A — Robin convective BC (Meier's own Ch. 7 choice):**
```
q = h · (T_flame − T_surface)
```
- **T_flame:** adiabatic CH₄/air at λ = 1.2 ≈ **1900 K**. Effective jet T at impingement *after entrainment* over SOD = 7 D will be substantially lower; Meier uses **T_ref = 1000 K** in Ch. 7 as a "conservative" choice for hydrothermal flames, which here is too low. Bracket **T_ref ∈ [1200, 1700] K**.
- **h:** Meier uses **h = 10 kW/(m²·K)** in Ch. 7 (taken from Potter Drilling, conservative). For a methane diffusion flame at SOD = 7 D with no swirl reaching the surface, plausible range is **h ∈ [5, 30] kW/(m²·K)**.

**Option B — Energy-balance constraint (preferred for ROP/TSE validation):**
- Set q such that **total integrated input over the patch and duration equals (38 kW × 1383 s) = 52.6 MJ**. This is the energy-input budget that the experimental TSE = 15.4 J/mm³ is derived against.
- **Fraction reaching the rock:** Meier's TSE = 15.4 J/mm³ is computed using *combustion enthalpy*, not measured heat flux. Significant energy leaves with the exhaust gases (~50–80 % depending on combustion efficiency and SOD). So either:
  - feed all 38 kW into the model and accept that the model TSE will inherit the same total-enthalpy convention (apples-to-apples with Meier's reported number); or
  - feed an efficiency-corrected fraction (e.g. 0.3 × 38 = 11 kW into the rock) and compare against **enthalpy-to-rock TSE ≈ 4–6 J/mm³**.

**Recommended for first run:** Option A with `h = 10 kW/(m²·K)`, `T_ref = 1400 K`, uniform over patch Ø 50 mm. Then sensitivity-sweep both knobs.

### Patch shape

- Uniform circular disk of Ø 50 mm centred on the active drilling face. (Meier's Ch. 7 uses a uniform disk of Ø 40 mm = r < 20 mm for a different geometry — adapt to this pilot's nozzle.)
- A Gaussian profile with σ ≈ 10 mm is a defensible alternative; **run a sensitivity check.**

### Moving boundary

- **The bit advances at the measured ROP = 1.5 m/h = 4.17 × 10⁻⁴ m/s.** A model that does not include material removal will not run for 23 min — the steady-state hole position is the validation target, not an assumption.
- A *non-moving-boundary* model can still partially validate against this test by (i) running for a fraction of a second to predict spall onset T and time, (ii) integrating spall-volume rate to predict ROP from first-spall dynamics, and (iii) comparing TSE = power / (volume rate × patch area).

## 5. Boundary conditions — mechanical / thermal far field

- **Side faces of rock domain (x, y boundaries):** weak convective to T_ext = 20 °C (h_ext ≈ 10 W/(m²·K) for natural convection in lab air) or adiabatic — heat diffusion length over 23 min at κ ≈ 7 × 10⁻⁷ m²/s is ~5 cm, so adiabatic is acceptable if domain is ≥ 0.2 m wide.
- **Bottom face:** roller / no-slip (rests on lab floor) and adiabatic.
- **Lateral mechanical confinement:** distributed normal stress σ_lat = **1 MPa** on x- and y-facing side boundaries (this matches Meier's "symmetric pads" lateral load).
- **Top face outside the heated patch:** traction-free; convective to lab air (T_ext = 20 °C, weak h_ext) for thermal.
- **Initial state:** uniform T_i = 20 °C, zero displacement, zero stress.

## 6. Operating timeline (Meier Fig. 8.7)

| Event | Time (s from ignition) | Notes |
|---|---|---|
| Ignite (hot-surface coil, ~100 W ramp, 0.2 A/15 s steps) | 0 | mass flows ramping |
| Steady combustion plateau | ~150 | T_chamber ≈ 1163 °C; ṁ_CH₄ = 2.47 kg/h, ṁ_air = 52 kg/h |
| Drilling phase 1 begins | ~470 | bit starts on rock surface |
| **Phase 1 drilling: 647 s** | 470 → 1117 | constant ROP = 1.5 m/h, advances ~16 cm |
| Break (install exhaust vacuum) | ~1117 → 1450 | ~5 min idle, combustion stays on or cycles |
| Drilling phase 2 begins | ~1450 | |
| **Phase 2 drilling: 736 s** | 1450 → 2186 | constant ROP = 1.5 m/h, advances ~18 cm |
| Bottom thermocouple triggers stop | ~2186 | full penetration through 50 cm block |

**Total drilling time: 1383 s (23.05 min). Total bit advance: 50 cm (through the block).**

For a model: simulate **continuous drilling for 1383 s** with constant input; the break is operational, not physical.

## 7. Validation targets

Ranked by tightness of the experimental measurement.

### Binding targets

1. **ROP = 1.5 m/h, constant.** Model ROP must (i) average within ±20 % of 1.5 m/h over the full 1383 s, and (ii) **not decay** during the run (this directly refutes Silva et al. 2018's claim of rate decline from rock heating). A monotonically decaying ROP is *wrong*, not just imprecise.

2. **TSE = 15.4 J/mm³.** Model TSE (= input thermal energy / excavated volume) must agree within ±30 %. Note the convention: Meier uses **HHV of methane × total fuel mass / measured cavity volume**. A model fed the same total enthalpy and producing 3.4 L of excavated volume passes; a model fed an efficiency-corrected smaller energy must report a correspondingly smaller TSE.

3. **Total excavated volume V = 3.42 L (3420 cm³)** over 1383 s. Equivalent to a constraint on the volumetric removal rate `dV/dt ≈ 2.47 cm³/s` averaged. A model that predicts the right ROP but a different hole shape will miss this; tighter than ROP alone.

### Geometric targets

4. **Hole shape: cylindrical, Ø ≈ 85 mm, slightly underreamed.** A model that predicts a half-spheroid cavity (= heat spreading without sweeping action) is wrong. The cylinder + underream comes from the bit advancing under gravity while the spall ejecta sweep the side walls. **Tolerance: aspect ratio Ø_hole / Ø_nozzle = 85 / 7.1 ≈ 12.** Models should produce 8–15×.

5. **No axial fractures through the block.** The 1 MPa lateral confinement was sufficient to prevent splitting at this geometry. Any model predicting macroscopic axial cracks at SHmin = 1 MPa, 660 kg block + 38 kW + 23 min is over-stressing the failure criterion.

### Particle / fracture targets

6. **Particle size distribution:**
   - **Average < 100 µm** (collected fraction sieved < 500 µm then laser-diffracted).
   - **90 % cumulative < 300 µm.**
   - **Particles are smaller than the ~2 mm grain size** → diagnostic of *intra-crystalline* fracture, not just inter-crystalline along grain boundaries → confirms thermal-stress (not pore-pressure) spallation mechanism.
   - Grain-resolved models (Walsh GEODYN, Vogler MOOSE, Liu damage-FEM) should reproduce the cumulative distribution shape; bulk-continuum models cannot match this directly but can predict avg spall thickness from λ = α/(mv) ≈ 50–100 µm at ROP = 1.5 m/h (Walsh-2014 length scale).

### Combustion / system targets (model-side sanity)

7. **Combustion chamber T plateau ~1163 °C** (steady for both phases) — a sanity check that the model's flame T input is consistent with Meier's instrumentation. Not a rock-side target.

8. **Flame detachment at ṁ_CH₄ ≈ 1 kg/h** during ramp (flame lifts from injector) — operational, not used for model validation.

## 8. Implementer assumptions / under-specified items

- **`h` (impingement heat transfer coefficient) not measured.** This is the single largest uncertainty. Bracket [5, 30] kW/(m²·K) on first sweep.
- **`T_ref` (effective jet temperature at the rock surface) not measured.** Adiabatic CH₄/air ~1900 K; entrainment over 7 D drops it to ~1200–1700 K. Bracket and sweep.
- **Patch shape and size** not directly measured — Meier reports nozzle Ø 7.1 mm and SOD = 7 D; the actual heated footprint at the surface is inferred. Use 30–50 mm Ø; sweep.
- **Energy-into-rock fraction.** Meier's TSE uses combustion enthalpy, not heat-flux measurements — for code-to-code comparison this is fine; for predicting TSE *from* first-principles heat transfer, you need to model the combustion efficiency loss.
- **Lateral confinement profile.** Meier reports SHmin = SHmax ≈ 1 MPa but does not give the spatial distribution of the side-pad load. Uniform on x- and y-faces is the only defensible default.
- **Drawworks mechanics** (the cable that lowers the bit at constant rate while the hole opens beneath it). Meier reports the bit is held in continuous contact with the rock; do not impose a weight-on-bit — it is gravity-driven self-feeding.
- **Initial bit position.** Wellhead glued to top face; bit tip starts at SOD = 7 D = 50 mm above the surface. As the hole opens, the bit follows by gravity. Most models will skip this and start with the bit at the rock surface.
- **Cuttings transport.** Cyclone separates particles ≥ 10 µm. The < 10 µm fraction is lost; PSD targets above are for the **collected fraction (60 % of total mass).** A model predicting the full PSD should adjust accordingly.

## 9. Model-class recommendations

Map the choice of kernel onto the validation targets:

| Kernel | Targets it can hit | Targets it cannot hit directly |
|---|---|---|
| **Meier-style uncoupled thermoelastic (Ch. 7)** | onset stress field; subsurface T; tension-amplification curve | ROP (#1), TSE (#2), V (#3), shape (#4), PSD (#6) |
| **Hu 2018 von Mises FEM** | onset T (#7), first-spall location | sustained ROP, PSD |
| **Liu 2024 continuous damage** | damage volume, modified MSE, mineralogy/heterogeneity sweeps | exact ROP (depends on removal rule), PSD |
| **Lyu / Rauenzahn-Tester analytical** | ROP (#1), TSE (#2) order of magnitude | shape (#4), PSD (#6), no-decay (#1b) |
| **Walsh GEODYN grain-scale** | PSD (#6), pore-fluid amplification, wave effects | runs too short to simulate 1383 s |
| **Vogler MOOSE regime indicator** | onset locations, no-splitting check (#5) | ROP, TSE, sustained drilling |
| **Kant 2017 Sp ≥ 1 analytical** | onset T, depth-amplification | ROP, TSE, PSD |

A complete validation against this test requires a kernel that includes (i) thermomechanical stress, (ii) a failure criterion, **and** (iii) material removal / moving boundary. The full set of targets (#1–#6) is not hittable by any single kernel in the wiki today — that is the gap this validation case maps out.

## 10. Source

[[sources/meier-2017]] — Meier, T. (2017), *Assessment of a contactless drilling tool and its development to access deep underground resources*, Doctoral thesis, ETH Zurich, Diss. No. 24021, https://doi.org/10.3929/ethz-b-000182030.

Primary chapters: **§ 8** (the demonstration itself), § 8.1.1–8.1.5 (process lines, burner, wellhead, control, safety), § 8.2.1–8.2.3 (data reduction, ignition + T evolution, drilling demonstration), Figs. 8.1–8.10, Table 8.2. Property set from Table 7.1 (§ 7.1).

Companion wiki pages:
- [[validation/meier-2017__pilot-scale-grimsel-660kg]] — half-page validation case (this packet's source).
- [[validation/meier-2017__grimsel-granite-properties]] — property table with cross-references.
- [[sources/meier-2017]] — full source page with conflicts and open questions.
