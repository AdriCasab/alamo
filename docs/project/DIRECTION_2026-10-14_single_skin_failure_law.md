# Direction: one skin-failure law for spallation and vaporisation

Agreed between the user and the advisor on 2026-10-14. This document records the decision and the
constraints that follow from it. The planning files (`ROADMAP.md`, `ACTIVE_STEP.md`) are updated through
`/plan`; this document is the reference they point to. Background and the arithmetic behind every number
here: `Claude_markdowns/2026-10-14.md` §6–§13; literature: `tests/MMWSpalling/validation/hoeser/
literature_minimum_flux_notes.md`.

## 1. The decision

**One physics model, not two.** A single failure law for the heated skin and a single enthalpy field.
The heat source — a convective flame (Robin surface patch) or a millimetre-wave beam (Beer–Lambert volumetric
deposition) — is the only thing that changes between cases. The regime at any point of the hole (nothing,
flake, melt, vapour) is an **outcome the model computes per column, never an input or a switch**.

**Why the law must read the temperature profile, not the surface flux.** Under a flame the surface flux and
the surface gradient are the same number (Fourier). Under a beam they are not: the heat is deposited over the
absorption depth, the surface loses heat by radiation and gas, the temperature maximum can sit below the
surface, and the net surface flux is zero or negative. A rule written on surface flux has no meaning under a
beam. A rule written on the temperature difference across a material length reads either profile.

## 2. The law

A flake forms when the thermal strain mismatch between the skin and the rock a fixed distance below it,
reduced by a softening function of the skin's own temperature, drives the stress intensity at a flaw past the
toughness:

```
Sp = R0 · s(T_s) · E·β·[T_s − T(z_s + d_ref)] / (1 − ν) · 1.12·sqrt(π·a) / K_Ic(T_s)   ≥ 1  →  flake
```

| symbol | meaning | status |
|---|---|---|
| `T_s` | the face's reconstructed surface temperature (P0) | exists |
| `T(z_s + d_ref)` | the temperature at the **constraint depth** below the face, read from the resolved field along the face normal (vertical column for top faces, the row behind a lateral face) | **new (step A)** |
| `d_ref` | the constraint depth — a fixed material length, §3 | **new (step A)** |
| `R0` | the measured strength ratio `UCS/σ_el(onset)` ≈ 0.53 (P2-2b/P2-2d) | exists, kept |
| `s(T_s)` | softening: 1 below ~900 K, falling to 0 toward the solidus (~1,200–1,450 K) | **new (step B)**; the P1 `melt_aware` skeleton is the hook |
| `E, β, ν` | at-temperature Grimsel tables | exist |
| `a` | the Weibull-drawn flaw, scale 8 µm | exists, unchanged |
| `K_Ic(T)` | the at-temperature toughness table | exists |

If the skin never reaches `Sp = 1`, the enthalpy field carries the column on: melting at the solidus with the
latent heat in `H`, and removal by the vapour front when `H ≥ H_vap` (Step 11; spall has priority over vapour
at a column). Nothing else decides the regime.

**Limits the law must reproduce** (each is a gate in step A or B):

- thin hot skin on cold rock → today's absolute criterion; the banked Kant onset 843.76 K stands;
- slowly soaked face → drive ∝ the conducted flux; in steady state fires only above a **minimum flux**
  `q_min = k·ΔT_c/d_ref` (Kant 2016's measured 0.3–0.5 MW/m²), so the Meier bore wall at 0.11–0.17 stops;
- continuous drilling → the floor runs hotter than first onset (the moving-boundary profile `ℓ = κ/v`):
  ~905 K at Meier's 1.5 m/h, ~870 K at 2.1 m/h — Wilkinson's steady-state rise;
- uniformly heated slab → zero drive (`unit/gradient_blind`);
- interior-peaked profile (deep absorption) → the differential is negative → **no surface spall**;
- skin approaching softening → `s → 0`, the brittle law switches itself off, the enthalpy path takes over.

## 3. The material length — fixed

`d_ref` is a property of the rock, in the same class as the flaw size, the toughness and the compressive
strength. Three rules, to be written as guardrails in every packet that touches it:

1. **Set once, on the calibration set.** Calibrated with the exact transient (erfc) solution against Kant
   2016's spall / no-spall standoffs (Gotthard spalls at 130 mm, not at 140). ⚠ Revised 2026-10-19 (user review): since d_ref ≈ k·ΔT_c/q_min the
   calibration must use Kant's rock's conductivity (Gotthard λ 3.27) ⇒ offline **~8.5 mm**; the property choice
   alone spans ~4.5–9 mm, so the value is set by the in-binary Gotthard ladder (P2-2e b2), not by this note. Kant is the programme's
   calibration set; Meier remains the test set (standing rule in `ROADMAP.md`).
2. **Never adjusted per case.** If Höser's periphery, the Meier twin or a beam column disagree with the
   prediction, that is a finding about the law, not a reason to move the length.
3. **A length in metres, never in cells.** `q_min ∝ 1/d_ref`, so a cell count would make the physics
   mesh-dependent. Converted per mesh (3 cells at 2 mm, 6 at 1 mm) and interpolated; the mesh twins check it.

Two independent routes bracket the same value, which is the reason to believe it is physical: half the
buckle wavelength of the thickest flake, `(C_L/2)·h_f,max = 5 × 1.5 mm = 7.5 mm` with Rauenzahn's aspect ratio
≈ 10; and grain size (1–5 mm) as the lower bound.

## 4. Why this form (and not the earlier ones)

- **A reference depth tied to the heated-layer thickness** (P1, 2026-10-01) is scale-free: the drive fraction is
  the same at every flux, so no minimum flux can exist; it needed a tenfold flaw refit and stranded half the
  Meier floor. P1's rejection of a *fixed* depth rested on a leg at 0.12 MW/m², three times below Kant's
  minimum, which the experiment says must not fire.
- **The column minimum** (P2-2a rev 6) is non-local and switches the model off under volumetric heating.
- **Isotropic averages** of the surrounding rock (ball, shell, smoothed field) are softer than a read behind
  the face: they put the skin's own neighbours into the reference and lose the Kant level for the same wall
  cure (§10 table). The constraint for a face's skin is the rock behind that face.
- **A flux cut-off** is a threshold, not a differential, and cannot be defined under a beam.
- Literature: every analytical model references stress to the far field; the published minimum flux is Dey &
  Kranz's plate-buckling bound; wherever a model has flux dependence a fixed physical length carries it
  (Walsh & Lomov 2014, Seagraves 2019); the thermal-shock stress-reduction factor `ψ(Bi) ≈ 0.31·Bi` is the
  formal cousin of `d_ref·q/k`.

## 5. The order of work

Each item is one packet in the `/plan → /implement → /verify` loop. A and B are validated on flame data
because only there does each end of the law act alone; C and D are the application.

| step | change | calibrated on | validated on | delivers |
|---|---|---|---|---|
| **P2-5** (running) | none — the Höser lab deck at a measured flux | — | Höser 2016/2018 | the flame-end anchor at 7 MW/m²; the periphery row is the first minimum-flux datum |
| **A — P2-2e** | the constraint depth `d_ref` at the four stress sites (lateral along the face normal); key default 0 = today | Kant 2016 spall/no-spall | Kant 2017 onsets (845–847 K predicted); the h = 500 leg marginal (952 K) and SOD 150 a **negative control**; Höser periphery; Meier S0 twin (wall share, mouth, rate, frozen-column census); `gradient_blind`; 1-D volumetric columns at two absorption depths (interior peak ⇒ no spall) | the low-flux end; the first beam prediction |
| **B — softening** | `s(T_s)` replacing the constant relief's physics; `R0` kept | literature band (quartz α–β 846 K → solidus) | Wilkinson's limiting rate at melting; Höser's 4–7 MW/m² onsets (no suppression); Kant level unchanged | the high-flux end; the flame's spall→melt boundary |
| **C — α(T)** | temperature-dependent absorption with energy closure (`MMWSpalling.H:6166–6169`: the deposition sum must telescope) | Zhang/Oglesby loss data | Zhang/Oglesby thermal case, energy ledger | the mechanism that thins the heated layer |
| **D — the regime map** | the deliverable: power density 1–100 MW/m² × absorption depth 0.5 mm–3 cm (frequency) × spot size, spall + vapour on | — | vapour end: ablation rates (data still to be found); spall end: the same law on the flame data | regime per column, rate, specific energy (GJ/m³), hole shape (spalled annulus, vaporised core) |

**Parked:** P2-5b (Meier floor re-levelling), the cold body and wall march, the Meier shape. Meier becomes a
regression that must not break, not a target. The flame work ends when B is accepted.

## 6. Standing constraints

- One physics change per packet; the tree identity holds at the key's default.
- Nothing is tuned to Meier or to a removed volume; the strength side of the law (`R0`, `a`, `K_Ic`) is one
  calibration knob and is not moved by A or B.
- Every prediction is pre-registered with its failure meaning; a wrong prediction is a finding.
- The relief factor is in use today (`spallation.sp.stress_relief = 0.531294` on every current leg). It is not
  removed by A; B replaces its physics with a temperature dependence and keeps the measured level.

## 7. Known unknowns

- A hot subsurface layer under deep absorption may fail by tensile cracking (Guo & Dyskin 2024), a mode
  outside a skin law; the model will say "no spall" there.
- No quantitative spallation data exist under millimetre waves; the map's spall side is a prediction whose
  physics is validated on the flame.
- The softening band is a new parameter, bounded by the literature, with one validation point (Wilkinson).
- The ablation-rate data for the vapour end are not yet on disk.
- Buckling onset is a bifurcation (LLNL grain-scale simulations): sensitive to roughness and heterogeneity;
  the frozen-column census guards the model-side version of that sensitivity.

## Addendum 2026-11-03 — step A superseded; the constraint depth is retired (P2-16)

Agreed between the user and the advisor on 2026-11-03. **§2's reference term `T(z_s + d_ref)`, §3 ("the material
length — fixed") and §5's step A no longer describe the code.** The rest of this document stands.

- **Why A is retired, on form rather than value.** Under any constraint depth the onset climbs steeply as the flux
  approaches the minimum. Re-fitting d only slides the climb along the flux axis (P2-14a c1: 1001 → 879 K over
  0.34–1.59 MW/m² at d = 2 mm). Kant 2016's own onset-by-standoff (Fig. 11), in the dataset that defines the minimum
  flux, shows no climb toward it. Half-space thermoelasticity has no reference depth: the in-plane stress is local in
  temperature (literature study 2026-11-02, `Claude_markdowns/reports/Granite spallation and melt physics.md`).
  P2-16 deletes `spallation.sp.constraint_depth`; a deck naming it aborts.
- **What replaces A's two jobs:**
  - *The onset* is today's criterion, local in temperature and referenced to the initial field (§2 with
    `T(z_s + d_ref)` read as `T_0`). It is flat with flux, as measured.
  - *The minimum flux* is a separate buckling thickness condition (P2-15, `spall.buckle_span`, face list only,
    default off). A flake of normal thickness x = c_f·ℓ fires only if C_L·x ≤ L_span, with
    C_L = π/√(6(1+ν)ε_th) (Rauenzahn & Tester 1989 Eqs. 15–16, after Preston, Dey & Kranz). The onset stays flat
    above the minimum and nothing fires below it: a step, as Kant measured.
  - L_span is a rock property in metres, set once. Whether it follows a measured property (e.g. ∝ 1/UCS) is P2-16's
    test. ⚠ L_span and c_f are one calibration.
- **What A's beam argument becomes (§1).** The case for "a rule written on the temperature profile" stands, but the
  quantity read from the profile is now the heated-layer depth ℓ, not a temperature difference across d. Under an
  interior peak ℓ is undefined today (the reconstruction reads "flat"). Defining it is the first design question of
  the MMW slice (P2-14c).
  - Expected regimes: a shallow absorption depth gives flame-like spalls; a deep one demands a span no granite has,
    so the brittle law switches itself off and the column goes on to melt.
  - A block heated through (~10 GHz) is a finite body whose cooled skin can be in tension, a mode the local criterion
    cannot represent (keep the full elastic stress solve for it).
- **§5's order now reads:** A done and retired → buckling condition (P2-15, P2-16) → merge slices (feet, MMW with
  the interior-peak ℓ) → B softening, C α(T), D the regime map, unchanged in intent.
