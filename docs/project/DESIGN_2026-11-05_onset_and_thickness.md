# Design note — where the heat flux enters a spall, and what the beam needs (2026-11-05)

Written at the user's request before any packet ("think more before any packet"). It settles the
physics behind the three "skin shortcuts" identified for the beam port (the HTML explainer
`Claude_markdowns/htmls/2026-11-05_three_shortcuts.html`, artifact EsgCShDVP4UVdYRkiZwsYv) after two
findings of the same day: the relief diagnostic (`tests/MMWSpalling/studies/relief_diag/README.md`)
and the discovery that an energy-release flake thickness would remove the minimum-flux mechanism.

The goal it serves: ONE source code that runs a flame jet and a millimetre-wave beam, predicts
spallation versus melting/vaporisation on unseen tests, with the fewest typed parameters and no
shortcut that only holds for a flame.

---

## 1. The chain, and the one question

```
heat source ──► T(x,t) ──► stress field σ(x,t) ──► ONSET: K_I(σ, flaw) ≥ K_Ic(T)   [condition A]
                                 │                      │ fires
                                 │                      ▼
                                 └──► FLAKE THICKNESS x ──► BUCKLING: C_L(T_s)·x ≤ L_span  [condition B]
                                                                          │ passes
                                                                          ▼
                                                                removal, enthalpy leaves
```

The question the user asked: *is the minimum flux not a function of the temperature gradient — no
spall if not enough heat is put in to form a strong gradient?* Yes. The precise statement is **where**
in the chain the gradient acts. There are exactly two physical places, and one place it must not be.

---

## 2. The onset does not see the gradient, and the data say it must not

**Physics.** Thermal stress is not produced by a temperature slope. It is produced by hot rock being
held by cold rock. A hot skin on a cold block cannot expand in-plane, so it is in biaxial compression
σ = E·ε_th(T_s − T_0)/(1 − ν). The stress depends on how hot the skin is relative to the body that
constrains it, not on the slope beneath the surface. For a thin skin on a wide cold block the
constraint is complete and σ is a function of T_s alone.

**Data.** First-spall temperature is a material property band, flat with flux:
- Wilkinson & Tester 1993 (Barre, Westerly): onset 400–600 °C, "relatively insensitive to the applied
  heat flux level", 0.2–2.9 MW/m²; the *steady-state* surface temperature rises with flux
  (500→930 °C), which is a different observable.
- Kant 2017 and Höser 2016 (Grimsel): 430–600 °C across 4.3–7.4 MW/m², "any trend was not recorded".
- Across six granites and two methods: 430–630 °C over 0.2–7.4 MW/m² (literature study 2026-11-02, Q1).

A criterion in which the gradient raised or lowered the onset temperature contradicts this.

**History — three attempts, two refuted by exactly that.**

| step | how the gradient entered the onset | what happened |
|---|---|---|
| P1 (2026-10-01), `driving_stress = differential`: stress jump across the flake, T(z) − T(z + h_f) | the form was exact (drive fraction 1 − e^{−c_f} to 16 digits) and passed the uniform-slab test, but it needed a 10× larger flaw to land the onset, the two cancelled as ONE calibration, and on Meier it **froze 48 % of the columns** under the burner: a column that falls behind soaks heat to depth, its gradient flattens, its drive collapses, it stops. Keying the onset on the gradient makes flatness self-reinforcing (`ARCHIVE_DONE.md` "P1 — the gradient-aware firing criterion"). |
| P2-2e (2026-10-14→24), `constraint_depth`: reference T at a fixed depth d | the onset **climbed as flux fell** for every d (Q = q/q_min 1.25 → ×2.02 … 4 → ×1.00; −122 K from h 2000 to h 500; 1001 K at 0.34 MW/m²). The measured onset is flat. Retired in P2-16. |
| P2-15 → P2-17: the gradient in the **thickness**, not the onset | onset stays local; the layer depth sets the flake; the span blocks thick flakes. Reproduces the flat onset and supplies Kant's minimum flux as a separate condition. |

**The "gradient-blind defect".** P1's uniform-slab test (a block preheated uniformly to the flame
temperature, zero flux) showed the absolute criterion firing 80 times where nothing should fire. That
is real, but it is not an onset defect: a uniformly hot block is **unconstrained** — nothing cold holds
it, so its stress is zero. The 1-D law assumes complete constraint, which that block lacks. The defect
is in the stress source's constraint assumption, i.e. the finite-body problem of section 3, not in the
criterion. Putting a gradient into the onset formula fixed the symptom and broke Meier.

**Conclusion.** The onset criterion stays local: K_I from the local stress at the surface flaw against
K_Ic at the local temperature. No slope term. This is also the literature's form (Kant 2017's LEFM
criterion; Wilkinson's "weakest flaw" onset).

---

## 3. First place the gradient acts: confinement, i.e. the stress field

The stress is σ = E/(1 − ν)·(ε_th − ε_body): the thermal strain of the hot material minus the strain
the surrounding body permits. The gradient's role is through the **ratio of heated depth to the body**:

| heating | heated depth vs body | constraint | σ_s / σ_full |
|---|---|---|---|
| flame skin, 1 mm on a 100 mm block | ≪ | complete | 1 (measured: spot centre = uniform slab within ±15 %, relief diagnostic) |
| 28 GHz beam, 1/e depth 3–7 cm on a 25 mm slab | ≳ | partial, slab bends, rim in tension | < 1, geometry-dependent |
| uniformly preheated block | = | none | 0 |

**Alternatives for the stress source.**

| form | parameters | flame | beam | status |
|---|---|---|---|---|
| S1 · 1-D law × R (production): σ = R·E(T_s)·[ε(T_s) − ε(T_0)]/(1 − ν) | R = 0.531 CAL-ONCE | right form; R unexplained | wrong: assumes full constraint, no rim tension | keep for flame (user, 2026-11-05) |
| S2 · the elastic operator (`spall.stress_source = mechanics`) | none for the stress | **cannot resolve a 1 mm skin**: reads 0.26 / 0.49 / 0.65 of exact at 2 / 1 / 0.5 mm, plateau at 0.65 (relief diagnostic, uniform ladder) | right: confinement and rim tension computed, provided the heated region spans many cells | use for beam legs only |
| S3 · 1-D law with an explicit confinement factor f(δ/W) from the operator's own ladder | one derived function | reduces to S1 | an interpolation of S2 | not needed if S2 runs the beam legs |

**What R is.** R = UCS/σ_full(850 K) was set so the full-constraint stress at the target onset equals
the compressive strength. The relief diagnostic showed it is **not geometry**: the spot centre is as
constrained as an infinite slab. It is a material-level factor — candidates: stiffness loss from
thermal microcracking (the E(T) table already halves by 600 °C; is it applied at T_s?), grain-scale
expansion mismatch (the Voronoi study the user kept), the flaw-size statistics. It stays declared and
is never re-fitted on a beam case. ⚠ With R on the flame path and the operator on the beam path, the
two legs use different stress levels at the same temperature by a factor ~1/R; a beam prediction's
onset is therefore stated with and without R, and the gap is a declared uncertainty until R is
explained.

---

## 4. Second place the gradient acts: the flake thickness, hence the minimum flux

**Physics (Preston 1934; Dey & Kranz 1985/87; Rauenzahn & Tester 1989; Dmitriyev).** A flaw under
thermal compression extends parallel to the surface, freeing a plate; the plate buckles if it is thin
relative to its span. The plate thickness follows the compressed layer, and the compressed layer is
the heated layer, δ = k·ΔT_s/q. Low flux → deep layer → thick plate → no buckling → no spall. That is
the only derived minimum-flux argument in the literature (study Q4), and it is the structure the
model has since P2-15: x = c_f·ℓ, bk = C_L(T_s)·x/L_span, blocked iff bk > 1, with
C_L = π/√(6(1 + ν)ε_th(T_s)) (R&T Eq. 16) and L_span = 12.0 mm a rock property derived on Kant's
Gotthard bracket (P2-17).

**Measured thicknesses** (study Q6): Canadian Mines Bureau field average 0.4 mm; MIT "one to two
grains thick", aspect ratio 8–15; LLNL shadow survey 4:1–10:1; Höser's Grimsel fines peak at 200–300 µm
with 3–4 mm flakes present; Rauenzahn's Weibull median 2.67·ℓ/m. Thinner at higher flux is the
reported trend. The model's 0.4·ℓ gives 360–530 µm on Höser's leg and 700–1300 µm on the calibration
deck — the right order, and `c_f` was chosen for that (model card row 20).

**Alternatives for the thickness selector.**

| form | thickness | parameters | flux dependence | flame | beam | minimum flux |
|---|---|---|---|---|---|---|
| T1 · thermal-skin fraction (production) | c_f·ℓ, ℓ from the two-cell fit | c_f = 0.4 CHOICE | thicker at low flux ✓ | calibrated | **undefined** (ℓ flat/infinite) | ✓ via B |
| T2 · stress-fraction depth | depth where σ(z)/σ(0) = e^{−c_f} (= 0.67) | same c_f, re-read as a stress fraction | identical to T1 on an exponential skin | identical by construction | defined; deep under volumetric heating ⇒ B blocks ⇒ melt before spall (the 28 GHz record) | ✓ via B |
| T3 · energy-release thin limit (Hutchinson–Suo) | smallest t with ∫₀ᵗ(1 − ν)σ²/E dz ≥ K_Ic²(1 − ν²)/E; uniform: K_Ic²(1 + ν)/σ² | none | **none at onset** (σ_s fixed by A) | 0.07–0.10 mm; 3–5× thinner than today | defined | **✗ lost**: window identical at every flux |
| T4 · buckling-selected crossing | depth where σ(z) = buckling stress of a plate of that thickness over L_span | none beyond L_span | thicker at low flux | always has a solution | defined | **✗ lost**: a crossing always exists |
| T5 · Weibull/Rauenzahn | 2.67·ℓ/m median | m (= 13 already) | as T1 | 0.2·ℓ: rate halves | as T1 | as T1 |

T3 and T4 are parameter-free and both delete the floor, for the same reason: they select the
thickness from the stress level, which condition A fixes, so the layer depth — the only carrier of
flux — drops out. **A parameter-free selector that ignores the heated depth cannot reproduce Kant's
minimum flux.** The fraction in T1/T2 is not a shortcut; it is the one declared constant of the
delamination-depth physics (where the parallel crack runs relative to the compressed layer), and it
is a measurable quantity: flake thickness over heated-layer depth, which Grimsel data could pin.

**T3 is still useful** as a lower bound — a flake thinner than t_min cannot free itself. It costs
nothing to compute per line and log; it becomes a third condition only if data ask for it.

**Recommendation.** T2: the same rule for flame and beam, identical to today on flame legs (no
recalibration), defined under the beam, one declared parameter. Replace ℓ by the depth at which the
stress along the line has decayed to e^{−c_f} of its surface value, measured on the same stress the
criterion reads (R·law on flame legs, the operator on beam legs), with the top half-cell's sub-cell
profile as the scan already uses it.

---

## 5. The surface temperature (the second shortcut) — unchanged by the above

The criterion fires on T_s; cells store averages; at 1 mm cells on a 1.4 mm skin the cell reads 72 %
of the surface rise. Some face value is mandatory (dropping it reintroduces a ~100 K mesh-dependent
onset). The two-cell exponential fit assumes heating from above and reads "flat" under a beam, where
the face is cooler than the cell below it. The face heat balance the face list already uses as its
fallback (`FaceRemoval.H` ~526, the Newton root of flux in minus losses minus half-cell conduction)
is sign-agnostic and is not a fit. Recommendation: face balance primary; the exponential kept as the
sub-cell shape only when the net face flux is into the rock; one refined-mesh twin as the arbiter.
Expected effect on flame legs: a few kelvin on flat columns only; a declared re-pin of the face-list
gate.

---

## 6. What the beam needs from this, and the regime prediction

With A local, the stress from the operator, T2 for the thickness and B unchanged, the model makes
one sharp, falsifiable claim under a beam: the stressed layer thickens over time toward the
penetration depth; B blocks thick layers; so spall is possible only **early and at high intensity**,
and otherwise heating runs to melt. Exposure time becomes an axis of the regime map. Tests: no spall
on Oglesby's granite slab at any power (28 GHz, melt at 50–150 s); Liu's 45 GHz no-melt/melt pair;
the laser record's spall at short exposure and 2–3 MW/m²; Oglesby's basalt spall as the one
millimetre-wave spall datum (needs basalt properties, a user decision).

---

## 7. Measurements that would close the open constants

| constant | what pins it |
|---|---|
| c_f (delamination fraction) | flake thickness and the heated depth on the same Grimsel test (shadow method on sieved flakes + IR surface T and flux) |
| R (material relief) | hot stiffness E(T) above 500 °C on Grimsel; or a Voronoi mismatch study |
| L_span | Grimsel's own minimum flux (the production form predicts **0.432 MW/m²**, P2-17 z2) |
| onset band | first-spall T under a beam (the user's 10 GHz rig, later) |

---

## 8. Decision record

- 2026-11-05, user: keep R = 0.531 declared for now; do not re-fit on a beam case.
- 2026-11-05, this note: the onset stays local (section 2); the thickness keeps its one fraction and
  moves to the stress profile (T2); the energy limit is logged, not applied; the operator is the beam
  stress source only. No packet written yet — the user asked to think first.
