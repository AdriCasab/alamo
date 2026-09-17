# Meier 2017 Grimsel pilot — validation case: status summary

Updated **2026-09-17**. Sections 1–5 below the line are the pre-D2 history
(June 2026) and are kept for the record; several of their numbers and
statements are stale and are flagged in the status block. The physics work
program and development workflow live in [TODO.md](TODO.md). Reference data:
`meier-2017-pilot-660kg-validation-reference.md` (read its **Errata** block
first). Page-cited thesis notes: `Claude_markdowns/meier_thesis_notes.md`.
Reader-facing write-up: `Claude_markdowns/meier_validation_report.html` / `.pdf`.

## Current status (D2b, review-accepted 2026-09-16; corrections 2026-09-17)

**Scored configuration:** `sp_meier_pilot/input_feet` + `sp_meier_pilot/test_feet`
(check-only scorer over `tests/MMWSpalling/studies/d2b_feet_rop/output/B_JM_A2`).
The uniform-beam configuration (`input_drilling` + `test`) is superseded and
kept unchanged for history.

**Model as scored (all inputs from the experiment, the burner drawings or a
published correlation; nothing tuned):**
- 3-D quarter domain 0.12 × 0.12 × 0.30 m, 2 mm cells (1 mm check), pinned
  Robin closure, Weibull flaws (a₀ 4 / 20 µm, m = 20, V0 = V_cell → firing
  T ≈ 821 K).
- Jet: Martin (1977) local h(r, s) at the measured flow (1449 W/m² K under the
  jet, 690 at the foot ring; `jet.py`/`walljet.py`), stagnation decay beyond
  5 D, energy-conserving outward march of T_gas (`jet_closure = enthalpy`),
  entrained gas = the jet's own exhaust (`jet_T_ent_mode = exhaust`),
  T_nozzle = 1900 K (adiabatic).
- Burner on feet: `nozzle_descent = feet`, three pads on r = 28–40 mm, nozzle
  50 mm above; no flux above the nozzle plane.

**Result (B_JM_A2, window 171–342 s):**

| criterion | measured | model | verdict |
|---|---|---|---|
| ROP | 1.3–1.6 m/h (band 1.04–1.92) | **1.36 m/h** (1 mm: 1.27–1.33) | pass |
| ΔT_fire | 500–560 K | 528 K | pass |
| hole Ø below 0.1 m | 85–93 mm | 92 → 85 mm | pass |
| min Ø | ≥ 80 mm | 84.9 (0.5 m projection); 78.2 over the drilled depth | see note |
| depth-mean Ø 0–0.5 m | 85–93 mm | 104 mm | **fail** (mouth funnel) |
| volume rate, whole excavation | 2.0–3.0 cm³/s | 5.68 cm³/s (3.27 inside Ø 93) | **fail** (funnel) |
| steady state | yes | no (centre pit runs to the box bottom) | **open** |
| mesh 2 → 1 mm | ≤ 5 % | ring −4.8 % (trimmed box, not like-for-like under exhaust mode) | soft |

Sensitivities: T_nozzle 1436 K → 0.67 m/h; ambient (293 K) entrainment → stall
after 12 mm; h × 3.5 / × 7 (the old literature values) → 2.43 / 2.82 m/h.

**Reading.** The rate of penetration is reproduced without tuning; the hole
shape is right below the mouth. Two model errors remain: (1) a 200 mm funnel
at the mouth (thesis Fig. 8.8 shows none) from a free-surface wall jet that
never entrains ambient air and exhaust recirculation switched on from t = 0;
(2) no steady state, because the h shape is clamped at 12 D stand-off so the
centre pit never equilibrates. The next packet (D2c, see
`Claude_markdowns/2026-09-17.md`) addresses both and re-scores.

**Corrections from the thesis reread (2026-09-17) that change how this case
is read:**
- The "1163 °C chamber temperature" is the igniter thermocouple upstream of
  the lifted flame (≈ 550 °C while drilling). **The nozzle gas temperature is
  unmeasured**; 1900 K adiabatic is an upper bound; cooling water 160 L/h with
  unreported ΔT. 1436 K is a low case, not a measurement.
- Volume 3.42 L was measured by filling the hole with water → score the whole
  excavation.
- Hole Ø 85 mm is Table 8.2 (p. 222). Nozzle 7.5 mm (drawings) vs 7.1 mm (text);
  Laval nozzle with 6 mm throat → supersonic jet; Martin is subsonic.
- The only measured impingement h in the thesis is 0.4–1.6 kW/m² K (Ch. 3
  calibration jets); 5–10 kW/m² K were literature assumptions.
- The burner was lowered by the crane with intermittent drawworks pushes; a
  model stall = where the operator pushed.
- Hole width vs depth read from Fig. 8.8: `meier_fig8_8_hole_profile.csv` (±5 mm;
  collar ≈ 95 mm in the top 3 cm, 85–88 mm through the middle, ≈ 80 mm at the
  bottom). This is the hole-shape target for the next packet.

**Stale in the history below:** §1's "30–50 mm impingement zone" and the
wall-jet/cuttings explanation of the hole width; §3a–3c numbers (pre-E1 leak
fix and pre-D2); §3e resolution figures (S1 showed ROP under a prescribed flux
is mesh-converged by construction; the open mesh question is the face-form
uptake of the jet march); §4's "energy partition is a free parameter" (the
jet march fixes the budget; the unknowns are T_nozzle and h) and P4 debris
shielding (withdrawn).

## Run notes (2026-09-15)

- **Overshoot time-step rule (binding for removal runs; not the conduction
  limit).** A step deposits `ΔT_step = q·dt/(ρ·Cp·dz)` into the top cell,
  and that is the firing-temperature overshoot (ρ·Cp = 2.1725e6 J/m³K). For
  ≤ 5 K: dt ≤ 11 ms at dz = 2 mm, q = 2 MW/m²; 1.4 ms at 15 MW/m²; 0.7 ms at
  dz = 0.125 mm, 2 MW/m². D2b used 4 ms (J-M) and 2 ms (J-10).
- **`h_col` is not a flake size** in self-heated drilling. It is a removal
  increment ≈ dz/2 − a_f. See
  `tests/MMWSpalling/studies/s1_surface_resolution/README.md`.

---

# History (pre-D2, June–September 2026; partly stale, see above)

## 1. The experiment (Meier 2017, ETH Diss. 24021, Ch. 8)

A **contactless** flame-jet drilling tool (no mechanical cutting — the burner
descends by gravity feed) drills a 660 kg Grimsel granite block:
38 kW combustion power (HHV), ROP = 1.5 m/h **constant** over 23 min,
V = 3.42 L through 50 cm of block, TSE = 15.4 J/mm³ (combustion-enthalpy
convention), hole ~85 mm Ø cylindrical, cuttings avg < 100 µm
(intra-crystalline spallation), lateral confinement ~1 MPa.

The ~85 mm hole exceeds the ~30–50 mm flame impingement zone because the
deflected hot gases (wall-jet) and entrained cuttings spall/scour the
perimeter — **not** because anything cuts mechanically.

**Caveat for this project:** Meier's source is a convective flame; an MMW
drill is contactless *and volumetric*. This case benchmarks the shared
spallation / removal / thermomechanics kernel with the MMW beam used as a
surface-flux surrogate. Native MMW validation (Beer-Lambert depth, α(T)
feedback) is separate and still ahead.

## 2. What was built (chronological)

1. **Uniform (flat-top) beam profile** — `beam.profile = uniform`,
   `beam.radius` (`src/Numeric/MMWBeam.H`, `BeamSource()` in
   `src/Integrator/MMWSpalling.H`). A Gaussian's central hot-spot pins the
   surface at the 5000 K enthalpy-table ceiling; the flat profile matches a
   flame's stagnation-zone flux.
2. **Sp-gated mechanical bit feed (Step 19)** — `bit.*` keys. A flat floor
   descending at the drawworks rate, removing rock only where the spall
   criterion fired. Used in the first scored drilling case; **superseded by
   Step 20** (kept for prescribed-rate studies).
3. **Plot-cadence fix** — `Integrator.cpp` `Evolve()`: unset `plot_dt`
   (default −1.0) wrote a phantom plotfile every whole second of sim time.
4. **Step 20: surface-normal (slope-aware) spall kernel** —
   `spall.surface_normal = 1`. The K_I criterion is evaluated along the local
   cavity-wall normal (vertical T profile sampled at depth/cosθ) and a flake
   of normal thickness h recedes its column by h/cosθ — the graph form of
   level-set normal motion. The `flake_coherence_length` cap is bypassed
   under the flag (the slope-aware law is self-regularizing; the cap is what
   froze the old runs into a 63° cone). Flag off = bit-identical to all
   prior behaviour (verified against Kant onset output to machine precision).
5. **Beer-Lambert energy-conservation fix** — `BeamSource` previously sampled
   the volumetric deposition pointwise at cell centres, silently losing
   energy once 1/α < dz: at the Meier surrogate's α·dz = 4 only **54%** of
   nominal beam power reached the rock. Now integrated analytically per cell
   (exact at any α·dz). Kant (no beam) unaffected; Hu family (α·dz ≈ 0.05)
   shifts < 0.1%; **the Meier scored baselines below are invalidated and
   await re-baselining** (TODO item R1).

## 3. Result history — read in order

### 3a. Committed scored test (pre-source-fix; numbers now stale)
`sp_meier_pilot/test` (case A energy, case B gated-bit drilling): all 8
checks passed — dV/dt −6%, TSE(HHV) −3%, emergent ROP 1.44 vs 1.5 m/h,
no decay, Sp ≥ 1, no runaway. **Caveat discovered later:** these runs
physically received ~54% of nominal beam power (item 5 above), so part of
that agreement rode on the bug. The test is retained but must be
re-baselined after the TODO physics program (it will currently fail).

### 3b. Stationary-beam cone stall → Step 20 fix
Without a bit, the vertical-only kernel froze into a self-limiting 63° cone
(the slope cap's maximum) at ~48 mm and baked to the 5000 K ceiling, while
rim rock sat at Sp ≈ 6 unremoved — a strictly-vertical-removal artifact (the
hot rim would physically spall hardest). With Step 20: **emergent cylinder**
at beam diameter (r = 24.7 mm constant over 80+ mm of depth), no stall, no
bake, the floor self-organizes onto Sp = 1.00 (marginal-stability attractor)
at ~850 K ≈ 570 °C — inside Kant's measured 553–694 °C spall band. **No bit
needed for sustained drilling.**

### 3c. Pure no-bit prediction (post-source-fix, current physics state)
3D, 42 mm footprint, 30% of 38 kW into rock (packet's efficiency estimate),
no pacing device, 56 s:

| Quantity | Model | Meier | ratio |
|---|---|---|---|
| ROP (16–56 s) | 2.70 m/h | 1.5 | 1.8× |
| dV/dt | 4.62 cm³/s | 2.47 | 1.9× |
| TSE (to-rock) | 2.3 J/mm³ | 4–6 | ~0.5× |
| Surface | Sp = 1.00, ~854 K | spall regime | ✓ |

The three discrepancies are one: the model removes rock ~1.9× too
efficiently per delivered joule — equivalently, matching 1.5 m/h implies a
delivery efficiency of **~17%** (packet's plausible range 20–50%). An
uncalibrated first-principles prediction within a factor of 2, with the
spallation regime, threshold-riding surface, steady advance, and cylindrical
shape reproduced. Candidate explanations for the factor: delivery efficiency
genuinely below 20%; missing cuttings-shielding physics (TODO P4); flux
profile structure (TODO P10); resolution (below). Mild within-run
deceleration (3.6 → 2.3 m/h) is unexplained — watch in re-baseline.

### 3d. 2D slab surrogate — validated
A 4-cell-thick slab (70×4×60) reproduces the 3D ROP within ~13–20% with
identical physics signatures at **~17× less compute** (~1 min vs ~10 min).
Geometry-bound quantities (V, TSE vs Meier's cylinder numbers) do NOT carry
over — trench ≠ hole. This is the development harness
(`sp_meier_pilot/input_2d_dev`); 3D is reserved for final scored runs.

### 3e. Resolution sensitivity — quantified, unconverged
1-D columns at fixed 2.06 MW/m²: sustained ROP 1.80 / 1.62 / 0.99 m/h at
dz = 2 / 1 / 0.5 mm. In the 2D Meier config the fine grid came out *faster*
instead (+30%). Sensitivity ≈ ±30–40%, **non-monotonic across
configurations, not converged at 0.5 mm**. Treat dz = 2 mm ROP as carrying a
factor ~1.4 numerical uncertainty (TODO N1).

## 4. Current limitations (each maps to a TODO item)

- **No melt mechanics**: the criterion would "spall" liquid; nothing
  evacuates melt (P1, P2). Blocks native MMW validation.
- **No cuttings/debris shielding**: removal is instantaneous and clean;
  real flakes blanket the surface (P4) — prime suspect for the 1.9×
  over-efficiency.
- **Hard binary threshold**: no subcritical (stress-corrosion) crack growth;
  marginal states (Sp ≈ 0.9) persist indefinitely where real granite fails
  in seconds (P5).
- **Pointwise local-thermoelastic stress**: correct for the advancing floor
  (it *is* Kant's half-space), overestimates confinement at cavity walls,
  no stress concentrations, cannot test block splitting (P6, P7).
- **Column (height-function) surface**: walls can be vertical but never
  overhang — no lateral hole widening below the rim (Meier's 85 mm vs
  50 mm footprint is out of scope by construction).
- **Energy partition is a free parameter** (~17% inferred vs 20–50%
  estimated); the model validates the product (efficiency × spallation
  kernel) (P10).
- **Resolution unconverged** at production dz (N1).
- **Scored test stale** post-source-fix (R1).

## 5. Files

| File | Role |
|---|---|
| `sp_meier_pilot/input_2d_dev` | **fast development harness** (~1 min) — see TODO.md workflow |
| `sp_meier_pilot/input_energy` | scored case A (38 kW, free spallation) — stale, awaiting R1 |
| `sp_meier_pilot/input_drilling` | scored case B (gated bit) — superseded by no-bit prediction, awaiting R1 |
| `sp_meier_pilot/input_drilling_25mm` | footprint sensitivity (pre-Step-20 funnel documented) |
| `sp_meier_pilot/test` | scoring harness (8 checks) — re-baseline under R1 |
| `TODO.md` | physics program + development workflow |
