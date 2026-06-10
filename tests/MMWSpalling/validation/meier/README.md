# Meier 2017 Grimsel pilot — validation case

Status: **closed out — scored regression, passing** (2026-06-10). Two scored
cases run by `sp_meier_pilot/test`; known kernel limitations documented below
and in the input headers.

Reference data packet: `meier-2017-pilot-660kg-validation-reference.md` (repo
root). Inputs and the scoring script live under `sp_meier_pilot/`.

---

## 1. What this case is — and an important caveat

Meier 2017 (ETH Diss. 24021, Ch. 8) is the meter-scale **thermal-spallation
drilling demonstration**: a 38 kW methane/air **flame jet** drills a 660 kg
Grimsel granite block at a steady **ROP = 1.5 m/h**, removing **3.42 L** at a
thermal specific energy **TSE = 15.4 J/mm³**, leaving an ~85 mm cylindrical
hole.

**Caveat (the big one).** Meier's heat source is a **convective flame jet with a
mechanical bit** — contact drilling, surface heating. An **MMW drill is the
opposite: contactless and volumetric** (the wave penetrates and deposits energy
over a Beer-Lambert depth, no bit). So this case does **not** directly validate
the MMW beam. Its value is as an **end-to-end benchmark of the shared
spallation/removal/thermomechanics kernel**, with the MMW beam used as a
*surface-flux surrogate* for the flame. The native MMW-beam validation belongs
to the Zhang/Oglesby (MMW thermal) and Hu (spallation) cases.

## 2. Key Meier targets (from the reference packet)

| # | Target | Value | Scored here? |
|---|---|---|---|
| 1 | ROP, sustained, no decay | 1.5 m/h | **yes** (case B: emergent, ±20%, no-decay) |
| 2 | TSE | 15.4 J/mm³ (HHV) / 4–6 (to-rock) | **yes** (case A: HHV ±30%; case B: to-rock band) |
| 3 | Excavated volume rate | ≈2.47 cm³/s | **yes** (case A early window, ±30%) |
| 4 | Hole shape | cylindrical, ~85 mm Ø | **no — kernel limitation** (see §6) |
| 5 | No axial splitting at SHmin ≈ 1 MPa | none | implicitly (no splitting observed) |

Grimsel granite (Meier Table 7.1): k=1.5 W/mK, ρ=2750, Cp=790, E=30 GPa, ν=0.3,
α=8e-6 /K. Failure: LEFM `Sp = K_I/K_Ic(T)`, a0=20 µm, K_Ic = Nasseri-2007
Westerly table (no Grimsel-specific value exists — documented substitution),
Weibull m=20.

## 3. Modelling approach and decisions

Heat source: the **MMW beam used as a surface-flux surrogate** — high
Beer-Lambert absorption (`beam.alpha_expr = 2000` → 0.5 mm e-folding) so ~all
power lands at the surface, and the beam is **surface-following** (it tracks
the receding floor). Chosen over the `surface_patch` convective BC because that
BC is pinned to the fixed top face and cannot follow a deep advance.

Decisions taken during bring-up:
- **Profile: uniform flat-top disk**, not Gaussian. A flame jet's stagnation
  heat flux is flat; a Gaussian's central hot-spot (~60 MW/m² peak) overheats
  the centre and pins the surface at the enthalpy-table ceiling (5000 K, see
  `src/Numeric/Material/Table.H`). → added `beam.profile = uniform`.
- **Footprint radius is inferred, not measured** (Meier gives only nozzle
  Ø=7.1 mm, SOD=7D). 25 mm = thermal impingement footprint; **42 mm = the
  as-drilled 85 mm hole** (the swept heat-delivery area). 42 mm is the scored
  drilling case; 25 mm is kept as a sensitivity (`input_drilling_25mm`).
- **Energy partition: two conventions, scored separately.** Case A feeds the
  full 38 kW and scores in Meier's HHV convention (valid only in the early
  window, before the cone-stall). Case B feeds the packet's
  efficiency-corrected **30 % (11.4 kW)** — full power with only spallation as
  a sink runs away, because ~50–80 % of Meier's combustion enthalpy leaves
  with the exhaust (a loss sink the model lacks) — and scores against the
  **enthalpy-to-rock TSE ≈ 4–6 J/mm³**.
- **Mechanical bit feed, gated on spallation**: a flat bit floor descends at
  the measured feed (1.5 m/h) and removes rock only where `Sp ≥ 1` fired that
  step, so **ROP is an emergent prediction** (capped at the feed), not an
  assumption. → added `bit.*` (Step 19).

## 4. Code added (all default-off; no effect on existing tests)

1. **Uniform beam profile** — `src/Numeric/MMWBeam.H` and the surface-following
   `BeamSource()` in `src/Integrator/MMWSpalling.H`.
   Keys: `beam.profile = gaussian|uniform`, `beam.radius`.
   Uniform: `I = P0/(π r²)` inside the disk, 0 outside, no Gaussian waist.
2. **Mechanical bit feed (Step 19)** — `ParseBitSettings` in
   `src/Integrator/MMWSpalling.H`, logic in `src/Integrator/MMWSpalling/Removal.H`.
   Keys: `bit.enabled`, `bit.feed_rate` (m/s), `bit.radius`, `bit.x0/y0/z0`,
   `bit.start_time`. A flat floor `z_bit(t)=bit_z0-feed·(t-start)` descends at
   the feed; within `bit.radius` the surface is forced down to it **gated on
   the column firing spall this step**; the bit is the sole remover when
   enabled (regime code 3). Requires `spall.enabled` or `vapor.enabled`.
3. **Plot-cadence fix** — `src/Integrator/Integrator.cpp` `Evolve()`: the
   `plot_dt` branch lacked a `plot_dt > 0` guard; with the default −1.0 it
   wrote a phantom plotfile every whole second of sim time.

## 5. Files (`sp_meier_pilot/`)

| File | What it is |
|---|---|
| `input_energy` | **scored case A** — 38 kW HHV, free spallation, no bit, 12 s |
| `input_drilling` | **scored case B** — 11.4 kW, Sp-gated bit at 1.5 m/h, r=42 mm, 56 s |
| `input_drilling_25mm` | unscored sensitivity — 25 mm thermal footprint (funnels; header documents why) |
| `test` | runs + scores both cases (`--no-run` to score existing output) |

```bash
/Users/tzetze20/Desktop/code/.venv/bin/python \
  tests/MMWSpalling/validation/meier/sp_meier_pilot/test
```

## 6. Results (2026-06-10, all checks PASS)

**Case A — energy/TSE (38 kW, free spallation, t ≤ 10 s):**
dV/dt = 2.32 cm³/s vs 2.47 (−6 %); TSE = 15.0 J/mm³ vs 15.4 (−3 %); hole
diameter 49.6 mm ≈ the 50 mm heated disk.

**Case B — sustained drilling (11.4 kW, gated bit, r = 42 mm, 56 s):**
emergent ROP = 1.44 m/h vs 1.5 (−4 %) with no decay (2nd/1st-half rate ratio
0.90); floor spall-ready throughout (Sp 6.3–9.7); TSE = 6.35 J/mm³ (band
2.8–7.8); Tmax = 1747 K — no runaway, thin hot skin maintained all the way
down.

**Known kernel limitation (documented, not scored).** Without the bit, the
stationary beam bores a self-limiting cone that stalls (~48 mm) and then bakes:
the **column-based removal scheme cannot spall the sloped cavity wall** — rim
columns sit at Sp ≈ 6 and bake instead of flaking off, a code artifact of
strictly-vertical removal (real spallation pops flakes along the surface
normal; the hot, overstressed rim would physically spall *hardest*). The same
mechanism makes the 25 mm case taper into a funnel. Cylindrical shape (#4) and
full-duration volume (#3 over 1383 s) are therefore out of scope until a
**surface-normal removal kernel** exists (planned next physics step; `phi` is
already a level set, so the natural formulation is `∂φ/∂t = −v_n|∇φ|`).

## 7. Recommended next steps

1. **Surface-normal spall removal** — the single change that unblocks targets
   #3/#4, retires the `flake_coherence_length` regularizer as a load-bearing
   constraint, and likely removes the need for the bit stand-in. Then re-run
   Meier as a true prediction (no bit) and re-score.
2. **Native MMW validation (the real project goal):** volumetric Beer-Lambert
   heating at a realistic absorption depth, contactless, no bit — against
   Zhang/Oglesby and Hu, then an MMW experimental anchor (e.g. Woskov).
