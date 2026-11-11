# Meier case — physics work program and development workflow

Updated 2026-06-10. Companion to [README.md](README.md) (status + limitations).

## Development workflow

The loop, per physics item below:

1. **Implement** the item (default-off input flag unless stated otherwise, so
   every existing test stays bit-identical).
2. **Run the 2D dev harness** — `sp_meier_pilot/input_2d_dev`, ~1 min on
   4 ranks:
   ```bash
   mpirun --oversubscribe --bind-to none -np 4 \
     bin/mmwspalling-3d-g++ tests/MMWSpalling/validation/meier/sp_meier_pilot/input_2d_dev
   ```
3. **Compare against the harness baseline** (recorded in the input header):
   central-column ROP, max Sp at the floor (should ride ≈ 1.00), Tmax
   (no runaway), trench wall shape. Decide whether the change moved the
   physics in the expected direction; update the baseline table in the input
   header when a change is accepted.
4. When the program is complete and the 2D results are satisfactory,
   **re-baseline the 3D scored test** (item R1) and optionally run the hero
   case (R2).

The 2D slab is validated against 3D to ~13–20% on ROP at ~17× less compute
(README §3d). Geometry-bound quantities (V, TSE vs Meier's cylinder numbers)
must NOT be scored in 2D — trench ≠ hole.

---

## P. Physics items (priority order)

### P1. Melt-aware spallation  — blocks native MMW work
Spallation requires a solid: a melt cannot store elastic energy or carry a
crack. Currently the criterion would "spall" liquid (K_Ic table clamps small
and finite; phase fractions Λ_S/Λ_L exist thermally but the criterion ignores
them).
- Scale the driving stress by a load-bearing skeleton function f(Λ_L): 1 in
  solid, → 0 near the rheologically critical melt fraction (Λ_L ≈ 0.4,
  rigidity percolation), 0 in full melt. (Equivalent to E(Λ_L) in the
  closed form; interpretation: no skeleton, no stored stress.)
- Anchor the depth scan at the topmost load-bearing (majority-solid) cell:
  melt above is overburden fluid; the solid–liquid interface is the LEFM
  free surface.
- Acceptance (2D): with melting enabled and flux pushed up, drilling
  transitions from spallation to a melt-limited regime instead of spalling
  liquid; with melting disabled, bit-identical.

### P2. Melt evacuation
With P1 alone, melt accumulates as an insulating spall-proof blanket and the
hole stalls — physically real (glassy plugs) but real MMW systems remove melt
by vaporization (modelled) and purge-gas shear (not modelled). Add a
purge/shear removal channel for liquid-fraction surface cells (rate
parameterized; calibration deferred).
- Acceptance (2D): melt-regime drilling reaches a steady ROP set by the
  purge parameter instead of stalling.

### P3. Temperature-dependent MMW absorption α(T)
The defining MMW feedback: granite's loss tangent rises steeply with T →
absorption runaway localizes deposition. `beam.alpha_expr(T)` plumbing
exists; needs a literature α(T) for granite/basalt at λ ≈ 2 mm and a
validation anchor (Woskov-style penetration data). This is where the model
stops being a flame surrogate and becomes an MMW model.

### P4. Cuttings / debris shielding  — promoted by the 1.9× over-prediction
Removal is currently instantaneous and clean; real flakes blanket the surface
and intercept the incoming flux until evacuated (Meier had to add a vacuum
cleaner mid-test). Prime suspect for the model removing rock ~1.9× too
efficiently per delivered joule (README §3c). Simplest model: a shielding
efficiency factor on the incident flux, parameterized by the recent removal
rate; or an explicit debris-residence time.
- Acceptance (2D): predicted ROP at 30% power moves from ~2.9 toward
  ~1.5 m/h with a physically defensible shielding parameter.

### P5. Subcritical crack growth (stress corrosion)
The criterion is a hard binary threshold; marginal states (Sp ≈ 0.9) persist
indefinitely where real granite fails in seconds (v ∝ (K/K_Ic)^n, n ≈ 30–50).
Softens the threshold from below; matters most near onset and at low flux.
- Acceptance (2D + 1-D): onset time at sub-threshold flux becomes finite and
  n-dependent; threshold-riding states fire at Sp slightly below 1.

### P6. Wall stress anisotropy
`local_thermoelastic` keeps full biaxial confinement everywhere; at a cavity
wall the surface-normal component is relieved (closer to uniaxial
tangential). Step 20 already computes the local normal per column — degrade
the biaxial form toward uniaxial with slope (exact in flat and vertical
limits). Reduces over-driving of wall spallation.

### P7. Elastic solve A/B and revival
The pointwise stress form cannot see stress concentrations, block-scale
splitting (Meier target #5 currently passes vacuously), confinement coupling
beyond the additive offset, or Voronoi grain-scale heterogeneity. The MLMG
solve was shelved for convergence reasons (see
mlmg-shelved-local-thermoelastic-notes.md).
- First step is cheap: A/B a flat-onset case (solve vs closed form) and a
  pre-drilled-hole wall-stress comparison to quantify the error being
  carried. Revival decision afterwards.

### P8. Radiation in the melt regime
Negligible at spall temperatures (~24 kW/m² at 850 K) but ~MW/m² at melt
temperatures — same order as the incident flux. Bundle with P1–P3; the
existing `losses.*` framework needs a cavity view-factor story eventually.

---

## N. Numerics items

### N1. Resolution convergence study  — standing TODO, quantified
Sustained ROP varies ±30–40% over dz = 2 → 0.5 mm and is **non-monotonic
across configurations** (1-D: fine slower; 2D Meier config: fine faster).
Not converged at 0.5 mm. Production dz = 2 mm numbers carry ~factor-1.4
numerical uncertainty. Needed: a proper sequence (2D, dz = 2/1/0.5/0.25 mm,
fixed physics) + understanding of the flake-quantization mechanism; consider
a sub-cell surface-temperature reconstruction in the criterion if bias
persists.

### N2. Enthalpy-table ceiling as input
`Numeric/Material/Table.H` hard-codes T_max = 5000 K; saturation silently
violates energy conservation. Make it an input and abort-or-warn on
saturation.

### N3. Implicit/long-horizon pathway
1383 s at explicit dt = 1 ms is ~10⁶ steps. The implicit thermal path exists
but is not compatible with the removal machinery (CLAUDE.md lesson 5).
Needed only for hero runs (R2); revisit then.

---

## R. Re-baselining (after the physics program)

### R1. Re-baseline the 3D scored test
The Beer-Lambert fix invalidated the committed bands (runs received ~54% of
nominal power — README §2.5). After P-items land: re-run `input_energy`,
replace the gated-bit case with the **no-bit prediction** scored in the
implied-efficiency framing (option 2: the efficiency at which the model
matches 1.5 m/h must fall in the plausible 20–50% window), bands widened per
N1. Add the emergent-cylindricality check (Step 20 makes shape scoreable).
Investigate the mild within-run deceleration (3.6 → 2.3 m/h) seen in §3c.

### R2. Hero run (optional)
Full 1383 s / 50 cm advance at coarse dx in 3D (or 2D first): tests the
long-horizon no-decay claim that 56 s windows cannot. ~10–20 h wall-clock;
overnight job; requires N3 or patience.

---

## Done (for the record)
- Uniform beam profile; Sp-gated bit feed (Step 19); plot-cadence fix.
- Step 20 surface-normal kernel + coherence-cap bypass (emergent cylinder,
  no stall, threshold-riding floor).
- Beer-Lambert per-cell energy-conserving deposition (54%-loss bug fixed).
- 2D slab surrogate validated (~13–20% on ROP, ~17× cheaper).
- First-principles no-bit prediction: factor ~1.8 high on ROP ⇒ implied
  delivery efficiency ~17%.
