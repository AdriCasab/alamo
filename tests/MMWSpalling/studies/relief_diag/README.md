# The relief diagnostic — is R = 0.531 geometry or material? (2026-11-05, advisor-run)

Run at the user's request ("can you just run it yourself right now"), no packet. Question: the criterion's stress is
the fully constrained 1-D thermal stress times a relief factor R = 0.531294 (`spallation.sp.stress_relief`,
CAL-ONCE on Kant). Does the elastic operator, solving equilibrium on the real finite geometry, deliver that
relief by itself (then R is geometry and retires), or not (then R hides a material mechanism)?

## What ran (`run.py`, `run_uniform.py`; argv in each leg's `argv.sh`; binary in `BINARY.txt`)

**Spot case** — the Höser lab flame (P2-16 HOS_R argv) on a FULL 100 × 100 mm block (100 mm deep), spot centred
(r < 20 mm), free sides, bottom roller (`zlo_roller_321`, the only elastic BC preset), removal on through the
column path, to 1.0 s (nothing fires before 1.0 s on these legs):

| leg | stress | surface T | Robin form | mesh |
|---|---|---|---|---|
| L | local 1-D law × R | reconstructed | surface (production form) | 2.08 mm |
| Lc, Lc1 | local 1-D law × R | cell average | cell | 2.08 / 1.04 mm |
| M, M1 | **the operator** (no relief at that site), el.interval 10 | cell average | cell | 2.08 / 1.04 mm |

The interlock chain forced the cell form on the comparison legs: `robin_form = surface` requires
`surface_T = reconstructed`, the flake energy form requires it too, and the operator refuses it (`MMWSpalling.H`
:791, :5680, :5701). Both comparison legs therefore see the same thermal BC and the same cell temperature, and
**ratio = σ_op / σ_1D = R · Sp_M / Sp_L** column by column (same flaw, same K_Ic(T), same T).

**Resolution ladder** (`U/`) — the `unit/face_removal` fixture (uniform 1866 K gas, h 700, 24 × 24 × 64 mm block,
free sides) heated uniformly over its whole top, at 2 / 1 / 0.5 mm. Exact answer at the centre, twelve skin
depths from any edge: full constraint, ratio = 1. Whatever the operator reads below 1 is its own discretisation
error for a skin of ℓ/dx cells.

Solver notes: a 50³ mesh with 25-cell boxes coarsens only once and **MLMG failed** on the first heated solve;
48³ / 96³ with 24-cell boxes (coarsenable four times) converges in 17 iterations, 0.5 s (48³) / 3.9 s (96³) per
solve. `el.interval = 1` had cost > 1 min per step on the bad mesh; at interval 10 the stress is quasi-static.

## Results (`RESULTS.txt`, `RESULTS1.txt`, `RESULTS_uniform.txt`)

Operator / 1-D law at the spot centre (r = 0) and at the centre of the uniformly heated block, against the
skin depth in cells (two-cell estimate on the centre column):

| ℓ/dx | spot centre, 2 mm | spot centre, 1 mm | uniform, 2 mm | uniform, 1 mm | uniform, 0.5 mm |
|---|---|---|---|---|---|
| 0.30 | 0.31 | — | 0.26 | — | — |
| 0.45–0.5 | 0.39–0.42 | 0.39 | 0.29–0.30 | 0.45 | — |
| 0.63–0.67 | — | 0.44 | — | 0.45 | — |
| 0.8–0.9 | — | 0.51–0.55 | — | 0.49–0.51 | — |
| 0.96–1.0 | — | — | — | 0.51 | 0.59 |
| 1.3–2.0 | — | — | — | — | 0.64–0.65 |

Rim: the operator puts the spot's rim (r ≈ 24 mm) in mild tension, Sp_M = −0.03 at 1 s (K_I 3 % of K_Ic).

## Reading

1. **The operator's surface stress is not converged at any affordable mesh.** On the uniform block, where the
   answer is 1, it reads 0.26–0.30 at 2 mm, 0.45–0.51 at 1 mm and 0.59–0.65 at 0.5 mm, and the 0.5 mm value
   plateaus at 0.65 while ℓ/dx grows from 1.6 to 2.0. The error is a function of the mesh, not only of ℓ/dx.
   ⇒ **for a flame skin (ℓ ≈ 1 mm) the operator is unusable as an absolute stress source at production meshes.**
   The 1-D law remains the only usable stress on the flame path. The operator is for geometries whose heated
   region is resolved — centimetre-scale beam heating, which is exactly where it is needed.
2. **Against its own uniform baseline, the spot centre shows no geometric relief.** At matched mesh and ℓ/dx the
   spot centre reads 1.2–1.4 × the uniform value at 2 mm and 0.9–1.1 × at 1 mm. A 40 mm spot on a 100 mm block
   is as constrained at its centre as an infinite uniformly heated slab, to within ±15 %. The discretisation
   error cancels to first order in that comparison.
3. ⇒ **Decision-tree outcome 2: the 47 % missing stress is NOT geometry.** R is not constraint the half-space
   law assumed and the block does not provide. It stays a DECLARED, unexplained level factor until a material
   mechanism replaces it: thermal microcracking / stiffness loss (check the E(T) table is applied at T_s),
   grain-scale expansion mismatch (the Voronoi study the user kept), or the flaw statistics. It must not be
   re-fitted on any beam case.
4. The rim tension exists and is small at 1 s (K_I ≈ 3 % of K_Ic): radial cracking is not predicted at this
   exposure on this block; Woskov's banding was for far hotter spots.

## Caveats
- The comparison legs use the cell Robin form and the cell temperature (forced by the interlocks); the
  production form L is reported only for its onset (fires at 1.00 s, 817.5 K at 2 mm).
- The operator's own convergence plateau (0.65 at 0.5 mm) is unexplained: node/cell averaging across a one- to
  two-cell skin plus the free-surface stress recovery are the suspects; a 0.25 mm uniform leg (48 × 48 × 128 mm
  at 0.25 mm = 1.9 M cells) would say whether it converges at all. Not run.
- `el.interval = 10` (one solve per 80 ms); the criterion reads an up-to-80-ms-old stress between solves.
