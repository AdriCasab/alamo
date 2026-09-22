# D2j-0b — the hot disk under a descending exclusion plane (written before launch; never edit)

Reviewer diagnostic, 2026-09-21 evening. Key-only: binary `6d047b50…` (the D2j-0 / D2i
binary), no source edit, no build. Purpose: test the mechanism read from the D2i cascade
timing (`Claude_markdowns/2026-09-21.md`, "After D2j-0") **without the feet and without the
jet**: does the `s <= 0` exclusion (no flux to rock that has risen above the nozzle plane)
turn the bounded D2j-0 edge offset into an inward, time-growing, mesh-dependent stall front?

## Configuration (all runs)

D2j-0 `input_hotdisk` unchanged below the patch block. Quarter domain 0.08 × 0.08 × **0.15** m
(deeper than D2j-0 so 250 s of centre recession, ≈ 108 mm, leaves ≈ 40 mm of rock), axis at the
xlo/ylo corner, surface at z = 0.15 m. Disk radius 0.030 m, h_in = 700 W/m²K, T_in = 1600 K
(the D2j-0 sharp set). Pinned closure, idle clock off, V0 = dz³, no jet closure, no feet.
A **prescribed plane** `nozzle_descent = prescribed`, `nozzle_z0 = 0.170` m (20 mm above the
surface), `nozzle_feed = 4.3056e-4` m/s (1.55 m/h, just under the D2j-0 centre rate 1.58 m/h
so the centre never reaches the plane), `nozzle_collision_radius = 0.002` (the axis column(s)
abort if the plane ever reaches them; nothing else is guarded). s = z_plane − z_face per column.
dt 16 ms at 2 mm, 8 ms at 1 mm (the D2i ladder, 2 K/step), `blocking_factor = 1`, mgs 20 / 40,
stop 250 s, plotfiles every 50 s, removal_events.csv on.

## Variants (h and T_gas as parser expressions in r-independent form; s is the only variable)

| case | above the plane (s ≤ 0) | below (s > 0) | what it isolates |
|---|---|---|---|
| **C0** | h 700, T 1600 (no exclusion) | h 700, T 1600 | the D2j-0 edge with a plane that does nothing: the kinematic baseline |
| **C1** | h **1e-3**, T 1600 (exclusion emulated; the parser forbids h ≤ 0 in the patch) | h 700, T 1600 | the exclusion alone: exposed rock cools to the bulk and to lab air |
| **C2** | h **100**, T **1000** (exposed rock held warm, equilibrium ≈ 780 K, below firing) | h 700, T 1600 | exclusion + a warm wall: does keeping the exposed rock warm remove the cascade? |
| C3 (2 mm only, exploratory) | h 100, T 1600 (exhaust-like wall; expected to spall the wall slowly) | h 700, T 1600 | what a physically heated wall does to the hole; not scored |
| C1d (2 mm only) | as C1 at dt 8 ms | | dt check on the scored variant |

Runs: C0/C1/C2/C3/C1d at 2 mm (minutes each); C0/C1/C2 at 1 mm (≈ 1–2 h each).

## Kinematic baseline (hand, from the D2j-0 S2/S1 profiles; not a prediction under test)

With the plane at the centre rate, a band with v/v_c = f crosses the plane after
t = s0/((1 − f) v_c) = 20 mm/((1 − f)·0.439 mm/s): 28–30 mm ≈ 51 s; 26–28 ≈ 71–73 s;
24–26 ≈ 152 s (2 mm) / 111 s (1 mm); 22–24 ≈ 268 / 163 s; 20–22 ≈ 351 / 217 s. So **a purely
kinematic sweep already exposes bands inward, faster at 1 mm** because the 1 mm profile is
slower at every radius. The cascade signature is therefore NOT "the front moves inward" but
**"a band stops firing before its own face crosses the plane, and its stop follows the crossing
of the band outside it"**, as in D2i.

## Metrics (fixed now)

From `<pf>_removal_events.csv` (every removal with column and `h_applied`) and the known plane
z_n(t) = 0.170 − 4.3056e-4·t. Per column: cumulative recession c(t), face z_f(t) = 0.15 − c(t)
(cell-quantised as the code does: the top face of the top solid cell), s(t) = z_n(t) − z_f(t).
Bands: 2 mm annuli in column-centre radius, 18–30 mm.

- **v_c**: mean column rate over r < 6 mm, 50–250 s.
- **t_stop(band)**: the median over the band's columns of each column's last removal time
  (columns still firing in the last 5 s count as t_end). **t_cross(band)**: the median over the
  band's columns of the first time s(t) ≤ 0 (∞ if never). **Δ_own = t_cross − t_stop**;
  **Δ_outer = t_stop(band) − t_cross(next band outward)**. **depth-below-plane at stop** = −s at
  t_stop.
- **Disk-mean rate** R(block) over r < 30 mm per 50 s block, and the mesh gap
  g(block) = R_1mm/R_2mm − 1.
- **r_front(t)**: inner edge of the contiguous run of 2 mm bands reaching 30 mm whose mean rate
  over the trailing 25 s is < 0.1 v_c, at t = 100, 150, 200, 250 s.

## Predictions

- **P1 (sanity)**: v_c within 5 % of 1.58 m/h in every run; C0 at each mesh reproduces the
  D2j-0 S2/S1 profile (disk-mean 50–150 s within 3 % of 1.133 / 1.014 m/h). If P1 fails the
  configuration is wrong; stop.
- **P2 (the cascade)**: in **C1**, at both meshes, at least two bands inside 28 mm stop with
  Δ_own > 10 s (they stop while still below the plane) and Δ_outer in 5–60 s; in **C0** no band
  inside 28 mm stops before 250 s (Δ_own undefined; t_stop = t_end). Failure of the C1 half
  means the exclusion alone does not cascade and the feet/jet feedback is essential (→ run the
  D2h OL pair, as the D2j-0 RESULTS said).
- **P3 (the growing mesh gap)**: in C1 the gap g grows monotonically in magnitude across the
  50 s blocks and exceeds −20 % in the last block (D2i reached −36 %); in C0 g stays within
  ±4 percentage points of its 50–100 s value (the standing ≈ −10 % D2j-0 offset).
- **P4 (the remedy)**: **C2** behaves like C0 under P2 and P3 (no band stops with Δ_own > 10 s;
  g standing). If P4 fails while P2/P3 hold, warming the exposed rock is not the fix and the
  corner needs the geometric remedy (the continuous front, `2026-09-21b.md`).
- **dt**: C1d matches C1 (2 mm) in every band's t_stop to within 5 s and in R(block) to 2 %.

## Decision rule

P1 ∧ P2 ∧ P3 ∧ P4 → the cascade is the exclusion cooling the exposed rock, and a wall-heating
closure above the plane is the justified next packet, with C2's outcome as its acceptance
prediction. P2 fails → feedback packet (D2h OL pair). P4 fails → front redesign, not wall
heating.
