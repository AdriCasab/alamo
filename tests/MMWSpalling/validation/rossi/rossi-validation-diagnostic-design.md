# Rossi 2018 damage-profile validation — diagnostic design

Step 16 shipped the `h_spall` mechanics under `sp_weibull` and deferred the
Rossi 2018 quantitative validation. This note captures the design
decisions for the follow-up packet that wires the Rossi-comparable
diagnostic. It belongs in the repo root next to
[kant-2017-central-aare-validation-reference.md](kant-2017-central-aare-validation-reference.md)
because, like Kant, this defines how a specific external dataset maps
onto an ALAMO observable.

## 1. Diagnostic gap (the conceptual error in Step 16's Rossi test)

Rossi 2018 measures, post-cooling, the **spatial distribution of distinct
cracks** as a function of depth (peak at 100–200 µm, falloff to baseline
by ~520 µm, density 4.4–7.6 mm/mm² at peak). Each crack contributes once
to one depth bin.

Step 16's `h_spall_field` records, per surface column, the **deepest**
`K_I ≥ K_Ic` crossing observed across all steps (running max anchored at
the immovable initial top-z slice). Under a monotonically intensifying
ramp, this is the LATE-ramp value, not the integrated population. A
column whose surface fires early (h_col ≈ 0) and grows by the end
(h_col ≈ 1.5 mm) registers only the late value.

These are different physical quantities. Step 16's test reads it
correctly; the packet's claim that this maps to Rossi's depth profile was
the design error.

## 2. Design decision: event semantics — edge-triggered, not residence-time

Under a monotonic ramp, once `K_I ≥ K_Ic` at a cell, it stays above. A
naive "count every step that K_I ≥ K_Ic" gives a residence-time map
peaked at the surface (the cell that crosses earliest and stays
super-critical longest). That's not a crack count.

**Spec**: count only the **first crossing** per cell. After that one
event, the cell's count is 1 and never increments again.

Implementation sketch — new field `crack_event_count_mf`, cell-centred
Set::Scalar (binary 0/1 at end of run), registered only when
`spallation_model == SpWeibull`. Each step, for each (i, j, k) in the
heated zone with `crack_event_count_mf(i, j, k) == 0`, evaluate
K_I(z_crack = depth_of_k; a_f at column top) and K_Ic(T at (i, j, k));
set the count to 1 if K_I ≥ K_Ic. Bin by k at end of run → depth
profile directly comparable to Rossi.

## 3. Design decision: stress relief after firing — three options

Rossi's measured distribution implicitly assumes a fired crack releases
local stress (the energy goes into crack opening; the rock around the
crack relaxes). Without stress relief, even the time-integrated count
overcounts deep events — every cell eventually crosses K_Ic threshold
as the ramp progresses.

| Option | Description | Effort | Physics fidelity |
|---|---|---|---|
| (i) Zero-out σ in fired cells | When event count flips 0→1, mark cell "fired" and zero σ going forward. Implement via an eigenstrain or a damage-like stiffness scaling (DON'T zero σ post-solve — breaks force balance). | High (mechanics coupling) | Closest to physical picture |
| (ii) Freeze from re-firing only | Don't touch σ; just clamp event count at 1 per cell. K_I deeper in the column still sees unrelieved σ_xx, so overcounts deep cracks — but the SHAPE of the depth profile is informative. | Low (post-process only) | Shape OK, magnitude inflated |
| (iii) Document as known limitation | No stress relief, count every first crossing, accept the deep-end overcount in the Rossi comparison. | None | Lowest |

**Recommendation**: start with **(ii)**. Simplest, gives a clean depth-
profile shape, zero mechanics coupling. If the resulting profile looks
like a step function (all cells eventually fire) with the
peak-at-100-200-µm structure washed out, escalate to (i). The mechanics
integration of (i) is the substantial cost — but the existing
`damage_law` D-field already provides a stiffness-scaling channel, so
(i) could be implemented as a soft damage that ramps stiffness to 0 in
fired cells, reusing existing plumbing.

## 4. Addendum — cell-quantization (a complement to §3, not a substitute)

Step 16 takeaway #3: at 1.5625 mm cells, sub-cell `h_col` events
(K_I ≥ K_Ic only at z_crack = 0, fails at s = 1) get suppressed because
the depth scan walks in integer-cell strides. This biases the GB/IG
firing ratio to **100 % GB** rather than Rossi's ρ_B/ρ_I ≈ 5–6.

A sub-cell binary search in the depth scan (between s and s+1 once the
crossing is bracketed) would let `h_col` be reported at fractional-cell
granularity, and lets the per-cell event count in §2 fire on cells where
the crossing is sub-dz.

This is a sharpening of the diagnostic, not a replacement of §2. The
follow-up packet should bundle it with the event-count field so the two
share a single mechanics pass.

## 5. Open questions for the next `/plan`

1. **Where does the event-count loop live?** Inside
   `UpdateRemovalAfterCohesive` (alongside the existing column scan,
   reusing the K_I helper) or a new dedicated method? The former is
   minimal-diff; the latter is cleaner if we want different sampling
   cadence (every step vs every Nth step).

2. **Which `a_f` does each cell use?** Three candidates:
   - (a) the surface cell's `flaw_a` (matches Step 16 column scan;
     treats the crack as originating from the surface flaw),
   - (b) the cell's own `flaw_a` (per-cell microstructure flaw — more
     physical for "where do cracks initiate"),
   - (c) the column's running-max `flaw_a`.
   Recommend **(a)** — consistent with Step 16 and the K_I integrator's
   "edge crack from the new surface" interpretation.

3. **Stop_time tuning.** Rossi reports cracks at 100–200 µm peak with
   tail to 520 µm at "~650 °C surface". Under §3 option (ii), the depth
   profile saturates as time progresses. Run to a fixed T_surf rather
   than fixed clock time, e.g. stop_time = `(650 - 20) / 5 = 126 s` for
   the canonical 5 °C/s ramp. (Step 16 used 130 s, close enough.)

4. **Cell-centred or node-centred event field?** The stress field is
   node-centred and gets `NodeToCellAverage`d for the K_I evaluation.
   Cell-centred is consistent with `flaw_a_mf` and matches the depth-
   scan stride. Recommend **cell-centred**.

5. **Geometry — keep Step 16's `50×50×50 mm + lateral roller box +
   full-face heating`, or revisit the patch + cold rim?** Step 16's
   takeaway #1 documents why the patch + free lateral config can't fire
   Sp on this machinery (σ_xx at patch center ≈ 0.3× of 1-D ideal). For
   the follow-up packet, recommend **keep Step 16's geometry** — the
   "central column of an infinite heated zone" approximation is a
   defensible reading of Rossi's far-from-rim measurement region, and
   reopening the patch + free lateral question costs an MLMG
   conditioning investigation that's out of scope. Note this as a
   reading of Rossi rather than a literal reproduction.

## 6. Status reference — Step 16 mechanics (DO NOT re-litigate)

The K_I-vs-K_Ic depth scan, `Sp_cluster_id_mf > 0` detach gate,
`h_spall_field_mf` running-max diagnostic, and the column-sweep material
removal under `sp_weibull` are correct and tested:

- 358 / 1024 top-z columns fired with `h_spall_field > 0` at end of ramp
- 358 cells removed (column-sweep wiring works)
- GB fraction of firings: 100 % (cell-quantization artifact — §4)
- `damage_law` regression sweep byte-identical
- §15c regression bundle PASS

The follow-up packet's scope is the new diagnostic field + the Rossi
PASS criteria using that field. The mechanics half is closed.

## 7. Step 16c findings — the option (L) σ_xx artefact

Step 16c shipped path (b) from §1: per-step h_col CSV log + time-binned
histogram, with σ_xx reconstructed linear cell-to-cell (option (L)).
Result: 528 first-firing events recorded; histogram clusters tightly
at **h_col ∈ [757, 787] µm** (mean 769 µm, std 5.8 µm — all in the bin
centred at 799 µm). R1 sub-criterion (1) FAILs (peak outside [100, 200]
µm); sub-criteria (2) and (3) both PASS (unimodal, sharp falloff).

**Mechanism — option (L) creates a sampling artefact at z = dz/2 − a_f.**
Under linear cell-to-cell σ_xx, the top cell has no adjacent cell ABOVE
to interpolate from, so σ_xx is **flat** in z ∈ [0, dz/2] (extrapolated
from the top-cell centre value). K_I integration over the crack length
(z_crack, z_crack + a_f) is then constant for `z_crack ∈ [0, dz/2 − a_f]`
— all integration samples fall in the flat region. K_I starts dropping
only past dz/2 − a_f, where the integral begins to pull from the
linear-interp region. The bisection brackets (0, dz) and converges to
the START of the K_I drop, i.e. z ≈ dz/2 − a_f for every firing column.
At dz = 1.5625 mm and a_f ∈ [4, 22] µm, this gives h_col ∈ [759, 779]
µm — matching the observed 757-787 µm range to within bisection
tolerance.

**Stop_time tuning does not move the peak.** Truncation analysis on the
captured CSV: cutting at any `stop_time` from 85 s through 126 s leaves
the peak at the same bin (centre 799 µm). The peak is determined by
the σ_xx reconstruction, not by the ramp endpoint. Per the Step 16c
packet's "Expected outcomes" protocol for "peak > 200 µm", step 1
(tighten stop_time) does not apply; step 2 (flag for reviewer) is the
active path.

**Why this is option (L)'s flaw specifically.** At first-firing time
t ≈ 80 s under prescribed-T heating (`T_expr = 293.15 + 5·t` K at the
face), T_face = 693 K but T_top_cell_centre ≈ 625 K (thermal-diffusion
gap across dz/2 = 781 µm). σ_xx scales linearly with (T − T_ref) under
1-D confinement, so σ_face ≈ 1.2 × σ_top_cell. Option (L) discards this
20 % surface-stress jump by treating σ_xx as flat from cell centre to
face. The K_I curve consequently has the (artificial) flat-then-drop
shape that pins the crossing at dz/2 − a_f.

## 8. Step 16d scope — option (a), with a pre-declared decision tree

**Step 16d should run option (a) only — T-derived σ_xx with the
prescribed-T face boundary condition.** Implementation: one lambda
swap. The lambda takes `T_expr(t)` at the face, linearly interpolates T
between the face value and the top-cell-centre value (and between
successive cell centres deeper), then computes
`σ_xx(z) = −E · α_th · (T(z) − T_ref) / (1 − ν)`. Implementation cost
is barely above the (c) "smooth K_Ic" alternative.

**Why NOT path (c) — smooth K_Ic between cells — as the next move.** K_I
is determined by the σ_xx profile it integrates over. In the flat
[0, dz/2] region, K_I is constant regardless of what K_Ic does. The
crossing's depth is set by where K_I starts to drop (i.e. by the σ_xx
profile's transition), not by K_Ic's continuity. K_Ic typically varies
only 5–15 % across the top half-cell under Nasseri — far too weak to
move the peak by ~600 µm. Path (c) would burn a packet on the wrong
variable; the mechanism analysis already tells us σ_xx is the
load-bearing knob.

**Step 16d pre-declared decision tree:**

| Outcome | Action |
|---|---|
| R1 PASSes (peak in [100, 200] µm + no deeper lobes + falloff ≤ 25 % of peak) | Ship 16d. Surface BC was the missing piece. Advance to Step 17. |
| R1 FAILs with peak still >> 200 µm | The σ_xx error has been fixed and the remaining gap is mesh-resolution-limited. Accept the mesh-limit conclusion (§9). Document and defer R1; advance to Step 17. **Do NOT iterate.** |
| R1 FAILs with peak in (20, 100) µm or (200, 500) µm (close but not in band) | Surface BC matters but residual physics is missing. Investigate (single round only) — could be `a_f` variation, K_Ic profile, or thermal-diffusion-rate sensitivity. Document findings; ship 16d if a clear physical reason for the offset is found, else defer per §9. **No sub-cell-reconstruction tweaks beyond option (a).** |

Path (b) — node-based interp on `stress_mf` — stays in the design note
as a possible alternative but is NOT the recommended 16d path. (a)
addresses the root cause (surface T jump) more directly and does not
depend on whether the mechanics solve produces a face-resolved σ_xx
under prescribed-T patch eigenstrain.

## 9. Mesh-resolution limit — what is and isn't recoverable

ALAMO needs mesh points (cell centres or nodes) to reconstruct σ_xx at
a given depth. Available z values under the working 16c geometry
(32³ in 50×50×50 mm, dz = 1.5625 mm):

- face: z = 0
- top-cell centre: z = dz/2 = 781 µm
- second-cell centre: z = 3·dz/2 = 2344 µm
- third-cell centre: z = 5·dz/2 = 3906 µm
- ...

Rossi's measured peak at 100–200 µm sits **between** the face (z = 0)
and the top-cell centre (z = 781 µm). A bisection peak in [100, 200]
µm requires K_I = K_Ic to cross somewhere in that range. With **only
two reconstruction points** in the relevant interval (face + top-cell
centre), any sub-cell σ_xx reconstruction strategy — option (L), (T),
(N), or any linear-style interp — has the same fundamental property:
**σ_xx between the two points is parameterized by one shape function**
(a linear ramp, or whatever the chosen scheme produces). K_I integrated
over a 22 µm-long crack window then varies smoothly between the two
end-point values; the K_I-vs-K_Ic crossing lands somewhere in
[0, dz/2] determined by that smooth interpolation.

Whether the crossing lands at 100 µm vs 600 µm depends on the relative
σ_face vs σ_top_cell values. Option (a) uses the prescribed T_face
which gives the correct surface σ_xx; that's the most physically
meaningful sub-cell reconstruction available with two points. **If the
crossing under (a) lands at, say, 350 µm rather than 100–200 µm, no
amount of further sub-cell tweaking will move it to 150 µm** — the
mesh simply doesn't have the resolution to manufacture a peak feature
in that depth range. A Rossi-faithful 100–200 µm peak under this mesh
geometry requires either:

a. **Anisotropic z-refinement near the surface** — `dz_top ≈ 20 µm`,
   giving ~75 reconstruction points across Rossi's 100–500 µm damage
   zone. Reopens MLMG conditioning constraints documented in Step 16's
   takeaway #1; needs careful elastic-solver engineering. Step 18+
   scope.
b. **A non-linear σ_xx model** at the surface that injects physics
   features beyond what the mesh resolves (e.g. an analytical
   thermal-stress correction layer). Out of scope for the 16-series.
c. **A different validation experiment** whose observable lives at
   cell resolution. Tabled in the original design note §1; remains
   tabled.

## 10. Anti-patterns — what NOT to do

Step 16 → 16b → 16c spent three packets on sub-cell reconstruction
refinements. The history is informative; the lesson is short:

- **Don't chase the peak with successive sub-cell reconstruction
  tweaks.** After (a), if R1 still FAILs, the mesh-resolution argument
  in §9 applies; accept it. Sub-cell tweaks beyond (a) tinker with the
  shape of an interpolation that doesn't have enough information to
  match Rossi at this mesh.
- **Don't smooth K_Ic without smoothing σ_xx.** K_I is determined by
  σ_xx integration; K_Ic continuity is second-order. Touching K_Ic
  alone won't move the peak (§7 + §8 above).
- **Don't extend `stop_time` past Rossi's experimental conditions.**
  Past 130 s = 670 °C surface, the simulation is no longer the Rossi
  experiment; the result is no longer comparable. The packet's
  protocol already forbids this; the lesson sticks for future runs.
- **Don't promote R3 (ρ_B/ρ_I) to a binding gate.** It reflects
  late-ramp σ_xx skew toward large-a_f columns, not the mesh-
  resolution issue. Path to 5:1 is orthogonal: higher peak σ_xx,
  re-examined a0_ig, or finer mesh (the §9 (a) path). All out of
  scope for the 16-series.

If Step 16d's (a) gives R1 PASS, ship and move on. If R1 still FAILs,
accept the mesh-limit conclusion, advance to Step 17 (per-cell `v_n`,
revised-approach Stage 4 — re-runs Zhang/Oglesby under `sp_weibull`).
The Step 18 confining-pressure sweep against Kant 2017 is a stronger
validation than further Rossi work on this mesh.

## 11. Final status (post-16d) — DEFERRED

Step 16d ran. Option (a) (T-derived σ_xx with prescribed-T face BC)
was implemented per spec and produced a broad h_col distribution
(range [573, 1556] µm, std 325 µm, peak at 705 µm) — successfully
removed the option (L) flat-extrapolation artefact from §7. But the
peak landed at 705 µm, not in Rossi's [100, 200] µm band.

Per §8 decision tree (binding outcome): R1 FAIL with peak >> 200 µm
→ accept mesh-resolution conclusion → advance to Step 17. Per §10
anti-patterns: do NOT try another reconstruction strategy.

**The Rossi quantitative R1 validation is DEFERRED.** Three packets
(16b, 16c, 16d) tried three different σ_xx reconstructions; all
three peaked outside [100, 200] µm:

| Packet | σ_xx model | Peak | Distribution |
|---|---|---|---|
| 16b | cell-centred event count | 47 µm (cell-quantized) | 1 bin |
| 16c | option (L) flat + linear | 799 µm (artefact at dz/2 − a_f) | 1 bin |
| 16d | option (a) T-derived face BC | 705 µm (broad, std 325 µm) | 10 bins |

The §9 mesh-resolution argument is empirically confirmed: at
dz = 1.5625 mm with two reconstruction points (face + top-cell-
centre) spanning 781 µm, no smooth interpolation scheme can
manufacture a K_I = K_Ic crossing in 100-200 µm.

### Conditions for revisiting Rossi quantitative R1

Path (a) — **anisotropic z-refinement near the surface** (Step 18+
scope). Refine `dz_top ≈ 20 µm` to give ~75 reconstruction points
across Rossi's 100-500 µm damage zone. Reopens the MLMG
conditioning constraints documented in Step 16 takeaway #1; needs
careful elastic-solver engineering. Defer until either
(i) the elastic MLMG handles anisotropic cells reliably (unlikely
on the current solver — Step 15b documented MLMG fails on even 2:1
aspect ratio), OR (ii) a different elastic solver is wired in.

Path (b) — **different validation experiment** with cell-resolution-
appropriate observables. Candidate: Friedrich-Wong slow-rate damage
(`revised-model-approach.md` §6). Its observable is total damage
volume `S_v(T)` from the quadratic damage law `k = E·(Δα)²/[8·(1−ν)·G_Ic]`,
not depth-resolved. Cell resolution becomes non-binding because
ALAMO would report an integral over the cooled sample, not a peak
location. The slow-rate regime is ALSO physically different from
Rossi's gradient-driven flame regime — so it's a different
validation question, not a substitute. May be worth a separate
packet (post-Step 18) regardless of Rossi.

Path (c) — **accept the depth-distribution shape, abandon the band
location**. R1 sub-criteria (2) + (3) PASS under 16d (unimodal,
sharp falloff to baseline). The shape IS Rossi-like; only the
specific peak depth differs. Argue that ALAMO captures the
physical mechanism (early-ramp shallow events → late-ramp deeper
events, GB-biased Weibull distribution) and that the specific
100-200 µm number is a mesh-resolution artifact of where ALAMO
puts cell centres rather than a model prediction. Report this in
the production V&V matrix as "Rossi shape: PASS; Rossi peak depth:
mesh-limited at current resolution". This is a documentation
choice, not a code change.

### Recommended ongoing state of the Rossi test

`tests/MMWSpalling/sp_rossi_damage_profile/`:
- Stays as an ongoing **mechanics regression**. P1-P4 binding
  (h_spall fires, material removes, GB-biased, h_spall bounded).
- R1/R2/R3 are reported but **not gated at the project level**.
  The test prints them; future changes to the diagnostic chain
  (event-count semantics, freeze-guard, etc.) should still see
  consistent values.
- The Rossi peak-location numbers stay informational. Future
  packets touching σ_xx reconstruction (e.g. a hypothetical
  anisotropic-mesh Step 18+) can re-promote R1 to a binding gate
  if the mesh supports it.

### What NOT to do (re-emphasized from §10)

- **Do not** chase the Rossi peak with successive sub-cell
  reconstruction tweaks on the current mesh. 16b → 16c → 16d
  exhausted that path.
- **Do not** smooth K_Ic-vs-z to chase the peak — K_I is
  determined by σ_xx integration; K_Ic continuity is second-order.
- **Do not** extend `stop_time` past 130 s to shift the peak —
  past 130 s the simulation no longer matches Rossi's experimental
  conditions.
- **Do not** modify Weibull parameters (`a0_gb`, `a0_ig`, m) to
  shift R3 toward 5:1 — those are calibrated against Kant; the
  24:1 ratio observed under sp_weibull's late-ramp σ_xx skew is
  a known finding, not a bug.
- **Do not** revisit the patch + free-lateral geometry to seek a
  smaller σ_xx that might shift the peak shallower — this reopens
  Step 16's MLMG conditioning failures.

The `sp_weibull` branch's validation story under the current mesh
ends with: **Kant onset (15b, validated), Kant sensitivity (15c,
validated), Rossi mechanics + shape (16d, validated; depth
location mesh-limited), Step 17 Zhang/Oglesby ROP, Step 18 Kant
confining-pressure sweep**.
