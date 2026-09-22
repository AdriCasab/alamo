# D2j-0b results — the hot disk under a descending exclusion plane (2026-09-21)

Reviewer diagnostic (key-only, binary `6d047b50…`, no source edit, no build). Predictions in
`CRITERION.md` (sha256 `2e6595eb…`, locked 15:04, before launch); metrics and scoring in
`analyze.py` → `analysis.md`, `d2j0b_bands.png`; field probes `step_probe.py`. Eight scored runs
(C0/C1/C2 at 2 and 1 mm, C3 and C1d at 2 mm), 2–3 min each at 2 mm, 26 min at 1 mm, all to
250 s, energy ledger closed. An exploratory surface-normal set (`CRITERION_addendum.md`) is
reported in §6 when finished.

## 0. Verdict in one paragraph

**The D2i ring stall is reproduced on a static hot disk with no feet, no jet and no gas closure,
and it is not caused by the exclusion rule, nor curable by heating the exposed wall.** A firing
column beside a taller, colder neighbour loses heat through the step into that neighbour's
interior at a rate ∝ 1/dx; when the neighbour is cold enough the pinned budget cannot pay it and
the column stalls, becoming the next cold wall. The exclusion (`s ≤ 0` → no flux) only makes the
first stall hard and cold early; the same cascade runs without it, three times slower at 2 mm
and unmistakably at 1 mm. **D2j-0's "bounded ≈ 10 % offset" was a 150 s window artefact**: the
plain disk's 1 mm/2 mm gap grows −10 → −15 → −20 % over 50 s blocks and its stalled ring walks
inward after 200 s. Heating the exposed surface (h 100 at 1000 or 1600 K) changes nothing, band
for band, at either mesh: the sink is 10–20 mm below the exposed surface. **This is a structural
property of the per-column surface (flux and losses on the top face only, vertical step faces
adiabatic and immobile), so no key-only configuration converges at a lateral hot/cold edge.**

## 1. Runs

| run | above plane (h, T) | dz | dt | v_c [m/h] | disk-mean 50–150 s | D2j-0 ref |
|---|---|---|---|---|---|---|
| C0_2mm | 700, 1600 (none) | 2 mm | 16 ms | 1.514 | 1.165 | 1.133 (+2.8 %) |
| C1_2mm | 1e-3, 1600 (exclusion) | 2 mm | 16 ms | 1.522 | 1.060 | |
| C2_2mm | 100, 1000 (warm wall) | 2 mm | 16 ms | 1.520 | 1.063 | |
| C3_2mm | 100, 1600 (exhaust-like) | 2 mm | 16 ms | 1.518 | 1.063 | |
| C1d_2mm | as C1 | 2 mm | 8 ms | 1.521 | 1.057 | |
| C0_1mm | none | 1 mm | 8 ms | 1.486 | 1.017 | 1.014 (+0.2 %) |
| C1_1mm | exclusion | 1 mm | 8 ms | 1.459 | 0.843 | |
| C2_1mm | warm wall | 1 mm | 8 ms | 1.459 | 0.843 | |

## 2. Band timing (median over the band's columns)

t_stop = last removal; t_cross = the band's own face rising above the plane; Δ_own = t_cross −
t_stop (> 0: stopped while still under the flux); Δ_outer = t_stop − t_cross of the band outside.

| band | C0 2 mm | C1 2 mm (Δ_own / Δ_outer) | C2 2 mm | C0 1 mm | C1 1 mm (Δ_own / Δ_outer) | C2 1 mm |
|---|---|---|---|---|---|---|
| 28–30 | 110 | 41 (+14 / —) | 41 | 198 | 40 (+13 / —) | 40 |
| 26–28 | 216 | 92 (+10 / +36) | 91 | 239 | 68 (+12 / +15) | 68 |
| 24–26 | firing | 131 (+8 / +29) | 133 | 239 | 91 (+13 / +11) | 91 |
| 22–24 | firing | 170 (+11 / +31) | 170 | firing | 118 (+12 / +13) | 118 |
| 20–22 | firing | 203 (+10 / +22) | 203 | firing | 139 (+14 / +9) | 139 |
| 18–20 | firing | 234 (+7 / +21) | 235 | firing | 162 (+15 / +8) | 161 |

C3 and C1d equal C1 at 2 mm to within 1–2 s in every band. In C1/C2 every band stops 4–6 mm
below the plane; D2i's bands stopped 5–50 mm below it, 33–59 s (2 mm) / 9–23 s (1 mm) after the
outer band crossed — the same pattern and the same mesh ratio as here (21–36 s / 8–15 s).

## 3. Disk-mean recession (r < 30 mm) and the mesh gap

| block | C0 2/1 mm [m/h] | gap | C1 2/1 mm | gap | C2 2/1 mm | gap |
|---|---|---|---|---|---|---|
| 50–100 s | 1.211 / 1.085 | −10.4 % | 1.151 / 0.996 | −13.5 % | 1.154 / 0.996 | −13.7 % |
| 100–150 s | 1.119 / 0.948 | −15.3 % | 0.968 / 0.690 | −28.7 % | 0.972 / 0.689 | −29.1 % |
| 150–200 s | 1.048 / 0.838 | −20.1 % | 0.759 / 0.414 | −45.5 % | 0.762 / 0.412 | −45.9 % |

r_front (inner edge of the stalled run, trailing 25 s) at 100 / 150 / 200 s: C0 none / 28 / 28 mm
at both meshes; C1 and C2 28 / 26 / 22 mm at 2 mm, 26 / 22 / 18 mm at 1 mm.

## 4. Predictions scored (as written in CRITERION.md; no re-windowing)

- **P1 — fails as written, sanity intent met.** v_c is within 5 % of 1.58 m/h at 2 mm
  (1.51–1.52) but 6–8 % low at 1 mm (1.459–1.486). The C0 disk-mean 50–150 s matches D2j-0 within
  0.2 / 2.8 %, so the configuration reproduces D2j-0; the 1 mm v_c shortfall is the cascade
  reaching r < 18 mm by 162 s and starving the centre columns laterally late in the 50–250 s
  window. This is a post-hoc reading and is flagged as such; it does not rescue the threshold.
- **P2 — C1 half holds, C0 half fails.** C1 cascades at both meshes (all five bands inside 28 mm
  stop with Δ_own > 10 s at 1 mm; three of five at 2 mm). **C0 also cascades**: 26–28 stops at
  216 s (2 mm); 26–28 and 24–26 at 239 s (1 mm). This was predicted in the daily notes at 15:26,
  after the 2 mm set and before the 1 mm data, on the grounds given in §5.
- **P3 — C1 half holds, C0 half fails.** C1's gap grows monotonically to −45.5 % (D2i: −36 %). C0's
  gap is not standing: −10.4 → −15.3 → −20.1 %.
- **P4 — fails.** C2 equals C1 in every band, block and front position at both meshes.
- **dt — holds.** C1d vs C1: t_stop within 0.5 s, block rates within 0.4 %.

By the decision rule: P2 fails in its C0 half, but not in the direction the rule anticipated
(the exclusion is *sufficient* and *not necessary*), and P4 fails, so **neither the feedback packet
nor the wall-heating packet follows.** What follows is §5.

## 5. The mechanism, from the fields

`step_probe.py` (row y = dz/2) at 150 s (2 mm) and 100 s (1 mm):

- Firing top cells sit at 700–815 K; their outward neighbour column is taller by a step of
  4–20 mm; the neighbour's cell *at the same k* (the wall cell the top cell conducts into) is at
  480–560 K in C0, C2, C3 alike and 400–500 K in C1. The exposed surfaces of the stalled columns
  differ enormously between variants — 804–809 K in C0 (held by the imposed pinned flux; D2i's
  "hot silent pad columns"), 430–500 K in C2/C3, 377–407 K in C1 — while the *wall cells* differ
  by 5–60 K, because they lie 10–20 mm below the exposed surface, beyond what 50 kW/m² of surface
  heating can reach in 30 s. That is why C1 = C2 = C3 and why C0 is merely slower.
- The loss through the step is k(T_top − T_wall)/dx per unit top area (top cell to one lateral
  cell): 1.5·(821 − 500)/dx = 240 kW/m² at 2 mm, 480 at 1 mm, against a pinned budget of
  q_pin ≈ 520 kW/m² minus downward conduction. The threshold T_wall at which the firing column
  starves therefore rises with 1/dx; at 1 mm even the warm interior of an un-excluded stalled ring
  (560–650 K) is enough, which is why C0 cascades at 1 mm and only slowly at 2 mm. An
  instantaneous per-column budget test (`wall_probe.py`) does not separate firing from stalled
  columns because the top-cell temperature cycles through 500–821 K between the two firings per
  cell; the band-timing evidence carries the conclusion, not the budget arithmetic.
- The inner side of a firing column is void at its k_top (the inner column is deeper), so it has
  exactly one cold lateral sink, the step. The exposed vertical face of the step receives no flux
  and no loss and cannot recede. In the continuum the re-entrant corner would be the hottest
  stress concentration and would spall away; in the per-column model it is a permanent,
  unheated, dx-sharp sink.

## 6. Exploratory: `spall.surface_normal = 1` (C1n / C0n)

Prediction (addendum, written before these ran): the kernel changes how deep a column cuts per
event, not the step flux, so the cascade stays. **Result at 2 mm: inert.** The key was applied
(command line logged; the removal logs differ from line 17 on, so the kernel ran), yet C1n equals
C1 in every band's t_stop (41 / 92 / 131 / 170 / 203 / 234 s), in r_front and in the block rates
to 0.3 %; C0n equals C0 likewise (26–28 stops at 217 vs 216 s). D2b's "inert under the pinned
closure" holds for the cascade too: the Step 20 kernel acts on the criterion and the cut depth,
not on the step flux, and is not a candidate remedy. **C1n_1mm (15 min): inert as well** — t_stop 40 / 68 / 91 / 118 / 139 / 162 s, r_front 26 / 22 /
18 mm, block rates 0.996 / 0.689 / 0.411 m/h, all equal to C1_1mm to the last digit shown.

## 7. What this means for the programme

1. **No key-only fix exists** for the D2i mesh failure. Every configuration tested since D2d
   (clip off, idle off, rec_fixed, hot disk, plane, warm wall) shares the structural cause.
2. **The D2j-0 conclusion "the edge is bounded, do not redesign the front" is withdrawn.** The
   150 s window hid the growth. The redesign question is open again, and the D2h open-loop pair
   would only re-measure the same cascade with the feet removed.
3. **The remedy must act on the step**: either the exposed vertical faces of the column surface
   become part of the heated, receding boundary (a lateral flux on exposed side faces plus
   lateral recession — the surface-normal idea done on the flux side, not only the criterion), or
   the front is represented continuously (the isotherm recession of `2026-09-21b.md`, whose slope
   term is now motivated by exactly this corner). A cheap arbiter for any candidate is this
   study: C1 at 2 and 1 mm, band table and gap — a converged fix gives Δ_outer independent of dx
   and a standing gap.
4. Until then the Meier statement of `2026-09-21a.md` (window 150–350 s) stands as written; its
   limitation 1 should be reworded from "thermal knife-edge at the pad rim" to the mechanism above.
