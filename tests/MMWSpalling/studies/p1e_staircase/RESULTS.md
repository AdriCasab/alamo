# P1e — the staircase side-face correction

**Verdict: the correction is RIGHT as a discretisation fix and it largely cures
the stall, but it FAILS both of its gates — (f) got WORSE, so the packet's
attribution of the azimuthal asymmetry to side-face area is REFUTED by its own
primary test.**

## Scorecard

| # | test | bar | measured | |
|---|---|---|---|---|
| (a) | key OFF vs `p1b.md5` | 0 of 2639 | **0** (and 0 vs `p1.md5`) | PASS |
| (b) | estimator vs analytic, rev 7 wall | within 0.05 | **0.025** | PASS |
| (c) | wall and pillar exactly 1 | exact | exact, plus the range invariant `[1/√2, 1]` | PASS |
| (d) | unit test, 6 legs | 1 / 1 / 0.707107 / identical | **19/19 checks** | PASS |
| (e) | `ledger_err_solid` | ≤ 5e-3 | **7.6055e-04** | PASS |
| **(f)** | per-column firing asymmetry, 35–44 mm | ≤ 1.25 | ⚠⚠ **2.04** (P1c 1.84, rev 7 1.01) | **FAIL** |
| **(g)** | roundness at 150/200/250 mm | ≤ 1.10 | ⚠ **1.257 / 1.393 / 1.619** (P1c 1.475/1.566/1.708) | **FAIL** |
| (h) | `side_sc_mean` | expect 0.78–0.84 | **0.860** (min 0.776) — ⚠ above the band | report |
| (i) | `side_P_in35` per 100 s | — | **48 / 37 / 52 / 67 / 80 %** | report |
| (j) | axis depth | −3 to −10 % vs P1b | **388 mm** (P1b 390, P1c 374) = **−0.5 %** | report |
| (k) | cells removed, one quadrant | −8 to −12 % vs P1b | **45,673** (P1b 72,091 = **−37 %**; P1c 36,667) | report |
| (l) | frozen at P1's cut | — | **3 of 316** (P1b 0, P1c 65) | report |
| (m) | axis columns frozen, 35–40 mm | pre-registered 10–12 % | **0 of 28 = 0.0 %** | report |
| (n) | D2z floor-step populations | re-read | `delta_mean` 5.891, `step_g_mean` 0.973, `riser_area_ratio` 0.461, `wall_step_f_multi` 0.829 | report |
| (o) | 35–44 mm firings at `h_f` = 1.5 mm | — | **61.9 %** (P1b 64.4, P1c 43.7) | report |

⚠⚠ **The control is P1c, not P1b.** P1c is `min_sides = 3` with the key OFF, so it
differs by exactly ONE key. P1b differs by two and cannot separate the correction
from the support rule. Every delta above names its reference.

## ⚠⚠ (f) FAILED IN THE WRONG DIRECTION — the attribution is refuted

Criterion firings per column (grid census), 35–44 mm ring, all t:

| | 0–15° | 15–30° | 30–45° | 45–60° | 60–75° | 75–90° | max/min |
|---|---|---|---|---|---|---|---|
| rev 7 (control, no rule) | 179.4 | 180.4 | 181.2 | 181.7 | 180.1 | 180.4 | **1.01** |
| P1b (`sides = 2`) | 82.0 | 42.8 | 33.1 | 31.9 | 43.0 | 79.6 | 1.81 |
| P1c (`sides = 3`, key off) | 84.8 | 57.1 | 47.8 | 53.2 | 46.2 | 85.2 | 1.84 |
| **P1e (key ON)** | **123.2** | 77.0 | **70.7** | 78.8 | 73.8 | **144.5** | ⚠⚠ **2.04** |

**The packet's mechanism predicted this would fall toward 1.01. It rose.**

The packet's thesis was that a round wall's diagonal columns receive two faces'
heat where an axis column receives one, so the diagonals run *hotter*, their
profiles flatten, and the gradient-aware criterion starves them. The correction
removes area in proportion to that over-count — **0.72 at the diagonals against
0.92 at the axes** — so it cools the diagonals *more* than the axes. The diagonals
then fired **even less relative to the axes**.

⇒ The local relationship is the ordinary one, **more heat means more firing**, so
the diagonal firing deficit cannot be caused by diagonal over-heating. Something
else makes the diagonals fire less, and the area over-count was **partially
offsetting** it. Witness 6 names where to look: the `h_expr` map and the jet
binning are already proven azimuthally symmetric, so **the reconstruction stencil
is the next suspect, not the heat map.**

⚠ The correction is still *physically right*: it is exact on every canonical case,
the ledger closes to 7.6e-4, and the laboratory is byte-identical. It removes real
over-delivery — **4/π − 1 = 27.32 %** on average and **√2 − 1 = 41.42 %** of
diagonal-to-axis contrast. It is simply **not the cause of the azimuthal firing
asymmetry**, which is what (f) was built to decide.

## What the correction DID fix: the stall

| | P1c (control) | **P1e** |
|---|---|---|
| run ends | **abort 475.2 s**, collision at r = 27.0 mm | **500.0 s, rc = 0** |
| frozen at P1's cut | 65 of 316 | **3 of 316** |
| axis columns frozen, 35–40 mm | **39.3 %** (11 of 28) | **0.0 %** (0 of 28) |
| axis depth | 374 mm | **388 mm** |
| cells removed (one quadrant) | 36,667 | **45,673** |

⚠⚠ **P5 FAILED in the favourable direction.** It pre-registered that 10–12 % of
axis columns at 35–40 mm would still freeze, with P1c measured at 10.6 %. **I
measure P1c at 39.3 %, not 10.6 %**, and P1e at **0.0 %**. So the starvation P5
named as residual is *gone*, and the prediction's own baseline does not reproduce.

⇒ Removing the over-delivery **did** fix the freezing — which is consistent with
the over-delivery being real — while leaving the firing asymmetry and the shape
anisotropy untouched. **Those are two different problems, and this step separates
them.**

## (g) roundness: improved at every depth, nowhere near the bar

| depth | P1c | **P1e** | rev 7 |
|---|---|---|---|
| 150 mm | 1.475 | **1.257** | 1.015 |
| 200 mm | 1.566 | **1.393** | 1.015 |
| 250 mm | 1.708 | **1.619** | 1.027 |

⚠ The hole is **lobed along the axes**: at 200 mm the axis sectors reach 43.0 mm
while the diagonals reach 30.9–33.6. That is the same sense as P1c and the
*opposite* of P1b, whose rule over-drilled the diagonals. Consistent with (f).

## The rate predictions, both missed in magnitude

P4 predicted axis depth **−3 to −10 %** and volume **−8 to −12 %** against P1b.
Measured: depth **−0.5 %** (388 vs 390 mm) and volume **−37 %** (45,673 vs
72,091). ⚠ But the volume figure is **not attributable to the correction**: P1b
ran `min_sides = 2` and its support rule alone removed 22,800 cells' worth of
extra firings. Against the one-key control P1c the volume **rose** (36,667 →
45,673, **+25 %**), because the hole no longer stalls. **Only the P1c comparison
is a statement about the correction.**

`side_P_in35` rises **48 → 80 %** across the run: side heat starts on the wall and
migrates onto the floor as the floor roughens, which is why a single "the
correction removes ~10 % of total heat, so the rate falls" prediction could not
have been right in either direction.

## Reproduction and the reference values

`factor.py` reproduces the packet's wall table to the digit — **0.920 / 0.793 /
0.720 / 0.719 / 0.789 / 0.913**, mean 0.809, max deviation 0.025 — and the floor
figures exactly: **mean 0.867, median 1.000**, histogram 1-face 225 / 2-face 180 /
3-face 15 of 420 (54 % single-face).

⚠⚠ **Three packet reference values did NOT reproduce, and one of them changed the
gate's definition:**

1. **(f)'s values only reproduce over ALL t.** rev 7 1.01, P1b 1.81, P1c 1.84 come
   out exactly with no time cut. Under a `t > 100 s` cut P1b reads **2.57** and
   P1c reads **infinity — a dead sector** with zero firings at 30–45°. The 1.25
   bar belongs to the all-t statistic, and that is what is scored above.
2. **(m)'s P1c baseline is 39.3 %, not 10.6 %.**
3. **`side_sc_mean` lands at 0.860, above P3's 0.78–0.84** — consistent with the
   wall mean 0.809 only if the floor (median 1.000) is included, which it is.

## Carried forward

* ⚠⚠ **The azimuthal firing asymmetry is NOT the side-face area count.** The
  successor is the **reconstruction stencil** (witness 6), not the heat map and not
  another support-rule round. The flux estimator — a three-parameter fit with a
  substrate base, replacing two fits that both assume decay to ambient and recover
  only **0.38** of the true flux at 30–50 mm — remains named, but (f) now says the
  *stencil's* azimuthal behaviour must be checked first.
* ⚠ `side_sc_mean` 0.860 means the shipped correction removes ~14 % of side-face
  area, i.e. ~7 % of total delivered heat (side faces carry ~half).
* ⚠⚠ **D2z's floor-step conclusions must be re-derived**, not assumed:
  `riser_area_ratio` 0.461 and `wall_step_f_multi` 0.829 are measured with the key
  on. The factor falls on the 46 % of floor steps with two or more open faces; the
  54 % single-face majority reads exactly 1, which is why the floor median is 1.000.
* ⚠⚠ `a0` and the drive fraction remain **ONE calibration** (refit ×2.46,
  criterion ×2.49). This step does not change how the profile reads, so none is owed.
* ⚠⚠ **No Meier rate, volume, width or standoff is quotable** from P1, P1b, P1c or
  P1e. The bore is reported per sector only — never as a single diameter, which is
  what produced the retracted "the hole is cylindrical" claim.
