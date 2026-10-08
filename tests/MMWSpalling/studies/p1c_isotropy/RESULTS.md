# P1c — a grid-isotropic support rule

**Verdict: `spall.support_min_sides = 3` FAILS both primary gates.** The named
successor **P1d** applies. No source change, no rebuild: the binary is P1b's
`9e35d750c931f6ab…` and `find src -newer bin/mmwspalling-3d-g++` is empty.

## Scorecard

| # | test | bar | measured | |
|---|---|---|---|---|
| (a) | same binary as P1b | `9e35d750…` | **exact match**; argv diff vs P1b is `plot_file` + `min_sides=2→3` only | PASS |
| (b) | rev 7 healthy, `min_sides=3`, 3 bands × 6 sectors | 0.0 %, 0 pillar-only | **0/316, 0/178, 0/213**, every sector | PASS |
| (c) | P1 stalled floor reach | ≥ 20 %, sector ≤ 2.0 | **24.1 %** (76/316), **1.84**, 6 pillar-only | PASS |
| **(d)** | leg reaches `stop_time` | 500.0 s, rc = 0, no collision | ⚠⚠ **ABORT at 475.184 s**, nozzle collision | **FAIL** |
| **(e)** | stall stays cured | frozen 0; step med ≤ 6, p90 ≤ 8 mm | ⚠⚠ **65 of 316 frozen**; step **4.0 / 26.0 mm** | **FAIL** |
| (f) | roundness, radius per sector at 150/200/250 mm | ≤ 1.05 | **1.475 / 1.566 / 1.708** — ⚠ *reported, not scored* (witness 6) | — |
| (g) | `ledger_err_solid` | ≤ 5e-3 | **7.5973e-04** | PASS |
| (h) | support-firing temperature | none below 700 K | **min 842.2 K**, median 1097.5 K | PASS |
| (i) | the decay | report | **tail/peak 0.000** — vacuous, as pre-registered | report |
| (m) | firings, histogram, pillars | report | **729** total; `{2:68, 3:657, 4:4}`; 68 pillar-only, all legitimate | report |

## ⚠⚠ (d) and (e): the rule is switched off and the stall comes back

```
ABORT: MMWSpalling.H:2231 (BuildFlameColumns) surface_patch nozzle collision:
  stand-off s = -8.87068e-14 <= 0 at column (8,10), time = 475.2
```

r = 27.0 mm — within a cell of P1's own collision at r = 26.4 mm. ⚠ The
stand-off is **−8.9e-14 m**, i.e. zero to round-off, not a deep burial: the
burner touched rock that had not receded, 95 % of the way through the run.

| | P1 | P1b (`min_sides=2`) | **P1c (`min_sides=3`)** |
|---|---|---|---|
| run ends | abort 177.6 s | **500.0 s, rc = 0** | ⚠⚠ **abort 475.2 s** |
| frozen columns of 316 | 152 | **0** | **65** |
| step at P1's cut, median / p90 | 21 / 70 mm | 2.0 / 8.0 mm | **4.0 / 26.0 mm** |
| support firings | — | 22,800 | **729** |
| cells removed (one quadrant) | — | 72,091 | **36,667** |

**729 firings against 22,800 is the whole story: 3-of-4 does not restrict the
rule, it switches it off.** The packet pre-registered the static floor at 111 and
warned that at that reach "the rule is very nearly switched off"; the dynamic
answer is 729, and it is not enough. Removal halves and the hole narrows
drastically — P5's predicted direction, far past its predicted magnitude.

(i) confirms the pre-registration exactly: **tail/peak = 0.000**, because support
firings stop entirely after t ≈ 350 s. **That is not evidence of repair** and was
correctly held out of the gate list.

## ⚠⚠ The finding: P1b's roundness was a CANCELLATION, not a round hole

Per-sector firing counts in the 35–44 mm ring, both legs at their final state:

| | 0° | 15° | 30° | 45° | 60° | 75° | max/min |
|---|---|---|---|---|---|---|---|
| P1b criterion (reg 1) | 3573 | 2178 | 1683 | 1970 | 2244 | 3533 | **2.12** |
| **P1c criterion (reg 1)** | 2121 | 1370 | 1004 | 1330 | 1109 | 2129 | **2.12** |
| P1b support (reg 5) | 1232 | 2140 | 2242 | 2639 | 2038 | 1260 | 2.14 |
| P1c support (reg 5) | 36 | 24 | 13 | 23 | 17 | 34 | 2.77 |

On the floor (r < 25 mm) the criterion reads **1.58 in both legs** as well. The
criterion's per-sector structure is **identical whatever `min_sides` is**, and it
favours the **axes**; P1b's support rule favoured the **diagonals** (11.27× in the
44–50 mm wall). In P1b the two partly **opposed** each other and the net roundness
came out 1.095. With the rule switched off, the criterion's own bias is exposed
and roundness **degrades to 1.475–1.708 with the sense inverted** — axes widest
(41–43 mm), diagonals narrowest (23–29 mm).

⇒ **An isotropic support definition alone will not produce a round hole.** P1d is
necessary but, on this evidence, not sufficient.

⚠ **Caveat, stated because it limits the claim.** Per-sector *counts* are
confounded by the grid: a 15° sector holds a different number of columns by
azimuth, and radial spacing along a diagonal is √2 larger, so some of the 1.58 /
2.12 is a property of the statistic rather than of the physics. The **depth and
radius** statistics are the unconfounded ones, and they carry the shape claim:
the floor (r < 25 mm) is round in both legs (**390 mm in every sector, 1.00**, and
366–370 mm, 1.01), while the outer bands are not.

## ⚠⚠ (j) would have passed vacuously — the packet called this in advance

Depth per sector in the 44–50 mm ring:

| | 0° | 15° | 30° | 45° | 60° | 75° | max/min |
|---|---|---|---|---|---|---|---|
| P1b | 14 | 106 | 255 | 245 | 134 | 16 mm | **18.21** |
| P1c | 10 | 6 | 7 | 6 | 6 | 11 mm | **1.83** |

P1c's ring looks ten times more even — **and it is 6–11 mm deep, i.e. essentially
undrilled.** Read alone, (j) would have reported a near-cure while the hole got
worse. This is precisely the vacuity the packet named when it replaced
depth-in-a-ring with radius-at-fixed-depth, and the replacement is what caught
the real answer. **(f) is why this step has a verdict at all.**

Support share in that ring is **0.0 % in every sector** — confirmation that
`min_sides = 3` never fires on the wall, which was the intent and is also why the
wall stopped advancing.

## Reproduction

`exposure2.py` reproduces both tables in the packet's `## Sources` exactly,
including rev 7's wall at `min_sides = 2` reading **21 / 43 / 89 / 91 / 53 / 21 %**
by sector — the grid-diagonal signature as a measured number. `roundness.py
--selftest` reproduces rev 7 **1.015 / 1.015 / 1.027** and P1b **1.095 / 1.095 /
1.108**, so the statistic is not blind to the defect. `analyze.py` run against
P1b's own log reproduces its t_max, frozen 0/316, ledger 7.5973e-04, min T
738.7 K, tail/peak **0.619**, 22,800 firings and the `{2: 22784, 3: 16}`
histogram.

The `k_top` cross-check in `exposure.py` now prints (deliverable 6): **0
mismatches** on both frozen logs, so the replayed geometry is the one the runs had.

## Carried forward

* **P1d is the successor**, as the packet named it: contiguous void spanning more
  than 180° of the 8-neighbour perimeter. ⚠ But it must be paired with the
  criterion's own azimuthal structure, which this step shows is **independent of
  the support rule** — otherwise P1d buys an isotropic rule on top of an
  anisotropic criterion.
* ⚠ `MMWSpalling.H:5536`'s abort message still tells the user a smooth cavity wall
  "reads exactly one exposed side, inward". P1b disproved that; no source change
  was permitted here, so it stands as a user-visible falsehood.
* ⚠⚠ `a0` and the drive fraction remain **ONE calibration** (refit ×2.46,
  criterion ×2.49). Nothing here changed the profile's sharpness, so none is owed.
* ⚠⚠ No Meier rate, volume, width or standoff is quotable from P1, P1b or P1c.
