# P1b -- unsupported hot rock comes off

Revision 2, implemented. One binary for every scored run; nothing rebuilt after
the first scored run of the final binary. `BINARY.sha256` pins it,
`MANIFEST.md` has a row for every leg.

## Scorecard

| | test | bar | measured | |
|---|---|---|---|---|
| **(a)** | `support_rule = 0` vs `p1.md5` | 0 differing / 2639 | **0** (and 0 vs `p0b7.md5`) | PASS |
| **(a)** | parse guards | abort outside 1-4 | aborts; also refuses the absolute form | PASS |
| **(b)** | unit test, 4 legs | tower fires, wall spared, flat identical, slab zero | all four | PASS |
| **(c)** | ⚠⚠ the rule is self-extinguishing | tail/peak < 0.20 | **0.657** — cause is GRID ANISOTROPY | **FAIL** |
| **(d)** | the run completes | 500 s, rc = 0, no collision | **500.0 s, rc = 0** | PASS |
| **(e)** | frozen columns at P1's cut | 152 -> < 30 | **0** | PASS |
| **(f)** | adjacent-column depth step | median 21 -> < 6 mm | **2.0 mm** (p90 8.0) ⚠ the median hides an 18x azimuthal spread | PASS, but see (c) |
| **(g)** | two-axis sweep | report all | `EXPOSURE_SWEEP.txt` | done |
| **(h)** | `ledger_err_solid` | <= 5e-3 on the column; report Meier | **7.60e-4** on Meier after fixing a defect this criterion caught | PASS |
| **(i)** | manifest row per scored leg | every leg | `MANIFEST.md` | PASS |
| **(j)** | bore | report only | see below | reported |
| **(k)** | ⚠⚠ inert on a HEALTHY floor | < 1 % vs > 30 % | **0.0 %** vs **35.4-40.8 %** | PASS |

## ⚠⚠ (c) FAILS — AND THE CAUSE IS GRID ANISOTROPY, NOT THE WALL CORNER

Support firings per 25 s rise to a peak of **1716** at 100-125 s and then settle
to a **plateau near 1000-1100**: tail/peak = **0.657** against a 0.20 bar.
**It was not tuned** — no setting was changed after seeing this number.

⚠⚠ **CORRECTION. An earlier version of this document attributed the plateau to
the 25-50 mm "wall corner" and to D2ad's radial fall-off. That was wrong, and it
was a grid artefact misread as physics.** The real cause is that the exposure
test is **not grid-isotropic**:

> On a Cartesian grid a **round** inward-descending step presents **two**
> inward-lower face neighbours along its **diagonals** (both `i-1` and `j-1` sit
> at smaller radius, hence deeper) and only **one** along its **axes** (there,
> the two tangential neighbours are at the same radius and the same depth).

So `support_min_sides = 2` fires on the diagonals of a smooth round wall and
spares its axes. The packet's own justification — *"a smooth hole wall has
exactly 1 exposed side (inward), so it is spared"* — holds **only for a
grid-aligned wall**, and I repeated it without testing it.

**Measured, final depth of the 44-50 mm ring by azimuth:**

| | 0-15° | 15-30° | 30-45° | 45-60° | 60-75° | 75-90° | spread |
|---|---|---|---|---|---|---|---|
| **P1b** | **14 mm** | 106 | **255** | 245 | 134 | **16 mm** | **18.2x** |
| rev 7 (healthy) | 298 | 310 | 313 | 313 | 310 | 300 | **1.1x** |

rev 7's ring is round to ±2.5 %. P1b's is gouged to 255 mm on the diagonals
while the axes sit at 14 mm. Support firings in 25-50 mm run **3.1:1
diagonal-to-axis** (5737 in 45-60° against 1727 in 0-15°).

⇒ **The rule is a removal law with a square-grid signature.** It re-drilled the
wall band the P1 criterion was meant to suppress, and it did so in a pattern the
mesh chose. The (c) plateau is this gouging sustaining itself, not floor
maintenance and not the D2ad corner.

⚠⚠ **Consequences for the other criteria.** (e) and (f) still pass and the stall
is genuinely cured — but **(f) passed on a MEDIAN over r < 40 mm, which hides an
18x azimuthal spread just outside it.** A median is the wrong statistic for a
defect with angular structure. The witness this step needed and did not have is
an **azimuthal-isotropy check**: the 44-50 mm ring median depth per 15° sector,
within ±3 %. rev 7 passes it; P1b fails it by 18x.

⚠ **Nothing about the hole's width or shape is quotable from this run**, and the
support share rise (4.8 -> 25.4 %) is partly a denominator effect: criterion
firings fall ~4x over the run as the standoff opens, so a flat support count
raises the share on its own.

### The fix is already keyed, and was never run on Meier

`support_min_sides = 3` is grid-isotropic for this geometry — three of four face
neighbours cannot be lower on a round wall — and it clears (k) at 0.0 % on rev 7
(`EXPOSURE_SWEEP.txt`). Its reach on P1's stalled floor is ~22 % against ~40 %,
so it may not clear the stall; **its Meier leg was never run, only the offline
exposure test.** If 3-of-4 misses the stall, the principled form is a contiguous
void span over **more than 180° of the 8-neighbour perimeter** (>= 5 of 8, or
opposing faces), which is isotropic by construction rather than by luck.

## (e) and (f): the stall is gone

| | P1 | **P1b** |
|---|---|---|
| frozen columns (last fired < 120 s, r < 40 mm, t <= 177.6 s) | 152 of 316 | **0 of 316** |
| adjacent-column depth step, median | 21 mm | **2.0 mm** |
| ... p90 | 70 mm | **8.0 mm** |
| run ends | **ABORT at 177.616 s**, nozzle collision at r = 26.4 mm | **500.0 s, rc = 0** |

P2 is confirmed: the colliding column had 4 of 4 sides exposed and the rule
removes that population before it can strand.

⚠ P5 predicted the rule would miss the 22 % of frozen columns with exactly one
exposed side, and that (e) might land between 30 and 60. It landed at **0**.
The reason is visible in the mechanism: the rule does not have to catch a
laggard directly. Removing its proud neighbours lowers the step the laggard
sits in, and the laggard's own profile re-steepens, so the ORDINARY criterion
takes it. The rule levels the floor; the criterion does the rest.

## ⚠⚠ (h): a real defect, found by the criterion, fixed, and re-run

The first complete 500 s run read `ledger_err_solid` = **-0.258** against P1's
**7.6e-4** -- a 300x regression of an invariant the P0b line spent a whole
revision closing.

**The cause, measured not guessed:** `ledger_E_flake` -- what left the SOLID set
with the rock -- was accumulated for regime 1 only. Every regime-5 void took a
whole cell's enthalpy out of the books unrecorded. The arithmetic pinned it
exactly: **-8.95 J per voided cell** against the **8.69 J** one 2 mm cell holds
at 500 K of superheat, over 3167 firings, every one of them `n_voided == 1`.

⚠ The error first departs from zero at t = 19 s and the first support firing is
at t = 16.7 s. The correlation and the per-cell arithmetic together are what
make this a diagnosis rather than a suspicion.

Fixed by accumulating the voided cell's own enthalpy (`E_sup_local`), which is
the sink the packet specifies. The run was redone on the fixed binary.

**After the fix, at the same cut (t <= 100 s): 1.39e-1 -> 7.60e-4**, the same
magnitude as P1's own 7.60e-4, with the running value at 9.39e-5. So the Meier
leg now CLEARS the 5e-3 bar, rather than only being reported.

⚠⚠ **And the fix is provably bookkeeping, not behaviour.** `output_fh_ledgerbug/`
is kept with its event-log sha256, and the fixed binary reproduces it **row for
row: 71,119 rows compared, 0 differing** -- across a full rebuild on a toolchain
that is NOT byte-reproducible. Every (c)/(d)/(e)/(f) number in this document
therefore stands unchanged by the fix. All eight unit-test legs also reproduce
byte for byte.

## (g) the two-axis sweep -- a sensitivity statement, not a search

Fraction of columns at r < 40 mm that the rule would qualify, offline, on both
geometries (`EXPOSURE_SWEEP.txt`, full table):

| setting | rev 7 (HEALTHY) | P1 (stalled) | (k) |
|---|---|---|---|
| sides >= 2, 2 mm | 30.7-32.9 % | 44.0-47.5 % | FAIL |
| sides >= 2, 4 mm | 12.0-12.7 % | 38.9-43.4 % | FAIL |
| **sides >= 2, 6 mm (shipped)** | **0.0 %** | **35.4-40.8 %** | **PASS** |
| sides >= 2, 8 mm | 0.0 % | 34.2-40.2 % | PASS |
| sides >= 3, 4 mm | 0.0 % | 22.2-22.8 % | PASS |
| sides >= 3, 6 mm | 0.0 % | 21.2-22.2 % | PASS |

⚠⚠ The sweep says more than "6 mm is a safe choice". On the HEALTHY floor the
qualifying fraction collapses 32.9 -> 12.0 -> **0.0 %** as the threshold goes
2 -> 4 -> 6 mm; on the STALLED floor it barely moves, 47.5 -> 43.4 -> 40.8 %.
**The two floors are not rough in the same way.** A healthy floor's roughness is
a one-to-two-cell phase lag; a stalled floor's is a four-plus-cell chasm. 6 mm
is the gap BETWEEN two distinct populations, not a cut through one -- which is
why the threshold is not a tuning dial and why (k) is the criterion that
matters.

⚠ `support_min_sides = 3` also clears (k) and would even allow 4 mm, but reaches
only ~22 % of P1's stalled columns against ~40 %. The shipped (2, 6 mm) clears
(k) with the widest reach.

⚠ One number does not reproduce: the packet's rev 7 row at 4 mm reads 8.5-9.8 %;
I measure **12.0-12.7 %**. The shipped 6 mm row reproduces (0.0 against the
packet's 0.0/0.3/0.0) and every P1 row reproduces within a point, so the
disagreement is confined to the one setting the packet already excludes -- but
4 mm fails (k) by MORE than the packet expected, not less.

## The witnesses

| | witness | result |
|---|---|---|
| 1 | `n_exposed` read from `k_top`, not a recession sum | from `k_top_global`, the ReduceIntMax'd array; every rank holds the same one |
| 2 | no decomposition dependence | the input is already reduced, so the count is identical on every rank **by construction**; (a) is 0/2639 at the tree's layouts |
| 3 | `n_voided == 1` on every support firing | **3167 of 3167** on the Meier leg; the assertion is extended to regime 5 and stayed green |
| 4 | no new temperature key | none. `1 - exp(-spall.flake_c_f)` is COMPUTED, never written as a literal |
| 5 | not silently inert | 19 support firings by t = 18.5 s, in the first 25 s window; (e) moved 152 -> 0 |
| 6 | the rule must stop | **IT DOES NOT.** (c) fails. See above |
| 7 | threshold left at one cell | it is not: 6 mm, and (k) run BEFORE the Meier leg, as the packet requires |
| 8 | firing temperature not 525 K | unit test brackets it at **745-750 K**; on Meier the median `T_surf` at a support firing is **861.4 K** (p10 753.3) |

⚠⚠ Witness 8 is what separates P6's expected rate rise from defect 8, and it is
clean: no support firing happened below ~745 K anywhere.

## (j) the bore, REPORTED, NOT SCORED

P4 pre-registered that the bore cannot be judged from this run and would widen
and overshoot before settling, because the 40-50 mm band is rugged NOW. The
firing table above confirms the premise directly: that band carries the largest
share of support firings at every time. **No bore number from this run is
quotable**, and none is quoted here.

⚠⚠ No Meier rate, volume, width or standoff from P1 is quotable, and the collar
(Ø 93.8 / 88.2 mm) must never appear without the cone (72 mm at 100 mm, 56 at
150) and the stall.

## Carried forward, unresolved

* ⚠⚠ **(c) failed because the exposure test is not grid-isotropic** — it drills
  the diagonals of a round wall (44-50 mm ring: 14 mm on the axes, 255 mm at
  30-60°, against rev 7's ±2.5 %). The successor is an isotropic support
  definition plus an **azimuthal-isotropy witness**, not the wall corner.
* ⚠⚠ **`a0` and the drive fraction are ONE calibration.** Unchanged here, so not
  owed here, but it must not be lost: the flaw enlargement alone multiplies the
  column rate by 2.46x and the differential criterion divides it back by 2.49x.
* `K_Ic`'s tip temperature still reads `T_at_zv`; `a0_gb` = 400 um still needs a
  grain-boundary justification; `ell >= dz/4` is still a mesh length setting a
  physical one; per-grain flaw draws; 7b (floor risers) still untested.
* ⚠ The estimator recovers 0.92 of the floor flux and **0.38 at 30-50 mm**. That
  is the same band, but after the anisotropy diagnosis it is NOT evidence the two
  are one problem -- the band was selected by the mesh here, not by the physics.
