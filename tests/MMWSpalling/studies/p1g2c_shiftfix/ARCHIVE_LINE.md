# P1g-2c `(z)` — the archive line

Deliverable 8. Lift verbatim into `docs/project/ARCHIVE_DONE.md`. Numbers from
`LEAK_SIZE.txt`, `FIXTURE_pre/postfix.log`, `ANCHORS_fixed.txt`, `SCORE_p1g2c.txt`.

---

**P1g-2c — the thermal shift's energy leak, closed (2026-10-06).** P1g-2b's `(b3)`
failed on the thermal leg and was correctly not waived. The cause was one
character of scope: at `Removal.H:3407` the row write-back was an **assignment**,
`Hw(ii,jj,kf) = rNew[base+m]`. A face exposed on two sides runs two rows under the
cos² split; the rows share the face cell at `m = 0`, and **cross-face** rows
collide at `m > 0` as well — a face at `(i,j)` shifting +x and one at `(i+2,j−2)`
shifting +y both write `(i+2,j)`. Both shares were booked into `E_shift_local`,
but the second assignment **discarded the first**, so the field lost only one of
them.

**Sized before it was touched.** The unclosed energy was **7,616 J**, **0.89 % of
the budget** and **15.3 %** of the shift's booked **49.8 kJ**. ⚠ `ledger_err_solid`
is **dimensionless** (`MMWSpalling.H:4587` divides by `max(|E_src|,1e-30)`), so
P1g-2b's "absolute 5e-3 J bar" misnamed it and its mitigation — that relocating
more enthalpy inflated the bar — was void, because the bar is already normalised
by the cumulative source. In joules the residual was **7,616 J against the three
clean legs' 45 J**.

**Why two reviews missed it: both existing ratios were tautologies.**
`lat_shift_ratio` read **1.0000** at p05/p50/p95 and `lat_shift_E_est/lat_shift_E`
read **0.947** — and both are built from the **intended** shift (`vol_disp`,
`E_shift_local`), never from the field after the write. The new diagnostic
`lat_shift_real` sums `H` over the touched cells **before and after** the
write-back, which is the only form that can see an overwrite.

**⚠⚠ And the two-axis path had never executed in any test.** The fixture *has* a
45° staircase (`diag_seed`), but only `run_sweep` used it and that hard-sets
`spall.lateral_removal=0` — no removal, no banking, no shift. Every leg that ran
the shift used an axis-aligned riser where *every open face is +x*. Worse,
`spall.lateral_thermal` appeared **zero times** in the committed test, so as
shipped it ran **no shift at all**. Two legs were added: `G1_FLAT_SHIFT`
(single-touch control) and `G1_DIAG45` (two rows sharing cells). ⚠ Both need
**n = 20**, not the fixture's 8: `ROWMIN = 6` lets a corner shift only while
`5 ≤ j ≤ n−6`, which is **empty at n = 8**, and `axis_seed(n)` leaves only `n/2`
cells behind the wall — four at n = 8, so the flat control shifted **zero** times
on the first run.

**The fix, and why it is not a blanket `+=`.** Per-cell deltas are accumulated
with a touch count, then applied: a cell touched **once** is **assigned** exactly
as before, a cell touched more gets `rH + Σ(rNew_ax − rH)`. A blanket increment is
algebraically identity-preserving on a flat wall but `old + (new − old) != new` in
floating point, so it would have drifted every frozen anchor at round-off and the
identity gate would have had to be loosened — which is how a real change hides.
**No key gates the fix**: a conservation bug is corrected, not made optional.

**Results.** On the fixture the 45° leg goes **0.669 → 1.000000** realised/booked
with `ledger_err_solid` at **2.4e-15**; the flat control holds **1.000000** and its
realised energy is **identical** before and after. **68 of 70** anchor files are
**byte-identical** across the fix — the only two that move are `G1_DIAG45`'s,
exactly where the fix acts. Identity **0 of 2639** with the lateral keys off. At
Meier scale the leak closes: **8.94e-3 → 7.78e-4**, matching the three non-thermal
legs to **5 significant figures**, unclosed **7,616 J → 54 J**, and
**realised/booked = 1.000000**. Measured two-axis coverage on the Meier wall is
**0.891** (P2 inferred ~0.80).

**⚠⚠ P3 IS REFUTED, AND THAT IS THE STEP'S REAL RESULT.** The honest thermal-half
figure is **−3.18 mm**, not the −4 to −5 mm predicted: **the fix's own bore effect
is −0.16 mm.** So a 15.3 % energy leak was worth **0.16 mm of bore**. P1g-2b's
−3.02 mm happened to sit within 0.16 mm of the truth — but it was still right to
refuse to quote it, because nothing short of this step could have known that. ⇒
**the energy books and the bore are only weakly coupled here**, which is a
general caution for this programme: a conservation failure is worth fixing on its
own terms, and its size is not a guide to how much geometry it was distorting.

**⚠ One prediction partially met.** P1's windowed bars wanted `≤1.5e-4` for
t > 50 s and `≤1.0e-4` for t > 150 s. The second is met (**9.87e-5**); the first
**misses by 2.5 %** (**1.54e-4**). Reported as a near-miss on a prediction, not a
gate failure — `(d1)`'s bar is 5e-3 and the write-back is now exact. ⚠ The two
comparison legs are **different trajectories**, not repeats: `CTRL4_FIXED` does
0.84× `CTRL2`'s lateral firings, so its residual is not required to equal
`CTRL2`'s, and a few per cent between runs that remove different amounts of rock
is not evidence of a leak.

**Also corrected: the refusal counter was mislabelled, not broken.**
`shifts + refusals == banks` **exactly**, and `refusals/voids = 1.01`. Every bank
attempts a shift; on a step where the bank **voided** its own face cell the row
gathers length 0 and the attempt landed in the `runlen < ROWMIN` branch, sharing
one counter named as a failure. Split three ways (`lat_ref_ell`, `lat_ref_void`,
`lat_ref_short`), the fixture shows it is **all void** on the flat leg (140/140)
and **96 of 122** on the 45° leg. ⇒ "the shift acts on only a quarter of events"
is **wrong**.

**⇒ P2-1 is the floor taper at 40–56 mm.** The depth-matched bore error stands at
**+16.2 mm** with the thermal half on and honest, so the remaining 13 mm is the
floor, not the wall.

⚠ `(c1b)` is **unscoreable**: the `evidence_p1g2` legs have **no recorded argv** —
`PROVENANCE.txt` carries binaries and compilation stamps only, and their launcher
was a session scratch script. `(c1)` was scored on a strictly larger replayable
set instead (70 files, 36 legs). Same provenance gap as Leg B's missing binary in
P1g-2b; the lesson is to record the **argv**, not only the binary hash.
