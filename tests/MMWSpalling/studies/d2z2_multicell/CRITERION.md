# D2z-2 CRITERION — hashed before the first leg, zero outputs on disk

Repair the three measurement defects D2z left, fix `StepConductance`'s
normaliser, and settle the multi-cell branch on Meier.

```
  g(delta, n, a) = (1 - e^{-delta a}) / (1 - e^{-n a}) ,  a = dz*v/alpha
```

No free parameter, no clamp, no key: the grid has already resolved the step as
`n` cells, so its own drop across that face goes as `1 - e^{-n a}`. D2z
normalised by `1 - e^{-a}` regardless of `n`. **`n = 1` is exact**, so nothing
D2z measured inside r < 25 mm (0.0 % multi-cell) moves.

Metric: `col_rate` and the metadata reader are **imported unchanged** from
`studies/d2j0b_plane/analyze.py`; the equal-burner-depth machinery (`shape`,
`rop`, `rings`, the volume bases) **unchanged** from
`studies/d2x_areasplit/analyze.py`. Frozen studies are read by path, never
re-run: `W2`, `W1`, `RB_D`, `S1`, `S1w`, `S2`, `S2w`, `S2f`, `A`, `D`, `M1`.

⚠ **TWO BINARIES.** The conductance is a source change, not a key.
`R2`/`R2w`/`SW0` ran the **D2z form** (sha `f9085c65…`); `R2d`/`R1d`/`SW`/`SWs`/
`S0` ran the **derived** form (sha `2dde9448…`). Both are asserted per leg in
`run.py` and recorded in each `.done`.

---

## (a) SCORED — the diagnostic repairs are diagnostics-only
`R2_removal_events.csv` **byte-identical** to
`d2z_stepheight/output/S2_removal_events.csv`. **Bar: identical.** This is what
licenses scoring `R2d`/`R1d` against the frozen `S2`/`S1`.

## (b) SCORED — the estimator, with BOTH populations specified identically
Live `delta_mean` against `delta_offline.py` on that leg's own event log.
**Bar: within 10 %**, on `R2w` at r < 25 and r < 30, and on `SW` at r < 50.

⚠⚠ **Both sides specified the same way**, or this fails for population reasons a
third time: never-fired and once-fired columns included at `phi = 1` on both
sides, and the radius applied to the **TALLER** column on both sides.
`d2z2_multicell/delta_offline.py` is a rewrite of D2z's, not a copy.

⚠ **The live `delta_mean` is accumulated ONCE PER PAIR**, at the taller column's
top exposed cell, not once per stepped face. The packet's wording says "once per
stepped face", but the branch runs for every exposed cell of the riser, so a
per-face accumulation would weight an `n`-cell pair `n` times while the offline
estimator weights every pair equally — and (b) compares them directly. The
per-face count is **not** lost: it is the new `delta_faces` column, and
`delta_faces / delta_n` is the mean riser height in cells.

## (c) SCORED — THE RIM TEST: does the derived form reopen what D2z cured?
`R1d` vs `R2d`, disk-mean `col_rate`, single **50–200 s** window.

* **r < 25 mm: bar |gap| ≤ 5 %, and it must reproduce `S1`/`S2`'s +1.39 % to
  within 0.5 pts.** The interior is 0.0 % multi-cell at both meshes, so `n = 1`
  everywhere and the derived form is exact there. A shift means something else
  moved.
* **r < 30 mm: SCORED against the frozen −0.71 %.** **|gap| ≤ 5 % ⇒ the rim cure
  was the population and drain correction and is real. |gap| > 5 % ⇒ part of the
  cure was the double count and the rim needs physics** — report it, do not
  patch it.

## (d) SCORED — the drain meter agrees with the model it reports
`wall_step_q_eff_in / wall_step_q_in` against `step_g_mean`, every closure leg.
**Bar: within 5 %.**

⚠⚠ **Recorded BEFORE launch, because it changes how a failure reads.** These are
two *different* means of the same `g`: the ratio is **flux-weighted**
(`Σ q g / Σ q`) and `step_g_mean` is **face-weighted** (`Σ g / N`). They are
equal only when `g` and `q` are uncorrelated, and they are strongly correlated —
a tall step has both a large `g` and a large `ΔT`. Worse, `q_raw` is **negative**
wherever the neighbour is hotter, so the flux-weighted mean is not even bracketed
by the range of `g`. The unit test measured **0.38** between them on a 4×4 grid.
**The 5 % bar is therefore a statement about that correlation, not about
self-consistency, and it is scored exactly as written.** The genuine
self-consistency check — which does catch the two accumulators drifting onto
different face sets — is the **exact identity** with the gate off:
`wall_step_q_in == wall_step_q` and `wall_step_q_eff_in == wall_step_q_eff`,
bit for bit. **Both are reported.**

**REPORTED alongside:** the ungated ratio and the gated one, so the
0.843-vs-0.335 contamination is visible rather than inferred.

## (e) SCORED — the power identity, on the legs that can carry it
`jet_P_face − (patch_P_robin + side_P_face + side_P_loss)` on `SW`, `SWs`, `SW0`
and `S0`. **Bar: ≤ 0.02 W.** The thermo log prints 6 significant figures, so
three ~10³ W columns quantise at ±0.005 W each ⇒ a 0.015 W floor; `S0` measured
0.0141 W against D2z's ±0.01 W bar, which was written below the instrument.
⚠ **The only check that catches a kernel↔march split** — `side_face_factor`
defaults to 1.0, numerically identical to the derived branch, so a missed call
site is invisible to everything else. ⚠ `jet_P_face` exists only under
`jet_closure`, so this is evaluable **only on the Meier legs**.

## (f) SCORED — identity
`S0` byte-identical to `d2z1_delete/output/RB_D`: thermo, `_removal_events.csv`,
`_jet_profile.csv`. **The header delta must be empty.** ⚠ `metadata` and
`diff.patch` are run provenance (compile date, plot_file path, the ParmParse
table — `query_default` ADDS a key to it — and a start stamp) and are excluded,
exactly as the 2639-file tree list excludes them. The tree gate is run separately
against `witness/d2z1.md5` (**not** `tree0922.md5`).

## (g) SCORED — Meier shape at equal burner depth
`SW` against `A` and `D`, depths capped 50 mm above the burner, **scored only on
rock the nozzle plane has passed**. Report mean |ΔØ|, Ø at 100 and 125 mm, burner
descent, pit lead, volume on **both** bases, and the 2 mm ring table 30–50 mm
**printed per ring, never averaged**.
⚠ **State which rate definition is used.** `A` reads **1.251 m/h** over the
100–150 mm *descent* stretch and **1.296** over the *time* window 150–400 s.

## (h) SCORED — the law A/B, with a floor that is a SYSTEMATIC not noise
`SW` vs `SWs` on (g)'s mean |ΔØ|. **Distinguished only when they differ by more
than 0.25 mm; otherwise UNDECIDED and reported as such.** The floor is a
scoring-depth systematic, from D2x's two readings of the *same* leg (`A` scores
2.27 mm at 150 mm descent, 2.51 at 180). **In a tie the preferred law remains
`wetted`.**

## (i) REPORTED — the conductance regression on Meier
`SW` vs `SW0` (derived vs D2z `g`) on everything in (g). **Reported, not
scored** — the derived form is right on its derivation; this measures how much
the error was worth.

## (j) REQUIRED WITNESS — the multi-cell population, per leg
`wall_step_f_multi`, `side_P_face_multi`, `step_g_max`, and the offline step
populations for every leg. ⚠ **If `wall_step_f_multi` is small on `SW`, (i) is
uninformative and must be reported as such.**

⚠ **`step_g_max` is bounded by `(n+1)/n`, NOT by 1.** `delta = n + phi_S - phi_T`
with both phases in [0, 1], so `delta` spans `[n-1, n+1]` and `g > 1` whenever
the shorter column is due to fire and the taller has just fired — the riser
really is momentarily taller than the `n` cells the grid resolved. The bound is
**2.00 at n = 1 falling to 1.14 at n = 7**; the D2z form's bound **rose** with
`n` (1.99 → 7.73 at a = 0.01). That direction is the repair.

---

## Pre-registered predictions

* **P1** `R2` byte-identical to `S2`. Held / refuted.
* **P2** live-vs-offline **≤ 10 %** on `R2w` (both radii) and `SW`. Held / refuted.
* **P3** the derived form is **exactly inert at 2 mm**: `R2d_removal_events.csv`
  byte-identical to `S2`'s. ⚠ **`R1d` must NOT be identical to `S1`** — the 1 mm
  rim is 26.6 % multi-cell. **If `R1d` is also identical, the change never
  reached a multi-cell face and the test is inert — plumbing, not physics.** The
  two halves are scored separately.
* **P4** the rim **HOLDS**: `R1d` vs `R2d` at r < 30 mm lands **|gap| ≤ 5 %**.
  Refuted if it reopens past 5 % toward the tuned constant's −9.11 %.
* **P5** `wall_step_f_multi` **≥ 20 %** on `SW`. Held ≥ 20 %, refuted < 20 % —
  and if refuted, (j) applies and (i) is reported uninformative.
* **P6** `SW` does **not** recover `A`'s wall band: Ø at 100 mm stays **< 87 mm**.
  At that depth `A` **84.4**, `D` **97.1**, Meier **87**.
* **P7** `SW` burner descent over the 100–150 mm stretch lands **1.3–1.6 m/h**.
  `A` reads **1.251** there. Above 1.6 ⇒ the side path is over-delivering.

## Refuting evidence on disk

* **`A` at 2.27 mm mean |ΔØ| with no side path at all is the number to beat.**
* **Per population.** On **one-cell** steps `SW` must draw **less** side power
  than `D`'s **1234.7 W**; more there means `delta` is inverted — check
  `delta_mean` first. On **multi-cell** steps more than `D` is **expected**.
  (D2x side power: `C` 223.2, `D` 1234.7, `E` 453.8 W.)
* **If `SW` comes out faster than `A`**, check `side_faces > 0` first.
* **No Meier agreement claim tighter than ±5 mm.** Fig. 8.8 is a visible width =
  a LOWER BOUND on diameter; band 85–95 ± 5 mm, rate band 1.3–1.6 m/h.
