# D2z CRITERION — hashed before the first leg, zero outputs on disk

The step-height closure: pay each riser its true sub-cell height
`delta = |n + phi_S - phi_T|` and correct the lateral face conductance by
`g(delta) = (1 - e^{-delta a})/(1 - e^{-a})`, `a = dz*v/alpha`. No free
parameter.

Metric: `col_rate` and the metadata reader are **imported unchanged** from
`studies/d2j0b_plane/analyze.py`. Frozen and re-baselined studies are read by
path and never re-run. **Reference is `studies/d2z1_delete/output/{W2,W1}`** —
the whole-cell pair. ⚠ The frozen `d2l_scan/A_f02{,_1mm}` pair is **half-cell**
and is no longer the reference.

Legs: `S2`/`S1` (`step`, 2/1 mm), `S2w`/`S1w` (`wetted`), `S2f` (`f = 0` with
the drain registered), `S0` (identity), `SW` (Meier, only if the gate passes).

---

## (a) SCORED — THE GATE: the interior mesh gap

Disk-mean rate, **r < 25 mm**, `S1` vs `S2`, on a **SINGLE 50–200 s window**.

* **Bar: |gap| ≤ 5 % = "cured".**
* **Separately pre-registered: |gap| ≤ 2 % = "as good as the fudge it
  replaces."** The whole-cell reference is **+0.38 %**, so a −4 % result is
  *cured but worse than `f = 0.2`* and must not read as a clean pass.

⚠⚠ **AGGREGATE, NOT PER-BLOCK.** A per-block gap inherits block noise from
*both* meshes: D2z-1 measured per-block shifts of +2.59 / −3.42 / +1.00 pts
against an aggregate of +0.085, and block-length sensitivity ran
6.90 → 1.95 → 0.85 → 0.00 pts as blocks lengthened 25 → 150 s **with the mean
invariant**. Per-block is **REPORTED**.

⚠ **NOT r < 30 mm** — that disc includes the 25–30 mm rim stall, where the same
reference reads −6.4 / −8.0 / −14.0 %. Report r < 30 mm alongside so the rim
stays visible.

## (b) SCORED — the drain falls

Flat-floor lateral drain in **watts** (`wall_step_q * wall_step_n * dx*dz`),
`S2` against `S2f`. **Bar: ≥ 8x reduction.**

⚠⚠ **The witness is fixed in this packet.** `WallStepDiagnostic` computed
`kc*(Tc - T_nbr)/dx` from the **raw** cell conductivity and knew nothing about a
factor applied inside `k_face`, so correct code would have failed this. It now
reports **`wall_step_q` (raw) and `wall_step_q_eff` (under the same `g(delta)`
face rule) side by side**; (b) is scored on the **effective** pair.

**Decompose the expected ~14x**, so a failure can be attributed:
* per-face `g(delta)` at the item-0a re-measured `delta = 0.157` (2 mm) is
  **g = 0.25**, i.e. **~4.0x**;
* the rest is population — steps exist **7.9 %** of the time here against
  30–44 % on the `f = 0` cascade, i.e. **~4x**.

Report both terms.

## (c) SCORED — the power identity

`jet_P_face − (patch_P_robin + side_P_face + side_P_loss)` to **±0.01 W** on
every leg. ⚠ **This is the only check that catches a kernel↔march split.**
`:2507`'s `side_f_m` defaults to 1.0 — numerically the same as the derived
branch — so a missed site is invisible to everything else.

## (d) SCORED — THE ESTIMATOR, NOT THE PHYSICS

⚠⚠ **Split deliberately.** Barring the live `delta` against a band measured on
the `f = 0.2` legs would repeat the first draft's mistake: the closure exists to
change how steps evolve, so if it heals *better* the step height drops **below**
the band and **correct code would fail.**

* **SCORED (estimator):** live `delta_mean` against the **offline both-column
  estimator run on that same leg's own event log**
  (`delta_offline.py`). **Bar: within ~10 %** — the gap is
  trailing-interval versus hindsight-interval, nothing else.
* **REPORTED (physics):** the same leg's `delta` against item 0a's re-measured
  whole-cell values, **0.137–0.171 (2 mm)** and **0.197–0.319 (1 mm)**,
  aggregate **0.157** and **0.270**. **A drop below that is a RESULT, not a
  failure.**

⚠ **`delta_mean` is a thermo variable over the LIVE face population and is
therefore TIME-WEIGHTED, not at-birth.** Stated here in words because that is
the exact mistake the first draft made.

`delta_fallback_cols`: **bar < 5 %** of side columns. Item 0a measured the
offline fallback fraction at **0.00 %** in every block at both meshes.

## (e) SCORED — WHICH AREA LAW, WITH A TIE RULE

`S1` vs `S1w` on (a)'s gap, and `S2` vs `S2w`. The law whose gap is inside ±5 %
wins.

⚠⚠ **"Smaller gap wins" is meaningless below the aggregate's own floor**, which
is roughly **0.5–1 %**. **Pre-registered: the laws are distinguished only when
their gaps differ by MORE than 1 pt; otherwise the outcome is UNDECIDED and
reported as such.** **In a tie the preferred law is `wetted`**, on physical
grounds stated in advance — a phase lag on a smooth floor has no vertical face.
If neither law is inside ±5 %, both are refuted and the answer is D2j's
continuous front. Report `riser_area_ratio` alongside.

## (f) SCORED — identity

`S0` byte-identical to `studies/d2z1_delete/output/RB_D` over the interval both
completed. ⚠ **A truncated run is a PREFIX** — compare only where both have
data. ⚠ **The 2639-file tree gate cannot replace this**: no input in that set
turns `side_face_flux` on. The tree gate is run separately against
`witness/d2z1.md5` (**not** `tree0922.md5`; see `witness/REFERENCE.md`).

## (f2) REQUIRED WITNESS — the step populations, across every leg

Fraction of time a step exists, fraction ≥2 cells deep, step height at birth,
and lifetime — for the closure legs, the whole-cell reference pair and the
`f = 0` control (`S2f`). Read-only on the event logs.
⚠ **A gate pass with a growing multi-cell fraction is a TRANSIENT pass**, and
(a) alone cannot tell the difference.
⚠ Item 0a measured **≥2-cell = 0.0 %** on the arbiter at both meshes, so the
multi-cell population the packet worries about (53.2 % on Meier `A`) **does not
exist on this geometry**. (f2) is therefore decisive only on `SW`.

## (g) REPORTED — the front-proxy error

The stated **1.35–1.75** factor against a true continuous front, recomputed from
`S2`/`S1` rather than carried as a claim.

## (h) REPORTED — Meier (`SW` only)

Against `A` and `D`, **at equal burner depth**, depths capped 50 mm above the
burner: shape, Ø at 100/125 mm, burner descent, pit lead, volume on **both**
bases, and the ring table at 2 mm rings 30–50 mm **printed per ring, never
averaged**. ⚠ **State which rate definition is used** — `A` reads **1.251 m/h**
over the 100–150 mm *descent* stretch and **1.296** over the *time* window
150–400 s.

---

## Pre-registered predictions

* **P1** the gate holds: aggregate |gap| **≤ 5 %** held, **> 5 %** refuted.
* **P2** `delta_mean` inside item 0a's re-measured band (0.137–0.171 at 2 mm,
  0.197–0.319 at 1 mm) held, outside refuted. ⚠ Scored as an **estimator**
  check per (d); a drop below the band from better healing is a RESULT.
* **P3** drain reduction **≥ 8x** held, **< 8x** refuted.
* **P4** the **`step`** law passes the gate and **`wetted`** does not held; the
  reverse, both, or neither refutes it.
  ⚠⚠ **Read the outcome correctly either way.** `step` delivers about the same
  heat as `f = 0.2` did while cutting the per-face drain ~4x, so it should beat
  the `f = 0.2` legs. `wetted` delivers almost no heat with the **same** drain
  cut, so its heal-to-drain ratio is about a third of `f = 0.2`'s and it is the
  more likely to reopen. ⇒ **A `wetted` PASS is the informative surprise** — it
  would mean the **drain correction alone cures the cascade**, the physically
  clean result. A `step` pass mostly **relabels `f` as a measured quantity**.
  Say which happened; do not report a `step` pass as though it settled the
  physics.
* **P5** `SW` does **not** recover `A`'s wall band: Ø at 100 mm stays **< 87 mm**
  held, **≥ 87** refuted. For the record at that depth: `A` **84.4**, `D`
  **97.1**, Meier **87**.

**Refuting evidence on disk — per population.** On **one-cell steps**
(`delta` ≈ 0.05–0.20) `SW` must draw **less** side power than `D`; if it draws
more there, `delta` is inverted — check `delta_mean` first. On **multi-cell
steps** more power than `D` is **expected**, because `step` charges the full
vertical face where `D` paid `RiserWeight`'s ramp. And if `S2` comes out
*faster* than `W2` overall, the key is not reaching the code — check
`side_faces` first.
