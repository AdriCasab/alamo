# D2m — predictions and decision rule
# (AMENDED 2026-09-22 12:22:10 after a mid-flight packet change; never edit again)

## AMENDMENT RECORD (read first)

**What happened.** `ACTIVE_STEP.md` was edited at **12:16:08**, after this file was first written
and hashed (**12:15:53**, sha256
`6d1daf69895a0ba52c17a33774d528383309b5bcf4818d2da91f123fc8641b6b`) and as the five runs were
launching. The packet grew 160 → 182 lines.

**What the packet changed.** Nothing about the runs; everything about what is *deciding*:
- item 4 gained a bullet making **the transition radius deciding alongside min Ø**: "the radius at
  which the recession rate falls to half its central value tracks `foot_r_outer` within one cell
  (2 mm) under H-feet, and does not move under H-thermal";
- a new item **5a** specifies the recession-rate profile as "the measurement everything else hinges
  on": 2 mm annular bins from the axis to **70 mm**, D2j-0's estimators (central rate, half-value
  radius, 10–90 % fall-off width), taken from the **removal log** and not from plotfile
  differencing, plotted as `d2m_vprofile.png` with each case's `foot_r_outer` marked and **Meier's
  implied floor half-width 42.5–46.5 mm** shaded;
- item 7 gained the corresponding `RESULTS.md` §4 and the second figure.

**What this file changed.** §5 previously carried a transition prediction as a *non-deciding*
extra, binned at 4 mm and defined against half of *peak*. It is now **deciding**, binned at 2 mm
out to 70 mm, and defined against half of the **central** value, with D2j-0's estimators — §5
below. Everything else in this file is **byte-unchanged**: the cases, the hypotheses H-feet /
H-thermal / H-asym and their numbers, the ROP direction and its sizes, the `S3648` prediction, the
scoring definitions, and the decision table.

**What was done about the runs.** The amendment changes no input key, so the five runs were
physically unaffected. They were **killed and deleted anyway**, per the standing rule ("if amended,
re-hash, kill and delete every affected run, and record it"). They were killed at t ≈ 22–28 s of
250, **no `.done` was written, no analysis was run, and no output was inspected**; `output/` was
removed before this file was rewritten. The relaunch therefore post-dates this hash entirely.

---

Key-only. No source change, no build, no `make`. Binary `bin/mmwspalling-3d-g++` =
`62dea451e6b1d1908e0ea12f730cdc3f7c01790f19d6d63b8c2b88befe0f67fa`, to be unchanged at the end.

**The question.** D2l's cut diameter was **78.0 / 77.1 / 77.7 / 77.5 / 78.0 mm** — flat to 0.9 mm
while absorbed power doubled and the rate changed 2.7×. A thermal diameter would have moved. The
burner's feet span r ∈ [28, 40] mm, i.e. **stance Ø 80**, within 2 mm of the observed cut. Does the
cut diameter follow the stance?

**Control:** `R40` = the D2i `LI0` keys at 2 mm, side flux off, stop 250 s — ROP 1.491 m/h,
V̇ 4.58 cm³/s, min Ø 78.0 mm, depth-mean Ø 107.6 mm, absorbed 1639 W (quarter). These are an
**identity check, not a prediction**.

## 1. Two couplings declared before launch (found by reading the source, not by running)

ACTIVE_STEP says "change only `foot_r_inner` / `foot_r_outer`" and "not the gas closure". Neither
key is confined to the support rule:

1. **`foot_r_inner` silently sets `nozzle_collision_radius`** (`MMWSpalling.H:5695`,
   `feet ? foot_r_inner : radius`), and the D2i keys never set it. Uncorrected, `S3648` would also
   move the collision guard 28 → 36 mm. **Every case pins
   `surface_patch.nozzle_collision_radius = 0.028`** — the control's own effective value, so the
   control is unchanged and the coupling is removed from the ladder. The identity guard proves the
   pin is inert.
2. **`foot_r_outer` also defines the free-surface boundary** (`MMWSpalling.H:2062`: free-surface
   when `r > foot_r_outer` and inside the aspect cone measured from `r − foot_r_outer`). There is
   no separate key, so this cannot be pinned away. **Bounded from D2l rather than by extra runs:**
   `P_fsoff0` removed the free surface *entirely* and moved min Ø by **+0.3 mm** and ROP by
   **−2.3 %**. Moving `foot_r_outer` only modulates it partially, so the effect on **min Ø is well
   under 0.3 mm** — negligible against the 64 → 96 mm swing H-feet predicts — but **material at the
   few-per-cent level for ROP**, which is why the ROP predictions below are stated with that
   caveat and `jet_fs_cols` is reported per case.

## 2. The hypotheses

- **H-feet** (ACTIVE_STEP): min Ø ≈ **2 × `foot_r_outer` within 5 mm** — ≈ 64 / 78 / 88 / 96 mm at
  r_outer 32 / 40 / 44 / 48.
- **H-thermal** (ACTIVE_STEP): min Ø stays **78 ± 3 mm** at every r_outer.
- **H-asym (mine, and what I actually expect).** Neither, because the mechanism is not symmetric:
  - **Outward (R44, R48): the stance binds.** The burner rests on the ring, so it cannot descend
    faster than the ring's 0.9 quantile erodes. Rock at larger r sees lower `h`, erodes more slowly,
    holds the burner higher, and *must* be cut for descent to continue. **min Ø should track
    2 × r_outer outward, to within a cell or two: ≈ 86–88 at R44 and ≈ 94–96 at R48.**
  - **Inward (R32): the stance does NOT bind.** `h = h_expr(r, s)` does not depend on the foot
    radii, so rock at r ∈ [32, 40] is heated *exactly as before* and should still be removed — it
    simply stops being load-bearing. **min Ø should therefore NOT collapse to 64; predict
    74–80 mm**, i.e. close to the control.
  - **So the prediction is asymmetric tracking: min Ø ≈ 2 × r_outer for r_outer ≥ 40, and roughly
    flat below it.** This is falsified by R32 coming out near 64 (pure H-feet) or by R44/R48 staying
    near 78 (pure H-thermal).

### 2a. Why min Ø alone is weak evidence — stated before seeing it

`z_foot` is the **0.9 nearest-rank quantile** of the annulus face heights, i.e. the burner rests
near the *top* of the rock under its pads, so ~90 % of annulus columns out to `foot_r_outer` lie at
or below `z_foot`. The hole at the feet depth is therefore **at least 2 × `foot_r_outer` wide
almost by construction**, and min Ø over [0, z_m] sits at the narrowest (deepest) point, which is
the feet depth. **So H-feet's outward branch is close to a structural identity of the support rule,
and confirming it carries little information.** The two informative tests are:
- **R32 (inward)**, which the support rule does *not* force; and
- **the transition radius (§5)**, which is the flame's own signature and is not constrained by the
  support rule at all.

This is why the amendment's promotion of the transition radius to a deciding metric is right, and
it is recorded here so that a confirmed H-feet is not over-read.

## 3. ROP direction — the reasoning's own sign test

A larger ring reaches into cooler, less-eroded rock, so its 0.9 quantile sits higher, the burner
rides higher, stand-off grows and heating falls.

- **Predict ROP falls monotonically as `foot_r_outer` rises.** Sizes, stated so they can be wrong:
  R32 **+10 … +25 %** over the control; R44 **−10 … −25 %**; R48 **−25 … −45 %**.
- **A rising ROP at larger r_outer refutes the reasoning even if the diameter tracks**, and will be
  reported as such rather than absorbed into a "H-feet holds" verdict.
- Because of coupling 2 above, **ROP differences smaller than ≈ 3 % are not attributable to the
  feet** and will not be read as signal.

## 4. `S3648` vs `R48` — does the outer edge alone set it, or the whole ring?

Both have r_outer = 48; `S3648` drops the inner 28–36 mm band. The inner rock is the **fastest
eroding** (highest `h`), so it contributes the **lowest** z_face values, and a 0.9 nearest-rank
quantile is taken near the **top** of each pad. Removing low values should barely move a high
quantile.

- **Predict `S3648` ≈ `R48`: min Ø within 3 mm and ROP within 5 %.** That is the "outer edge sets
  it" reading. A material difference means the whole ring matters, and the two must be reported
  separately rather than averaged.

## 5. The recession-rate profile — DECIDING (amended; was non-deciding)

**Estimators, D2j-0's, reused unchanged** (`d2j0_hotdisk/analyze.py`, copied into this study rather
than run in the frozen folder):
- per-column recession rate in **m/h** from the **removal log** (Σ `h_applied` over the window ÷
  window), never plotfile differencing;
- **2 mm annuli**, `EDGES = arange(0, 0.0701, 0.002)` — the axis out to **70 mm**;
- **central value `v_c`** = mean per-column rate over **r < 6 mm**;
- **`r_half`** = `cross(prof, 0.5·v_c)`, the **outermost** radius where the binned profile falls
  through the level, linearly interpolated between bin centres;
- **`w_edge`** = `r(0.1·v_c) − r(0.9·v_c)`, the 10–90 % fall-off width.

**Window.** The profile is read at **matched geometry**: the 100 s window ending at each case's
`t_m`, the time its feet reach the common depth `z_m`. The fixed [150, 250] s window is reported
alongside as a cross-check; **if the two disagree on `r_half` by more than one cell (2 mm), both
are reported and no single reading is forced.**

**The deciding prediction:**
- **H-feet:** `r_half` tracks `foot_r_outer` **within one cell (2 mm)** — ≈ 32 / 40 / 44 / 48 mm.
- **H-thermal:** `r_half` does not move — fixed within 2 mm across the ladder.
- **H-asym (mine):** `r_half` is set by `h(r, s)`, not by the support rule, so it **stays near the
  control's value in every case** — I expect **H-thermal to hold on `r_half` while H-feet holds on
  min Ø outward**. That combination is the sharpest form of H-asym: the stance sets a *floor* on
  the cut by forcing removal of rock it stands on, while the flame's own erosion edge does not
  move. If `r_half` does track `foot_r_outer`, H-asym is refuted on its central claim.

**Meier's implied floor half-width is 42.5–46.5 mm** (Ø 85–93) and is shaded for orientation. It is
**not** scored here — no Meier number comes out of this packet.

## 6. Scoring — fixed here, before any result

- **Burner ROP** = least-squares slope of `nozzle_z` over thermo rows in [150, 250] s, m/h. The
  D2f/D2g/D2h/D2i/D2k/D2l estimator, unchanged.
- **Ø(z)** = 2 × azimuth-mean `r_wall`, 9 sectors (the `d2c_steady/score.py` definition, reused
  through `d2l_scan/analyze.py` by path), reported **as a full profile at 20 mm intervals**, never
  at a single depth (D2l: the profiles cross).
- **min Ø** = the minimum of that profile over [0, z_m], **and the depth at which it occurs**.
- **Matched geometry:** every case is read at the time its feet reach **z_m**, the deepest feet
  depth every case reaches. **Never at matched time** (D2l: Ø at matched time tracks dwell).
- **V̇** = 4 × Σ `h_applied`·dx·dy over removal rows in [150, 250] s ÷ 100 s.
- **J/mm³ to rock** = 4 × `jet_P_face` / V̇, reported **against the thermodynamic floor
  1.147 J/mm³ only**. The 4–6 band is an assumed 30 % delivery, not a measurement, and is **not**
  used as a target.
- **No Meier score comes out of this packet** — case-to-case differences at fixed mesh only.
- Nothing is re-windowed, re-binned, or re-chosen after results are seen.

## 7. The decision

| outcome | reading |
|---|---|
| min Ø **and** `r_half` track 2 × r_outer in **both** directions | **H-feet.** The hole is the width of its own tool stance; Meier's was wider than his. The next packet must target what removes rock **beyond** the feet, and **no Meier diameter may be quoted until that exists**. |
| min Ø flat at 78 ± 3 and `r_half` fixed | **H-thermal.** The feet hypothesis is dead, the width is thermal, the gas-path packet proceeds with the wall as its target. |
| **min Ø tracks outward but `r_half` does not move** — H-asym | **Mixed, and reported as mixed.** The stance sets a *floor* on the cut by forcing removal of the rock it stands on, but the flame's erosion edge is where the width really comes from. Both the support rule and the gas path are implicated, and the next packet needs the **radial heat distribution**, not just the wall. |
| `S3648` ≠ `R48` | the whole ring matters, not the outer edge; report separately. |

Reported, not deciding: absorbed power, depth-mean Ø, `jet_s_c`, `patch_min_standoff`,
`jet_fs_cols`, ledger closure, cap ratio, `v_c`, `w_edge`.
