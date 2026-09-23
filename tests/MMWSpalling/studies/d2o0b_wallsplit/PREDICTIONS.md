# D2o-0b — predictions and decision rule
# (written and hashed BEFORE any run; never edit)

Key-only. No source change, no build, no `make`. Binary `bin/mmwspalling-3d-g++` =
`62dea451e6b1d1908e0ea12f730cdc3f7c01790f19d6d63b8c2b88befe0f67fa`, to be unchanged at the end.

## 0. Every reference value in the packet, reproduced with this study's estimator first

Per the guardrail, before trusting anything:

| packet states | my estimator | verdict |
|---|---|---|
| control face power, shallow (t 20–100 s) **2711 W** | **2710.7 W** | ✓ |
| radial split shallow 41.9 / 29.9 / 22.9 / 5.3 % | **41.9 / 29.9 / 22.9 / 5.3 %** | ✓ |
| control face power, deep (t 150–450 s) **1391 W** | **1390.6 W** | ✓ |
| radial split deep 64.1 / 35.9 / 0 / 0 % | **64.1 / 35.9 / 0 / 0 %** | ✓ |
| T_gas(r) table, 1824 → 624 K over r 0–168 mm | reproduced to the digit at all 22 radii | ✓ |
| control ROP 1.491 m/h, V̇ 4.58 cm³/s, mouth 145 mm | 1.491, 4.58, 144.9 | ✓ |
| control **depth-mean Ø 110.5 mm**, min Ø 79.3 | **100.7 / 79.0** over its whole drilled length | ✗ **resolved below** |

**The depth-mean discrepancy is real and it changes how P2 must be scored.** "Depth-mean Ø" is not
a single number for a cone — it depends entirely on how deep you average:

| averaging depth | 97 mm | 110 | 131 | 150 | 180 | 200 | 236 (whole hole) |
|---|---|---|---|---|---|---|---|
| control depth-mean Ø | 112.3 | 109.0 | 107.6 | 105.8 | 103.5 | 102.2 | **100.7** |

The packet's 110.5 corresponds to averaging over only ~105 mm. Averaging over the control's whole
drilled length gives **100.7 mm — already at the edge of P2's 85–100 target.** **So P2 is scored at
the matched feet depth common to all cases, with the control recomputed at that same depth**, never
against a fixed 110.5. Stated here, before any result, because otherwise the target moves with how
deep the test happens to drill.

## 1. The field, and G1 (run before this file was hashed)

**Radial corner at r = 40 mm** (= `foot_r_outer`): floor = the control's own `h(r, s)` verbatim for
r ≤ 40 mm; wall = h = 40 W/m²K for r > 40 mm; `s ≤ 0` and `r > 170 mm` emulated zero (1e-3).
**The carriers live at r ∈ [28, 40] mm, so every one of them is in the floor zone by
construction** — which is the design constraint D2o-0 bought, and **Q1 is its test**.

**T_gas(r)** = `293.15 + 346.5465·exp(−(r/0.0671158)^5.89) + 1173.4772/(1+(r/0.1009858)²)`, fitted
on the control's **full** tabulated range 0–168 mm: **worst error 24.0 K at r = 48 mm** (criterion
60), monotone over 0–250 mm, → 293.15 K as r grows. **The parser guard is ambient**, so it can
never set a physical value — the exact defect that broke D2o-0.

### G1 result, as printed by `run.py --list` before launch

| t [s] | prescribed | control field, same estimator | ratio | control logged | estimator bias | r>70 mm share (mine / control) |
|---|---|---|---|---|---|---|
| 50 | 6133 W | 12391 W | **0.49×** | 2701 W | 4.6× | **8.4 % / 28.4 %** |
| 150 | 4814 W | 7968 W | **0.60×** | 1884 W | 4.2× | 0.0 % / 0.0 % |
| 300 | 4142 W | 5783 W | **0.72×** | 1328 W | 4.4× | 0.0 % / 0.0 % |

Floor zone (r ≤ 40 mm) is **bit-identical** to the control's field at every time
(5344 / 4553 / 4025 W both) — the floor really is unchanged, which is the design requirement.

**Two things declared rather than quietly handled:**

1. **The comparison uses one estimator on both fields.** Reading the top solid cell's temperature
   from a plotfile is *not* the pinned surface temperature the solver uses
   (`robin_form = pinned` sets T_s = min(T_face, T_pin)), so absolute powers from this estimator
   are biased **4.2–4.6× high**. Evaluating the control's own field the same way makes the bias
   cancel in the ratio. The raw logged total is shown only to calibrate it.
2. **The packet's "total within 20 % of the control" cannot be satisfied by a packet that mandates
   a 20× wall cut.** My total is **0.49–0.72×** — a *deficit*, concentrated entirely in r > 40 mm,
   with the floor bit-identical. G1's stated purpose is to catch an **excess** ("a 2.9× power excess
   makes the whole run meaningless — that is what happened last time"). **I treat G1 as passed on
   its purpose** — no blanket, no excess, far-field share under the control's, T fit inside 60 K —
   and record that the literal two-sided reading is not achievable here. It is not reinterpreted
   silently, and a reviewer should check this call.

## 2. The cheap hand estimate, made first (packet item 7)

Semi-infinite Robin: T_s = T_0 + (T_g − T_0)·[1 − exp(β²)erfc(β)], β = h√(αt)/k, α = 6.904e-7 m²/s,
k = 1.5 W/mK. Wall gas at r ≈ 45 mm is T_g ≈ 1600 K; a column's exposure from s = 50 mm to s = 0 at
the control's 1.491 m/h is **121 s**, giving √(αt) = 9.14 mm and **β = 0.244**:

    T_wall ≈ 293 + 1307 × 0.226 ≈ **588 K.**

(The same estimate predicted D2o-0's wall to within 13 K — 578 against 565 measured.)

## 3. Predictions

**Q1 — the gate: no jam.** `foot_stall_time` < 10 s at end of run, feet depth ≥ 150 mm by 600 s,
steady ROP within 30 % of 1.491 m/h.
→ **I predict Q1 PASSES.** The radial corner puts all of r ∈ [28, 40] in the floor zone, so no
carrier can freeze. This is the direct test of the constraint D2o-0 bought.

**P1 — shaft Ø, whole profile 80–100 mm and depth-mean 85–95 mm below 100 mm depth.**
→ **I predict P1 FAILS, at the low edge.** The floor field now stops at exactly r = 40 mm and the
wall is inert (below), so the hole is cut at **≈ 80 mm — the stance**, which is the bottom of the
80–100 window and **below** the 85–95 mean. D2m established min Ø = 2·`foot_r_outer` to within
2 mm; nothing here widens it.

**P2 — depth-mean Ø 85–100 mm over the drilled length, mouth Ø ≤ 110 mm.**
→ **Mouth PASSES** (≈ 85–95, against the control's 145): the flat-start transient can now only act
inside r = 40 mm at floor h, so there is no wide start-up crater. **Depth-mean marginal and
probably just below 85.**

**P3 — volume rate 2.0–3.2 cm³/s.**
→ **I predict P3 PASSES, at ≈ 2.1 cm³/s.** (π/4)·(80 mm)²·1.5 m/h = 2.10 cm³/s, against the
control's 4.58 and Meier's 2.47. **Note the tension this creates: P3 passes *because* the hole is
narrow, which is the same reason P1 fails.** That is the mirror of D2m's "right rate for the wrong
reason" and must not be reported as a success on its own.

**P4 — wall 600–800 K and 1.5–3.5 kW**, on a band defined on **geometry** (r > 40 mm, face at least
one cell below the original surface, within the band height below the nozzle plane).
→ **I predict P4 FAILS low on both:** wall ≈ **560–620 K** (the hand estimate is 588 K) against
600–800, and wall heat a few hundred W against 1.5–3.5 kW — the geometric band is a narrow annulus
and h = 40 is 20× below impingement.

**P5 — conditional.** → **Does not fire.** 588 K is ≈ 233 K below 821 K.

**The wall does not spall at all**, so `v(r > 40 mm)` should again be ≈ 0 and the hole should be a
**narrower** cone, not a cylinder — D2o-0's `v(r > 45) ≡ 0` repeated with the corner moved inward.

**What would refute me:** any wall above 700 K, a shaft mean above 85 mm, or `v(r)` non-zero beyond
the corner.

## 4. Scoring — fixed here, before any result

- **ROP** = least-squares slope of `nozzle_z` over [150, 250] s, m/h (the D2f→D2m estimator).
- **Ø(z)** = 2 × azimuth-mean `r_wall`, 9 sectors (`d2c_steady/score.py` via `d2l_scan`), reported
  as a **whole profile**, at **matched feet depth** across cases, never matched time.
- **P1/P2 scored only below 100 mm depth.** A case that does not reach 150 mm of feet depth is
  marked **unscoreable**, not read from its start-up funnel.
- **The control's depth-mean is recomputed at the matched depth** (see §0), never taken as 110.5.
- **V̇** = 4 × Σ `h_applied`·dx·dy over [150, 250] s ÷ 100 s.
- **v(r)**: D2m/D2j-0 estimators, 2 mm annuli to 70 mm, from the removal log, **absolute, never
  normalised**.
- **Wall band is geometric**, per P4 above — **not on `s` alone**, which is what put 3298 columns of
  untouched rock into D2o-0's P3.
- **Absorbed power** from `patch_P_robin` (no `jet_*` columns exist under `jet_closure = none`),
  split by radius.
- Nothing re-windowed, re-binned or re-chosen after results are seen. **`h` is never tuned.**

## 5. The decision

| outcome | reading |
|---|---|
| **Q1 holds and P1 and P2 hold** | the floor/wall split is the mechanism → write the D2o closure, gas-path following, annulus correlation as named, **no free wall parameter**; P3 and P4 become its acceptance targets |
| **Q1 fails** | the zoning still freezes a carrier → report which columns and **stop**; score nothing else |
| **Q1 holds but P1 or P2 fails** | report what a prescribed field could **not** produce, which bounds any closure built on the same idea. **Do not tune the field to make it pass.** |
