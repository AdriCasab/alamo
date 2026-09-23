# D2o-0b results — prescribed floor/wall heating, field bounded in radius (2026-09-22)

**Q1 FAILS, so by the packet's rule nothing else is scored. But the failure is a different one
from D2o-0's, and the diagnosis is the deliverable.**

- **The packet's radial-corner constraint was right, and it worked.** Frozen carriers fell from
  **15 of 162 to 2 of 162**, `foot_stall_time` from **219 s to 11 s**, ROP from **0.108 to
  0.771 m/h**, feet depth at the end from 40 to 78 mm. The tool no longer locks.
- **But a sharp radial corner is itself a step-face cascade generator.** The h discontinuity
  (700 → 40, a factor 17) produces a **+402 K lateral temperature step** exactly at the
  load-bearing annulus edge. The conduction sink `k·ΔT/dx` takes **302 kW/m²** of a q_pin ≈
  520 kW/m² input at 2 mm — **58 %** — which is why the burner descends at half rate. The control's
  step is +77 K, taking 11 %.
- **The constraint the next packet inherits is therefore stronger than the one D2o-0 bought: *any*
  sharp corner in `h` does this, radial or stand-off. The transition must be smooth over several
  cells.** Making the corner radial traded a carrier-freeze for a lateral-conduction drain.
- **P5's trigger is aimed at the wrong quantity.** P5 watches for the wall within 50 K of 821 K;
  the wall here is at 417 K so P5 never fires — yet the mesh-divergence risk is real and comes from
  the **lateral gradient at the corner**, not from the wall's proximity to threshold.

## 0. Hashes, binary, provenance

- **Binary** `bin/mmwspalling-3d-g++` = `62dea451e6b1d1908e0ea12f730cdc3f7c01790f19d6d63b8c2b88befe0f67fa`,
  **unchanged at the end**. No source edit, no build, no `make`, nothing committed.
- **`PREDICTIONS.md`** sha256 `d2ee12218b4d84f3c1c249f7df05bacbef6c5c9525907d4bee8bd0d3b513ec0e`,
  read-only, **15:49:55**, written after G1 passed and before `output/` existed. The three 2 mm
  legs completed 16:05.
- **Frozen, imported by path, none edited:** `d2i_gate/` (control and LI0 keys), `d2j0b_plane/`,
  `d2l_scan/analyze.py`, `d2m_feet/analyze.py`, `d2b_feet_rop/analyze.py`, and
  **`d2o0_wallsplit/`** — copied and fixed, never run in place. **The control was not re-run.**

### Every packet reference reproduced with this study's estimator first

Gas power **2710.7 W** shallow / **1390.6 W** deep, radial splits **41.9/29.9/22.9/5.3 %** and
**64.1/35.9/0/0 %**, and all 22 `T_gas(r)` values — all match the packet exactly.
**One disagreed and was resolved before launch:** the packet's control depth-mean Ø of 110.5 mm.
My estimator gives **100.7 mm** over the control's whole drilled length, because depth-mean is not
a single number for a cone — it runs 112.3 mm at a 97 mm averaging depth down to 100.7 at 236 mm,
and 110.5 corresponds to ~105 mm. **P2 is therefore scored at the matched depth with the control
recomputed there**, never against a fixed 110.5 (`PREDICTIONS.md` §0).

### G1 — field integrity, passed before launch

| t [s] | prescribed | control field, same estimator | ratio | r>70 mm share (mine / control) |
|---|---|---|---|---|
| 50 | 6133 W | 12391 W | **0.49×** | **8.4 % / 28.4 %** |
| 150 | 4814 W | 7968 W | 0.60× | 0.0 % / 0.0 % |
| 300 | 4142 W | 5783 W | 0.72× | 0.0 % / 0.0 % |

Floor zone (r ≤ 40 mm) **bit-identical** to the control's field at every time (5344 / 4553 /
4025 W both). `T_gas(r)` fit worst error **24.0 K** at r = 48 mm (criterion 60), monotone over
0–250 mm, → 293.15 K, **guard at ambient** so it can never set a physical value. D2o-0's 1400 K
blanket is gone.

**Two calls declared rather than handled quietly:**
1. **G1 uses one estimator on both fields.** Reading the top solid cell's temperature from a
   plotfile is *not* the pinned surface temperature the solver uses, so absolute powers are biased
   **4.2–4.6× high**; evaluating the control's own field identically cancels the bias in the ratio.
   This was found *by* G1 — the first version compared my estimator against the control's logged
   march and read a spurious 2.3–3.1× excess.
2. **The packet's "total within 20 % of the control" cannot be met by a packet that mandates a 20×
   wall cut.** My total is 0.49–0.72× — a *deficit*, entirely in r > 40 mm, floor bit-identical.
   G1's stated purpose is catching an **excess** (D2o-0's 2.9×). **Treated as passed on its
   purpose**, and flagged for the reviewer.

### The 1 mm leg was ended early, on the user's instruction

Stopped at **t = 56.9 s of 600** (first cut from 600 s to ~260 s, then ended). No `.done` was
written and none was fabricated; see `output/W_1mm_TRUNCATED.md`. **Q1 already fails at 2 mm and
refinement can only make the step-face sink worse**, the packet requires both meshes only for a
**positive** claim, and shape is unscoreable regardless. **This study is single-mesh and nothing in
it is converged.**

## 1. Q1 — the gate

| case | h_wall | ROP 150–250 s [m/h] | vs control | feet depth at 600 s | `foot_stall_time` end | Q1 |
|---|---|---|---|---|---|---|
| control | — | **1.491** | — | 236 mm | 1 s | — |
| `W_2mm` | 40 | **0.771** | 0.52× | 78 mm | 11 s | **FAIL** |
| `Wsens_2mm` | 20 | **0.522** | 0.35× | 77 mm | 14 s | **FAIL** |
| `Wz_2mm` | 0 | **0.340** | 0.23× | 39 mm | 47 s | **FAIL** |

Q1 needed stall < 10 s, feet ≥ 150 mm by 600 s, ROP within 30 % of 1.491. **Every case misses on
depth and rate; `W_2mm` misses the stall bound by 1 s.** Against D2o-0's 0.108 m/h and 219 s stall
this is a large improvement that still does not clear the gate.

## 2. Which columns — the diagnosis the rule asks for

**(a) The carrier freeze is essentially cured.** Annulus columns (r ∈ [28, 40] mm) with s ≤ 0,
frozen at the emulated zero: **2 of 162** at 435–438 K, against **15 of 162** at 444–482 K in
D2o-0. The radial corner does what it was designed to do.

**(b) The corner is now the problem.** At t = 600 s in `W_2mm`:

| ring | mean surface T | |
|---|---|---|
| r 36–40 mm (firing, inside the corner) | **~820 K** | 55 of 162 annulus columns at or above 815 K |
| r 40–44 mm (just outside the corner) | **~418 K** | |
| **step across the corner** | **+402 K** | control: **+77 K** |

| mesh | k·ΔT/dx | as a fraction of q_pin ≈ 520 kW/m² |
|---|---|---|
| 2 mm | **302 kW/m²** | **58 %** |
| 1 mm (projected) | **604 kW/m²** | **116 % — exceeds the input** |

**Partial mesh confirmation, indicative only.** From the truncated 1 mm leg, at the common time
t = 50 s:

| mesh | T(r 36–40) | T(r 40–44) | step | k·ΔT/dx |
|---|---|---|---|---|
| 2 mm | 753.7 K | 598.6 K | +155.1 K | 116 kW/m² |
| 1 mm | 796.8 K | 643.0 K | **+153.7 K** | **231 kW/m²** |

**The temperature step is mesh-independent (155.1 vs 153.7 K) while the sink it drives doubles.**
That is the D2j-0b step-face signature exactly: a physical ΔT, a conduction drain that scales as
1/dx. t = 50 s is early transient and **nothing here is scored** — it is reported as support for
the diagnosis, not as a mesh result.

## 3. Shape — observed, NOT scored (Q1 failed)

Per the decision rule, these are recorded for the next planner and are **not** scored. No case
reached the 150 mm feet-depth floor, so P1 is unscoreable in any event; the matched depth collapsed
to 39 mm (set by `Wz_2mm`).

| case | depth-mean Ø | mouth Ø | V̇ [cm³/s] | `patch_P_robin` |
|---|---|---|---|---|
| control | 109.8 | **128.2** | **4.58** | 1639 W |
| `W_2mm` | 79.1 | **79.3** | **2.00** | 1078 W |
| `Wsens_2mm` | 78.2 | 78.6 | 1.88 | 941 W |
| `Wz_2mm` | 78.1 | 78.4 | 1.63 | 735 W |

**The funnel is gone** — mouth 79 mm against the control's 128 — and V̇ 2.00 sits in the band Meier
occupies (2.47). **But this must not be read as success:** the hole is 79 mm because it is exactly
the skirt bore, and the rate is halved. It is the same "right number for the wrong reason" as D2m's
`R44`. A cylinder of the correct 85–95 mm was not produced.

**P4 is unscoreable for a reason worth recording: the geometric wall band is EMPTY.** No rock
outside r = 40 mm was ever excavated, so there is no wall to measure — the hole is the skirt bore
and everything beyond it is untouched original surface. **This vindicates the geometric band
definition**: D2o-0's `s`-only band would have counted those thousands of untouched columns as
"wall" (it counted 3298 and reported a meaningless 333 K).

Ledger ≤ 5.6e-15 on every run.

## 4. Pre-registration scorecard

| prediction | outcome |
|---|---|
| **Q1 passes** (radial corner stops the jam) | **REFUTED on the gate, upheld on the mechanism.** Carriers 15 → 2, stall 219 → 11 s, ROP 0.108 → 0.771 — the freeze is cured, but a new drain halves the rate |
| P1 fails low, hole cut at ≈ 80 mm = the stance | **HELD** — 78–79 mm, though formally unscoreable |
| P2 mouth passes (≈ 85–95 vs control 145) | **HELD, and better than predicted** — 79 mm |
| P3 passes at ≈ 2.1 cm³/s | **HELD** — 2.00 cm³/s (not scored: Q1 failed) |
| P4 fails low; wall ≈ 560–620 K | **unscoreable** — the band is empty, no wall exists |
| P5 does not fire | **HELD by its letter**, but **the trigger is aimed at the wrong quantity** |
| the wall does not spall, v(r > corner) ≈ 0 | **HELD** |
| **the corner-driven cascade** | **NOT PREDICTED.** I predicted an inert wall and a narrow hole; I did not predict that the corner's own lateral step would drain the load-bearing ring |

The hand estimate (588 K) could not be checked: the geometric band is empty.

## 5. What the next packet inherits

1. **The design constraint, now in its general form: no sharp corner in `h`, anywhere.** A
   stand-off corner freezes the carriers (D2o-0); a radial corner drains them laterally (here).
   **The floor→wall transition must be smooth over several cells**, so that no cell-to-cell ΔT
   large enough to drive a 1/dx conduction sink ever forms at the load-bearing radius. This is the
   single most reusable result of D2o-0 and D2o-0b together.
2. **Any future wall closure must be checked for the lateral step, not only for the wall
   temperature.** P5's 50 K-from-threshold trigger would have passed this run silently. The right
   diagnostic is `k·ΔT/dx` at the transition against q_pin.
3. **A grazing wall cannot widen the hole** — confirmed again: the hole is exactly the skirt bore,
   and the wall band is empty because nothing outside r = 40 mm was ever removed.
4. **Owed:** a 1 mm leg for any positive claim from this configuration; and a re-test once the
   transition is smoothed, which is the first thing to try — it is key-only and cheap.
