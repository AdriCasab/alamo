# D2v — the body-in-the-hole h map: separate the funnel from the wall

**Verdict: (a) PASSES, and by a wide margin.** The funnel is removed:
`M1`'s mean |Ø_model − Ø_Meier| over 10 / 25 / 50 mm falls from **37.2 mm to
3.8 mm** against a bar of 20. Every gated criterion passes.

| gate | result |
|---|---|
| **(a) SCORED — the funnel** | `M1` **3.8 mm** (bar 20); `M0` 37.2, `M2` 42.8 — **PASS** |
| **(c) item 0's witness** | `side_P_face_above` **2363 W** on `M2`, **exactly 0.00** on `M1` (max *and* min) — **PASS** |
| **(d) no blanket** | mean recession beyond r = 76 mm: `M0` 0.216, `M1` **0.000**, `M2` **0.000** mm (bar 2) — **PASS** |
| **P1 byte identity** | `M0` == frozen `d2i_gate/LI0_2mm`: **114 599 event rows and 4 184 thermo rows BYTE-IDENTICAL** |
| domain watch | clearance 38–120 mm (bar 20) — **PASS** |
| tree identity | **2639/2639**, zero differing lines, keys unset, binary pinned |
| regressions | sweep's 8 unit + `jet_d2c` + `dev2d` + `s1_smoke` all rc=0; plus `spall_event`, `side_face_flux`, `convective_patch`, `jet_enthalpy`, `robin_feet` **PASS** |
| `ledger_err` | ≤ **1.24e-14** on all three legs |
| item 3 trace | **PASSED T1–T10 on both branches** before launch |

## 1. The shape

Diameter vs depth at t = 401.6 s [mm]:

| depth | 10 | 25 | 50 | 75 | 100 | 125 | 150 | 175 | 200 |
|---|---|---|---|---|---|---|---|---|---|
| **Meier\*** | 96 | 92 | 88 | 87 | 87 | 87 | 86 | 85 | 85 |
| `M0` stale map | 144.3 | 129.1 | 114.3 | 105.3 | 99.5 | 94.4 | 85.1 | 75.3 | 66.3 |
| **`M1` confined** | **100.4** | **95.4** | **91.7** | **89.6** | **88.7** | **87.8** | **87.0** | **86.1** | **85.2** |
| `M2` + wall heating | 136.9 | 135.1 | 132.4 | 129.9 | 125.5 | 118.8 | 111.7 | 105.4 | 99.4 |

\* a **visible width** = a **lower bound** on Ø, ±5 mm; shape only.

`M0` is a **funnel** — 144 mm at the mouth narrowing to 66 mm at 200 mm depth,
crossing the reference. `M1` is a **near-cylinder** that tracks the reference
within **4.4 mm at the mouth and 0.2 mm at 200 mm**, i.e. inside the reference's
own ±5 mm at every depth. **This is not a validation claim** — see §5.

## 2. Attribution, which was the packet's purpose

**`M1 − M0` is the funnel fix**: −43.9 / −33.7 / −22.6 mm at 10 / 25 / 50 mm.
The cause was located before launch and removed at the boundary condition: item
3's trace measured `max dT = 0.013 K` over the 2463 columns beyond r = 76 mm,
where the stale map runs `h = 362 W/m²K`.

**`M2 − M1` is the wall heating**: **+36.7 / +24.8 / +14.2 mm** at
100 / 150 / 200 mm. P3 held. But the effect is far larger than a wall correction
— see §3.

## 3. ⚠ `M2` widens the hole enormously, and (d) passes for a precise reason

`M2` strips **76.7 mm** from the r = 50–76 mm ring, against `M0`'s 32.5 and
`M1`'s 0.5. Yet its mean recession **beyond** r = 76 mm is **exactly 0.000 mm**.

**So the confinement works exactly where it was built to work, and the widening
happens entirely inside the mouth.** (d) is passed on its own terms, not evaded:
the above-plane columns carrying `h > 10 W/m²K` number 670 and reach only to
**r = 73 mm**, inside `r_out`. That is the packet's design doing its job.

⚠ **But the mechanism behind `M2` is unbounded, and this is not in the packet.**
`JetEnthalpyMarch` skips above-plane columns at both `:2126` and `:2212`. That
is *why* `jet_T_above` survives the march — and it also means **the above-plane
heat is never deducted from the jet's enthalpy budget**. The march's invariants
(`P_face <= P_cap`) cover below-plane columns only, and the global
`ledger_err` cannot see it because the flux is still counted in `P_robin`.
`M2`'s above-plane source is therefore capped only by the radial confinement,
never by the gas's own enthalpy. **A future packet that widens the annulus band
must fix this first.**

## 4. ⚠ P6 refuted on its own stated test, and the step

**P6 predicted a numerical floor of ~96–100 mm at 10–50 mm depth, refuted if
`M1` or `M2` cut below 92 mm. `M1` cuts 91.69 mm at 50 mm depth.** Refuted —
narrowly, and in the *favourable* direction: the 16 mm blend does not impose as
hard a floor as item 1b argued. The floor argument was sound in form but ~5 mm
pessimistic in magnitude.

**The lateral step at the wall**, which the packet asked to be watched:

| leg | largest step in `z(r)`, r < 76 mm | median step | lateral ΔT, max | `k·ΔT/dx` |
|---|---|---|---|---|
| `M0` | 20.3 mm (10 cells) at 49→51 mm | 6.9 mm | 312 K | 234 kW/m² |
| `M1` | **112.3 mm (56 cells)** at 43→45 mm | **1.3 mm** | 357 K | 268 kW/m² |
| `M2` | 40.6 mm (20 cells) at 65→67 mm | 6.3 mm | 359 K | 269 kW/m² |

**The large geometric step is not a pathology here — it is what a cylinder looks
like in `z(r)`.** `M1`'s *median* step is 1.3 mm against `M0`'s 6.9: a flat floor
plus one wall, where `M0` is a cone. The thermal step that actually drives
D2o-0b's cascade is **comparable on all three legs** and below D2o-0b's
+402 K / 302 kW/m² on every one, with the D2k/D2r-1 side-face remedy on
(`f = 0.2`).

⚠ **This is nevertheless exactly where the 1 mm run must look.** D2t established
that `z(r)` does **not** converge in the taper band even where `D(z)` does, and
here the wall sits at r ≈ 44 mm, **inside the 40–56 mm w1 blend** — so the
wall's radius is set by the taper. Per the packet: if it steps, **widen the
taper, do not raise `h`.**

## 5. ⚠ What is NOT claimed

* **No Meier validation claim, and no agreement quoted tighter than ±5 mm.** The
  reference is a **lower bound** on diameter. `M1` landing within 0.2–4.4 mm of
  it is a **shape** result on a single mesh, not a validation.
* **No convergence claim.** One mesh. D2t showed the side-face lever's effect
  appears as a *mesh gap*, which a single mesh cannot see. This is a screening run.
* **A pass on (a) is "the funnel is removed down to the numerical floor", NOT
  "the model reproduces Meier's diameter".**
* **Volume (reported):** `M0` 0.97×, `M1` **0.82×**, `M2` 1.29× of Meier
  pro-rata. `M1` now removes *less* than Meier, where D2l found the model
  removing 1.9× too much. Reported, not scored.
* `side_face_factor = 0.2` remains **underived** and sets the wall's magnitude.

## 6. Deviations from `CRITERION.md`, disclosed

1. **Scored at 401.6 s, not the fixed 500 s.** `M1`/`M2` ended at ~402 s on the
   **frozen harness's own bottom watchdog** (`z − jet_s_c <= BOTTOM_M`), which is
   part of the d2i configuration `M0` must keep, not something added here — they
   drill deeper (362 mm vs `M0`'s 280 mm at the same time). 401.6 s is the
   **only** common time; it is forced by the data, not chosen. **Cost measured
   and reported: (a) on `M0` is 37.2 mm at *both* times (Δ = 0.00 mm)** — the
   shallow ladder (a) scores is set early and is insensitive to the window.
2. **The metric offset** already recorded in `PREFLIGHT.md` §1: my interpolated
   Ø(z) reads +0.3 to +3.3 mm above the packet's quoted table, so (a)'s `M0` is
   **37.2 mm, not 35.3**. One metric is applied to all three legs, so `M1 − M0`
   and `M2 − M1` are untouched.
3. **`M0` has no `side_P_face_above` column at all** (registered only when
   `side_face_factor_above >= 0`), so (c) records it as **absent**, not 0 — the
   same deviation D2r-1 recorded.
4. **`input_drilling`, not `input_feet_d2c`** — the packet's configuration table
   and its item 5 disagree; item 5 wins and the frozen argv is imported verbatim.
