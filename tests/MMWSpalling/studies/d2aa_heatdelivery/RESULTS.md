# D2aa RESULTS — heat delivery under the burner body

`CRITERION.md` sha `d81932d18019649bc4e7901855aed6d92c91e53e5178a2a4d9a27042a0e494e3`,
hashed with **zero outputs on disk**. Full scoring in `analysis.txt`.

**Headline: the `WL` leg discriminated. `A`'s 1.251 m/h WAS the step drain, so
removing wall heat un-throttles the ring to ~2.1 m/h, not to Meier's band. The
wall lever is real but smaller than predicted, the pads leg JAMMED exactly as
pre-registered, and specific energy is invariant — so heat cannot be moved
between regions to buy rate.**

## Scorecard

| criterion | result |
|---|---|
| (a) rate, shape, volume | reported below; `A` at 2.27 mm still the best shape on record |
| (b) false-pass test | **not triggered** — no leg lands in the 1.3–1.6 m/h band |
| (c) power identity | **PASS on every leg** — 6.0e-3 to 9.9e-3 W against 0.02 |
| (d) region split closes | **PASS** — 1e-8 to 5e-8 relative, bar 1e-4 |
| (e) identity | **PASS** — `S0` byte-identical to `RB_D`, header delta EMPTY; tree **2639/2639** vs `d2z1.md5` |
| (f) specific energy | **1.31–1.33 J/mm³ on every leg** |
| (g) item 0 audits | reported below; **the packet's decay-clip figure is corrected** |

| prediction | verdict |
|---|---|
| **P1** (reviewer, PRIMARY) `WL` 1.8–2.2 m/h | **HELD** — 2.073 |
| **P2** (planner alt) 1.3–1.6 m/h | **REFUTED** |
| **P3** side path 0.3–0.7 kW | **REFUTED** — 1.32 kW; 0.83 at the band's other end |
| P3b `PD` ring power ≥ 5× down | held by construction, **VACUOUS** |
| P4 `BO` below both single edges | numerically held, **contaminated by the jam** |
| P5 false-pass guard | **VACUOUS** — no leg in band |
| **P6** > half the charged side AREA under the pads | **REFUTED** — 49.5 % (58.8 % of the POWER) |
| **P7** pads-leg jam signature | **CONFIRMED** |

---

## 1. The discrimination — P1 held, P2 refuted

| leg | power to rock | rate | mean \|ΔØ\| | Ø@100 | Ø@125 | pit lead | vol |
|---|---|---|---|---|---|---|---|
| `SW` frozen | 8.05 kW | 3.694 | 2.80 | 87.6 | 85.8 | 38 | 1.08 |
| **`WL`** h=142 | **5.48 kW** | **2.073** | 3.81 | 86.7 | 83.6 | 88 | 1.23 |
| **`WL81`** h=81 | **4.94 kW** | **1.780** | 2.73 | 85.2 | 82.7 | 98 | 1.22 |
| `PD` shielded | 5.09 kW | **jammed** | — | — | — | — | — |
| `BO` both | 5.43 kW | **0.411** | — | — | — | — | — |
| `A` no side path | 4.18 kW | 1.251 | 2.27 | 83.8 | 81.3 | 128 | 1.27 |
| `D` f = 0.2 | 10.35 kW | 3.799 | 7.56 | 94.2 | 90.6 | 32 | 1.16 |

`WL` lands at **2.073 m/h** — inside P1's 1.8–2.2 and on D2x's un-throttled
pad-ring rate of **2.09–2.11**. `WL81` reads **1.780 m/h** with Ø@100 = **85.2**,
which is exactly P1's predicted 84–85 wall band.

**⇒ `A`'s 1.251 m/h was the step drain.** The planner had already withdrawn the
premise P2 rested on; the measurement confirms the withdrawal. Neither wall-`h`
leg reaches Meier's band, and both pass through it only if the power falls
further.

Ring table, 2 mm rings, **per ring, never averaged** [mm recession]:

| leg | 30 | 32 | 34 | 36 | 38 | 40 | 42 | 44 | 46 | 48 |
|---|---|---|---|---|---|---|---|---|---|---|
| `WL` | 167.6 | 166.7 | 161.4 | 155.7 | 147.7 | 132.2 | 66.1 | 20.6 | 5.2 | 0.4 |
| `WL81` | 151.1 | 149.9 | 142.4 | 135.6 | 124.4 | 106.7 | 48.1 | 14.2 | 3.4 | 0.1 |
| `PD` | **0.0** | **0.0** | **0.0** | **0.0** | **0.0** | 49.7 | 50.0 | 50.6 | 50.8 | 51.5 |
| `BO` | 72.3 | 74.7 | 77.3 | 80.8 | 82.6 | 67.3 | 66.7 | 64.5 | 57.8 | 49.2 |
| `SW` | 256.3 | 259.4 | 258.4 | 258.1 | 258.3 | 247.8 | 62.1 | 1.7 | 0.0 | 0.0 |
| `A` | 124.9 | 119.2 | 112.7 | 104.3 | 90.1 | 69.1 | 37.9 | 11.7 | 3.4 | 0.2 |

⚠ **No Meier agreement claim tighter than ±5 mm** — Fig. 8.8 is a visible width,
i.e. a LOWER BOUND on diameter; band 85–95 ± 5 mm.
⚠ **Shape:** `A` at 2.27 mm is still the best on record. `WL81` at 2.73 is the
closest any side-path leg has come; `WL` at 3.81 is worse than `SW`'s 2.80. The
refuting-evidence note stands: the side-face path still costs shape.
⚠ Pit leads 88 and 98 mm are inside the registered 150 mm ceiling, so neither
`WL` leg is disqualified.

## 2. P3 REFUTED at BOTH ends of the band

Side path, full-hole, 150–400 s:

| leg | total | bore | pads | wall |
|---|---|---|---|---|
| `SW` | 3540.9 W | — | — | — |
| `WL` h=142 | **1321.6** | 183.5 | 669.6 | 468.6 |
| `WL81` h=81 | **825.5** | 133.7 | 416.2 | 275.6 |
| `PD` | 2243.7 | 623.2 | **−53.4** | 1674.0 |
| `BO` | 2455.7 | 210.3 | 1204.1 | 1041.3 |

The predicted **7–9x** cut is really **2.7x** at h = 142 and **4.3x** at h = 81.
Naive scaling (`142/685 x 3.54 = 0.73 kW`) assumes a fixed driving ΔT; it is not
fixed, because less heating leaves the wall cooler and raises `T_gas − T_wall`,
clawing part of the reduction back. **Running both ends of the band is what made
that a mechanism rather than a missed number** — and it was the reason for
running both, since choosing h = 81 alone would have landed P3 and would have
been choosing `h` to hit a rate.

⚠ **`PD`'s `side_P_pads` is NEGATIVE (−53.4 W).** With the ring shielded the gas
deposits nothing there but the faces still radiate. "Adiabatic" is adiabatic to
the gas only; it is not adiabatic.

## 3. P7 CONFIRMED — the pads leg jammed, and that is a result

`PD`: rings 30–38 mm at **0.0 mm recession**, never reached 150 mm descent,
`wall_step_f_multi_pads` pinned at **1.0000**. That is the pre-registered
signature verbatim: **the model's lateral-recession limit** — undercut is
impossible in this code, so a shielded support ring simply stops receding and the
burner stops — **not physics**, and to be read as D2o-0's 0.11 m/h jam.

`BO` at **0.411 m/h** is in the same territory, so P4 "holds" numerically while
measuring the jam rather than the enthalpy coupling it was written for. Its
fitted power-rate exponent is **5.59** against 1.50 for both `WL` legs — not a
physical exponent, but what burning near-`WL` power and going nowhere looks like.

**⇒ The pads lever cannot be tested in this code as written.** It is blocked on
the same vertical-recession limitation the line has carried since D2p, not on a
missing coefficient.

## 4. Specific energy is invariant — the constraint this packet establishes

| leg | P to rock | vol | J/mm³ |
|---|---|---|---|
| `WL` | 5.48 kW | 1.23 | **1.32** |
| `WL81` | 4.94 | 1.22 | **1.31** |
| `SW` | 8.05 | 1.08 | **1.33** |
| `A` | 4.18 | 1.27 | **1.33** |

**The ratio does not move.** It is the sensible heat to the pinned surface
temperature, as energy-set recession requires. **So heat cannot be redistributed
between bore, pads and wall to buy rate: the rate follows the POWER.**

The packet's own `rate ∝ P^1.59` predicts every leg within ±4.5 % (fitted
exponent 1.50 on both `WL` legs, 1.65 on `A`), and inverting it reproduces the
packet's target:

| target | power needed |
|---|---|
| 1.30 m/h | 4.17 kW |
| 1.45 m/h | 4.47 kW |
| 1.60 m/h | 4.75 kW |

`WL81` sits at **4.94 kW** ⇒ the residual is **~0.2–0.5 kW**, now a number rather
than a hypothesis.

⚠ If rock really takes 4–6 J/mm³ the model's spall energetics are wrong by ~3x —
a larger, **separate** finding, testable against Wilkinson and Tester's measured
surface temperatures. **Not a licence to keep the power**, and not acted on here.

## 5. P6 — refuted on area, held on power

Charged side area, cm² full-hole: `WL` bore 18.4 / pads 65.7 / wall 48.7 ⇒ pads
**49.5 %**, just under half. On **power** pads (669.6 W) vs wall (468.6 W) is
**58.8 %**. Pads is the largest single band either way, so the packet's conclusion
— `WL` alone cannot account for the 3.54 kW — stands; only the ">half by area"
form of it fails.

⚠ **Item 0(a)'s `bore 1 / pads 10 / wall 7 cm²` per quarter could NOT be
reproduced**, because it is a claim about `SW`, which is frozen and has no region
columns. On `WL` the split is 4.6 / 16.4 / 12.2 per quarter — same **ordering**
and a similar pads:wall ratio (1.35 against the packet's 1.43), but ~1.8x larger
in total on a leg with different physics. A bridge leg (SW's argv + the new
columns) would close this; it was not run.

## 6. Item 0 audits (g)

**(a) The refuted-premise table reproduces EXACTLY** on `SW`'s own log, 150–400 s:
`side_faces` 6030.7 (packet 6031), `side_cols_above` 5081.0 (5081, 84.3 %),
`delta_faces` 949.7 (950), `side_P_face` 885.2 W quarter (885.2). The identity
`side_faces − side_cols_above = delta_faces` holds exactly. Charged area after the
wetted law is 15.7 cm² (packet "≈16–18") ⇒ **564 kW/m²** (packet "≈560").
**The packet's central premise — the wall lever is `h`, there is no area excess —
stands.**

**(c) THE PACKET'S DECAY-CLIP FIGURE IS WRONG, and the correction was registered
before any leg ran.** `SW` and `A` both set `jet_core_length = 8.0`,
`jet_D = 0.0075`, so the core length is **8 D = 60 mm, not "twelve nozzle
diameters = 90 mm"**. It is not occasionally in play: the stand-off exceeds it in
**100 % of rows on every leg** (`A` 157.3 mm mean, `SW` 87.8, `D` 81.5), so the
decay is always engaged and never clipped at phi = 1. **But it is not a
discriminator**: under `jet_T_ent_mode = exhaust` the mixing temperature is
already 1463 K, and measured `T_stag` is **1776–1780 K across all three legs**
despite `s_c` spanning 81–157 mm. The pit-lead ceiling was registered as asked and
is not load-bearing.

**(b) The D2l factor-of-four audit: NOT DONE.** It requires re-deriving D2l's
absorbed-power fraction from its own logs and was displaced by the parse defects
in §7. Carried as debt; it is a reporting item and changes nothing here.

## 7. What went wrong, and what it cost

**Three parse defects, all mine, all the same shape.** The pad-ring guard had
**three** wrong homes before it worked:

1. beside `pad_shield_open`'s own range check it ran **before** the feet branch
   reads `foot_r_inner/outer`, so it aborted exactly the legs that DID set the
   ring — `PD` and `BO` died at 0.1 min;
2. moved inside the feet branch it stopped firing at all, leaving the key
   **silently inert** without feet — worse than the abort;
3. at the end of `Parse` it was unreachable, because the expression path returns
   early at the `use_expr = true` exit — the path every study leg takes.

It is now `CheckPadRing`, called on **both** exits. Separately,
`pad_shield_open` read `pp.contains` **after** `query_default` had added the key —
the trap documented two lines above it in the same file.

**One guardrail broken.** Three rebuilds happened while the first `S0` was
running. The project rule says a mid-sweep rebuild **invalidates the identity
leg**, so that `S0` was killed and re-run clean against the final binary
(`27beeeae…`). The edits in question touched only the pad-shield guard, which
`S0` does not set — but byte-identity is the only oracle this line has, and that
argument is not a substitute for the re-run.
