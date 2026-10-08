# P1 — the gradient-aware firing criterion

`spall.driving_stress = absolute | differential` and
`spall.driving_ref = flake | fixed`, both default to the shipped behaviour.
Under `differential` the crack is driven on the stress **jump across the
candidate flake** rather than on the absolute superheat.

## The defect, demonstrated rather than asserted

`tests/MMWSpalling/unit/gradient_blind/` holds a column uniformly at 820 K with
**no net surface flux** — the flame gas sits at the rock's own temperature, so
the Robin flux is identically zero and no gradient ever forms.

| leg | firings | |
|---|---|---|
| `absolute` | **80** | the defect: a soaked slab fires exactly as a thermal skin does |
| `differential` + `flake` | **0** | |
| `differential` + `fixed` 0.5 / 1.0 / 2.0 mm | **0 / 0 / 0** | |
| `fixed` with no depth | parse abort | the depth has no default, by design |

A uniform slab has **zero** jump, not a small one: the two-cell reconstruction
of a uniform column correctly returns "no decay", so the reference temperature
*is* the surface temperature and the subtraction cancels identically.
**(b) PASSES** — and the packet notes nothing in this model had ever passed it.

## (c) THE DISCRIMINATOR — the flake scale is the physical length

`drive_frac`, the fraction of the absolute drive the differential form applies,
reads **0.330 = 1 − e^{−c_f} on both Kant legs, at every mesh and every flaw
size**. The two legs are physically very different — `R150` is slow and deep
(`ell_eff` 3.7–4.8 mm), `R500` fast and thin (0.9–1.3 mm) — and the flake
reference holds them together, which is what matches Kant's measured
390–560 K band.

**One refit factor puts both inside the band.** Measured sweep at 2 mm,
onset superheat [K]:

| `a0` [µm] | 8 | 40 | 60 | **80** | 100 |
|---|---|---|---|---|---|
| R150 | no fire | 522.3 | 487.0 | **461.3** | 441.1 |
| R500 | 650.4 | 526.8 | 492.1 | **465.8** | 441.6 |

**80 µm — a factor of 10 on the shipped 8 µm — is chosen on the 150 kW/m² leg,
where it reproduces the shipped onset temperature**: 463.3 K absolute vs
**461.76 K** differential (`S150_2.0`), −1.5 K. `kic.scale` stays 1.0 and nothing
is fitted against Meier.

⚠⚠ **Corrected (revision 5). "Where it fires on the anchor does not change"
holds on ONE leg of the pair, not both.** The sweep row above reads 465.8 for the
500 leg, but **465.8 is rev 7's `R500`, not a P1 leg at all.** Measured on the P1
binary at 2 mm: `S500` fires at **478.33 K** and `R500` at **482.28 K** — the
500 kW/m² leg moves **+12.5 to +16.5 K**. Both are still inside Kant's 390–560 K
band, so (c) is unaffected, but the refit is **not** onset-neutral and must not
be described as if it were.

⚠ `a0_gb` → **400 µm** at the shipped ratio 5. Below the packet's 0.5 mm stop
threshold, but **approaching grain scale** (Grimsel grains 1–3 mm) and reported
as a result rather than buried.

### `fixed` fails, and harder than pre-registered

The packet predicted `fixed` would fail R150 at ≈635 K. Measured: **there is no
flaw size at any reference depth that puts both legs in band.**

| `d_ref` | `a0` 20 µm | 40 µm | 80 µm | `drive_frac` 150 / 500 |
|---|---|---|---|---|
| 0.5 mm | no fire / 563.7 OUT | no fire / 486.3 IN | no fire / 389.4 OUT | 0.000 / 0.369–0.464 |
| 1.0 mm | no fire / 457.2 IN | no fire / 364.3 OUT | 583.5 OUT / 282.4 OUT | 0.164 / 0.777 |
| 2.0 mm | 597.9 OUT / 389.4 OUT | **493.1 IN** / 307.4 OUT | 316.4 OUT / 243.2 OUT | 0.388 / 0.943 |

The `drive_frac` column is the mechanism: a fixed length gives the two legs
drives differing by 2–4×, so a factor that fixes one throws the other out.
**P2 confirmed — this is the packet's main claim.**

## (c′) the anchor was validating a boundary condition we do not run

Neither Kant input sets `surface_patch.robin_form` and the harness never passed
it, so every shipped leg fell to the parser default `cell` while the Meier
production leg runs `surface`. Surface-form legs (`S150`, `S500`) are added; the
cell-form legs (`R150`, `R500`) are kept as the control, and `C150`/`C500` remain
the diagnostic controls.

| | surface form `S` (production) | cell form `R` (control) |
|---|---|---|
| 150 kW/m² firing time | 23.41 / 23.95 / 24.05 s, **2.68 %** | 18.93 / 21.60 / 22.86 s, **18.64 %** |
| 500 kW/m² firing time | 2.15 / 2.21 / 2.34 s, 8.40 % | 1.54 / 1.61 / 1.98 s, 25.61 % |
| 150 kW/m² `dT_fire` | **0.218 %** | 0.392 % |
| 500 kW/m² `dT_fire` | **0.301 %** | 0.974 % |

**(c′) PASSES at 150 kW/m²** (bar 5 %). The 500 kW/m² residual 8.40 % is the cold-start
estimator error the packet pre-registered as "cleanly isolated for the
successor", and it is halved from the cell form's 25.61 %.

⚠⚠ **(a) is therefore scored on the surface-form legs** — 0.218 % / 0.301 %,
both inside 0.5 %. That is the packet's own instruction in (c′) ("keep the
cell-form legs as the control"), not a choice made after seeing that the
cell-form R500 misses the bar at 0.974 %.

⚠ **There cannot be a surface-form (w) control.** `robin_form = surface` aborts
unless `spall.surface_T = reconstructed` (MMWSpalling.H:626) — which is exactly
the key the diagnostic control turns off. So (w) is necessarily a cell-form
comparison: it isolates the surface-T key, not the Robin form. Those legs were
written, the parser refused them, and they were withdrawn. (w) still passes:
the control's true surface temperature moves **30.29 %** at 150 kW/m² and
**46.67 %** at 500 kW/m² against the reconstruction's 0.39 / 0.97 %.

## ⚠⚠ (d) FAILS at 2.85 %, and the cause is the refit, not the criterion

| configuration | column rate spread |
|---|---|
| rev 7: `absolute`, `a0` = 8 µm | **1.51 %** |
| `absolute`, `a0` = 80 µm (the mandated refit) | **4.30 %** |
| **P1: `differential`, `a0` = 80 µm** | **2.85 %** |

⇒ **at the same flaw size the differential criterion IMPROVES mesh-independence
(4.30 → 2.85 %).** The bar was set against a measurement taken with a flaw ten
times smaller. An 80 µm flaw samples the near-surface profile where the shipped
`T_at_zv` lerp is mesh-dependent, and it does so for `K_Ic` at the crack tip as
well as for σ over the flaw.

**(d′) PASSES**: the level is **2.476 m/h against rev 7's 2.506, −1.2 %** — the
refit restored it rather than distorting it. The energy balance (c′) *improved*
to +0.06 / +0.17 / +0.12 % (rev 7: −1.06 / −0.35 / +0.71), and P2's level
against the continuum limit is now 1.214 / 1.215 / 1.215 against an analytic
1.213.

⚠ The fix would be to read `K_Ic`'s tip temperature from the exponential too.
That changes behaviour on **every** reconstructing leg, which is outside what
this packet asked for, so it is named for the successor rather than taken.

## ⚠⚠ Reading the two ends from different profiles does not work

The packet mandates the exponential for the **reference**. Leaving the crack
side on `T_at_zv` mixes two profiles, and that lambda's top-half-cell lerp is
mesh-dependent: an 80 µm flaw tip sits at **α = 0.08** of the (T_surf, T_top)
gap at `dz` = 2 mm and at **α = 0.32** at 0.5 mm. Measured cost: **2.4 % mesh
spread at R500**, which fell to 0.97 % once both ends read one profile. A
difference is only well defined on **one** profile.

## (a) identity half — the keys are inert when off

`tree_sweep_0922/witness/p1.md5`, 2639 files, **0 differing vs `p0b7.md5` and 0
differing vs `d2z1.md5`.** Both new keys default to the shipped behaviour and the
guard confirms the binary did not change across the sweep.

⚠ The monitor that announced this gate reported "0 files" because it sampled
`p1.md5` one second before the hash finished writing. **A watcher keyed on a
file's existence races the writer**; key it on the writer's exit instead.

## ⚠⚠ Meier: the run ABORTED at t = 177.6 s on a nozzle collision

```
ABORT: MMWSpalling.H:2214 (BuildFlameColumns) surface_patch nozzle collision:
stand-off s = -6.67e-06 <= 0 at column (i,j) = (6,11), time = 177.616
(z_nozzle = 0.475993, z_face = 0.476)
```

That column sits at r = √(13² + 23²) = **26.4 mm**, inside the Ø56 mm collision
circle. The nozzle descended 74 mm on its prescribed 1.5 m/h feed while that
column had given up only **24 mm** of rock, so the standoff ran to zero.

**The 60–160 s scoring window is complete**, so (e), (e′), (h) and (j) are all
measurable; (f) is not, and (g)'s bore half was never reached. The abort is not
a harness failure — it is the result.

### (e) PASSES, and by two orders of magnitude

Firings in the 60–160 s window, by radius, both runs truncated to t ≤ 177.616 s
so the comparison is like for like:

| band [mm] | rev7 (`absolute`, 8 µm) | **P1** (`differential`, 80 µm) | change |
|---|---|---|---|
| 0–20 | 15 742 | 13 850 | −12 % |
| 20–30 | 12 869 | 4 569 | −64 % |
| 30–40 | 10 259 | 1 253 | −88 % |
| **40–50** | 7 841 | **38** | **−99.5 %** |
| **50–60** | 2 573 | **0** | **−100 %** |
| 60+ | 395 | **0** | −100 % |
| **r ≥ 40** | **10 809** (21.8 %) | **38** (0.19 %) | **−99.6 %** |
| **r ≥ 50** | **2 968** (6.0 %) | **0** (0.00 %) | **−100 %** |

The bar was "fall by ≥ half". **Split at r = 50 mm as the packet requires**: the
50 mm+ half is a clean zero, and the 40–50 mm half is −99.5 %. ⚠ Per the packet
the 40–50 mm band is floor under floor-intensity heat, not side-heated wall, so
**those gains are artefact-assisted and are labelled so** — and (j) below now
puts a number on exactly how much of it is artefact.

### ⚠ (e′) FAILS at 94.7 % — but the denominator collapsed

Cap grip = per-firing `rc_flag` bit 16, the bar < 20 % beyond r = 40 mm.

| | rev7 | P1 |
|---|---|---|
| grip, all radii | **30.6 %** | **8.6 %** |
| grip, r ≥ 40 mm | 99.8 % | **94.7 %** |
| grip, r ≥ 50 mm | 100.0 % | *no firings* |

My estimator **reproduces the packet's stated 30.6 %** on rev7 to the digit
before being applied to P1, so the 8.6 % is trustworthy.

⚠ **(e′) fails as written and passes in substance.** Over the whole hole the
cap's grip fell **3.6×, from 30.6 % to 8.6 %**, well inside the 20 % bar — the
stopgap is largely retired. But in the wall band the gradient criterion does not
loosen the cap, it **deletes the firings**: 38 survive out of 10 809, and on
those 38 the cap still binds 94.7 % of the time. A criterion that empties the
population it was meant to relieve has not satisfied (e′); it has made (e′)
unmeasurable there. ⚠ **The grip no longer tightens with depth** (rev7's
18.7 → 38.6 % ladder), which was P0b's reason for calling the stopgap degrading.

### ⚠⚠ (f) WITHDRAWN (revision 5), and the shape error has reversed sign

The registered 150–300 mm section was never drilled — the run stalled at 186 mm.
What *is* comparable is the **collar**, which is converged: rev7's 10 mm and
25 mm diameters are **identical at 177.6 s and at 500 s**.

Diameter [mm] vs Meier Fig. 8.8 (±5 mm read precision, a LOWER bound on Ø):

| depth [mm] | Meier | rev7 @177.6 s | rev7 @500 s | **P1 @177.6 s** |
|---|---|---|---|---|
| 10 | 96 | 122.8 (+26.8) | 122.8 (+26.8) | **92.6 (−3.4)** |
| 25 | 92 | 110.9 (+18.9) | 110.9 (+18.9) | **86.2 (−5.8)** |
| 50 | 88 | 100.4 (+12.4) | 103.0 (+15.0) | 78.0 (−10.0) |
| 75 | 87 | 94.9 (+7.9) | 101.9 (+14.9) | 76.4 (−10.6) |
| 100 | 87 | 91.6 (+4.6) | 99.6 (+12.6) | 70.0 (−17.0) |
| 150 | 86 | 82.2 (−3.8) | 99.2 (+13.2) | 54.0 (−32.0) |

**At the collar — the one converged section — P1 cuts the error from +26.8 mm to
−3.4 mm, inside Fig. 8.8's own ±5 mm.** That is the first time the mouth has
come out right. ⚠⚠ But below ~50 mm the hole becomes a **cone**: P1 under-cuts
progressively where rev7 over-cut. The deep rows are not a fair test for either
run (that rock was drilled seconds earlier and has had no time to widen), so the
honest claim is the collar only.

Depth-normalised, which does not depend on run length — removed volume per mm of
axis advance, ×4 for the quarter domain:

| | Meier | rev7 | **P1** |
|---|---|---|---|
| volume / axis depth [mL/mm] | **5.93** | 6.7 (**+13 %**) | **3.1 (−47 %)** |
| window volume rate [m³/s] | 6.84e-6 | 7.17e-6 | **1.96e-6** |

⚠⚠ **Corrected (revision 5): Meier's figure is 5.93 mL/mm, not 7.1.** The
denominator is the **576 mm the burner descended** (1.5 m/h × 1383 s), not the
482 mm of finished hole. The check is self-consistency: 3.42 L over 576 mm
back-solves to an **86.9 mm** bore, which matches the reported 85–87 mm, whereas
482 mm would imply **95 mm**, which does not. ⇒ **rev7 was +13 % OVER Meier on
this measure, not −6 % under**, and P1 is **−47 %**, not −56 %. The sign of
rev7's error on removed volume flips with this correction.

⚠ **P1 still removes roughly half the rock per millimetre of hole that Meier
did.** The volume rate reproduces P0b's tabulated 1.794e-6 m³/s exactly as the
per-quadrant value (7.173/4), which validates the estimator.

The stall is **azimuthally ragged, not a ring-wide slowdown** — at r = 24–28 mm
the column depths run 24 / 98 / 160 mm (min/median/max). The colliding column
fired **58 times and voided only 12 cells**: whichever column gets ahead stays
ahead, and the laggards never catch up. Standoff: rev7 133.1 mm (above Fig. 8.8's
~100 mm ceiling), **P1 28.2 mm mean** and falling to zero — the right direction,
far past the target.

### (j) the two-cell fit under-reads the true flux by 2.7× in the wall band

Median `q_fit / q_real` and `drive_frac` in the 60–160 s window:

| band [mm] | `q_ratio` | `drive_frac` |
|---|---|---|
| 0–20 | **0.915** | 0.330 |
| 20–30 | 0.729 | 0.330 |
| 30–40 | **0.377** | 0.279 |
| 40–50 | **0.377** | 0.191 |

⚠⚠ **This is the quantitative confirmation of the packet's own warning.** On the
floor the two-cell exponential fit recovers 92 % of the real surface flux. In the
30–50 mm band it recovers **38 %** — and it is exactly there that `drive_frac`
drops below the self-similar 0.330 and the firings vanish. **A measurable
fraction of (e)'s suppression is the fit misreading the flux, not the physics
discriminating.** `drive_frac` = 0.330 on the floor is the self-similar value *by
construction*, so the floor is not being discriminated at all.

### (h) PASSES

| | bar | rev7 | **P1** |
|---|---|---|---|
| `(side_P_face + side_P_loss) − (rin + rout)` | ≤ 0.02 W | 0.0014 | **0.0100 W** |
| `(rin + rout) − (top + buried)` | ≤ 0.02 W | 0.0010 | **0.0099 W** |
| `ledger_err_solid` (Meier) | ≤ 5e-3 | 9.998e-5 | **1.054e-4** |
| `ledger_err_solid` (column, worst leg) | ≤ 5e-3 | — | **5.547e-4** (`S_2.0`) |

⚠ **P5 moved inboard**: `rin / rout` = 1898 / 128 W, **93.7 % inside 40 mm**
against rev7's 37.8 %. The narrow hole puts the burner's heat on the floor.
⚠ `ledger_P_loss` rose 135 → 1331 W because the un-removed wall stays hot.

### (g) the `fixed` arm — Kant half measured, bore half NOT RUN

The Kant half is in §(c): `fixed` has no flaw size at any `d_ref` ∈ {0.5, 1.0,
2.0} mm that holds both legs in band (drive fractions split 0.000/0.369,
0.164/0.777, 0.388/0.943 against `flake`'s 0.330/0.330). **The bore half was not
run.** With `flake` itself colliding at 177 s, three more 500 s Meier legs would
cost ~6 h to compare a form already refuted on the anchor. Stated as not run
rather than inferred.

---

# Revision 5 — the missing control, the withdrawals, and the census

Revision 4 was review-rejected. The criterion itself was **not** at issue and was
not re-done. Same binary `e485a62912d14bda…`, nothing rebuilt.

## 1. ⚠⚠ (d)'s attribution now has a run behind it — and it holds

Revision 4 excused (d)'s failure with *"at the same flaw size the differential
criterion improves mesh-independence (4.30 → 2.85 %)"*. **That 4.30 % had no
preserved run.** `p0_steady_column/test` hardcoded `driving_stress=differential`
and `sp.a0=80.0e-6`, `VARIANTS` had no flaw axis, and the harness `rmtree`s each
leg directory — so the number could not be recovered from anything on disk. Two
control variants (6 legs) were added and run:

| leg | configuration | rate at 2.0 / 1.0 / 0.5 mm [m/h] | spread |
|---|---|---|---|
| `A08` | `absolute`, **8 µm** (= P0b rev7 leg `S`) | 2.506 / 2.506 / 2.544 | **1.51 %** |
| `A80` | `absolute`, **80 µm** (the mandated refit) | 6.158 / 6.006 / 5.900 | **4.30 %** |
| `S` | `differential`, 80 µm (**P1**) | 2.476 / 2.433 / 2.406 | **2.85 %** |

* ⚠⚠ **`A08` reproduces P0b's archived 2.506 / 2.506 / 2.544 m/h to the digit**,
  spread 1.51 %, deviation **+0.1 %**. The harness has not moved since P0b, so
  the comparison is valid rather than merely asserted.
* **The attribution holds: `A80` 4.30 % > P1 2.85 %.** At the same flaw the
  `absolute` criterion is the *worse* of the two. The 4.30 % figure quoted
  without a run turns out to have been right, but it is only now evidence.
* ⚠ The flaw key here is the **deterministic scalar `spallation.sp.a0`**, not
  `weibull.a0_ig`, because `weibull.enabled = 0` on this input. `a0_ig` is the
  Kant/Meier key and is inert on this harness.

### ⚠⚠ An unasked-for finding: the refit and the criterion very nearly cancel

`A80`'s rate **level** is **6.0 m/h against `A08`'s 2.5 — the 80 µm flaw alone
multiplies the removal rate by 2.4.** The differential criterion pulls it back to
2.476, which is `A08`'s level to within **−1.2 %**. So the two halves of this
packet are a **matched pair**: enlarging the flaw makes firing much easier, and
applying only the 0.330 stress jump makes it correspondingly harder.

⇒ **(d′)'s "the refit restored the level" was not luck, and it was not neutral
either.** Neither change is separately defensible at these settings; they are
defensible only together, and **that coupling was never pre-registered.** The
successor should treat `a0` and the drive fraction as one calibration, not two.

## 2. (f) and (g)'s bore half are WITHDRAWN

Not "not scorable" — **withdrawn**. The converged 150–300 mm section was never
drilled (the run stalled at 186 mm) and **no P1 leg can reach it until the stall
is fixed**, so carrying them as unscored would imply a measurement that cannot be
taken. What is reportable, because it is converged:

| | P1 | Meier | rev 7 |
|---|---|---|---|
| Ø at 10 mm | **93.8 mm** | 96 | 124.8 |
| Ø at 25 mm | **88.2 mm** | 92 | 112.9 |

⚠⚠ **It must never be reported alone.** Beside it: the hole tapers to **72 mm at
100 mm** and **56 mm at 150 mm**, and **half the floor under the burner is
stranded** (§5). ⚠⚠ **No Meier rate, volume, width or standoff is quotable from
P1.**

## 3. (e′) restated — score the claim the treatment did not destroy

The banded form ("cap grip < 20 % beyond r = 40 mm") is **void**: the criterion
deletes 10 771 of the 10 809 firings there, so the denominator is gone. Scored
instead on the substantive claim:

| | bar | measured |
|---|---|---|
| **cap grip over ALL radii** | **< 20 %** | **30.6 % → 8.6 %** — **PASS**, a 3.6× retirement of the stopgap |
| (reason the banded form is void) | — | 38 survivors of 10 809; cap binds on 94.7 % of them |

⚠ **Deliverable 5's second witness has the same shape and is the same
observation, said once**: it asks for `drive_frac` ≪ 0.1 on the wall band, the
lowest band median is **0.191** at 40–50 mm, and nothing beyond 50 mm fires at
all. Both are one fact — *the treatment removed the population it was to be
measured on* — not two independent misses.

## 4. Corrections folded in

| | revision 4 wrote | measured (2 mm, P1 binary) |
|---|---|---|
| (c) `S150` onset | 461.3 K | **461.76 K** (reproduces rev 7's 463.3, −1.5 K) |
| (c) `S500` onset | "465.8 K" | ⚠⚠ **478.33 K**. 465.8 is **rev 7's `R500`**, not a P1 leg; P1's `R500` is **482.28**. ⇒ the 500 leg moves **+12.5 to +16.5 K** — onset-neutrality holds on **one leg, not both** |
| (c′) other half | absent | **`S500` 8.40 %**, above the 5 % bar — pre-registered as cold-start estimator residual, but it belongs in the scorecard |
| (h) column half | 1.054e-4 (Meier) | worst **column** leg **5.547e-4** (`S_2.0`), still far inside 5e-3 |
| (j) reference | "could not reproduce 0.78/0.53/0.26/0.10/0.20" | ✅ those were a **planner estimate**, never a measurement — `q_ratio` is this packet's own new column. The correction travels with the numbers |
| Meier volume/mm | 7.1 mL/mm | **5.93** (3.42 L over the **576 mm the burner descended**) ⇒ rev 7 **+13 %**, not −6 %; P1 **−47 %**, not −56 % |

## 5. ⚠⚠ The exposed-face census — the collision is a ROUGHNESS instability

Reduced from the shipped event log by `census.py` (no new run). t ≤ 177.616 s,
r < 40 mm, "frozen" = last fired before t = 120 s.

| | P1 | rev 7 |
|---|---|---|
| columns under the burner that fired at all | 316 | — |
| ... that **stopped firing** before t = 120 s | **152 (48 %)** | **0** |
| median depth, frozen vs live | **48 mm vs 180 mm** | — |
| frozen columns with **≥ 2 exposed sides** | **78 %** (1:22, 2:31, 3:26, 4:21 %) | — |
| live columns with **0 exposed sides** | ⚠ **51 %** (see below) | — |
| neighbour-to-neighbour step at a frozen column | median **22 mm**, p90 **76 mm** (cell = 2 mm) | — |
| **the colliding column (6,11)**, r = 26.4 mm | **4 of 4 sides exposed**, 24 mm deep, **122 mm below its deepest neighbour**, last fired **19.1 s**, stood **159 s** | — |

**Every row reproduces the planner's census except one.** `r < 40 mm` with
`t < 120 s` is the unique cut that returns both 152/316 and 78 %, and at that
same cut the live-column figure is **51 %, not the 69 % stated**. Six exposure
definitions and a sweep over the radius and freeze thresholds were tried; none
returns 69 % while also returning the two headline rows. ⚠ **The correction
strengthens the finding**: even the *live* columns stand proud of a neighbour
about half the time, so the roughness is not confined to the frozen set.

**The mechanism.** A column that fires sheds its skin, exposes cool rock, gets a
fresh steep profile and fires again; a column that lags soaks heat to depth,
flattens, its drive collapses, and it stops firing and lags further. Because the
differential criterion keys on the *gradient*, flatness is self-reinforcing —
**this feedback is introduced by the criterion.** The frozen region is **not a
plateau**: median neighbour step 22 mm with a p90 of 76 mm, against a 2 mm cell,
is a field of towers and chasms.

⇒ ⚠⚠ **The collision is at r = 26.4 mm — INSIDE the resolved floor, not in the
wall band P1 was aimed at. It is a new failure mode, not the one P1 set out to
fix**, and it needs its own unit test in the successor. The axis rate is almost
unchanged (6.30/5.00/3.83/3.19/2.93/2.70/2.48 vs rev 7's
6.57/5.10/4.09/3.49/3.08/2.94/2.75 m/h) while the minimum standoff collapses
56 → 0.5 mm against rev 7's 59 → 164 mm. The gap is a **minimum over columns**,
so it tracks the frozen ones. ⇒ **a uniformity failure, not a rate failure.**
