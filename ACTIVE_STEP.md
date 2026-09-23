# Active Step: D2v - The body-in-the-hole h map: separate the funnel from the wall

Status: ready for implementation

**The dominant shape error is a FUNNEL, not a narrow hole.** Measured on frozen
runs with the interpolated metric: the model cuts **141 mm at 10 mm depth against
Meier's 96** — ~45 mm too wide at the mouth — and is about right at 150–200 mm.
The cause is located: the impinging wall-jet `h_expr` **never falls off in
radius**, so the whole block top is excavated at 160–360 W/m²K while `s > 0`.

This packet gives the burner a **radially confined h map** — impingement under
the jet, annulus on the wall, nothing beyond the mouth — and attributes the
result in three legs so the **funnel fix** and the **wall heating** are separable.

D2u is withdrawn before launch (scalar above-plane `h` would have blanketed 89 %
of columns; its P2 was refuted by data on disk). Record in `ARCHIVE_DONE.md`.

## Context

**The error, measured** (`d2i_gate/LI0_2mm` at 572 s; `d2m_feet/R40` agrees at
250 s):

| depth [mm] | 10 | 25 | 50 | 75 | 100 | 150 | 200 |
|---|---|---|---|---|---|---|---|
| **model** | **141** | 128 | 113 | 105 | 99 | 92 | 86 |
| Meier, digitised | 96 | 92 | 88 | 87 | 87 | 86 | 85 |

**The cause.** `input_feet:163`'s `h_expr` evaluated directly:

| r [mm] | 5–18 | 25 | 30 | **40** | 50 | **76** | 120 | 168 |
|---|---|---|---|---|---|---|---|---|
| h [W/m²K] at s = 50 mm | 1449 | 1092 | 912 | **686** | 549 | **362** | 229 | **164** |

It decays as ~1/r and never reaches a floor. `surface_patch.radius = 0.2` and the
quarter domain's corner is at r = 168 mm, so **every column is in the patch** and
every column with `s > 0` gets this. Physically r > 40 mm is under the skirt and
in the exhaust annulus, **not** under a wall jet (impingement ~700–1500 vs
grazing ~40 W/m²K).

**The collar is made early, below the plane, then freezes** — because those
columns pass above the plane where `h` is hard-zeroed:

| t [s] | 100 | 200 | 300 | 340 | 450 | 572 |
|---|---|---|---|---|---|---|
| Ø @ 10 mm | 126 | **141** | 141 | 141 | 141 | 141 |
| Ø @ 50 mm | 67 | 106 | **113** | 113 | 113 | 113 |
| mouth ring r = 50–76 mm, recession [mm] | 11.0 | 27.9 | 33.3 | **33.6** | 33.6 | 33.6 |

**This is the floor/wall split D2o-0 attempted.** It failed twice — on a sharp
radial corner, then on a step-face cascade. **Both failure modes now have
remedies:** D2t's smooth taper (>= 8 coarse cells) and D2k/D2r-1's side-face
heating. That is why it is worth attempting again, and why the taper width and
the side-face keys are not optional.

### ⚠ PRE-FLIGHT ITEM 0 — THE STANDING RULE, six packets bought it

**Before writing code: name the call site of every new or changed behaviour and
demonstrate, from the literal key set of each leg, that it is reached.** Record in
`PREFLIGHT.md`. **If any leg cannot reach it, STOP and report.**

### ⚠ THE BLANKET — the way this packet fails catastrophically

Measured at 572 s (`nozzle_z` 214 mm, `foot_z` 164 mm): **3210 of 3600 columns
are above the plane (89 %), 2294 of them never touched.** Giving those a
*uniform* `h` at ~1690 K delivers ~110 kW/m² against ~24 kW/m² of losses and
fires undisturbed block top — **D2o-0's 1400 K blanket again.** The above-plane
state **must be the same radially confined expression**, never a scalar.
Criterion (d) is the witness and it can stop the packet.

## Sources — every line and number opened this session

| site | what |
|---|---|
| `MMWSpalling.H:1773` | `hj = (s > 0.0) ? surface_patch.h_f(...) : 0.0;` — **above-plane `h` hard-zeroed under `jet_enthalpy`; `h_expr` is not consulted.** The line to change. |
| `:1781` | `flame_Tg[col] = surface_patch.jet_T_ent;` — **293.15 K** above the plane. |
| `:1776` | aborts unless `hj` is finite and `>= 0` (enthalpy branch). |
| `:1653-1656` | side-face term and the `side_f_above` selector. **Sentinel < 0 falls back to `side_face_factor`.** |
| `:2754-2885` | the **ledger twin** — mirror any side-face change or `ledger_err` breaks. |
| `:622`, `:639` | `jet_T_exhaust`, `jet_T_rec` are registered ⇒ in `thermo.dat`. |
| `:2691-2693` | losses are gated on `surf_c`, **not** on the patch. |

**Frozen runs, import by path, never re-run:** `studies/d2i_gate/output/LI0_2mm`
(0.4 m domain, `timestep 0.016`, 700 s), `studies/d2m_feet/output/R40`,
`studies/d2t_radial_taper/`, `studies/d2r1_aboveplane_sides/`.

### The Meier configuration — use `input_feet_d2c`, NOT `input_feet`

| key | `input_feet` (STALE) | `input_feet_d2c` (USE THIS) |
|---|---|---|
| `geometry.prob_hi` | `0.12 0.12 0.3` | **`0.12 0.12 0.4`** |
| `timestep` | `0.004` | **`0.004` -> set `0.016`** (proven at 2 mm; the frozen D2i/D2m runs use it) |
| `stop_time` | 340 | set **500** |

`amr.n_cell 60 60 200` ⇒ **dx = dy = dz = 2 mm**. `surface_patch.radius = 0.2`,
`foot_r_outer = 0.040`, `jet_closure = enthalpy`, `jet_T_ent = 293.15`,
`jet_T_ent_mode = exhaust`, `nozzle_descent = feet`. **`side_face` appears ZERO
times** ⇒ `side_face_flux` defaults off.

⚠ **The 4 ms timestep makes any runtime estimate ~4x too low.** Set 16 ms.

### The scoring reference — `validation/meier/meier_fig8_8_hole_profile.csv`

| depth [mm] | 10 | 25 | 50 | 75 | 100 | 125 | 150 | 175 | 200 | 250 |
|---|---|---|---|---|---|---|---|---|---|---|
| visible width [mm] | 96 | 92 | 88 | 87 | 87 | 87 | 86 | 85 | 85 | 84 |

⚠ **A VISIBLE WIDTH ⇒ a LOWER BOUND on Ø**, at **±5 mm** (the saw cut need not
pass through the axis). Thesis D = 85 mm (Tab. 8.2, p. 222); **3.42 L / 0.5 m
implies a mean of 93 mm**. ⇒ **score the magnitude against a band of 85–95 mm
±5, and use the csv for SHAPE only.** Never quote agreement tighter than ±5 mm.

### The gas temperature above the plane

Measured on `LI0_2mm` at 572 s: **`jet_T_exhaust` = 1690 K, `jet_T_rec` = 1689 K**
— **uncooled**. Meier's ~186 W/K water jacket is **not in the model**.
⇒ **declare `jet_T_above` with the jacket drop stated as an assumption**, and
report the run's own exhaust range beside it. **Do not adjust it to fit a result.**

## Goal

### 1. The h map — key-only, ONE expression that BRANCHES ON `s`

⚠ **The map must NOT be the same blend above and below the plane.** **REJECTED
FORM, computed for the discarded 36 -> 52 mm ramp** (kept because it is the
argument, not the design): the impingement term is clamped to its 15 mm stand-off
value at the wall, so above-plane columns at r = 42–47 mm would receive

| r [mm] | 40 | 42 | 45 | 47 | 50 | 52 |
|---|---|---|---|---|---|---|
| blended h [W/m²K] | 557 | 467 | 350 | 282 | 193 | 142 |
| q vs 1690 K at a pinned 822 K wall [kW/m²] | **483** | **405** | **303** | **245** | 168 | 123 |

against a floor `q_pin` of ~520 kW/m². **That is ~1 m/h of wall recession above
the plane and would drive the hole out to the hand-over radius regardless of
physics — the packet would decide its own answer.** The **adopted** 40 -> 56 mm
ramp is no better above the plane (w1 = 0.625 at r = 46 mm gives ~430 W/m²K), which
is why the fix is the `s` branch, **not** a change of radii.

**Structure (key-only, `if(s>0.0, …, …)` as D2j-0b used):**

```
  floor branch  (s > 0) :  h_floor(r,s) = h_imp(r,s)*w1(r) + h_ann*(1 - w1(r))
  wall  branch  (s <= 0):  h_wall(r)    = h_ann
  both then x the mouth taper:  h = branch * w2(r) + h_out*(1 - w2(r))
```

- `h_imp` = the **existing** `input_feet:163` expression, unchanged.
- `h_ann` = **142 W/m²K**. Defensible band **81** Dittus-Boelter / **142** area
  expansion / **250** entry length at **Re = 5568**. ⚠ Provenance: **D2r-0
  pre-flight arithmetic, not any run**, and that pre-flight is where a *laminar*
  correlation was wrongly used on a turbulent flow — quote the corrected band only.
- `h_out` = **1.0e-3** (outside the mouth, vented by Meier's sealed wellhead; the
  parser forbids `h <= 0`).
- **`w1`: 1 for `r <= r_imp = 0.040`, smooth monotone to 0 by `r_ann_r = 0.056`.**
  ⚠ **`r_imp` MUST be the skirt radius `foot_r_outer = 0.040`, not less.** The pads
  sit at 28–40 mm and the burner rides the highest rock under them; at
  `r_imp = 0.036` the outer pad rock loses **10 % at 38 mm and 20 % at 40 mm**,
  recedes slower and slows the feed — **D2o-0's jam mechanism in mild form**.
- `w2`: 1 out to `r_mouth = 0.060`, smooth to 0 by `r_out = 0.076`.
- **Both blends >= 8 coarse cells (>= 16 mm at dz = 2 mm).** D2o-0b failed on a
  step; D2t passed on 8 cells.
- **No `min`/`max` in the blends** (numeric-first rule) and **never
  `a/max(x,num)`** (AMReX parser bug).

### ⚠ 1b. PRE-REGISTER THE NUMERICAL FLOOR THIS CREATES

A 16 mm blend from 40 to 56 mm hands over at a centre of **48 mm**, so **this
configuration cannot cut narrower than Ø ~ 96–100 mm at 2 mm.** Meier reads
**96 at 10 mm depth but 88 at 50 mm and 85 at 200 mm**. ⇒ **the model will remain
wider than Meier at depth BY CONSTRUCTION, and a pass on (a) must not be read as
"the model reproduces Meier's diameter".** **Taper width is a NUMERICAL
parameter** — to be tested at 8 mm on a 1 mm mesh later, **never tuned to move a
result.**

### 2. The source change — evaluate the SAME expression above the plane

At `:1773`/`:1781`, gated on two new `surface_patch` keys, **both defaulting to
sentinels so unset runs are byte-identical**:

- `jet_above_h` (bool, default 0): when 1, `hj = h_f(x, y, r, s, time)` for **all**
  `s` — i.e. the confined map applies above the plane too.
- `jet_T_above` (scalar, default −1): when `>= 0`,
  `flame_Tg[col] = (s > 0.0) ? jet_T_ent : jet_T_above`.

**The parser MUST abort if one is set without the other, and if
`jet_T_above <= jet_T_ent`.** Above-plane `flame_Tg` is 293.15 K, so `h` without a
hot `T_gas` **refrigerates the wall** — a code invariant, not a guardrail.

⚠ `h_imp` clamps `s` to `[0.015, 0.09]`, so above the plane it would return the
**strongest** near-field profile. **That is exactly why item 1 puts a flat
`h_ann` on the wall branch instead of a blend.** Do **not** reinstate a blended
above-plane map — see the rejected-form table in item 1.

### 3. FIRST, AND IT CAN STOP THE PACKET — trace the COMPILED map, BOTH branches

⚠ `flame_h` is a plain `std::vector` (`:5825`), **never written to a plotfile,
`thermo.dat` or any csv**, and a numpy re-evaluation is only a filter because the
D2o-0 bug lived **inside the AMReX parser**. Use D2t's method: **a few-second
2 mm run stopped before any spall**, reading `temp` from the first plotfile and
binning by radius.

**⚠ TWO CORRECTIONS TO D2t's RECIPE, both verified:**

1. **The march makes `T_gas` radius-dependent**, so the raw rise is not `h(r)`.
   `LI0_2mm_jet_profile.csv` carries `time, bin, r_lo, n_cols, T_gas, …` — **per
   bin `T_gas` and `r_lo`, but NO `h` column.** So **normalise the measured rise
   by `(T_gas(r) − T_0)` read from that csv.** (`T_gas` starts near **1207 K** and
   varies with `r`.)
2. **⚠ A FEET-DESCENT RUN NEVER EXECUTES THE ABOVE-PLANE BRANCH.** The nozzle
   rides 50 mm above the feet, so at t = 0 **every column has `s = +50 mm > 0`** —
   all below the plane. A few-second run therefore **cannot test the wall branch
   at all.** ⇒ **a SECOND trace leg is required**: `nozzle_descent = prescribed`
   with the plane already **below** the original surface, the same expression,
   `jet_closure = enthalpy`, and the collision guard off.

**The prescribed-descent leg's key set, spelled out.** Start from the frozen
`d2i_gate` set and **strip every feet-only key**, then add the prescribed pair:

| strip | why |
|---|---|
| `nozzle_descent = feet` -> **`prescribed`** | the mode switch |
| `foot_body_clearance` | **aborts** without feet (`:6287`) |
| `foot_r_inner`, `foot_r_outer`, `foot_standoff`, `foot_rule`, `foot_npads`, `foot_pad_quantile` | `query_required` **only** under feet (`:6317`), so ParmParse strict-checks them as unused (CLAUDE.md lesson 3) |
| `nozzle_collision_radius`, `nozzle_feed_max` | feet-path guards |
| **`jet_free_surface = 1`** | ⚠ see below |
| add `nozzle_z0`, `nozzle_feed` | **forbidden with feet** (`:6312`), **required** without |

⚠ **`jet_free_surface` does NOT abort under prescribed descent — it silently
misbehaves, which is worse.** The free-surface test is
`r > surface_patch.foot_r_outer` (`:2223`), and `foot_r_outer` **defaults to 0.0**
(`:5962`) and is never parsed outside the feet branch. ⇒ **every column would be
treated as free surface.** **Strip it explicitly**; do not rely on an abort.

⚠ **Expect a HIGHER exhaust temperature on `M1`/`M2` than the frozen 1690 K.**
The confined map extracts less power from the block top, so more enthalpy stays
in the gas. **The declared `jet_T_above` will therefore sit on the LOW side —
state the direction and do NOT adjust it to match.**

Require, across both legs: monotone decreasing; `h_imp` inside `r_imp`;
**`h_ann` on the WALL branch at every `r` inside `r_mouth`**; the floor by
`r_out` with **no rise anywhere beyond it**. Hash as its **own** pre-flight,
separate from `CRITERION.md`. **If either trace is not the intended map, STOP.**

### 4. The metric — interpolate, and reuse D2t's analyzer

Final `k_top` per column from `_removal_events.csv` -> `z_face`; bin on **fixed
physical 2 mm radial bins**; **linearly interpolate the crossing radius** between
bracketing bins for each depth. ⚠ D2t's `D = 2*r_centre` snapped to a 4 mm
lattice — **4 mm is 4.7 % of an 85 mm hole and equals the reference's own ±5 mm**,
the whole discrepancy budget. Score Ø(z) and volume. **Never rate** (D2l: rate
agreement was a coincidence of two errors).

### 5. The campaign — three legs, 2 mm, 500 s

| leg | literal keys |
|---|---|
| `M0` | **the frozen `d2i_gate` key set, imported from its harness** (`pinned_idle_cycles=1.0e9`, `foot_body_clearance=0`), `stop_time=500.0` |
| `M1` | `M0` + the new `h_expr` (item 1) + `surface_patch.side_face_flux=1` + `surface_patch.side_face_factor=0.2` + **`surface_patch.side_face_factor_above=0.0`** |
| `M2` | `M1` but `surface_patch.side_face_factor_above=1.0` + `surface_patch.jet_above_h=1` + `surface_patch.jet_T_above=<declared>` |

**`M1 − M0` is the funnel fix. `M2 − M1` is the wall heating.** That separation is
the packet's purpose.

⚠ **`M0` MUST BE BUILT FROM THE FROZEN RUN'S EXACT KEYS, NOT FROM
`input_feet_d2c`.** They differ in **two** keys — `pinned_idle_cycles` **2.0 vs
1.0e9** and `foot_body_clearance` **1 vs 0** — so an `input_feet_d2c` baseline
would miss P1's bar for **configuration** reasons and the miss would be
misattributed. Import the `d2i_gate` harness key set; `M1`/`M2` then differ from
it by the listed keys **only**.

⚠ **Above-plane TOPS inside the mouth are heated as well as sides**, exactly as in
D2t. State it; it is not a side-only mechanism.

⚠ **`side_face_factor_above=0.0` on `M1` is MANDATORY, not cosmetic.** Unset, the
sentinel falls back to `side_face_factor` = 0.2, and with `h = 0` against 293 K the
Newton solve returns `q_in < 0` — **`M1` would COOL its above-plane faces** and
stop being a control.

⚠ **DOMAIN WATCH:** assert the deepest column stays **>= 20 mm** above the domain
floor; shorten `stop_time` if not.

⚠ **NO 1 mm RUN IN THIS PACKET.** D2t showed the lever's effect appears as a
*mesh gap*, so **a single mesh cannot see it**: this is a **screening** run and
carries **no convergence claim**.

### 6. Reachability ledger — carry into `PREFLIGHT.md`

| behaviour | call site | reached because |
|---|---|---|
| confined `h` below the plane | `:1773` true branch | every leg sets `h_expr`; `M0` uses the old one, `M1`/`M2` the new |
| confined `h` **above** the plane | `:1773`, gated `jet_above_h` | `M2` sets it; `M0`/`M1` do not |
| above-plane `T_gas` | `:1781`, gated `jet_T_above >= 0` | `M2` only; parser aborts if set without `jet_above_h` |
| side-face term | `:1655-1656` inside `col_in_patch` -> `nfx+nfy>0` | `M1`/`M2` set `side_face_flux=1`; `radius = 0.2` covers the domain |
| per-region factor | `:1653-1654` on `fs[col_c] > 0.0` | `M1` sets 0.0, `M2` sets 1.0 |
| ledger twin | `:2754-2885` | witnessed by `ledger_err` |

## Criteria — write into `CRITERION.md`, hash read-only **BEFORE THE FIRST LEG**

- **(a) THE SCORED ONE — THE FUNNEL.** Mean `|Ø_model − Ø_Meier|` over depths
  **10 / 25 / 50 mm** on `M1`. `M0` measures **35.3 mm** (excess 45 / 36 / 25).
  **Bar: <= 20 mm.** Report the whole ladder and `M2` as well.
- **(b) THE WALL.** `M2 − M1` on Ø at **100 / 150 / 200 mm**, and against the
  85–95 ± 5 band. **Reported, not gated** — one mesh cannot settle it.
- **(c) ITEM 0's WITNESS.** `side_P_face_above` **> 0 on `M2`, exactly 0 on
  `M0`/`M1`**. **Zero on `M2` ⇒ inert ⇒ STOP and report.**
- **(d) NO BLANKET — this one can kill the packet.** ⚠ **Score the MEAN, not the
  max.** Measured on the frozen run beyond r = 76 mm: **mean 0.08 mm but max
  4.00 mm** (a stray two-cell event), and **exactly 0 beyond 90 mm** — so a max-
  based bar would fail on the baseline itself. **Bar: mean recession beyond
  r = 76 mm <= 2 mm on every leg**; report the max and the beyond-90 mm figure
  alongside. Also report the count of above-plane columns with `h > 10 W/m²K`,
  which must be confined to the annulus band. `M0` already strips **33.6 mm** from
  the r = 50–76 mm ring — report that too.
- **(e) VOLUME.** Total excavated x4 against Meier's 3.42 L / 0.5 m pro-rated to
  the depth reached. **Reported.** D2l found the model removes ~1.9x too much.
- **(f) THE LATCH.** Above-plane column count and plane crossings per block, as
  D2t. **Reported.**

### Pre-register, with numbers, before launch

- **P1 — `M0` reproduces the frozen `d2i_gate/LI0_2mm` event log BYTE FOR BYTE up
  to 500 s**, and therefore the measured table 141 / 128 / 113 / 105 / 99 / 92 /
  86. A ±5 mm bar was too weak: `input_feet_d2c` differs from the frozen run in
  two keys, so a miss would have been **configuration, not physics**. **Byte
  identity removes that ambiguity. Any difference ⇒ STOP** — the baseline this
  packet is attributed against is not what it claims to be.
- **P2 — `M1` reduces the collar.** Direction only. **Refuted if Ø @ 10 mm on
  `M1` exceeds `M0`'s.**
- **P3 — `M2` is wider than `M1` at 100–200 mm depth.** **Refuted if `M2 <= M1`
  there.**
- **P4 — (d) HOLDS: no blanket.** **This is the one that kills the packet.**
  Refuted if any leg strips more than 2 mm beyond r = 76 mm.
- **P5 — (a)'s MAGNITUDE IS DECLARED OPEN.** No number predicted. Declaring it
  open is the discipline D2t honoured and D2s failed.
- **P6 — THE FLOOR APPLIES AT SHALLOW DEPTH ONLY, on `M1`/`M2`.** Item 1b's
  ~96–100 mm floor holds where the mouth taper and the early formation set Ø, i.e.
  **10–50 mm depth on `M1` and `M2`**. **Refuted if `M1` or `M2` cuts below 92 mm
  at 10–50 mm depth.**
  ⚠ **It does NOT hold at depth, and must not be scored there.** A wall column at
  r = 46 mm still sees ~430–455 W/m²K **below** the plane and recedes at ~1.16 m/h
  against a 1.5 m/h feed, so it stays below the plane for most of the run and cuts
  ~160 mm deep; a column at 50 mm is overtaken after ~210 s having cut only ~38 mm.
  Expect `M1` ~100 mm at 50 mm depth, 92–96 mm at 100–160 mm, 88–92 mm deeper.
  **`M0` measures 86 mm at 200 mm depth and has no blend at all**, so a blanket
  "below 92 mm at depth" clause would fire on the control. **Deep Ø is REPORTED,
  not gated.**
  **A pass on (a) is "the funnel is removed down to the numerical floor", NOT
  "the model reproduces Meier's diameter".**

### What would make this inert

| # | cause | witness |
|---|---|---|
| N1 | compiled map != intended | item 3's thermal trace — **STOP before launch** |
| N2 | `jet_above_h` not wired at `:1773` | (c): `side_P_face_above` = 0 on `M2` |
| N3 | above-plane state unconfined ⇒ blanket | (d): recession beyond r = 76 mm |
| N4 | `side_face_factor_above` unset on `M1` ⇒ it cools instead of controlling | `side_P_face_above` **negative** on `M1` |
| N5 | stale input: 0.3 m domain or 4 ms | domain watch; campaign log runtime |
| N7 | trace run uses feet descent ⇒ **wall branch never executes** | the second trace leg (prescribed descent, plane below the surface) |
| N8 | `M0` built from `input_feet_d2c` ⇒ two-key config drift | P1 byte identity against the frozen log |
| N6 | `side_face_flux` off | `side_cols_above` = 0 |

### Refuting evidence already on disk

- **The measured Ø(z) table governs P1.** If `M0` disagrees, the measurement wins
  and the discrepancy is the finding.
- **D2l:** with the side-face closure off the model removed **1.9x too much rock
  while ROP read in-band** ⇒ score shape and volume, **never rate**.
- **D2t:** `z(r)` does not converge in the band even where `D(z)` does. A clean
  Ø(z) here does **not** mean the wall is resolved.
- **D2o-0:** a radial corner in `h` caused a step-face cascade; **D2o-0b:** a
  sharp corner gave a +402 K step. **If (a) passes but the wall shows a large
  lateral step, the taper is too narrow** — widen it, do not raise `h`.

## Guardrails

- **NO TUNING TO HIT THE BAND.** `h_ann = 142` declared and never varied;
  `side_face_factor` 0.2 / 1.0 are **frozen stand-ins**; `jet_T_above` declared
  from the run's own exhaust with the jacket drop **stated**, never fitted. The
  taper radii are set by geometry (`foot_r_outer`, the hole wall, the measured
  mouth ring) and confirmed by item 3 — **not adjusted to move Ø(z)**.
- **The above-plane state is an EXPRESSION, never a scalar.** 89 % of columns are
  above the plane and all are in-patch; a scalar is a blanket. This is the fault
  that withdrew D2u.
- **Both tapers >= 8 coarse cells.** D2o-0b failed on a step. **Taper width is a
  NUMERICAL parameter** — narrow it only with the mesh (8 mm on 1 mm), **never to
  move a result.**
- **`r_imp` is the skirt radius `foot_r_outer = 0.040`.** Starting the blend
  inside the pad annulus (28–40 mm) costs the outer pad rock 10–20 % of its `h`,
  slows its recession and slows the feed — **D2o-0's jam mechanism.**
- **The wall branch is `h_ann`, flat.** Do not let the impingement term reach
  above-plane columns: at r = 42–47 mm a blend would deliver **245–405 kW/m²** and
  drive the hole out to the hand-over radius regardless of physics.
- **Score shape and volume, never rate.**
- **Never quote agreement tighter than ±5 mm**, and always state the reference is
  a **lower bound** on Ø.
- **No 1 mm run, no convergence claim, no shape verdict from 2 mm alone.**
- New keys default to sentinels; unset runs **byte-identical**; tree identity
  must stay **2639/2639**.
- `spall.per_cell_removal` **off** (`:4490-4494`, `:4905-4909` unfixed).
  `side_face_h_max` never set. **Do not use `sn_cos_col`** — empty unless
  `spall.surface_normal` is on.
- No AMR, no third mesh, no effective-area weight (D2s withdrawn), no shield, no
  `body_radius`, and **do not smooth the s-switch in this packet** — that is two
  variables; measure it in (f).
- ⚠ **`f = 0.2` is still UNDERIVED** and sets the wall's magnitude. Retire it by
  derivation, never a scan — **not here**.
- Do not touch `ext/`, `bin/`, `obj/`, `build/`, `compile_commands.json`,
  `configure`, `LICENSE`, or unrelated integrators. Frozen studies are read-only.
- **Do not add a guardrail that forbids a cheap verification.**

## Commands

```bash
# Build
EIGEN=$PWD/ext \
  CPLUS_INCLUDE_PATH=/opt/homebrew/include \
  LIBRARY_PATH=/opt/homebrew/lib:/opt/homebrew/Cellar/gcc/15.2.0_1/lib/gcc/current \
  make -j8

# Identity: unset keys must be byte-identical
bash tests/MMWSpalling/studies/tree_sweep_0922/witness/guarded_post_all.sh

# Item 3 -- the compiled-map thermal trace, BEFORE any leg
caffeinate -dis python3 tests/MMWSpalling/studies/d2v_hmap/hmap_trace.py

# Campaign (3 legs, 2 mm, 16 ms, 500 s)
caffeinate -dis python3 tests/MMWSpalling/studies/d2v_hmap/run.py --jobs 3

# Score
/Users/tzetze20/Desktop/code/.venv/bin/python \
  tests/MMWSpalling/studies/d2v_hmap/analyze.py \
  | tee tests/MMWSpalling/studies/d2v_hmap/analysis.txt

# Regressions
python3 tests/MMWSpalling/spall_event/test tests/MMWSpalling/spall_event/output
```

## Follow-on

- **If (a) passes:** 1 mm / 500 s (~3 h, overnight) for the convergence claim,
  then the D3 write-up.
- **If (a) passes but the wall steps:** widen the taper. **Do not raise `h`.**
- **If (a) fails:** (d) and the trace say whether the map reached the rock at all;
  the collar may also be set by the feet/support rule rather than by `h(r)`.
- **Named successors, in order:** per-bin `h(r)` from an annulus closure;
  smoothing the s-switch; deriving `f`.
- **Parked by the user:** 0.5 mm third mesh, AMR, the Fig. 8.8 colour original.
- **Withdrawn:** D2u (scalar blanket), D2s (effective area — any revival needs a
  **per-edge** weight from the drop across each face).
