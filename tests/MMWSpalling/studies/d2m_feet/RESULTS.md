# D2m results — is the cut diameter set by the feet or by the flame? (2026-09-22)

**Verdict: MIXED, and reported as mixed — both mechanisms are active, on different quantities.**

- **The cut diameter is a support-rule artefact.** min Ø tracks `2 × foot_r_outer` to within
  **0.5–2.1 mm** across a 32 mm swing of stance, in **both** directions: 63.5 / 79.3 / 87.0 / 93.9
  mm against stances of 64 / 80 / 88 / 96. The model's hole is the width of its own tool stance.
- **The flame's erosion is stance-invariant where it matters.** At r = 45–50 mm the **absolute**
  recession rate is the same to **5–7 %** across the whole ladder, while the centre moves **39 %**.
  The far field is set by the flame; only the centre responds to the stance.
- **Meier's hole was wider than his tool** (Ø 85–93 on Ø 80 feet). Ours cannot be, by construction.
  **The mechanism that removes rock beyond the stance is missing**, and no wall-`h` work supplies
  it.

**The measurement 5a was added to get, answered:** the model's floor is **not flat**. At the
control, v(r)/v_c is 0.83 at 29 mm and 0.62 at 39 mm — **already falling well inside 40 mm** —
and it never reaches 0.1·v_c within 70 mm, so `w_edge` is undefined for four of the five cases.
**The model drills a cone with a long tail, not a cylinder with an edge.** Meier's hole is
cylindrical, so this is a shape failure the support rule alone does not explain.

## 0. Hashes, binary, provenance

- **Binary** `bin/mmwspalling-3d-g++` = `62dea451e6b1d1908e0ea12f730cdc3f7c01790f19d6d63b8c2b88befe0f67fa`,
  **verified unchanged at the end**. No source edit, no build, no `make`, nothing committed.
- **`PREDICTIONS.md`** sha256 `d735f95fef9e12e9502cb37fe8c18e9f784bca8a1f08ec90843ad16af6654c9a`,
  read-only, **12:22:10**. All five runs completed **12:40–12:42**, so every run post-dates the hash.
- **Frozen, imported by path, none edited:** `d2l_scan/` (Ø(z), volume, matched-depth estimators),
  `d2i_gate/` (the LI0 keys), `d2j0b_plane/` (the runner), `d2b_feet_rop/analyze.py`. D2j-0's
  recession estimators and `d2c_steady`'s `az_rwall` were **copied** into this study rather than
  run in their frozen folders.

### The mid-flight amendment — read this before the numbers

`ACTIVE_STEP.md` was edited at **12:16:08**, after `PREDICTIONS.md` was first hashed (12:15:53,
`6d1daf69…`) and while the first launch was starting. It **changed the deciding metric**: item 4
promoted the transition radius to deciding, and a new item 5a specified the recession profile
(2 mm bins to 70 mm, D2j-0 estimators, half-value against the **central** value). My hashed file
had that prediction as a *non-deciding* extra at 4 mm bins against half of *peak*.

Per the standing rule, the five in-flight runs were **killed at t ≈ 22–28 s of 250 and deleted** —
**no `.done` was written, no analysis was run, and no output was inspected** — `PREDICTIONS.md` was
amended and re-hashed, and everything was relaunched from t = 0. The amendment changed **no input
key**, so the runs were physically unaffected and were discarded anyway. Full record at the top of
`PREDICTIONS.md`; the killed campaign log is kept as `campaign_killed_amendment.log`.

### Two source couplings the packet did not anticipate, declared before launch

1. **`foot_r_inner` silently sets `nozzle_collision_radius`** (`MMWSpalling.H:5695`,
   `feet ? foot_r_inner : radius`) and the D2i keys never set it, so `S3648` would have moved the
   collision guard 28 → 36 mm. **Pinned to 0.028 in every case** — the control's own effective
   value. **The identity guard proves the pin is inert.**
2. **`foot_r_outer` also defines the free-surface boundary** (`MMWSpalling.H:2062`). It cannot be
   pinned away, and **it turned out to be a much larger confound than the pre-launch bound
   allowed** — see §6.

## 1. Identity

| run | vs | rows (t ≤ 250 s) | verdict |
|---|---|---|---|
| `R40` | `d2i_gate/LI0_2mm` | 2605 | **byte-identical** |

**PASS.** ACTIVE_STEP names the partner "D2l's `M_f0`"; no such run exists — D2l's f = 0 row was
read from `d2i_gate/output/LI0_2mm`, and that is what `R40` reproduces. Named explicitly rather
than substituted silently.

## 2. The ladder

| case | ring [mm] | stance Ø | **min Ø** | 2×r_outer | Δ | ROP [m/h] | Δ ROP | absorbed [W] | V̇ [cm³/s] | depth-mean Ø |
|---|---|---|---|---|---|---|---|---|---|---|
| R32 | 28–32 | 64 | **63.5** | 64 | −0.5 | 1.925 | +29.1 % | 1258 | 3.59 | 91.8 |
| R40 | 28–40 | 80 | **79.3** | 80 | −0.7 | 1.491 | — | 1639 | 4.58 | 110.5 |
| R44 | 28–44 | 88 | **87.0** | 88 | −1.0 | 1.301 | −12.8 % | 1841 | 5.11 | 120.0 |
| R48 | 28–48 | 96 | **93.9** | 96 | −2.1 | 1.179 | −21.0 % | 2021 | 5.62 | 128.1 |
| S3648 | 36–48 | 96 | **94.8** | 96 | −1.2 | 1.154 | −22.6 % | 2052 | 5.70 | 129.4 |

- **H-feet holds on min Ø, in both directions.** My pre-registered **H-asym is refuted on its
  inward branch**: I predicted R32 would stay at 74–80 mm because rock at 32–40 mm is heated
  identically regardless of the foot radii. It collapsed to 63.5. Being heated is not enough —
  rock outside the stance is never brought under the nozzle, its stand-off grows, and it lags.
- **But this was pre-registered as weak evidence** (`PREDICTIONS.md` §2a): `z_foot` is the 0.9
  quantile of annulus heights, so ~90 % of rock out to `foot_r_outer` lies at or below the feet and
  the hole is ≥ 2·r_outer at the feet depth *almost by construction*. **Confirming it carries
  little information, and it is not read as more than it is.**
- **ROP direction PASSES**: monotone fall, 1.925 → 1.179 m/h. Predicted sizes R32 +10…+25 %
  (observed **+29.1 %**, slightly over), R44 −10…−25 % (**−12.8 %**, in range), R48 −25…−45 %
  (**−21.0 %**, slightly under). Direction right, magnitudes roughly right at both ends.
- **Absorbed power *rises* with stance (1258 → 2021 W) while ROP *falls*** — a wider hole absorbs
  more and descends more slowly.

## 3. `S3648` vs `R48` — the outer edge alone sets it

min Ø **94.8** vs **93.9** (Δ +0.9, predicted within 3); ROP **1.154** vs **1.179** (Δ −2.1 %,
predicted within 5 %); r_half 54.3 vs 53.5. **Prediction holds.** Dropping the inner 28–36 mm band
changes almost nothing, as expected from a 0.9 nearest-rank quantile taken near the top of each pad.

## 4. The erosion edge — and a retracted reading

**The pre-registered test FAILS as written.** H-feet on `r_half` required
|r_half − r_outer| ≤ 2 mm; the offsets are **+5.4 / +6.6 / +6.6 / +5.5 mm**.

**A slope reading I added and then retracted, in full.** Seeing that r_half spans 16.0 mm for
16 mm of r_outer (`d r_half/d r_outer` = **1.02**), I first wrote that the erosion edge moves
one-for-one with the stance and that H-feet therefore held in substance. **That was wrong.**
`r_half` is a crossing of `0.5·v_c`, and `v_c` is **not constant** across the ladder — it falls
3.317 → 1.968 m/h. The level moves, so the crossing moves, even if the profile does not.

**Absolute rates settle it:**

| case | r_outer | v_c | 25 mm | 35 mm | 45 mm | 50 mm | 60 mm |
|---|---|---|---|---|---|---|---|
| R32 | 32 | 3.317 | 2.564 | 1.801 | 1.303 | 1.155 | 0.453 |
| R40 | 40 | 2.443 | 2.249 | 1.692 | 1.278 | 1.137 | 0.882 |
| R44 | 44 | 2.153 | 1.984 | 1.649 | 1.254 | 1.127 | 0.880 |
| R48 | 48 | 1.968 | 1.777 | 1.589 | 1.221 | 1.106 | 0.862 |
| S3648 | 48 | 1.939 | 1.757 | 1.575 | 1.219 | 1.096 | 0.860 |

Ladder spread as % of mean: **25 mm: 39 %**, 35 mm: 14 %, **45 mm: 7 %**, **50 mm: 5 %**,
60 mm: 54 % (the last driven by R32 alone, whose narrow hole puts 60 mm out of reach; excluding
R32 the four remaining agree to **2.5 %**).

**At 45–50 mm the absolute erosion rate is invariant to 5–7 % across a 16 mm change of stance,
while the centre moves 39 %.** The flame's far-field erosion is not set by the stance — only the
centre is. **H-thermal holds on the erosion edge; H-feet holds on min Ø.** That is the packet's
"mixed" row, reached for a different reason than H-asym predicted.

**Window cross-check, as pre-registered.** Matched-geometry and fixed-[150, 250] s `r_half` differ
by up to 5.2 mm (R32), more than one cell, so **both are reported and no single reading is
forced**. The fixed window compresses the spread (slope ≈ 0.72) because the faster cases have
already drilled past the comparison depth — which is exactly why matched geometry is the deciding
window.

**`w_edge` is undefined for R40, R44, R48 and S3648**: the profile never falls to 0.1·v_c within
70 mm. That is not a measurement failure — **the erosion has no edge.**

## 5. The floor is a cone, not a cylinder — the answer 5a was added to get

v(r)/v_c by 2 mm annulus at the control (R40): 1.00 at the axis, 0.95 at 21 mm, **0.83 at 29 mm**,
**0.62 at 39 mm**, 0.52 at 45 mm, and still ≈ 0.36 at 69 mm.

5a asked whether the floor recedes **flat out to ~45 mm**, as Meier's cylindrical hole requires,
or is **already falling at 40 mm**. **It is already falling far inside 40 mm, and it never stops
falling.** Meier's implied floor half-width is 42.5–46.5 mm; at that radius our model is removing
rock at roughly **half** the centre rate, so the floor is a cone with a long skirt. This is a shape
failure that the support rule does **not** explain and that the next packet must address on the
heat-distribution side.

## 6. The free-surface confound — larger than the pre-launch bound allowed

`jet_fs_cols` over the window: **0.0 / 17.9 / 87.2 / 202.7 / 247.0** for R32 / R40 / R44 / R48 /
S3648. The pre-launch bound (`PREDICTIONS.md` §1) took D2l's `P_fsoff0`, which removed the free
surface entirely **when it had only ~18 columns**, and moved min Ø by +0.3 mm and ROP by −2.3 %.
At R48 there are **11× more** free-surface columns, so **that bound does not cover the outer end of
the ladder**, and I am flagging it rather than leaning on it:

- **min Ø survives easily.** The confound would have to be ~100× its measured size to produce a
  30 mm tracking. The H-feet conclusion on min Ø is robust.
- **The ROP percentages carry this uncertainty** and should be read as direction-plus-rough-size,
  not as calibrated numbers.
- **A clean ladder would need `jet_free_surface = 0` throughout** — out of scope here ("not the
  gas closure"), and named for the next packet.

Ledger closure ≤ 4.8e-15 and cap ratio ≤ 0.21 on every run. J/mm³ to rock is **1.40–1.44**
(1.22–1.26 × the floor 1.147) and **flat across the ladder** — the same fixed removal cost D2l
found, unchanged by stance.

## 7. What the next packet inherits

1. **No Meier diameter may be quoted until something removes rock beyond the stance.** Our hole is
   the width of the tool; Meier's was 5–13 mm wider than his. This is now measured, not inferred.
2. **The real shape failure is the cone.** The floor is already at 0.62·v_c by 39 mm and has no
   edge. The next packet belongs on the **radial heat distribution**, which D2l showed no existing
   key can reach — `jet_free_surface` and the far-field exponent move axial decay and mouth
   dilution, never the radial profile.
3. **Orientation only, explicitly not a score** (single mesh, 250 s, relative study; the packet
   forbids Meier scoring): R44 happens to land at ROP 1.301 m/h and min Ø 87.0 mm, both inside
   Meier's ranges, while still removing **5.11 cm³/s against 2.47** — because the depth-mean Ø is
   120 mm. **Matching rate and cut diameter together still leaves the volume 2× high.** The funnel
   above the cut, not the cut, is where the volume goes.
4. **A stance ladder with `jet_free_surface = 0`** would make the ROP numbers clean (§6).
5. **Untouched and still owed:** the 2639-file sweep on the working tree before any rebuild.
