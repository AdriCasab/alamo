# D2q-2c results — the arbiter: does the corner drain spall the wall?

## Verdict: **NO. P4 (the inert branch) is CONFIRMED.**

**`pcr_lateral_cells = 0` at every step of all four runs**, over 250 s, with up
to **5224** exposed lateral faces present. Per-cell removal never fired once.
Under a cap derived from the physics rather than fitted, **the corner drain does
not carry the wall to its firing point.**

This was pre-registered as P4 and predicted, not rationalised afterwards.
**The cap must not be raised to rescue it** — the next move is a gas closure
that raises the wall temperature honestly.

## 0. Provenance

| item | value |
|---|---|
| `CRITERION.md` sha256 | `6f837a707a5167b8e293a88d90ae161e2aa7c400c41c16025731bdd9b0aa0258`, read-only, hashed **23:41:26** before `output/` existed |
| binary | `20ab657331d72dd25b4b34eddb8698295c90306c63e9400b4b2569adebab18fc` (provenance only, not an oracle) |
| identity, all keys off | **2639/2639 PASS**, via `guarded_post_all.sh`, which verified the binary did not move mid-sweep |
| runs | C1/C0 × 2 mm/1 mm, all reached t ≈ 250 s (2 mm 3.4 min, 1 mm 21.7 min) |

**Frozen, imported by path, never re-run:** `d2j0b_plane/` (metrics *and* the
key-off anchors), `d2l_scan/` (attribution baseline).

## 1. Scorecard

| criterion | result |
|---|---|
| **(a)** mesh gap ≤ 5 % every block, both pairs | **FAIL** — C1 −10.9 / −22.2 / −30.3 % |
| **(e)** ≤ C0 key-off + 15 % | **FAIL** — 5 of 6 blocks over |
| **(f)** key-on does nothing lateral until firing | **vacuously satisfied** — never fired |
| **P4** inert branch | **CONFIRMED** |

**But neither (a) nor (e) failed *because of per-cell removal*, and saying so
would be the central error available here.** The mechanism did not act. Both
numbers are the effect of the capped side flux alone.

### Disk-mean rate [m/h], r < 30 mm, key ON

| block | C1_2mm | C1_1mm | gap | C0_2mm | C0_1mm | gap |
|---|---|---|---|---|---|---|
| 50–100 | 1.4400 | 1.2827 | **−10.9 %** | 1.4722 | 1.4180 | −3.7 % |
| 100–150 | 1.4260 | 1.1091 | **−22.2 %** | 1.5056 | 1.4375 | −4.5 % |
| 150–200 | 1.3592 | 0.9472 | **−30.3 %** | 1.4973 | 1.4668 | −2.0 % |

## 2. Attribution — the pre-registered baseline does its job

D2l's `f01` is side-face flux at `f = 0.1`, i.e. h_side ≈ 70 W/m²K against my
cap of 80, **with no per-cell removal in the code at all**:

| block | C1_2mm (cap 80, per-cell ON) | `f01_2mm` (h≈70, per-cell ABSENT) |
|---|---|---|
| 50–100 | 1.4400 | 1.4638 |
| 100–150 | 1.4260 | 1.4450 |
| 150–200 | 1.3592 | 1.3748 |
| gaps | −10.9 / −22.2 / −30.3 % | −10.6 / −20.2 / −27.2 % |

**Within 2–3 % on rate and 1.5–3 points on gap.** The campaign reproduces a pure
side-heating run. Per-cell removal is not merely unhelpful here — it is
**undetectable**, which is exactly what a zero firing count predicts.

## 3. Why (e) "fails", stated honestly

(e) was pre-registered against **C0 key-off**, following the packet's own worked
example (D2k's 1.651 vs C0's 1.085 = +52 %). By that anchor C1 exceeds the
ceiling in 5 of 6 blocks. **That verdict stands as pre-registered.**

But the same data read against **C0 key-ON** (1.47 / 1.51 / 1.50 at 2 mm) would
**pass comfortably** — because capped side heating raises *every* configuration,
including the no-plane control where nothing anomalous happens. Both readings
are recorded; the pre-registered one governs; and the attribution in §2 shows
the excess belongs to side heating, not to per-cell removal.

## 4. A genuine finding that is not about the mechanism

**Capped side-face heating largely fixes the cascade when there is no exclusion
plane.** C0's mesh gap goes from **−10.4 / −15.3 / −20.1 %** (key off) to
**−3.7 / −4.5 / −2.0 %** at h ≤ 80. With the plane, C1 stays bad
(−10.9 / −22.2 / −30.3 %).

Consistent with D2l's `f_min = 0.5`: cascade suppression scales with side-heating
strength, and 80 W/m²K is far below the 350 W/m²K that `f = 0.5` implies. The
exclusion plane is what keeps C1 cascading, by zeroing `h` above it — visible in
`side_P_face`, only **35 W** on C1_2mm against **132 W** on C0_2mm.

## 5. ⚠ LATENT DEFECT, verified in code, not triggered by this campaign

**The cluster labeller scans the whole column for `Sp >= 1`**
(`MMWSpalling.H:4490–4494`):

```cpp
for (int k = klo; k <= khi; ++k)
    if (Sp_a(i, j, k) >= 1.0) fire[col_idx(i, j)] = 1;
```

That was safe **only because `Sp_field_mf` was previously written solely at
`k == ksurf`**. D2q-2b item 6 writes `Sp` at non-top exposed cells, so the moment
a **lateral** cell reaches `Sp >= 1` its column is marked firing, `clu(i,j,k_top)
> 0.0` opens the detach gate (`Removal.H:652`), and **the column TOP spalls
vertically because a SIDE face was ready.**

**It did not bite here** — lateral `Sp` never exceeded 0.4533, so these numbers
are unaffected. **But it would fire precisely when the mechanism starts
working**, converting the intended lateral removal into a spurious vertical one.
Any future packet that raises the wall temperature must fix this first, or its
first positive result will be an artefact.

## 6. Pre-flight record

- **Stop check:** `T_fire` = **822.2 ± 9.5 K** (measured at Sp ≈ 1, not
  extrapolated); cap needs **149.1 s** of continuous exposure to reach it;
  run is 250 s ⇒ proceed. Independent reproduction of the packet's `h*`: at
  h = 88.5 the wall reaches 822 K in **121.9 s** against the stated 121 s dwell.
  **The check was right to pass the packet and right about the risk** — the
  caveat recorded with it (a lateral face is exposed for well under 250 s) is
  exactly what happened.
- **Lateral and top faces lie on ONE `Sp(T)` curve** (−0.89 / −1.05 / +3.03 % in
  the well-populated bins) — the "same physics" claim, measured.
- **Contiguity grep:** no third bug. Two live `rem(k+1)` reads, both audited as
  correct. One stale comment corrected (`BeamEnergyBalance` still documented
  invariant (a) in its retired contiguity form).
- **Harness guard** `guarded_post_all.sh` built, self-tested, and used.

## 7. Decision, per the pre-registered rule

**Inert (P4) ⇒ report it.** The next move is the gas closure that raises the wall
temperature honestly — **not a bigger cap, not a larger `side_face_factor`.**
Neither was touched.

**What this campaign does and does not establish.** It establishes that under a
cap set below `h*`, an exposed wall in this arbiter never reaches its firing
point in 250 s, so per-cell removal cannot act. It does **not** establish that
the mechanism is wrong, untestable, or unable to spall a wall that *is* hot
enough — that remains untested, and §5 must be fixed before it is tested.
