# D2l Gate 0 results — why the model drills 2.7× too fast (2026-09-22)

**Verdict: the window is EMPTY. STOP — no code was written, Gate 1 was not run, and the
`side_face_h_scale` closure was not started.**

`f_min` = **0.5** (below it the D2j-0b cascade returns) and `f_rate` = **0** (the burner rate
re-enters Meier's 1.3–1.6 m/h only with the side term switched off entirely). `f_rate < f_min`
with no overlap at all, so **no value of `f` — and therefore no single scale on the side-face
heat transfer — can suppress the mesh cascade and produce Meier's rate at the same time.** That
is the pre-registered STOP condition.

**The larger finding is that the closure was aimed at the wrong quantity.** Across the whole
ladder the model removes rock at an almost constant **1.32–1.43 J/mm³ of absorbed energy**, about
1.2× the thermodynamic floor. Volume rate is therefore just proportional to absorbed power, and
**at f = 0 — the side term completely off — the model already absorbs 6.56 kW (17.3 % of Meier's
38 kW HHV) and removes 4.58 cm³/s against Meier's 2.47.** Matching Meier's volume at the model's
conversion needs an absorbed **3.46 kW, i.e. 9.1 % delivery**. The factor ≈ 2 discrepancy is
present with the side-face term *off*; D2k's term made an existing error larger, it did not create
it. **A wall-`h` closure cannot reach it.**

## 0. Hashes, binaries, provenance

- **Binary:** `bin/mmwspalling-3d-g++` = `62dea451e6b1d1908e0ea12f730cdc3f7c01790f19d6d63b8c2b88befe0f67fa`
  (the D2k campaign build). **No source change, no build, no commit.** Gate 0 was key-only
  throughout, as the packet requires.
- **`PREDICTIONS.md`** sha256 `b1b6686f992754f69bdea78c0ace87335ccba3e29ad2e63cf3fe73ecd9b2d613`,
  written 2026-09-22 10:26:57, made read-only **before `output/` existed**.
- **Gate 1 was NOT run** (it is conditional on a Gate 0 "proceed"). `bin/mmwspalling-3d-g++-amr` =
  `2eeda1a5f95cdf547de9d4a26da03df61835a4c96d56fb556045e6402dabb693` is recorded here only so the
  next packet can identify the build; the working-tree AMR delta remains unverified and **no
  rebuild happened**.
- **Frozen studies imported by path, none edited:** `d2j0b_plane/` (arbiter keys and the
  `metrics`/`band_timing`/`r_front` estimators), `d2i_gate/` (the LI0 Meier keys, which pull
  d2h → d2f → d2e → d2c), `d2k_sideface/` (identity reference), `d2b_feet_rop/analyze.py` (`Run`).
  `az_rwall` was **copied** into `analyze.py` from `d2c_steady/score.py` rather than running that
  script, which writes into its own frozen `output/`.

### Harness identity — the guard on reaching through two frozen studies

| run | vs D2k | rows | verdict |
|---|---|---|---|
| `A_f10` | `C1_2mm` | 625 | **byte-identical** |
| `M_f10` | `LI0_2mm` (t ≤ 250 s) | 2605 | **byte-identical** |

`run.py --identity` **PASS**. The f = 1.0 legs were deliberately re-run rather than read from
D2k's output, so that this check exists.

## 1. Deviation from the packet, declared before launch

ACTIVE_STEP puts 0c at f = 1. **Both probes are inert there**, for reasons that are properties of
the geometry the side flux produces, and this was measured from the existing D2k run and written
into `PREDICTIONS.md` §1 *before* anything was launched:

- the far law is `pow(12D/max(12D, s), n)` with 12D = **90 mm**, and at f = 1 the hole is
  self-similar at `jet_s_c` 67–69 mm, so every column sits in the clamp region where the factor is
  exactly 1 for any `n`;
- `jet_fs_cols` = 0.0 through the whole scoring window, so free-surface dilution is already off.

Each probe was therefore run **at f = 1 as specified** (to settle inertness by measurement rather
than by argument) **and at f = 0** (`P_fsoff0`, `P_flat0`), where both are live and which is the
configuration the competing explanation actually describes. No new key was used; nothing in the
packet was dropped.

A second deviation, forced by the data and recorded here: ACTIVE_STEP says to run the 1 mm legs
"at whichever f first fails at 2 mm and the one above it". **Nothing failed at 2 mm at any f** —
the single-mesh band-stall signature never fired, not even at f = 0.1 (§2). Rather than pick a
boundary by eye after seeing results, **the full 1 mm ladder was run** — four legs instead of two,
no discretion about which. This was the right call: the mesh gap degrades steeply at f ≤ 0.2 while
the 2 mm stall detector stays silent, so the packet's cheap trigger would have missed `f_min`
entirely.

## 2. Gate 0a — `f_min` = 0.5

**The 2 mm single-mesh signature does not work.** No band stalls inside 28 mm and no `r_front`
walks inward at *any* f, including f = 0.1, although the key-off C1 at 2 mm does stall. The
declining 2 mm block rates at f = 0.1 (1.464 / 1.445 / 1.375) are a hint, but the criterion does
not fire. **The cascade is visible only in the mesh gap.**

| f | 50–100 s | 100–150 s | 150–200 s | verdict |
|---|---|---|---|---|
| 1.0 | **+3.4 %** | **+1.5 %** | **−3.2 %** | pass |
| 0.5 | **+3.7 %** | **+1.4 %** | **−2.1 %** | pass |
| 0.2 | −6.4 % | −8.0 % | −14.0 % | **FAIL** |
| 0.1 | −10.6 % | −20.2 % | −27.2 % | **FAIL** |
| 0 (key off, D2j-0b C1) | −13.5 % | −28.7 % | −45.5 % | FAIL (reference) |

**`f_min` = 0.5.** The f = 0.1 gaps are already most of the way to the key-off cascade, and both
failing rows grow with time in the way D2j-0b identified. **P1 (predicted `f_min` = 0.5) holds.**

**Caution, not a finding:** both passing rows drift negative in the last block (−2.1 %, −3.2 %)
over only 200 s. That is the week's recurring failure mode — "converged at two meshes over a short
window" — and **P2 (that `f_min` rises with refinement, because the step-face sink scales as 1/dx
while the side supply does not) is UNTESTED here.** The queued 0.5 mm arbiter leg is its test.

## 3. Gate 0b — `f_rate` = 0, and the cost is flat

Compared at **matched geometry**: each run at the time its feet reach the common depth **131 mm**,
the deepest depth every run on the ladder reached. Never at matched time.

| f | ROP 150–250 s [m/h] | absorbed (qtr) [W] | side/face | V̇ [cm³/s] | depth-mean Ø | min Ø | TSE(HHV) | TSE(to-rock) |
|---|---|---|---|---|---|---|---|---|
| 1 | **4.066** | 3204 | 53 % | 9.73 | 101.9 | 78.0 | 3.9 | 1.32 |
| 0.5 | **3.759** | 2524 | 47 % | 7.58 | 97.8 | 77.1 | 5.0 | 1.33 |
| 0.2 | **3.087** | 1980 | 37 % | 5.75 | 100.4 | 77.7 | 6.6 | 1.38 |
| 0.1 | **2.137** | 1837 | 24 % | 5.25 | 104.1 | 77.5 | 7.2 | 1.40 |
| 0 | **1.491** | 1639 | 0 % | 4.58 | 107.6 | 78.0 | 8.3 | 1.43 |

Bands: ROP 1.3–1.6 m/h · Ø 85–93 mm · TSE(HHV) 15.4 · TSE(to-rock) 4–6 J/mm³. Ledger ≤ 8.1e-15
and cap ratio ≤ 0.39 on every run, so nothing here is a budget artefact.

- **`f_rate` = 0.** Even f = 0.1 gives 2.14 m/h, well above the band. **P3 (predicted `f_rate`
  ≈ 0.1) is REFUTED** — the rate response is steeper at small f than predicted.
- **Ø does not respond to f.** Minimum Ø is **77.1–78.0 mm at every f including f = 0** — flat to
  0.9 mm across the whole ladder and below Meier's 85–93 throughout. **P4 holds.** Cutting `f`
  buys rate and buys no shape: the D2b failure mode restated.
- **TSE(to-rock) is 1.32–1.43 J/mm³ at every f**, against the 4–6 band and a thermodynamic floor
  ρc_pΔT_fire = 1.147. **P5's core claim holds: the to-rock band fails at every f, including the
  one that fixes the rate.** TSE(HHV) never reaches 15.4 anywhere on the ladder — it is 8.3 even at
  f = 0. *(P5's predicted numbers were partly wrong: TSE(HHV) at f = 0.1 was predicted 14–15 and
  came out 7.2, because V̇ at low f was underestimated.)*

### Two traps in the shape metric, both found by checking rather than by results

1. **Ø at matched time is meaningless here.** Ø(t_end) at the matched depth runs
   104 / 95 / 88 / 78 / 94 mm across f = 1 → 0 and tracks **dwell** (119 / 98 / 45 / 0 / 258 s),
   not f: the slowest run (f = 0) shows the second-widest hole. D2k's caveat is confirmed
   quantitatively.
2. **A single depth is fragile.** The Ø(z) profiles of different runs **cross**, and the matched
   feet depth landed near a crossing. Reported at one depth, every probe in §4 reads "Ø unchanged";
   the full profiles differ completely. The profile is therefore reported in full:

| f | 10 mm | 30 mm | 50 mm | 70 mm | 90 mm | 110 mm | 130 mm |
|---|---|---|---|---|---|---|---|
| 1 | 121 | 107 | 104 | 104 | 103 | 94 | 78 |
| 0.5 | 126 | 107 | 99 | 96 | 94 | 88 | 77 |
| 0.2 | 136 | 117 | 104 | 97 | 92 | 85 | 78 |
| 0.1 | 141 | 122 | 109 | 102 | 96 | 86 | 78 |
| 0 | 145 | 126 | 113 | 106 | 100 | 88 | 78 |

Even at matched geometry the mouth is **wider at lower f**, because a slower run spends far longer
with the upper hole exposed. **Freshly cut rock is ~78 mm Ø regardless of f**; everything above is
time-exposure. The funnel is not an `f` effect.

## 4. Gate 0c — no probe produces the discriminating signature

Deciding metric: **Ø at the matched depth z_m**, as fixed in `PREDICTIONS.md` §3 before any run.

| probe | f | change | ROP [m/h] | Δ rate | Δ Ø (deciding) | Δ depth-mean Ø | signature? |
|---|---|---|---|---|---|---|---|
| `P_flat` | 1 | n 1.0 → 0.5 | 4.066 | +0.0 % | +0.0 | +0.0 | no |
| `P_fsoff` | 1 | free surface off | 4.066 | −0.0 % | +0.0 | +3.0 | no |
| `P_flat0` | 0 | n 1.0 → 0.5 | 1.460 | −2.1 % | +0.0 | +0.1 | no |
| `P_fsoff0` | 0 | free surface off | 1.457 | −2.3 % | +0.3 | +17.1 | no |

- **`P_flat` at f = 1 is inert to round-off, as predicted.** Every physical thermo column
  (`nozzle_z`, `foot_z`, `ledger_P_robin`, `jet_P_face`, `patch_P_robin`, `jet_s_c`, `jet_T_stag`,
  `side_P_face`) differs by **exactly zero**; only `ledger_err` (~1e-16) and `jet_P_decay`
  (~1e-12 W on a ~1e-12 base) move. **P6 holds in substance**, though strictly it predicted
  byte-identity and one round-off-level diagnostic differs.
- **The same probe at f = 0 is live** (`nozzle_z` 2 mm, `P_robin` 615 W, `s_c` 18 mm), confirming
  that the inertness at f = 1 is caused by the geometry, exactly as diagnosed.
- **P8 is partly refuted:** `P_fsoff0` was predicted to raise the rate; it lowers it 2.3 %, because
  removing the free surface removes entrained mass, so the wall jet cools faster along r.
- **P7 and P9 hold.** No probe gives rate down *with* Ø up.

**A reviewer should check this one.** Scored on the **depth-mean** instead, `P_fsoff0` reads
**+17.1 mm with rate −2.3 % — the signature**, and the verdict flips. It is not treated as one for
two reasons: it is not the pre-registered metric, and it moves the **wrong way** — the depth-mean
is already too *wide* (112 mm against 85–93), so +17 mm is a worse funnel, not a better hole. The
quantity that fails *low* is the minimum Ø, and that moves +0.3 mm. I switched the deciding column
to the depth-mean while exploring, saw the verdict flip, and reverted to the pre-registered metric.

**0c is a weak test of the competing explanation, and does not refute it.** The two available keys
move the far-field *axial* decay and the mouth dilution. Neither reshapes the **radial** profile,
which is what "h(r, s) is too peaked" actually means. **Gate 0c shows that the existing keys cannot
test explanation 2, not that explanation 2 is wrong.**

## 5. The decision, and what the next packet inherits

**`f_rate` = 0 < `f_min` = 0.5 — the window is EMPTY.** Per ACTIVE_STEP: **STOP. Write no code.**
Gate 1 was not run and no closure was written. **`f` was not tuned, and no value from this scan may
be kept as a setting** — `f` corrects area, and the scan was a diagnostic.

**The wall-`h` closure is dead on its own terms, and for a sharper reason than the empty window.**
The volume-per-absorbed-joule is flat at 0.70–0.76 cm³/kJ across the ladder, so:

| f | absorbed, full hole | % of 38 kW HHV | V̇ [cm³/s] vs Meier's 2.47 |
|---|---|---|---|
| 1 | 12.82 kW | 33.7 % | 9.73 (3.9×) |
| 0 | 6.56 kW | 17.3 % | 4.58 (1.9×) |
| *needed* | *3.46 kW* | *9.1 %* | *2.47* |

- **With the side term entirely off the model is still ≈ 1.9× too fast by volume.** The error D2l
  set out to explain is not created by D2k's term; D2k doubled an error that was already there.
- **Matching ROP at f = 0 matches one number while volume runs 1.9× high**, because the hole is the
  wrong shape — 78 mm where it cuts, 145 mm at the mouth. That is the "right rate for the wrong
  reason" the packet warned about, and it applies to **D2i/D2d themselves**, not only to a tuned f.
- Either delivery to the rock is ≈ 9 % (below every estimate on record: the packet's 30 %, the
  README's implied ~17 %), or the model removes rock ≈ 2× too cheaply per absorbed joule, or both.
  **Deciding between those is the next packet**, and neither is a wall-`h` question.

**Named for the next planner, not started here:**

1. **The delivery/cost split.** An independent anchor on absorbed power is needed — Meier's
   unreported 160 L/h cooling-water ΔT is the one measurement that would settle it (already on the
   user's list in ROADMAP).
2. **Radial distribution, untested.** The funnel (145 mm mouth against 78 mm at the cut) is the
   shape failure, and no existing key reaches it. A radial-profile probe is key-only and cheap.
3. **`f_min` under refinement (P2), untested.** If it rises, D2k's arbiter PASS is itself only a
   two-mesh statement. The queued 0.5 mm leg is the test, and it should run before any further
   claim about the side-face remedy.
4. **The working-tree AMR delta is still unverified** — Gate 1 was not reached, so the 2639-file
   sweep on the current tree remains owed before *any* rebuild.
