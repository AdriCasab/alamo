# D2g Stage 1: long mesh check of the D2f support rule (q = 0.9) at 0.40 m

## 0. Criterion (fixed before launch)

The criterion is `CRITERION.md` (read-only), sha256 `ba914d27bd88f1b7c3c046defecd2517acf52e29a5348f6c4f43e8af5f99a43a`, mtime 2026-09-19 18:33:21,
written before either run was launched. Its content:

**Pair.** `L09_2mm` and `L09_1mm`: D2f Q09 keys (pads, q 0.9, no clip),
0.40 m full domain, stop cap 700 s, ended by the bottom watchdog.

**Burner ROP.** The least-squares slope of `nozzle_z` against time over the
thermo rows with a ≤ t ≤ b (and t > 0), in m/h. This is `d2e run.blocks()`
and `analyze_pairs.window_rop`.

**t_m** = min(t_end of the two runs), where t_end is each run's last thermo
row.

**Gap** = ROP(1 mm) / ROP(2 mm) − 1.

**PASS requires both:**
- **(a)** the matched-window gap over [150 s, t_m] is within **±5 %**;
- **(b)** every **100 s** block is within **±5 %**. The blocks start at
  150 s: [150, 250], [250, 350], …, and the last one is [150 + 100k, t_m].
  That last block is dropped if it is shorter than 50 s.

**Otherwise FAIL.** No other window, block length or offset is evaluated for
the verdict.

**Reported, not deciding:**
- the 50 s blocks from 0 s ([0, 50], [50, 100], …, the last ending at t_m if
  ≥ 25 s), for continuity with D2f;
- the growth diagnostic: the least-squares slope of the 100 s block gaps
  against their block mid-times, in percentage points per 100 s.

**Outcomes:**
- **PASS:** D2g Stage 2 (scored runs) is planned, quoting "mesh-checked to
  t_m at 0.40 m".
- **FAIL:** stop. The next packet diagnoses the late drift. The leads are
  named from the logs (rim spall rate, flank removal, T_rec, s_c,
  `jet_r_reach`, pinned share) and not tested here.

## 1. Identity and binary

- **Binary** `6d047b50bdcde619952aff907d88f89ed0cae8f9d3ed645b43ba0ac2cd647005`,
  verified at the start (18:33) and at the end. No source edit, no build.
- **Identity PASS.** Rows to t ≤ 299.9 s of `thermo.dat` are **byte-identical**
  to the D2f Q09 runs, at both meshes (3124 rows each). The only key that
  differs from Q09 is the stop cap (700 s against 300 s), plus the 1 mm plot
  interval.

## 2. Runs

| run | status | t_end | wall | note |
|---|---|---|---|---|
| L09_2mm | bottom watchdog | 573.2 s | 51 min | |
| L09_1mm | bottom watchdog | 621.4 s | 1202 min reported | includes a user pause (SIGSTOP 19:44 → SIGCONT 20:06+) and long idle/loaded periods overnight; the clean rate is ≈ 19 wall-s per simulated s, so the computing time is ≈ 3 h |

- Both runs ended on the bottom watchdog (the centre reaching 40 mm above the
  domain floor), not on the 700 s cap, and neither stalled.
- **The 1 mm run outlived the 2 mm one** (621 against 573 s) because its hole
  is shallower, which is the failure itself.

## 3. The pair (t_m = 573.2 s)

**Matched window 150–573 s: 2 mm 1.493 m/h, 1 mm 1.067 m/h, gap −28.5 %.**

| 100 s block (criterion) | 2 mm ROP | 1 mm ROP | gap |
|---|---|---|---|
| 150–250 s | 1.492 | 1.466 | −1.8 % |
| 250–350 s | 1.511 | 1.322 | **−12.5 %** |
| 350–450 s | 1.510 | 0.903 | **−40.2 %** |
| 450–550 s | 1.424 | 0.502 | **−64.8 %** |

- **Growth diagnostic: −21.7 points per 100 s.**
- 50 s blocks (continuity with D2f): +2.0, +0.1, −1.9, +0.1, −3.0, −4.3,
  −17.5, −27.7, −51.9, −54.8, −77.8 %. The first six reproduce D2f Q09
  exactly.
- The 2 mm run holds 1.38–1.54 m/h across the whole run. **All of the
  divergence is the 1 mm run slowing down.**

## 4. Per-block diagnostics

| block | rim spall 38–40 mm [mm/min] | annulus mean | 40–60 mm removal [cm³/s] | mean T_rec [K] | mean s_c [mm] | mean reach [mm] | pinned share | max abs ledger_err |
|---|---|---|---|---|---|---|---|---|
| 150–250 s | 24.8 / 24.3 | 28.3 / 27.7 | 1.76 / 1.55 | 1559 / 1564 | 136.4 / 136.0 | 62.0 / 61.0 | 0.21 / 0.19 | 1.8e-15 / 5.4e-15 |
| 250–350 s | 25.2 / 22.4 | 28.3 / 25.9 | 1.03 / 0.59 | 1623 / 1639 | 151.5 / 152.2 | 52.5 / 51.0 | 0.14 / 0.10 | 8.3e-15 / 7.9e-15 |
| 350–450 s | 24.5 / 15.4 | 27.8 / 21.1 | 0.55 / 0.17 | 1663 / 1686 | 161.5 / 170.0 | 48.1 / 46.2 | 0.12 / 0.07 | 7.0e-15 / 6.4e-15 |
| 450–550 s | 24.8 / **8.4** | 27.1 / 15.6 | 0.29 / 0.04 | 1688 / 1716 | 170.3 / **191.0** | 45.7 / 44.0 | 0.10 / **0.05** | 3.4e-15 / 8.5e-15 |

(2 mm / 1 mm. `rim.py` reconstructs `foot_z` exactly, 0.000 mm, in both runs.)

**Feet against centre:**

| t | feet depth 2 mm / 1 mm | centre depth 2 mm / 1 mm | s_c 2 mm / 1 mm | `jet_P_face` 2 mm / 1 mm |
|---|---|---|---|---|
| 250 s | 104.7 / 105.0 mm | 200 / 200 mm | 145 / 145 mm | 1442 / 1388 W |
| 350 s | 146.7 / 142.7 mm | 254 / 252 mm | 157 / 159 mm | 1202 / 1108 W |
| 450 s | 188.0 / 167.5 mm | 304 / 298 mm | 166 / 181 mm | 1052 / 933 W |
| 550 s | 228.7 / **183.0** mm | 352 / 337 mm | 173 / **204** mm | 945 / 807 W |

- **The centre keeps pace at both meshes** (352 against 337 mm at 550 s).
- **The ring stalls at 1 mm.** Its rim spall falls from 24.3 to 8.4 mm/min,
  while the 2 mm rim holds ~25 throughout.
- **So the stand-off runs away** (204 mm at 1 mm against 173), the jet decays
  over a longer distance, less power reaches the face, and the ring slows
  further. It is a feedback, and within this domain it runs away at 1 mm and
  does not at 2 mm.
- `foot_carry_cols` is 67–73 at 1 mm against 20–23 at 2 mm: three times as
  many annulus columns are standing at or above their pad height, i.e. the
  pad is riding on rock that is not firing.
- The ledger closes to ≤ 8.5e-15 everywhere, so this is not an energy-accounting
  error.

## 5. Verdict: FAIL

Both criteria fail:
- **(a)** the matched window is **−28.5 %** (limit ±5 %);
- **(b)** three of the four 100 s blocks are outside, at −12.5, −40.2 and
  −64.8 %.

**Consequence, per the pre-registered rule: stop. D2g Stage 2 (the scored
runs) is not planned on this rule.** The support rule is not the problem
D2f left it as: q = 0.9 with no clip is converged only to ~300 s, and the
pair diverges completely after that.

**Leads for the next packet** (named from the logs, not tested here):
1. **The ring stops firing at 1 mm** (rim spall 24 → 8 mm/min, pinned share
   0.19 → 0.05). Why a finer mesh cools the ring is the question. Candidates
   the logs point to: the `s ≤ 0` rule that excludes columns above the nozzle
   plane, the per-column pinned criterion, and the half-cell removal
   increment.
2. **The stand-off feedback:** s_c 204 mm at 1 mm against 173 at 2 mm, with
   `jet_P_face` down to 807 W. Whether a centre equilibrium exists at 1 mm in
   a 0.40 m domain is open; D2b saw the same runaway with the clipped h law.
3. **`foot_carry_cols` is 3× higher at 1 mm.** The pad rides on non-firing
   rock, which is where the quantile and the (retired) clip used to act.
4. **The flank dries up at 1 mm** (40–60 mm removal 0.04 against 0.29 cm³/s),
   so the hole stops widening while the centre keeps deepening.
5. **T_rec is 28 K hotter at 1 mm** with less face power: the gas leaves
   hotter because the face takes less from it.

**What stands from the earlier steps:** the D2d mesh failure (the clip) and
the D2f q-insensitivity are unaffected. What this step removes is the
assumption that the no-clip rule is mesh-converged at scored times.

**Figure:** `d2g1_blocks.png` (50 s and 100 s block gaps, with D2f Q09
overlaid).
