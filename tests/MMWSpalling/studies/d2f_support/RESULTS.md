# D2f results: the contact support rule as Meier ran it (mesh gate and q-sensitivity)

This covers the pad rule with **no clip** at `foot_pad_quantile` q = 1.0,
0.9 and 0.8, on D2e's mesh-pair geometry: the 0.40 m full quarter domain,
stop 300 s. These are **diagnostics only.** No ROP here is read against
Meier's band, and every hole number is direction only.

Full tables: `analysis.md` (`analyze.py`).
Figures: `d2f_blocks.png`, `d2f_q.png`, `d2f_rim.png`, `d2f_profile.png`.

## 0. Hashes

| item | sha256 | note |
|---|---|---|
| binary | `6d047b50bdcde619952aff907d88f89ed0cae8f9d3ed645b43ba0ac2cd647005` | the D2e build, verified at the start (17:24) and the end; no build or relink in D2f |
| `PREDICTIONS.md` | `db37f84b366cb3c3d9c5e87f904104db61e0195d085175570d2b15a3709c2012` | mtime 2026-09-18 17:27:07, written by `predict.py` before any D2f run; read-only |

## 1. Harness identity (E6: PASS)

`run.py` loads d2e's `run.py` by path and adds three switches. Each replaces
`foot_pad_quantile` and `foot_body_clearance` in d2c's `keys()` list.
Q09_2mm's command differs from D2e Pa_2mm's only in `stop_time` (300 against
200 s).

| file | rows (t ≤ 199.968 s, Pa_2mm's end) | Q09_2mm vs d2e Pa_2mm |
|---|---|---|
| thermo.dat | 2084 | byte-identical |
| removal_events.csv | 68 447 | byte-identical |
| jet_profile.csv | 1303 | byte-identical |

## 2. Runs and wall time

| run | status | t_end | wall | note |
|---|---|---|---|---|
| Q10_2mm | ok | 299.9 s | 7.8 min | three at a time |
| Q09_2mm | ok | 299.9 s | 7.8 min | three at a time |
| Q08_2mm | ok | 299.9 s | 7.8 min | three at a time |
| Q10_1mm | ok | 299.9 s | ≈ 178 min of computing | two at a time, with Q09_1mm |
| Q09_1mm | ok | 299.9 s | ≈ 178 min of computing | the optional leg, run so that outcome (b) needs no D2g condition |

- **No stall.** The longest `foot_stall_time` in any run is 8.2 s (limit 60
  s), and no watchdog fired.
- **The 1 mm pause.** The `.done` wall_s of 1135.6 min includes a
  **user-requested pause** (SIGSTOP at t ≈ 98 s, 2026-09-18 18:35 → SIGCONT
  2026-09-19 10:33, 958 min; `PAUSED.md`). The processes were suspended, not
  killed, so the runs are not restarts.
- **Machine load.** After the resume, heavy host load (load average ≈ 33,
  `kernel_task` 55 %) halved the rate for part of the time. That puts the
  computing time between 20 and about 50 wall-s per simulated s. It is not a
  clean cost figure; D2e's clean rate was 19.5.
- **Output:** 9.0 GB, git-ignored.

## 3. Gate (i): the Q10 pair (strict per-pad maximum)

| block | 2 mm ROP [m/h] | 1 mm ROP [m/h] | gap |
|---|---|---|---|
| 0–50 s | 1.556 | 1.577 | +1.3 % |
| 50–100 s | 1.543 | 1.528 | −1.0 % |
| 100–150 s | 1.474 | 1.462 | −0.8 % |
| 150–200 s | 1.437 | 1.398 | −2.7 % |
| 200–250 s | 1.426 | 1.409 | −1.2 % |
| 250–300 s | 1.502 | 1.357 | **−9.7 %** |

**The matched window 150–300 s passes:** 1.448 against 1.403 m/h, **−3.1 %**
(`score.mesh` PASS).
- At 300 s the other quantities agree: s_c 153.3 / 153.7 mm (+0.2 %), T_rec
  1618 / 1630 K (+0.7 %), centre depth 226 / 224 mm, pinned share 0.20 /
  0.18.
- Reach, T_rec, s_c and `jet_P_face` agree to 0.5 % at 100 and 200 s.
- There are no regime-4 rows, and the ledger is ≤ 7.3e-15.

**The every-block criterion fails in the last block (−9.7 %), so gate (i)
fails.** Two things make up that block, and both are reported, not explained
away:
1. **The descent is quantised.**
   - The pad rule averages three per-pad heights, so the burner steps down in
     dz/3: 0.67 mm at 2 mm, 0.33 mm at 1 mm.
   - A 50 s block holds about 20 mm of descent, so one extra 2 mm step is
     worth about 3.3 % of a block's ROP.
   - Q10_2mm took 28–29 steps (20.0 mm) in each block from 150 s, then 30
     (20.7 mm) in 250–300 s. That is 1.43 → 1.50 m/h, while the 1 mm run
     showed nothing similar.
2. **A late 1 mm slowdown, at both quantiles.**
   - After ~250 s both 1 mm runs descend about 1 mm less per block (19.0
     against 20.0 mm).
   - Rim spall at 1 mm is 23.3–23.8 mm/min, against 24.4–25.0 at 2 mm.
   - Flank removal (40–60 mm) falls faster at 1 mm: 0.75–0.93 against
     1.19–1.28 cm³/s in the last block. T_rec runs 12–17 K hotter and
     `jet_P_face` 3–6 % lower at 1 mm.
   - **This is not D2e's amplifier.** There the hotter recirculation made the
     1 mm burner faster; here the 1 mm burner is slower despite the hotter
     gas.
   - It shows up in Q09 as well, with the gap drifting +0.1 → −3.0 → −4.3 %
     over the last three blocks. At 300 s it is not known whether this
     divergence grows or settles.

**The Q09 pair (q 0.9, optional leg): converged to 300 s.**
- Matched window 150–300 s: 1.491 against 1.459 m/h, **−2.1 %**.
- Blocks: +2.0, +0.1, −1.9, +0.1, −3.0, −4.3 %, all within 5 %. The first
  three are identical to D2e P-a, as they must be.
- At 300 s: s_c −0.2 %, T_rec +1.1 %, centre depth −0.4 %, pinned share
  0.19 / 0.16.

## 4. Gate (ii): q-sensitivity (2 mm, 150–300 s)

| q | R_Q (hand) | ROP [m/h] | hand transient | s_c / T_rec at 300 s | hand s_c / T_rec | mean `foot_carry_cols` | burner / annulus mean [mm/min] | pad-setting band |
|---|---|---|---|---|---|---|---|---|
| 1.0 | 40.00 mm | **1.448** | 1.678 | 153.3 mm / 1618 K | 152.8 / 1604 | 8.6 | 24.1 / 28.1 | 38–40 (all pads) |
| 0.9 | 38.97 mm | **1.491** | 1.741 | 153.3 mm / 1624 K | 150.7 / 1614 | 22.3 | 24.9 / 28.4 | 38–40 (all pads) |
| 0.8 | 37.91 mm | **1.535** | 1.809 | 150.7 mm / 1633 K | 148.5 / 1624 | 36.6 | 25.6 / 28.6 | 38–40 (all pads) |

**Spread ±2.9 % about the mean (1.491 m/h); −0.044 m/h per +0.1 of q.** The
hand model gives ±3.8 % and −0.066. **Gate (ii) holds.**

Rim-zone recession by band (spall only, no clip), mm/min:

| q | 28–30 | 30–32 | 32–34 | 34–36 | 36–38 | 38–40 |
|---|---|---|---|---|---|---|
| 1.0 | 32.5 | 30.8 | 29.0 | 27.3 | 26.2 | 24.6 |
| 0.9 | 32.9 | 31.0 | 29.3 | 27.5 | 26.4 | 24.8 |
| 0.8 | 33.4 | 31.2 | 29.5 | 27.7 | 26.6 | 25.1 |

- The pad is a cone: the inner annulus recedes faster than the rim.
- Under any q the 38–40 mm band sets the pad height, and q only picks how
  deep into that band's slow tail the burner rides. That is why q moves the
  rate by only ~3 %.
- `rim.py` reconstructs `foot_z` exactly (0.000 mm) at every q. The
  reconstruction goes through a wrapper in `analyze.py` that sets `rim.QUANT`
  from each run's switch, so `rim.py` itself is unchanged.

**Hole at 300 s, direction only** (azimuth-mean Ø, `score.py`'s method;
`d2f_profile.png`):

| depth | q 1.0 | q 0.9 | q 0.8 | D2e B0 (q 0.9 + clip) | Fig. 8.8 visible width (finished hole) |
|---|---|---|---|---|---|
| 10 mm | 147 | 145 | 143 | 147 | 96 |
| 25 mm | 130 | 129 | 126 | 130 | 92 |
| 50 mm | 114 | 113 | 110 | 114 | 88 |
| 75 mm | 105 | 104 | 102 | 105 | 87 |
| 100 mm | 92 | 92 | 92 | 92 | 87 |
| min to the feet depth | 78.7 | 77.7 | 75.3 | 79.3 | — |

- **Walls above the nozzle plane (73–79 mm at 300 s) are final:** no gas
  reaches them.
- **Below it the wall is still being cut,** and the centre pit reaches
  226–230 mm.
- **q moves the wall by 2–4 mm,** and the clip made no difference.
- **The mouth funnel (+38 / +26 mm at 25 / 50 mm) is unchanged by the
  support rule.** It is a gas-side item for D2g, and the wellhead caveat
  applies.

## 5. PREDICTIONS rows

| id | expectation | measurement | verdict |
|---|---|---|---|
| E1 (gate i) | Q10 1 mm / 2 mm within 5 % over 150–300 s **and** every block after 50 s | window −3.1 %; blocks +1.3, −1.0, −0.8, −2.7, −1.2, **−9.7** % | **REFUTED** (the last block) |
| E2 (gate ii) | the 2 mm ROP for q 0.8 / 0.9 / 1.0 within ±5 % of the mean | 1.535 / 1.491 / 1.448 m/h, ±2.9 % | **held** |
| E3 | the ROP falls as q rises | 1.535 > 1.491 > 1.448 | **held** |
| E4 | the ROP is within ±27 % of the hand value | q 1.0 −13.7 %, q 0.9 −14.4 %, q 0.8 −15.1 % | **held** (×3) |
| E5 | s_c at 300 s is within ±20 mm of the hand value | 153.3 / 153.3 / 150.7 against 152.8 / 150.7 / 148.5 mm | **held** (×3) |
| E6 | Q09_2mm byte-identical to D2e Pa_2mm to 200 s | every row of all three files | **held** |

- **The hand model runs ≈ 14 % fast at every q,** as it did in D2c (the
  same model).
- **The hand model gets the q-slope's sign and size right:** −0.066 against
  −0.044 m/h per 0.1.
- **The hand model had no prediction for E1,** by construction: it has no
  column scatter.

## 6. Decision outcome: (b)

The pre-registered outcome for **(ii) holds, (i) fails** is **(b): D2g scores
q = 0.9**, justified as the jolt allowance.
- **Q09 is converged to 300 s** (matched window −2.1 %, every block
  ≤ 4.3 %), so the adoption carries no D2g condition.
- **The maximum's failure is a single block,** made up of the 2 mm step
  quantisation (one 0.67 mm step ≈ 3.3 %) and the late 1 mm slowdown. Its
  matched window passed (−3.1 %).
- **Q09 passes everything,** but its last blocks trend the same way (−3.0,
  −4.3 %).
- **So the difference between q = 1.0 and q = 0.9 on mesh behaviour is
  smaller than the decision makes it look.** The outcome rule is applied as
  written.

## 7. What D2g inherits

- **The support rule:** `nozzle_descent = feet`, `foot_rule = pads`,
  `foot_pad_quantile = 0.9`, `foot_body_clearance = 0`.
  - This is contact with the highest rock under each foot, within the
    axial-motion allowance.
  - q moves the rate by only ~3 % per 0.2, so the allowance barely matters.
  - The clip is retired: it has no physical counterpart and is the mesh
    seed.
- **The late 1 mm drift must be checked beyond 300 s.**
  - After ~250 s both pairs show the 1 mm run slowing slightly, with its
    flank drying faster and T_rec hotter. The gap is −4.3 % (Q09) in the last
    block.
  - D2g's scored configuration needs a mesh check over a window long enough
    to see whether that gap grows (for example 1 mm to 450–600 s at 0.40 m,
    3–4 h at the clean rate).
- **The ROP is one-sided.** The rule models the contact limit. Meier's
  operator advanced by crane and read the hook load, so his time-averaged
  0.5 m / 1383 s = **1.30 m/h** (the text says 1.5) is at or below what the
  rule gives.
- **The gas over the surface:** D2c's free-surface dilution assumes 293 K
  room air. In phase 1 the mouth was inside a sealed, water-cooled wellhead
  with confined, cooled exhaust, and the wellhead also caps any funnel. The
  mouth comparison is therefore weaker evidence than D2c/D2d treated it, and
  the funnel (+38 mm at 25 mm) is unchanged here.
- **Two phases** (647 s + 736 s): phase 2 restarted in a cooled, part-drilled
  hole under rough vacuum. The model is one continuous run at 1 atm.
- **Still open:**
  - the outer-flank gas reach;
  - the rim's three slots;
  - lifting the burner (Meier's holds and jolts shape the hole locally, for
    example the 10 cm underream);
  - P-b (the recirculation amplifier), not needed while no seed is present.
