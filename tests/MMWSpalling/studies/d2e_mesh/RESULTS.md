# D2e results: speed-up acceptance and the mesh-seed diagnosis

Controlled 2 mm / 1 mm pairs on the D2d mesh-pair geometry (0.40 m full
quarter domain). The switches are **diagnostics, not scored
configurations**, and no pair's ROP is read against Meier's band.

## 0. Hashes

| item | sha256 | note |
|---|---|---|
| `PREDICTIONS.md` | `7245201c1fb11b68c5219f026dc294db1c82ec7de72e1e21f6b74f836e564a79` | mtime 2026-09-18 15:05:21, written by `predict.py` before any pair run; read-only |
| binary before D2e (Perf) | `1f4d5a011feef2102075f20eee2996751910972f7b9d3682602da1e96d4ac89f` | relinked 09-18 11:27:47 by the user's speed-up session |
| binary for every D2e run | `6d047b50bdcde619952aff907d88f89ed0cae8f9d3ed645b43ba0ac2cd647005` | the single D2e build (rec_fixed), 09-18 14:41:45; recorded per run in `.done` |
| witness reference | `witness/ref_d2c.md5` (= `post.md5`, 2639 files) | D2c final binary; copied from the volatile scratchpad |

## 1. Perf acceptance and `max_grid_size`

**Perf: ACCEPTED.** The speed-up binary (`1f4d5a01…`, no rebuild) reproduces
**2639 / 2639 witness files byte-identical** to the D2c reference
(`witness/perf.md5` vs `witness/ref_d2c.md5`, empty diff). Every witness test
exits with rc 0:
- `jet_d2c`;
- the Meier `dev2d`;
- `robin_face`, `robin_pinned`, `robin_jet`, `jet_enthalpy`, `robin_feet`;
- `beam_void_closure`, `scalar_flaw`, `spall_event`;
- the S1 smoke.

**After the D2e build (rec_fixed)** the witness set is again **2639 / 2639
identical** (`witness/post_d2e.md5`), and `unit/jet_d2c` passes (30 checks:
26 + (k) + the 3 new aborts).

**`max_grid_size` 20 vs 60** (B0 at 2 mm, 16 ms, 30 s, 4 ranks,
`run.py --mgs-check`):

| file | rows | 20 vs 60 |
|---|---|---|
| thermo.dat | 313 | differs only in `ledger_err` (max abs 4.4e-16) |
| removal_events.csv | 7772 | byte-identical |
| jet_profile.csv | 255 | byte-identical |

- **The two agree, so the step uses 60.** The only differing column is the
  energy-ledger residual. It is itself a relative round-off quantity
  (~1e-15), and its difference, 4.4e-16, is a reordered sum over the boxes.
  Every physics column and both event logs are byte-identical, and a rerun
  reproduced this exactly.
- **B0_2mm stays at 20,** because its identity reference (R7b_2mm_dt16) ran
  at 20.
- **Wall time,** stepping only (from the t = 0 plotfile header to the last
  thermo row, 1 s file-time resolution): 38 s at 20, 35 s at 60, a ratio of
  0.92. A 30 s run cannot resolve the claimed 1.15× (1.58/1.37); the realised
  speed-up is measured on the full runs in §4.

## 2. Timestep and output policy (`run.py`)

- **dt:** a fixed step per mesh, **16 ms at 2 mm and 8 ms at 1 mm**, both
  5.66 K/step by the harness formula (limit 6 K, so no halving). The 1 mm
  8 ms step is gated by B0-1 mm (§4).
- **No steady-stop.** The stall (60 s) and bottom watchdogs are on. Stops are
  P-a at 200 s and everything else at 300 s.
- **Plotfiles.** They are written every 100 s at 2 mm. At 1 mm they are
  written only at t = 0 and at the stop (`amr.plot_int` = stop/dt).
  - **This departs from the packet's 150 s.** A 1 mm plotfile is 1.55 GB,
    and ALAMO also writes one at `stop_time`, so 150 s would write three per
    run (4.7 GB), against the packet's other limit of 3 GB per 1 mm run.
  - Two plotfiles (about 3.1 GB) is the minimum.
  - The analysis reads only `thermo.dat`, `removal_events.csv` and
    `jet_profile.csv`.
- **Keys** are `studies/d2c_steady/run.py keys()` for `MATRIX["R7b_2mm"]` /
  `["R7_1mm"]`. The pair switch replaces the scored key in place, so no key is
  given twice. Against the reference commands, B0 differs only in
  `stop_time` and `amr.plot_int`, plus, at 1 mm, `timestep` (8 ms) and
  `amr.thermo.plot_int`.

## 3. Rim tool validation on the D2d pair (before `PREDICTIONS.md`)

`rim.py` rebuilds the mask-top staircase from `removal_events.csv`, where each
event sets k = k_top − n_voided. It then applies `FeetDescent`'s pad rule.
Full output: `rim_d2d_validation.md`.

- **`foot_z` agrees exactly (worst |Δ| = 0.000 mm)** at every 10 s check
  time in R7b_2mm (56 checks) and R7_1mm (30), and no event's `k_top`
  disagrees with the running staircase (limit 0.5 mm; 09-18b matched to
  0.3 mm).
- **Rim (38–40 mm band) recession, spall + clip, mm/min, by 50 s block:**

  | block | 0–50 | 50–100 | 100–150 | 150–200 | 200–250 | 250–300 |
  |---|---|---|---|---|---|---|
  | 2 mm | 14.4 + 10.3 | 15.1 + 11.3 | 14.6 + 10.6 | 14.1 + 10.6 | 14.4 + 10.6 | 14.2 + 10.3 |
  | 1 mm | 12.4 + 13.7 | 13.8 + 15.5 | 14.8 + 16.7 | 15.7 + 17.6 | 16.7 + 19.1 | 17.4 + 19.3 |

- **Against 09-18b §2.3:**
  - it reproduces **2 mm 14 + 11 ≈ 25 mm/min, flat**, and **1 mm 12→17
    spall + 14→19 clip = 26→37**;
  - **the pad height is set by the 38–40 mm band** on all three pads in
    every block at both meshes. Inner-band columns at the pad height number
    zero, except once: 1 % of the 36–38 mm band, at 1 mm, at 250 s;
  - the 40–60 mm band removal (1 mm 1.37 → 0.04, 2 mm 1.88 → 1.22 cm³/s over
    150–300 s) reproduces §2.4.
- **Not reproduced: "86–100 % of that band sits at H".** Here 57–71 % (2 mm)
  and 66–87 % (1 mm) of the 38–40 mm annulus columns stand exactly at their
  pad height at the block ends. 09-18b does not state its definition, so this
  number is reported and not used.
- **The 0–50 s block.** The packet's 1.25 / 1.31 m/h is not reproduced by
  either the block fit (1.592 / 1.656) or the endpoint difference (1.44 /
  1.56). The burner stands still for the first ~3 s. From 50 s on the fit
  agrees with the packet to 0.01 m/h. `PREDICTIONS.md` fixes the fit
  definition.

**`rim.py` fixes after `PREDICTIONS.md` (15:20; no expectation changed).**
The first pass on the 2 mm pair runs showed three tool problems:
1. It checked `foot_z` at nominal 10 s marks against the *nearest* thermo row.
   At 16 ms those rows fall every 0.096 s, so a pad step inside that gap
   read as 0.67 mm in Pb_2mm.
2. It dropped the last block when a run ends one step short of its stop
   (299.904 s).
3. It applied the pad quantile to P-c, which runs `foot_rule = mean`.

Now the reconstruction runs at the thermo rows' own times, keeps edges
within 0.5 s of t_end (clipped), and applies the run's rule (the annulus mean
for P-c). The D2d validation tables are byte-unchanged, and `foot_z` is exact
(0.000 mm) in R7b, R7 and every D2e 2 mm run.

## 4. Runs

### B0 at 2 mm: harness identity (PASS)

`B0_2mm` (16 ms, `max_grid_size` 20, stop 300 s; 6.2 min alone) against
`studies/d2c_steady/output/R7b_2mm_dt16`:

| file | rows | result |
|---|---|---|
| thermo.dat | 3125 (every B0 row) | byte-identical |
| removal_events.csv | 93 099 (every B0 row) | byte-identical |
| jet_profile.csv | 1602 | byte-identical |

B0's files are a byte-for-byte prefix of the reference's. The first `--identity`
call reported FAIL because it compared rows to a nominal 300 s: B0 stops at
299.904 s, the last step before `stop_time`, while the reference ran on to
568 s and has a few rows in 299.92–300 s. The comparison now uses B0's own end
time. The harness was not changed in any other way, and the reference was not
touched.

### B0 at 1 mm: the step gate (PASS)

The gate was read at 16:08:43 by `gate_watch.py` from B0_1mm's first 150 s
(8 ms, `max_grid_size` 60), against R7_1mm (4 ms, the D2c binary,
`max_grid_size` 20):

| block | B0_1mm (8 ms) | R7_1mm (4 ms) | diff |
|---|---|---|---|
| 0–50 s | 1.652 m/h | 1.656 m/h | −0.22 % |
| 50–100 s | 1.759 m/h | 1.758 m/h | +0.06 % |
| 100–150 s | 1.875 m/h | 1.888 m/h | −0.67 % |
| s_c at 150 s | 120.0 mm | 120.67 mm | −0.55 % |
| T_rec at 150 s | 1533.7 K | 1533.8 K | −0.01 % |

**The 8 ms step at 1 mm is accepted** (limits 2 % / 1 %). P-a, P-b and P-c
run at 1 mm at 8 ms.

### Run log

| run | status | t_end | wall | note |
|---|---|---|---|---|
| B0_2mm | ok | 299.9 s | 6.2 min | alone, `max_grid_size` 20 |
| Pa_2mm | ok | 200.0 s | 5.5 min | three at a time |
| Pb_2mm | ok | 299.9 s | 7.7 min | three at a time |
| Pc_2mm | ok | 299.9 s | 7.7 min | three at a time |
| B0_1mm | ok | 299.9 s | 97.4 min | two at a time |
| Pa_1mm | ok | 200.0 s | 64.6 min | two at a time |
| Pb_1mm | **stopped at 79.5 s** | — | 28.4 min | the user's decision (below) |
| Pc_1mm | **not started** | — | — | the user's decision (below) |

**Stopped by the user, 2026-09-18 16:53.** After P-a, the user reread
Meier Ch. 8 with me and decided to move to a physical support rule for the
burner: it rests on the **highest rock under its feet** (pp. 209–223; see
`note_support_rule.md`). The 0.9-quantile pad rule and the clip are the
workarounds that rule replaces, so P-b and P-c at 1 mm were no longer needed.
- Pb_1mm was stopped with `run.py --kill` (binary + `plot_file` match), and
  the chain that would have started Pc_1mm was stopped before it did.
- B0_1mm was 13 s from its end and was allowed to finish.
- No other run was affected, and nothing is scored from the stopped run.

**Speed-up realised:**
- 1 mm, two runs at a time: **19.5 wall-s per simulated s**, against 54 for
  R7 in D2d (4 ms, D2c binary, `max_grid_size` 20). That is **2.8×**, the
  code speed-up, the 8 ms step and `max_grid_size` 60 together.
- 2 mm alone: 1.24 wall-s per simulated s, against about 6.7 for the D2d
  probe at 4 ms.

### Pairs with both meshes

#### B0: none (the scored configuration)

Matched window 150–300 s: 2 mm 1.483 m/h, 1 mm 2.116 m/h, **gap +42.7 %**.

| block | 2 mm ROP | 1 mm ROP | gap |
|---|---|---|---|
| 0–50 s | 1.588 | 1.652 | +4.1 % |
| 50–100 s | 1.549 | 1.759 | +13.6 % |
| 100–150 s | 1.494 | 1.875 | +25.5 % |
| 150–200 s | 1.438 | 1.995 | +38.8 % |
| 200–250 s | 1.499 | 2.131 | +42.2 % |
| 250–300 s | 1.499 | 2.203 | +46.9 % |

| quantity | t | 2 mm | 1 mm |
|---|---|---|---|
| reach [mm] | 100 s | 168.3 | 169.0 |
| reach [mm] | 200 s | 63.1 | 52.6 |
| T_rec [K] | 100 s | 1556.8 | 1552.9 |
| T_rec [K] | 200 s | 1552.8 | 1601.4 |
| s_c [mm] | 100 s | 112.0 | 109.0 |
| s_c [mm] | 200 s | 137.3 | 126.0 |
| P_face [W] | 100 s | 2747.4 | 2749.0 |
| P_face [W] | 200 s | 1641.7 | 1411.9 |

| mesh | rim spall + clip by block [mm/min] | annulus mean | burner | 40–60 mm [cm³/s] | regime-4 rows | max abs ledger_err |
|---|---|---|---|---|---|---|
| 2mm | 13.8+10.3, 15.3+11.3, 14.8+10.6, 14.0+10.3, 14.3+10.3, 14.5+11.0 | 29.5, 31.2, 29.3, 28.0, 28.2, 28.4 | 26.5, 25.8, 24.9, 24.0, 25.0, 25.0 | 1.24, 1.96, 1.93, 1.90, 1.70, 1.26 | 930 | 9.5e-15 |
| 1mm | 12.3+13.7, 14.0+15.4, 14.7+16.5, 15.6+17.8, 16.7+18.6, 17.2+19.5 | 31.6, 32.1, 32.3, 33.8, 35.7, 36.9 | 27.5, 29.3, 31.3, 33.3, 35.5, 36.7 | 1.41, 1.94, 1.88, 1.40, 0.50, 0.04 | 10078 | 5.0e-15 |

`score.mesh`:

Matched window 150–300 s (t_end 300 s at 1 mm, 300 s at 2 mm; steady windows are NOT required, D2d decision 2).

| quantity | 2 mm | 1 mm | 1 mm vs 2 mm |
|---|---|---|---|
| burner ROP over the window | 1.483 m/h | 2.116 m/h | +42.7 % |
| s_c at the window end | 152.0 mm | 133.0 mm | -12.5 % |
| T_rec at the window end | 1621 K | 1686 K | +4.0 % |
| centre depth at the window end | 226 mm | 243 mm | +7.5 % |
| pinned share | 0.19 | 0.12 | -0.07 |

**Mesh verdict: FAIL** (burner ROP within 5 %).

#### P-a: foot_body_clearance = 0

Matched window 150–200 s: 2 mm 1.455 m/h, 1 mm 1.456 m/h, **gap +0.1 %**.

| block | 2 mm ROP | 1 mm ROP | gap |
|---|---|---|---|
| 0–50 s | 1.608 | 1.640 | +2.0 % |
| 50–100 s | 1.574 | 1.576 | +0.1 % |
| 100–150 s | 1.526 | 1.498 | -1.9 % |
| 150–200 s | 1.455 | 1.456 | +0.1 % |

| quantity | t | 2 mm | 1 mm |
|---|---|---|---|
| reach [mm] | 100 s | 168.3 | 169.0 |
| T_rec [K] | 100 s | 1556.0 | 1555.1 |
| s_c [mm] | 100 s | 112.0 | 112.0 |
| P_face [W] | 100 s | 2750.7 | 2745.1 |

| mesh | rim spall + clip by block [mm/min] | annulus mean | burner | 40–60 mm [cm³/s] | regime-4 rows | max abs ledger_err |
|---|---|---|---|---|---|---|
| 2mm | 23.7+0.0, 26.9+0.0, 25.5+0.0, 24.5+0.0 | 29.5, 31.5, 29.2, 28.2 | 26.8, 26.2, 25.4, 24.3 | 1.24, 1.97, 1.94, 1.89 | 0 | 3.0e-15 |
| 1mm | 25.4+0.0, 26.3+0.0, 25.0+0.0, 24.2+0.0 | 31.4, 30.7, 28.9, 27.7 | 27.3, 26.3, 25.0, 24.3 | 1.41, 1.96, 1.93, 1.82 | 0 | 5.1e-15 |

`score.mesh`:

Matched window 150–200 s (t_end 200 s at 1 mm, 200 s at 2 mm; steady windows are NOT required, D2d decision 2).

| quantity | 2 mm | 1 mm | 1 mm vs 2 mm |
|---|---|---|---|
| burner ROP over the window | 1.455 m/h | 1.456 m/h | +0.1 % |
| s_c at the window end | 136.0 mm | 136.3 mm | +0.2 % |
| T_rec at the window end | 1557 K | 1562 K | +0.4 % |
| centre depth at the window end | 170 mm | 171 mm | +0.6 % |
| pinned share | 0.23 | 0.22 | -0.01 |

**Mesh verdict: PASS** (burner ROP within 5 %).


### Two-mm-only runs (no 1 mm partner; reported, no verdict)

| block ROP [m/h] | 0–50 | 50–100 | 100–150 | 150–200 | 200–250 | 250–300 | T_rec / s_c at 300 s |
|---|---|---|---|---|---|---|---|
| B0_2mm | 1.588 | 1.549 | 1.494 | 1.438 | 1.499 | 1.499 | 1621 K / 152.0 mm |
| Pb_2mm (rec_fixed 1550 K) | 1.581 | 1.540 | 1.496 | 1.465 | 1.443 | 1.431 | 1550 K / 152.0 mm |
| Pc_2mm (mean rule, no clip) | 1.972 | 1.868 | 1.775 | 1.766 | 1.777 | 1.772 | 1658 K / 142.2 mm |

Pb_1mm's 0–50 s block (1.644 m/h, before it was stopped) is +4.0 % over
Pb_2mm, as B0's is. At 0–50 s the recirculation is still at its start value,
so this says nothing about the amplifier.

## 5. PREDICTIONS rows

| pair | expectation | measurement | verdict |
|---|---|---|---|
| B0 | 2 mm byte-identical to R7b_2mm_dt16 to 300 s | every B0 row of thermo, removal events and jet profile | **held** |
| B0 | 1 mm step gate: blocks within 2 %, s_c and T_rec at 150 s within 1 % | −0.22 / +0.06 / −0.67 %; −0.55 %, −0.01 % | **held** |
| B0 | gap over 150–300 s ≥ +30 % | +42.7 % (D2d +42.8 %) | **held** |
| P-a | early gap (0–50, 50–100 s) ≤ 3 % | +2.0 %, +0.1 % | **held** |
| P-a | burner slower than B0 at 2 mm | *faster* in every block: 1.608 vs 1.588, 1.574 vs 1.549, 1.526 vs 1.494, 1.455 vs 1.438 m/h | **REFUTED** |
| P-a | rim = spall only | 0 regime-4 rows at both meshes | **held** |
| P-b | gap(250–300) ≤ gap(50–100) + 3 points | Pb_1mm stopped at 79.5 s | **not tested** |
| P-c | 150–300 s gap ≤ 5 %; burner follows the annulus mean | Pc_1mm not run | **not tested** |

**The P-a refutation, as it stands.** With the clip off, the rim recedes by
spall alone at 24–27 mm/min. That is as fast as spall plus clip did at 2 mm
(about 25). The clip replaced spall on the rim rather than adding to it: a
clipped cell restarts cold and spalls later. The expectation assumed the two
add, and the numbers show they do not.

## 6. Decision-rule verdict

- **The seed is identified: the clip.** P-a removes the early gap (+2.0 %,
  +0.1 %) and the whole mesh gap, **+0.1 % over the matched window
  150–200 s**, against +38.8 % for B0 in the same block. Every paired
  quantity agrees to 0.6 % or better: s_c 136.0 / 136.3 mm, T_rec
  1557 / 1562 K, centre depth, pinned share and reach.
  - The rim rate converges: spall at 24–27 mm/min at both meshes, and zero
    clip rows.
  - The flank band removal matches: 1.2–2.0 cm³/s at both.
- **The amplifier half of the rule is not tested,** because P-b was stopped.
  So "understood" in the pre-registered sense (P-a removes the early gap
  **and** P-b stops the growth) is **not reached**.
  - What P-a does show: **without the seed the gap does not grow at all** over
    200 s. The recirculation loop has nothing to amplify.
  - Whether it would amplify another seed of the same size is P-b's question,
    and it stays open.
- **The named fallbacks** (the pin idle timeout, the half-cell increment) are
  **not needed**. P-a shrank the early gap.

## 7. What the pairs support for D2f (not a design)

- **The burner-rate closure is the mesh-sensitive element**: the 0.9
  quantile plus the whole-cell clip.
  - The quantile alone, with no clip, converges over 200 s: P-a.
  - The clip alone is what diverges: B0 against P-a.
- **The user's decision (2026-09-18):** replace the rule with the physical
  support Meier describes. The burner rests on the **highest rock under each
  foot** and descends as that rock is removed.
  - The operator jolts keep it free (pp. 219–221), and the clearance is set by
    the advance rate (p. 223).
  - The 0.9 quantile (a guard against one never-fired 2 mm column) and the
    clip (the cleanup of rock the quantile leaves inside the feet) are the
    workarounds it retires.
- **What the pairs imply for that rule:**
  1. **A per-pad maximum is set by the single slowest column.** Its mesh
     behaviour is not tested here. P-a converges because the 0.9 quantile
     averages over a tail.
  2. **Rock under the steel must be shielded** (09-18b §4): no jet flux from
     above. It then goes only by undercutting from the bore side or by the
     push, and **that rule must be shown converged by this pair method
     before any scored run.**
  3. The mouth result is unaffected: the flank removal is the same at both
     meshes and in every pair.
- **Every Meier number from D2b to D2d was made under the quantile-plus-clip
  rule, and is superseded once D2f changes it.** The gas-side closures and
  every non-Meier validation are unaffected.
