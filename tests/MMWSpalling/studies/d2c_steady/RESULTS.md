# D2c results: steady-state and no-funnel closures (pre-registration; no scored runs)

## 0. Pre-registration record

- `PREDICTIONS.md` sha256 `a28f19875a7240f73b40eef0d2d70eda050e6147bf3e7d0d9e80a4274b5bd651`
- `PREDICTIONS.md` mtime 2026-09-17 13:54:01
- Written by `predict.py` before any D2d run and before the D2c cost probes.
- Scored configuration: ACTIVE_STEP §Scored configuration, unchanged.

<!-- NOZZLE -->

## 1. Nozzle state (`nozzle.py`, deliverable 7; output verbatim)

```
Laval state (throat Ø 6, exit Ø 7.5 mm, mdot 0.01513 kg/s, R 290, p_amb 1.013 bar)
 T0 [K] gamma p0 [bar]    M_e p_e [bar] T_e [K] u_e [m/s]   J [N] (p_e-p_a)A_e [N]
   1600  1.25    5.539  1.838     0.952    1125    1173.6  17.488          -0.2693
   1600  1.30    5.463  1.860     0.892    1053    1172.2  17.203          -0.5329
   1600  1.35    5.391  1.883     0.838     987    1170.7  16.938          -0.7747
   1750  1.25    5.793  1.838     0.996    1231    1227.4  18.494          -0.0766
   1750  1.30    5.713  1.860     0.933    1152    1226.0  18.196          -0.3523
   1750  1.35    5.638  1.883     0.876    1080    1224.3  17.919          -0.6051
   1900  1.25    6.036  1.838     1.037    1336    1278.9  19.458           0.1081
   1900  1.30    5.953  1.860     0.972    1251    1277.4  19.148          -0.1792
   1900  1.35    5.875  1.883     0.913    1173    1275.7  18.859          -0.4426

Momentum diameter D_e = 2 mdot/sqrt(pi rho J) (gamma = 1.3) and core lengths
 T0 [K] T_surr [K] D_e [mm] 5 D_e [mm] 8 D_e [mm]
   1600     293.15    3.751      18.76      30.01
   1600     800.00    6.197      30.98      49.58
   1600    1500.00    8.485      42.43      67.88
   1750     293.15    3.647      18.24      29.18
   1750     800.00    6.025      30.13      48.20
   1750    1500.00    8.251      41.25      66.00
   1900     293.15    3.556      17.78      28.44
   1900     800.00    5.874      29.37      46.99
   1900    1500.00    8.043      40.21      64.34

jet_De_ref (full jet, at jet_T_ent = 293.15 K) for ALAMO:
  T_nozzle 1600 K: jet_De_ref = 0.0037512240761077564 m (gamma 1.25: 3.721 mm, 1.35: 3.780 mm)
  T_nozzle 1750 K: jet_De_ref = 0.003647390490008811 m (gamma 1.25: 3.618 mm, 1.35: 3.675 mm)
  T_nozzle 1900 K: jet_De_ref = 0.0035555966839916735 m (gamma 1.25: 3.527 mm, 1.35: 3.583 mm)

Check against the ACTIVE_STEP §Sources table (1900 K, gamma 1.3):
  p0 [bar]                 5.953  (packet 5.95 ± 0.01)  ok
  M_e                      1.860  (packet 1.86 ± 0.005)  ok
  p_e [bar]                0.972  (packet 0.97 ± 0.01)  ok
  T_e [K]               1250.748  (packet 1250 ± 10)  ok
  u_e [m/s]             1277.414  (packet 1277 ± 5)  ok
  mdot*u_e [N]            19.327  (packet 19.3 ± 0.1)  ok
  D_e 293 K [mm]           3.555  (packet 3.55 ± 0.06)  ok
  D_e 800 K [mm]           5.874  (packet 5.8 ± 0.1)  ok
  D_e 1500 K [mm]          8.043  (packet 8 ± 0.06)  ok
  8 D_e 293 K [mm]        28.437  (packet 28.5 ± 0.6)  ok
  8 D_e 1500 K [mm]       64.343  (packet 64 ± 0.6)  ok
  J used (incl. pressure thrust -0.179 N) = 19.148 N, -0.93 % vs mdot*u_e -> D_e +0.47 %
  all reproduced
```

Notes:
- **J includes the pressure thrust**, as the packet asks: J = 19.15 N at
  1900 K. The packet's 19.3 N is ṁ·u_e alone (19.33 N). The difference
  (−0.9 %) makes D_e larger by 0.47 %.
- **D_e at 800 K** is 5.87 mm (5.85 with ṁ·u_e); the packet's 5.8 is rounded.
- **Every other packet value is reproduced,** including 8·D_e = 28.4 mm in
  ambient air and 64 mm in 1500 K exhaust.
- **`jet_De_ref` per T_nozzle** (γ 1.3; used by `run.py` and `predict.py`):
  1600 K 3.751 mm, 1750 K 3.647 mm, 1900 K 3.556 mm.
- **γ sensitivity:** ±0.05 in γ moves D_e by about ±0.8 %.

<!-- DRYRUN -->

## 2. Dry run: D2c scoring on D2b's B_JM_A2 (deliverable 10)

The D2b scored run, re-scored with the D2c rules. No D2b table is rewritten; the D2b
configuration has no clearance and no free-surface dilution.

### B_JM_A2 (bottom at 341.6 s)

| criterion | target | simulation | verdict |
|---|---|---|---|
| steady window | ≥ 150 s, |ds_c/dt| < 0.02 mm/s, 50 s sub-fits ±10 % | none (fallback 171–342 s; s_c drift +0.39 mm/s; sub-fits 1.39, 1.39, 1.34 m/h) | FAIL: no steady window |
| ROP | [1.04, 1.92] m/h | 1.36 m/h (ring columns 1.52) | FAIL (no window) |
| whole volume rate | [1.98, 2.96] cm³/s | 5.68 cm³/s (inside Ø 93, reported: 3.27) | FAIL (no window) |
| mechanical share | ≤ 5 % (above = finding) | no clearance in this run | — |
| wall profile (mouth) | 0 ≤ Ø − width ≤ 10 mm at 25 and 50 mm (Fig. 8.8) | z 25: Ø 164.5 vs 92 (+72.5); z 50: Ø 130.8 vs 88 (+42.8); RMS over 3 rows inside 71 mm: 83.9 mm; rows below the width: none | FAIL (no window) |
| depth-mean Ø (drilled depth) | 85–93 mm | 161.1 mm over 0–71 mm (max-r Ø 163.0) | FAIL (no window) |
| minimum Ø (drilled depth) | ≥ 80 mm; 76–80 unresolved at 2 mm | 117.1 mm at 70 mm (max-r Ø 118.8) | FAIL (no window) |
| ΔT_fire | 500–560 K | 528 K (T_fire p10–p90 815–828 K) | FAIL (no window) |
| removal power / R | vs 2.83 kW; R ≤ 1.05 | 6.52 kW (2.30× Meier); R = 0.828 (face power ×4 7.87 kW) | FAIL (no window) |
| ledger | round-off | 1.8e-14 | PASS |
| mesh | full domain 1 mm vs 2 mm, ROP within 5 % | score.py --mesh (D2d run 7) | — |

Profile at the window end (Ø in mm; model azimuth-mean / max-r / Fig. 8.8 width):

| depth | 10 mm | 25 mm | 50 mm | 75 mm | 100 mm | 150 mm | 200 mm | 250 mm |
|---|---|---|---|---|---|---|---|---|
| model (az-mean) | 214 | 165 | 131 | — | — | — | — | — |
| model (max r) | 216 | 165 | 132 | — | — | — | — | — |
| Fig. 8.8 | 96 | 92 | 88 | 87 | 87 | 86 | 85 | 84 |

Reading:
- **Steady window:** none, as expected. The D2b run hit the bottom with s_c still drifting (+0.39 mm/s), so every scored row reads FAIL (no window). The numbers are shown over the fallback window.
- **Whole-excavation volume:** 5.68 cm³/s, 1.9× the top of the band. D2b's in-Ø 93 score was 3.27 cm³/s and is now only reported.
- **Profile metric:** the mouth reads 165 / 131 mm at 25 / 50 mm against 92 / 88 mm (Fig. 8.8), and the RMS over the drilled depth is 83.9 mm. This is the D2b funnel, the target of the D2c free-surface closure.
- **Drilled-depth Ø:** depth-mean 161.1 mm, minimum 117.1 mm, over 0–71 mm (nozzle-plane depth at the window end). The minimum applies the 2 mm rule.
- **Profile CSV:** `output/B_JM_A2_fig88_profile.csv` in this study; no D2b output was changed.

<!-- PROBES -->

## 3. Cost probes and the D2d schedule (deliverable 12)

The probes were run alone, after `PREDICTIONS.md` was hashed (13:54), using the
scored configuration (`run.py keys()`, 0.40 m domain). They are not scored.
The ACTIVE_STEP mtime (12:46:45) was re-checked first and was unchanged.

| probe | cells | dt | simulated | wall | cell-steps/s | wall s per simulated s |
|---|---|---|---|---|---|---|
| probe_2mm | 60² × 200 | 4 ms | 29.9 s | 200 s | 2.69e+07 | 6.7 |
| probe_1mm | 120² × 400 | 4 ms | 9.9 s | 550 s | 2.59e+07 | 55.6 |

Projection: to each run's expected end (the hand model's bottom time at 0.40 m, PREDICTIONS.md; the cap if the hand has no bottom). Alone = the probe rate; three at a time = × 0.72 (D2b Stage B slowest / lone).

| run | dz | dt | expected end (hand) | alone | three at a time | to the cap alone |
|---|---|---|---|---|---|---|
| R1_scored | 2 mm | 4 ms | 498 s (bottom) | 0.93 h | 1.29 h | 2.23 h |
| R2_clamp | 2 mm | 4 ms | 354 s (bottom) | 0.66 h | 0.91 h | 2.23 h |
| R3_n05 | 2 mm | 4 ms | 433 s (bottom) | 0.81 h | 1.12 h | 2.23 h |
| R4_1600 | 2 mm | 4 ms | 716 s (bottom) | 1.33 h | 1.85 h | 2.23 h |
| R5_1750 | 2 mm | 4 ms | 588 s (bottom) | 1.09 h | 1.52 h | 2.23 h |
| R6_d2b_n1 | 2 mm | 4 ms | 557 s (bottom) | 1.04 h | 1.44 h | 2.23 h |
| R7_1mm | 1 mm | 4 ms | 498 s (bottom) | 7.69 h | 10.67 h | 18.53 h |
| Oa_core5 | 2 mm | 4 ms | 518 s (bottom) | 0.96 h | 1.34 h | 2.23 h |
| Ob_ricou | 2 mm | 4 ms | 532 s (bottom) | 0.99 h | 1.37 h | 2.23 h |

Run 7 (1 mm, full domain): 907 s simulated in 14 h alone.

**Probe sanity** (not a score; `jet_*` identities from thermo; ledger 2e-15):
- **State at 30 s (2 mm):**
  - s_c 74 mm, T_stag 1870 K, T_rec (= `jet_T_mouth`) 1604 K, D_e 8.32 mm;
  - the hand transient gives 74 mm, 1865 K, 1578 K and 8.25 mm at the same time;
  - descent 12.0 mm (2 mm steps), against 14.4 mm in the hand model.
- **Free surface is active:** 3,284 → 3,176 marching free-surface columns,
  with m_exit/m_stag = 4.25 → 3.66. The exhaust over the free surface is at
  519–585 K, while the recirculated mouth gas is at about 1600 K.
- **Clearance fired:** 90 regime-4 cells (7.2e-7 m³ in the sector) in 30 s.
  `jet_P_neg` stayed 0: no rock was hotter than its gas.
- **1 mm at 10 s:** s_c 59 mm (2 mm: 58), T_rec 1626 K.
- **Cost:** 2.7e7 and 2.6e7 cell-steps/s. D2b's lone probes gave 1.9e7 and
  2.4e7.

**Proposed D2d schedule** (15 cores; three 4-rank runs at a time):
- **Night 1**, all 2 mm; each batch takes about as long as its slowest run.
  - **Batch 1 (≈ 1.5 h):** R1 scored, R2 clamp, R6 (D2b closures + n = 1).
  - **Batch 2 (≈ 1.9 h):** R4 1600 K, R5 1750 K, Oa core 5.
  - **Batch 3 (≈ 1.4 h):** Ob Ricou.
  - R3 (n = 0.5) is **not run**: PREDICTIONS P1 needs a ≥ 1.55 m domain.
  - Total ≈ 5 h, with the 2 mm cap (1200 s, 2.2 h alone) as the backstop.
- **Night 2: R7 (1 mm, full domain) alone.**
  - The hand model puts the bottom stop at 498 s, about 7.7 h. The hand runs
    about 13 % fast, so the simulated bottom may come nearer 560 s, about
    8.7 h.
  - **Proposed cap: 750 s simulated** (≈ 11.6 h alone), which keeps it
    within 14 h.
- **Run 7 cannot reach a steady window within about 14 h alone.** The hand
  window needs a ≥ 0.65 m domain (1.6× the cells) and about 950 s simulated.
  That is ≈ 950 × 55.6 × 1.625 s ≈ 24 h alone. In a 0.40 m domain, 14 h
  covers 907 s, but the centre reaches the bottom at about 500 s. **D2d
  decides.** This packet does not trim the domain.
- **The full-domain mesh criterion needs both runs on steady windows.** Per
  P1, neither R1 nor R7 has one at 0.40 m, so the mesh check as defined
  cannot pass at 0.40 m. This is also a D2d decision.

<!-- NOTES -->

## 4. Notes for the reviewer and D2d

- **Clearance order (verified in `Advance`).** Removal runs first, then the
  thermal advance, which calls `BuildFlameColumns` → `FeetDescent`. The cut
  armed at step n acts before PASS 1 of step n + 1, so the thermal removal of
  that step sees the cleared surface. Unit test (h) checks this order: foot_z
  and the pad heights are rebuilt from the event log, and each cut lands on
  the pad height of the tops before it.
- **Literature check of `jet_fs_aspect` / `jet_fs_exponent`: not done in
  D2c** beyond the packet's citations. The values are unchanged (aspect 1,
  exponent 1), and the hand sensitivities are in PREDICTIONS P5. On the mouth
  they give: aspect 0.5 → Ø25 142 mm, aspect 2 → 126 mm, exponent 0.8 →
  138 mm; scored 134 mm.
- **Fig. 8.8 was not re-digitised,** as the packet asked, and no overlay
  was rendered. `score.py` reads the CSV as given, and no disagreement is
  recorded.
- **`score.py` interpretation of the profile rule.**
  - The mouth passes iff 0 ≤ Ø − width ≤ 10 mm at 25 and 50 mm.
  - Any CSV row inside the drilled depth with Ø below the visible width
    fails the profile, because the width is a lower bound.
  - Ø(z) is the azimuth-mean (9 sectors) of the per-sector max r of
    columns deeper than z.
  - The drilled depth is the nozzle-plane depth at the window end. For the
    D2b dry run that is 71 mm, so the drilled-depth Ø metrics there cover
    the mouth only.
- **The D2b-configuration hand model runs fast and hot** against B_JM_A2:
  ROP +13 %, centre depth +11 %, s_c +6 %, T_rec +1 %. The P6 bands (±27 %
  on rates) are twice that ROP error.
- **What the predictions tell D2d** (pre-registered, not results):
  1. **Domain height.** At 0.40 m no steady window is predicted before the
     bottom stop; a steady window needs ≥ 0.65 m.
  2. **ROP.** The steady ring ROP is 1.93 m/h, just above the band. The
     transient ROP from 150 s to the bottom is 1.79 m/h.
  3. **Mouth.** It stays too wide (134 / 118 mm at 25 / 50 mm), so
     free-surface dilution alone does not close the funnel.
  4. **n.** The ring barely depends on n (+1.2 % for n = 0.5, −5.7 % for
     the clamp).
  5. **T_nozzle.** 1.3–1.6 m/h would need T_nozzle from below 1600 K up to ≈ 1730 K (steady),
     conditional on the scored closures.
- **Hand-model limits** (see the `handmodel.py` docstring):
  - axisymmetric to 0.12 m, whereas the quarter's corners reach 0.17 m (the
    probe's m_exit/m_stag is 3.7–4.3, against the hand's 2.4);
  - the face is at T_fire everywhere;
  - a single clearance band, so the mechanical share (0.04 %) is a lower
    bound.

## 5. D2d pre-run record (written before the first campaign run)

Timestamp **2026-09-17 21:46:49**. Nothing below may change once a run has started.

### Frozen-input checks
- `PREDICTIONS.md` sha256 `a28f19875a7240f73b40eef0d2d70eda050e6147bf3e7d0d9e80a4274b5bd651`
  — **matches** the hash in ACTIVE_STEP (`a28f1987...5bd651`).
- Binary `bin/mmwspalling-3d-g++` sha256
  `4b4f9c51f71f5ee8ac7522c54082d7656b1076f45966bbc65a7a46135178cdd4`,
  built 2026-09-17 13:42:57 (the D2c final build). **No rebuild during the
  campaign.**
- `validation/meier/sp_meier_pilot/input_feet_d2c`: regenerated by
  `run.py --write-input` and **byte-identical** (it is pinned to the frozen
  D2c height and cap; the campaign heights live in `MATRIX`).
- Harness after the amendments: `run.py` sha256 `240a09158f965faf52e45f9290b19098d77e34ebb845aff42a55bfd0677fa8cf`,
  `score.py` sha256 `874e7fee611cf3fd08db6792ebd283b7660d23883a75f34e656375d36d131f1a`.

### Decisions (ACTIVE_STEP §Context 1-6, verbatim in substance)
1. **Domain height 0.70 m** for every 2 mm matrix run except the mesh
   partner. The hand model needs 0.65 m and runs 13 % fast, so 0.70 m leaves
   margin. Window (a) stays measurable, from 150 s to the time the centre is
   0.36 m deep, which is where a 0.40 m domain would have stopped.
2. **The mesh check is a matched transient window on a 0.40 m pair, full
   domain.** The pair is `R7_1mm` and the new `R7b_2mm`: same keys, same
   height, only dz differs. The criterion is the burner ROP over the same
   window [150 s, min(t_end)] within 5 %. A steady 1 mm run would need about
   24 h, and mesh convergence does not require a steady state. Also reported:
   s_c and T_rec at common times, the pinned share, and the centre depth at
   the window end. **This replaces the D2c steady-window mesh criterion.**
3. **Run 7 is capped at 700 s**, and `R7b_2mm` gets the same cap.
   `CAP_1MM = 1200` is replaced by 700.
4. **The minimum Ø is scored to the feet depth** (`foot_z` at the window
   end), not only to the nozzle plane; both are reported. The depth-mean Ø
   and the mouth keep the nozzle-plane range.
5. **Run 6 is `R6_d2bjet_n1`:** "D2b jet closures + n = 1"; only the momentum
   diameter, the core length and the free surface revert. Exponential, count
   and clearance stay on.
6. **The T sweep stays at 1600, 1750 and 1900 K.** The 1.3 m/h end is
   reported as "at or below 1600 K", open-ended. `Oc_1450` is optional and
   only extends the sweep; it does not change the scored reading.

### Amended scoring definitions (`score.py`)
- **(i)** The minimum Ø is taken over [0, feet depth at the window end]; the
  nozzle-plane minimum is also printed. On the D2b dry run this changes the
  one row and recovers D2b's 78.2 mm minimum (117.1 mm over the drilled
  depth). Every other row of the dry-run table is unchanged (diffed).
- **(ii)** `--mesh` uses the matched fixed window [150 s, min(t_end)].
- **(iii)** New `--window-a`: [150 s, t_a], with t_a the first time the
  centre depth `lz − (nozzle_z − jet_s_c)` reaches 0.36 m (t_a = t_end,
  flagged, if never reached). The PREDICTIONS P6 bands are stated on this
  window, so `--p6` reads it.
- `--p6` evaluates every P6 line mechanically: held / refuted / not tested.

### Matrix and watchdogs
- `R3_n05` is **not run** (P1: it needs ≥ 1.55 m). It is reported from the
  hand model only.
- Watchdogs: stalled (60 s without descent), bottom (centre ≤ 40 mm above the
  domain bottom), steady-stop (50 s after a steady window is found), and the
  per-run cap (2 mm 1300 s, mesh pair 700 s).
- Plotfiles: 200 s at 2 mm, 350 s at 1 mm. R7 writes about 4.7 GB (3
  plotfiles of 1.55 GB), and each 0.70 m run about 2.5 GB.

### Parse smokes (0.5 s simulated, before the campaign)
- `R1_scored` at 0.70 m: rc 0, 60 × 60 × 350 cells, dt 4 ms; the thermo
  header is byte-identical to the D2c 2 mm probe's (47 columns).
- `R7b_2mm` at 0.40 m: rc 0, 60 × 60 × 200 cells, dt 4 ms.

### Projection and planned order

| probe | cells | dt | simulated | wall | cell-steps/s | wall s per simulated s |
|---|---|---|---|---|---|---|
| probe_2mm | 60² × 200 | 4 ms | 29.9 s | 200 s | 2.69e+07 | 6.7 |
| probe_1mm | 120² × 400 | 4 ms | 9.9 s | 550 s | 2.59e+07 | 55.6 |

Projection: to each run's expected end, the earliest of the hand model's bottom time at that run's height, its steady-stop and its cap. Alone = the probe rate scaled by the cell count; three at a time = × 0.72 (D2b Stage B slowest / lone).

| run | dz, height | dt | expected end (hand) | alone | three at a time | to the cap alone |
|---|---|---|---|---|---|---|
| R1_scored | 2 mm, 0.7 m | 4 ms | 950 s (steady-stop) | 3.09 h | 4.29 h | 4.23 h |
| R2_clamp | 2 mm, 0.7 m | 4 ms | 666 s (bottom) | 2.16 h | 3.01 h | 4.23 h |
| R4_1600 | 2 mm, 0.7 m | 4 ms | 1115 s (steady-stop) | 3.63 h | 5.04 h | 4.23 h |
| R5_1750 | 2 mm, 0.7 m | 4 ms | 1025 s (steady-stop) | 3.33 h | 4.63 h | 4.23 h |
| R6_d2bjet_n1 | 2 mm, 0.7 m | 4 ms | 1000 s (steady-stop) | 3.25 h | 4.52 h | 4.23 h |
| R7_1mm | 1 mm, 0.4 m | 4 ms | 498 s (bottom) | 7.69 h | 10.67 h | 10.81 h |
| R7b_2mm | 2 mm, 0.4 m | 4 ms | 498 s (bottom) | 0.93 h | 1.29 h | 1.30 h |
| Oa_core5 | 2 mm, 0.7 m | 4 ms | 965 s (steady-stop) | 3.14 h | 4.36 h | 4.23 h |
| Ob_ricou | 2 mm, 0.7 m | 4 ms | 980 s (steady-stop) | 3.19 h | 4.43 h | 4.23 h |
| Oc_1450 | 2 mm, 0.7 m | 4 ms | 1235 s (steady-stop) | 4.02 h | 5.58 h | 4.23 h |

Run 7 (1 mm, 0.4 m): 907 s simulated in 14 h alone; its cap is 700 s (decision 3).

Planned order (stated before launching; it may be reordered for machine load):
1. **batch A** — R1_scored, R4_1600, R5_1750 (3 jobs);
2. **batch B** — R2_clamp, R6_d2bjet_n1, R7b_2mm (3 jobs);
3. **R7_1mm alone**;
4. the optional runs (Oa, Ob, Oc) only if the budget allows.

**Addendum 2026-09-18 (D2e Goal 7; post-start amendment).** `score.py` was
edited on 2026-09-18 at 08:32, after R1–R6 had finished. The change was in
`sweep()` only: the single-basis rule, so that no fit mixes a steady ROP with
transient ones. Its sha256 went from `874e7fee…` to
`2de24904572a4aec25165f57492cc97365a40bcf33e81f4ca9b2ef8d17f2ee88`. No run
used `score.py`, and every other function is unchanged.

<!-- D2D -->

## 6. Pre-registration check (PREDICTIONS.md P6)

`score.py --p6` on R1, R2, R4, R5 and R6. The bands are the hashed P6 bands and
are evaluated on window (a), [150 s, the time the centre reaches 0.36 m].

| P6 item | prediction | measurement | verdict |
|---|---|---|---|
| P1 scored | no steady window before the 0.40 m bottom depth; s_c settles in 160–200 mm | steady window: yes; mean s_c 156 mm, drift +0.103 mm/s | held |
| P1 clamp | no centre equilibrium (|ds_c/dt| stays ≥ 0.02 mm/s) | drift +0.476 mm/s, steady window no | held |
| P1 n = 0.5 | needs a ≥ 1.55 m domain | R3_n05 not run (decision) | not tested |
| P2 scored ROP | 1.31–2.27 m/h | 1.50 m/h (window (a) 150–566 s) | held |
| P2 clamp ROP | 1.21–2.1 m/h | 1.41 m/h (window (a) 150–384 s) | held |
| P3 1600 K ROP | 0.91–1.57 m/h | 0.75 m/h (window (a) 150–884 s) | REFUTED |
| P3 1750 K ROP | 1.11–1.91 m/h | 1.29 m/h (window (a) 150–682 s) | held |
| P3 ordering | ROP increases with T_nozzle | 0.75 → 1.29 → 1.50 m/h (1600 → 1750 → 1900 K) | held |
| P4 mouth | stays too wide: 134 / 118 mm ±20 at 25 / 50 mm | 130 / 114 mm vs Fig. 8.8 92 / 88 | held |
| P4 depth-mean Ø | 91–131 mm | 106 mm | held |
| P4 whole volume | 3.17–5.48 cm³/s | 3.47 cm³/s (window (a) 150–566 s) | held |
| P4 mechanical share | ≤ 5 % | 2.91 % | held |
| P5 run 6 mouth | wider than run 1 (hand 174 vs 134 mm at 25 mm) | 167 vs 130 mm | held |
| P5 run 6 ROP | 1.71 m/h ±27 % | 1.46 m/h | held |
| P5 core 5 / Ricou ordering | Ricou < core 5 < scored | optional runs not available | not tested |
| P5 aspect / exponent | hand model only | no D2d run (by design) | not tested |

**11 held, 1 refuted, 4 not tested.**

### The refutation: P3 at 1600 K

The hand model predicted 0.91–1.57 m/h; the run gives **0.75 m/h**, 18 % below
the band. It is not explained away, and the run is not repeated.

- **What the run shows.** R4 never reaches a steady window, and its 50 s
  sub-fits swing between 0.38 and 0.88 m/h over 442–884 s. The descent is
  intermittent: the ring fires, stalls and fires again.
- **Hypothesis for D3** (stated, not tested here). The hand model assumes a
  face that fires continuously at T_fire, so its ROP is the pinned closed form
  at the local gas temperature. At 1600 K the gas over the ring is only about
  150–200 K above the firing threshold, so parts of the face fall below it and
  wait for conduction to catch up. That intermittency lowers the mean rate and
  is outside the hand model by construction. The same mechanism would explain
  why R5 (1750 K) is still transient at the 1300 s cap while R1 (1900 K)
  settles by 470 s.
- The **ordering** prediction still holds: 0.75 → 1.29 → 1.50 m/h.

### Window (a) table (the P6 basis)

| run | window (a) | centre 0.36 m reached | ROP | whole volume | mech share | Ø 25 | Ø 50 | depth-mean Ø | min Ø (feet) | mean s_c | s_c drift | mean T_rec |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| R1_scored | 150–566 s | yes | 1.50 m/h | 3.47 cm³/s | 2.91 % | 130 | 114 | 106 | 78 | 156 mm | +0.103 mm/s | 1632 K |
| R4_1600 | 150–884 s | NO (t_end) | 0.75 m/h | 1.84 cm³/s | 2.79 % | 124 | 106 | 101 | 77 | 163 mm | +0.113 mm/s | 1431 K |
| R5_1750 | 150–682 s | yes | 1.29 m/h | 2.71 cm³/s | 3.08 % | 127 | 111 | 101 | 78 | 153 mm | +0.084 mm/s | 1532 K |

| run | window (a) | ROP | note |
|---|---|---|---|
| R2_clamp | 150–384 s | 1.41 m/h | the centre reaches 0.36 m early (no far-field decay) |
| R6_d2bjet_n1 | 150–633 s | 1.46 m/h | see §7 (b) |

### T_nozzle sweep

Basis: window (a), TRANSIENT (not every run has a steady window, so no steady ROP enters the fit)

| T_nozzle | steady-window ROP | window (a) ROP | window used for the fit |
|---|---|---|---|
| 1600 K | 0.64 m/h (no window: fallback) | 0.75 m/h | 150–884 s |
| 1750 K | 0.91 m/h (no window: fallback) | 1.29 m/h | 150–682 s |
| 1900 K | 1.65 m/h (steady) | 1.50 m/h | 150–566 s |

Meier's 1.3–1.6 m/h corresponds to T_nozzle 1760 K to above 1900 K on the window (a), TRANSIENT basis (linear interpolation between the swept points). This is an inference conditional on J-M, n = 1, free-surface dilution and core 8; it is never reused as an input.

The 1.3 m/h end is **at or below 1600 K** only on the steady basis of the hand
model; on the run basis the whole band sits above 1760 K. Both readings are
conditional and neither is reused as an input.

**Addendum 2026-09-18 (D2e Goal 7), "P1 scored" row.** The settled s_c in
R1's steady window (470–694 s) is **172 mm** (mean 172.0, 174 at the window
end). The table's 156 mm is the window (a) mean. The P1 band (160–200 mm)
applies to the settled value, so **"held" stands**.

<!-- VERDICT -->

## 7. Verdict for D3

Campaign 2026-09-17 21:47 → 2026-09-18 12:07, binary `4b4f9c51…` throughout
(no rebuild). Order used: batch A (R1, R4, R5), batch B (R2, R6, R7b), R7
alone. The optional runs Oa, Ob and Oc were **not run** (budget).

| run | status | t_end | wall |
|---|---|---|---|
| R1_scored | steady-stop | 694.1 s | 3.85 h |
| R2_clamp | bottom | 727.8 s | 2.64 h |
| R4_1600 | steady-stop (watchdog; see the caveat in (d)) | 883.7 s | 4.50 h |
| R5_1750 | cap | 1299.9 s | 5.65 h |
| R6_d2bjet_n1 | steady-stop | 1097.5 s | 3.67 h |
| R7b_2mm | bottom | 561.2 s | 1.39 h |
| R7_1mm | **stopped early at 300.1 s** (user decision, 2026-09-18) | 300.1 s | 5.01 h |

No run died abnormally, and nothing is scored from a restart. The R7 stop
was the user's call to save about 6 h (400 s more at about 54 wall-s per simulated s). It was taken after the 150 s matched
window was already valid, and it is recorded in `R7_1mm.done` (`status
stopped_early`, `stop_reason`). Wall times include periods when other load on the
host slowed the runs; they are not clean cost figures.

**The headline for D3 is (e): the 2 mm results are not mesh-converged.** Every
2 mm number below, the ROP PASS included, is conditional on that.

### (a) The scored table (R1)

**A steady window exists: 470–694 s.** The steady ROP is **1.65 m/h** (ring
columns 1.69), a PASS against [1.04, 1.92], window type **steady**. On the
transient P6 basis, window (a) = 150–566 s, the ROP is 1.50 m/h.

| criterion | R1 (steady 470–694 s) | verdict |
|---|---|---|
| ROP | 1.65 m/h | PASS |
| whole volume rate | 2.59 cm³/s | PASS |
| mechanical share | 4.22 % | ok (≤ 5 %) |
| mouth Ø 25 / 50 mm | 129.9 / 114.0 vs 92 / 88 (+37.9 / +26.0) | FAIL |
| depth-mean Ø | 100.9 mm (0–246 mm) | FAIL |
| minimum Ø (feet depth) | 77.6 mm at 282 mm | unresolved at 2 mm |
| ΔT_fire | 528 K | PASS |
| removal power / R | 2.97 kW (1.05×), R = 0.779 | PASS |
| ledger | 8.1e-15 | PASS |
| mesh | matched window, +42.8 % | **FAIL** |

The rate criteria pass untuned. The shape criteria fail (mouth and
depth-mean Ø), and the mesh criterion fails.

### (b) The mouth / funnel (R1 vs R6)

| | R6 (D2b jet + n = 1) | R1 (momentum D_e, core 8, free surface) | Fig. 8.8 |
|---|---|---|---|
| Ø at 25 mm | 167.4 | 129.9 | 92 |
| Ø at 50 mm | 132.9 | 114.0 | 88 |
| depth-mean Ø | 102.6 | 100.9 | 85–93 |
| steady ROP | 1.50 m/h | 1.65 m/h | 1.04–1.92 |

**Free-surface dilution is a funnel lever, but not a sufficient one.** It
narrows the mouth by 37 mm at 25 mm and by 19 mm at 50 mm, closing about half
the excess over Fig. 8.8 (from +75 to +38 mm and from +45 to +26 mm). About
**50–58 % of the gap remains**, and since Fig. 8.8's visible width is a lower
bound, the real residual may be somewhat smaller. R6 changes three switches
at once (D_e, core, free surface), so this is their joint effect. Oa and Ob,
which would have separated them, were not run. The depth-mean Ø barely moves
(−1.7 mm): the excess sits in the top 100 mm.

### (c) The far-field law (R1 vs R2)

- **Centre.** With the 12 D clip, the centre never equilibrates. R2's s_c
  drifts at +0.45 mm/s, there is no steady window, and the centre reaches
  0.36 m by 384 s against 566 s for R1. Its mean s_c over window (a) is
  201 mm against 156 mm. The n = 1 law is what gives R1 a steady window.
- **Ring.** The law barely changes the burner. The window (a) ROP is 1.41 vs
  1.50 m/h, and the mouth Ø is identical (130 / 114 mm). The extra heat under
  the clip goes into the centre pit instead: whole volume 4.28 vs
  3.47 cm³/s on window (a), and removal power 3.65 kW (1.29× Meier) against
  2.97 kW.

The far law is a centre closure, not a ROP or shape lever.

### (d) ROP(T_nozzle), an inference

Window (a), transient basis (one basis for all three, because R4 and R5 have
no steady window): **0.75 → 1.29 → 1.50 m/h at 1600 → 1750 → 1900 K**.
Meier's 1.3–1.6 m/h maps to **1760 K to above 1900 K** by linear
interpolation. The lower end of the band (1.04 m/h) falls between 1600 and
1750 K. Below 1600 K the curve is open, because Oc (1450 K) was not run.

This is conditional on J-M, n = 1, free-surface dilution and core 8, **and on
the 2 mm mesh**. Given (e), the 1 mm curve would sit higher and the inferred T
lower. It is never reused as an input. P3 at 1600 K is refuted (§6).

**Caveat on R4.** The run-time watchdog found a window from 680 s at
t = 833.5 s and stopped the run 50 s later. The scorer, judging at t_end, finds
no window, because the added 50 s held stalled sub-fits (0.38 and 0.49 m/h).
So R4 stopped at 884 s without its centre reaching 0.36 m, and its window (a)
runs to t_end (flagged in §6). The watchdog's window was transient, and the
scorer's verdict (no window) stands. The run was not repeated. For D3: the
steady-stop watchdog should re-check the window at the stop time, or stop
later.

### (e) Mesh (matched window) and domain height

**Mesh: FAIL.** R7_1mm vs R7b_2mm, both 0.40 m, same keys, window 150–300 s:

| quantity | 2 mm | 1 mm | 1 mm vs 2 mm |
|---|---|---|---|
| burner ROP | 1.486 m/h | 2.122 m/h | **+42.8 %** (limit 5 %) |
| s_c at 300 s | 152.0 mm | 132.3 mm | −12.9 % |
| T_rec at 300 s | 1623 K | 1686 K | +3.9 % |
| centre depth at 300 s | 226 mm | 243 mm | +7.5 % |
| pinned share | 0.19 | 0.12 | −0.07 |

The gap is **systematic and growing**, not a single event. In 50 s blocks
the 2 mm ROP stays between 1.46 and 1.59 m/h, while the 1 mm ROP climbs
1.66 → 1.76 → 1.89 → 2.00 → 2.13 → 2.22 m/h. Both start together (0–50 s:
1.59 vs 1.66). At 1 mm the centre stand-off stalls near 130 mm instead of
opening towards 150 mm, so the ring sees hotter gas (T_rec +63 K) and a lower
pinned share. The 1 mm burner is already **above the Meier band (2.12 m/h)**
and still accelerating when it stops. Where it would settle, and whether it
would settle at all, is unknown. The 300 s stop keeps the matched window
valid and does not change this verdict: the gap is +13 % by 50–100 s
and +26 % by 100–150 s, both before the window opens.

What this means: D2b–D2d's 2 mm ROP agreement with Meier is **not a
converged result**. The jet march, h and clearance closures all read
cell-scale geometry (bin radii, s_c, the solid top), and at least one of them
is mesh-sensitive enough to shift the burner ROP by about 40 %. The run does
not identify which one. That diagnosis is D3's first task, before any physics
item in (h).

**Domain height (R1 at 0.70 m vs R7b at 0.40 m, common window 150–561 s):**
burner ROP 1.50 vs 1.56 m/h, **+4.08 % at 0.40 m, a finding (> 2 %)**. s_c
−0.8 %, T_rec +0.15 %, centre depth +1.1 %, mean exhaust T +0.13 %. The
shallower domain runs slightly fast; the cause (the bottom boundary
holding heat under the pit is the obvious candidate) is not tested. This is an order of magnitude below
the mesh effect.

### (f) Mechanical share and minimum Ø to the feet depth

- **Mechanical share: 4.22 %** (R1 steady), ≤ 5 %, so no finding. It is
  2.2–3.1 % over window (a) in every run and 4.56 % in R6.
- **Minimum Ø: 77.6 mm at 282 mm depth** (R1, over 0–296 mm, which includes
  the feet zone), in the 76–80 mm "unresolved at 2 mm" band. Over the drilled
  depth alone it is 81.1 mm. R2 and R7b give 77.6 and 79 mm, R6 76.1 mm at
  432 mm. **The clearance does not demonstrably hold the hole at ≥ 80 mm.**
  The shortfall (≤ 4 mm in Ø) is within one 2 mm cell per side, and the
  max-radius Ø at the same depth is 79.8 mm, so the 2 mm runs can neither
  confirm nor refute it. The 1 mm run stopped before its feet reached a
  comparable depth, so it does not resolve this.

### (g) Energy

R1 steady window: removal power **2.97 kW** (1.05× Meier's 2.83 kW), face
power (×4 quadrants) **3.81 kW**, removal efficiency **R = 0.779**, and ΔT_fire
528 K inside 500–560 K. The mean T_rec over window (a) is **1632 K**, rising
to about 1690 K by the centre's 0.36 m depth. The ledger closes to 8e-15.
Across runs R sits at 0.77–0.78. The removal power scales with the ROP (R6
2.54 kW, R2 3.65 kW), and R2's excess is centre-pit volume, not burner
advance. Energy is not the residual, at 2 mm.

### (h) The five open physics items, ranked by the residuals (not tested)

1. **Side-wall heating.** This is the largest shape residual: the mouth
   excess is +38 / +26 mm after dilution, and the depth-mean Ø is 101 against
   85–93 mm. Dilution removed about half of it, and the rest lies in the top
   100 mm, where the wall sees exhaust gas for the longest time. The wall-jet
   h along the side wall is the most likely over-prediction.
2. **Recovery temperature.** T_rec is 1630–1690 K and climbs with depth. It
   sets the gas temperature along the upper wall (see 1), and it moves with
   the mesh (+63 K at 1 mm). A lower recovery temperature would narrow the
   mouth and slow the burner together.
3. **Supersonic Martin.** The centre closure governs s_c, and s_c is where
   the mesh gap concentrates (−13 %). The far law (c) shows that the centre
   stagnation h controls whether any steady state exists. A supersonic
   correction to the Martin stagnation h acts directly on that.
4. **Pad shielding.** The minimum Ø at the feet (77.6 mm) is 2.4 mm short of
   80 and unresolved at 2 mm. This is a small residual, but the only one
   located at the pad.
5. **cp(T).** Energy already closes: removal power 1.05×, R stable at 0.78,
   ΔT_fire in band. cp(T) would shift the gas temperatures by a few per cent,
   below every other residual.

**Precondition for all five:** the mesh effect in (e) (+43 % ROP) is larger
than any residual above. D3 should resolve it first, finding which closure
reads cell-scale geometry and converging it. Until then, no physics item can
be ranked on 2 mm residuals with confidence, including this ranking.

**Addendum 2026-09-18 (D2e Goal 7; from the D2d review and the 09-18b
diagnosis). The tables above are unchanged.**
- **R1's steady window is marginal:** its s_c drift is 0.0199 mm/s against
  the 0.02 mm/s limit.
- **The 2 mm ROP PASS is a property of the 2 mm rim cycle** (09-18b §2.5):
  - the burner rate is the pad rim's spall plus whole-cell clip,
    14 + 11 ≈ 25 mm/min at 2 mm (reproduced by `studies/d2e_mesh/rim.py`);
  - it is not a converged property of the physics.
- **The T_nozzle inference in (d) is not to be quoted** until the mesh
  problem is resolved (D2e).
- **The mouth result (b) is unaffected by the seed.** The funnel forms in the
  first ~120 s, while the nozzle is above the surface, with the same flank
  removal at both meshes: the 40–60 mm band removes 1.3–2.0 cm³/s at either
  dz over 0–150 s (`rim.py`).
- **For the sweep, only the run basis (window (a)) is to be quoted.**

## Timestep ladder on R7b_2mm (2026-09-18, `dt_cap.py`; reference = R7b_2mm at DT_CAP 4 ms)
Same case and binary-identical physics, only `timestep` changed (the 5 K halving bypassed). Scored
with score.py (fallback window, no steady window on the 0.40 m domain for any of them) plus a fixed
150–550 s nozzle-descent fit. Pre-stated acceptance: ROP within 2 % of 4 ms, min Ø / mouth Ø within
one cell (2 mm), ledger at round-off.

| dt | K/step (rule) | ROP scorer | ROP fit 150–550 | t_bottom | min Ø feet depth | mouth Ø 25/50 | <h_applied> | <T_fire> | n_events | ledger | wall (new binary) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 4 ms | 1.4 | 1.60 | 1.562 | 561.2 s | 79.3 | 129.6 / 113.8 | 1.010 mm | 810.6 K | 142631 | 2.8e-14 | 83 min (old binary) |
| 8 ms | 2.8 | 1.60 | 1.561 | 563.4 s | 78.8 | 129.7 / 113.8 | 1.014 mm | 810.8 K | 142828 | 9.2e-15 | 22.0 min |
| 16 ms | 5.7 | 1.58 | 1.546 | 567.9 s | 79.3 | 129.9 / 114.1 | 1.021 mm | 811.0 K | 141161 | 9.5e-15 | 12.5 min |
| 32 ms | 11.3 | 1.67 | 1.591 | 574.1 s | 76.9 | 129.9 / 114.0 | 1.036 mm | 811.2 K | 140810 | 4.8e-15 | 5.5 min |
| 64 ms | 22.7 | 1.66 | 1.583 | 578.9 s | 77.2 | 129.9 / 114.1 | 1.064 mm | 812.1 K | 137383 | 4.5e-15 | 2.8 min |
| 128 ms | 45.3 | 1.55 | 1.522 | 588.2 s | 78.6 | 130.3 / 114.1 | 1.113 mm | 814.2 K | 131919 | 5.3e-15 | 1.5 min |

Verdict: 8 ms and 16 ms PASS (ROP −0.1 % / −1.0 %, t_bottom +0.4 % / +1.2 %, Ø within 0.5 mm).
32 ms FAILS (scorer ROP +4.4 %, min Ø −2.4 mm, a wall row drops below the Fig. 8.8 width); 64 and
128 ms fail likewise. The failure mode is not the firing temperature (+0.6 K at 32 ms, +3.6 K at
128 ms; the pinned closure caps the surface) but the flake: <h_applied> grows 1.01 → 1.11 mm and the
event count falls 7.5 %, because the depth scan K_I(z) ≥ K_Ic(T(z)) runs on a profile that has
heated past the crossing. This is exactly what the event-location proposal (interpolate the column
profile to the Sp = 1 crossing) removes; the ladder is the evidence for it.
Rule adopted for D2d: fixed step with overshoot ≤ ~6 K/step by the harness formula, i.e. 16 ms at
2 mm and 8 ms at 1 mm (DT_CAP → per-mesh; 1 mm cost 5.0 h → ~1.8 h per 300 s with the 09-18 binary).
