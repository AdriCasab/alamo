# D2l Gate 0 analysis — key-only scan of the side-face factor

Binary `62dea451…`; predictions hashed before launch (sha256 in RESULTS.md §0).

## 1. Gate 0a — the cascade floor (D2j-0b C1 arbiter, side flux on)

### Single-mesh signature at 2 mm (band stalling — what the scan detects)

| f | disk-mean rate by block [mm/s] | bands stalling inside 28 mm | r_front walks in |
|---|---|---|---|
| 1.0 | 1.5971 / 1.6207 / 1.6761 | none | no |
| 0.5 | 1.5832 / 1.6164 / 1.6610 | none | no |
| 0.2 | 1.5985 / 1.5563 / 1.6178 | none | no |
| 0.1 | 1.4638 / 1.4450 / 1.3748 | none | no |

### Mesh gap where both legs were run (the deciding condition)

| f | block | 2 mm | 1 mm | gap | verdict |
|---|---|---|---|---|---|
| 1.0 | 50–100 s | 1.5971 | 1.6509 | **+3.4 %** | ok |
| 1.0 | 100–150 s | 1.6207 | 1.6455 | **+1.5 %** | ok |
| 1.0 | 150–200 s | 1.6761 | 1.6233 | **-3.2 %** | ok |
| 0.5 | 50–100 s | 1.5832 | 1.6415 | **+3.7 %** | ok |
| 0.5 | 100–150 s | 1.6164 | 1.6382 | **+1.4 %** | ok |
| 0.5 | 150–200 s | 1.6610 | 1.6257 | **-2.1 %** | ok |
| 0.2 | 50–100 s | 1.5985 | 1.4957 | **-6.4 %** | **> 5 %** |
| 0.2 | 100–150 s | 1.5563 | 1.4314 | **-8.0 %** | **> 5 %** |
| 0.2 | 150–200 s | 1.6178 | 1.3920 | **-14.0 %** | **> 5 %** |
| 0.1 | 50–100 s | 1.4638 | 1.3084 | **-10.6 %** | **> 5 %** |
| 0.1 | 100–150 s | 1.4450 | 1.1531 | **-20.2 %** | **> 5 %** |
| 0.1 | 150–200 s | 1.3748 | 1.0007 | **-27.2 %** | **> 5 %** |

**f_min = 0.5** — the smallest f on the ladder meeting all three conditions.

**f is a diagnostic dial here, not a calibration.** It scales q_side at exactly the place an h scale would, which is why it maps the window; its physical meaning is still AREA, and no value in this table may be kept as a calibrated setting.

## 2. Gate 0b — rate, shape and cost vs f (Meier LI0 keys, 2 mm, 250 s)

Every run is compared at **matched geometry**: the time its feet reach the common depth **131 mm**, the deepest depth every run on the ladder reached. Never at matched time.

| f | ROP 150–250 s [m/h] | absorbed (quarter) [W] | side / face | V̇ [cm³/s] | depth-mean Ø [mm] | min Ø [mm] | TSE(HHV) [J/mm³] | TSE(to-rock) [J/mm³] |
|---|---|---|---|---|---|---|---|---|
| 1 | **4.066** | 3204 | 53 % | 9.73 | **101.9** | **78.0** | 3.9 | 1.32 |
| 0.5 | **3.759** | 2524 | 47 % | 7.58 | **97.8** | **77.1** | 5.0 | 1.33 |
| 0.2 | **3.087** | 1980 | 37 % | 5.75 | **100.4** | **77.7** | 6.6 | 1.38 |
| 0.1 | **2.137** | 1837 | 24 % | 5.25 | **104.1** | **77.5** | 7.2 | 1.40 |
| 0 | **1.491** | 1639 | 0 % | 4.58 | **107.6** | **78.0** | 8.3 | 1.43 |

Ø(z) profile at matched geometry [mm], z = depth below the original top:

| f | 10 mm | 30 mm | 50 mm | 70 mm | 90 mm | 110 mm | 130 mm |
|---|---|---|---|---|---|---|---|
| 1 | 121 | 107 | 104 | 104 | 103 | 94 | 78 |
| 0.5 | 126 | 107 | 99 | 96 | 94 | 88 | 77 |
| 0.2 | 136 | 117 | 104 | 97 | 92 | 85 | 78 |
| 0.1 | 141 | 122 | 109 | 102 | 96 | 86 | 78 |
| 0 | 145 | 126 | 113 | 106 | 100 | 88 | 78 |

Bands: ROP 1.3–1.6 m/h · Ø 85.0–93.0 mm · TSE(HHV) 15.4 · TSE(to-rock) 4.0–6.0. Thermodynamic floor ρc_pΔT_fire = 1.147 J/mm³.

**f_rate = 0.000** (run directly in the band).

Ø at matched feet depth is not confound-free either: a faster run **reaches that depth earlier**, so its wall has had longer to widen by the end of the run. Dwell = t_end − t(feet reach the matched depth) makes the size of that bias visible; Ø(t_end) is the same wall measured at the end of the run instead.

| run | t at matched depth [s] | dwell [s] | Ø at match [mm] | Ø at t_end [mm] |
|---|---|---|---|---|
| f = 1 | 131 | 119 | 78.0 | 104.1 |
| f = 0.5 | 152 | 98 | 77.1 | 95.3 |
| f = 0.2 | 205 | 45 | 77.7 | 87.6 |
| f = 0.1 | 250 | 0 | 77.5 | 77.5 |
| f = 0 | 314 | 258 | 78.0 | 94.0 |

| run | cap max | max abs ledger_err | jet_s_c [mm] | jet_fs_cols | feet depth end [mm] |
|---|---|---|---|---|---|
| f = 1 | 0.39 | 4.1e-15 | 68.5 | 0.0 | 267 |
| f = 0.5 | 0.27 | 5.4e-15 | 83.4 | 0.0 | 233 |
| f = 0.2 | 0.19 | 8.1e-15 | 108.5 | 0.0 | 171 |
| f = 0.1 | 0.16 | 3.7e-15 | 126.0 | 0.0 | 131 |
| f = 0 | 0.14 | 3.4e-15 | 136.4 | 17.9 | 236 |

## 3. Gate 0c — the competing explanation (distribution)

Discriminating signature: **rate DOWN and Ø UP together.** No cut to side flux can produce that, so either probe showing it implicates the distribution.

Ø at the matched feet depth **97 mm**.

| probe | f | change | ROP [m/h] | Δ rate | **Ø at z_m [mm]** | **Δ Ø (deciding)** | Δ depth-mean Ø | Δ min Ø | identical to base? | signature? |
|---|---|---|---|---|---|---|---|---|---|---|
| `P_flat` | 1 | far-field exponent n 1.0 -> 0.5 | 4.066 | +0.0 % | **83.5** | **+0.0** | +0.0 | +0.0 | no | no |
| `P_fsoff` | 1 | jet_free_surface off | 4.066 | -0.0 % | **83.5** | **+0.0** | +3.0 | +0.0 | no | no |
| `P_flat0` | 0 | far-field exponent n 1.0 -> 0.5 | 1.460 | -2.1 % | **82.3** | **+0.0** | +0.1 | +0.0 | no | no |
| `P_fsoff0` | 0 | jet_free_surface off | 1.457 | -2.3 % | **82.6** | **+0.3** | +17.1 | +0.3 | no | no |

**Deciding column is Δ Ø at z_m**, the metric fixed in PREDICTIONS.md §3 before any run. The depth-mean is shown beside it because a single depth is fragile — the Ø(z) profiles of different runs cross, and z_m can land near a crossing. **A mean-based rule would have flipped this verdict**: `P_fsoff0` moves the depth-mean by +17.1 mm and would read as the signature. It is not treated as one, for two reasons: it is not the pre-registered metric, and it moves the wrong way — the depth-mean is already too WIDE (112 mm against Meier's 85–93), so +17 mm is a worse funnel, not a better hole. The quantity that fails LOW is the minimum Ø (D2d: 77.6 against 85–93), and that moves +0.3 mm.

Base rows for reference:

| base | f | ROP [m/h] | depth-mean Ø [mm] | min Ø [mm] |
|---|---|---|---|---|
| `M_f10` | 1 | 4.066 | 102.3 | 83.5 |
| `LI0_2mm` | 0 | 1.491 | 112.0 | 82.3 |

## 4. The decision (PREDICTIONS.md §2, applied)

- `f_min` = **0.5** (cascade suppression floor)
- `f_rate` = **0.000** (rate re-enters 1.3–1.6 m/h)

**Window is EMPTY: f_rate < f_min.**

Per ACTIVE_STEP: **STOP. Write no code.** The wall closure cannot satisfy both gates — the same parameter controls cascade suppression and the power excess and they pull opposite ways. Do not run Gate 1, do not write `side_face_h_scale`.

