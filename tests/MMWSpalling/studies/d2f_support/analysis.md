## Harness identity (E6)

| file | rows (t <= 199.968 s, Pa_2mm's end) | Q09_2mm vs d2e Pa_2mm |
|---|---|---|
| thermo.dat | 2084 | byte-identical |
| removal_events.csv | 68447 | byte-identical |
| jet_profile.csv | 1303 | byte-identical |

Identity: PASS

### Q10: foot_pad_quantile = 1.0, foot_body_clearance = 0

Matched window 150–300 s: 2 mm 1.448 m/h, 1 mm 1.403 m/h, **gap -3.1 %**.

| block | 2 mm ROP | 1 mm ROP | gap |
|---|---|---|---|
| 0–50 s | 1.556 | 1.577 | +1.3 % |
| 50–100 s | 1.543 | 1.528 | -1.0 % |
| 100–150 s | 1.474 | 1.462 | -0.8 % |
| 150–200 s | 1.437 | 1.398 | -2.7 % |
| 200–250 s | 1.426 | 1.409 | -1.2 % |
| 250–300 s | 1.502 | 1.357 | -9.7 % |

| quantity | t | 2 mm | 1 mm |
|---|---|---|---|
| reach [mm] | 100 s | 168.3 | 169.0 |
| reach [mm] | 200 s | 63.1 | 62.8 |
| T_rec [K] | 100 s | 1555.8 | 1555.5 |
| T_rec [K] | 200 s | 1550.3 | 1552.1 |
| s_c [mm] | 100 s | 112.2 | 112.0 |
| s_c [mm] | 200 s | 138.0 | 137.7 |
| P_face [W] | 100 s | 2750.0 | 2742.8 |
| P_face [W] | 200 s | 1653.5 | 1644.9 |

| mesh | rim spall + clip by block [mm/min] | annulus mean | burner | 40–60 mm [cm³/s] | regime-4 rows | max abs ledger_err |
|---|---|---|---|---|---|---|
| 2mm | 23.7+0.0, 26.9+0.0, 25.5+0.0, 24.5+0.0, 25.1+0.0, 24.4+0.0 | 29.5, 31.4, 29.3, 28.0, 28.2, 28.1 | 25.9, 25.7, 24.6, 24.0, 23.8, 25.0 | 1.24, 1.97, 1.94, 1.90, 1.73, 1.28 | 0 | 4.9e-15 |
| 1mm | 25.4+0.0, 26.3+0.0, 25.1+0.0, 24.0+0.0, 24.3+0.0, 23.8+0.0 | 31.4, 30.7, 28.9, 27.5, 27.3, 26.6 | 26.3, 25.5, 24.4, 23.3, 23.5, 22.6 | 1.42, 1.96, 1.93, 1.86, 1.47, 0.93 | 0 | 7.3e-15 |

`score.mesh`:

Matched window 150–300 s (t_end 300 s at 1 mm, 300 s at 2 mm; steady windows are NOT required, D2d decision 2).

| quantity | 2 mm | 1 mm | 1 mm vs 2 mm |
|---|---|---|---|
| burner ROP over the window | 1.448 m/h | 1.403 m/h | -3.1 % |
| s_c at the window end | 153.3 mm | 153.7 mm | +0.2 % |
| T_rec at the window end | 1618 K | 1630 K | +0.7 % |
| centre depth at the window end | 226 mm | 224 mm | -0.9 % |
| pinned share | 0.20 | 0.18 | -0.02 |

**Mesh verdict: PASS** (burner ROP within 5 %).

### Q09: foot_pad_quantile = 0.9, foot_body_clearance = 0

Matched window 150–300 s: 2 mm 1.491 m/h, 1 mm 1.459 m/h, **gap -2.1 %**.

| block | 2 mm ROP | 1 mm ROP | gap |
|---|---|---|---|
| 0–50 s | 1.608 | 1.640 | +2.0 % |
| 50–100 s | 1.574 | 1.576 | +0.1 % |
| 100–150 s | 1.526 | 1.498 | -1.9 % |
| 150–200 s | 1.455 | 1.456 | +0.1 % |
| 200–250 s | 1.521 | 1.475 | -3.0 % |
| 250–300 s | 1.462 | 1.400 | -4.3 % |

| quantity | t | 2 mm | 1 mm |
|---|---|---|---|
| reach [mm] | 100 s | 168.3 | 169.0 |
| reach [mm] | 200 s | 61.7 | 60.6 |
| T_rec [K] | 100 s | 1556.0 | 1555.1 |
| T_rec [K] | 200 s | 1556.7 | 1562.6 |
| s_c [mm] | 100 s | 112.0 | 112.0 |
| s_c [mm] | 200 s | 136.0 | 136.2 |
| P_face [W] | 100 s | 2750.7 | 2745.1 |
| P_face [W] | 200 s | 1623.0 | 1594.6 |

| mesh | rim spall + clip by block [mm/min] | annulus mean | burner | 40–60 mm [cm³/s] | regime-4 rows | max abs ledger_err |
|---|---|---|---|---|---|---|
| 2mm | 23.7+0.0, 26.9+0.0, 25.5+0.0, 24.5+0.0, 25.1+0.0, 25.0+0.0 | 29.5, 31.5, 29.2, 28.2, 28.5, 28.4 | 26.8, 26.2, 25.4, 24.3, 25.3, 24.4 | 1.24, 1.97, 1.94, 1.89, 1.62, 1.19 | 0 | 5.6e-15 |
| 1mm | 25.4+0.0, 26.3+0.0, 25.0+0.0, 24.2+0.0, 24.5+0.0, 23.3+0.0 | 31.4, 30.7, 28.9, 27.7, 27.6, 26.8 | 27.3, 26.3, 25.0, 24.3, 24.6, 23.3 | 1.41, 1.96, 1.93, 1.82, 1.29, 0.75 | 0 | 7.9e-15 |

`score.mesh`:

Matched window 150–300 s (t_end 300 s at 1 mm, 300 s at 2 mm; steady windows are NOT required, D2d decision 2).

| quantity | 2 mm | 1 mm | 1 mm vs 2 mm |
|---|---|---|---|
| burner ROP over the window | 1.491 m/h | 1.459 m/h | -2.1 % |
| s_c at the window end | 153.3 mm | 153.0 mm | -0.2 % |
| T_rec at the window end | 1624 K | 1641 K | +1.1 % |
| centre depth at the window end | 228 mm | 227 mm | -0.4 % |
| pinned share | 0.19 | 0.16 | -0.03 |

**Mesh verdict: PASS** (burner ROP within 5 %).

## q-sensitivity (2 mm, 150 s – end)

| q | ROP [m/h] | hand | s_c / T_rec at end | hand s_c / T_rec | mean foot_carry_cols | max stall [s] | burner / annulus mean [mm/min] | pad-setting bands | foot_z recon. |
|---|---|---|---|---|---|---|---|---|---|
| 1.0 | 1.448 | 1.678 | 153.3 / 1618 | 152.8 / 1604 | 8.6 | 8.2 | 24.1 / 28.1 | 38–40, 38–40, 38–40 | 0.000 mm |
| 0.9 | 1.491 | 1.741 | 153.3 / 1624 | 150.7 / 1614 | 22.3 | 7.8 | 24.9 / 28.4 | 38–40, 38–40, 38–40 | 0.000 mm |
| 0.8 | 1.535 | 1.809 | 150.7 / 1633 | 148.5 / 1624 | 36.6 | 7.4 | 25.6 / 28.6 | 38–40, 38–40, 38–40 | 0.000 mm |

Rim-zone spall by 2 mm band over 150 s – end [mm/min] (no clip in any D2f run):

| q | 28–30 | 30–32 | 32–34 | 34–36 | 36–38 | 38–40 |
|---|---|---|---|---|---|---|
| 1.0 | 32.5 | 30.8 | 29.0 | 27.3 | 26.2 | 24.6 |
| 0.9 | 32.9 | 31.0 | 29.3 | 27.5 | 26.4 | 24.8 |
| 0.8 | 33.4 | 31.2 | 29.5 | 27.7 | 26.6 | 25.1 |

Hole at the run end (direction only; azimuth-mean Ø, score.py):

| q | t | Ø 25 mm | Ø 50 mm | min Ø to feet depth (at depth) | feet depth |
|---|---|---|---|---|---|
| 1.0 | 300 s | 130.3 | 114.1 | 78.7 (122 mm) | 123 mm |
| 0.9 | 300 s | 128.6 | 113.1 | 77.7 (124 mm) | 125 mm |
| 0.8 | 300 s | 126.4 | 110.5 | 75.3 (128 mm) | 129 mm |

## PREDICTIONS rows

| id | expectation | measurement | verdict |
|---|---|---|---|
| E6 | Q09_2mm byte-identical to d2e Pa_2mm to 200 s | see identity table | held |
| E1 (gate i) | Q10 1 mm / 2 mm within 5 % over 150–300 s and every block after 50 s | window -3.1 %; blocks +1.3, -1.0, -0.8, -2.7, -1.2, -9.7 % | REFUTED |
| E2 (gate ii) | 2 mm ROP for q 0.8 / 0.9 / 1.0 within ±5 % of the mean | 1.535 / 1.491 / 1.448 m/h, mean 1.491, spread ±2.9 %; -0.044 m/h per +0.1 q | held |
| E3 | ROP falls as q rises | 1.535 > 1.491 > 1.448? | held |
| E4 | q 1.0 ROP within ±27 % of hand 1.678 | 1.448 m/h (-13.7 %) | held |
| E5 | q 1.0 s_c at end within ±20 mm of hand 152.8 | 153.3 mm | held |
| E4 | q 0.9 ROP within ±27 % of hand 1.741 | 1.491 m/h (-14.4 %) | held |
| E5 | q 0.9 s_c at end within ±20 mm of hand 150.7 | 153.3 mm | held |
| E4 | q 0.8 ROP within ±27 % of hand 1.809 | 1.535 m/h (-15.1 %) | held |
| E5 | q 0.8 s_c at end within ±20 mm of hand 148.5 | 150.7 mm | held |

**Decision outcome: (b): score q = 0.9 (jolt allowance)**
