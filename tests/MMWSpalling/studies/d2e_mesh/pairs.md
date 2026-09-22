## B0 identity (2 mm) and 1 mm step gate

| file | rows (t <= 299.904 s, B0's end) | B0_2mm vs R7b_2mm_dt16 |
|---|---|---|
| thermo.dat | 3125 | byte-identical |
| removal_events.csv | 93077 | byte-identical |
| jet_profile.csv | 1602 | byte-identical |

Identity: PASS

| block | B0_1mm (8 ms) | R7_1mm (4 ms) | diff |
|---|---|---|---|
| 0–50 s | 1.652 m/h | 1.656 m/h | -0.22 % |
| 50–100 s | 1.759 m/h | 1.758 m/h | +0.06 % |
| 100–150 s | 1.875 m/h | 1.888 m/h | -0.67 % |
| s_c at 150 s | 0.12 | 0.12067 | -0.55 % |
| T_rec at 150 s | 1533.7 | 1533.8 | -0.01 % |

Step gate (blocks within 2 %, s_c and T_rec within 1 %): PASS

### B0: none (the scored configuration)

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

### Pa: foot_body_clearance = 0

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

## PREDICTIONS rows

| pair | expectation | measurement | verdict |
|---|---|---|---|
| B0 | 2 mm byte-identical to R7b_2mm_dt16 to 300 s | see §4 identity table | held |
| B0 | 1 mm step gate (blocks 2 %, s_c/T_rec 1 %) | see §4 gate table | held |
| B0 | gap over 150–300 s ≥ +30 % | +42.7 % | held |
| P-a | early gap (0–50, 50–100 s) ≤ 3 % | +2.0 %, +0.1 % | held |
| P-a | burner slower than B0 at 2 mm (every block to 200 s) | 1.608 vs 1.588, 1.574 vs 1.549, 1.526 vs 1.494, 1.455 vs 1.438 | REFUTED |
| P-a | rim = spall only (no regime-4 rows) | 0 / 0 rows | held |
| P-b | all | not run | not tested |
| P-c | all | not run | not tested |
