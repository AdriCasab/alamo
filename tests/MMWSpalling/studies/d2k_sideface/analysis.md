# D2k arbiter analysis (side-face flux on the D2j-0b C1 cascade)

## 1. Disk-mean rate (r < 30 mm) by block, 2 mm / 1 mm

| block | 2 mm | 1 mm | gap | D2j-0b C1 (key off) | D2j-0b C0 |
|---|---|---|---|---|---|
| 50–100 s | 1.5971 | 1.6509 | **+3.4 %** | -13.5 % | -10.4 % |
| 100–150 s | 1.6207 | 1.6455 | **+1.5 %** | -28.7 % | -15.3 % |
| 150–200 s | 1.6761 | 1.6233 | **-3.2 %** | -45.5 % | -20.1 % |

**(a) PASS** — every deciding block within 5 %.

## 2. Band timing (D2j-0b format): t_stop / t_cross / Δ_own / Δ_outer

A band **stalls** only if it STOPS FIRING before its own face crosses the plane.
D2j-0b marks a band still firing at the end with t_stop = t_end (`250 (firing)`),
so t_stop >= t_end - 5 s is NOT a stall however large t_cross - t_stop looks.

| run | band [mm] | n | t_stop | t_cross | Δ_own | below plane at stop [mm] | Δ_outer |
|---|---|---|---|---|---|---|---|
| C1_2mm | 18–20 | 15 | 250 (firing) | ∞ | — | 22.5 | ∞ |
| C1_2mm | 20–22 | 17 | 250 (firing) | ∞ | — | 22.5 | ∞ |
| C1_2mm | 22–24 | 16 | 250 (firing) | ∞ | — | 22.5 | ∞ |
| C1_2mm | 24–26 | 23 | 250 (firing) | ∞ | — | 22.5 | ∞ |
| C1_2mm | 26–28 | 19 | 250 (firing) | ∞ | — | 22.5 | ∞ |
| C1_2mm | 28–30 | 25 | 250 (firing) | ∞ | — | 22.5 | ∞ |
| C1_1mm | 18–20 | 61 | 250 (firing) | ∞ | — | 24.4 | ∞ |
| C1_1mm | 20–22 | 66 | 250 (firing) | ∞ | — | 24.4 | ∞ |
| C1_1mm | 22–24 | 69 | 250 (firing) | ∞ | — | 24.4 | ∞ |
| C1_1mm | 24–26 | 81 | 250 (firing) | ∞ | — | 24.4 | ∞ |
| C1_1mm | 26–28 | 86 | 250 (firing) | ∞ | — | 24.4 | ∞ |
| C1_1mm | 28–30 | 89 | 250 (firing) | ∞ | — | 24.4 | ∞ |

**(b) PASS** — no band inside 28 mm stalls with Δ_own > 10 s at either mesh (the better outcome).

## 3. r_front(t) at 1 mm [mm] (inner edge of the stalled run reaching 30 mm)

| run | 100 s | 150 s | 200 s | 250 s |
|---|---|---|---|---|
| C1_2mm | — | — | — | — |
| C1_1mm | — | — | — | — |

**(c) PASS** — no inward-walking front at 1 mm.

## 4. Side power as a fraction of top-face power (reported, not deciding)

| run | block | side_P_face [W] | patch_P_robin [W] | side / top | side_faces | max abs ledger_err |
|---|---|---|---|---|---|---|
| C1_2mm | 50–100 s | 26.7 | 387.8 | 6.9 % | 13 | 1.5e-15 |
| C1_2mm | 100–150 s | 26.9 | 388.0 | 6.9 % | 13 | 3.2e-15 |
| C1_2mm | 150–200 s | 29.0 | 387.8 | 7.5 % | 14 | 2.9e-15 |
| C1_1mm | 50–100 s | 29.5 | 387.1 | 7.6 % | 57 | 4.4e-15 |
| C1_1mm | 100–150 s | 29.2 | 387.1 | 7.6 % | 57 | 4.0e-15 |
| C1_1mm | 150–200 s | 29.6 | 387.1 | 7.6 % | 57 | 5.5e-15 |
| C1f07_2mm | 50–100 s | 27.3 | 388.5 | 7.0 % | 19 | 2.6e-15 |
| C1f07_2mm | 100–150 s | 27.4 | 388.6 | 7.1 % | 19 | 3.8e-15 |
| C1f07_2mm | 150–200 s | 28.4 | 388.4 | 7.3 % | 20 | 3.8e-15 |

## 5. f sensitivity at 2 mm (pre-registered, never fitted)

| f | run | 50–100 s | 100–150 s | 150–200 s |
|---|---|---|---|---|
| 1.0 | C1_2mm | 1.5971 | 1.6207 | 1.6761 |
| 0.7 | C1f07_2mm | 1.5839 | 1.6205 | 1.6634 |

f = 0.7 vs f = 1.0: -0.8 %, -0.0 %, -0.8 %

## 6. Verdict: PASS

(a) blocks +3.4, +1.5, -3.2 % -> ok; (b) ok; (c) ok.

Per CRITERION.md: run the conditional Meier pair under a separately hashed D2i criterion.
