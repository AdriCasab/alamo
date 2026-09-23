# D2o-0b analysis — prescribed floor/wall heating, field bounded in radius

Binary `62dea451…`; predictions hashed before launch (RESULTS §0). Compared at **matched feet depth 39 mm**, never matched time. **The control's depth-mean is recomputed at that depth, not taken as the packet's 110.5** (PREDICTIONS §0: it is 112.3 at 97 mm and 100.7 over the whole hole).

## 1. Q1 — the gate: did the radial corner stop the jam?

| case | h_wall | ROP 150–250 s [m/h] | vs control | feet depth at 600 s | `foot_stall_time` end [s] | Q1 |
|---|---|---|---|---|---|---|
| control | — | **1.491** | — | 236 mm | 1 | — |
| W_2mm | 40 | **0.771** | 0.52× | 78 mm | 11 | **FAIL** |
| Wz_2mm | 0 | **0.340** | 0.23× | 39 mm | 47 | **FAIL** |
| Wsens_2mm | 20 | **0.522** | 0.35× | 77 mm | 14 | **FAIL** |

Q1 requires `foot_stall_time` < 10 s, feet depth ≥ 150 mm by 600 s, and ROP within 30 % of 1.491 m/h. D2o-0 (stand-off corner) gave ROP 0.108 m/h and stall 219 s.

## 2. Shape — P1 and P2

| case | shaft Ø min/mean/max (z>100mm) | depth-mean Ø | mouth Ø | V̇ [cm³/s] | P1 whole | P1 mean | P2 mean | P2 mouth | P3 |
|---|---|---|---|---|---|---|---|---|---|
| control | — unscoreable | **109.8** | **128.2** | **4.58** | — | — | FAIL | FAIL | FAIL |
| W_2mm | — unscoreable | **79.1** | **79.3** | **2.00** | — | — | FAIL | PASS | FAIL |
| Wz_2mm | — unscoreable | **78.1** | **78.4** | **1.63** | — | — | FAIL | PASS | FAIL |
| Wsens_2mm | — unscoreable | **78.2** | **78.6** | **1.88** | — | — | FAIL | PASS | FAIL |

P1 whole profile 80–100 mm and mean 85–95; P2 depth-mean 85–100 and mouth ≤ 110; P3 V̇ 2.0–3.2 cm³/s (Meier 2.47, control 4.58).

## 3. Ø(z) at matched geometry [mm], 20 mm intervals

| case | 20 |
|---|---|
| control | 109 |
| W_2mm | 79 |
| Wz_2mm | 78 |
| Wsens_2mm | 78 |

## 4. v(r) — absolute, 2 mm annuli (never normalised)

| case | v_c [m/h] | 25 mm | 35 mm | 45 mm | 50 mm | 55 mm | 60 mm |
|---|---|---|---|---|---|---|---|
| control | 3.757 | 2.623 | 1.768 | 1.163 | 0.966 | 0.694 | 0.524 |
| W_2mm | 2.161 | 1.822 | 1.378 | 0.000 | 0.000 | 0.000 | 0.000 |
| Wz_2mm | 0.921 | 0.627 | 0.276 | 0.000 | 0.000 | 0.000 | 0.000 |
| Wsens_2mm | 1.870 | 1.541 | 1.216 | 0.000 | 0.000 | 0.000 | 0.000 |

## 5. Wall — P4, P5, on a GEOMETRIC band

A wall column has r > 40 mm, its face **at least one cell below the original rock surface**, and 0 < s ≤ 45 mm. **Not `s` alone** — that is what put 3298 columns of untouched rock into D2o-0's P3.

| case | t [s] | band cols | T_wall min/mean/max [K] | 821 − T_max | wall heat [W] | P4 T | P4 P | P5 |
|---|---|---|---|---|---|---|---|---|
| W_2mm | — | 0 | — | — | — | — | — | — |
| Wz_2mm | — | 0 | — | — | — | — | — | — |
| Wsens_2mm | — | 0 | — | — | — | — | — | — |
| control | 550 | 79 | 589 / **748** / 840 | -19 K | — (impinging field) | PASS | — | **FIRES** |

Wall temperature by depth [K], 50 mm bins:

| case | 0–50 | 50–100 | 100–150 | 150–200 | 200–250 | 250–300 | 300–350 | 350–400 |
|---|---|---|---|---|---|---|---|---|
| control | — | — | — | 767 | 736 | — | — | — |

## 6. Mesh — 2 mm vs 1 mm at matched feet depth

**The 1 mm leg was stopped at t = 56.9 s of 600** and never reached the 150–250 s scoring window, so it produces no row here. It was ended deliberately: Q1 already fails at 2 mm and refinement can only make the step-face sink worse (it scales as 1/dx), the packet requires both meshes only for a **positive** claim, and shape is unscoreable regardless. Its one usable datum — the corner step at the common time t = 50 s — is in RESULTS.md §2. **This study is single-mesh and nothing in it is converged.**

## 7. Ledger and absorbed power

| case | max abs ledger_err | patch_P_robin [W] |
|---|---|---|
| control | 3.4e-15 | 1639 |
| W_2mm | 5.6e-15 | 1078 |
| Wz_2mm | 4.2e-15 | 735 |
| Wsens_2mm | 5.2e-15 | 941 |
