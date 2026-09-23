# D2o-0 analysis — prescribed floor/wall heating

Binary `62dea451…`; predictions hashed before launch (RESULTS §0). Ø compared at **matched feet depth 48 mm**, never matched time. P1/P2 are scored **below 10 cm depth**, where no rock saw the flat-start transient (PREDICTIONS §4 Q2).

## 1. THE RESULT: the prescribed split jams the burner

| case | h_wall | skirt | ROP 150–250 s [m/h] | vs control | nozzle descent in 450 s | foot_stall_time at t_end [s] | patch_P_robin [W] |
|---|---|---|---|---|---|---|---|
| control | — | — | **1.491** | — | 236.0 mm | 1 | 1639 |
| W_2mm | 40 | grazing | **0.108** | 0.07× | 55.3 mm | 23 | 1137 |
| Wz_2mm | 40 | zero | **0.053** | 0.04× | 48.0 mm | 219 | 775 |
| Wsens_2mm | 20 | grazing | **0.025** | 0.02× | 50.7 mm | 147 | 1067 |

**The burner descends ~27 mm and then stops.** `patch_min_standoff` grows from 60 to 181 mm over the run: the centre pit keeps deepening (it is floor-zone and stays impingement-heated) while the feet are held up. The hole never reaches 10 cm depth, so **P1 and P2 are unscoreable as defined** — there is no shaft.

### Why: a locking ratchet between the 0.9 pad quantile and the s = 45 mm switch

The feet rest on the **0.9 nearest-rank quantile** of the annulus face heights, so ~10 % of annulus rock sits **above** the pad plane, i.e. at s < the 50 mm stand-off. **The floor/wall switch at s = 45 mm cuts straight through that population** — the switch is only 5 mm below the plane the feet themselves define.

| run | annulus cols in the WALL zone (s ≤ 45 mm) | their surface T | fate |
|---|---|---|---|
| control | 5 of 162 | **680–814 K** | at the 821 K threshold — keep firing, tool descends |
| W_2mm | 15 of 162 | **444–482 K** | 340 K below firing — **frozen permanently, tool jams** |

`foot_carry_cols` = 18 in W_2mm against 20–26 in the control, and those carriers cannot spall at h = 40. **Any annulus column that drifts below the switch freezes and holds the burner up for the rest of the run.** This is a property of where the transition was placed, not of the grazing coefficient: all three variants (h = 40 grazing, h = 40 with a zero skirt, h = 20) stall, and their hole shapes are identical to 0.1 mm.

## 2. Shape — such as it is

| case | h_wall | skirt | shaft Ø min/mean/max (z>10cm) | **mouth gain** | **steady band gain** | shaft gain | mouth Ø | ROP [m/h] | V̇ [cm³/s] |
|---|---|---|---|---|---|---|---|---|---|
| control | — | — | — (no shaft: hole < 10 cm) | **+28.4** | **+24.7** | — | 136.7 | 1.491 | 4.58 |
| W_2mm | 40 | grazing | — (no shaft: hole < 10 cm) | **+7.9** | **+8.8** | — | 95.8 | 0.108 | 1.93 |
| Wz_2mm | 40 | zero | — (no shaft: hole < 10 cm) | **+7.9** | **+8.8** | — | 95.8 | 0.053 | 1.90 |
| Wsens_2mm | 20 | grazing | — (no shaft: hole < 10 cm) | **+7.9** | **+8.8** | — | 95.8 | 0.025 | 1.91 |

P1 target shaft Ø **85–95 mm**; P2 target gain **2.5–7.5 mm/side**; skirt = 80 mm.


Control values, which is how the definitions were checked: mouth gain **+28.4**, steady band gain **+24.7** mm/side. **These are not the control's own values** (+32.4 / +5.0, which reproduce ACTIVE_STEP's ~33): the comparison depth collapsed to 48 mm because the test runs stalled, so every case is read inside its own start-up funnel. That is why the shape table decides nothing here.

- **W_2mm**: **P1 and P2 UNSCOREABLE** — the hole never reached 10 cm depth (feet at 40 mm at end of run). Reporting them from the start-up funnel would be meaningless.
- **Wz_2mm**: **P1 and P2 UNSCOREABLE** — the hole never reached 10 cm depth (feet at 40 mm at end of run). Reporting them from the start-up funnel would be meaningless.
- **Wsens_2mm**: **P1 and P2 UNSCOREABLE** — the hole never reached 10 cm depth (feet at 40 mm at end of run). Reporting them from the start-up funnel would be meaningless.

## 3. Ø(z) at matched geometry [mm], 20 mm intervals

| case | 20 | 40 |
|---|---|---|
| control | 118 | 87 |
| W_2mm | 85 | 78 |
| Wz_2mm | 85 | 78 |
| Wsens_2mm | 85 | 78 |

## 4. v(r) — absolute, 2 mm annuli (D2m estimators, never normalised)

| case | v_c [m/h] | 25 mm | 35 mm | 45 mm | 50 mm | 55 mm | 60 mm |
|---|---|---|---|---|---|---|---|
| control | 3.595 | 2.696 | 1.848 | 1.318 | 1.118 | 0.835 | 0.658 |
| W_2mm | 2.813 | 2.357 | 1.700 | 0.000 | 0.000 | 0.000 | 0.000 |
| Wz_2mm | 1.970 | 1.666 | 1.395 | 0.000 | 0.000 | 0.000 | 0.000 |
| Wsens_2mm | 2.669 | 2.301 | 1.680 | 0.000 | 0.000 | 0.000 | 0.000 |

**Q6 asks whether v(r) becomes flatter (cylinder) or steeper (narrower cone).** Ratio v(45 mm)/v_c:

| case | v(45)/v_c |
|---|---|
| control | 0.367 |
| W_2mm | 0.000 |
| Wz_2mm | 0.000 |
| Wsens_2mm | 0.000 |

## 5. Wall temperature and wall heat (P3, P4, P5)

Surface temperature of the band columns (0 < s ≤ 45 mm) at end of run, from the plotfile. **Wall heat is computed from the prescribed field and the wall area** — with `jet_closure = none` there is no march and therefore no in-code wall-power diagnostic.

| case | t [s] | band cols | T_wall min/mean/max [K] | 821 − T_max | **wall heat [W]** | P3 | P4 | P5 |
|---|---|---|---|---|---|---|---|---|
| W_2mm | 450 | 239 | 444 / **565** / 612 | +209 K | **160** | FAIL | FAIL | no |
| Wz_2mm | 450 | 3298 | 320 / **333** / 441 | +380 K | **11** | FAIL | FAIL | no |
| Wsens_2mm | 450 | 400 | 417 / **496** / 526 | +295 K | **143** | FAIL | FAIL | no |
| control | 550 | 84 | 589 / **749** / 840 | -19 K | — (impinging field) | PASS | — | **FIRES** |

Wall temperature by depth [K], 50 mm bins (end of run):

| case | 0–50 | 50–100 | 100–150 | 150–200 | 200–250 | 250–300 | 300–350 | 350–400 |
|---|---|---|---|---|---|---|---|---|
| W_2mm | 565 | — | — | — | — | — | — | — |
| Wz_2mm | 333 | — | — | — | — | — | — | — |
| Wsens_2mm | 496 | — | — | — | — | — | — | — |
| control | — | — | — | 767 | 739 | — | — | — |

## 6. Mesh — 2 mm vs 1 mm by 100 s block (a single mesh is not quotable)

_1 mm leg not complete_

## 7. Ledger and absorbed power

| case | max abs ledger_err | patch_P_robin [W] | J/mm³ to rock (floor 1.147) |
|---|---|---|---|
| control | 3.4e-15 | 1639 | 1.43 |
| W_2mm | 1.0e-14 | 1137 | 2.35 |
| Wz_2mm | 2.9e-15 | 775 | 1.64 |
| Wsens_2mm | 5.0e-15 | 1067 | 2.24 |
