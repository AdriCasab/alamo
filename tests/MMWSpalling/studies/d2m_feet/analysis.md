# D2m analysis — is the cut diameter set by the feet or by the flame?

Binary `62dea451…`; predictions hashed before launch (RESULTS.md §0). All cases read at **matched geometry**: the time each case's feet reach the common depth **81 mm**. Never at matched time.

## 1. The ladder

| case | ring [mm] | stance Ø | **min Ø** | 2×r_outer | min Ø − stance | at depth | ROP [m/h] | Δ ROP | absorbed [W] | V̇ [cm³/s] | depth-mean Ø |
|---|---|---|---|---|---|---|---|---|---|---|---|
| R32 | 28–32 | 64 | **63.5** | 64 | -0.5 | 80 mm | **1.925** | +29.1 % | 1258 | 3.59 | 91.8 |
| R40 | 28–40 | 80 | **79.3** | 80 | -0.7 | 80 mm | **1.491** | +0.0 % | 1639 | 4.58 | 110.5 |
| R44 | 28–44 | 88 | **87.0** | 88 | -1.0 | 80 mm | **1.301** | -12.8 % | 1841 | 5.11 | 120.0 |
| R48 | 28–48 | 96 | **93.9** | 96 | -2.1 | 80 mm | **1.179** | -21.0 % | 2021 | 5.62 | 128.1 |
| S3648 | 36–48 | 96 | **94.8** | 96 | -1.2 | 80 mm | **1.154** | -22.6 % | 2052 | 5.70 | 129.4 |

## 2. The two deciding metrics

**min Ø** (weak evidence outward — see PREDICTIONS §2a: the 0.9 quantile makes Ø ≥ 2·r_outer at the feet depth close to a structural identity) **and r_half**, the flame's own erosion edge, which the support rule does not constrain.

| case | r_outer | **min Ø** | H-feet (2×r_outer ±5) | H-thermal (78±3) | **r_half [mm]** | r_half − r_outer | H-feet on r_half (±2 mm) |
|---|---|---|---|---|---|---|---|
| R32 | 32 | **63.5** | yes | no | **37.4** | +5.4 | no |
| R40 | 40 | **79.3** | yes | yes | **46.6** | +6.6 | no |
| R44 | 44 | **87.0** | yes | no | **50.6** | +6.6 | no |
| R48 | 48 | **93.9** | yes | no | **53.5** | +5.5 | no |
| S3648 | 48 | **94.8** | yes | no | **54.3** | +6.3 | no |

**The pre-registered test, as written:** H-feet on r_half requires |r_half − r_outer| ≤ 2 mm at every case. Offsets are +5.4, +6.6, +6.6, +5.5 mm, so the literal test **FAILS**.

**Slope — computed after seeing the result, reported as an addition, not a substitution.** The literal test conflates offset with slope, and here the offset is nearly constant (+5.4 to +6.6 mm) while r_half spans **16.0 mm** for **16 mm** of r_outer: **d r_half / d r_outer = 1.02** (intercept +5.2 mm). _At face value this reads as the edge moving one-for-one with the stance. **It is an artefact — see the absolute table below, which retracts it.**

#### Absolute v(r) at fixed radii — this retracts the slope reading above

`r_half` is a level crossing of `0.5·v_c`, and **`v_c` itself falls from 3.32 to 1.97 m/h across the ladder**. So r_half can march outward while the profile stands still. The absolute rates settle it:

| case | r_outer | v_c | 25 mm | 35 mm | 45 mm | 50 mm | 60 mm |
|---|---|---|---|---|---|---|---|
| R32 | 32 | 3.317 | 2.564 | 1.801 | 1.303 | 1.155 | 0.453 |
| R40 | 40 | 2.443 | 2.249 | 1.692 | 1.278 | 1.137 | 0.882 |
| R44 | 44 | 2.153 | 1.984 | 1.649 | 1.254 | 1.127 | 0.880 |
| R48 | 48 | 1.968 | 1.777 | 1.589 | 1.221 | 1.106 | 0.862 |
| S3648 | 48 | 1.939 | 1.757 | 1.575 | 1.219 | 1.096 | 0.860 |

Ladder spread (max−min as % of mean): **25 mm: 39 %**, **35 mm: 14 %**, **45 mm: 7 %**, **50 mm: 5 %**, **60 mm: 54 %**.

**At 45–50 mm the absolute erosion rate is invariant to 5–7 % across a 16 mm change of stance, while the centre moves 39 %.** The flame's far-field erosion is therefore *not* set by the stance; only the centre is. `r_half` tracked `foot_r_outer` purely through its own normalisation. **H-thermal holds on the erosion edge; H-feet holds on min Ø, which PREDICTIONS §2a flagged in advance as close to a structural identity of the support rule.**

### Cross-check: the fixed [150, 250] s window

| case | r_half matched-geometry | r_half fixed window | Δ |
|---|---|---|---|
| R32 | 37.4 | 42.7 | -5.2 |
| R40 | 46.6 | 50.7 | -4.1 |
| R44 | 50.6 | 53.1 | -2.5 |
| R48 | 53.5 | 54.2 | -0.7 |
| S3648 | 54.3 | 54.2 | +0.0 |

**The two windows disagree by more than one cell somewhere — both are reported and no single reading is forced** (PREDICTIONS §5).

## 3. `S3648` vs `R48` — outer edge, or the whole ring?

- min Ø **94.8** vs **93.9** mm (Δ +0.9, predicted within 3)
- ROP **1.154** vs **1.179** m/h (Δ -2.1 %, predicted within 5 %)
- r_half **54.3** vs **53.5** mm

**The outer edge alone sets it.**

## 4. Ø(z) at matched geometry [mm], 20 mm intervals

| case | 20 | 40 | 60 | 80 |
|---|---|---|---|---|
| R32 | 112 | 98 | 77 | 64 |
| R40 | 134 | 118 | 94 | 79 |
| R44 | 143 | 128 | 103 | 87 |
| R48 | 153 | 136 | 111 | 94 |
| S3648 | 154 | 138 | 112 | 95 |

(Meier 85–93 mm shown in the figure for orientation only — **no Meier score comes from this packet**.)

## 5. Recession-rate profile (ACTIVE_STEP 5a) — D2j-0 estimators

Matched-geometry window, 2 mm annuli, from the removal log.

| case | v_c [m/h] | r(0.9 v_c) | **r_half** | r(0.1 v_c) | **w_edge** | r_outer |
|---|---|---|---|---|---|---|
| R32 | 3.317 | 20.8 | **37.4** | 61.6 | **40.8** | 32 |
| R40 | 2.443 | 27.0 | **46.6** | nan | **nan** | 40 |
| R44 | 2.153 | 28.7 | **50.6** | nan | **nan** | 44 |
| R48 | 1.968 | 29.8 | **53.5** | nan | **nan** | 48 |
| S3648 | 1.939 | 26.4 | **54.3** | nan | **nan** | 48 |

Meier's implied floor half-width is 42.5–46.5 mm. **Does the model's floor recede flat out to ~45 mm, or is it already falling at 40?** v(r)/v_c by 2 mm annulus:

| case | 1 | 3 | 5 | 7 | 9 | 11 | 13 | 15 | 17 | 19 | 21 | 23 | 25 | 27 | 29 | 31 | 33 | 35 | 37 | 39 | 41 | 43 | 45 | 47 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| R32 | 1.01 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 0.98 | 0.98 | 0.97 | 0.94 | 0.90 | 0.82 | 0.77 | 0.71 | 0.66 | 0.61 | 0.58 | 0.54 | 0.51 | 0.47 | 0.45 | 0.42 | 0.39 | 0.37 |
| R40 | 1.01 | 0.99 | 1.00 | 1.00 | 1.01 | 1.00 | 1.00 | 0.99 | 0.97 | 0.95 | 0.95 | 0.95 | 0.92 | 0.90 | 0.83 | 0.77 | 0.73 | 0.69 | 0.66 | 0.62 | 0.58 | 0.55 | 0.52 | 0.49 |
| R44 | 0.99 | 1.01 | 1.00 | 1.00 | 1.00 | 0.98 | 0.99 | 0.98 | 0.96 | 0.95 | 0.94 | 0.94 | 0.92 | 0.93 | 0.89 | 0.85 | 0.79 | 0.77 | 0.72 | 0.68 | 0.64 | 0.60 | 0.58 | 0.55 |
| R48 | 0.98 | 1.00 | 1.00 | 0.99 | 0.98 | 0.98 | 0.97 | 0.97 | 0.95 | 0.93 | 0.93 | 0.91 | 0.90 | 0.91 | 0.91 | 0.89 | 0.84 | 0.81 | 0.76 | 0.72 | 0.68 | 0.65 | 0.62 | 0.59 |
| S3648 | 1.00 | 1.02 | 0.99 | 0.96 | 0.99 | 0.96 | 0.97 | 0.97 | 0.94 | 0.92 | 0.93 | 0.90 | 0.91 | 0.90 | 0.90 | 0.88 | 0.84 | 0.81 | 0.77 | 0.73 | 0.69 | 0.66 | 0.63 | 0.60 |

## 6. Reported, not deciding

| case | J/mm³ to rock | × floor 1.147 | jet_s_c | patch_min_standoff | jet_fs_cols | ledger_err | cap max |
|---|---|---|---|---|---|---|---|
| R32 | 1.40 | 1.22 | 125.5 | 67.2 | 0.0 | 2.8e-15 | 0.12 |
| R40 | 1.43 | 1.25 | 136.4 | 89.2 | 17.9 | 3.4e-15 | 0.14 |
| R44 | 1.44 | 1.26 | 140.9 | 96.6 | 87.2 | 4.8e-15 | 0.16 |
| R48 | 1.44 | 1.25 | 143.9 | 101.3 | 202.7 | 3.9e-15 | 0.17 |
| S3648 | 1.44 | 1.25 | 144.7 | 102.3 | 247.0 | 4.8e-15 | 0.21 |

The 4–6 J/mm³ band is an assumed 30 % delivery, not a measurement, and is not used as a target; the floor 1.147 J/mm³ (ρc_pΔT_fire) is the only independent number here.
