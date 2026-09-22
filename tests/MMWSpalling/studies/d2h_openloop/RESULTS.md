# D2h results — three key-only pairs on the D2f support rule (reviewer diagnostic, 2026-09-20)

Criterion: `CRITERION.md`, sha256 `0e23e117…da0c3`, written 15:02:54 before launch (first run
directory 15:03). Binary `6d047b50…`, no build. Harness `run.py` imports d2f → d2e → d2c by
path; nothing frozen was edited. Runs under `caffeinate -dis`; paused 16:25–16:5x by
SIGSTOP/SIGCONT (pause, not restart; `PAUSED.md`). 1 mm plotfiles every 50 s (8 × 1.58 GB per run).

## Runs
| run | switch | status | t_end | wall |
|---|---|---|---|---|
| I0_2mm / I0_1mm | pinned_idle_cycles 2 → 1e9 | ok / ok | 349.9 s | 9 min / 225 min (incl. pause) |
| RF_2mm / RF_1mm | jet_T_ent_mode exhaust → rec_fixed 1600 K | ok / ok | 349.9 s | 9 min / 226 min |
| OL_2mm / OL_1mm | RF + prescribed nozzle 1.508 m/h, feet + free-surface keys dropped | **aborted** | 119.3 s | — |

OL aborted on the nozzle collision guard: under `prescribed` its radius defaults to the patch
radius (0.2 m), under `feet` to `foot_r_inner` (0.028 m). The guard fired when the nozzle plane
reached the untouched far surface. Fixed in `run.py` (`nozzle_collision_radius = 0.028`); not rerun.

## 1. Burner ROP and pad-band (28–40 mm) spall recession, 2 mm / 1 mm (gap)
| pair | ROP 150–250 s | ROP 250–350 s | pad v 150–250 | pad v 250–350 |
|---|---|---|---|---|
| L09 (D2g ref) | 1.492 / 1.466 (−1.8 %) | 1.511 / 1.322 (**−12.5 %**) | 28.3 / 27.7 (−2.4 %) | 28.3 / 25.9 (**−8.7 %**) |
| **I0 idle off** | 1.491 / 1.461 (−2.1 %) | 1.505 / 1.444 (**−4.0 %**) | 28.3 / 27.7 (−2.4 %) | 28.2 / 27.1 (**−3.8 %**) |
| RF rec fixed | 1.528 / 1.488 (−2.6 %) | 1.470 / 1.249 (−15.0 %) | 29.0 / 28.1 (−3.2 %) | 27.4 / 24.3 (−11.2 %) |

## 2. Silent fraction (since last firing > 2 t_cell) by 4 mm band at 350 s; pad firings/col/100 s
| run | 24–28 | 28–32 | 32–36 | 36–40 | 40–44 | 44–48 | 48–52 | fir/col |
|---|---|---|---|---|---|---|---|---|
| L09 2 mm | 0 | 0 | 0 | 0 | 0 | 9 | 75 | 46.9 |
| L09 1 mm | 2 | 4 | 9 | 9 | 42 | 81 | 100 | 85.6 |
| I0 2 mm | 0 | 0 | 0 | 0 | 0 | 0 | 52 | 46.7 |
| **I0 1 mm** | 0 | 0 | 0 | **0** | **1** | 22 | 99 | 90.3 |
| RF 2 mm | 0 | 0 | 0 | 0 | 0 | 26 | 94 | 45.4 |
| RF 1 mm | 5 | 7 | 9 | 24 | 40 | 79 | 100 | 80.5 |

## 3. State at 350 s (2 mm / 1 mm)
| pair | T_rec K | s_c mm | reach mm | pinned share | idle cols | P_face W |
|---|---|---|---|---|---|---|
| L09 | 1646 / 1666 | 157 / 159 | 49.8 / 48.5 | 0.126 / 0.080 | 25 / 594 | 1202 / 1108 |
| I0 | 1640 / 1649 | 157 / 157 | 50.0 / 48.6 | 0.136 / 0.127 | 0 / 0 | 1231 / 1188 |
| RF | 1600 / 1600 | 159 / 161 | 49.6 / 48.3 | 0.118 / 0.077 | 41 / 658 | 1139 / 1059 |

## 4. Fields (plotfiles; T_top by band, silent = above the nozzle plane)
- Q09/L09 1 mm at 300 s (idle on): band 50–52 mm silent (T_top 618 K, 92 % < 700 K, Sp 0.39).
- I0 1 mm at 300 s (idle off): 50–52 mm half firing (737 K, 41 % < 700 K, Sp 0.71); 48–50 firing at 804 K.
- I0 1 mm at 350 s: bands 34–48 mm all firing at 786–809 K; 48–50 silent (614 K) — its face
  (89.7 mm) is above the plane (feet 145 → plane 95 mm). I0 2 mm at 350 s: 48–50 firing (face
  103 mm, below the plane). With the idle clock off, silence = above the plane, both meshes.
- Lateral T at the rim-floor level is the same in all cases (40–42 mm: 585–592 K; 44–46: 428–438 K):
  the cold rock beside the rim is unchanged; what changed is whether the edge columns keep firing.

## 5. Verdict (per CRITERION.md)
- **I0 removes the front**: 1 mm pad recession 250–350 s within 5 % of 2 mm (−3.8 %) AND
  36–40 mm silent fraction 0 % at 1 mm. Burner ROP gap −4.0 % (L09: −12.5 %).
- **RF does not** (−11.2 %, 24 %): the recirculation ratchet is not the amplifier.
- OL not tested (abort).
- **The reviewer's prior ("I0 does not remove the front") was wrong**; the other agent's
  reading of the idle clock was right.

## 6. Mechanism (from the code paths and the fields)
The pin state is read only by `SurfaceCellFlux` (via `pin_T_eff`). Under the pin the cell
absorbs h(T_gas − T_pin) − losses **regardless of its own temperature**, so it always reaches the
firing condition. After the idle clock drops the pin, the cell gets the face form
G(T_f − T_c), which is self-limiting: the cell approaches T_f and, when lateral loss into the
cold rock beside it (∝ 1/dx) holds T_f near or below the firing threshold, it settles at a
sub-critical plateau (the 665 K columns seen in the fields) and never fires again. Its cold
top then lengthens the next column's firing interval past 2 t_cell → idle → same fate: the
front walks inward. Both ingredients are cell-scaled: t_cell ∝ dz (the clock) and the lateral
loss ∝ 1/dx. Removing the clock leaves the state rule T_s = min(T_f, T_pin), which is enough.

## 7. What this does and does not establish
- Established: the D2g late divergence to 350 s is the idle clock acting at a lateral edge.
- Not established: convergence beyond 350 s (the I0 gap went −2.1 → −4.0 %; L09 went −1.8 →
  −12.5 → −40 %). A long I0 pair to the 0.40 m bottom (~550 s; 1 mm ≈ 3 h alone) is the
  decisive test of the support rule. The wall-side silence (44–52 mm) is still the `s ≤ 0` rule
  and is not addressed here.
- Changing the default of `pinned_idle_cycles` rebases every pinned-closure result since C1;
  do it as a new default only after re-checking the reference set, or keep the key and set it
  in the Meier inputs.
