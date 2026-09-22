# D2e pre-registration: controlled 2 mm / 1 mm pairs

Written 2026-09-18 15:05:21 by `predict.py`, before any pair run (the harness identity and
`max_grid_size` checks are not pair runs). **Frozen: never edit.** Its sha256
and mtime are in RESULTS §0.

The switches are **diagnostics, not candidate scored configurations**. No
pair's ROP is compared with Meier's band.

## Definitions

- **Block ROP:** the least-squares slope of `nozzle_z` over the thermo rows
  with a ≤ t ≤ b and t > 0 (`run.py blocks()`), in m/h, for the blocks 0–50,
  50–100, …, 250–300 s.
- **Gap:** the 1 mm / 2 mm block-ROP ratio − 1.
- **Matched window:** 150–300 s, `score.mesh` (for P-a, 150–200 s, its
  stop).
- **Rim, pad and annulus quantities:** from `rim.py`. The rim is the
  38–40 mm band. spall = regimes 1–3, clip = regime 4. The annulus mean
  covers every column in 28–40 mm, and the burner rate is the block ROP in
  mm/min.
- **"Burner follows the annulus mean"** means that over 150–300 s the burner
  rate is closer to the annulus-mean recession than to the rim-band
  recession. It is a binary test, and no tolerance is added.

## Reference numbers from the D2d mesh pair (same definitions)

| block | 0–50 | 50–100 | 100–150 | 150–200 | 200–250 | 250–300 |
|---|---|---|---|---|---|---|
| R7b_2mm ROP [m/h] | 1.592 | 1.557 | 1.497 | 1.457 | 1.494 | 1.513 |
| R7_1mm ROP [m/h] | 1.656 | 1.758 | 1.888 | 2.000 | 2.132 | 2.221 |
| gap | +4.0 % | +12.9 % | +26.1 % | +37.3 % | +42.7 % | +46.8 % |

- **1 mm at 150 s (R7):** s_c = 120.67 mm, T_rec = 1533.85 K.
- **Matched window 150–300 s (D2d):** 2.122 vs 1.486 m/h, +42.8 %.
- **Rim, by block, at 2 mm:**
  - spall 14.4, 15.1, 14.6, 14.1, 14.4, 14.2 mm/min;
  - clip 10.3, 11.3, 10.6, 10.6, 10.6, 10.3 mm/min;
  - annulus mean 30.2, 31.1, 29.3, 28.1, 28.4, 28.3 mm/min;
  - burner 26.5, 25.9, 25.0, 24.3, 24.9, 25.2 mm/min.
- **Rim, by block, at 1 mm:**
  - spall 12.4, 13.8, 14.8, 15.7, 16.7, 17.4 mm/min;
  - clip 13.7, 15.5, 16.7, 17.6, 19.1, 19.3 mm/min;
  - annulus mean 31.4, 32.2, 32.5, 34.0, 35.8, 36.9 mm/min;
  - burner 27.6, 29.3, 31.5, 33.3, 35.5, 37.0 mm/min.
- **Pad height:** set by the 38–40 mm band on every pad in every block, at
  both meshes.
- **`rim.py` validation:** it reproduces `foot_z` to
  0.000 mm at every check time.
- **The packet quotes 1.25 / 1.31 m/h for 0–50 s.** That value comes from a
  definition not reproduced here. The gaps from 50 s on agree with it to
  0.01 m/h.

## Expectations

| pair | switch (both meshes) | isolates | expectation |
|---|---|---|---|
| B0 | none (the scored configuration) | harness (2 mm) and 1 mm step gate | **2 mm:** byte-identical to R7b_2mm_dt16 to 300 s. **1 mm gate:** 0–50, 50–100, 100–150 s ROP blocks within **2 %** of R7 (1.656, 1.758, 1.888 m/h), and s_c and T_rec at 150 s within **1 %** (120.67 mm, 1533.85 K). **Gap:** 1 mm / 2 mm over 150–300 s ≥ +30 % (D2d +42.8 %) |
| P-a | `foot_body_clearance = 0` | the seed (clip cycle) | early gap (0–50, 50–100 s blocks) ≤ 3 % (D2d +4.0 %, +12.9 %); burner slower than B0 at 2 mm; rim = spall only (no regime-4 rows). **Prior:** D2b had no clearance, and its trimmed mesh check gave ring −4.8 % / burner −6.4 % at 1 mm, small and of the opposite sign. It also differed in the far law (clamp), decay diameter (nozzle, core 5), free surface (off), bin update and domain (trimmed r ≤ 60 mm), so it is consistent with clearance as the seed but not a clean test |
| P-b | `jet_T_ent_mode = rec_fixed`, `jet_T_rec_fixed = 1550` | the amplifier (recirculation loop) | gap in the 250–300 s block ≤ its 50–100 s value + 3 points (no growth); D2d grew +13 → +47 % |
| P-c | `foot_rule = mean`, `foot_body_clearance = 0` (the clearance/pads abort forces this) | slow-tail selection (read **against P-a**; lowest value) | 150–300 s gap ≤ 5 %; burner follows the annulus mean (definition above) |

## Decision rule

- **The mesh problem is understood if P-a removes the early gap and P-b
  stops the growth.**
- **Named fallback.** If P-a does **not** shrink the early gap, the next
  suspects are:
  1. the pin idle timeout: `PinRule` unpins a column idle for
     `pinned_idle_cycles · t_cell`, with
     `t_cell = pin_rhoCp·dz·(pin_T − T_amb)/q_pin`, which scales with cell
     height;
  2. the half-cell removal increment on rim columns (`h_col = dz/2 − a_f`).

  They are recorded in RESULTS §6 as the start of the next plan, and are
  **not tested in this step.**
- **If P-b does not stop the growth,** the amplifier is not recirculation
  alone. RESULTS says what the logs point to (flank exclusion by s ≤ 0 is the
  other link).
- **If the 1 mm step gate fails,** the 8 ms / 1 mm step is not acceptable.
  Any 1 mm run already started at 8 ms is stopped and recorded. P-a, P-b and
  P-c then run at 1 mm at 4 ms, and the 1 mm baseline is **R7 itself** (4 ms,
  bit-exact binary), so B0-1 mm is not rerun.
- **What the result means for D2f** (stated, not acted on): if P-a/P-c
  converge, the burner-rate closure (feet quantile + whole-cell clip) is the
  mesh-sensitive element. The tube-wall packet must then replace it with a
  rim rule shown converged by this pair method.
