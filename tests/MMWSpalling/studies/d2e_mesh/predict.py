#!/usr/bin/env python3
"""Write and hash d2e_mesh/PREDICTIONS.md (D2e Goal 4), before any pair run.

The expectations and the decision rule are the packet's (ACTIVE_STEP D2e §4),
verbatim in substance. The only numbers added are read here from the D2d mesh
pair's logs (studies/d2c_steady/output/{R7b_2mm, R7_1mm}) with the D2e block
definition and rim.py; no new modelling. Refuses to overwrite an existing file.
"""
import datetime
import hashlib
import importlib.util
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rim  # noqa: E402

_spec = importlib.util.spec_from_file_location("d2e_run", os.path.join(HERE, "run.py"))
E = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(E)

OUTF = os.path.join(HERE, "PREDICTIONS.md")
D2 = os.path.join(E.D2C, "output", "R7b_2mm")
D1 = os.path.join(E.D2C, "output", "R7_1mm")
EDGES = [0.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0]


def block_rops(pf):
    b, th = E.blocks(pf, EDGES)
    return b, th


def main():
    if os.path.exists(OUTF):
        sys.exit(f"{OUTF} exists; PREDICTIONS.md is written once")
    b2, th2 = block_rops(D2)
    b1, th1 = block_rops(D1)
    gap = [x1 / x2 - 1.0 for x1, x2 in zip(b1, b2)]
    sc150 = float(np.interp(150.0, th1[:, 0], th1[:, 2]))
    tr150 = float(np.interp(150.0, th1[:, 0], th1[:, 3]))
    r2 = rim.analyse(D2, EDGES)
    r1 = rim.analyse(D1, EDGES)
    ri = int(np.nonzero(np.isclose(rim.BANDS[:-1], 0.038))[0][0])
    blk = lambda r, key: ", ".join(f"{bl[key]:.1f}" for bl in r["blocks"])
    rimtxt = lambda r, key: ", ".join(f"{bl[key][ri]:.1f}" for bl in r["blocks"])
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    L = f"""# D2e pre-registration: controlled 2 mm / 1 mm pairs

Written {now} by `predict.py`, before any pair run (the harness identity and
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
| R7b_2mm ROP [m/h] | {' | '.join(f'{x:.3f}' for x in b2)} |
| R7_1mm ROP [m/h] | {' | '.join(f'{x:.3f}' for x in b1)} |
| gap | {' | '.join(f'{x * 100:+.1f} %' for x in gap)} |

- **1 mm at 150 s (R7):** s_c = {sc150 * 1e3:.2f} mm, T_rec = {tr150:.2f} K.
- **Matched window 150–300 s (D2d):** 2.122 vs 1.486 m/h, +42.8 %.
- **Rim, by block, at 2 mm:**
  - spall {rimtxt(r2, 'sp')} mm/min;
  - clip {rimtxt(r2, 'cl')} mm/min;
  - annulus mean {blk(r2, 'ann_mean')} mm/min;
  - burner {blk(r2, 'burner')} mm/min.
- **Rim, by block, at 1 mm:**
  - spall {rimtxt(r1, 'sp')} mm/min;
  - clip {rimtxt(r1, 'cl')} mm/min;
  - annulus mean {blk(r1, 'ann_mean')} mm/min;
  - burner {blk(r1, 'burner')} mm/min.
- **Pad height:** set by the 38–40 mm band on every pad in every block, at
  both meshes.
- **`rim.py` validation:** it reproduces `foot_z` to
  {max(r2['worst'], r1['worst']) * 1e3:.3f} mm at every check time.
- **The packet quotes 1.25 / 1.31 m/h for 0–50 s.** That value comes from a
  definition not reproduced here. The gaps from 50 s on agree with it to
  0.01 m/h.

## Expectations

| pair | switch (both meshes) | isolates | expectation |
|---|---|---|---|
| B0 | none (the scored configuration) | harness (2 mm) and 1 mm step gate | **2 mm:** byte-identical to R7b_2mm_dt16 to 300 s. **1 mm gate:** 0–50, 50–100, 100–150 s ROP blocks within **2 %** of R7 ({b1[0]:.3f}, {b1[1]:.3f}, {b1[2]:.3f} m/h), and s_c and T_rec at 150 s within **1 %** ({sc150 * 1e3:.2f} mm, {tr150:.2f} K). **Gap:** 1 mm / 2 mm over 150–300 s ≥ +30 % (D2d +42.8 %) |
| P-a | `foot_body_clearance = 0` | the seed (clip cycle) | early gap (0–50, 50–100 s blocks) ≤ 3 % (D2d {gap[0] * 100:+.1f} %, {gap[1] * 100:+.1f} %); burner slower than B0 at 2 mm; rim = spall only (no regime-4 rows). **Prior:** D2b had no clearance, and its trimmed mesh check gave ring −4.8 % / burner −6.4 % at 1 mm, small and of the opposite sign. It also differed in the far law (clamp), decay diameter (nozzle, core 5), free surface (off), bin update and domain (trimmed r ≤ 60 mm), so it is consistent with clearance as the seed but not a clean test |
| P-b | `jet_T_ent_mode = rec_fixed`, `jet_T_rec_fixed = 1550` | the amplifier (recirculation loop) | gap in the 250–300 s block ≤ its 50–100 s value + 3 points (no growth); D2d grew {gap[1] * 100:+.0f} → {gap[5] * 100:+.0f} % |
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
"""
    with open(OUTF, "w") as f:
        f.write(L)
    sha = hashlib.sha256(open(OUTF, "rb").read()).hexdigest()
    mt = datetime.datetime.fromtimestamp(os.path.getmtime(OUTF)).strftime("%Y-%m-%d %H:%M:%S")
    print(f"PREDICTIONS.md sha256 {sha}\nmtime {mt}")


if __name__ == "__main__":
    main()
