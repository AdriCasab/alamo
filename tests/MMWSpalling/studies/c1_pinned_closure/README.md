# C1: pinned-surface Robin closure study

**Purpose.** C1 checks `surface_patch.robin_form = pinned`, the flame Robin
closure. S1 showed that neither Robin form converges at dz = 2 mm: `cell`
over-absorbs by 42 % and `face` under-absorbs by 82 %. Every run fires at a
fixed cell-average temperature T_fire, so the resolved surface rides at
T_fire.

After a column's first spall firing the closure pins the surface at that
firing's top-cell T. The flame flux and the losses are both evaluated there:

`q = h(T_gas − T_pin) − εσ(T_pin⁴ − T_a⁴) − h_c(T_pin − T_a)`

The face form is kept as a min-rule fallback, T_s = min(T_face, T_pin), and it
also covers idle columns.

This is a study, not a regression; there is no `test`. The regression is
`tests/MMWSpalling/unit/robin_pinned`. Results are in
[RESULTS.md](RESULTS.md).

## How to run

```bash
VENV=/Users/tzetze20/Desktop/code/.venv/bin/python
$VENV run.py --part 1 --jobs 6     # 1-D A/B vs S1 (1 rank per run)  -> output/p1_*
$VENV run.py --part 2 --jobs 2     # 2-D Meier-harness check (4 ranks per run) -> output/p2_*
$VENV analyze.py                   # regenerates RESULTS.md above <!-- DISCUSSION -->
```

Nothing is copied from other directories:
- **Part 1** uses the S1 base input (`../s1_surface_resolution/input`) and S1's
  `overrides()` with the `R-face` / `R-face-1400` sources, then
  `robin_form=pinned`.
- **Part 2** overrides `validation/meier/sp_meier_pilot/input_2d_dev` on the
  command line. The input itself is not edited.

`analyze.py` loads S1's helpers via importlib, because both files are named
`analyze.py`.

## Settings to carry

- **`weibull.V0 = V_cell`** in every Part 2 run (user decision 2026-09-15):
  `8e-9` at 2 mm and `1e-9` at 1 mm. With the default 1e-9, T_fire drifts
  ~+10 K per halving (S1b). This applies to C1 overrides only; defaults and
  existing inputs are unchanged. Part 1 keeps S1's degenerate a = 20 µm draw,
  so it compares like-for-like with S1.
- **Overshoot time-step rule** (binding for removal runs):
  ΔT_step = q·dt/(ρ·Cp·dz) ≤ 5 K, with q the closed-form flux at the expected
  T_fire. `run.py` halves dt from 1 ms until the rule holds and prints dt and
  K/step per run.
- **`surface_patch.pinned_idle_cycles`** (default 2) is reported as a
  sensitivity, never tuned. It is a pinned-only key, so face or cell runs must
  not carry it; the strict parser aborts on it.
- **Stop times.**
  - Part 1: min(60 s, 150 mm of closed-form recession).
  - Part 2: min(40 s, 80 mm at the closed form with T_fire = 822 K).
- **Not checkpointed.** The pin state and the energy ledger restart from face
  / re-baseline after a restart.
