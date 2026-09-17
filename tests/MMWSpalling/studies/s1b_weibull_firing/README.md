# S1b: firing temperature under the production Weibull block

**Purpose.** Packet C's pinned-surface Robin closure uses
`q = h·(T_gas − T_fire)`, so T_fire has to be known under the Meier flaw
settings. This study measures it with those settings (`a0_gb = 20 µm`,
`a0_ig = 4 µm`, `m = 20`, `seed = 12345`) on the S1 1-D column under the F1
beam. It then checks whether T_fire drifts with dz through the Weibull
volume factor (V_cell/V0)^(1/m), and compares V0 choices. It is a study,
not a regression: there is no `test`. Results and the recommendation are in
[RESULTS.md](RESULTS.md).

## How to run

```bash
VENV=/Users/tzetze20/Desktop/code/.venv/bin/python
$VENV run.py --jobs 6                 # W-def, W-cell, W-seed, S-4p4 -> output/
$VENV run.py --meier-keys             # + W-def 2 mm with one Meier 2D key each
$VENV analyze.py                      # regenerates RESULTS.md above <!-- DISCUSSION -->
```

Nothing is copied from S1:
- `run.py` uses `../s1_surface_resolution/input` and S1's `overrides()` for
  the dz / F1 keys;
- `analyze.py` loads S1's plotfile helpers.

**Common settings:**
- 1 rank, dt = 1 ms, 60 s;
- `spall.removal_events_csv = 1`, the complete removal log (A1), which
  supplies T_top and a_f at every firing;
- overshoot q·dt/(ρ·Cp·dz) = 0.47 / 0.95 / 1.9 / 3.8 K per step at
  dz = 2 / 1 / 0.5 / 0.25 mm, all within the ≤ 5 K rule.

**Runtime:** each 0.25 mm run takes ~1–2 min on an idle machine (the
recorded timings were under heavy load).

**Cases:**

| Case | Setting |
|---|---|
| W-def | default V0 = 1e-9 m³ |
| W-cell | V0 = dz³ (vol_factor ≡ 1); at 1 mm it must equal W-def byte-for-byte, and `analyze.py` checks this |
| W-seed | seed 777 |
| S-4p4 | scalar a0 = 4.436 µm (no scatter; needs A1's scalar-flaw path) |
