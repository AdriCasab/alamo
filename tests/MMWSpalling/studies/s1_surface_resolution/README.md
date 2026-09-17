# S1: 1-D surface-resolution study

**Purpose.** This asks how the mesh size dz at the drilled surface changes
four things: ROP, the temperature at which the spallation criterion fires,
the absorbed flux, and melt. It covers fixed-flux (beam) and Robin (flame)
sources on 1-D columns: 1 × 1 × N cubic cells, 250 mm deep, Neumann-0 on
every face. It is a study, not a regression: there is no `test`, and
`scripts/runtests.py` does not pick it up. Results and the discussion are in
[RESULTS.md](RESULTS.md) (accepted 2026-09-15).

## How to run

```bash
VENV=/Users/tzetze20/Desktop/code/.venv/bin/python
$VENV run_sweep.py --smoke      # one dz per case (<= 10 s) + 1-vs-4-rank identity check
$VENV run_sweep.py --jobs 8     # full matrix, 1 rank per run  -> output/
$VENV analyze.py                # regenerates RESULTS.md above <!-- DISCUSSION -->
```

- `input` is the base; `run_sweep.py` sets dz, source, stop time and
  plot_file by CLI `key=value` overrides.
- `STOP_OVERRIDE` in `run_sweep.py` holds measured time-to-200-mm stop times
  for the fast cases.

## Caveats to carry

- **Degenerate Weibull draw.** The base input uses `weibull.enabled = 1`,
  `a0_gb = a0_ig = 20e-6`, `m = 1e12` as a scalar flaw. This worked around a
  segfault with `weibull.enabled = 0` + removal. Since A1 (2026-09-15) the
  scalar path works (`weibull.enabled = 0` + `spallation.sp.a0`), and
  `unit/scalar_flaw` shows the two give identical T/H. `run_sweep.py` is left
  as is, so the results stay reproducible.
- **`h_col_events.csv` is edge-triggered.** It logs one row per column's
  *first* firing at each new top cell; repeat firings at the same top cell
  are not logged. Never sum `h_col` for recession: it undercounts about 2×.
  For a complete record use `spall.removal_events_csv = 1` (A1), which
  writes one row per (step, column) removal with `h_applied` and `n_voided`.
- **`h_col` is a removal increment, not a flake size.** In self-heated runs
  the criterion (`top_cell`, `local_thermoelastic`) is a cell-average
  temperature threshold. The depth scan returns `h_col ≈ dz/2 − a_f`, and
  two firings on consecutive steps void one cell. Flake and PSD claims are
  retired for the drilling regime. Energy and ROP are unaffected: every
  voided cell leaves at ρ·Cp·ΔT_fire.
- **Overshoot time-step rule (binding for removal runs).** One step deposits
  `ΔT_step = q·dt/(ρ·Cp·dz)` into the top cell, and that is the
  firing-temperature overshoot (ρ·Cp = 2.1725e6 J/m³K). For ≤ 5 K:
  - dt ≤ 11 ms at dz = 2 mm, q = 2 MW/m²;
  - dt ≤ 1.4 ms at 2 mm, 15 MW/m²;
  - dt ≤ 0.7 ms at 0.125 mm, 2 MW/m².

  The S1 M-cases (15 MW/m², dt = 1 ms) overshoot 50–110 K per step at fine
  dz. Their firing-T spread is that overshoot, not physics.
