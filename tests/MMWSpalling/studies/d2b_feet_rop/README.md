# D2b: burner on feet and the scored Meier ROP / hole-diameter test

**Purpose.** Meier's burner rests on three feet (Ø80 mm, nozzle 50 mm above
the feet, nozzle Ø7.5 mm). It descends when the rock under the feet recedes,
so the measured ROP is the recession at the burner radius. D2a2 used a
prescribed feed instead, and the flank lagged. D2b places the burner on the
simulated rock and scores the result against Meier. It uses the D2a2 wall-jet
enthalpy closure, with one change: the stagnation mix entrains the hole's own
exhaust rather than 293 K air.

This is a study, not a regression, so there is no `test`. The regressions are
`tests/MMWSpalling/unit/robin_feet` (code) and
`validation/meier/sp_meier_pilot/test_feet` (the scored Meier test, a
check-only script over this study's Stage B J-M output). Results and
pre-registration are in [RESULTS.md](RESULTS.md). **No parameter was adjusted
toward Meier's ROP, hole diameter or volume rate.**

## How to run

```bash
VENV=/Users/tzetze20/Desktop/code/.venv/bin/python
$VENV feetmodel.py                     # hand model (numbers only)
$VENV run.py --predict                 # RESULTS.md §0: hand prediction + Stage A rule (refuses to overwrite)
$VENV run.py --stage A --jobs 3        # wall-treatment decision, J-M 150 s x 3
$VENV analyze.py --stage A             # applies the rule, prints the choice
$VENV run.py --probe --wall A2         # cost probe, 10 s of J-M at 2 mm
$VENV run.py --probe-1mm --wall A2     # cost probe at 1 mm (trimmed domain)
$VENV run.py --stage B --wall A2       # J-M (scored), J-5, J-10
$VENV run.py --stage D --wall A2       # mesh check: 1 mm and 2 mm, r <= 60 mm
$VENV run.py --stage C --wall A2       # fixed 293 K; fixed + nozzle mass; T_nozzle 1436 K
$VENV analyze.py                       # RESULTS.md tables between the markers + PNGs
$VENV run.py --write-input-feet --wall A2 --stop-feet <s>   # validation input
```

Runs use 4 MPI ranks, with `--jobs 3` in parallel by default. Every run is
`validation/meier/sp_meier_pilot/input_drilling` plus CLI overrides; the full
key list is in `run.py` `keys()`. Output goes to `output/<name>/`, with a
`<name>.done` file holding the run metadata and its stop status.

## The code (in `src/Integrator/MMWSpalling.H`; default-off)

- **`surface_patch.nozzle_descent = feet`** (`FeetDescent`, once per step at
  level 0):
  - Live columns in the annulus [`foot_r_inner`, `foot_r_outer`] are split
    into `foot_npads` azimuth sectors of the observed azimuth range.
  - Each pad's height is the nearest-rank `foot_pad_quantile` of its column
    top faces (`foot_rule = pads`). Alternatives: `max` or `mean` over the
    whole annulus.
  - z_foot = mean pad height; z_target = z_foot + `foot_standoff`.
  - z_n = min(z_prev, max(z_target, z_prev − `nozzle_feed_max`·dt)). The
    burner never rises; the first step sets z_n = z_target.
  - `nozzle_z0` and `nozzle_feed` are rejected under `feet`. The
    `nozzle_collision_radius` default becomes `foot_r_inner`.
- **`surface_patch.jet_T_ent_mode = exhaust`** (requires
  `jet_entrained_mass = 1`):
  - The stagnation mix uses T_rec = the previous step's `jet_T_exhaust`
    (T_ent on the first step).
  - The march floor stays at `jet_T_ent`.
  - The budget becomes P_nozzle + P_recirc, with
    P_recirc = (m − mdot)·cp·(T_rec − T_ent).
- **Thermo columns:**
  - `nozzle_z`, `foot_z`, `foot_cols`, `foot_carry_cols` (columns at or
    above their pad height), `foot_stall_time` (simulated time
    since z_n last decreased; replaced the packet's per-step
    `foot_stalled_steps`, which counted almost every step at 2 mm);
  - `jet_T_rec`, `jet_P_recirc`.
- **Not checkpointed.** `nozzle_z` and T_rec restart from the first-step rule.

## Harness details

- **Watchdog (`run.py`).** It polls `thermo.dat` every 10 s of wall time and
  stops the run, keeping the data, when either:
  - `nozzle_z` has not descended for 60 s of simulated time (`stalled`); or
  - the centre face is 40 mm above the domain bottom (`bottom`).
- **dt** follows the overshoot rule (≤ 5 K per step; 4 ms cap), evaluated at
  the hand model's stagnation flux.
- **Stall signal.** The watchdog status, plus the longest hold computed
  from `nozzle_z` by `analyze.py`. The study runs were made before
  `foot_stall_time` existed.
- **Stage A window** is 50–150 s. Stages B, C and D use the steady window
  found by `analyze.py`: descent-rate sub-fits within ±10 % and
  `jet_s_c` drift under 0.02 mm/s. When no such window exists, the report
  uses the second half of the run (at least 150 s; 60 s for Stage D) and
  marks it "not steady".

## Hand model (`feetmodel.py`)

- The D2a2 march with a firing face (T_s = T_fire) everywhere inside the
  quantile radius r_q = 38.97 mm.
- s rises linearly from s_c at the core edge to 50 mm at r_q.
- T_rec is solved as a fixed point.
- The wall profile comes from the freeze depth
  z(r) = 50·v(r)/(v − v(r)), measured below the feet.
  - **Correction to the packet:** the packet's d(r) = 50·v/(v − v(r)) is the
    feet's depth when column r freezes. It is deeper than the column's own
    depth by the 50 mm stand-off.
- Under exhaust recirculation the model has **no centre equilibrium**: the
  centre stays ahead of the ring, and the deep-pit limit (`S_C_MAX`) gives
  the ring ROP. Numbers are in RESULTS §0.
  - **Cause (found after §0 was written):** the h expression clips the
    stand-off at 12 D. Beyond 90 mm the centre keeps its 12 D coefficient,
    1.86× the ring's, while recirculation stops T_stag from decaying
    toward T_fire.
  - With the unclipped Martin shape, the J-M exhaust case does reach a
    centre equilibrium, at s_c ≈ 0.40 m and 1.93 m/h (T_rec 1732 K),
    against 1.82 m/h clipped.
  - So the simulated pit depth is an artefact of extrapolating beyond
    Martin's range, and the wall treatment is its limiter by default. The
    ring feels the pit only through the face power, so report face power
    against `jet_s_c`.
- **Stage A's A3 is effectively A2 (no coherence cap).** Under
  the pinned closure, the surface-normal kernel's cos θ thickness is divided
  back out by the applied recession.

## Scoring (see the RESULTS §0 addendum)

- Volume rate and rock-side power are scored inside Ø 93 (r ≤ 46.5 mm).
- Hole Ø is scored over the 0.5 m block. Columns above the nozzle keep
  their depth. Active columns are projected with the run's v(r) until the
  nozzle reaches them.
