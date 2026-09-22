# D2c: steady-state and no-funnel closures, with pre-registration

**Purpose.** D2b scored the Meier test on a transient.
- The 12 D clip in `h_expr` kept the centre pit running to the domain bottom,
  so no run had a steady window.
- The mouth grew a funnel, so the whole-excavation volume was 2.3× Meier's.

D2c adds four closures, all default-off:
- a far-field h law;
- a momentum decay diameter;
- free-surface dilution;
- burner-body clearance.

It also adds two ride-alongs, fixes the scoring, and pre-registers hand
predictions.

**D2c runs no scored simulation.** The campaign is D2d. It runs from
`run.py --stage D2d` and is scored by `score.py`.

## How to run (the order matters)

```bash
VENV=/Users/tzetze20/Desktop/code/.venv/bin/python
$VENV nozzle.py                  # Laval state, D_e, jet_De_ref (checks the packet table)
$VENV predict.py                 # PREDICTIONS.md + RESULTS.md §0 hash (refuses to overwrite)
$VENV score.py --dry-run ../d2b_feet_rop/output/B_JM_A2   # RESULTS.md §Dry run
$VENV run.py --probe 2mm         # 30 s, scored configuration, 0.40 m domain
$VENV run.py --probe 1mm         # 10 s, 1 mm, full domain
$VENV run.py --project           # wall-clock projection (RESULTS.md §Probes)
$VENV run.py --write-input       # validation/meier/sp_meier_pilot/input_feet_d2c
# D2d only:
$VENV run.py --stage D2d [--cases R1_scored ...]
$VENV score.py output/R1_scored
$VENV score.py --mesh output/R1_scored output/R7_1mm
$VENV score.py --sweep output/R4_1600 output/R5_1750 output/R1_scored
```

## The code (`src/Integrator/MMWSpalling.H`, `MMWSpalling/Removal.H`; lev 0)

Every key below defaults to D2b behaviour. The thermo columns appear only when
their key is non-default, and are appended after all older columns.

### `surface_patch.jet_decay_diameter = nozzle | momentum`

- Under `momentum`:
  - φ = min(1, `jet_core_length`·D_e/s_c), with
    D_e = `jet_De_ref`·√(T_mix/`jet_T_ent`);
  - T_mix is the lagged T_rec under exhaust, else `jet_T_ent`.
- `jet_De_ref` is the **full** jet's 2ṁ/√(πρJ) at `jet_T_ent`, from
  `nozzle.py`. Never scale it by the sector's mdot share.
- `jet_core_length` must be given explicitly.
- Requires `jet_stagnation = decay`.
- Thermo: `jet_D_e`.
- **Ricou–Spalding** entrainment = `momentum` + `jet_core_length = 3.125` +
  `jet_entrained_mass = 1`. No separate key.
- With `jet_De_ref = jet_D`, core 5 and fixed T_ent, the result is
  bit-identical to `nozzle` (unit test (c)).

### `surface_patch.jet_free_surface = 0 | 1`

Keys: `jet_fs_aspect` (1) and `jet_fs_exponent` (1). Requires
`nozzle_descent = feet`.

- **Classification.**
  - z0 = each column's face height at the first march call. It is host
    state and **not checkpointed**.
  - A marching column is free-surface if r > `foot_r_outer` and
    z0 − z_face ≤ aspect·(r − `foot_r_outer`).
- **Per bin, after its heat is removed:**
  - m_out = m·(r_hi/r_lo)^(e·n_fs/n);
  - T_gas ← T_ent + (T_gas − T_ent)·m/m_out;
  - the jet's excess enthalpy flux is conserved.
- Empty bins do not dilute.
- A free-surface column in the axis bin (r_lo = 0) aborts.
- **T_rec** (exhaust mode, next step) is the inlet T_gas of the first bin
  holding a free-surface column (`jet_T_mouth`), else T_exhaust.
  - Modelling choice: the gas that recirculates is the gas in the hole, not
    the diluted gas leaving over the surface.
  - In a deep hole the nozzle plane is below the original surface, so no
    free-surface column marches, and the closure acts only while the feet are
    shallower than about 50 mm.
- **Invariants** are unchanged: P_face + m_exit·cp·(T_exhaust − T_ent) = P_cap.
- Thermo: `jet_fs_cols`, `jet_m_ratio` (m_exit/m_stag), `jet_T_mouth`.
- Profile CSV extra columns: `n_fs`, `m_ratio`.

### `surface_patch.jet_bin_update = explicit | exponential`

- Per bin, q is evaluated at T_in and at T_in − 1 K through the same
  `PinRule` / `SurfaceCellFlux` calls:
  - G = Σ(dq)·dA, T_eq = T_in − P_in/G, N = G/(m·cp);
  - the bin's columns see T_eff = T_eq + (T_in − T_eq)(1 − e^−N)/N;
  - P_bin is re-evaluated at T_eff, and T_out = T_in − P_bin/(m·cp).
- For affine q this gives T_out = T_eq + (T_in − T_eq)·e^−N exactly.
- G ≤ 0 falls back to explicit.
- Profile CSV extra column: `T_eff`.

### `surface_patch.jet_negative_flux = drop | count`

- Under `count`, P_bin is the signed sum, so rock hotter than the gas heats
  the gas.
- The clamp acts only when P_bin > 0.
- Invariants: P_face is signed, and P_face − `jet_P_neg` ≥ 0 replaces
  P_face ≥ 0. The rest are unchanged.
- Thermo: `jet_P_neg`. Profile CSV extra column: `P_neg`.

### `surface_patch.foot_body_clearance = 0 | 1`

Requires `nozzle_descent = feet`, `foot_rule = pads` and removal.

- **Arming.** `FeetDescent` arms every annulus column whose top cell is above
  its pad's nearest-rank cell. It uses the same pick on the integer tops,
  because z_face is monotone in k.
- **Order within a step.** Removal runs first, then the thermal advance
  (which calls `BuildFlameColumns` → `FeetDescent`). The cut therefore acts
  at the **next** step's removal pass, before PASS 1, so the thermal passes
  see the cleared surface.
- **The cut.**
  - Cells above the pad's k are flipped to removed **by index**:
    regime 4, D = 0, thermal state kept.
  - The column's φ is shifted by n·dz, so the new top keeps its sub-cell φ.
  - A tie-flipped cell gets φ = −tiny, and the new top is held at φ ≥ 0.
  - **Pin state is untouched.**
- **No ratchet:** clipping values above a nearest-rank quantile down to it
  leaves the quantile unchanged (unit test (h)).
- **Ledger:** unchanged. Void H is frozen and still counted, so the cut is
  ledger-neutral.
- **Logging.**
  - `removal_events.csv` gets one row per cleared column: regime 4,
    `k_top` = the top before the cut, `h_applied` = n·dz, `n_voided` = cells
    flipped, the other fields 0.
  - Thermo: `foot_mech_vol` (cumulative m³ in the sector) and
    `foot_mech_cols` (per step).
- **Guard order.** The removal guard in the code is a backstop:
  `surface.follow_mask = 1` already requires removal and aborts first.

### Not checkpointed

Like D2b's `nozzle_z` and T_rec: `jet_fs_z0`, `foot_clear_k`, and the T_rec
source.

### Far-field h law (helpers only)

- `walljet.h_expr(h_ref, far, n)` / `h_py(...)`:
  - `far="clamp"` (the default) is D2b's string, byte-identical to
    `sp_meier_pilot/input_feet` (unit test (j));
  - `far="power"` multiplies by `pow(0.09/max(0.09,s),n)`, with numeric
    arguments first.
- `feetmodel.Cfg` carries the far law, the momentum D_e and the core length.
  The default `Cfg()` reproduces `feetmodel.py`'s D2b output byte for byte.

## Study files

- **`nozzle.py`:** the isentropic Laval state and D_e at T0 = 1600/1750/1900 K,
  with γ = 1.25/1.3/1.35. J includes the pressure thrust, −0.18 N at 1900 K;
  the packet's 19.3 N is ṁ·u_e alone.
- **`handmodel.py`:** steady state (`feetmodel` with the D2c closures) and a
  transient burner-on-feet model. The transient covers time to a steady
  window, the mouth with dilution, the whole excavation and the mechanical
  share. See the docstring.
- **`predict.py` → `PREDICTIONS.md`:** P1–P6, the model check against D2b,
  and the domain-height rule. The sha256 is in `RESULTS.md` §0.
- **`score.py`:** the ACTIVE_STEP scoring table. It reuses
  `d2b_feet_rop/analyze.py` and writes `output/<run>_fig88_profile.csv`. The
  definitions are in its docstring. Two interpretations:
  - the mouth passes iff 0 ≤ Ø − width ≤ 10 mm at 25 and 50 mm;
  - any Fig. 8.8 row inside the drilled depth with Ø < width fails the
    profile, because the width is a lower bound.
- **`run.py`:** the D2d matrix, watchdogs (stalled, bottom, steady-stop,
  cap), probes, projection and `--write-input`. Metadata uses `wall_s` for
  wall-clock seconds; D2b's `.done` files overwrote `wall` with it.

## Scored configuration (fixed; not to be changed after any output)

- J-M, T_nozzle 1900 K (1600 and 1750 K are the sweep).
- Exhaust recirculation, with entrained mass.
- Momentum D_e, core 8.
- Far law power with n = 1.
- Free surface with aspect 1 and exponent 1.
- Exponential bin update; rock-to-gas heat counted.
- Feet with pads at the 0.9 quantile, plus clearance.
- Wall treatment A2, `pinned_idle_cycles = 2`.
- Quarter domain 0.12 × 0.12 × 0.40 m at 2 mm (the domain-height rule; see
  PREDICTIONS P1).

**No parameter was adjusted toward Meier's ROP, hole diameter or volume
rate.**
