# D2a2: energy-conserving wall-jet T_gas(r, s)

**Purpose.** D2a's study profile made T_gas a function of stand-off only, so
face power grew with the domain (20.9–29.0 kW against a ~20 kW jet cap), and
the J-5 / J-10 hole widths came from the s ≤ 0 clamp. D2a2 replaces the gas
temperature with a conservation law. The new key is
`surface_patch.jet_closure = enthalpy`.

This is a study, not a regression, so there is no `test`. The regression is
`tests/MMWSpalling/unit/jet_enthalpy`. Results and pre-registration are in
[RESULTS.md](RESULTS.md). **No parameter was adjusted toward Meier's ROP, hole
diameter or volume rate.**

## How to run

```bash
VENV=/Users/tzetze20/Desktop/code/.venv/bin/python
$VENV walljet.py                  # arithmetic below + pre-registration numbers
$VENV run.py --part 1 --jobs 3    # 2-D slab, 18 runs, 4 ranks each -> output/p1_*
$VENV run.py --probe              # 3-D cost probe (10 s of J10_19)
$VENV run.py --part 2 --jobs 3    # 3-D quarter, 9 runs -> output/p2_*
$VENV analyze.py                  # RESULTS.md tables between the markers + PNGs
```

Everything is a CLI override of `validation/meier/sp_meier_pilot/input_2d_dev`
(Part 1) or `input_drilling` (Part 2).

The first run set implemented the first packet version, which put the nozzle
temperature at stagnation (now `jet_stagnation = nozzle`). It is kept in
[v1_nozzle/](v1_nozzle/RESULTS.md) for its runaway-centre and coherence-cap
findings.

## The closure (code: `JetEnthalpyMarch` in `src/Integrator/MMWSpalling.H`)

Steady enthalpy balance of the wall jet, marching outward in r:

    mdot cp dT_gas/dr = −q(r) · 2πr,    q = h(r, s) (T_gas − T_s)

- **Discrete form, once per step at level 0.** Live in-patch columns with
  s > 0 are binned by r (width `jet_dr`, default max(dx, dy)).
- **Stagnation (free-jet centreline decay, default `jet_stagnation = decay`).**
  - s_c is the mean stand-off of the innermost non-empty bin.
  - φ = min(1, `jet_core_length` · `jet_D` / s_c); φ = 1 for
    `jet_stagnation = nozzle`, or when no column marches (s_c = 0).
  - T_stag = T_ent + (T_nozzle − T_ent) φ.
  - The march mass is m = `jet_mdot`, or `jet_mdot`/φ with
    `jet_entrained_mass = 1`, which carries the diluted mass instead of
    discarding the excess.
  - This line is what bounds the central pit: as the centre deepens, s_c
    grows and T_stag falls. Without it, the centre runs away from any feed
    (v1).
- **March.** From the innermost bin, with T_gas = T_stag:
  - `P_bin = Σ max(0, q_robin) dx dy`, where `q_robin` is `SurfaceCellFlux`'s
    h(T_gas − T_s) evaluated with the H update's inputs (top-cell T_old and k,
    the robin form, and the C1 pin from the same pin rule at this h and T_gas);
  - the bin's columns get `flame_Tg = T_gas`;
  - `T_gas ← T_gas − min(P_bin, avail)/(m cp)`, where
    avail = m cp (T_gas − `jet_T_ent`).
- **Columns with s ≤ 0** (above the nozzle plane) get h = 0 and are skipped;
  they keep their radiation and convection losses.
- **Columns hotter than the gas** (q < 0) neither absorb nor return enthalpy
  in the march. The kernel still applies their (negative) q.
- **Invariant (abort, to 1e-9 of the nozzle budget):**
  - jet_P_face ≤ jet_P_cap = m cp (T_stag − T_ent);
  - jet_P_face + jet_P_exhaust = jet_P_cap;
  - jet_P_decay = jet_mdot cp (T_nozzle − T_ent) − jet_P_cap ≥ 0, and = 0
    unless the decay runs without entrained mass. The nozzle budget closes as
    face + exhaust + decay.
- **Overdraw.** If a bin would take more than the enthalpy left above T_ent,
  the excess is supplied by the T_ent floor, the recirculation reservoir. It is
  not counted in `jet_P_face`; it appears as `patch_P_robin − jet_P_face`
  ("floor supply" in RESULTS). This happens only when T_ent exceeds the rock
  temperature; it is zero in every run here.

**Closed forms (predictions only).** With h flat (h_ref) inside r_core,
h = h_ref r_core/r beyond, and a uniform T_s:

- the core depletes the excess by exp(−π r_core² h_ref/(m cp));
- beyond the core the excess decays as exp(−(r − r_core)/L), with
  L = m cp/(2π h_ref r_core).

`unit/jet_enthalpy` (b) verifies the discrete march against its recursion to
6e-12, and first-order convergence to this continuous solution.

### Keys

| key | meaning | default |
|---|---|---|
| `surface_patch.jet_closure` | `none` (D2a) or `enthalpy` | `none` |
| `surface_patch.jet_mdot` | jet mass flow [kg/s], **the share of the modelled sector** | required |
| `surface_patch.jet_D` | nozzle bore [m] (7.5e-3 here) | required |
| `surface_patch.jet_core_length` | potential core [nozzle diameters] | 5 |
| `surface_patch.jet_stagnation` | `decay` (centreline decay) or `nozzle` (T_nozzle at stagnation) | `decay` |
| `surface_patch.jet_entrained_mass` | 1: march with mdot/φ (diluted mass) | 0 |
| `surface_patch.jet_T_nozzle` | gas temperature entering the stagnation bin [K] | required |
| `surface_patch.jet_cp` | T-averaged cp [J/kgK] | required |
| `surface_patch.jet_T_ent` | floor / recirculation temperature [K] | 293.15 |
| `surface_patch.jet_dr` | radial bin width [m] | max(dx, dy) |
| `surface_patch.jet_profile_interval` | writes `<plot_file>_jet_profile.csv` every this many seconds (0 = off) | 0 |

- **Requirements.** `enthalpy` requires `h_expr`, forbids `T_flame_expr`, and
  is single-level.
- **New thermo columns** (after all older ones): `jet_P_face`, `jet_P_cap`,
  `jet_T_stag`, `jet_s_c`, `jet_P_decay`, `jet_T_exhaust`, `jet_P_exhaust`,
  `jet_r_reach`.
  - `jet_r_reach` is the largest column r over bins with P_bin > 0.
- **Profile CSV** columns: `time,bin,r_lo,n_cols,T_gas,P_bin,P_share`. The
  row's T_gas is the bin inlet value, which is the value its columns receive.

**The march sees only the modelled face, so `jet_mdot` must be the modelled
sector's share of the jet.**
- **Part 2** (quarter domain, axis at the corner, Neumann-0 symmetry):
  mdot/4, and every power is reported ×4.
- **Part 1** (8 mm slab through the axis): mdot · w/(π r_core) = 2.055e-3
  kg/s.
  - This is the slab's share at the core edge, where the strip area 2w·dr
    equals the annulus share.
  - Beyond r_core the slab keeps too much gas, so Part 1 is machinery and
    trend only.

## Constants and their sources

| quantity | value | source |
|---|---|---|
| nozzle bore D | **7.5 mm** | burner drawing VT5 sheet 14/17 (thesis text: 7.1 mm; the drawing wins, ROADMAP "Burner drawings") |
| stand-off SOD | 50 mm = 6.67 D | drawing sheet 1/17 (nozzle exit to foot tip) |
| r_core | 2.5 D = 18.75 mm | Martin's lower validity bound, as in D2a |
| mdot | 0.01513 kg/s | thesis measured (52 kg/h air + 2.47 kg/h CH4), not the drawing's design flows |
| T_nozzle | 1900 K / 1436 K | thesis adiabatic flame / chamber thermocouple 1163 °C, uncorrected (a lower bound) |
| cp | 1250 J/kgK | planner's T-averaged value for the products; **not checked against a table here**; ±8 % is the stated sensitivity |
| T_fire | 821 K | C1 with weibull.V0 = V_cell |
| ρCp | 2.1725e6 J/m³K | Grimsel block (2750 × 790) |
| Meier rock-side removal power | 2.83 kW | 2.47e-6 m³/s × ρCp × 528 K |
| h_ref bracket | 5e3, 1e4 W/m²K | reference bracket; 1e4 from Meier Ch. 7 (Potter Drilling, not measured) |

## Martin h at D = 7.5 mm

`walljet.py` reuses D2a's `jet.py` with `jet.D = 7.5e-3`. That covers Martin
(1977): G, F, and the local coefficient h_loc = (1/2r) d(r² h_avg)/dr, whose
constants D2a verified against Zuckerman & Lior (2006). It also covers the
Sutherland air properties.

- **`verify_hloc()` at D = 7.5 mm:** G is reproduced over 2.5–7.5 D (H/D = 2,
  7, 12) to **1.0e-12**.
- **Film temperature per case:** (T_nozzle + T_fire)/2.

| T_nozzle | T_film | μ [Pa s] | k [W/mK] | Pr | Re | F | G(2.5D, SOD) | Nu_avg | h_avg | **h_loc = h_Martin** |
|---|---|---|---|---|---|---|---|---|---|---|
| 1900 K | 1360.5 K | 4.974e-5 | 0.0805 | 0.726 | 5.16e4 | 781.3 | 0.2182 | 149.0 | 1599 | **1449 W/m²K** |
| 1436 K | 1128.5 K | 4.461e-5 | 0.0715 | 0.733 | 5.76e4 | 841.6 | 0.2182 | 161.2 | 1536 | **1392 W/m²K** |
| (1500 K, D2a's film) | 1160.5 K | 4.535e-5 | 0.0728 | 0.732 | 5.66e4 | 832.2 | 0.2182 | 159.3 | 1546 | 1400 |

- **D = 7.5 vs 7.1 mm:** D2a's h_Martin was 1527 W/m²K (D = 7.1 mm, H = 7 D,
  1500 K film). At the same film, D = 7.5 mm and SOD = 50 mm give 1400
  (−8.3 %), as the ROADMAP note predicted (h ∝ D^−1.5 at fixed mdot).
  Re = 4 mdot/(π D μ) falls from 5.98e4 to 5.66e4.
- **Film choice:** 1449 and 1392 differ by 4 %. Using each case's own film is
  consistency, not an extra knob.

**h(r, s)** keeps D2a's form at the new D:
h_ref · h_loc(max(r, 2.5D), clamp(s, 2D, 12D)) / h_loc(2.5D, SOD).
- It is flat inside r_core and follows Martin's decay beyond.
- Every `min`/`max` takes its numeric argument first: the AMReX 25.12
  `a/F2(x, number)` rewrite trap, guarded by `unit/jet_enthalpy` (g).

## Settings

- **Common:**
  - `robin_form = pinned`, `pinned_idle_cycles = 2`;
  - `surface.follow_mask = 1`, `energy_ledger.enabled = 1`,
    `spall.removal_events_csv = 1`;
  - `beam.P0 = 0`, `weibull.V0 = dz³`;
  - `nozzle_collision_radius = 2.5 D`, feed 1.5 m/h (Part 1 also feed 0),
    `nozzle_z0 = z_top + 50 mm`.
- **Time step:** the overshoot rule q dt/(ρCp dz) ≤ 5 K with the stagnation
  maximum q. dt = 4 ms halved until the rule holds: 4 / 2 / 1 ms, 0.8–3.0 K
  per step.
- **Anchor grid:** h_ref ∈ {h_Martin, 5e3, 1e4} × T_nozzle ∈ {1900, 1436} K,
  named JM_19, JM_14, J5_19, J5_14, J10_19, J10_14.
  - Baseline: `jet_stagnation = decay`, `jet_core_length = 5`, nozzle mass.
- **Part 1** (2-D slab, 4 ranks, 60 s cap): the 6 anchors at feed 0 and
  1.5 m/h. JM_19 at 1.5 m/h adds cp ×0.92 / ×1.08, T_ent = 600 K,
  jet_dr = 1 mm, `jet_stagnation = nozzle`, and `jet_entrained_mass = 1`.
- **Part 2** (3-D quarter, 60 × 60 × 100 at 2 mm, 0.12 × 0.12 × 0.20 m, axis
  at the xlo/ylo corner, 4 ranks, 300 s cap): the 6 anchors at 1.5 m/h,
  mdot/4.
  - Treatment runs at 1900 K: J5_19_ent and J10_19_ent (entrained mass), and
    J5_19_noz (nozzle stagnation).
- **Stop rule** (the packet's "300 s or 40 mm above the bottom"):
  - **Decay:** the centre leads the nozzle by at most s_eq + 10 mm, where s_eq
    is the stand-off at which the centre's closed-form ROP equals the feed.
    For feed 0, s_eq is where T_stag = T_fire and the centre stalls.
    - Part 2 stops: 300 s, except J5_19 at 253 s and J10_19 at 232 s.
  - **Nozzle:** the centre outruns any feed (v1), so the run stops at the
    centre's fastest closed-form rate: J5_19_noz at 31 s.
- **Time step (decay):** the stagnation q is evaluated at s ≥ SOD − 2 mm if
  the centre outruns the feed at SOD, else at s ≥ 5 D. That gives dt = 4 ms
  everywhere except J10_19 (2 ms).
- **Cost probe:** 10 s of J10_19 took 97 s, loaded by the witness run. That
  gives 20–38 min per Part 2 run alone, so no run was shortened.
- **Known weaknesses (stated, not fixed):**
  - no cold entrainment, which bounds the flux, and therefore the removal
    reach, from above;
  - the decay line with nozzle mass discards (1 − φ) of the excess
    (`jet_P_decay`); entrained mass is the other bracket;
  - the exhaust enthalpy is discarded, which biases the hole narrower;
  - axisymmetric radial binning;
  - column-top patch only: no side-wall flux, no 1/cos θ;
  - pin state and ledger are not checkpointed.
