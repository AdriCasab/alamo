# Which tests cover which site — the per-build check

**The full gate runs ONCE PER PACKET, at the end, on the binary that ships.** Its
claim is a property of that binary: *with every new key off, it reproduces the
frozen tree*. That claim is not made per edit, so it does not need proving per
edit.

During development, run the two or three tests that cover the mechanism the
packet touches. The packet's reachability ledger already names the sites; this
table turns those into a **~2-minute check per build** instead of a 5–13 minute
one.

⚠ **This does not replace the gate and must not be used as a "fast tier".** The
gate's value is catching what nobody anticipated; this table only covers what
somebody did. Tiering the gate itself by cost would skip its slow legs (HOS_S, CAL_S,
FEET_S) -- exactly the jet, closure and feet paths -- so it selects against the gate's
whole purpose.

**Every kept test, from its own deck: `run_kept.py` (P2-29) -- the USER's periodic check, every few days, not a packet
step** (the user, 2026-11-10). `python run_kept.py run <label>` archives each output, launches every deck, scores, and
fails on a stale output or on a test neither in its manifest nor in its EXCLUDED list (`kept_<label>.txt`, one line per
test). A packet proves the runner on `run <label> --subset --plant-stale` (latent_heat, beam, stefan, face_dirichlet,
~10 min). ⚠ P2-28's b4 scored two tests against outputs written days earlier -- a test that only reads an output is not
a test of the binary that ships.

## Site → test

P2-10 deleted the side-face family and its six tests (`side_face_flux`, `side_face_area`,
`side_face_step_height`, `side_staircase`, `lateral_drive`, `wall_march`); their row is gone.
P2-16 deleted the constraint depth; its two tests (`constraint_ref_lock`, `face_cref`) moved to
`tests/MMWSpalling/retired_units/` and their row is replaced by the retired-key row.
P2-26 deleted the closed switches (census B3, C2, C3): `unit/sp_melt_skeleton` (the `melt_aware` A/B) moved to
`retired_units/` (README there); `removal_increment` lost its SN leg and checks (e)–(f); `jet_d2c` lost the
`rec_fixed` legs KT, K15, A17–A19 and check (k); `face_heat`'s R_pad / R_above are refused by the strict
`surface_patch` scope now (the keys are gone), not by P2-12's message. Restated in the same packet because they
were already RED when it ran them (each pre-dated it): `lateral_flux_lock`'s `FaceListHeat` count (2 since P2-23's
flameless list, exclusive branches) and `buckle_lock` B2 (P2-23's log-only `bk_lb` call); `face_removal --check`
reads FAIL (4) on its 2026-10-27 outputs, byte for byte the P2-18/P2-19 lock logs -- the bar is "unchanged".
P2-27 (inputs + census C1a/C1b + G, `studies/p2_27_cleanup/`): NEW `unit/interlocks` (the cross-key aborts moved into
`CheckInterlocks`, one-rule legs I1-I20 compared with the P2-26 binary's messages under `--pre`, and the per-level
explicit dt check D1-D5), `unit/beam_schedule` (`beam.schedule.file` == the arrays, monotone t_end both forms) and
`unit/face_feet` (the pad clearance on the face list). RETIRED to `retired_units/`: `dp_yield`, `gb_cohesive`,
`spall_event` (left the tree list: 8 paths), `regime_low_high_power`, `heterogeneous_kappa`; the Hu damage-law
validation (`hu_end_to_end`, `hu_spall_onset`'s damage decks + `diag/`, the two damage-law extensions) to
`retired_validation/hu_damage_law/`; `hu_spall_onset/test_lefm` + `input_granite2_lefm` re-express Hu on the LEFM rule
(a RECORD). `amr_microstructure_regrid` lost its D / kappa_eff / k_eff checks (fields gone); it still reads rc 1 on
`is_grain_boundary` at level 0 (pre-existing, recorded before and after). `surface_T_recon` passes again (`,q_ratio`
added to its header constant). The tree list gained `unit/jet_d2c/output` (270 paths) and the face gate `FEET_S` /
`FEET_0` (the feet deck, clearance on / off).

**P2-26: outside-the-gate decks whose reference changed BY DESIGN** (the deleted line was live). Measured before/after
-- the pre-image deck on the P2-25 code vs the edited deck on the P2-26 binary, same argv
(`studies/p2_26_closed_switches/AB.txt`, field by field in `AB_FIELDS.txt`):
- `unit/seed_void` (surface_normal 1, coherence 0.004): only the scan's thickness report `h_spall_field` moves (39
  cells, ≤ 0.34 mm); the removal state is identical (1,938 removed cells both sides). Its test PASSES.
- `validation/meier/sp_meier_pilot/input_energy` (coherence 0.004): the cap was BINDING -- 4,648 → 12,100 removed cells
  at 4 s. `input_drilling` / `input_drilling_25mm` identical over 1 s / 2 s prefixes; `input_2d_dev_melt` thickness
  reports only over 15 s. The test, scored on both sides (`meier_ab.py`, `MEIER_AB.txt`): case A (energy) FAILS on BOTH binaries (pre-existing: A1 dV/dt 0.00 cm³/s both; A2 TSE at 10 s 10.22 → 3.93 J/mm³ vs 15.4; the hole saturates by 3 s at 37.2 → 96.8 cm³); case B (drilling) PASSES on both: ROP 1.44 m/h both (the bit sets it), V(56 s) 100.6 → 122.9 cm³, B4 TSE 6.35 → 5.20 J/mm³, Tmax 2268 → 2078 K.
- `extensions/hu_spall_onset_mmwbeam_lefm` (coherence 0.002): binding on `input_pulsed_brittle_400w` (721 → 1,526
  removed cells, full 30 s), `_hires` (218 → 495 at 8 s), `input_wellbore_3d` (351 → 498 at 8 s); `_fine`,
  `_smallgrain` identical over 0.5 s. No test runs these decks (its `test` uses `input_granite2`): record only.
- INERT, identical over the full run (the verdict is unchanged by identity): `unit/gradient_blind`,
  `unit/surface_T_recon`, `validation/p0_steady_column`, `validation/kant/sp_kant_onset_meiercfg` (both decks).
  ⚠ `unit/surface_T_recon` FAILS its check (3) on BOTH binaries: the expected removal-log header predates P1's
  `q_ratio` column (pre-existing, not P2-26's; its deck's run is byte-identical before/after).

| you changed | run | cost |
|---|---|---|
| **removal** (`FaceRemoval.H`: `RemovalStep`, the line rule, the face list, the buckling span, the evaluate-only mode; `Shared.H`: the scan pieces, the remap, the step prologue) | `witness/face_gate.py run <label>` + `cmp p2_28 <label>` (THE gate since P2-28: 19 short legs, every one on the face list), `unit/face_removal --check`, `unit/face_list`, `unit/lateral_flux_lock` + `unit/buckle_lock` (source locks) | ~5 min + 60 s + ~10 s |
| the enthalpy **march** (`JetEnthalpyMarch`, `jet_P_face`, per-bin `T_gas`; since P2-28 it bills the list's faces only) | `unit/face_heat` (v)/(vi), the face gate (M0S, HOS_S, CAL_S) | ~1 min + the gate |
| `BuildFlameColumns` (`flame_h`, `flame_Tg`, `flame_s`, the h/T maps; the constant-patch branch's list half, P2-28) | `unit/face_heat`, `unit/face_dirichlet`, the face gate | ~1 min + 10 s + the gate |
| the **feet rule** / burner descent (`FeetDescent`, `foot_*`) and the **pad clearance** (`FrFootClearance`) | `unit/face_feet` (P / N / P3 + PT: pad rule, booking, freeze, ranks), the face gate's `FEET_S` / `FEET_0` | ~3 min |
| the **prescribed patch** (P2-28: `FaceEntry::T_presc`, the half-cell Dirichlet term `2 kc (T - Tc)/dz^2` in the kernel / ledger / diagnostics, the criterion's T_s on a prescribed face) | `unit/face_dirichlet` (the exact solution at two meshes, F == C, the ledger), `validation/hu/hu_spall_onset/test_lefm` (T_s == T_B1(time) on every row), the face gate's `HU_S` | 10 s + ~10 s + the gate |
| the **evaluate-only mode** (P2-28: the list with spall and vapour off; `list_state()` / `evaluate_only()`) | `validation/kant/sp_kant_onset/test_face` (Eq. 15, h 150 / 250), `validation/kant/sp_kant_pressure_sweep/test_face` (R1, R2, R4), the face gate's `KV_S` | ~20 s + ~30 s |
| the **top-face Robin** / `SurfaceCellFlux` (cell / face / surface forms; the pinned form retired in P2-28) | `unit/robin_face` (a), `unit/face_beam` z1 | 7 + 25 s |
| the **flake depth of the fallback branch** and the **profile helpers** (`ProfileEFoldDepth` & co., `Shared.H`) | `unit/profile_depth`, the face gate | 3 s + the gate |
| the **beam** / Beer-Lambert deposition (`BeamSource`, `BeamEnergyBalance`, `MMWBeam.H`, the loss table, `ColumnOpticalDepth`) | `unit/beam`, `unit/beam_tau`, `unit/face_beam`, `unit/face_melt`, `unit/face_vapor`, `unit/retired_keys` (R5, R6, K3, K4) | 5 + 25 + 25 + 30 + 15 s |
| **decomposition / AMR** plumbing | `face_gate.py rankcheck <label>`, `unit/face_feet` (g), `unit/amr_microstructure_regrid` (rc 1 recorded since P2-27) | ~1 min |
| the **cross-key interlocks** (`CheckInterlocks`) and the **explicit time-step check** (`CheckTimeStepLimit`; 2d+1 on a prescribed top) | `unit/interlocks` (re-based at P2-28 on the face-list fixture; 15 one-rule legs, D1-D5) | ~40 s |
| the **beam power schedule** (`MMWBeam.H`: `schedule.file` / the arrays) | `unit/beam_schedule` (+ the face gate's `BEAM_S`) | ~17 s |
| the **retired-key lists** in `MMWSpalling::Parse` (P2-10, P2-16, P2-27, P2-28's presence lists; the `retired_at_default` table) | `unit/retired_keys` (source half R1-R11 + one-step K legs; P2-28: R11 + K32-K44, `--pre <stage-2 pre-image> --red R11`), `unit/lateral_flux_lock` (R11's twin) | 3 + 90 s |
| the **buckling condition** (`spall.buckle_span`, `BuckleRatio`, the withdrawal) or `ElasticConstants`' nu | `unit/buckle_lock`, `unit/lateral_flux_lock`, `unit/face_buckle` | 1 + 340 s |
| the **face list WITHOUT a flame** (P2-23) and the **blocked-line record** `spall.buckle_log_csv` | `unit/face_beam`, the face gate's `BEAM_S`, `validation/oglesby/test` on a short B0 | 25 s + the gate |
| the **melt and the vapour** (P2-24) | `unit/face_melt`, `unit/face_vapor`, `unit/latent_heat`, the face gate (`BEAM_S`) | 30 + 15 s + the gate |
| the **energy ledger** (`EnergyLedger`; `ledger_err_raw` -- renamed P2-28, not a closure on a spalling deck -- and `ledger_err_solid`, the closing row) | `unit/face_feet` (e), `unit/face_vapor`, `unit/face_heat` (v), `unit/face_dirichlet` (d), `unit/robin_face` (a) | ~3 min |

⚠ **P2-28 (census H): THE COLUMN REMOVAL PATH IS DELETED** -- one removal path, the face list. The tree gate
(`hash.sh`, `post_all_par.sh`, `guarded_post_all_par.sh`, `tree_compare.py`) is a RECORD: `p2_27.md5` (2,901 paths) is the
column era's last list, and 2,225 of its paths died with the column path (the 676 `robin_face/form_*` that survived are
no-removal cell/face Robin controls, covered by `unit/robin_face` (a)). The face gate is THE gate (`REFERENCE.md`).
Retired with the column path (`retired_units/README.md`, `retired_validation/column_path/README.md`): `robin_pinned`,
`robin_jet`, `jet_enthalpy`, `robin_feet`, `jet_d2c`, `pad_clear_ledger`, `pin_clock`, `removal_increment`, `scalar_flaw`,
`beam_void_closure`, `seed_void`, `surface_T_recon`, `sp_v_n_regime`, `sp_lefm`, `robin_face`'s ledger legs, the column
`gradient_blind` deck; site A's Kant onset decks and tests, the pressure-sweep site-A deck, `sp_kant_onset_meiercfg`,
`p0_steady_column`, `rossi`. COVERAGE GAPS NAMED, NOT CLOSED: the jet march's per-bin recursion vs the strip closed form
(was `jet_enthalpy`); the feet descent's feed cap / stall / foot_rule modes (was `jet_d2c`); a beam path across an
interior void (was `beam_void_closure`).

**P2-27: the decks whose t = 0 plotfile changed BY DESIGN (stage 4).** The six `microstructure.phase0.{kappa, beta,
E, mu, rho, Cp}` are refused under `law = grimsel` (`E, beta, mu` only when `spall.stress_law = grimsel` too), so a law
phase's t = 0 property fields come from `UpdatePropertyFields` (the law at the start temperature) instead of the typed
scalars: at 293.15 K `kappa_phase` 1.5 -> 3.43053 W/mK, `beta_phase` 8.0e-6 -> 8.4e-6 1/K, `E_phase` 30 -> 53.3 GPa,
`mu_phase` 11.54 -> 20.5 GPa (at 297.5 K, BEAM_S: 3.40512, 8.41631e-6, the same E and mu); `rho_phase`, `Cp_phase`
were already law-filled and did not move; every later plotfile and every thermo column is identical (the field audit,
`AUDIT_4_fields.txt`). Out of the gate, every deck that includes the rock file or the two edited fixtures moves the same
way at t = 0 only; `face_list`, `face_melt`, `face_vapor`, `face_beam` pass on the new binary.
⚠ `validation/kant/sp_kant_onset_meiercfg/input_weibull:103` sets `weibull.V0 = 8e-9` on a 1 mm mesh (V0 != V_cell, a
legacy record, not edited): the only deck where the Weibull size correction is not neutralised (MODEL_CARD §6 row 11).

⚠⚠ **EVERY beam fixture and gate deck declares `beam.eps_r = 1`, so NEITHER gate covers anything that reads the
refractive index.** P2-22 proved it the hard way: a build computing the incident power with a 0.55 m wide beam
passed every unit test and would have passed both gates (the `-1` "no void path" sentinel had been fed into the
in-rock spread formula); only a smoke run on a granite deck caught it. A packet touching an index-dependent
quantity owes a fixture with `eps_r != 1`. Related: byte identity on an `eps_r != 1` deck does not cross the P2-22
boundary -- see `REFERENCE.md`'s P2-22 block before re-running P2-20's or P2-21's a3 as a regression check.

## The other saving, which is larger than the gate

Every study carries an `S0` identity leg on the Meier geometry at **500 s**
(12–13 min) — as long as the whole serial gate. **A 100 s prefix proves the same
thing**, and the prefix rule is already accepted practice (D2x, and D2z's
criterion (f): *"a truncated run is a PREFIX — compare only where both have
data"*).

Measured on `d2aa_heatdelivery/output/S0`, so this is not an argument from
principle:

| mechanism | first non-zero |
|---|---|
| `patch_P_robin`, `jet_P_face`, `pinned_cols_face`, `nozzle_z`, `foot_cols`, `jet_s_c`, `ledger_E_H` | **0.10 s** (first step) |
| `side_P_face`, `side_faces` | **3.26 s** |

By 100 s there are already **2332 side faces** and **20,535 removal rows — 22.9 %
of the whole 500 s log**. And byte-identity diverges at the *first* differing
floating-point operation: D2z-2's two conductance forms parted at **row 343,
t = 16.08 s** of a 250 s leg, 6.4 % in.

⇒ **Run `S0` to `stop_time = 100` by default** (one key on the leg's argv; the
comparison is already prefix-aware). **Reserve the full 500 s for a packet that
touches something which first activates late** — anything keyed to crater depth,
to the pit passing the decay core length, or to a late free-surface transition.

⚠ Not applied retroactively to D2aa, whose `S0` ran the full 500 s and is
already scored.
