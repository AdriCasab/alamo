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
somebody did. Tiering the gate itself by cost would skip `jet_enthalpy`,
`robin_feet` and `jet_d2c` — which are exactly the jet, feet and closure paths —
so it selects against the gate's whole purpose.

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
| the enthalpy **march** (`JetEnthalpyMarch`, `jet_P_face`, per-bin `T_gas`, `side_P_col`) | `unit/jet_enthalpy`, `unit/jet_d2c` | 239 + 189 s |
| `BuildFlameColumns` (`flame_h`, `flame_Tg`, `flame_s`, the h/T maps) | `unit/jet_enthalpy`, `unit/robin_jet` | 239 + 25 s |
| the **feet rule** / burner descent (`FeetDescent`, `foot_*`, `nozzle_descent`) | `unit/robin_feet` | 217 s |
| the **top-face Robin** / `SurfaceCellFlux` | `unit/robin_face`, `unit/robin_pinned` | 7 + 19 s |
| **removal, the production path** (`FaceRemoval.H`: the line rule, the face list, the buckling span; `Shared.H`: the scan pieces, the remap, the step prologue) | `witness/face_gate.py run <label>` + `cmp p2_24 <label>` (the face-list gate, P2-18; re-pinned at P2-19 and P2-24; P2-23 adds the flameless `BEAM_S`: 14 short legs, every one on the face list), `unit/face_removal --check`, `unit/face_list --check`, and (P2-27) `unit/lateral_flux_lock` + `unit/buckle_lock`: both source locks read `FaceRemoval.H` and went RED unseen from P2-23 to P2-26 because no face-list edit was routed to them | ~150 s + 60 s + ~5 s |
| the **flake depth of the fallback branch** and the **profile helpers** (`ProfileEFoldDepth` / `ProfileFallDepth` / `ProfilePeak`, `Shared.H`; the `nodepth` refusal, `FaceRemoval.H`) | `unit/profile_depth` (compiles the helpers' own text and runs z0a–z0d), the face gate | 3 s + ~150 s |
| **removal, the column path** (`Removal.H`, retiring; `h_c`, the flip test, `removed_mf`) | `unit/removal_increment` (W, PRESC; the SN leg retired at P2-26), `unit/scalar_flaw` (`unit/spall_event` is the DAMAGE lineage's test) | ~15 s |
| the **beam** / Beer-Lambert deposition (`BeamSource`, `BeamEnergyBalance`, `MMWBeam.H`: the dielectric description eps_r / tan_delta -> alpha, R, P2-20; the loss table `beam.tan_delta_T` and `ColumnOpticalDepth`, P2-21) | `unit/beam_void_closure`, `unit/beam` (eps_r 1 fixture, alpha 10), `unit/beam_tau` (U1 the two-layer optical-depth sum, U2 sum path = closed form on a uniform column, U3 unit/beam byte-identical), `unit/retired_keys` (R5, K3, K4: the typed alpha_expr / R abort; R6 the alpha association lock) | 13 + 5 + 25 + 3 s |
| **mesh / AMR / decomposition** plumbing | `studies/s1_surface_resolution --smoke` | 9 s |
| the **cross-key interlocks** (`CheckInterlocks`, P2-27: every cross-key abort of Parse, called LAST) and the **explicit time-step check** (`CheckTimeStepLimit`, per level; face list ABORTS, elsewhere WARNS) | `unit/interlocks` (with `--pre <old binary>` the messages are compared) | ~40 s (80 with `--pre`) |
| the **beam power schedule** (`MMWBeam.H`: `schedule.file` / the arrays, monotone `t_end`) | `unit/beam_schedule` (+ the face gate's `BEAM_S`, which reads oglesby's file) | ~17 s |
| the **pad clearance** (regime 4: `FrFootClearance` on the face list, the column block in `Removal.H`; `FeetDescent` arms both) | `unit/face_feet` (P / N / P3 + a PT re-run: pad rule, booking, freeze, ranks) and the face gate's `FEET_S` / `FEET_0`; the column path: `unit/pad_clear_ledger`, `unit/jet_d2c` (CLD) | ~3 min + 1 min |
| the **retired-key lists** in `MMWSpalling::Parse` (P2-10's flame-wall keys; P2-16's `spallation.sp.constraint_depth`; P2-18's `retired_at_default` table, which P2-26 extended with eight keys) | `unit/retired_keys` (source half + one-step abort legs; P2-26: R7 the B3 keys, R8 the C2/C3 keys, K5–K13), `unit/lateral_flux_lock` (P2-26: no B1/B3/C2/C3 identifier left in `src/`) | 3 + 60 s |
| the **buckling condition** (`spall.buckle_span`, `BuckleRatio`, the face list's withdrawal of a fired line) or `ElasticConstants`' ν (the law's, P2-14a) | `unit/buckle_lock` (source; B7 locks the ν), `unit/lateral_flux_lock` (the one declared ν read), `unit/face_buckle` (P2-15: A unmoved, N blocked, witnesses, bk in metres) | 1 + 340 s |
| the **face list WITHOUT a flame** (P2-23: the surface-patch parse before the `enabled` return, the switch `fl_on`, `FaceListHeat`/`FaceListTg` at h = 0 / T_gas = 0, the flameless list half in `AdvanceMicrostructure`, the ledger's switch) and the **blocked-line record** `spall.buckle_log_csv` | `unit/face_beam` (z1: the listed faces' loss at their own T_s against the analytic; U3 the P2-22 binary refuses), `face_gate.py` (its `BEAM_S` leg is this deck), `validation/oglesby/test` on a short B0 | 25 s + the gate |
| the **melt and the vapour** (P2-24: `MeltFraction`; the beam's `beam.melt_factor` / `beam.R_melt` in `ColumnOpticalDepth` and the per-column `R_col` at every `(1 - R)` site; the molten / lid exclusions before the crack scan; the face list's vapour rule and `ledger_E_vapor`; the `bc_*` melt record; the phase table's boxcars, `Table.H`) | `unit/face_melt` (U1 the Beer–Lambert integral with R_col, the molten count; U1b the lid count; U2 the no-op keys bit-identical), `unit/face_vapor` (U3 the void at t_v and its booking; U4/U5 no void; U6 the buried count; U3 on the P2-23 binary refuses), `unit/latent_heat`, the face gate (`BEAM_S` carries the keys), `validation/oglesby/test --onset` on a short M0 | 30 + 15 s + the gate |
| the **energy ledger** (`EnergyLedger`, `ledger_err`) | `unit/robin_face`, `unit/scalar_flaw` (both check `ledger_err`; so do `robin_jet`, `robin_pinned`, `jet_d2c`, `pad_clear_ledger`) | 7 + 9 s |

⚠ **The tree gate never runs the face list** (no gate deck sets `surface_patch.face_list = 1`;
census 2026-11-05 finding 1). The production removal path is gated by `face_gate.py` (P2-18),
whose reference list is declared in `REFERENCE.md` beside the tree list. Run BOTH at the end of
a packet that touches removal. `unit/jet_d2c/output` is in the oracle since P2-27 (270 paths; its retired legs'
stale outputs were moved to `studies/p2_27_cleanup/jet_d2c_retired_legs/` first, or they would have been hashed as
constants).

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

⚠ `unit/beam_void_closure` **exits 1 by design**: its T3 is D2z-1's documented
debt. `rc = 1` from that test alone is the expected state, not a regression.

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
