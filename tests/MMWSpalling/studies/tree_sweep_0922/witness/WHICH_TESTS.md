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

| you changed | run | cost |
|---|---|---|
| the enthalpy **march** (`JetEnthalpyMarch`, `jet_P_face`, per-bin `T_gas`, `side_P_col`) | `unit/jet_enthalpy`, `unit/jet_d2c` | 239 + 189 s |
| `BuildFlameColumns` (`flame_h`, `flame_Tg`, `flame_s`, the h/T maps) | `unit/jet_enthalpy`, `unit/robin_jet` | 239 + 25 s |
| the **feet rule** / burner descent (`FeetDescent`, `foot_*`, `nozzle_descent`) | `unit/robin_feet` | 217 s |
| the **top-face Robin** / `SurfaceCellFlux` | `unit/robin_face`, `unit/robin_pinned` | 7 + 19 s |
| **removal, the production path** (`FaceRemoval.H`: the line rule, the face list, the buckling span; `Shared.H`: the scan pieces, the remap, the step prologue) | `witness/face_gate.py run <label>` + `cmp p2_24 <label>` (the face-list gate, P2-18; re-pinned at P2-19 and P2-24; P2-23 adds the flameless `BEAM_S`: 14 short legs, every one on the face list), `unit/face_removal --check`, `unit/face_list --check` | ~150 s + 60 s |
| the **flake depth of the fallback branch** and the **profile helpers** (`ProfileEFoldDepth` / `ProfileFallDepth` / `ProfilePeak`, `Shared.H`; the `nodepth` refusal, `FaceRemoval.H`) | `unit/profile_depth` (compiles the helpers' own text and runs z0a–z0d), the face gate | 3 s + ~150 s |
| **removal, the column path** (`Removal.H`, retiring; `h_c`, the flip test, `removed_mf`) | `unit/removal_increment`, `unit/scalar_flaw` (`unit/spall_event` is the DAMAGE lineage's test) | ~15 s |
| the **beam** / Beer-Lambert deposition (`BeamSource`, `BeamEnergyBalance`, `MMWBeam.H`: the dielectric description eps_r / tan_delta -> alpha, R, P2-20; the loss table `beam.tan_delta_T` and `ColumnOpticalDepth`, P2-21) | `unit/beam_void_closure`, `unit/beam` (eps_r 1 fixture, alpha 10), `unit/beam_tau` (U1 the two-layer optical-depth sum, U2 sum path = closed form on a uniform column, U3 unit/beam byte-identical), `unit/retired_keys` (R5, K3, K4: the typed alpha_expr / R abort; R6 the alpha association lock) | 13 + 5 + 25 + 3 s |
| **mesh / AMR / decomposition** plumbing | `studies/s1_surface_resolution --smoke` | 9 s |
| the **retired-key lists** in `MMWSpalling::Parse` (P2-10's flame-wall keys; P2-16's `spallation.sp.constraint_depth`; P2-18's `retired_at_default` table) | `unit/retired_keys` (source half + one-step abort legs), `unit/lateral_flux_lock` | 3 + 60 s |
| the **buckling condition** (`spall.buckle_span`, `BuckleRatio`, the face list's withdrawal of a fired line) or `ElasticConstants`' ν (the law's, P2-14a) | `unit/buckle_lock` (source; B7 locks the ν), `unit/lateral_flux_lock` (the one declared ν read), `unit/face_buckle` (P2-15: A unmoved, N blocked, witnesses, bk in metres) | 1 + 340 s |
| the **face list WITHOUT a flame** (P2-23: the surface-patch parse before the `enabled` return, the switch `fl_on`, `FaceListHeat`/`FaceListTg` at h = 0 / T_gas = 0, the flameless list half in `AdvanceMicrostructure`, the ledger's switch) and the **blocked-line record** `spall.buckle_log_csv` | `unit/face_beam` (z1: the listed faces' loss at their own T_s against the analytic; U3 the P2-22 binary refuses), `face_gate.py` (its `BEAM_S` leg is this deck), `validation/oglesby/test` on a short B0 | 25 s + the gate |
| the **melt and the vapour** (P2-24: `MeltFraction`; the beam's `beam.melt_factor` / `beam.R_melt` in `ColumnOpticalDepth` and the per-column `R_col` at every `(1 - R)` site; the molten / lid exclusions before the crack scan; the face list's vapour rule and `ledger_E_vapor`; the `bc_*` melt record; the phase table's boxcars, `Table.H`) | `unit/face_melt` (U1 the Beer–Lambert integral with R_col, the molten count; U1b the lid count; U2 the no-op keys bit-identical), `unit/face_vapor` (U3 the void at t_v and its booking; U4/U5 no void; U6 the buried count; U3 on the P2-23 binary refuses), `unit/latent_heat`, the face gate (`BEAM_S` carries the keys), `validation/oglesby/test --onset` on a short M0 | 30 + 15 s + the gate |
| the **energy ledger** (`EnergyLedger`, `ledger_err`) | `unit/robin_face`, `unit/scalar_flaw` (both check `ledger_err`; so do `robin_jet`, `robin_pinned`, `jet_d2c`, `pad_clear_ledger`) | 7 + 9 s |

⚠ **The tree gate never runs the face list** (no gate deck sets `surface_patch.face_list = 1`;
census 2026-11-05 finding 1). The production removal path is gated by `face_gate.py` (P2-18),
whose reference list is declared in `REFERENCE.md` beside the tree list. Run BOTH at the end of
a packet that touches removal. `unit/jet_d2c` runs in the tree sweep but its output is not
hashed (`hash.sh`): its legs KT/K15/A17–A19 use `jet_T_ent_mode = rec_fixed`, which retires
(census B3), so it was left run-only rather than added to the oracle.

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
