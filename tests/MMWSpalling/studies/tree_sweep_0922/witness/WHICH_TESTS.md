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

| you changed | run | cost |
|---|---|---|
| the enthalpy **march** (`JetEnthalpyMarch`, `jet_P_face`, per-bin `T_gas`, `side_P_col`) | `unit/jet_enthalpy`, `unit/jet_d2c` | 239 + 189 s |
| `BuildFlameColumns` (`flame_h`, `flame_Tg`, `flame_s`, the h/T maps) | `unit/jet_enthalpy`, `unit/robin_jet` | 239 + 25 s |
| the **feet rule** / burner descent (`FeetDescent`, `foot_*`, `nozzle_descent`) | `unit/robin_feet` | 217 s |
| **side faces** (`side_face_*`, `k_face`'s lateral branch, `WallStepDiagnostic`) | `unit/side_face_flux`, `unit/side_face_area`, `unit/side_face_step_height` | ~60 s |
| the **top-face Robin** / `SurfaceCellFlux` | `unit/robin_face`, `unit/robin_pinned` | 7 + 19 s |
| **removal** (`Removal.H`, `h_c`, the flip test, `removed_mf`) | `unit/spall_event`, `unit/removal_increment`, `unit/scalar_flaw` | ~15 s |
| the **beam** / Beer-Lambert deposition | `unit/beam_void_closure` | 13 s |
| **mesh / AMR / decomposition** plumbing | `studies/s1_surface_resolution --smoke` | 9 s |
| the **energy ledger** (`EnergyLedger`, `ledger_err`) | any side-face test — `ledger_err` is checked in all of them | — |

⚠ `unit/beam_void_closure` **exits 1 by design**: its T3 is D2z-1's documented
debt. `rc = 1` from that test alone is the expected state, not a regression.

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
