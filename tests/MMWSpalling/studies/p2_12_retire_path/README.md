# P2-12 — retire P2-11b's surface-path heat mode (2026-10-30)

Done by the advisor at the user's instruction ("delete the lines we aren't using"; no `/implement`). P2-11b
(review-accepted 2026-10-29) closed the bounded heat-model check with X1 NOT joined to the Höser deck, so its mode
(`surface_patch.face_heat = path`, `surface_patch.barrier_expr`) had no user.

## What changed

| file | change |
|---|---|
| `src/Integrator/MMWSpalling.H` | restored byte-for-byte to `studies/p2_11b_path/src_pre/MMWSpalling.H`, P2-11b's own pre-edit snapshot (the P2-10 C2 source); −435 / +12 lines |
| `src/Integrator/MMWSpalling/FaceRemoval.H` | restored to `src_pre/FaceRemoval.H` (the fallback flake reads `flame_Tg[c.col]` again) |
| `tests/MMWSpalling/unit/lateral_flux_lock/test` | the two lines P2-11b restated put back: `FrFaceFlux` called 4 times, not 5 (recovered from a run-time `diff.patch` record; identical to the reconstruction, `gates.py` docstring) |
| `tests/MMWSpalling/unit/face_heat/` | moved to `retired_units/face_heat/` (6.6 GB, mostly outputs), as P2-10 kept its retired tests |

The pre-deletion files are in `src_p2_11b/`. `src/Integrator/MMWSpalling/Removal.H` was untouched by P2-11b and is
unchanged.

## Gates (`chain.sh`, log `CHAIN.log`): ALL PASS on binary `9278aa02d8036a74`

| gate | claim | result |
|---|---|---|
| d1 | both headers byte-identical to `src_pre` | PASS |
| d1b | the binary carries no P2-11b string (`face_heat`, `barrier_expr`, …) | PASS (0; the stale binary had 5) |
| d2 | tree sweep 0 of 2639 vs `p2_10_c2` and vs `p2_11b`; controls `d2aa` 26, `d2t` 863 | PASS (`GATE_tree.txt`) |
| d3 | witness set W (M0T, 6 face-removal legs, 4 face-list legs) vs P2-10's `W_c2`, P2-10's comparator | PASS (`W_p2_12_vs_c2.txt`) |
| d4 | every kept unit test vs P2-10's `UNIT_c2`: same rc, identical verdict lines, no exception | PASS (`UNITS_p2_12_vs_c2.txt`) |
| d5 | X1's own argv aborts and names `surface_patch.face_heat` (strict `surface_patch` parse) | PASS (`RETIRED.txt`) |

Five unit tests exit rc 1 on this binary **and on P2-10's C2** (unchanged):
- `beam_void_closure` (T3, a documented debt);
- `amr_microstructure_regrid`;
- `heterogeneous_kappa`;
- `sp_melt_skeleton` (the molten-cap criterion and OFF-leg removal checks fail);
- `surface_T_recon` (event-log header gating).

They are recorded in the model card, not fixed here.

## ⚠ The stale-binary incident (`stale_binary/`)

The first gate run used a binary that still carried P2-11b:
- `cp -p` kept `src_pre`'s Oct-29 mtimes, so make judged the object up to date and did not recompile.
- d2 passed on it, because position mode was already identical.
- d5 caught it: X1's argv ran with rc 0, and `barrier_expr` on its own aborted at P2-11b's own check.

The fix was to touch both headers and rebuild, and to add gate d1b. ⚠ **Restoring a source file with `cp -p` (or any
copy that keeps an old mtime) needs a `touch` before `make`.**

## Not done here

- (done 2026-10-30 by `/plan`: P2-11b and P2-12 archived in `docs/project/ARCHIVE_DONE.md`; P2-13 is active.)
- P2-11b's records and run outputs (`studies/p2_11b_path/`, ~8.5 GB) are kept.

## ⚠ REVERSED the same day: X1 restored (user, 2026-10-30, "realistic over matching")

X1's coordinate (heat by distance along the surface from the stagnation point) is judged the more realistic of the
two flame heat rules (`docs/project/MODEL_CARD.md` §11 item 5), so it stays in the code, default-off, as the basis for
future wall-heat work. Restored from `src_p2_11b/` (3 files, byte-identical); `unit/face_heat` moved back; rebuilt with
`touch` first. Binary `0a96a095ad30f45e`.

Gates (`restore/chain_restore.sh` + `chain_restore2.sh`, log `restore/CHAIN.log`): **ALL PASS**

| gate | claim | result |
|---|---|---|
| r1 | the three files byte-identical to `src_p2_11b/`; `unit/face_heat` present | PASS |
| r1b | the binary carries the X1 strings | PASS (5) |
| r2 | tree 0 of 2639 vs `p2_11b` and `p2_10_c2`; controls 26 / 863 | PASS |
| r5 | X1's own argv to 15 s vs the original X1 leg: plotfiles 00000/01875 byte-identical, thermo 313 rows × 92 non-timer columns identical, event log (8,539 rows) and jet profile exact prefixes | PASS |
| r3 | W vs `W_c2` | PASS |
| r4 | all 38 kept unit tests vs P2-11b's `UNIT_p2_11b` (face_heat included): same rc, identical verdict lines | PASS |

⚠ r5's first two scorings failed on MY window definitions, not on the restore, and both fixes are declared in
`gates_restore.py`:
1. "time < 15 s" counted 2 original rows stamped 14.999999999999234: the step that STARTS at 15 s, which a run
   stopped at 15.0 never executes.
2. The jet profile is logged once per second, so its last row before 15 s is 14.008 s.

The claim is exact-prefix identity.

**Net effect of P2-12: zero source change.** What remains is the record, two lessons (touch after a `cp -p`
restore; compare logs as exact prefixes), and the model card built in the same session.
