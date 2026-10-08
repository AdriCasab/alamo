# P1 — MANIFEST (criterion (i))

Binary under test: `e485a62912d14bda4368c0b21beb756c437518946003d005647b537276ae39b6`
(`bin/mmwspalling-3d-g++`, built 2026-10-01 10:42:57). **Verified unchanged after
every scored run below.**

⚠ The toolchain is NOT byte-reproducible: rebuilding from identical sources gives
a different sha256, so the only evidence a leg ran on this binary is that nothing
was rebuilt after it. Staleness is decided by comparing each event log's mtime
against the binary's.

`rows` counts lines including the header, so a 1-row file is a header-only log —
that is the **expected** result for the `differential` legs of the slab test.

## (b) the uniform-slab test — `tests/MMWSpalling/unit/gradient_blind/`

harness: `tests/MMWSpalling/unit/gradient_blind/test`, output under `.../output/`

| leg | event-log sha256 (16) | rows | vs final binary |
|---|---|---|---|
| `abort` | `(no events — abort leg)` | 0 | ok |
| `absolute` | `48ba07a03b90ad6b` | 81 | ok |
| `diff_fixed05` | `0a14c7924128448e` | 1 | ok |
| `diff_fixed10` | `0a14c7924128448e` | 1 | ok |
| `diff_fixed20` | `0a14c7924128448e` | 1 | ok |
| `diff_flake` | `0a14c7924128448e` | 1 | ok |

`abort` is the parse-guard leg: `driving_ref = fixed` with no `driving_ref_depth`
must abort, so it has no event log by construction.

## (a) (d) (d′) — the 1-D column matrix, `tests/MMWSpalling/validation/p0_steady_column/`

harness: `tests/MMWSpalling/validation/p0_steady_column/test`

| leg | event-log sha256 (16) | rows | vs final binary |
|---|---|---|---|
| `B3_0.5` | `f9c7f6588fce4214` | 18317 | ok |
| `B3_1.0` | `7937dbe48e493e32` | 4429 | ok |
| `B3_2.0` | `e8f271b51d8edd49` | 1101 | ok |
| `B5_0.5` | `16e4dc337f30dc47` | 23297 | ok |
| `B5_1.0` | `2e55e3d3dfd18ddc` | 5121 | ok |
| `B5_2.0` | `24db69496d01ac47` | 1361 | ok |
| `C_0.5` | `0b5c814cf8b3cf7b` | 23381 | ok |
| `C_1.0` | `1fc2f776c83f1b9c` | 7025 | ok |
| `C_2.0` | `1ab9b17c0dd3ae65` | 1757 | ok |
| `D_0.5` | `5c4ce6da85fea466` | 17409 | ⚠ **STALE — predates the final binary** |
| `D_1.0` | `3b64f99ee07a28df` | 2241 | ⚠ **STALE — predates the final binary** |
| `D_2.0` | `58538e1c3840155c` | 289 | ⚠ **STALE — predates the final binary** |
| `F_0.5` | `fe66d0a1ada55d10` | 20949 | ok |
| `F_1.0` | `0b887cd2b3223910` | 4285 | ok |
| `F_2.0` | `5fa7a9717457fab1` | 833 | ok |
| `P_0.5` | `ef7773302d58f4a5` | 21453 | ok |
| `P_1.0` | `1aa31a7adeccc84e` | 5785 | ok |
| `P_2.0` | `bc9e2870fdea5780` | 1697 | ok |
| `S_0.5` | `b2a53d0d4934cfac` | 21249 | ok |
| `S_1.0` | `f084344ad7dc6483` | 4741 | ok |
| `S_2.0` | `25d92ccc198406ab` | 1281 | ok |

### ⚠⚠ Revision 5 — the (d) control arm, 6 new legs

Added to `VARIANTS` as `("A80"/"A08", "surface", 0.4, "absolute", a0)`. argv is
the matrix argv with `spall.driving_stress=absolute` and `spallation.sp.a0` set
as shown; everything else identical to leg `S`.

| leg | `spallation.sp.a0` | event-log sha256 (16) | rows | vs final binary |
|---|---|---|---|---|
| `A08_2.0` | 8.0e-6 | `47068feac77c7389` | 1221 | ok |
| `A08_1.0` | 8.0e-6 | `f6ce2a95ec341ebc` | 5021 | ok |
| `A08_0.5` | 8.0e-6 | `c89784823a6d73b6` | 21465 | ok |
| `A80_2.0` | 80.0e-6 | `258c617b2349f7e6` | 5401 | ok |
| `A80_1.0` | 80.0e-6 | `7476bda79d7e8dab` | 25853 | ok |
| `A80_0.5` | 80.0e-6 | `c280f9fb9db84352` | 108073 | ok |

⚠ `A08_2.0`'s hash is **byte-identical to the standalone smoke run** made before
the harness edit, which confirms the refactored `VARIANTS` emits the argv that
was smoked.

⚠⚠ **The 18 pre-existing legs were NOT re-run and their hashes are unchanged**
(spot-checked `S_2.0` `25d92ccc198406ab`, `S_1.0` `f084344ad7dc6483`, `S_0.5`
`b2a53d0d4934cfac`, `C_2.0`, `F_2.0`, `P_2.0`, `B3_2.0`, `B5_2.0` — all match the
rows above). The harness gained `--only <PREFIX>` for exactly this reason:
`run()` `rmtree`s what it runs, so adding a control arm without it would have
**destroyed the event logs this manifest hashes**.

⚠⚠ **`D_0.5` / `D_1.0` / `D_2.0` are STALE** — their hashes are character-for-
character the P0b rev7 values (`5c4ce6da85fea466`, `3b64f99ee07a28df`,
`58538e1c3840155c`; see `../p0b_one_surface_T/MANIFEST.md`), and they were
already flagged stale there. `D` is the `surface_T = cell_average` diagnostic
control. **No P1 number is scored on it**; it is quoted only as the printed
control, and the P1 keys are inert on it in any case because the reconstruction
it would modify is switched off. Do not quote these three rows as P1 evidence.

## (c) (c′) (w) — the Kant anchor, 22 legs, `tests/MMWSpalling/validation/kant/sp_kant_onset_meiercfg/`

harness: `tests/MMWSpalling/validation/kant/sp_kant_onset_meiercfg/test`
(7 families × 3 meshes + 1 parse guard; `R` reconstructed/cell-form,
`S` surface-form, `C` diagnostic control, `W` Weibull)

| leg | event-log sha256 (16) | rows | vs final binary |
|---|---|---|---|
| `C150_0.5` | `b25154f2a0c15488` | 69121 | ok |
| `C150_1.0` | `a22894f3bba49bc8` | 8641 | ok |
| `C150_2.0` | `e78f62d97749b005` | 1009 | ok |
| `C500_0.5` | `4194071a47a152bc` | 108289 | ok |
| `C500_1.0` | `555b42ab553c7313` | 13249 | ok |
| `C500_2.0` | `20aea95321295941` | 1585 | ok |
| `D150_2.0` | `(no events — abort leg)` | 0 | ok |
| `R150_0.5` | `d7d637c0a6b1e599` | 13825 | ok |
| `R150_1.0` | `53755c0ba277cc7b` | 2305 | ok |
| `R150_2.0` | `c72ecdffc30133b4` | 433 | ok |
| `R500_0.5` | `7c1dd6a1b425ef87` | 108289 | ok |
| `R500_1.0` | `0afbebf20736835f` | 13249 | ok |
| `R500_2.0` | `a035a2f3d01e4b1e` | 1585 | ok |
| `S150_0.5` | `d42a0c16f12a3ecb` | 13825 | ok |
| `S150_1.0` | `e66c417195fc56f9` | 1729 | ok |
| `S150_2.0` | `5d3a2b9666bc5ad9` | 289 | ok |
| `S500_0.5` | `22e5f55fcef1a086` | 108289 | ok |
| `S500_1.0` | `3d73f90ab4237c1d` | 13249 | ok |
| `S500_2.0` | `9509030469f8d417` | 1585 | ok |
| `W150_0.5` | `43c5bf0aa67b542c` | 11511 | ok |
| `W150_1.0` | `a71fbda4a575f393` | 2182 | ok |
| `W150_2.0` | `afc66cbb5b22aeda` | 432 | ok |

`D150_2.0` is the `robin_form = surface` + `surface_T = diagnostic` parse-guard
leg (MMWSpalling.H:626) — the abort documented in §(c′) as the reason no
surface-form (w) control can exist. Its `run.log` is the witness.

## (e) (e′) (f) (h) (j) — the Meier leg

| item | value |
|---|---|
| output | `output_fh/` (`FH*`), argv recorded verbatim in `output_fh/cmd.txt` |
| event log | `output_fh/FH_removal_events.csv`, 67 031 rows, t = 0.32 … 177.616 s |
| thermo | `output_fh/FH/thermo.dat`, 1851 rows |
| ⚠ status | **ABORTED at t = 177.616 s** on a nozzle collision at r = 26.4 mm |
| scoring window | 60–160 s — **complete**, so (e), (e′), (h), (j) are valid |
| ⚠ (f) | **NOT SCORABLE**: the registered 150–300 mm section was never drilled |
| comparator | `../p0b_one_surface_T/output_fh_rev7/`, truncated to t ≤ 177.616 s |

⚠ The rev7 comparator was produced by the **rev7** binary
(`1c22018e9ad8e8bc93bf948b1324df3d5839600a3ec9e9b96eaebf64536408cf`), not this
one. That is intended — it is the baseline — but it means the rev7 column of
every table in `RESULTS.md` is a cross-binary comparison, not an identity test.

## (a) identity half — the tree gate

| item | value |
|---|---|
| witness | `../tree_sweep_0922/witness/p1.md5`, 2639 files, written 2026-10-01 11:58:16 |
| vs `p0b7.md5` | **0 differing** |
| vs `d2z1.md5` | **0 differing** |
| guard | "binary unchanged across the sweep; p1 identity result is valid" |
| invocation | `guarded_post_all.sh p1 --force` (the cache is keyed on the binary alone) |

## Reduction scripts

The numbers in `RESULTS.md` were produced by throwaway reductions over the files
above (radius = `hypot((i+0.5)·dx, (j+0.5)·dx)` with dx = dz = 2 mm; depth =
`(299 − k_top + n_voided)·dz`; cap grip = `rc_flag` bit 16; volume = Σ`n_voided`
× 8e-9 m³ × 4 for the quarter domain). Two independent estimators agree:
the depth-field sum and the Σ`n_voided` sum both give 1.3200 L for rev7 at
177.6 s, and the window volume rate reproduces P0b's tabulated 1.794e-6 m³/s
exactly as the per-quadrant value (7.173e-6 / 4).

## Deliverable 5 — the exposed-face census

| item | value |
|---|---|
| script | `census.py` (a reduction; it starts no run) |
| input | `output_fh/FH_removal_events.csv` |
| input sha256 | `26d8996f26ae2c90c57989a426cb58d763edd040596bb6b91d149b035f0209da` |
| cut | t ≤ 177.616 s, r < 40 mm, frozen = last fired before t = 120 s |

⚠ That cut is the **unique** one returning both of the planner census's headline
rows (152/316 frozen and 78 % with ≥ 2 exposed sides). At that same cut the
live-column figure is **51 %**, not the 69 % stated; see `RESULTS.md` §5.
