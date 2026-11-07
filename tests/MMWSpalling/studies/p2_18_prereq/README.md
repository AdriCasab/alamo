# P2-18 — prerequisites for the retire packets (census A1–A3), 2026-11-05

Done by the advisor at the user's instruction, without an `/implement` packet (P2-12's precedent). The
user's three asks: short witness runs; "remove from the input files inputs which are no longer in
use"; the shared helpers in a neutral file.

## What changed

| | |
|---|---|
| **A1 face-list gate** | `witness/face_gate.py` (13 short legs, every one on the face list; 1,172 hashed files; ~2.5 min). Baseline `face_p2_18_pre.md5` hashed on the pre-change binary BEFORE any source edit. Declared in `REFERENCE.md`; rows in `WHICH_TESTS.md`. |
| **A3 helper move** | `src/Integrator/MMWSpalling/Shared.H` (new, 314 lines, class-body fragment included before Removal.H): ElasticConstants, RefillTempGhosts, ScanTAtZv, ScanKIcPa, ScanSigmaT, DepthScanH, RemapProfile moved verbatim from the top of Removal.H; plus `RemovalStepBegins(lev)` (the guards) and `ResetRemovalStepFields(lev)` (the per-step reset), which Removal.H now calls at the same two points. The column path's PASS 1 (its own copy of the top-solid scan) replaced by the existing `ColumnTopSolid(lev)`. Removal.H 2,580 → 2,282 lines. |
| **A2 retired at default** | `MMWSpalling::Parse`: a `retired_at_default` table after the P2-16 abort list. Each entry is read once (marks a strict-scope key used) and compared with the value the deleted code accepted; equal or `*` → run; else abort naming key, value, note. Entries today: `surface_patch.face` = zhi, `surface_patch.T_ext` = `*`, `surface_patch.type` = prescribed_temperature. |
| **A2 dead parses deleted** | `microstructure.phase<i>.name`, `.alpha_attenuation` (REQUIRED, never read) + their members; `surface_patch.type`/`face`/`T_ext` + the `T_ext` member; `microstructure.phase<i>.grimsel.T_ref` (overwritten before use). |
| **A2 deck cleanup** | 602 lines in 91 kept decks (`A2_DECK_EDITS.txt`): every DEAD-key line, and every line that set a RETIRING key to its parse default (`damage.enabled 0`, `solver explicit`, `kappa_damage_alpha 0.0`, `spall.sample top_cell`, `surface_normal 0`, `flake_coherence_length 0.0`, `melt_aware 0`, `bit.enabled 0`, `weibull.size_volume cell`, `flake_length twocell`, `flake_fraction cf`). Live values stay. Records (`legs/*/argv.sh`, study `run.py`, `pre/`, `retired_units/`) untouched. |

## Gates (`chain_a3.sh` on the A3 binary, `chain_a2.sh` on the final binary) — ALL PASS

| gate | A3 binary | final binary |
|---|---|---|
| tree gate vs `p2_16.md5` (controls d2aa / d2t) | **0 of 2639** (26 / 863) | **0 of 2639** (26 / 863) |
| face gate vs `face_p2_18_pre.md5` | **0 of 1172** substantive (13 thermo files differ in the wall-clock timer columns only) | **0 of 1172** |
| `buckle_lock`, `retired_keys`, `lateral_flux_lock` | rc 0, 0 FAIL | rc 0, 0 FAIL |
| `face_removal --check` (a3 i: one copy of each helper) | 9/9 PASS | 9/9 PASS |
| retired-at-default legs (K1–K6, `grimsel_calibration`, stop 1e-9) | — | face=zhi RUNS; face=xlo ABORTS naming the key and "retired at its default (zhi)"; T_ext=293.15 RUNS; type=foo ABORTS likewise; phase0.alpha_attenuation=2000 RUNS (top-level, silently inert); grimsel.T_ref=293.15 ABORTS as "specified but not used" (strict scope; no deck sets it) |
| every kept deck starts (93 decks, stop 1e-9) vs census RUNNABLE_ALL | — | 91 unchanged; `hoeser/input_lab` is a template (needs argv, as before); 2 Hu `diag/` decks aborted once — see findings |
| rank check (−np 1 vs −np 4, CAV + T00_1_k1) | — | counts and CSVs identical, thermo columns agree to 1e-15 |

## Findings

1. **Plotfile bytes are rank-layout dependent**: a byte rank check is impossible by construction; the
   face gate is pinned at −np 4 like the tree gate, and `rankcheck` compares the physics.
2. **`hu_spall_onset/diag/input_sandstone2_clampall` crashes intermittently (bus error) at
   stop_time 1e-9**: 1 of 4 runs on the final binary, 0 of 2 on the pre-change binary with the
   original deck — not reproducible, not caused by the edits (its cleaned lines were name,
   alpha_attenuation, kappa_damage_alpha 0.0, type, face, T_ext, solver). The `diag/` records
   retire with Hu's end-to-end case (decision 5). `input_granite2_base` aborted once in the chain
   and ran cleanly on re-run.
3. `face_removal --check` still reports P2-9's four tilt-rate FAILs (TL20/TL35/TL45/TR, 4–8 %):
   P2-10 recorded them as the P2-9 record's per-firing deficit; unchanged here.
4. `unit/jet_d2c` runs in the tree sweep but is not hashed; left run-only because its legs use
   `rec_fixed`, which census B3 retires.

Files: `BINARY_pre.txt` / `BINARY_a3.txt` / `BINARY_a2.txt`, `bin_pre_mmwspalling-3d-g++`, `build_*.log`,
`tree_*.log`, `GATE_tree_*.txt`, `face_*.log`, `GATE_face_*.txt`, `lock_*_*.log`, `K_*.log`, `RUNNABLE_a2.txt`,
`face_a2_rankcheck.log`, `a2_apply.py` + `a2_dry.txt` + `A2_DECK_EDITS.txt`, `TEST_EDITS.txt`, `runnable.py`.
