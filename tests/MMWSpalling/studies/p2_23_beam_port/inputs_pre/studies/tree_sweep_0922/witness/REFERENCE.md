# Tree-sweep reference list

The sweep has **no `run.py` and no `--check`/`--repin`**. It runs the tests
against the binary at its hardcoded path (`bin/mmwspalling-3d-g++`) and writes
an md5 list; a *re-pin* is a **declaration** made here about which list is the
reference, not a script flag.

## Current reference

**Tree gate: `p2_18.md5`** (P2-18, 2026-11-05), byte-identical to `p2_16.md5`, which was
byte-identical to every list back to `p2_2a_r7.md5` (P2-2a was the last re-pin that moved a
hash: 26 `jet_enthalpy` strip files). The two fixed controls that prove the comparator is live:
`d2aa.md5` differs in exactly **26** files, `d2t.md5` in **863**. The D2z-1 declaration below is
history; its re-pin procedure ("How to re-pin") still applies.

**Face-list gate: `face_p2_19.md5`** since P2-19 (declared below); before it **`face_p2_18.md5`** (P2-18 A1, 1,172 files), produced by
`face_gate.py run p2_18` and byte-identical (but the three wall-clock timer columns) to
`face_p2_18_pre.md5`, the same legs on the pre-P2-18 binary. ⚠ The tree gate never runs the
face list (no gate deck sets `surface_patch.face_list = 1`; census 2026-11-05 finding 1), so a
packet that touches removal passes **both** lists: `guarded_post_all_par.sh <label>` against
`p2_18.md5` and `face_gate.py run <label>` + `face_gate.py cmp p2_18 <label>`. The face legs are
short on purpose (13 legs, ~2.5 min at -np 4): five `unit/face_removal` fixtures to 8 s, the two
`unit/face_list` builders, `unit/face_buckle` F_S to 8 s, two `unit/face_heat` path legs, and
three production decks truncated by `stop_time` (P2-9's Höser leg M0 minus the P2-10 tokens to
3 s, P2-16's HOS_R to 2 s, `kant/grimsel_calibration` to 10 s), every one past its first firing.
Hashes are **-np 4 hashes**: plotfile bytes depend on the box-to-rank layout, so
`face_gate.py rankcheck <label>` compares the physics (thermo columns numerically, CSVs byte
for byte) at -np 1 vs -np 4; on P2-18 the legs agree to round-off (`face_p2_18_RANKCHECK.txt`,
ledger closure columns differ at 1e-15). To re-pin the face list: run, `cmp`, classify every
differing path, declare it here, keep the old list.

Both lists were taken on binary `sha256 5561de02089a0e97…` recorded in
`studies/p2_18_prereq/BINARY_a2.txt`; the pre-change binary is kept at
`studies/p2_18_prereq/bin_pre_mmwspalling-3d-g++`.

**P2-19 re-pin of the face list: `face_p2_19.md5`** (2026-11-05, binary `sha256 c34e428b8eba6b9a…`,
`studies/p2_19_flake_depth/BINARY.txt`; tree gate `p2_19.md5` 0 of 2639 vs `p2_18.md5`, so the tree reference does
not move). `face_gate.py cmp p2_18 p2_19`: 1,172 files each side, **26 differ, 0 only on one side, and no `Cell_D`,
`Cell_H`, `Header` or `Checkpoint` differs** (the audit, `studies/p2_19_flake_depth/A1.txt`). Every differing path is
one of two declared classes:
- **(i) the 13 `thermo.dat`:** the header changes only by `fr_noskin_{floor,wall,roof}` → `fr_nodepth_{floor,wall,roof}`
  (the P2-19 refusal replaces the no-skin refusal, same site and position), and every other non-timer column is
  identical on every row. The renamed columns' values differ where the two refusals differ: `CAL_S` 0 → 1,111 (wall),
  `fr/BORE_1` 862,051 → 928,106 (+4 roof), `fr/SPH_1` 1,067,053 → 2,147,217 (floor); 0 → 0 on the rest.
- **(ii) the 13 `*_removal_events.csv`:** the header gains `,ell_p,z_peak,ell_pk,h_E` at the end, and the pre-existing
  columns of every row are byte-identical (row counts equal on every leg).
The fallback branch fires on no gate leg, so the removal and the field are untouched. The old list `face_p2_18.md5`
is kept. **From P2-19 on, a face-touching packet compares against `face_p2_19.md5`.**

**P2-21 (2026-11-06, binary `sha256 650d31436e506203…`, `studies/p2_21_loss_table/BINARY.txt`): NOT re-pinned.**
Tree gate `p2_21.md5` vs `p2_18.md5`: **0 of 2639** (controls d2aa 26, d2t 863): no gate deck carries a varying
`beam.tan_delta_T`, so the column optical-depth sum is reached by none and no header gains its nine columns. Face gate
`face_gate.py cmp p2_20 p2_21`: **IDENTICAL** (timers only) -- HOS_S and CAL_S now read the 17-pair table through
`rock_grimsel.inputs`' nested `FILE =` but deposit nothing (`beam.can_deposit()` false), so their headers stand.
Also byte-identical: `unit/beam` against `p2_20_beam_deposition/unit_beam_output_pre/` (28 files, `unit/beam_tau` U3)
and the Zhang deck with a one-pair table against P2-20's `legs/Z1/run` (112 files + thermo.dat, a3). The references
stay `p2_18.md5` and `face_p2_19.md5`.

**P2-22 (2026-11-06, binary `sha256 3604d38df733579e…`, `studies/p2_22_inrock_spread/BINARY.txt`): NOT re-pinned.**
The in-rock beam spread (the waist law now uses the wavelength inside the rock, so it depends on `beam.eps_r`).
Tree gate `p2_22.md5` vs `p2_18.md5`: **0 of 2639** (controls d2aa 26, d2t 863). Face gate
`face_gate.py cmp p2_21 p2_22`: **IDENTICAL** (timers only), and `face_p2_19.md5` -- the declared reference -- is
still identical to `face_p2_22.md5` (13 differing files, every one a `thermo.dat`, timers only). The references
stay `p2_18.md5` and `face_p2_19.md5`.

⚠⚠ **BYTE IDENTITY ON AN eps_r != 1 DECK DOES NOT CROSS THIS BOUNDARY, BY DESIGN.** Every gate fixture declares
`beam.eps_r = 1` (the in-rock wavelength then equals the free-space one), which is why both gates are clean. On a
deck with `eps_r = 5.02` the spread changes: `studies/p2_22_inrock_spread/C1_vs_Z1.txt` records the Zhang deck
against P2-20's `legs/Z1/run` as **75 of 112 plotfile files differing**, with the largest relative move
**6.4e-06** (`beam_E_dep`) and 4.3e-06 (`beam_E_trans`); `beam_closure_err` moves 83 % *relative* because it is
tiny and got better (-6.54e-07 -> -1.08e-07 on the cold first step). ⇒ **P2-20's a3 and P2-21's a3 ("the Zhang deck
byte-identical to `legs/Z1/run`") are SUPERSEDED from P2-22 on and must not be re-run against `bin/` as regression
checks** -- they would fail for this declared reason. The binary those rows were taken on is preserved at
`studies/p2_22_inrock_spread/bin_pre/mmwspalling-3d-g++` (`650d31436e506203…`) for anyone who needs the old bytes.

⚠⚠ **THE eps_r = 1 FIXTURE SET IS A GATE HOLE for anything that reads the refractive index.** Because every unit
fixture and every gate deck sets `eps_r = 1`, a P2-22 build that computed the incident power with a **0.55 m wide**
beam passed every unit test and would have passed both gates; only a smoke run on a granite deck caught it. The
`-1` sentinel for "no void path" had been fed into the spread formula. ⇒ any packet touching an index-dependent
quantity owes a fixture with `eps_r != 1`; until one exists, neither gate covers that path.

**P2-20 (2026-11-06, binary `sha256 0d071f167acfc08f…`, `studies/p2_20_beam_deposition/BINARY.txt`): NOT re-pinned.**
Tree gate `p2_20.md5` vs `p2_18.md5`: **0 of 2639** (controls d2aa 26, d2t 863), including the four beam-ON gate tests
(dev2d, beam_void_closure, scalar_flaw, s1_smoke), whose decks now carry `beam.eps_r = 1.0` and a `beam.tan_delta`
string that reproduces the old typed α bit for bit (the evaluation order 2π√ε′·(tan δ/λ) makes every gate pair exact;
`Z0.txt`). Face gate `face_gate.py cmp p2_19 p2_20`: **IDENTICAL** (timers only). The references stay `p2_18.md5` and
`face_p2_19.md5`.

## D2z-1 declaration (2026-09-25), superseded as the current reference by P2-18

**`d2z1.md5`** — pinned by **D2z-1** (2026-09-25), whole-cell removal made
unconditional and `spall.removal_increment` deleted.

Produced by `bash witness/guarded_post_all.sh d2z1`, whose guard fingerprints
the binary before and after and fails loudly if it moved mid-sweep (it did not).

Audit against the superseded list is kept at
`studies/d2z1_delete/tree_repin_audit.txt`. Summary: 2639 files both sides,
**no path appears or disappears**, 837 differing, and exactly **one** file gains
a column — `unit/beam_void_closure/output/thermo.dat` gains
`sn_recession_superseded`, the Step 20 witness registered only under
`spall.surface_normal = 1`. No other header changes anywhere.

## D2z-2 verification (2026-09-28) — reference UNCHANGED

`bash witness/guarded_post_all.sh d2z2` → `witness/d2z2.md5`. Compared against
`d2z1.md5`: **2639/2639 byte-identical, zero differing hashes, no path added or
removed.** The guard confirmed the binary never moved
(`GUARD: binary unchanged across the sweep`).

D2z-2 added one key (`surface_patch.step_diag_radius`) and eight thermo columns
(`delta_faces`, `side_P_face_multi`, `wall_step_n_in`, `wall_step_q_in`,
`wall_step_q_eff_in`, `step_g_mean`, `step_g_max`, `wall_step_f_multi`) and
changed `StepConductance`'s normaliser. **None of it reaches a run that does not
set `surface_patch.side_face_step_height`**, which no input in the 2639-file set
does — so `d2z1.md5` remains the reference and is NOT re-pinned.

⚠ `d2z2.md5` is a **verification** list, not a new reference.

## D2aa verification (2026-09-28) — reference UNCHANGED

`guarded_post_all.sh d2aa` → **2639/2639 byte-identical to `d2z1.md5`**, zero
differing hashes, no path added or removed. Guard confirmed the binary held.

D2aa added `surface_patch.pad_shield_open`, the `side_face_h` / `side_face_T`
pair and seven thermo columns (`side_A_bore/pads/wall`, `side_P_bore/pads/wall`,
`wall_step_f_multi_pads`), and changed the `SurfaceCellFlux` call at three sites.
**None of it reaches a run that does not set the new keys**, and no input in the
2639-file set does — so `d2z1.md5` remains the reference and is NOT re-pinned.

## The sweep is now PARALLEL — and it was validated by reproduction

`guarded_post_all_par.sh <label>` runs the same 11 tests in three static lanes
and calls the same `hash.sh`. Proven, not argued:

| | serial | parallel |
|---|---|---|
| wall time | 794 s (13.2 min) | **274 s (4.6 min)** |
| files | 2639 | 2639 |
| **differing hashes** | — | **0** |

`d2aa_par.md5` is byte-identical to `d2aa.md5` on the same binary. Lane totals
measured 254 / **272** / 208 s against predicted 265 / 293 / 236; lane B is the
critical path. Three lanes, not eleven: each test runs `mpirun -np 4`, and the
machine has 15 logical / 5 performance cores.

⚠ `post_all.sh` and `wit.sh` are **kept unchanged** as the serial reference. If
the parallel list ever disagrees with a serial one, the serial is right.

## Gate results are memoised on the BINARY HASH

Both guards read and write `witness/.gate_cache/<sha256>`. The gate's claim is a
property of the binary — *this binary, with every new key off, reproduces the
frozen tree* — so a second request for the same hash is a no-op. `--force`
overrides. The entry is written **only after the guard confirms the binary held**,
so a mixed-binary sweep can never seed it.
⚠ Keyed on the hash **only**, never the label: two labels against one binary is
exactly the case worth short-circuiting. Do not "improve" this into a label cache.

## Superseded

**`tree0922.md5`** — the pre-D2z-1 reference. **Preserved untouched.** Nothing
is deleted or overwritten; earlier labels (`d2x`, `d2v`, `d2t`, …) also remain.

## How to re-pin next time

1. `bash witness/guarded_post_all.sh <label>` → `witness/<label>.md5` and
   `witness/<label>_wit/`.
2. Diff against the current reference above and classify every differing path.
3. Update this file to name the new reference and what it audited.
4. Leave the old list in place.
