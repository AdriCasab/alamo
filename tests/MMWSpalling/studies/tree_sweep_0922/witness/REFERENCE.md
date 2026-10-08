# Tree-sweep reference list

The sweep has **no `run.py` and no `--check`/`--repin`**. It runs the tests
against the binary at its hardcoded path (`bin/mmwspalling-3d-g++`) and writes
an md5 list; a *re-pin* is a **declaration** made here about which list is the
reference, not a script flag.

## Current reference

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
