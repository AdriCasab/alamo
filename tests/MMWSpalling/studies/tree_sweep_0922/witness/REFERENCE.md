# Tree-sweep reference list

The sweep has **no `run.py` and no `--check`/`--repin`**. It runs the tests
against the binary at its hardcoded path (`bin/mmwspalling-3d-g++`) and writes
an md5 list; a *re-pin* is a **declaration** made here about which list is the
reference, not a script flag.

## Current reference

**P2-29 (2026-11-10): the reference is UNCHANGED -- `face_p2_28.md5`.** `face_gate.py cmp p2_28 p2_29a` (stage 1, the
zhang phase law + `Table.H`'s sharp band; `43d87dae…`) and `cmp p2_29a p2_29` (stage 2, the separate heat kernels
deleted; `ccd15b91…`, `studies/p2_29_one_heat_path/BINARY_2.txt`) both read IDENTICAL (2,010 / 2,010 files; the
17 differing are thermo wall-clock timers). `face_p2_29a.md5` and `face_p2_29.md5` are kept beside it.

**THE gate is the face gate: `face_p2_28.md5`** (2,010 files, 19 legs) since P2-28 (2026-11-09; binary `sha256
ee4b6b9b0ed947e2…`, `studies/p2_28_column_path/BINARY_2.txt`). The next packet runs `face_gate.py run <label>` and
`face_gate.py cmp p2_28 <label>`, plus `face_gate.py rankcheck <label>`. **The tree gate is a RECORD**: P2-28 deleted the
column removal path (`MMWSpalling/Removal.H`) that most tree legs drove; `hash.sh`, `post_all_par.sh`,
`guarded_post_all_par.sh` and `p2_9_face_removal/tree_compare.py` carry a record header and are not run. **`p2_27.md5`
(2,901 paths) is the column era's last list** (`p2_28a.md5`, stage 1, is byte-identical to it); it and every earlier list
stay on disk.

**P2-28 re-pin, two builds** (`studies/p2_28_column_path/`, `CHAIN.log`):
- **p2_28a** (stage 1, additive: the prescribed patch and the evaluate-only mode on the face list, the half-cell
  Dirichlet factor corrected at both sites; `e291bc04…`, `BINARY_1.txt`): tree vs `p2_27` **0 of 2,901 changed**, 0
  only-in (controls d2aa 1,230 / d2t 1,346); face `cmp p2_27 p2_28a`: the 16 legs IDENTICAL (thermo timers only), **481
  only-in-new = KV_S 295 + HU_S 157 + GB_S 29** (the three new legs: Kant evaluate-only, Hu prescribed + buckling,
  gradient_blind on the list).
- **p2_28** (stage 2, the column path deleted; `ee4b6b9b…`): face vs `p2_28a` read by
  `studies/p2_28_column_path/audit_components.py face p2_28a p2_28` (**PASS**, `K2b_audit.log`): 2,010 / 2,010 paths, none
  on one side only; 280 plotfiles = **140 removed-ok** (each lost exactly the column fields it carried, ⊆ {`Sp_field`,
  `Sp_cluster_id`, `h_spall_field`, `crack_event_count`, `T_surf_field`, `ell_eff_field`}; every kept component
  bit-identical FAB by FAB, Header and Cell_H min/max of the kept identical) + **140 identical**; 17 `thermo.dat` lose
  exactly `recon_*` (8), `surf_cols_robin`, `surf_q_mean`, `patch_P_robin` where present, `ledger_err` renamed
  `ledger_err_raw` with identical values, every other column identical (timers excluded); every CSV identical. The
  planted mutation (`--mutate`, `Sp_field` left out of the declaration) FAILS on 140 plotfiles. Rankcheck rel 1e-15.
  ⚠ The packet declared FEET_S's h_col CSV as the one path that leaves; it never existed on the face list (only the
  column path wrote it), so "nothing on one side only" is the bar met.
⇒ **the reference becomes `face_p2_28.md5`**; `face_p2_27.md5`, `face_p2_28a.md5`, `p2_27.md5`, `p2_28a.md5` are kept.

*(Superseded by P2-28:)* **Tree gate: `p2_27.md5`** (2,901 paths) and **face gate: `face_p2_27.md5`** (1,529 files) since P2-27 (2026-11-09,
declared in the P2-27 block below; binary `sha256 27433870a7c28591…`, `studies/p2_27_cleanup/BINARY_5.txt`). The
controls against `p2_27.md5` read **d2aa 1,230 / d2t 1,346 changed, 270 only in new, 8 only in ref** (P2-27 stage 3
removed three components from every cell plotfile, so the old 26 / 863 no longer apply).

**P2-27 re-pin, five builds** (`studies/p2_27_cleanup/`, `CHAIN.log`; every stage under its own bar, every bar PASS):
- **p2_27a** (inputs: numeric retired-key compare, CheckInterlocks, the dt check, beam.schedule.file; `ccfa2baa…`):
  tree vs `p2_26` **0 of 2,909 changed**, 0 only-in-ref, **270 only-in-new = `unit/jet_d2c/output`** entering the oracle
  (declared `MOVERS_1_jet_d2c.txt`; its retired legs' stale outputs moved out first); face `cmp p2_24` **IDENTICAL**
  (M0S now replays its record as written: `face_gate.py`'s DROP lost the two C2 tokens).
- **p2_27b** (C1a, the damage lineage; `7dc8bde5…`): tree vs `p2_27a` **0 changed, 8 only-in-ref** (`unit/spall_event`,
  retired; `hash.sh` dropped its two directories), 0 only-in-new; face `cmp p2_24` **IDENTICAL**.
- **p2_27c** (C1b, D / kappa_eff / k_eff deregistered; `eed9343c…`): tree vs `p2_27b` **1,314 of 2,901 changed = every
  Cell_D of the 708 CELL plotfiles**; `studies/p2_27_cleanup/audit_components.py` read each component by component:
  every one lost exactly {D, kappa_eff, k_eff}, every kept component bit-identical, Header and Cell_H min/max of the
  kept identical; the 708 node plotfiles, every CSV and thermo identical (old outputs archived first:
  `tree_archive_p2_27b/`). Face: 84 cell plotfiles (336 Cell_D, 168 Header/Cell_H) the same; rankcheck rel 1e-15.
- **p2_27d** (the six grimsel-inert scalars refused; rev 2 `5b575c85…`, rev 1 `be7bd689…` kept): tree **0 changed**;
  face vs `p2_27c` **exactly the 70 declared paths** (`MOVERS_4_face_t0.txt`: the t = 0 Cell_D + Cell_H of the 14
  grimsel legs), fields ⊆ {kappa, beta, E, mu}_phase (the law at the start temperature instead of the typed scalars;
  `audit_fields.py`); every later plotfile identical. (Its face run also ran the stage-5 FEET legs, already listed:
  FEET_S aborted on the then-present regime-4 refusal, FEET_0 ran.)
- **p2_27** (G, the pad clearance on the face list; `27433870…`, ships): tree vs `p2_27d` **0 changed**; face vs
  `p2_27d`: the 14 old legs and FEET_0 unchanged (timers only), **156 new files = FEET_S** (14 regime-4 rows; FEET_0 0);
  rankcheck rel 1e-15.
⇒ **the references become `p2_27.md5` and `face_p2_27.md5`**; every earlier list is kept.

*(Superseded as the current reference by P2-27:)* **Tree gate: `p2_26.md5`** since P2-26 (2026-11-08, declared below); before it **`p2_18.md5`** (P2-18, 2026-11-05),
byte-identical to `p2_16.md5`, which was byte-identical to every list back to `p2_2a_r7.md5` (P2-2a was the last re-pin
before P2-26 that moved a hash: 26 `jet_enthalpy` strip files). The two fixed controls that prove the comparator is
live: `d2aa.md5` differs from `p2_18.md5` in exactly **26** files and `d2t.md5` in **863**; against `p2_26.md5` they read
**72** and **863** (P2-26's 46 movers overlap d2aa's set nowhere and d2t's everywhere). The D2z-1 declaration below is
history; its re-pin procedure ("How to re-pin") still applies.

**P2-26 re-pin of the tree list: `p2_26.md5`** (2026-11-08, binary `sha256 e0551036cce9a4cc…`,
`studies/p2_26_closed_switches/BINARY_2.txt`; the closed-switch sweep, census B1–B3 + C2 + C3). Two builds:
- **Stage 1** (B1 parse order, B2 text, B3 the pad shield / above-plane pair / `rec_fixed` / flux length and fraction /
  flux stressed volume; binary `980de09e…`, `BINARY_1.txt`): tree `p2_26a.md5` vs `p2_18.md5` **0 of 2639** (controls 26 /
  863); face `cmp p2_24 p2_26a` **IDENTICAL**. `p2_26a.md5` is kept.
- **Stage 2** (C2 `spall.surface_normal` + `spall.flake_coherence_length`, C3 `spall.melt_aware` + `melt.lambda_crit`):
  `p2_26.md5` vs `p2_18.md5` **46 of 2639 differ, 0 on one side only**, every one inside the 233 paths pre-registered
  before any edit (`MOVERS_PRE.txt`: `unit/beam_void_closure` 169 + `dev2d` 64), classified path by path against the
  archived old bytes (`pre2_out/`, hashing to `p2_18`) in `GATE_b1.txt` and field by field in `B1_FIELDS.txt`:
  - **beam_void_closure (39):** `thermo.dat` loses exactly the column `sn_recession_superseded` (C2's conditional
    registration) and every shared column is identical on all 201 rows; 38 cell-plotfile fabs (steps 2000–20000, fabs
    2 and 3 of 4) differ.
  - **dev2d (7):** the cell-plotfile fab 0 of steps 8000–56000. dev2d writes no `thermo.dat` and no event log.
  - **The moved FIELDS** are only the column scan's thickness reports: `h_spall_field` (up to 0.34 mm of a 1 mm cell:
    the scan ran in the surface-normal frame, z = s·dz·cos θ; now vertical) and, in one cell of beam_void_closure's step
    5000 plotfile, `spall_thickness` / `RoP_spall` (h_col/cos θ → h_col, 0.8 %). **`removed`, `phi`, `Temp`, `H` are
    identical in every plotfile**: the removal is whole-cell (one dz per firing) whatever the scan reports. Node plotfiles
    and the t = 0 plotfiles do not move.
  - **C3 moved nothing**, as pre-registered: no path of the six melt-aware gate tests differs (their hottest cell is
    655–659 K against the 1448 K solidus, Λ_L 0 in every cell of every plotfile; `Z0.txt` §4).
  Face gate `cmp p2_24 p2_26`: 1,217 files each side, **IDENTICAL** (14 `thermo.dat` differ by the wall-clock timers
  only); `rankcheck p2_26` decomposition-stable (rel 1e-15). ⚠ `face_gate.py`'s `DROP` gained
  `spall.flake_coherence_length` / `spall.surface_normal`: M0S replays P2-9's M0 argv, which writes them at `0.0` / `0`,
  and the retired-at-default table compares STRINGS (`"0.0" != "0"`), so the first stage-2 run aborted M0S (archived:
  `studies/p2_26_closed_switches/face_2_attempt1/`). Dropping a default-valued token is inert (stage-2 cmp IDENTICAL).
⇒ **the tree reference becomes `p2_26.md5`** (`p2_18.md5` and `p2_25.md5` kept); the face reference stays
`face_p2_24.md5`. The next packet runs `tree_compare.py <label> p2_26 d2aa d2t` (controls 72 / 863) and
`face_gate.py cmp p2_24 <label>`.

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

**P2-23 (2026-11-06, binary `sha256 42a115a0004ea460…`, `studies/p2_23_beam_port/BINARY.txt`): tree NOT re-pinned;
the face gate GAINS a leg.** Tree gate `p2_23.md5` vs `p2_18.md5`: **0 of 2639** (controls d2aa 26, d2t 863): with the
surface patch on, the new face-list switch equals the old one, and no gate deck sets `spall.buckle_log_csv`. Face gate
`cmp p2_22 p2_23`: the 13 existing legs **identical** (0 substantive differences, thermo timers only), 0 paths only in
p2_22, and 45 paths only in p2_23, all under `BEAM_S/` (the new flameless leg, `unit/face_beam/input` to 1 s;
`studies/p2_23_beam_port/A1.txt`). ⇒ the tree reference stays `p2_18.md5`; **the face reference becomes
`face_p2_23.md5`** (the 13 legs byte-identical to `face_p2_19.md5`'s substance plus `BEAM_S`). The next packet runs
`face_gate.py cmp p2_23 <label>`.

**P2-25 (2026-11-07, binary `sha256 c052c348876d515b…` = the P2-24 addendum build with `grimsel.lambda.cap_theta`,
`studies/p2_25_regime_map/BINARY.txt`): NOT re-pinned.** Tree gate `p2_25.md5` vs `p2_18.md5`: **0 of 2639** (controls
d2aa 26, d2t 863); face gate `cmp p2_24 p2_25`: **IDENTICAL** (timers only), so the addendum key's default (−1) is inert
on both lists. ⚠ `face_gate.py cmp` now compares thermo columns with NaN EQUAL to NaN (`col_moved`, P2-25 DR-1): under
the old list `!=`, two parses of the same file disagreed on every NaN-holding column (BEAM_S: 11), so a1 could not read
IDENTICAL once a timer moved the file's md5 (`studies/p2_25_regime_map/Z1.txt`). `rankcheck` keeps the old test
(unchanged by the packet). ⚠ This sweep took 57 min instead of ~6: three 4-rank lanes plus Spotlight indexing on 5
performance cores tipped MPI's busy-waiting ranks into stalls (the long legs 5–128× slower, the short ones normal); the
hashes are unaffected (GUARD: binary unchanged). The references stay `p2_18.md5` and `face_p2_24.md5`.

**P2-24 (2026-11-06, binary `sha256 98a261d39bf217a3…`, `studies/p2_24_melt_vapour/BINARY.txt`): tree NOT re-pinned;
the face reference RE-PINNED to `face_p2_24.md5`.** Tree gate `p2_24.md5` vs `p2_18.md5`: **0 of 2639** (controls d2aa
26, d2t 863): the melt set now in the rock file changes the phase table only above 1410.5 K, and no gate deck gets there
(the hottest face-gate leg is CAL_S at 818 K; `studies/p2_24_melt_vapour/Z0.txt` §2). Face gate `cmp p2_23 p2_24`:
1,217 files each side, 0 on one side only, **32 differ**, every one in a declared class (`studies/p2_24_melt_vapour/A1.txt`,
scored by `a1.py`):
- **(i) the 13 flame legs' `thermo.dat`:** two appended columns `fr_n_molten`, `fr_n_lid` (registered on every face-list
  deck), 0 on every row; every other non-timer column identical on every row. Nothing else of theirs differs.
- **(ii) BEAM_S (`unit/face_beam/input`, which includes the #19 deck and so now carries the beam's melt keys and
  `vapor.enabled = 1`):** its `thermo.dat` keeps all 91 common columns identical and appends 14 (`beam_P_beam`, the
  `bc_*` melt record, the four counts, `ledger_E_vapor`; the counts and the vapour book 0); its three plotfiles' `Header`
  and `Cell_D_*` differ because they gain the vapour fields `RoP_vap`, `vapor_event`, `vapor_thickness`, while all 43
  common fields are identical in double precision. Its buckle log is byte-identical: the #19 deck turned the log off
  (P2-24 Goal 7) and the fixture turns it back on, because BEAM_S is the only gate leg that writes one (the first P2-24
  chain stopped at a1 on exactly that missing file; `studies/p2_24_melt_vapour/run1_a1stop/`).
⇒ the tree reference stays `p2_18.md5`; **the face reference becomes `face_p2_24.md5`**. The next packet runs
`face_gate.py cmp p2_24 <label>`; `face_p2_23.md5` is kept.

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

**P2-30 (2026-11-10): reference unchanged (`face_p2_28.md5`).** The header's split into slices, the comment moves and the
two code changes read IDENTICAL on `face_p2_30a`, `face_p2_30b` and `face_p2_30` (2,010 files; 17 differ in timer columns only).
⚠ `_buckle_log.csv` is HEADER-ONLY in both legs that write it (BEAM_S, HU_S end before any line fires), so this gate does not
reach the buckle ROW writer; and no leg writes a vapour row. A change to either writer needs `studies/p2_30_structure/covers.py`.

## How to re-pin next time

⚠ Since P2-28 the gate is the FACE gate: `face_gate.py run <label>` → `witness/face_<label>.md5`, then
`face_gate.py cmp p2_28 <label>` (and a component audit when components leave, as P2-27 / P2-28 did). The tree steps
below are the record of how the column era re-pinned.

1. `bash witness/guarded_post_all.sh <label>` → `witness/<label>.md5` and
   `witness/<label>_wit/`.
2. Diff against the current reference above and classify every differing path.
3. Update this file to name the new reference and what it audited.
4. Leave the old list in place.

**2026-11-10, after P2-30: THE REFERENCE BECOMES `face_p2_30c.md5` (21 legs, 3,051 files; binary `sha256 0db6623783c9961b…`).**
Two legs added, because P2-30's step zero found that NO leg wrote a vapour row or a buckle-log data row:
`VAP_S` (`unit/face_vapor/input` as it stands: 64 regime-2 rows in the event log) and `BKL_S` (`unit/face_buckle`'s F_S
+ `spall.buckle_log_csv=1`, 8 s: 3,456 firing rows and 3,456 buckle-log rows). Continuity: `face_gate.py cmp p2_28
p2_30c` reads **2,010 / 2,010 common files, 0 substantive differences (17 thermo files, timers only), 1,041 files only
in `p2_30c`** (the two new legs) — its VERDICT line says DIFFERENT only because of those extra paths. The new legs'
three logs are byte-identical to P2-30's covers (`studies/p2_30_structure/covers/post4a`). The binary is P2-30's source
after the comment review of the same day (comments only: `studies/p2_30_structure/comment_review_check.py`, the
pre-image `src_pre_5`). **From here: `face_gate.py run <label>` then `face_gate.py cmp p2_30c <label>`.**
`face_p2_28.md5` is kept as the 19-leg record.


**P2-31 (2026-11-10): THE REFERENCE BECOMES `face_p2_31.md5` (21 legs, 3,051 files; binary `sha256 e705b8591d083230…`,
`studies/p2_31_fixes_mesh/BINARY_p2_31.txt`).** Four builds, each gated (`studies/p2_31_fixes_mesh/RESULTS.txt`):
- **p2_31a** (stage 1: the reconstruction reads the start cell's initial temperature, two keys parsed above the spall-off
  return, three new aborts, the robin_form message, `kappa()` split into `diffusivity()` / `conductivity()`;
  `3f29e987…`): vs `p2_30c` **0 substantive differences on all 21 legs** (`GB_S`, the one leg with T_ref != IC, included).
- **p2_31b0** (stage 2a: the dead list out, the vapour keys retired at their defaults, `microstructure_enabled` and seven
  unreachable aborts gone; `ba1f6dd5…`): vs `p2_31a` and vs `p2_30c` **0 substantive differences**.
- **p2_31b** (stage 2b: `regime_field` removed -- it duplicated `regime` at every write site; `296fadf3…`): vs `p2_31b0`
  the plotfiles of the 19 removal legs lost exactly `regime_field`, every other component FAB by FAB, every thermo column
  and every CSV identical (`components.py`; KV_S / HU_S never carried it).
- **p2_31** (stage 3: the overlap rule -- the phi shift reads the +z line's own h, the plotted thickness / rate the largest
  h through the cell; `e705b859…`): vs `p2_31b` **0 substantive differences on all 21 legs**. The rule changes nothing
  the gate observes: the gate's overlaps (`fr/TL45_1` 792 cells, `SPH_1` 10, `BORE_1` 1) are interior cells, where the
  old shared slot fed nothing; the first observable difference is on `TD` at 1 mm from 3.304 s (`unit/face_removal
  --ovl`, `studies/p2_31_fixes_mesh/OVL_SEARCH.txt`).
The intermediate directories `face_p2_31a/`, `face_p2_31b0/`, `face_p2_31b/` were deleted for disk (16 GiB free at the
start); their md5 lists are kept and each is regenerable with `face_gate.py run <label>` on the binary kept in
`studies/p2_31_fixes_mesh/bins/` (`pre_2` = p2_31a, `pre_2b` = p2_31b0, `pre_3` = p2_31b; local only, `.gitignore`d).
**From here: `face_gate.py run <label>` then `face_gate.py cmp p2_31 <label>`.**
