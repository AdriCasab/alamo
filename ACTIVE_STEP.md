# Active Step: P2-26 — The closed-switch sweep: delete every removal switch the record closed, fix two parse defects, and re-pin the gate for the three switches that gate decks still set

Status: ready for implementation (planned 2026-11-07; P2-25 ACCEPTED and archived — `ARCHIVE_DONE.md` *"P2-25"*; the user's direction 2026-11-07: the merge of the removal paths and the cut of accumulated code come next, the stress-solver packet deferred because the regime map's depth boundary is geometric). This is census packets **B1, B2, B3, C2, C3** of `docs/project/CODE_CENSUS_2026-11-05.md` §2; the DAMAGE lineage (C1), the feet port (G), the column path (H) and Zhang (D) follow in their own packets.

## Context

The production model runs on the face list; the column path and a decade of study switches are still compiled, parsed
and guarded. The 2026-11-05 census ordered their removal and the user decided every open question that day (log
`Claude_markdowns/2026-11-05.md`, "Census §4 open decisions"): Q6 `melt_aware` DELETE, Q7 `surface_normal` RETIRE
with the coherence length in one re-pin. This packet takes every deletion that is off in production and either off
in every gate deck (B1–B3: zero differing files expected) or on in a named handful of gate decks (C2, C3: a declared
re-pin whose movers are listed BEFORE the run). P2-18 already removed the DEAD keys and cleaned 602 default lines
from the decks, so the closed switches below have **no live setter** except the three re-pin keys. Two defects ride
along because they sit in the same parse block: the flake-cap check that reads `h_spall_max` before it is parsed, and
a diagnostic key parsed inside a closure branch. Carried: every deletion is one commit's worth in a chain with the two
gates after each build (SR 9: `src_pre/` for EVERY touched file, P2-25 DR-4); nothing is tuned; no physics on the
production path changes — the face gate must read IDENTICAL throughout.

## Sources (opened 2026-11-07; `M:` = `src/Integrator/MMWSpalling.H`, `R:` = `MMWSpalling/Removal.H`, `F:` = `MMWSpalling/FaceRemoval.H`, today's lines)

- **The census** `docs/project/CODE_CENSUS_2026-11-05.md` §2 (B) :184–276 and (C2) :333–364, (C3) :365–391 — the item
  lists below are its, re-resolved to today's lines here; its `M:` numbers are the 11-05 snapshot and have drifted.
- **B1 parse defects.** `spall.flake_h_max` parsed `M:6156`, its check against `value.spall_h_spall_max` `M:6159–6162`;
  `spall.h_spall_max` parsed `M:6302` (AFTER the check ⇒ the check reads the member's initialiser, not the deck).
  `surface_patch.side_power_interval` parsed `M:8805–8809` inside the `jet_closure = enthalpy` branch (a `jet_closure
  = none` deck cannot set it). Dead members: `fl_wtop` (`M:3509`, written `M:4027, 4035`, never read);
  `sp_surface_T_mode` (set `M:6105, 6111, 6117`, declared `M:9568`, never read — the two bools are).
- **B2 stranded text / dead locals:** the census's list (B2) by content, not line: the P1b comment blocks, the
  `pc`/`p_n`/`q_rob q_los` unused bindings, the `(void)` pairs, the duplicated pad-ring block (`M:2328–2427` region,
  dies with B3's pad shield). Comments and locals only — bit identity.
- **B3 closed switches (0 live setters; test scripts and study `run.py` records name them):**
  `surface_patch.pad_shield_open` (parse `M:8471–8484`, `CheckPadRing` `M:8481, 8937, 8954`, uses `M:1061, 2328,
  2389, 2427, 3890, 3905, 8415`); `jet_above_h` / `jet_T_above` (parse `M:8446–8456`, uses `M:1064, 2363–2403,
  3467, 3897, 8778–8827`); `jet_T_ent_mode = rec_fixed` + `jet_T_rec_fixed` (parse `M:8694–8716`, uses `M:2648,
  2787, 2790, 3023`); `spall.flake_length = flux` + `spall.flake_fraction = weibull` (parse `M:6176–6185`; the
  `ell_flux` branch `M:3196–3240` in `BuildSurfaceTColumns`; thermo `len_fb_h0/q0` registered `M:300–306`; members
  `M:9373–9374, 9391–9392`); `weibull.size_volume = flux` + `weibull.V_ref` (parse `M:9934–9952`,
  `WeibullSizeFactor` `M:8347`, callers `M:6842`, `R:` one, members `M:9938–9946`). Evidence each is closed: census
  C.7.3 (D2aa, D2v, D2e, P1g-3 rejected, P2-2d "no target").
- **C2** `spall.surface_normal` (parse `M:6024–6050`, registration of `sn_recession_superseded` under the key
  `M:1009–1011`, `SP_SN_COS_MIN` `M:6506, 6537, 9219`, `ColumnSurfaceCos` `M:6508, 6637`, `R:256–257, 1711, 2109`,
  `sn_recession_superseded` `R:1710, 2105, 2124`; the `sn_cth`/`sn_icth` parameters `R:931–933, 1165, 1246, 1271,
  1704`; the literals `1.0, 1.0` at `F:591`); `spall.flake_coherence_length` (parse `M:6313–6316`, `spall_flake_coh_len`
  `M:1145, 9133`, `R:1510–1511`). **Live setters:** `surface_normal = 1` in `unit/beam_void_closure/input`,
  `unit/seed_void/input`, `validation/meier/sp_meier_pilot/input_2d_dev`, `…/input_2d_dev_melt`;
  `flake_coherence_length` 0.004 in `unit/{beam_void_closure,surface_T_recon,gradient_blind,seed_void}/input`,
  `validation/meier/sp_meier_pilot/input_{energy,drilling,drilling_25mm,2d_dev,2d_dev_melt}`,
  `validation/kant/sp_kant_onset_meiercfg/input_{model,weibull}`, `validation/p0_steady_column/input`; 0.002 in five
  `extensions/hu_spall_onset_mmwbeam_lefm/input_*`. **Gate members among them:** `beam_void_closure` (169 files) and
  `dev2d` = `meier/input_2d_dev` (64 files) — the C2 movers; the thermo column `sn_recession_superseded` leaves both
  and the criterion's `/cos θ` sampling changes `Cell_D`.
- **C3** `spall.melt_aware` + `melt.lambda_crit` (parse `M:6069–6076`, members + `MeltSkeleton` `M:9196–9211`, refusal
  `M:1127`, `R:723, 729, 1071–1079` (the Λ_L gather and `LL_at_didx`), `R:1174–1179, 1203, 1229` (`f_skel`)). **Live
  setters (all `= 1`, `lambda_crit 0.4`):** `unit/{robin_pinned/input_pinned, robin_jet/input_jet,
  jet_enthalpy/input_quarter, jet_enthalpy/input_strip, robin_feet/input_feet, robin_feet/input_column,
  scalar_flaw/input}`, `studies/s1_surface_resolution/input`, `validation/meier/sp_meier_pilot/input_2d_dev_melt`,
  plus `unit/sp_melt_skeleton/input`. Every one of these decks sets `phase0.T_m 1473.0`, `fus_width 50.0` (solidus
  1448 K). **Zero-compute reading (this session, the p2_25 sweep's outputs in `unit/*/output`):** the maximum `T_top`
  over every regime-1 event is **655 K** on robin_pinned, robin_jet, jet_enthalpy, robin_feet and scalar_flaw — no
  firing ever saw Λ_L > 0, so `f_skel ≡ 1` there and the deletion is bit-neutral on those 1,458 files. The **s1 smoke
  legs** (`s1_surface_resolution/output_smoke`: R-face, R-cell, F1s, F1, M-2000, M-100, R-face-1400; 142 files;
  15 MW/m² flame, no events CSV written) are the only candidates for a C3 move — `Lambda_L` is a plotted field there.
- **The retired-at-default table** `M:157–176` (`RetiredAtDefault {key, value, note}`; `"*"` = any value; the
  comparison is a STRING match on the joined tokens, so a numeric default must be written as decks write it).
  **ALAMO ignores an unknown top-level key silently** (SR 1): `spall.*`, `weibull.*` (`pp.query_default` at `M:9934`)
  and `melt.*` are top-level and NEED a table entry; `surface_patch.*` is a strict scope and aborts on an unread key.
- **Locks and tests that name the symbols:** `unit/lateral_flux_lock/test:308` (the depth-scan stress regex
  `f_skel \* \(sp_relief \* phase_E_loc`), `:405–436` (`WeibullSizeFactor`, `size_volume`, `V_ref` clauses);
  `unit/retired_keys/test` (clauses R1–R6, MUTATIONS, K legs on `unit/gradient_blind/input`);
  `unit/removal_increment/test:14, 71, 168–180` (the SN leg); `unit/face_heat/test:49–50` (BASE tokens
  `spall.flake_coherence_length=0.0`, `spall.surface_normal=0`); `unit/jet_d2c/test:139, 200–205` (RECF legs KT,
  K15, A17–A19, run but not hashed); `unit/sp_melt_skeleton/test` (the A/B test of the skeleton; `output_off`,
  `output_on`); `unit/seed_void/input`. Study `run.py` scripts naming the keys are RECORDS and are not edited.
- **Gates:** tree `witness/p2_18.md5` (controls d2aa 26, d2t 863; the P2-25 list `p2_25.md5` equals it at 0/2639),
  face `face_p2_24` (13 legs; `cmp p2_24 p2_25` IDENTICAL). Re-pin procedure: `witness/REFERENCE.md` "How to re-pin
  next time" (new list → diff → classify every path → declare → keep the old list). `retired_units/` holds retired
  tests (`constraint_ref_lock`, `face_cref`). `WHICH_TESTS.md:34, 46` name `removal_increment` and `jet_d2c`.
- **Records to touch:** `MODEL_CARD.md` §6 rows 22 (`flake_length`), 23 (`flake_coherence_length`), the §12 item 7
  dead/inert list (:231), `WHICH_TESTS.md`, `witness/REFERENCE.md`.

## Goal

Two builds, in order. Stage 1 (B1 + B2 + B3) must read 0/2639 and IDENTICAL; stage 2 (C2 + C3) is the declared re-pin.
Before each build: copy every file the stage will touch to `studies/p2_26_closed_switches/src_pre_<stage>/`.

1. **B1 — two parse defects, two dead members.** Move the `spall.h_spall_max` parse (`M:6302`) above the
   `spall.flake_h_max` check so the check reads the deck's value (bit-neutral: every production deck has
   `flake_h_max 0.0015 < h_spall_max 0.015`; a deck violating it now aborts loudly). Move `side_power_interval`'s
   parse out of the enthalpy branch (its member and consumers unchanged). Delete `fl_wtop` and `sp_surface_T_mode`
   (keep the two bools and the `surfT` string parse).
2. **B2 — stranded text and dead locals**, from the census list re-resolved; nothing that reaches an expression.
3. **B3 — the five closed switches, parse to consumer**, with their members, thermo columns (`len_fb_h0/q0`),
   `CheckPadRing`, `WeibullSizeFactor` (every caller returns to the bare flaw: under `cell` the factor was exactly 1),
   the `ell_flux` branch of `BuildSurfaceTColumns`, the `rec_fixed` branch of the march. `jet_T_ent_mode` keeps
   `fixed` / `exhaust`; `rec_fixed` aborts as an unknown value. Top-level keys join the retired-at-default table:
   `spall.flake_length` → `twocell`, `spall.flake_fraction` → `cf`, `weibull.size_volume` → `cell`, `weibull.V_ref` →
   `*` (inert number). Strict-scope keys (`pad_shield_open`, `jet_above_h`, `jet_T_above`, `jet_T_rec_fixed`) are
   simply unparsed: a setter aborts in the scope. Tests: `jet_d2c` drops KT, K15, A17–A19 (and its (k) checks);
   `lateral_flux_lock` drops the `WeibullSizeFactor` / `size_volume` / `V_ref` clauses and gains a name clause that
   none of the five switches' identifiers survive in `src/`.
4. **C2 — `surface_normal` and the coherence length.** Delete parse, message, members, the conditional registration,
   `ColumnSurfaceCos`, `SP_SN_COS_MIN`, the site-A `/cos θ` sampling, the `R:` uses, the `sn_cth`/`sn_icth`
   parameters of `DepthScanH` (constant 1.0 after the deletion: bit-safe) and `F:591`'s literals; delete the
   coherence cap (`R:1510–1511` block) and `spall_flake_coh_len`. Both keys join the table at `0` (and the
   face_heat BASE tokens are dropped so no argv writes `0.0`). Edit the 17 + 4 live setters (remove the line; the
   Meier, Kant-meiercfg, p0 and Hu decks are records whose next run changes by design — say so in each file's
   header line). Tests: `removal_increment` drops the SN leg and checks (e)–(e4); `seed_void/input` edited.
5. **C3 — `melt_aware` and `lambda_crit`.** Delete parse, members, `MeltSkeleton`, the refusal, the Λ_L gather,
   `LL_at_didx`, `f_skel` and the `ScanSigmaT` parameter (`R:128–146` region; the literal at `F:577`). Both keys join
   the table (`melt_aware` → `0`, `lambda_crit` → `*`). Edit the 10 setters. `unit/sp_melt_skeleton` moves to
   `retired_units/` with a one-line README entry (it tests the deleted feature). `lateral_flux_lock:308`'s stress
   regex loses `f_skel \* `.
6. **The retired-keys lock.** `unit/retired_keys/test` gains, for THIS packet's table entries, one R-clause (each key
   string appears in `src/` exactly once, in the table) and K legs on `gradient_blind/input`: each key at its default
   runs (rc 0), at a non-default aborts naming the key (`spall.surface_normal=1`, `spall.flake_coherence_length=0.004`,
   `spall.melt_aware=1`, `spall.flake_length=flux`, `weibull.size_volume=flux`), and one strict-scope probe
   (`surface_patch.pad_shield_open=0.5` aborts). MUTATIONS: one per new clause.
7. **Records.** `MODEL_CARD.md` §6 rows 22–23 marked deleted (this packet), §12 item 7 updated; `WHICH_TESTS.md`
   (`sp_melt_skeleton` retired, `removal_increment` SN leg gone, `jet_d2c` RECF legs gone); `witness/REFERENCE.md`
   re-pinned to `p2_26.md5` naming every classified path, `p2_18.md5` kept; the daily log. Study
   `studies/p2_26_closed_switches/` with `z0.py`/`Z0.txt`, `chain.sh`, `src_pre_1/`, `src_pre_2/`, `BINARY_1.txt`,
   `BINARY_2.txt`, `GATE_*.txt`, `RESULTS.txt`. Nothing committed.

## Reachability ledger (step zero resolves every row; stop on a failed row)

| # | change | read at | consumed at | guard | witness |
|---|---|---|---|---|---|
| L1 | B1: the `flake_h_max ≥ h_spall_max` check reads the DECK's `h_spall_max` | `M:6156–6162` vs `M:6302` (today the check runs first) | the abort | `sp_flake_h_max > 0` | z0 prints both parse lines' order before and after; K-leg: `spall.flake_h_max=0.02 spall.h_spall_max=0.015` aborts (today it does not) |
| L2 | B1: `side_power_interval` parsed on every `jet_closure` | `M:8805` (inside the enthalpy branch today) | `v.side_power_interval` consumers `M:4337, 4477, 4483` | — | a `jet_closure = none` deck (CAL) with `surface_patch.side_power_interval=1.0` runs (rc 0) after; aborts before |
| L3 | B3: `WeibullSizeFactor` ≡ 1 on every gate deck (none sets `size_volume = flux`) ⇒ deleting the factor and its callers is bit-exact | `M:8347` (`if (!sp_size_flux) return 1.0;`), callers `M:6842`, `R:` | the flaw read | `sp_size_flux` false everywhere | stage-1 tree 0/2639, face IDENTICAL |
| L4 | B3: the `ell_flux` branch and `len_fb_*` columns exist only under `flake_length = flux` (registered conditionally `M:300–306`) | `M:3196–3240` | thermo header | `sp_flake_len_flux` | no gate header has `len_fb_h0` (z0 greps the p2_25 outputs); stage-1 0/2639 |
| L5 | C2: `sn_recession_superseded` is registered ONLY under `surface_normal` (`M:1009–1011`) and the `/cos θ` sampling changes `Cell_D` only where the key is 1 | `M:1009`, `M:6506–6537`, `R:1165–1271` | thermo header, `Cell_D` | `spall_surface_normal` | the C2 movers are EXACTLY `beam_void_closure` (169) + `dev2d` (64) = 233 paths; every other path identical |
| L6 | C2: `DepthScanH`'s `sn_cth`/`sn_icth` are 1.0 on every call once the column uses go; `F:591` passes `1.0, 1.0` today | `R:931–933`, `F:591` | the scan's depth sampling | — | face IDENTICAL (the face list passed literals 1.0 already) |
| L7 | C3: `f_skel = MeltSkeleton(Λ_L(z_v))` differs from 1 only where a scanned cell has Λ_L > 0 | `R:1174–1179` | the driving stress | `melt_aware` | the five unit melt-aware legs: max firing `T_top` 655 K < solidus 1448 K ⇒ identical (1,458 files); the s1 smoke legs (142) are the pre-registered candidate movers — z0 reads `Lambda_L` from their kept plotfiles and names which legs ever had Λ_L > 0 in a cell that fired |
| L8 | the CONJUNCTION a retired key needs: the string appears ONCE in `src/` (the table) AND a default-valued setter runs AND a non-default setter aborts naming the key | `M:157–176` | Parse | `pp.contains(key)` | retired_keys R-clauses + K legs (Goal 6); the string compare means `0` not `0.0` — the face_heat BASE tokens are dropped, not rewritten |
| L9 | units: `flake_h_max`, `h_spall_max`, `flake_coherence_length` [m]; `lambda_crit` [–]; `V_ref` [m³]; `jet_T_rec_fixed`, `jet_T_above` [K]; `side_power_interval` [s] | — | — | — | z0 lists each with its unit beside the deck values it deletes |
| L10 | the deck edits change no gate output EXCEPT through C2/C3's code (a deleted line at a live value = the code change itself; a deleted line at the default = nothing) | the 31 edited decks | — | — | z0's edit list marks each line DEFAULT or LIVE; every LIVE line is in a C2/C3 deck |

## Legs and the chain (exact tokens)

`PY=/Users/tzetze20/Desktop/code/.venv/bin/python`; `W=tests/MMWSpalling/studies/tree_sweep_0922/witness`;
`S=tests/MMWSpalling/studies/p2_26_closed_switches`. Serial, under caffeinate (the tree sweep's three lanes ran 57 min
under Spotlight indexing on P2-25: exclude the repo from Spotlight first, or expect it). Run order:

| step | what | expected |
|---|---|---|
| z0 | offline: the census items re-resolved to today's lines (one table: item, file:line, action); the deck-edit list with DEFAULT/LIVE per line; the gate-header grep (`len_fb_h0`, `sn_recession_superseded` in the p2_25 outputs); the C3 `Lambda_L` read of the s1 smoke plotfiles; the pre-registered mover list; `src_pre_1/` copied; `BINARY_pre.txt` | movers = 233 (C2) + the named s1 legs (C3) |
| s1 | build stage 1 (B1 + B2 + B3), `BINARY_1.txt` | compiles; `find src -newer` lists only the touched files |
| g1 | `guarded_post_all_par.sh p2_26a`; `tree_compare.py p2_26a p2_18 d2aa d2t`; `face_gate.py run p2_26a`; `cmp p2_24 p2_26a` | 0/2639 (26 / 863); IDENTICAL |
| k1 | `unit/retired_keys/test` (the B3 clauses and K legs), `unit/lateral_flux_lock/test`, `unit/buckle_lock/test`, `unit/face_removal/test --check`, `unit/jet_d2c/test` (RECF legs gone), L1/L2's two probe legs | all PASS; the probes behave as L1/L2 say |
| s2 | `src_pre_2/`; build stage 2 (C2 + C3) + the 31 deck edits + the test edits, `BINARY_2.txt` | compiles |
| g2 | `guarded_post_all_par.sh p2_26`; `tree_compare.py p2_26 p2_18 d2aa d2t`; classify EVERY differing path; `face_gate.py run p2_26`; `cmp p2_24 p2_26` | differing = exactly the pre-registered movers; IDENTICAL |
| k2 | `retired_keys` (all clauses), `lateral_flux_lock`, `buckle_lock`, `face_removal --check`, `removal_increment` (SN gone), `face_heat`, `beam_void_closure` (re-pinned; T3 FAIL by design), the `rankcheck p2_26` | all PASS / as declared |
| r | `REFERENCE.md` re-pinned to `p2_26.md5`; `RESULTS.txt` with the line count before/after (`wc -l` of the four headers) | — |

## Criteria

| # | run order | test | bar | gate / report / verdict | what a failure of this row means |
|---|---|---|---|---|---|
| z0 | 1 | offline: every census item resolved to a present symbol at today's line; the deck-edit list (31 files) with DEFAULT/LIVE; the mover list; `src_pre_1/` complete | every item found or struck with a reason; every LIVE line in a C2/C3 deck; the s1 `Lambda_L` read done | **gate** | the census and the code have diverged more than drift — re-plan the item |
| a1 | 2 (after s1) | tree gate stage 1 vs `p2_18` | **0 of 2,639** (controls 26 / 863) | **gate** | a "closed" switch was live on a gate deck, or a B2 edit reached an expression |
| a2 | 2 | face gate `cmp p2_24 p2_26a` | IDENTICAL | **gate** | a switch sat on the face list's path |
| a3 | 2 | `retired_keys` B3 clauses + K legs; `lateral_flux_lock`; `buckle_lock`; `face_removal --check`; `jet_d2c` | PASS | **gate** | a lock still names a deleted symbol, or a key is silently ignored instead of refused |
| a4 | 2 | L1 probe (`flake_h_max=0.02 h_spall_max=0.015` on `gradient_blind/input`) aborts naming both keys; L2 probe (`side_power_interval=1.0` on the CAL deck, `stop_time` one step) runs | both | **gate** | the parse order did not change, or moved the wrong way |
| b1 | 3 (after s2) | tree gate stage 2 vs `p2_18`: the set of differing paths | = the pre-registered movers, every path classified (`sn_recession_superseded` gone / `Cell_D` moved / `Lambda_L`-driven), nothing else | **gate** | an unlisted mover = a deletion touched a deck the pre-registration missed |
| b2 | 3 | face gate `cmp p2_24 p2_26` | IDENTICAL | **gate** | the face list read `f_skel` or the scan parameters after all |
| b3 | 3 | `retired_keys` all clauses + MUTATIONS + K legs (incl. `surface_normal=1`, `coherence=0.004`, `melt_aware=1` ABORT; defaults run) | PASS | **gate** | a retired key is silently inert (SR 1) |
| b4 | 3 | `removal_increment` (SN leg gone), `face_heat`, `beam_void_closure` T1–T4 on the re-pinned output, `rankcheck p2_26` | PASS / T3 by design / physics identical across ranks | **gate** | the deletion changed decomposition behaviour |
| m1 | 3 | the line count of the four headers before/after; the number of keys parsed before/after | report (expected ≈ 600–700 lines, 12 keys) | report | — |
| m2 | 3 | the C3 movers: which s1 smoke legs moved and in which fields; the maximum Λ_L in a fired cell on each | report: the deletion changed those legs' physics by exactly the skeleton's factor — a RECORD of the closed study, not a regression | report | a C3 move on a unit leg whose max firing T was 655 K means Λ_L was read where no melt existed — a bug in what was deleted, or in the reading |
| m3 | 3 | the 31 edited decks: each starts (`stop_time` one step) on the new binary | rc 0 on all | report | a deck still carried a retired key at a live value elsewhere |
| e0 | 4 | verdict | — | **verdict** | a1–a4, b1–b4 hold ⇒ the record's closed switches are gone, the two defects fixed, the production path untouched, and the gate re-pinned to `p2_26` with every moved path named |

## Pre-registered predictions
- **P1.** Stage 1 is 0/2639 and IDENTICAL: no gate deck sets any B3 switch (verified: 0 live setters), the B1
  order fix changes no accepted deck, B2 is text.
- **P2.** Stage 2 moves exactly `beam_void_closure` (169 files) + `dev2d` (64) through C2, plus whichever s1 smoke legs
  z0 finds with Λ_L > 0 in a fired cell (M-2000 and M-100 the likely ones; R-face-1400 possible; the R/F legs not);
  the five unit melt-aware test dirs (1,458 files) do not move.
- **P3.** The face gate is IDENTICAL after both stages (no face-list deck sets any of the twelve keys; the face
  list already passed `1.0, 1.0` to the scan and never read `f_skel` — if it did, b2 is the finding).
- **P4.** About 600–700 source lines leave (B3 ≈ 300, C2 ≈ 200, C3 ≈ 80, B1/B2 ≈ 80) and twelve input keys stop
  being parsed; the production decks and the rock/model-form files are untouched.

## Refuting evidence on disk
- `witness/p2_25.md5` = `p2_18.md5` at 0/2639 on binary `c052c348` — the baseline every stage compares to.
- `unit/*/output` (the p2_25 sweep): max firing `T_top` 655 K on the five melt-aware unit legs — if a C3 move appears
  there, this reading was wrong, not the deletion.
- `studies/p2_18_prereq/` (A2): the retired-at-default pattern, its K legs K0–K6 and `a2_apply.py` (the deck-edit
  tool and its DEFAULT/LIVE rule) — reuse, do not reinvent.

## What would make this inert, and the witness for each
- A deleted key becomes silently ignored instead of refused (top-level keys are not strict): the table entry and
  its K leg (b3) — this is the whole point of Goal 6.
- A B2 "comment" edit removes a line with code on it: a1 (0/2639) or a compile error.
- The s1 `Lambda_L` read looks at the wrong cells (not the fired ones): m2's bug clause; cross-check with the
  plotfile's `spall_event` field.
- `src_pre_<stage>/` missing a touched file: `find src -newer BINARY_pre` vs the copy list (z0 and s1).
- The tree sweep's binary moves mid-sweep: `guarded_post_all_par.sh` exits 2 (SR 9) — rerun.

## Guardrails
- Nothing on the production path changes: `rock_grimsel.inputs`, `model_form_flame.inputs`, every face-list deck,
  the face gate's thirteen legs, the beam, B, c_f, the span, R, the flaw — untouched.
- No deletion from C1 (DAMAGE), G (feet), H (`Removal.H`, column-top state, `robin_form` values, pinned) or D
  (Zhang) — they are their own packets. `RefreshKappaEff` and the B4 dead face-list evaluations are left to H.
- Study `run.py` scripts, `legs/`, `*_wit/`, `src_pre/` copies and `inputs_pre/` are records: not edited.
- Each stage's build is preceded by its `src_pre_<stage>/` copy of EVERY touched file (P2-25 DR-4).
- The old lists (`p2_18.md5`, `p2_25.md5`, `face_p2_24`) stay on disk; `p2_26.md5` is declared beside them.
- Nothing committed by this packet.

## Commands
```bash
cd /Users/tzetze20/amr_tools/alamo
PY=/Users/tzetze20/Desktop/code/.venv/bin/python; S=tests/MMWSpalling/studies/p2_26_closed_switches
W=tests/MMWSpalling/studies/tree_sweep_0922/witness; U=tests/MMWSpalling/unit
$PY $S/z0.py > $S/Z0.txt                                              # z0: items, edits, movers, Lambda_L read
EIGEN=$PWD/ext CPLUS_INCLUDE_PATH=/opt/homebrew/include \
  LIBRARY_PATH=/opt/homebrew/lib:/opt/homebrew/Cellar/gcc/15.2.0_1/lib/gcc/current make -j8   # s1, later s2
bash $W/guarded_post_all_par.sh p2_26a; $PY tests/MMWSpalling/studies/p2_9_face_removal/tree_compare.py p2_26a p2_18 d2aa d2t
$PY $W/face_gate.py run p2_26a; $PY $W/face_gate.py cmp p2_24 p2_26a
for t in retired_keys lateral_flux_lock buckle_lock jet_d2c; do $PY $U/$t/test; done; $PY $U/face_removal/test --check
# ... stage 2 edits, src_pre_2, rebuild, then:
bash $W/guarded_post_all_par.sh p2_26; $PY tests/MMWSpalling/studies/p2_9_face_removal/tree_compare.py p2_26 p2_18 d2aa d2t
$PY $W/face_gate.py run p2_26; $PY $W/face_gate.py cmp p2_24 p2_26; $PY $W/face_gate.py rankcheck p2_26
for t in retired_keys lateral_flux_lock buckle_lock removal_increment face_heat beam_void_closure; do $PY $U/$t/test; done
wc -l src/Integrator/MMWSpalling.H src/Integrator/MMWSpalling/*.H src/Numeric/MMWBeam.H   # m1 (14,021 before)
```
