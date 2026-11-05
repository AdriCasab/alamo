# MMWSpalling code census, 2026-11-05: the deletion plan

Read-only. Nothing under `src/` or `tests/` was edited. This document turns the four census parts
into an ordered deletion plan under the five scope decisions the user made on 2026-11-05. It
supersedes `docs/project/CODE_CENSUS.md` Part D (2026-10-08) as the order of work.

**Sources and citation form.** Every claim cites the part that established it, then the code line.
- `[P1]` = `census_2026-11-05/1_keys.md` (keys); `[P2]` = `2_regions_main.md` (MMWSpalling.H by
  function); `[P3]` = `3_removal_paths.md` (the two removal paths); `[P4]` = `4_tests_and_gate.md`
  (tests and the tree gate). `R7`, `A.4`, `§3.3` name a report or section inside a part.
- `M:` = `src/Integrator/MMWSpalling.H`, `R:` = `src/Integrator/MMWSpalling/Removal.H`,
  `F:` = `src/Integrator/MMWSpalling/FaceRemoval.H`.
- "the gate" = the 2,639-file tree gate (`studies/tree_sweep_0922/witness/`, list `p2_16.md5`).
  "W′" = the face-list witness gate that packet A1 creates.
- Where parts disagree, both readings are given and marked **⚠ parts disagree**.

**The scope decisions (user, 2026-11-05).**
1. The LEFM criterion is THE model. The Drucker-Prager damage + cohesive lineage retires. Hu decks
   on it re-express on LEFM or retire; the Hu conduction/thermoelastic/breakage decks on the
   global-material path merge with Zhang.
2. Flame capability is KEPT. One source must run flame and beam. The face list's flame-specific
   pieces get a volumetric meaning in the beam-port packet.
3. Zhang MERGES onto the phase-law path; then the separate global-material path and the implicit
   solver are deleted (implicit: decision noted, section 4 Q1).
4. Production is single-phase, but Voronoi/grain-boundary machinery is KEPT as a declared future study.
5. Both stress sources are KEPT (mechanics `stress_mf` and `local_thermoelastic`).
6. Closed study switches (per the record) retire. The column path (`Removal.H`) retires AFTER the
   beam, the vapour/melt rebuild and the feet are on the face list.

---

## 1. Summary

**Sizes** [P2 header, P3 header]: `MMWSpalling.H` 9,375 lines (6,292 code / 2,674 comment / 409
blank); `Removal.H` 2,580 (1,008 comment-only); `FaceRemoval.H` 937. Total 12,892.

**`MMWSpalling.H` by lineage** [P2 R7]: PROD-SHARED 2,969 · FLAME 1,686 (170 CLOSED) · COLUMN 1,634
(223 CLOSED) · DAMAGE 640 · DIAG 564 · MICRO 543 · BEAM 526 (45 CLOSED) · ZHANG 414 · MECH 354 ·
stranded text 45.

**After the deletions.** `MMWSpalling.H` ≈ **6,430 ± 20** (DAMAGE + COLUMN + CLOSED + ZHANG + stranded
text out, merge code in) [P2 R7]; ≈ **6,215** if the pinned closure (~135) and `PatchDiagnostics`
(80) also go (my arithmetic on [P2 R7]). `Removal.H` → 0, with ~330 lines moved into the face-list
file (seven helpers `R:23-260`, prologue `R:262-287, 345-368, 508-547`; [P3 B.1, B.2]). Kept as
whole lineages: FLAME (minus CLOSED), BEAM, MICRO, MECH. The beam, vapour and feet ports add code
that is not estimated here.

**Findings that change how the deletions must be done**
1. **The gate never exercises the face list.** No gate deck sets `surface_patch.face_list = 1`;
   every gate deck is on the column path [P3 :443; P4 §3.1, §3.2]. Any edit to face-list code is
   invisible to the gate, and deleting `Removal.H` empties 1,963 of 2,639 files with no successor
   [P4 §3.3(b)]. A face-list witness gate must exist first (A1).
2. **DAMAGE fields are in every plotfile.** `D`, `kappa_eff`, `k_eff` are registered plotted on every
   microstructure run (`M:502-504`), so deregistering them changes all **2,466** `Cell_D` hashes with
   no physics change [P4 §3.3(a)]. The **default criterion is `damage_law`** (`M:8693`, `M:9140`)
   [P4 tags; P1 row `spallation.model`]. `regime` (`M:534`) is NOT damage-only: the face list writes
   it (`F:825-837`) [P3 A.1].
3. **Kept decks name retired-lineage keys at their defaults.** Of 95 kept deck files:
   `kappa_damage_alpha` 78, `solver` 76, `damage.enabled` 59, `spall.sample` 43, `surface_normal` 25,
   `flake_coherence_length` 24 [P4 §1.2]. P2-10's abort-on-`contains` list (`M:119-122`) would break
   them. Needs a "retired at default" treatment (A2).
4. **`melt_aware = 1` is set in gate decks covering 1,600 files** (robin_pinned, robin_jet,
   jet_enthalpy, robin_feet, scalar_flaw, s1 smoke). Its deletion is hash-neutral only if no firing
   cell ever has Λ_L > 0, which must be measured [P4 §3.3(d3)].
5. **The vapour front double-spends beam power.** `BeamSource` adds (1−R)P into H (`M:2048-2066`,
   `M:4918-4948`); the vapour front turns the same (1−R)P into recession at ρL_tot (`R:1706-1709`,
   `R:1722`) and freezes the voided cells' H (`R:2242-2257`); no ledger term (`R:2431-2435`) [P3 C.4].
6. **Parse-order bug:** `flake_h_max < h_spall_max` (`M:5610-5612`) runs before `h_spall_max` is
   parsed (`M:5752`, behind the return at `M:5681`), so it compares with the default 0.010 [P1
   Surprise 2; P2 R4].
7. **`hc.heat.alpha` is not inert on production:** `R:297-301` aborts every spall step if it is ≤ 0,
   before the face-list dispatch; the damage-law h_spall it guards is computed and discarded [P1
   Surprise 3; P3 B.2 B1].
8. **The Sp_field site runs every production step for nothing the trajectory reads:**
   `UpdateSpAfterMechanics`, `UpdateSpClusters`, `BuildSurfaceTColumns`, three `GatherColumnSurface`
   sweeps, each with an MPI reduce; output is plotted only (`M:5992-6343`, `M:6353-6557`) [P2 R2, R3].
   On the column path it is load-bearing (the detach gate, `R:1162-1165`) [P3 A.2].
9. **Required-but-dead keys sit on the production decks:** `phase0.alpha_attenuation`
   (`rock_grimsel.inputs:38`), `phase0.name` (`:22`), `phase0.A_D` (`:39`), `surface_patch.face`
   and `pinned_idle_cycles` (`grimsel_calibration/input:94-95`, HOS_R argv) [P1 DEAD table; deck lines
   read for this document].

**Order**: A (gate, retired-at-default, helper move) → B (0-of-2639 deletions) → C (re-pins: DAMAGE,
surface_normal+coherence, melt_aware) → D (Zhang merge, then delete) → E beam port → F vapour → G feet
→ H `Removal.H` and the column-top state. D is independent of C, E, F and G.

---

## 2. The deletion list

### 2.0 Dependency graph

```
A1 face-list gate W'  ─┐
A2 retired-at-default ─┼─> B1-B4 (0/2639) ─> C1 DAMAGE ─> C2 surface_normal+coherence ─┐
A3 move helpers        │                     C3 melt_aware (measure first) ─────────────┼─> E beam ─> F vapour ─┐
   + prologue          │                                                                │                       ├─> H Removal.H
                       │                     G feet (needs A1, A3 only) ────────────────┼───────────────────────┤   (C4, C5 fold in)
                       └─> D1 Zhang merge ─> D2 delete separate path (independent of C, E, F, G, H)
```
C1 before C2 is a preference (do DAMAGE first, [P4 §3.3]), not a code dependency. E needs C2 and C3
only for its in-place parameter deletions.

Gate treatments used below:
- **0/2639** — the gate runs as-is and must show 0 differing files (P2-10 template [P4 §3.1]).
- **0/W′** — the face-list witness set from A1 must show 0 differing files.
- **re-pin** — the packet names every path it moves before it runs, and re-pins by D2z-1's procedure
  (new list → diff → classify every path → declare in `REFERENCE.md` → keep the old list) [P4 §3.1].

### (A) Prerequisites

#### A1. Add face-list witness decks to the gate (W′)

- **Why first.** The gate cannot see face-list code (finding 1). P2-10 hit the same blindness and
  built a separate witness set W, hashed on the pre-deletion binary (`W_bB.md5`) [P4 §3.1].
- **What.** Build W′ and hash it on today's binary before any source edit:
  - P2-10's W without its key-0 members: `unit/face_removal` T00_1_k1, TL35_1, TL45_1, SPH_1, BORE_1;
    `unit/face_list` CAV, STAIR; M0T (P2-9's M0 Höser leg minus 13 retired tokens) [P4 §3.1, §3.3].
  - `unit/face_buckle` F_0, F_S [P4 §3.3].
  - **Recommendation (not from the parts):** the two production legs HOS_R and CAL from
    `studies/p2_16_span_rule/legs/` at a short `stop_time`. HOS_R covers `jet_closure = enthalpy`
    and `buckle_span = 0.012`; CAL covers `jet_closure = none` and `buckle_span = 0` [P1 production
    path]. And one `unit/face_heat` path leg, because `face_heat = path` is kept (section 4 note).
- **Also fix in the same packet.** `witness/REFERENCE.md` still declares `d2z1.md5` the reference
  (stale) [P4 §3.1]. `WHICH_TESTS.md` has no row for `face_removal` or `face_list`, and its
  "removal" row names a DAMAGE test [P4 §1.1]. `unit/jet_d2c` runs in the gate but is not hashed
  [P4 §3.1].
- **Gate:** none (no source change). **Risk:** low. The W′ legs must be decomposition-stable; P2-10's
  W already compared column by column [P4 §3.1].

#### A2. The "retired at default" treatment

- **Problem.** P2-10's retired list aborts on `pp.contains(key)`, even at the old default
  (`M:119-122`) [P4 §1.2]. Top-level keys are not strict: deleting a parse makes decks silently
  ignore the key [P1 fact 1]. Strict scopes (`surface_patch.*`, `beam.*`, `grimsel.*`) abort on an
  unread key [P1 fact 2].
- **What.** One helper and one list, read at the top of `Parse` beside `M:88-130`:
  - top-level key, single old default (`damage.enabled = 0`, `conductivity.kappa_damage_alpha = 0.0`,
    `solver = explicit`, `spall.sample = top_cell`, `spall.surface_normal = 0`,
    `spall.flake_coherence_length = 0`, `spall.melt_aware = 0`, `spall.removal_energy = flake`, …):
    query it; equal to the old default → continue; anything else → abort naming the packet;
  - top-level key with no single default (`phase<i>.alpha_attenuation`, `phase<i>.name`): delete
    the parse and the production-deck lines (section 3.9); other decks become silently inert, which
    is correct because the value was never read [P1 DEAD];
  - strict-scope key (`surface_patch.face`, `T_ext`, `type`, `pinned_idle_cycles`,
    `pad_shield_open`, `jet_above_h`, `jet_T_above`, `jet_T_rec_fixed`): deleting the parse already
    aborts every setter. Edit the decks in the same packet; the numbers are bounded (section 3).
- **Locks.** `unit/retired_keys` gains one R-clause (the name appears only in the list) and one
  K-leg pair per list (default runs, non-default aborts) [P4 §1.2]. `unit/lateral_flux_lock`
  asserts the list text, so it changes too [P1 Surprise 10].
- **Gate:** 0/2639 and 0/W′ (no key's behaviour changes). **Risk:** low; the danger is a default
  typed wrong in the list. Take each default from the parse line in section 3.

#### A3. Move the seven shared helpers and the `Removal.H` prologue out of `Removal.H`, no behaviour change

- **What moves** [P3 B.1, "What deleting Removal.H requires first"]:

  | helper | from | callers |
  |---|---|---|
  | `ElasticConsts` / `ElasticConstants()` | `R:23-47` | `R:692, 750`; `F:471`; site A `M:6022, 6059` |
  | `RefillTempGhosts` | `R:49-62` | `R:2571`; `F:926` |
  | `ScanTAtZv` | `R:64-103` | `R:1397`; `F:567` |
  | `ScanKIcPa` | `R:105-127` | `R:1411, 1423`; `F:582, 585` |
  | `ScanSigmaT` | `R:129-148` | `R:1527`; `F:577` |
  | `DepthScanH` | `R:150-227` | `R:1543`; `F:591` |
  | `RemapProfile` | `R:229-260` | `R:2108`; `F:288` |

- **Prologue pieces the face list silently depends on** [P3 B.2]: the guards (`R:272-276`, stage B0),
  the per-step reset of `regime`, `RoP`, `spall_*`, `vapor_*`, `regime_field` (`R:345-368`, B3; the
  face list writes these on voided cells `F:830-835` and relies on the reset), and `k_top_global`
  (`R:508-537`, B6; read only for the CSV `k_top` column `F:895`). Move them to a function the
  dispatch at `M:1454` calls before `UpdateRemovalAfterCohesive`.
- **Keep the signatures.** Do not drop `ScanSigmaT`'s `f_skel` or `DepthScanH`'s `sn_cth`/`sn_icth` on
  the move. [P3 B.1] says "drop on the move", but the column path still passes real values while
  `melt_aware = 1` and `surface_normal = 1` decks are in the gate (`R:1527`, `R:1544`). Dropping them is
  bit-safe only at the face-list call sites, which pass literal 1.0 (`F:577`, `F:591`) [P2 R8]. They
  go in C2/C3 (section 2E explains how the beam port depends on that).
- **Locks.** `unit/face_removal --offline` (a3 i) reads `Removal.H` to assert the helpers are shared
  (`test:312-332`); `unit/lateral_flux_lock` opens `Removal.H` unconditionally at `test:89`;
  `unit/buckle_lock` `read()` opens it [P4 §1, §1.2]. Re-point all three at the new file.
- **Gate:** 0/2639 and 0/W′. **Risk:** low if it is one textual move. The fragment must stay a
  class-body include like `Removal.H` (CODE_CENSUS D.1; CLAUDE.md lesson 1).

### (B) Deletions that are 0-of-2639 today

These delete code that is off in every gate deck. Each still needs **0/W′**, because several sit on
code the face list runs (`FaceHr`, the march, the ledger).

#### B1. Parse defects, DEAD keys and dead members

- **Fix first (defect, not deletion):** move the `flake_h_max < h_spall_max` check (`M:5610-5612`) after
  the `h_spall_max` parse (`M:5752`), or move that parse above the return at `M:5681` [P1 Surprise 2;
  P2 R4]. The second choice also closes the latent onset-deck trap at `BuildSurfaceTColumns:2979, 3019`
  [P1 report ii E]. Also move `surface_patch.side_power_interval`'s parse (`M:8206`) out of the
  `jet_closure = enthalpy` branch: a `jet_closure = none` deck such as CAL cannot set it [P1 Surprise 8].
- **DEAD keys** [P1 DEAD; P2 R6]:

  | key | parse | delete with it |
  |---|---|---|
  | `microstructure.phase<i>.alpha_attenuation` | `M:347-348` (REQUIRED) | member `MineralPhase::alpha_atten` `M:8375`; `rock_grimsel.inputs:38` |
  | `microstructure.phase<i>.name` | `M:338` | member `MineralPhase::name`; `rock_grimsel.inputs:22` |
  | `microstructure.phase<i>.grimsel.T_ref` | `Grimsel.H:115` (strict) | nothing sets it; the overwrite at `M:7060` stays |
  | `surface_patch.face` | `M:8014` (strict) | 2 production + 30 validation + 12 unit decks set `zhi` (incl. gate inputs) |
  | `surface_patch.T_ext` | `M:8041` (strict) | member `M:7822`; 27 Hu validation + 4 unit decks |
  | `surface_patch.type` | `M:7940`, check `M:7952-7955` (strict) | 17 Hu decks |

- **Dead members** [P2 R6]: `fl_wtop` (`M:3290`, written `M:3808, 3816`), `sp_surface_T_mode` (`M:8965`,
  written `M:5555-5567`; only the two bools are read).
- **Gate:** 0/2639 and 0/W′. Removing `surface_patch.face = zhi` from gate inputs is an input edit
  with no output change. **Risk:** the strict-scope deck edits (≈90 files) are the whole risk; a
  missed one aborts at parse, loudly.

#### B2. Stranded text, dead locals and parameters

- **Stranded text** (45 lines, DEAD lineage) [P2 R6]: `M:622`, `M:1401-1404`, `M:2846-2853`,
  `M:3035-3039`, `M:3064-3068`, `M:3198-3233` (36 lines, 694 above its function at `M:3928`), `M:3234`
  tombstone, `M:3939-3954`, `M:4479-4484`, `M:5588-5598`, `M:5899-5905`, `M:5942-5949`, `M:6574-6578`,
  `M:8360-8361`, `M:8764-8774`, `M:8796`; and 12 stale line references (`M:637-638`, `787`, `2888`,
  `2967`, `2991`, `4500`, `5478-5481`, `5669`, `6584`, `7215`, `7703`, `8932`). Removal.H's P1b blocks
  (`R:639-658`, `R:1610-1628` and five more) die with H; do not spend a packet on them.
- **Dead locals and parameters** [P2 R6]: `pc` (`M:1926-1928`), `p_n` (`M:1945-1947`), `q_rob`/`q_los`
  (`M:2063`), the unread `phi` bindings (`M:1916-1917`, `M:4779-4780`, `M:5920`), `face_top`
  (`M:3912, 3920`), `z_lo` (`M:6855-6856`), `DX` (`M:5175`), `(void) dt` (`M:2950`), `(void) z_n`
  (`M:3518`), the redundant `(void)` pair in `TimeStepComplete` (`M:1192-1193`). The duplicated pad-ring
  block (`M:2166-2184` = `M:2205-2222`, third copy `FaceHr :3686-3689`) goes in B3 with the pad shield.
- **Gate:** 0/2639, 0/W′. **Risk:** nil if comments and locals only.

#### B3. Closed study switches that are off in every gate deck

| switch | evidence it is closed | code | gate | tests |
|---|---|---|---|---|
| `surface_patch.pad_shield_open` | ARCHIVE L13601 "D2aa … Wall lever measured and exhausted"; C.7.3 [P1 CLOSED] | parse `M:8004-8005`, `CheckPadRing` `M:7887-7895` (called `M:8338, 8355`), `BuildFlameColumns :2116` and the twice-pasted block `M:2166-2184`/`M:2205-2222`, `FaceHr :3686-3689`, abort `M:986-991` [P1; P2 R4, R6] | 0 hashed files [P4 §3.3(d6)] | `face_heat` R_pad restated (abort becomes strict-scope "unused key") |
| `surface_patch.jet_above_h` / `jet_T_above` | ARCHIVE L12600 D2v; L18980 "set on M2 only"; C.7.3 [P1 CLOSED] | parse `M:8182, 8187-8199`, `BuildFlameColumns :2156, 2190-2191`, `FaceHr :3678` (jet zeroing above the plane) [P1] | 0 [P4 (d6)] | `face_heat` R_above restated |
| `jet_T_ent_mode = rec_fixed` + `jet_T_rec_fixed` | ARCHIVE L8974 D2e diagnostic; C.7.3 [P1 CLOSED] | parse `M:8094, 8110-8116`; `JetEnthalpyMarch :2569-2571` | 0 hashed: `jet_d2c` legs KT/K15/A17-A19 run but `jet_d2c/output` is not hashed [P4 §3.1, (d5)] | retire those `jet_d2c` legs |
| `spall.flake_length = flux`, `spall.flake_fraction = weibull` | ARCHIVE L18497 "P1g-3 … REJECTED ON PROVENANCE"; C.7.3 [P1 CLOSED] | `BuildSurfaceTColumns` ell_flux `M:2974-3014`; thermo `sp_len_fb_h0/q0` (reg. `M:227-229`, members `M:8794-8795`); `R:2488, 2540`; `fr_abort` `M:1052-1054` [P2 R2; P3 D] | 0 [P4 (d7)] | none. **⚠** C.7.3 calls `flux` closed at P1g-3, yet `studies/p2_2d_kic` (80 legs) and `p2_2e_cref` (204 legs) ran it later [P4 (d7)] |
| `weibull.size_volume = flux`, `weibull.V_ref` | ARCHIVE L19170, L19250-19252 P2-2d "unvalidated, per-ladder anchor, no target" [P1 CLOSED] | `WeibullSizeFactor` `M:7752-7759`, site A `M:6286`, `R:1197-1204`, abort `M:9333` [P2 R4; P3 B.2] | 0 [P1 gate column "—"] | `lateral_flux_lock` loses its `size_volume=flux` abort leg (`:414`) [P4 §1.2] |

- **Code budget:** FLAME CLOSED 170 lines (pad shield, above-plane pair, rec_fixed, `robin_form = face`
  aborts; the last goes in C4) [P2 R7].
- **Locks:** one `retired_keys` R/K pair per switch; `lateral_flux_lock` name clauses
  (`pad_shield_open`, `jet_above_h`, `rec_fixed`) [P4 §1.2 (d)].
- **Gate:** 0/2639. **0/W′ is load-bearing here:** `FaceHr` runs on production (`M:3673-3691`) [P2 T1
  #37]. **Risk:** low; strict-scope keys abort their setters, which are only `face_heat` refusal legs
  and the `jet_d2c` legs.
- **Not in B3, deliberately:**
  - `surface.refine_depth`/`refine_height`: CLOSED is the weakest verdict in part 1 (not in C.7.3; only
    `studies/amr_band`) [P1 CLOSED]; it is the AMR band of the column path (`M:692` requires it for AMR
    + removal) [P2 R4]. It goes with H's AMR decision.
  - `weibull.detach_mode = connected_cluster` (C.7.2, validation-only) [P4 (d8)]: it acts in
    `UpdateSpClusters`, a column-path consumer; it dies with H.
  - `face_heat = path`: **⚠ parts disagree.** [P2 header, R7] counts it (~385 lines) as deletable
    because "the P2-11b review recorded that the user OK'd deleting the path code". [P1 FLAME-CAP row]
    says it was "deleted and RESTORED by the user in P2-12 (ARCHIVE L24861) — NOT closed". The later
    record wins; it is kept under decision 2.
  - `robin_form = face` as a control: 1,003 gate files → C4.

#### B4. Dead computations on the face-list step

[P2 R6] and [P3 D]:
- `AdvanceMicrostructure` evaluates `SurfaceCellFlux` for every in-patch top cell (`M:2003-2009`),
  then zeroes it for every listed cell (`M:2020`). Every in-patch top cell is listed (same radius test
  in `BuildFaceList`).
- `JetEnthalpyMarch` under the face list: the `jTc`/`jkc` gather with two `ReduceRealSum`
  (`M:2485-2504`) feeds only the `SurfaceCellFlux` + `PinRule` call in `q_col` (`M:2603-2611`), whose
  result is replaced by the listed faces' sum (`M:2617-2629`).
- `EnergyLedger`: `SurfaceCellFlux` per listed cell, zeroed at `M:4179`.
- `RefreshKappaEff` at `F:929` (and `R:2576`) writes `kappa_eff_mf`, which the kernel overwrites
  in the same step at `M:2072` [P2 R6]. **Measure before deleting the `R:2576` call** (it is on the
  column path, so it is 0/2639 only if nothing reads `kappa_eff` in between).
- **Gate:** 0/2639 (no gate deck runs the face list) and **0/W′ is the real test**. **Risk:** a
  floating-point side effect would surface only on W′, which is why A1 comes first.
- The Sp_field site and its sweeps (finding 8) are NOT here: they are plotted, so removing them
  changes hashes; they go in H.

### (C) Deletions that need a declared re-pin

#### C1. DAMAGE (Drucker-Prager damage + cohesive), in two packets

Do DAMAGE before COLUMN: it removes only 8 gate files if the plot fields are kept [P4 §3.3].

**C1a — computations and keys out; fields kept registered.**
- **First:** A1, A2. Make `spallation.model` required or default it to `sp_weibull` (default is
  `damage_law`, `M:8693`, `M:9140`); every gate deck except `spall_event` sets it explicitly [P4 §3.3(a)].
- **Delete in `MMWSpalling.H`** [P2 T1, R7]:

  | piece | lines |
  |---|---|
  | `ParseDamageSettings` | `M:5013-5110` |
  | `DamageGrowthRHS`, `IntegrateDamageRK4`, `UpdateDamageAfterMechanics` | `M:5130-5142`, `5144-5160`, `5162-5271` |
  | `ParseCohesiveSettings`, `UpdateCohesiveAfterMechanics`, `ResetCohesiveDiagnostics` | `M:5273-5324`, `5326-5389`, `5391-5398` |
  | `PrincipalStressRatio`, `FiniteMatrix`, `SpallSampleMode` | `M:5883-5897`, `6917-6923`, `8525` |
  | members | `M:8450-8501`, `M:8507-8528` |
  | Parse DAMAGE lines (47), `ParseSpallSettings` damage-law keys after `M:5681` (65), `ParseSpallationModel` (10), `InitializeMicrostructure` damage/cohesive reset (28) | per [P2 T1 #3, #63, #106, #88] |
  | aborts | `M:627, 630` (damage/cohesive need microstructure), `M:5069`, `M:9170, 9173` (sp_weibull refuses cohesive/damage) [P2 R4] |
  | `UpdateModel` damage branches, `Regrid` read | `M:1264, 1290-1292`, `M:1211` [P1 DAMAGE rows] |

- **Delete in `Removal.H`** [P3 B.2, B.4]: B1 damage-law h_spall `R:289-321` (**runs every production
  step**, finding 7), B9 PASS 1.5 `R:560-626`, the damage branch `R:1632-1661`, the regime choice's
  damage half `R:1958-1964`.
- **The κ-damage factor.** `kc = kappa·exp(−a_D·D)` appears 8 times: `M:1923, 1943, 2500, 3152,
  3757, 4149, 5125` and `F:454` [P2 R1]. With `kappa_damage_alpha = 0.0` in every other gate deck it is
  `exp(0) = 1` exactly [P4 §3.3(a)], so replacing it by `kappa` is bit-safe there. `RefreshKappaEff`
  (`M:5112-5128`) becomes a copy. The face list writes `D = 0` on voided cells (`F:825-837`) [P3 A.1]:
  keep the write while `D` stays registered.
- **Keys** (35 + `damage.ic.*`) [P1 DAMAGE]: retired-at-default for `damage.enabled` (59 decks),
  `conductivity.kappa_damage_alpha` (78), `spall.sample` (43) [P4 §1.2]; `microstructure.phase<i>.A_D`
  is REQUIRED on every microstructure deck (`M:349`) and set 0.0 on production (`rock_grimsel.inputs:39`)
  [P1 row]: retire at default 0.0. `spallation.model = damage_law` and `spall.sample = column_max |
  shallow_max` refuse by value [P4 retirement table].
- **Tests** [P4 §1, retirement table (a)]: retire (move to `retired_units/`) `dp_yield`, `gb_cohesive`,
  `spall_event`, `regime_low_high_power`. Re-point `amr_microstructure_regrid` (drop `test:137-198`) and
  `heterogeneous_kappa` (keep checks 1 and 4); both already FAIL (rc 1). Validation: retire
  `hu_spall_onset` (+ `diag/`) and `hu_end_to_end` unless Q5 re-expresses them; extensions
  `hu_spall_onset_mmwbeam`, `hu_spall_onset_voronoi` become dead records [P4 §2, §4.1].
- **Locks:** `retired_keys` R-clause "no `DamageGrowthRHS`, `IntegrateDamageRK4`,
  `UpdateDamageAfterMechanics`, `UpdateCohesiveAfterMechanics`, `PrincipalStressRatio`,
  `damage_enabled`, `cohesive_enabled` in `src/` outside the list" + K-leg "`damage.enabled = 1` aborts
  naming the step"; `lateral_flux_lock` ADD the same name clause (P2-10's `DELETED` regex,
  `test:931-946`) [P4 §1.2].
- **Gate:** **re-pin 8 files** (`spall_event` output + output_no_event) dropped from `hash.sh`; the
  other 2,631 must be 0. 0/W′.
- **Risk:** medium. `hc.heat.alpha`'s only production reader is the abort at `R:297-301`; after C1a it
  is truly inert, but it is still REQUIRED by `HeatConduction.H:79` [P1 report i], so its deck lines
  stay. Grain-boundary flags (`is_grain_boundary`) are MICRO and stay; the cohesive law only read them.

**C1b — deregister `D`, `kappa_eff`, `k_eff`.**
- **What:** the registrations at `M:502-504` and their member fabs [P4 §3.3(a)]. `regime` (`M:534`)
  STAYS: the face list writes it (`F:825-837`); only its "damage_law back-compat" comment changes.
- **Gate:** **whole-tree re-pin of 2,466 `Cell_D` files**, with an audit showing only the component
  count moved (per D2z-1's procedure) [P4 §3.3(a)]; W′ re-pinned the same way.
- **Risk:** CODE_CENSUS B1, "registration order is a hidden contract". Do it alone, no other change.

#### C2. `spall.surface_normal` + `spall.flake_coherence_length` (one re-pin)

Both are set in the same two gate decks, so one re-pin covers both.
- **Gate legs:** `unit/beam_void_closure` (169 files) and `dev2d` (`meier/input_2d_dev`, 64) = **233
  files** [P4 §3.2, (d1), (d2)]. The thermo header loses `sn_recession_superseded` (`M:938-940`), and the
  criterion depth `/cos θ` changes `Cell_D` [P4 (d1)].
- **`surface_normal`, delete** [P2 R8; P3 A.4]:
  - parse + Message: `M:5448-5501` [P2] — **⚠ parts disagree** on the span: [P3] gives `M:5470-5499`;
  - members: `M:8593-8598` (`spall_surface_normal`), `M:8551-8557` (`sn_recession_superseded`),
    `M:8619-8622` (`SP_SN_COS_MIN`);
  - face-list refusal `M:1038`; `ColumnSurfaceCos` `M:5958-5990` (with the stranded Sp header
    `M:5942-5949`);
  - site A: `M:6081-6090`, `M:6189-6194`, `M:6220`, `M:6300-6304`;
  - column path: `R:548-558`, `R:1224-1231`, `R:1460-1463`, `R:1544`, `R:1558-1569`, `R:1807-1810`,
    `R:2000-2009`, `R:2403-2424`;
  - **`DepthScanH`'s `sn_cth`/`sn_icth` parameters** (`R:158-216`, used at `R:174-175, 179, 201-204,
    209-211, 215`) and the literals `1.0, 1.0` at `F:591`. After the column uses are gone the
    parameters are constant 1.0 everywhere, so this is bit-safe.
- **`flake_coherence_length`, delete:** parse `M:5763`, B23 `R:1788-1842`, `fr_abort` `M:1055-1056`
  [P1 CLOSED; P3 B.2].
- **Keys:** retired-at-default 0 for both. Decks: `surface_normal` 25 (`=0` ×21, `=1` ×4);
  `flake_coherence_length` 24 (`0.004` ×12, `0.0` ×7, `0.002` ×5) [P4 §1.2]. Production sets both
  at 0 (`model_form_flame.inputs:19, :41`). P2-17's face-list argv names both at 0 (14 of them) [P4 §4.2(iii)].
- **Tests:** `removal_increment` SN leg retires; `beam_void_closure` re-pins (it already FAILs T3 by
  design); `face_heat` BASE token dropped; `seed_void` and `meier/input_2d_dev*` re-pinned or
  re-keyed [P4 §1, (d1)]. Locks: name clause (`spall_surface_normal`, `sn_recession_superseded`,
  `icth`, `spall_flake_coh_len` absent) + K-legs.
- **Alternative:** fold C2 into H, since both legs leave the gate there anyway [P4 (d2)]. But E needs
  C2's `DepthScanH` parameter deletion (section 2E), so doing C2 before E is cheaper.
- **Risk:** low physics risk (both switches closed) but it is a declared re-pin of 233 files. The
  criterion half of `surface_normal` is an open decision (Q7); the recommendation is to retire it.

#### C3. `spall.melt_aware` — measure first

- **Measure (packet C3-0, no source change kept):** build a binary with the skeleton factor forced to
  1 and run the gate. Deletion is hash-neutral **iff no firing cell in any leg ever has Λ_L > 0**; the
  likely movers are s1 smoke M-2000 and M-100 at 15 MW/m² [P4 (d3)]. 1,600 files are at stake:
  robin_pinned 323 + robin_jet 448 + jet_enthalpy 288 + robin_feet 138 + scalar_flaw 261 + s1 smoke 142
  [P4 §3.2].
  - **⚠ parts disagree on the setter count:** [P1 Surprise 6] says "`spall.melt_aware = 1` by five gate
    inputs". [P3 gate table] names eight inputs (`robin_pinned/input_pinned`, `robin_jet/input_jet`,
    `jet_enthalpy/input_quarter`, `input_strip`, `robin_feet` both inputs, `scalar_flaw/input`,
    `studies/s1_surface_resolution/input`). [P4 §3.2] lists six test directories.
  - **⚠ parts disagree on whether it is closed:** [P1] classifies it BEAM, not CLOSED; [P2 T1 #96]
    flags it CLOSED; [P3 C.2, checklist] says "closed study, delete it"; [P4 (d3)] calls it a
    contradiction, because C.7.3 lists it closed while the plan ports "melt". Decision Q6.
- **Delete** [P2 R8; P3 A.4]: parse + `melt.lambda_crit` `M:5502-5526` (**⚠** [P3] gives
  `M:5510-5526`); members `M:8599-8605`; `MeltSkeleton` `M:8606-8618`; refusal `M:1037`; the Λ_L
  gather `R:1017-1057`; `LL_at_didx` `R:1366-1386`; the skeleton `R:1464-1479`, `R:1501-1502`;
  **`ScanSigmaT`'s `f_skel` parameter** (`R:128-146`, multiplies `R:138, 143`) and the literal at
  `F:577`.
- **Keys:** retired at default 0 (11 decks name it; `=1` ×9) [P4 §1.2]; `melt.lambda_crit` retired
  (10 decks). Production: `model_form_flame.inputs:22`.
- **Tests:** `sp_melt_skeleton` retires (already FAILs, rc 1) unless Q6 keeps it;
  `meier/input_2d_dev_melt` becomes a record [P4 (d3)].
- **Gate:** 0/2639 if the measurement says so; otherwise re-pin exactly the movers, named before the
  run. 0/W′. **Risk:** the measurement can show the skeleton DID fire in the s1 legs; then the deletion
  is a physics change on those legs and must be declared as one.

#### C4. `robin_form = cell / face` as key values — recommend folding into H

- **What it costs:** 1,003 gate files (`robin_face` 798, `robin_pinned` F 45, `jet_enthalpy` strip 118,
  s1 R-legs 42) [P4 (d4)]. `cell` is the parse default (`M:7962`), so retiring it needs a new default,
  which moves every Kant onset deck, `p0_steady_column`, `gradient_blind`, `surface_T_recon` and
  `convective_patch` [P4 (d4)].
- **What it saves: almost no code.** The cell form (`SurfaceCellFlux :4035-4044`) is production's
  loss path for out-of-patch column tops. The face form's Newton solve (`:3987-4022`) is `FrFaceFlux`'s
  fallback (`F:148`) and `RecedeAlongLines :523` [P2 R2]. Only the key values, their aborts
  (`M:704, 713`, `M:7970`) and the `robin_face` test retire.
- **Recommendation:** leave the values until H, when the gate is rebuilt anyway; then refuse them by
  value ("must be pinned or surface") and re-point `robin_face` to a `surface`-form replica [P4
  retirement table (d4)].

#### C5. The pinned closure (only if Q4 retires it) — fold into H

- **Code** (~135 lines) [P2 R7]: `M:717-740`, `M:828-838`, `M:1891-1897`, `BuildPinEffective`
  `M:3076-3093`, `PinRule` `M:3098-3115`, `SurfaceCellFlux` pinned branch `M:4023-4033`, members
  `M:7614-7632`, `SurfacePatch::Parse` `M:7976-7998`; the pin record `R:2119-2138` and pin state
  [P3 B.2 B28, D].
- **Production decks:** `surface_patch.pinned_idle_cycles = 1.0e9` is parsed under `surface` too, in a
  strict scope (`M:7976-7993`), so the production lines must go with the parse:
  `kant/grimsel_calibration/input:94` and the HOS_R argv token [P1 FLAME-CAP row]. Unit setters:
  `face_heat`, `face_list`, `face_removal`, `removal_increment`, `jet_d2c`, `jet_enthalpy`, `robin_feet`,
  `robin_jet`, `robin_pinned`; validation `meier/sp_meier_pilot` ×2 [P1 row].
- **Gate:** the pinned legs are most of the column gate (robin_pinned, robin_jet, jet_enthalpy quarter,
  robin_feet) [P4 §3.2]. Doing this before H means re-pinning files that H deletes anyway.
- **Tests:** `robin_pinned` retires; `pin_clock` retires with pinned [P4 §1].

### (D) The Zhang merge, then the separate-path deletion

#### D1. Merge Zhang onto the phase-law path

- **What** [P2 R1]: Zhang as `microstructure.phase0.law = zhang` in a single-phase expression
  microstructure. Its `Material` goes through the same polymorphic `Table::Build` the Grimsel law uses
  (`M:7051-7066`). `UpdatePropertyFields` (`M:4512-4587`) is generalised from "grimsel" to "any
  T-dependent law". Estimate: +40-80 lines [P2 R7].
- **Decks to re-express** [P4 §2]: `zhang/zhang_oglesby/input` (its test reads only `losses.*` and the
  plotfile, so it survives a deck rewrite); `hu/hu_conduction`, `hu/hu_thermoelastic`,
  `hu/hu_thermoelastic_amr`, `hu/hu_breakage_index` (all ZHANG/MAT + PATCH-T; the micro path already
  carries PATCH-T through `M:1996`, and AMR on the micro path exists in `unit/amr_microstructure_regrid`).
  Mechanics stays (decision 5). The breakage index is a stress indicator, not damage.
- **Unit tests to re-point** onto a one-phase microstructure: `beam`, `equilibrium`, `stefan`,
  `convective_patch` [P4 §1]. Their H is in kelvin today (`M:6680-6690`), so their references re-pin
  [P2 R1]; see defect 5.2.
- **Implicit:** decision Q1. If Zhang re-runs explicit, D1 includes a timing check on `zhang_oglesby`.
- **Gate:** 0/2639 (no gate deck is ZHANG) [P4 §3.3(c)], 0/W′. The witness is the re-pointed
  unit tests and the Zhang and Hu decks against their validation references. **Risk:** medium.
  Zhang validation is the thermal anchor (CLAUDE.md lesson 4). The re-expressed deck must reproduce the
  `zhang_oglesby` reference, not today's binary, wherever the two old kernels differed (cavity
  attenuation, defect 5.3).

#### D2. Delete the separate path

- **Delete** [P2 R1, T1]:

  | piece | lines |
  |---|---|
  | `AdvanceConstAlpha` | `M:1472-1511` |
  | `AdvanceMaterial` | `M:1513-1633` |
  | `AdvanceMaterialImplicit` (per Q1) | `M:1642-1771` |
  | `BuildEnthalpyTable`, `SetPhaseFractionsConstAlpha` | `M:6626-6630`, `M:6635-6642` |
  | global/closed-form branches of `ConvertTtoH`, `InvertHtoT` | `M:6676-6690`, `M:6724-6757` |
  | `material` + `microstructure` abort | `M:1136` |
  | `enum class SolverKind` and members `material`, `material_table`, `implicit_*`, `T_ref/T_m/T_v/lh_*` | `M:7575`, `M:7567-7580` |
  | Parse | `M:155-161`, `M:174-181`, `M:252-267` |
  | `UpdateModel`'s homogeneous `alpha[]/eta` branch (MECH, reachable only without a microstructure) | `M:1361-1398` |
  | the second `hc.heat.alpha` parse into `therm_alpha` (after C1a removed `R:297`) | `M:153` |
  | the four-way dispatch | `M:1457-1460` (becomes one call) |

- **Keys** [P1 ZHANG]: `material.type` (+ 30 grouped `material.zhang.*`/`material.constant.*`),
  `solver`, `implicit.{max_iter, tol_abs, tol_rel, verbose}`, `phase.{T_ref, T_m, T_v, lh_m, lh_v}`.
  `solver = explicit` is named by 75 decks → retired at default; refuse `implicit` by value [P4 §1.2,
  (c)]. `material.type` (11 decks) aborts on presence once D1 has rewritten them.
  **Also dies:** the top-level `alpha` eigenstrain key, whose only reader is `UpdateModel :1392`
  inside the deleted branch [P1 MECH row]. Production sets `alpha = 0.0 0.0`
  (`hoeser/input_lab:95`, `grimsel_calibration/input:57`) → retired at default.
- **Kept:** `surface_patch.T_expr` (classed ZHANG in [P1], but it is the PATCH-T temperature that Hu
  and Rossi use on the micro path); `losses.*`.
- **Locks:** `lateral_flux_lock` or `retired_keys`: the three Advance functions and
  `SolverKind::Implicit` are gone; `material.type` aborts [P4 retirement table (c)].
- **Gate:** 0/2639, 0/W′. **Risk:** low once D1 is accepted. 414 ZHANG lines + ~39 MECH lines [P2 R7].

### (E) The beam-port packet (separate; not a deletion)

The beam port adds behaviour. It is listed here because it is the hinge between the deletions and H.

- **First:** A1 (W′ must hold for flame decks through the port), A3 (recommended, so the helpers it
  edits are already out of `Removal.H`). Its in-place parameter deletions need **C2 and C3**: `f_skel`
  and `sn_cth`/`sn_icth` are bit-safe to delete only at the face-list call sites, which pass literal
  1.0 (`F:577`, `F:591`) [P2 R8]. The column path passes real values (`R:1527`, `R:1544`) while
  `melt_aware`/`surface_normal` decks remain. If E runs before C2/C3, it leaves the parameters and
  keeps passing 1.0.
- **Deletions E does NOT need:** C1 DAMAGE, D Zhang, B3 closed flame switches, C4, C5,
  `flake_coherence_length` (the face list already refuses it), and H.
- **Lines it must change or bypass** [P3 A.4]:

  | # | site | lines | why |
  |---|---|---|---|
  | 1 | parse gate | `M:994-1000` | requires `convective_flame` + `h_expr` + `robin_form = surface` |
  | 2 | parse aborts | `M:1033` (vapour → F), `M:1037` (melt_aware, deleted in C3, not lifted) | |
  | 3 | list-build gate | `M:1883-1886` | `fl_on` needs `surface_patch.enabled && !prescribed_T`; `BuildFlameColumns` is the only caller of `FaceListHeat`/`FrAssignTs`/`FaceListTg` (`M:2228-2235`) |
  | 4 | patch radius | `M:3734-3736` | faces outside `surface_patch.radius` are not listed; a beam deck has no patch |
  | 5 | BC face temperature | `F:138-156`, `F:162-253` | T_s from flame h/T_g or `ReconFromX`; under a beam only losses set the Newton |
  | 6 | R1 reconstruction | `F:229, 234-242, 469, 500-506`; `M:9000-9124` | gas cap, decaying-inward test, exponential T_s and h_f |
  | 7 | R1 fallback | `F:507-550` | ell_q = k(T_s−T0)/q_net, noskin: the flame-defined skin length |
  | 8 | R2 scan anchor and walk | `F:551-600`; `R:64-103`, `R:159-227` | walk stops at the first failure; search depth h_spall_max/(mu·dx) |
  | 9 | B thickness and strain | `F:607-635` | x = c_f·ell_eff; ε_th at T_s |
  | 10 | R3 heat remap | `F:701-717`; `R:236-260` | h = `rec.h_f`; exponential sub-cell blend with ell_eff |
  | 11 | window sizing | `F:368-369, 377-380` | Dmax, Lg from h_spall_max, flake_h_max |
  | 12 | beam's sight of a partial cell | `M:1924`, `M:2048-2061`, `SolidFacePath` `M:4328-4331` | a cell with f < 1 absorbs as a full cell; the kernel never reads `fr_f` |

- **The flame-specific pieces that need a volumetric meaning** (decision 2) [P3 A.3]:
  1. *Where h and T_gas enter:* `FrFaceFlux` (`F:148, 153`), the removal fallback (`F:518-525`, so the
     flake thickness is h-and-T_gas-defined on every fallback line), the gas-temperature cap
     `SurfaceGasTempCap` (`M:9117-9124`), and the list's existence only inside the patch
     (`M:3734-3736`, `M:1883-1885`).
  2. *A profile that decays inward:* `RawSurfaceX` returns x > 0 only if T_c > T_next > T0
     (`M:9000-9010`); `ReconFromX` extrapolates an exponential (`M:9036-9053`); `RemapProfile`'s blend
     (`R:240-247`); `DepthScanH` stops at the first sample with K_I < K_Ic and bisects assuming
     pass → fail with depth (`R:171-225`), so **a subsurface peak is never reached**.
  3. *The flat case:* under a beam-only deck h = 0, so q_net = −losses < 0 on every flat line, and
     every flat line is noskin (`F:544-549`). Volumetric heating makes T_next > T_c normal, so **on the
     beam today the face list would recede only through `through`** [P3 A.3 item 3].
  4. *mu:* /mu is applied to normal lengths only (`F:241, 505, 535-536, 377-380`). The cancellation
     holds only if T depends on normal distance alone; Beer-Lambert deposits along −z only
     (`M:2033-2061`), so wall lines see a radial profile [P3 A.3 item 4].
- **Defects E should fix:** alpha(T) energy conservation (5.6); the uniform-only closure check
  (`M:4859-4861`); refuse beam + face list until the port lands (a beam with the face list is accepted
  today and deposits through column tops only [P2 R4, R5]).
- **In-place deletions** (bit-safe at the face-list sites once C2/C3 have landed): `ScanSigmaT`'s
  `f_skel` + `F:577` literal; `DepthScanH`'s `sn_cth`/`sn_icth` + `F:591` literals; `MeltSkeleton` if
  C3 did not already remove it.
- **Tests:** re-point `beam_void_closure` and `scalar_flaw` (event-log completeness and rank identity
  are wanted), and `seed_void`'s S2 beam-path claim [P4 §1]. Meier `input_energy`, `input_2d_dev`
  re-point; `input_drilling*` wait on Q2 [P4 §2].
- **Gate:** W′ must be 0 for every flame leg (the port must not move flame behaviour); new beam
  face-list legs are added to W′ as the port's own witness, hashed when accepted.

### (F) The vapour/melt rebuild on the enthalpy path

- **What is wrong today** [P3 C.1, C.3, C.4]:
  1. **Double-spend** (finding 5). `BeamSource` puts (1−R)P into H on every non-void cell
     (`M:2048-2066`, `M:4918-4948`); `ledger_P_beam` books it (`M:4227, 4295`). The vapour front
     turns the same (1−R)P into a recession rate rop_v = (1−R)Q/(ρ·L_tot) (`R:1706-1709`) or
     (1−R)P/(ρ·L_tot·A_beam) (`R:1722`), charging ρ·L_tot = ρ[c_p(T_vap_lo−T0) + L_m + L_v] per
     removed volume. PASS 3 voids the cells and freezes their H (`R:2242-2257`); the remap and
     flake sink are regime 1 only (`R:2075, 2268, 2431-2435`); no `ledger_` vapour term exists.
  2. The gate tests the column-top cell's H only (`H_top ≥ H_vap`, `R:1725`), but h_v = rop·dt can
     exceed dz, so the phi shift voids deeper cells whatever their H (`R:2201-2212`).
  3. `Q_local` uses the Gaussian waist formula even under `beam.profile = uniform`
     (`R:1703-1708` vs `BeamIncidentIntensity`'s uniform branch `M:4964-4970`); it ignores spreading
     at the surface depth and radiation at the vapour temperature.
  4. The ledger: regime 2 and regime 3 book nothing; whole-cell regime 1 without remap books nothing
     (`R:2431`); `ledger_err` cannot see any removal (void H stays in `E_H`, `M:4134`);
     `ledger_err_solid` would drift by the H of every vapour-voided cell.
  5. Melt: latent heat is already on the enthalpy path (`L_m`, `L_v` are plateaus of the per-phase
     H(T) table, `M:390-393`, `M:7039-7042`, `M:7063-7066`; `InvertHtoT` returns Λ_S/Λ_L/Λ_V,
     `M:277-279`). **No melt removal or evacuation exists on either path** [P3 C.2].
- **What the rebuild must be** [P3 B.3, checklist]: a per-face-cell regime-2 test on H (not on a column
  top), a recession charged against H rather than against beam power, per-cell regime 2 in
  `regime`/`regime_field`/`vapor_*`, its own ledger sink, and the spall-versus-vapour priority (`R:1941-1965`).
  *Recommendation (mine):* because L_v is already a plateau of H(T), the removal can void a listed
  cell whose H has reached H(T_vap_hi) and book that H into a vapour sink. No separate rate law
  spends the beam a second time.
- **First:** E. Lift `M:1033`. **Keys kept:** `vapor.enabled`, `vapor.T0`; `vapor.A_beam` and
  `vapor.prescribed_rop.*` exist only for the rate law and are set by no input [P1 BEAM rows; P1 Counts].
- **Tests:** `sp_v_n_regime` re-points (it scores the local-Q radial variation, which the rebuild
  replaces, so its assertion changes); the vapour half of `regime_low_high_power` is covered by it [P4 §1].
- **Gate:** additive on W′ (no current W′ leg runs vapour). **Risk:** high (new physics). It must close
  `ledger_err_solid` with the new sink.

### (G) The feet port

- **What moves** [P3 B.2 B5, B.3]: the D2c mechanical clearance `R:387-506`. An annulus column proud of
  its pad has its cells above `foot_clear_k` voided as regime 4 with D = 0, H frozen; the column phi
  shifts; `foot_mech_E/vol/cols` are booked (`R:486, 495-496`). On the list: each cell above the pad
  index gets f = 0 (`fr_f_mf`), removed = 1, regime 4, `foot_mech_E += H·V`, H untouched, plus the face-list
  carry (`fr_carry`) for the new face cell and the CSV regime-4 row (`R:2501-2507`). D = 0 is dropped
  if C1 has landed.
- **Why the face list refuses it today:** the block runs before the dispatch but would not zero
  `fr_f_mf` (`M:1035-1036`) [P3 B.2]. Lift that abort.
- **Already on the shared side:** `FeetDescent` (`M:2302-2419`) produces `foot_clear_k` (`M:7681`,
  written `M:2411-2417`); only `R:442-505` consumes it [P3 D]. The ledger term `foot_mech_E`
  (`M:892-903`, `M:4310`) stays.
- **Decks:** Meier `input_feet`, `input_feet_d2c` also need `robin_form = surface` + `h_expr` (they
  are pinned today) [P4 §2].
- **Tests:** re-point `robin_feet`, `jet_d2c`, `pad_clear_ledger` (the face list refuses
  `foot_body_clearance`, `M:1033`), `jet_enthalpy`, `robin_jet` [P4 retirement table (b)].
- **First:** A1, A3. Independent of E and F. **Gate:** additive W′ legs; W′ flame legs 0.

### (H) Delete `Removal.H` and the column-top state

- **First:** A1-A3, C1, E, F, G; decisions Q2 (bit), Q3 (Rossi), Q4 (pinned), Q8 (site A); the new W′
  reference declared in `REFERENCE.md`, with the 2,639 list kept as the record of the column era
  [P4 §3.3].
- **Gate:** **the 2,639 gate cannot pass and has no successor inside itself**: 1,963 of 2,639 files
  leave (only the 676 `robin_face` form files survive, and they are `cell`/`face` controls)
  [P4 §3.3(b)]. 0/W′ for everything H does not deliberately remove; a declared W′ re-pin for the
  plotted fields below.
- **What moves** (done in A3) [P3 "What deleting Removal.H requires first"]: the seven helpers, the
  guards, the B3 reset, `k_top_global` (or `ColumnTopSolid`, `M:4589-4611`, which already computes it).
  `M:1454` calls `RecedeAlongLines` directly.
- **What dies** (all of `Removal.H`, 2,580 lines, plus in `MMWSpalling.H`) [P3 D; P2 R2, R3, R7]:

  | state | lines | note |
  |---|---|---|
  | `BuildSurfaceTColumns`, `surf_*` arrays, `surf_x_hold` | `M:2866-3033`, `M:8834-8852` | output feeds diagnostics only under the face list |
  | `SharedSurfaceRecon`, `SurfaceRobinT` | `M:3046-3062`, `M:3069-3074` | `SurfaceRobinT`'s four consumers are overridden on every listed cell (`M:2015-2020`, `M:2618-2626`, `M:4180-4182`); `PatchDiagnostics :3159` is the live one |
  | `GatherColumnSurface`, `ColumnSolidRunBelowTop`, `ColumnSurfaceCos` | `M:5907-5940`, `M:4677-4714`, `M:5958-5990` | |
  | site A `UpdateSpAfterMechanics`, `UpdateSpClusters` | `M:5992-6343`, `M:6353-6557` | **see Q8: the Kant onset decks need a criterion without removal** |
  | plotted fields `Sp_field`, `Sp_cluster_id`, `T_surf_field`, `ell_eff_field`, `h_spall_field`, `crack_event_count`; thermo `recon_*` | reg. incl. `M:576-577`, `M:583-591` | W′ re-pin (or keep registered as constants and remove in a follow-up) [P2 R7] |
  | `rc_ell_prev`, pin state, `h_col_events` CSV, `DetachMode` + `connected_cluster` | `M:8858`, `M:7621`, `M:5536, 8904, 8913`, `M:8888` | |
  | `PatchDiagnostics` | `M:3122-3196` | reports `patch_P_robin` from a flux the face list does not apply (`M:2083-2084`) [P2 R6] |
  | AMR column repair: `RepairRemovalStateAfterRegrid`, `CheckSurfaceCovered`, the surface band in `TagCellsForRefinement` | `M:4349-4436`, `M:4442-4477`, `M:4982-5010` | the face list refuses AMR (`M:1042`); the dispatch already skips `CheckSurfaceCovered`. Port only if AMR comes to the list |
  | dead `E_sup_local` | `R:2068, 2217` | regime 5 is never assigned |
  | `UpdateSurfaceMaskFromPhi` | `M:6793-6806` | **production caller `F:928`**; read only when `follow_mask = 0` and by the ledger's missing-column count [P2 T1 #78, R7]. Delete with that ledger column, or keep |

- **Keys** that die: `spall.removal_energy` (production requires `flake`; no face-list reader
  [P1 Surprise 4]) → retired at default `flake`; `spall.surface_T` → retired at default
  `reconstructed`; `surface_patch.face_list` → mandatory, retired at default 1; `spall.clusters_csv`,
  `spall.h_col_events_csv`, `spall.prescribed_ratio`, `weibull.detach_mode`/`A_crit`, `bit.*` (per Q2),
  `spallation.sp.prescribed_stress.*`/`prescribed_T_surface.*` if site A goes (Q8).
- **Tests** [P4 retirement table (b)]: retire `removal_increment`, `seed_void` (its four call sites are
  column-top scans; `spall.seed_void` stays as a fixture tool), `surface_T_recon` (re-point check 1),
  `robin_pinned` (per Q4). Re-point `gradient_blind` (also the runtime fixture of `lateral_flux_lock`
  `test:704` and `retired_keys`), `robin_face` ledger legs, `face_removal` (drop key-0 legs and (a3 i)'s
  column half), `face_list` (drop OFF0/CAVU, clause iv), `face_buckle` (drop F_S_k0), `face_heat`
  (R_nolist). Validation: re-point `kant/sp_kant_onset_meiercfg` (successor `grimsel_calibration`),
  `p0_steady_column`; `rossi` per Q3.
- **Locks** [P4 §1.2]: `lateral_flux_lock` stops opening `Removal.H` (`test:89`); its 48
  `Removal.H`-reading sites are deleted or restated (`:187`, `:189/:213`, `:222`, `:257`, `:282-293`,
  `:357`, `:398`, `:414`, `:435`, `:607`; counts `:311`, `:316`, `:354`, `:401`, `:857-904`); its fixture
  moves to a face-list deck. `buckle_lock` drops B2's `Removal.H` half and B3 (abort under
  `face_list = 0`). `retired_keys` drops the `Removal.H` `SUBST` and moves its fixture. NEW clause: "no
  `UpdateRemovalAfterCohesive`, `Removal.H` or column-only names in `src/`". `WHICH_TESTS.md` is rewritten.
- **Studies:** the 35 column-path studies stop running; they are records [P4 §4.2(ii)].
- **Risk:** high in volume, low in physics, if E-G landed first with W′ green.

---

## 3. The key table, condensed

From [P1] (242 keys + 3 grouped rows). "Fate" maps each key to a packet in section 2. Defaults and
parse lines are [P1]'s. Retired-at-default = A2 treatment.

### 3.1 PROD (80 keys): stay, with these exceptions

| family | keys | note [P1] |
|---|---|---|
| `losses.*` | `epsilon`, `epsilon_high`, `h_conv`, `sigma_SB`, `T_amb`, `T_eps_break` | `T_eps_break = 1.0e30` is the default on production: the two-emissivity switch is unused; could retire with `epsilon_high` (not on any list) |
| `hc.heat.alpha` | REQUIRED at `HeatConduction.H:79` and `M:153` | production reader is the abort `R:297-301` only → inert after C1a; `M:153` copy dies in D2 |
| `microstructure.*` | `mode`, `nphases`, `phase_expr`, `phase<i>.{beta, Cp, E, fraction, kappa, law, mu, rho, T_ref}` | `beta, Cp, kappa, rho` REQUIRED, value-inert under `law = grimsel`; `E > 0 && mu > 0` load-bearing via the abort `M:6023` [P1 Surprise 5]. Two Poisson ratios survive (`grimsel.nu` used; `E/(2mu)−1` only without the law, `R:43`) |
| `microstructure.phase<i>.grimsel.*` | `alpha_cap_theta`, `cp.{a,b,th0}`, `E0`, `K_Ic0`, `lambda.{a,th0}`, `nu`, `rho` | strict scope (`M:372`), parsed for every phase |
| `spall.*` | `buckle_span`, `enabled`, `flake_c_f`, `flake_h_max`, `h_spall_max`, `h_spall_n_segments`, `removal_energy`, `stress_law`, `stress_source`, `surface_T` | `flake_h_max`: parse-order bug (5.1). `h_spall_n_segments`: set by nothing, read at default 0 → hard-code. `removal_energy`, `surface_T`: parse-time requirements only on the face list → retired at default in H |
| `spallation.*` | `model`, `sp.{a0, kic.law, kic.scale, n_segments, stress_relief}`, `confining.p` | `model` default `damage_law` → C1a. `sp.a0` REQUIRED, inert under Weibull (`rock_grimsel.inputs:64`) but read when `weibull.enabled = 0` (`F:452`) |
| `surface.follow_mask` | | required by `surface` form and face list → mandatory after H |
| `surface_patch.*` | `enabled`, `face_list`, `h_expr`, `jet_*` (14), `mode`, `nozzle_*` (4), `radius`, `robin_form`, `T_flame_expr`, `x0`, `y0` | `jet_T_ent_mode` value `rec_fixed` → B3. `robin_form` values → C4/C5. `face_list` → mandatory after H |
| `weibull.*` | `a0_ig`, `detach_mode`, `enabled`, `m`, `seed`, `V0` | `detach_mode = connected_cluster` dies in H. `m` is behind `sp_weibull` (`M:9164`) |

### 3.2 FLAME-CAP (17): KEPT (decision 2), except

| key | fate |
|---|---|
| `face_heat`, `barrier_expr` | KEPT (restored P2-12; ⚠ [P2] disagrees, see B3) |
| `foot_body_clearance` | KEPT; ported in G (`R:432` consumer moves) |
| `foot_npads`, `foot_pad_quantile`, `foot_r_inner`, `foot_r_outer`, `foot_rule`, `foot_standoff`, `nozzle_feed_max` | KEPT (`FeetDescent`) |
| `jet_dr`, `jet_free_surface`, `jet_fs_aspect`, `jet_fs_exponent` | KEPT (march) |
| `h_conv`, `T_flame` (constant-coefficient flame) | KEPT; Hu decks and gate tests set them |
| `pinned_idle_cycles` | dies with the pinned form (C5/Q4); set 1.0e9 on production, inert under `surface` (Message `M:770-774`), strict scope |

### 3.3 BEAM (32): KEPT, except

| key | fate |
|---|---|
| `beam.*` (12) | KEPT. `alpha_expr` is the only attenuation (5.6) |
| `bit.*` (7) | per Q2 (`R:1894-1931`); face list refuses (`M:1034`); used by `meier/input_drilling*` |
| `spall.melt_aware`, `melt.lambda_crit` | C3 / Q6 |
| `microstructure.phase<i>.{fus_width, L_m, L_v, T_m, T_vap_hi, T_vap_lo}` | KEPT (latent plateaus in H(T); F needs them) |
| `vapor.enabled`, `vapor.T0` | KEPT for F |
| `vapor.A_beam`, `vapor.prescribed_rop.{enabled, expr}` | die with the rate law in F; `prescribed_rop` set by no input |

### 3.4 ZHANG (12 + 2 grouped): die in D2

`implicit.{max_iter, tol_abs, tol_rel, verbose}` (only `zhang_oglesby`), `material.type` +
`material.zhang.*` (20) + `material.constant.*` (10), `phase.{lh_m, lh_v, T_m, T_ref, T_v}` (only
`unit/stefan` + Hu thermoelastic decks), `solver` (75 decks at `explicit` → retired at default).
`surface_patch.T_expr` is KEPT (PATCH-T on the micro path, Hu/Rossi). Note: under `convective_flame` it
is consumed and ignored (`M:8034`).

### 3.5 MICRO (5): KEPT (decision 4)

`conductivity.h_gb0` (read every step, inert at 0; 78 decks), `microstructure.grain_expr`,
`number_of_grains`, `seed`, `weibull.a0_gb` (REQUIRED `M:9294`, inert on the grain-free rock,
`rock_grimsel.inputs:67`).

### 3.6 MECH (14): KEPT (decision 5), except

| key | fate |
|---|---|
| `alpha` (top-level eigenstrain) | dies in D2 (only reader `UpdateModel :1392`) |
| `phase<i>.alpha_beta.*` (5) | KEPT (MICRO/MECH study); `linear_strain` set by no input |
| `spallation.sp.kic.{K_Ic0, T_melt, T_ref}` | KEPT (linear law, Kant verification decks); silently ignored under `grimsel_table` [P1 report ii H] |
| `spallation.sp.prescribed_stress.*`, `prescribed_T_surface.*` | change only `Sp_field` (`M:6006-6007`, `6147`, `6221-6225`) → fate follows Q8 |
| `weibull.A_crit` | dies with `connected_cluster` in H |

### 3.7 DAMAGE (35 + grouped `damage.ic.*`): die in C1

`cohesive.{delta_c, delta_max, enabled, f_t, G_c, K_n, prescribed_delta.enabled, prescribed_delta.expr}`
(only `gb_cohesive`; `G_c`/`f_t` are set by no input and are silently ignored when both deltas are
given, `M:5293`); `damage.{cohesion, enabled, gb_multiplier, ic.*, m, n, phi_deg,
prescribed_stress.enabled, prescribed_stress.{xx,yy,zz,xy,yz,zx}_expr, psi0, stiffness_floor}`;
`conductivity.kappa_damage_alpha` (78 decks, `=0.0` ×75; production `hoeser/input_lab:105`,
`grimsel_calibration/input:81`); `microstructure.phase<i>.A_D` (REQUIRED, `rock_grimsel.inputs:39`);
`spall.{alpha, C_h, damage_threshold, prescribed_ratio.enabled, prescribed_ratio.expr, sample,
shallow_depth, t_spall, xi_s}`. `spall.alpha`/`C_h`/`t_spall` feed the h_spall computed and discarded
on every production step (`R:296-298`).

### 3.8 CLOSED (12), with evidence

| key | parse / default | evidence [P1 CLOSED] | gate setters | fate |
|---|---|---|---|---|
| `spall.flake_coherence_length` | `M:5763` / 0.0 | C.7.3 "superseded for the feet line (set to 0 there)"; face list refuses > 0 (`M:1055-1056`) | `beam_void_closure` 0.004, `meier/input_2d_dev` 0.004 | C2 |
| `spall.flake_fraction` | `M:5632` / cf | ARCHIVE L18497 (P1g-3 REJECTED); face list refuses `weibull` (`M:1052-1054`) | — | B3 |
| `spall.flake_length` | `M:5626` / twocell | ARCHIVE L18497 "REJECTED ON PROVENANCE"; refuses `flux` | — | B3 (⚠ two later studies ran `flux`, [P4 (d7)]) |
| `spall.surface_normal` | `M:5474` / 0 | C.7.3 (Step 20; recession half superseded by D2z-1, ARCHIVE L13103); own Message `M:5494-5501`; refusal `M:1038` | `beam_void_closure` 1, `meier/input_2d_dev` 1 | C2 (criterion half Q7) |
| `surface.refine_depth` | `M:687` / 0.0 | **weakest verdict**: not in C.7.3; `studies/amr_band/RESULTS.md` G1/G2 PASS, wall-time ratio 1.00; face list refuses AMR (`M:1042`) | — | H (AMR decision) |
| `surface.refine_height` | `M:688` / 0.0 | as `refine_depth` | — | H |
| `surface_patch.jet_above_h` | `M:8182` / 0 | ARCHIVE L12600 D2v; L18980 "set on M2 only" | — (only `face_heat` R_above) | B3 |
| `surface_patch.jet_T_above` | `M:8187` / −1.0 | as `jet_above_h`; pair both-or-neither (`M:8189-8196`) | — | B3 |
| `surface_patch.jet_T_rec_fixed` | `M:8110` / REQUIRED if rec_fixed | ARCHIVE L8974 D2e diagnostic; C.7.3 | `jet_d2c` KT/K15 (run, not hashed) | B3 |
| `surface_patch.pad_shield_open` | `M:8004` / −1.0 | ARCHIVE L13601 D2aa "Wall lever measured and exhausted"; C.7.3 "P1 confirmed, P2 refuted" | — (only `face_heat` R_pad) | B3 |
| `weibull.size_volume` | `M:9322` / cell | ARCHIVE L19170, L19252 P2-2d "unvalidated, per-ladder anchor, no target"; refuses `flux` | — | B3 |
| `weibull.V_ref` | `M:9329` / −1.0 | as `size_volume`; set by no input | — | B3 |

[P2 header] also treats `melt_aware`, `robin_form = face` (as a control) and `jet_T_ent_mode =
rec_fixed` as CLOSED; [P1] classes `melt_aware` as BEAM (see C3).

### 3.9 DEAD (6), with evidence and the production lines to remove

| key | parse | evidence | production lines | other setters |
|---|---|---|---|---|
| `microstructure.phase<i>.alpha_attenuation` | `M:347-348`, **REQUIRED** | zero readers in the three files; the beam reads `beam.alpha_expr`; MODEL_CARD §12 item 11 "never read" | `rock_grimsel.inputs:38` (`= 2000.0`) | 29 unit dirs, 26 validation files, every gate input (top-level: silently inert once the parse goes) |
| `microstructure.phase<i>.grimsel.T_ref` | `Grimsel.H:115` (strict) | Grimsel's methods never read it; `Table::Build` (`Table.H:57`) receives a copy overwritten at `M:7060` | none | none |
| `microstructure.phase<i>.name` | `M:338` | `grep '\.name\b'` finds only the parse | `rock_grimsel.inputs:22` (`= grimsel`) | 29 unit dirs, 26 validation files (top-level) |
| `surface_patch.face` | `M:8014` (strict) | consumed, checked `== zhi`, no reader | `grimsel_calibration/input:95`; HOS_R argv `surface_patch.face=zhi` | 30 validation + 12 unit decks: **must be edited in the same packet** |
| `surface_patch.T_ext` | `M:8041` (strict) | member `M:7822`, zero readers | none | 27 validation (Hu, Kant onset, p0, Rossi) + 4 unit decks: edit in the same packet |
| `surface_patch.type` | `M:7940` (strict) | compared only under `prescribed_T` (`M:7952-7955`) | none | 17 Hu decks: edit in the same packet |

**Other production-deck lines that go when their parse goes** (read from the decks for this document):

| line | key | packet |
|---|---|---|
| `rock_grimsel.inputs:39` | `phase0.A_D = 0.0` (REQUIRED, DAMAGE) | C1a |
| `grimsel_calibration/input:94`, HOS_R argv | `surface_patch.pinned_idle_cycles = 1.0e9` (strict) | C5 |
| `hoeser/input_lab:105`, `grimsel_calibration/input:81` | `conductivity.kappa_damage_alpha = 0.0` | C1a |
| `model_form_flame.inputs:51` | `damage.enabled = 0` | C1a |
| `model_form_flame.inputs:20` | `spall.sample = top_cell` | C1a |
| `model_form_flame.inputs:19`, `:41` | `spall.surface_normal = 0`, `spall.flake_coherence_length = 0.0` | C2 |
| `model_form_flame.inputs:22` | `spall.melt_aware = 0` | C3 |
| `model_form_flame.inputs:27`, `:39-40` | `weibull.size_volume = cell`, `spall.flake_length = twocell`, `flake_fraction = cf` | B3 |
| `hoeser/input_lab:156`, `grimsel_calibration/input:58` | `solver = explicit` | D2 |
| `hoeser/input_lab:95`, `grimsel_calibration/input:57` | `alpha = 0.0 0.0` | D2 |
| `model_form_flame.inputs:32` | `spall.removal_energy = flake` | H |
| `model_form_flame.inputs:50` | `bit.enabled = 0` | H (per Q2) |

Keep `hc.heat.alpha` (`hoeser/input_lab:133`, `grimsel_calibration/input:67`): REQUIRED by
`HeatConduction.H:79` whatever happens to `M:153`.

### 3.10 DIAG (7) and RETIRED-ABORT (22)

- DIAG: `energy_ledger.enabled`, `spall.removal_events_csv`, `surface_patch.jet_profile_interval`,
  `side_power_interval` (parse-scope defect, B1), `spall.seed_void` (fixture tool, kept) are KEPT;
  `spall.clusters_csv` and `spall.h_col_events_csv` die in H.
- RETIRED-ABORT: 21 P2-10 keys (`M:97-118`, abort `M:119-122`) and `spallation.sp.constraint_depth`
  (P2-16, `M:127-130`). The 16 `surface_patch.*` entries would abort as "unused key" anyway (strict
  scope); the list is load-bearing only for `spall.*` [P1 Surprise 10]. They stay.

---

## 4. Open decisions for the user

**Q1. The implicit solver: port it or re-run Zhang explicit?**
- *Evidence:* `AdvanceMaterialImplicit` (`M:1642-1771`) asserts `finest_level == 0` (`M:1645`),
  hard-codes Neumann on every face (`M:1739-1741`) and ignores voids [P2 R1]. Exactly one deck uses it
  (`zhang_oglesby`, `solver = implicit` ×1) [P4 §1.2]. With a microstructure, `solver = implicit` is
  silently ignored (`M:1457`) [P2 R1, R4]. It has no microstructure equivalent; it is "a decision, not
  a mechanical merge" [P2 R1].
- *Recommendation:* re-run `zhang_oglesby` explicit on the merged path and delete the implicit solver
  in D2. Port only if the explicit run's wall time is prohibitive. Measure that in D1, not before.

**Q2. The bit (regime 3): port or retire?**
- *Evidence:* `meier/input_drilling` and `input_drilling_25mm` need `bit.enabled`. The bit lives only in
  `R:1894-1931` and the B26 priority; the face list refuses it (`M:1034`); it is not on the port list
  [P3 B.2, B.3; P4 §2]. Those decks are the base of most D2* and P1* studies [P4 §2].
- *Recommendation:* retire with H. The burner-body mechanical removal that Meier's geometry needs is
  the feet clearance (G), which is ported. The drawworks bit is a separate Step 19 device that no
  current prediction uses. The drilling decks become records.

**Q3. Rossi: port a prescribed-T surface onto the face list, or retire?**
- *Evidence:* `rossi/sp_rossi_damage_profile` is LEFM + Weibull despite its name, uses PATCH-T and
  scores `h_col_events.csv`. The face list refuses PATCH-T (`fl_on` needs `!prescribed_T`,
  `M:1883-1885`) and writes no h_col CSV [P4 §2; P3 A.4 row 3].
- *Recommendation:* port PATCH-T onto the list only if Q5 re-expresses Hu onset on LEFM, because the
  same port serves both (`FrFaceFlux` with T_s := the patch T). Otherwise retire Rossi with H.

**Q4. The pinned robin form: keep or retire?**
- *Evidence:* the face list requires `surface` (`M:999`); production moved to `surface` at P0b
  (ARCHIVE L14061) [P1 FLAME-CAP row]. Pinned is about 135 lines [P2 R7] and dominates the column gate
  (robin_pinned, robin_jet, jet_enthalpy quarter, robin_feet) [P4 §3.2]. The Meier feet decks are
  pinned [P4 §2]. `pin_clock` tests it [P4 §1].
- *Recommendation:* retire it with H (C5). Re-key the Meier feet decks to `surface` + `h_expr` in G.
  Its `pinned_idle_cycles = 1.0e9` lines on the production decks go in the same packet.

**Q5. Hu: re-express on LEFM, or retire?**
- *Evidence:* `hu_spall_onset` (+ `diag/`) uses a D-threshold onset proxy; `hu_end_to_end` uses
  damage-law removal on the column path [P4 §2]. An LEFM onset would be Sp_field under PATCH-T with
  LTE and removal off, a mechanism the Kant onset decks already prove [P4 §2]. Sp_field is criterion
  site A, which H would otherwise retire (Q8).
- *Recommendation:* re-express `hu_spall_onset` as an LEFM onset (it is a Hu number the model has
  never been scored against on its real criterion). Retire `hu_end_to_end`, which validates a retired
  removal model, and the `diag/` records. This ties Q5 to Q8 (site A kept as the onset criterion) and to
  Q3 (PATCH-T on the list).

**Q6. `melt_aware`: closed, or the softening hook for the "ceiling" step?**
- *What the existing code does* [P3 C.2; P2 R8]: f(Λ_L) = clamp(1 − Λ_L/λ_crit) (`M:8611-8618`)
  multiplies both branches of the scan's driving stress (`ScanSigmaT`, `R:138, 143`). It can only
  **stall** spallation under a molten cap; it adds or removes no energy. Λ_L comes from a column-profile
  gather with its own `ReduceRealSum` (`R:1017-1057`). Site A never applied it, so the two criterion
  sites disagree under the key [P2 R8]. Its test `sp_melt_skeleton` already FAILs (rc 1) [P4 §1]. The
  face list never runs it (`M:1037`).
- *Is it a usable hook?* Only the multiplier slot (one parameter in `ScanSigmaT`) is reusable. The
  part that costs, the Λ_L gather, is column-only and must be rebuilt on lines anyway. On the face list,
  Λ_L per line cell would come from the gathered H (`F:400-465`, column `G[o+2]`) through the phase
  table. The skeleton also keys on liquid fraction only; softening below the melt is already in the
  Grimsel law's T dependence (`alpha_cap_theta`, E(T)).
- *Recommendation:* treat it as closed and delete it in C3. If the ceiling step needs softening, build
  it in that packet from the line's H, with its own evidence; re-adding a multiplier is a few lines.
  Measure C3-0 first in any case: if the skeleton fired in the s1 legs, that is evidence for the
  ceiling step.

**Q7. `surface_normal`'s criterion half: retire?**
- *Evidence:* the recession half is already superseded (D2z-1, ARCHIVE L13103; `sn_recession_superseded`)
  [P1 CLOSED]. The criterion half divides depth by cos θ at site A and in the column scan
  (`M:6189-6194, 6220, 6300-6304`; `R:1460-1463`; `DepthScanH`'s `sn_cth`/`sn_icth`) [P2 R8]. The face
  list measures along the normal through mu = n·e* and applies /mu to the normal lengths (`F:241, 505,
  535-536, 377-380`) [P3 A.3 item 4], which is the per-line equivalent. The key is a Message rather than
  an abort only because two gate decks set it (`M:5486-5490`) [P3 A.4].
- *Recommendation:* retire it in C2. The cost is the 233-file re-pin, which C2 pays once for both switches.

**Q8. (added) Criterion site A: keep as the no-removal onset criterion, or give the face list an
onset-only mode?**
- *Evidence:* without removal, site A is the only criterion [P2 R3]. The Kant onset decks (`sp_kant_onset` ×5), the
  pressure sweep, and a re-expressed Hu onset (Q5) run with removal off [P4 §2]. Site A's
  `sigma_at_depth` (`M:6219-6276`) repeats `ScanSigmaT` verbatim and its `K_Ic` switch
  (`M:6307-6316`) repeats `ScanKIcPa` [P2 R3]. On production it runs every step for plotted output only
  (finding 8).
- *Recommendation:* keep site A as the onset criterion, rewritten on the shared helpers that A3 moves,
  and run it only when removal is off. That removes the production waste and the second copy of the
  arithmetic, and keeps Kant/Hu onset alive. Retiring it outright breaks the six kept Kant onset/pressure-sweep decks.

**Note (not a question).** `face_heat = path` is KEPT: [P1] records that the user restored it in P2-12
after the P2-11b OK that [P2] cites. If that is wrong, its ~385 lines [P2 R7] are a B3 item.

---

## 5. Defects found that are not deletions

**5.1 Parse-order bug on `spall.flake_h_max`.**
- *Lines:* the check `flake_h_max < h_spall_max` at `M:5610-5612` reads `spall_h_spall_max` before
  `spall.h_spall_max` is parsed at `M:5752`, which sits after the early return at `M:5681`. It always
  compares with the member default 0.010 (`M:8570`) [P1 Surprise 2; P2 R4].
- *Effect:* `h_spall_max = 0.015, flake_h_max = 0.012` aborts wrongly; `h_spall_max = 0.001,
  flake_h_max = 0.005` passes. Production (0.0015 vs 0.015) passes either way.
- *Related trap:* `BuildSurfaceTColumns :2979, 3019` reads `h_spall_max` on onset decks with
  `spall.enabled = 0`, where the deck's value is never parsed [P1 report ii E].
- *Fix in:* B1.

**5.2 Kelvin enthalpy in `AdvanceConstAlpha`.**
- *Lines:* H = T − T_ref in kelvin (`ConvertTtoH :6680-6690`). The kernel adds Q [W/m³] and
  q_loss/dz to a_th·∇²T [K/s] with no 1/(ρc_p) (`M:1508`) [P2 R1].
- *Effect:* the kernel is consistent only if the input takes ρc_p = 1. `beam`, `equilibrium` and
  `stefan` run on it.
- *Fix in:* D1/D2. The path is deleted; the tests re-point to a one-phase microstructure with real
  ρ, c_p, and their references re-pin.

**5.3 Cavity attenuation in the Zhang kernels.**
- *Lines:* all three non-microstructure kernels pass a geometric beam path (`BeamSource(..., DX[2])`,
  z0 − z; `dz_abs = −1`). Solid cells below a cavity are attenuated through the void: the E1 leak
  survives on this path. `AdvanceMaterialImplicit` also ignores voids in conduction [P2 R1, R5].
- *Fix in:* D2 (the path is deleted). D1's merged Zhang uses `SolidFacePath` (`M:4327-4331`).
  `zhang_oglesby` has no removal, so the leak should not move its reference; D1 must confirm that.

**5.4 The vapour front double-spends beam power.**
- *Lines:* `M:2048-2066`, `M:4918-4948`, `M:4227, 4295` (beam into H); `R:1706-1709`, `R:1722`
  (the same power into recession); `R:2242-2257` (voided H frozen); `R:2075, 2268, 2431-2435` (no
  sink); `R:1725` (top-cell-only gate); `R:2201-2212` (deeper cells voided whatever their H)
  [P3 C.4].
- *Also:* the Gaussian waist under a uniform profile (`R:1703-1708`); regime 2/3 unbooked in the ledger
  [P3 C.3].
- *Fix in:* F. Do not port the rate law.

**5.5 The Sp_field site and its sweeps run every production step for diagnostics only.**
- *Lines:* `UpdateSpAfterMechanics` (`M:5992-6343`, called at `Advance:1446`), `UpdateSpClusters`
  (`M:6353-6557`), `BuildSurfaceTColumns` (`M:2866-3033`, called at `M:1434`), and three
  `GatherColumnSurface` sweeps (`M:2882`, `M:6074`, `M:6378`), each with an MPI reduce. A production
  step makes 5 top sweeps + 1 run sweep where 2 are needed [P2 R2, R3].
- *Effect:* the outputs are plotted and hashed only. `FaceRemoval.H` never reads `Sp_field_mf`
  [P2 R3]. `PatchDiagnostics` likewise reports a column-form flux the face list does not apply
  (`M:3122-3196`, `M:2083-2084`) [P2 R6].
- *Fix in:* H (with Q8: site A runs only with removal off). No interim gating, because the
  fields are hashed on W′.

**5.6 alpha(T) does not conserve beam energy.**
- *Lines:* `BeamSource` uses `a = alpha_f(T_cell)` (`M:4932`) both for the cell's own absorption and
  for the attenuation over the whole path above it (`M:4933-4943`). The per-cell integral telescopes to
  1 − e^{−a·s_N} only if every cell has the same `a` [P2 R5].
- *Not caught:* the closure check `P_dep ≤ P_inc(1 + 1e-9)` runs on the uniform profile only
  (`M:4859-4861`). The code itself says it "needs the cumulative optical depth … required before real
  alpha(T) MMW work" (`M:4723-4727`). The Gaussian profile adds an O((dz/z_R)²) error, because `Pi`
  uses the geometric depth (`M:4974-4975`).
- *Fix in:* E (cumulative optical depth Σ a_j·dz along the solid path; closure check on both profiles).

**5.7 Smaller defects.**

| defect | lines | fix in |
|---|---|---|
| `hc.heat.alpha` called inert but aborts every spall step via the discarded damage-law h_spall | `R:297-301` [P1 Surprise 3] | C1a |
| `solver = implicit` with a microstructure is silently ignored | `M:1457` [P2 R4] | D2 |
| beam + face list accepted; the beam deposits through column tops only, and partial cells (`fr_f < 1`) absorb as full | `M:2048-2061`, `M:4328-4331` [P2 R4, R5; P3 A.4 row 12] | E (refuse until ported) |
| `side_power_interval` parsed inside the enthalpy branch | `M:8206` [P1 Surprise 8] | B1 |
| `kic.K_Ic0`/`T_ref`/`T_melt` silently ignored under `grimsel_table`; cohesive `G_c`/`f_t` silently ignored when both deltas are given | `M:9221-9233`, `M:5293` [P1 report ii D, H] | C1a (cohesive); K_Ic: a value-check in B1 |
| k_eff bulk diagnostic counts phases rank-locally, no MPI reduction | `M:7368-7378` [P2 R6] | the MICRO study (inert for one phase) |
| `RefreshKappaEff` result overwritten in the same step | `F:929`, `R:2576` vs `M:2072` [P2 R6] | B4 |
| face-list CSV carries fixed placeholders (regime 1, Sp_top 0, d_ell_rel/delta/q_ratio 0) and reads `k_top` from the column path | `F:894-903`, `F:895` [P3 A.1] | A3 (k_top), H (header) |
| double microstructure refresh in `TimeStepBegin` and `TimeStepComplete` (AMR mechanics only) | `M:1176-1197` [P2 T1 #6-7] | MECH cleanup, any time |
| `SurfacePatch::Parse`: under `convective_flame`, `surface_patch.T_expr` is consumed and ignored | `M:8034` [P1 report ii I] | B1 (refuse) |

---

## Appendix: contrast with the October order (CODE_CENSUS D.2)

D.2 (2026-10-08) made "2,639 of 2,639 hashes match" the acceptance criterion of every packet. It ran
the work as a pure refactor: delete closed keys, remove dead code, build single copies, cut
`Removal.H` into stages, reorder members, split files. Three facts from this census change that:
1. The gate is blind to the production path (finding 1). An acceptance criterion of 2,639/2,639
   certifies nothing about the face list. A1 replaces it with 2,639 + W′.
2. The column path is now deleted, not refactored. D.2 items 5-6 (`ColumnState`, `ColumnCriterion`,
   cutting `Removal.H` at its stages) would spend effort on code H removes. Only the seven helpers and
   the prologue are worth moving (A3).
3. Several deletions are physics-visible on the gate (DAMAGE plot fields, `surface_normal`, possibly
   `melt_aware`). D.2's rule "no physics number, 0 of 2,639" holds only for B and D. C needs declared
   re-pins.

D.2 items 1-4 survive as A1, B1-B3. Item 7 (member order, `ValidateSettings`) and item 8 (file split)
come after H, when there is one removal path left to lay out.
