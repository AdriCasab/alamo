# Census part 2: MMWSpalling.H, function by function (2026-11-05)

Read-only map of `src/Integrator/MMWSpalling.H` as it stands on 2026-11-05: **9,375 lines,
6,292 code / 2,674 comment / 409 blank**. Re-anchored by function name against
`docs/project/CODE_CENSUS.md` (2026-10-08, 11,033 lines then). Nothing under `src/` or
`tests/` was edited. Every line number below was opened.

**How the numbers were made** (all in `work/`, re-runnable from the repo root):

| script | output | what it does |
|---|---|---|
| `p2_split.py` | `p2_split_raw.tsv` | brace-depth splitter (strings/comments stripped); every block opened at class-body depth, signature start, end, code/comment/blank, leading comment block. Boundaries spot-checked by reading (all 106 depth-2 blocks). |
| `p2_callers.py` | `p2_callers.tsv` | who calls each function, across `MMWSpalling.H`, `Removal.H`, `FaceRemoval.H` |
| `p2_aborts.py` | `p2_aborts.tsv` | all 185 `Util::Abort` sites with their guarding condition |
| `p2_ranges.py` | `p2_lineage_by_line.tsv` | assigns every one of the 9,375 lines a lineage and a CLOSED flag (function base + hand-checked sub-ranges); prints the budget of report 7 |
| `p2_table.py`, `p2_members.py` | `p2_table.md`, `p2_members.md` | the two tables below |

**What "production" means here.** The production configuration is now a file:
`tests/MMWSpalling/validation/shared/model_form_flame.inputs` + `rock_grimsel.inputs`
(both Meier and the Kant calibration `FILE =` them). It sets `surface_patch.face_list = 1`,
`robin_form = surface`, `spall.surface_T = reconstructed`, `removal_energy = flake`,
`stress_source = local_thermoelastic`, `stress_law = grimsel`, `kic.law = grimsel_table`,
`phase0.law = grimsel`, `weibull.enabled = 1`, `energy_ledger.enabled = 1`,
`spall.buckle_span = 0.012`, `surface.follow_mask = 1`, one phase in expression mode,
`melt_aware = 0`, `surface_normal = 0`, vapor/bit/damage off, and `max_level = 0`. The
face-list Parse block (`MMWSpalling.H:967-1112`) refuses most of the column path's options, so
"reached on production" below means "reached with that file".

**Lineage codes.** PROD = PROD-SHARED. The other nine are as in the brief. CLOSED is a separate
flag: a default-off key whose study closed per CODE_CENSUS C.7.3 (`surface_normal`,
`melt_aware`, `flake_length`/`flake_fraction`, `flake_coherence_length`, `pad_shield_open`,
`jet_above_h`/`jet_T_above`, `jet_T_ent_mode = rec_fixed`, `robin_form = face` as a control).
Two classification calls worth stating:
- `UpdateSpAfterMechanics` is COLUMN (the Sp_field criterion site) and only its
  `stress_mf` branch is MECH (`:6061-6062, :6125-6126, :6270-6276`).
- `face_heat = path` (`FaceListPath` and its helpers) is FLAME. It is a flame heat-map
  feature, and the P2-11b review recorded that the user OK'd deleting the path code.

---

## Table 1. Functions in file order

`lines` is the signature-to-closing-brace range. `+N` is the leading comment block attached
above it, which is counted in the lineage column. Lineage shows the main lineage, then the
other lineages it contains, by line count. Indented rows are lambdas or sub-blocks of the row
above.

| # | name | lines | count | code / comment / blank (+lead comment) | lineage | guard (runs only if) | one-line purpose |
|---|---|---|---|---|---|---|---|
| 1 | `MMWSpalling` | 79-82 | 4 | 4 / 0 / 0 (+0) | PROD | always | constructors; the ParmParse one calls Parse |
| 2 | `MMWSpalling` | 86-86 | 1 | 1 / 0 / 0 (+0) | PROD | always | ParmParse constructor -> Parse |
| 3 | `Parse` | 88-1126 | 1039 | 699 / 299 / 41 (+0) | PROD (+FLAME 151, COLUMN 111, MICRO 64, DAMAGE 47, BEAM 36, MECH 36, DIAG 32, ZHANG 31, DEAD 1); CLOSED 50 | always (queryclass) | every key: retired-key aborts, base parses, the Parse*Settings calls, microstructure/phase parse, 60 RegisterNewFab and 125 RegisterIntegratedVariable calls, 42 aborts, the cross-key abort matrix, face-list validation |
|  |   retired-key guards | 88-130 | 43 |  | PROD | pp.contains(key) | abort on any P2-10/P2-16 retired key |
|  |   microstructure block | 292-624 | 333 |  | PROD (+MICRO 56, COLUMN 54, DAMAGE 35, MECH 33, BEAM 6, DEAD 1) | microstructure.nphases | phase parse + field registrations (order = plotfile layout) |
|  |   cross-key matrix + thermo registrations | 626-941 | 316 |  | PROD (+FLAME 117, DIAG 32, COLUMN 31); CLOSED 26 | various | the late abort matrix and keyed thermo columns |
|  |   buckle_span block | 943-965 | 23 |  | PROD | spall.buckle_span>0 | P2-15 aborts + bk_* witnesses |
|  |   face-list block | 967-1112 | 146 |  | PROD (+FLAME 34, BEAM 1); CLOSED 2 | surface_patch.face_list=1 | path aborts, F-requirements, fr_abort refusals, fl_/fr_ registrations |
| 4 | `Initialize` | 1129-1158 | 30 | 16 / 13 / 1 (+0) | PROD | always | enthalpy table (material only), base inits, T_init copy, phi/seed void (removal), mask, microstructure, T->H |
| 5 | `StoreInitialTemperature` | 1166-1174 | 9 | 9 / 0 / 0 (+6) | PROD | microstructure (T_init_mf exists) | the one write of T_init_mf, the stress-free reference |
| 6 | `TimeStepBegin` | 1176-1186 | 11 | 8 / 3 / 0 (+0) | PROD | always; refresh only if microstructure && max_level>0 && el.type!=disable | base hooks + microstructure refresh for AMR mechanics |
| 7 | `TimeStepComplete` | 1188-1197 | 10 | 10 / 0 / 0 (+0) | PROD | always; same refresh guard | base hooks + the same refresh again (B15 double refresh) |
| 8 | `Regrid` | 1199-1217 | 19 | 13 / 6 / 0 (+0) | PROD | AMR regrid | removal repair (lev>0), mask, microstructure re-init, H->T, mechanics repair |
| 9 | `UpdateModel` | 1219-1399 | 181 | 152 / 25 / 4 (+0) | MECH | el.type != disable (Mechanics base skips it otherwise) | nodal elastic model from cell E/mu/free strain; Hu prescribed_T top-node override; homogeneous alpha[] branch |
| 10 | `struct SurfaceRecon` | 1408-1423 | 16 | 13 / 3 / 0 (+7) | PROD (+DEAD 4) | type | result of ReconFromX (T_surf, ell, h_f, T_bar, flags); preceded by a stranded, stale Advance header (1401-1404) |
| 11 | `Advance` | 1425-1470 | 46 | 24 / 22 / 0 (+0) | PROD | always | the step: surface T columns, property refresh, mechanics, damage, cohesive, Sp_field, clusters, removal, swap, one conduction kernel, H->T |
| 12 | `AdvanceConstAlpha` | 1472-1511 | 40 | 39 / 0 / 1 (+0) | ZHANG | no material block, no microstructure | constant-alpha Laplacian + beam + cell-average losses (Step 2-4 regressions) |
| 13 | `AdvanceMaterial` | 1513-1633 | 121 | 91 / 22 / 8 (+0) | ZHANG | material.type set, solver=explicit, no microstructure | explicit kappa(T) harmonic-face conduction, void-aware, prescribed_T/const flame patch, beam, losses |
| 14 | `AdvanceMaterialImplicit` | 1642-1771 | 130 | 99 / 20 / 11 (+7) | ZHANG | material.type set, solver=implicit | backward-Euler MLABecLaplacian, single level, Neumann-only, ignores voids |
| 15 | `AdvanceMicrostructure` | 1788-2085 | 298 | 199 / 83 / 16 (+15) | PROD (+BEAM 36, FLAME 7) | microstructure.nphases set (production (Meier/Kant face-list deck)) | THE production kernel: face conductivities, top-face SurfaceCellFlux, face-list FrFaceFlux, beam, H update; calls BuildFaceList/BuildFlameColumns/BuildPinEffective and the diagnostics |
|  |   lambda k_face | 1930-1966 | 37 |  | PROD | per cell | harmonic face conductivity, void zero, GB contact resistance |
| 16 | `BuildFlameColumns` | 2100-2236 | 137 | 89 / 48 / 0 (+13) | FLAME; CLOSED 50 | surface_patch.enabled && owner (microstructure) | per-column flame_h/Tg/s (constants or h_expr at stand-off from FeetDescent or prescribed feed); pad shield x2; calls FaceListHeat, FrAssignTs, JetEnthalpyMarch, FaceListTg |
| 17 | `FaceListTg` | 2243-2276 | 34 | 34 / 0 / 0 (+5) | PROD | face_list=1 (fl_ready) | copies column flame_Tg to each face (position); path: per-face Tg and the fl_s/fl_Tg fields |
| 18 | `FeetDescent` | 2302-2419 | 118 | 111 / 7 / 0 (+24) | FLAME | surface_patch.nozzle_descent=feet | burner-on-feet nozzle height from pad quantiles; arms foot_clear_k (foot_body_clearance) |
| 19 | `JetEnthalpyMarch` | 2470-2845 | 376 | 310 / 61 / 5 (+49) | FLAME; CLOSED 1 | surface_patch.jet_closure=enthalpy | serial radial wall-jet enthalpy march: stagnation state, per-bin draw, dilution, profile CSV, 1e-9 invariant abort |
|  |   lambda q_col | 2601-2629 | 29 |  | FLAME | per column | column Robin draw; under the face list the SurfaceCellFlux/PinRule part is discarded and the listed faces are billed |
|  |   lambda q_face | 2662-2668 | 7 |  | FLAME | face_heat=path | one face billed alone |
| 20 | `BuildSurfaceTColumns` | 2866-3033 | 168 | 107 / 55 / 6 (+19) | COLUMN (+DEAD 7) | spall.surface_T in {diagnostic, reconstructed} && owner | per-column two-cell exponential reconstruction (surf_T/ell/h_f/T_bar/flag) + flat-hold + flake_length=flux q_net length |
| 21 | `SharedSurfaceRecon` | 3046-3062 | 17 | 17 / 0 / 0 (+11) | COLUMN | callers: Sp_field site (recon_on), column removal path | rebuilds a SurfaceRecon from the per-column arrays |
| 22 | `SurfaceRobinT` | 3069-3074 | 6 | 6 / 0 / 0 (+5) | COLUMN (+DEAD 5) | returns 0 unless robin_form=surface | column-top reconstructed surface T for SurfaceCellFlux |
| 23 | `BuildPinEffective` | 3076-3093 | 18 | 17 / 1 / 0 (+0) | FLAME | robin_form=pinned && sp_on | pin_T_eff/pin_reason per column via PinRule |
| 24 | `PinRule` | 3098-3115 | 18 | 14 / 4 / 0 (+3) | FLAME | robin_form=pinned | the C1 idle-clock pin rule for one column |
| 25 | `PatchDiagnostics` | 3122-3196 | 75 | 75 / 0 / 0 (+5) | DIAG | owner && (pinned \|\| (sp_on && use_expr)) | pinned branch counters + patch_P_robin by the COLUMN form (not the flux the face list applies) |
| 26 | `struct FaceEntry` | 3261-3285 | 25 | 17 / 8 / 0 (+62) | PROD | type | one listed solid/void face; preceded by SurfaceCellFlux's stranded 37-line header (3198-3234) and the P2-8 banner |
| 27 | `FrPathOffsets` | 3331-3347 | 17 | 16 / 1 / 0 (+20) | FLAME | face_heat=path | Dijkstra neighbour offsets \|o\|<=4 |
| 28 | `FaceListPath` | 3349-3530 | 182 | 168 / 9 / 5 (+0) | FLAME | surface_patch.face_heat=path | P2-11b: surface-path distance from stagnation by Dijkstra, barrier, h at r:=s_path |
| 29 | `FaceCentre` | 3533-3546 | 14 | 14 / 0 / 0 (+1) | FLAME | callers only under face_heat=path | a face centre (FaceListHeat inlines the same arithmetic) |
| 30 | `FrGauss` | 3558-3564 | 7 | 7 / 0 / 0 (+10) | PROD | face_list=1 | Gaussian weight for the smoothed fraction |
| 31 | `FrSmoothAt` | 3565-3573 | 9 | 9 / 0 / 0 (+0) | PROD | face_list=1 | 5x5x5 smoothed solid fraction |
| 32 | `FrSmoothedGradient` | 3574-3582 | 9 | 9 / 0 / 0 (+0) | PROD | face_list=1 | central-difference gradient of the smoothed fraction (the normal) |
| 33 | `FrInitFraction` | 3585-3597 | 13 | 13 / 0 / 0 (+1) | PROD | face_list=1 (once) | solid fraction 1/0 from removed_mf |
| 34 | `FillFractionGhosts` | 3603-3649 | 47 | 47 / 0 / 0 (+4) | PROD | face_list=1 | fraction ghost rule: gas above, mirror at symmetry planes, clamp elsewhere |
| 35 | `FaceWeight` | 3652-3660 | 9 | 9 / 0 / 0 (+1) | PROD | face_list=1 | F3 area weight 1/(\|nx\|+\|ny\|+\|nz\|) |
| 36 | `FaceH` | 3663-3669 | 7 | 7 / 0 / 0 (+1) | PROD | face_list=1, face_heat=position | F2 h at a face centre (r from the patch centre) |
| 37 | `FaceHr` | 3673-3691 | 19 | 19 / 0 / 0 (+3) | PROD (+FLAME 6); CLOSED 6 | face_list=1 | h_expr at (x,y,r,s) with jet zeroing above the plane and the pad shield |
| 38 | `BuildFaceList` | 3694-3844 | 151 | 140 / 9 / 2 (+1) | PROD | face_list=1 && sp_on && removal && !prescribed_T | F1+F3: every solid/void face in the patch, gathered/sorted on every rank, per-cell index |
| 39 | `FaceListHeat` | 3848-3901 | 54 | 49 / 5 / 0 (+2) | PROD (+FLAME 7) | face_list=1 (inside BuildFlameColumns, use_expr) | F2 over the list (position) or FaceListPath (path); fl_h field; top-face mismatch witness |
| 40 | `FaceListDiagnostics` | 3906-3926 | 21 | 19 / 2 / 0 (+3) | DIAG | face_list=1 (fl_on) | per-class robin power fl_P_*, fl_id_march |
| 41 | `SurfaceCellFlux` | 3928-4045 | 118 | 74 / 44 / 0 (+0) | PROD (+FLAME 11) | every non-void owner cell (microstructure); FrFaceFlux on the face list | the one Robin/loss routine: surface, face-Newton, pinned and cell forms |
| 42 | `EnergyLedger` | 4064-4313 | 250 | 219 / 23 / 8 (+17) | DIAG | energy_ledger.enabled && lev==0 (on the production deck) | compensated H sums, Robin/loss re-evaluation (face-list mirror), side-power CSV, ledger_err/_solid |
|  |   lambda neumaier | 4092-4097 | 6 |  | DIAG | - | compensated sum |
| 43 | `SolidFacePath` | 4327-4331 | 5 | 5 / 0 / 0 (+12) | BEAM | beam on && void-aware | Beer-Lambert path (n_above+1/2)dz |
| 44 | `SurfaceOwner` | 4337-4337 | 1 | 1 / 0 / 0 (+4) | PROD | always | finest level owns the surface |
| 45 | `RepairRemovalStateAfterRegrid` | 4349-4436 | 88 | 71 / 14 / 3 (+10) | COLUMN | Regrid, lev>0 && removal (AMR; face list refuses max_level>0) | rebuild column-linear phi / clean removed after regrid |
| 46 | `CheckSurfaceCovered` | 4442-4477 | 36 | 36 / 0 / 0 (+4) | COLUMN | column removal path, max_level>0 | abort if the surface left the refined band |
| 47 | `UpdatePropertyFields` | 4512-4587 | 76 | 68 / 6 / 2 (+33) | PROD | some phase law=grimsel (production (Meier/Kant face-list deck)) | per-step Grimsel refresh of kappa, E, mu, rho, Cp, secant beta + prop_* witnesses; preceded by ColumnTopSolid's stranded header (4479-4484) |
| 48 | `ColumnTopSolid` | 4589-4611 | 23 | 23 / 0 / 0 (+0) | PROD | void-aware && (beam \|\| follow_mask \|\| ledger); BuildFaceList; mech clearance; tagging | top solid cell per column, ReduceIntMax |
| 49 | `ColumnSolidAboveCount` | 4627-4665 | 39 | 36 / 2 / 1 (+14) | BEAM | beam on && void-aware | solid cells above each (col,k), ncols*nz reduction |
| 50 | `ColumnSolidRunBelowTop` | 4677-4714 | 38 | 35 / 2 / 1 (+10) | COLUMN | BuildSurfaceTColumns; column removal path | contiguous solid run below each column top |
| 51 | `BeamEnergyBalance` | 4748-4887 | 140 | 101 / 33 / 6 (+32) | BEAM | beam power > 0 (microstructure kernel epilogue) | E1 beam closure P_inc/P_dep/P_trans, invariant (a), uniform double-count guard |
| 52 | `BeamSource` | 4918-4947 | 30 | 26 / 4 / 0 (+10) | BEAM | every kernel; returns 0 at P=0 | Beer-Lambert x Gaussian/uniform source, analytic per-cell integral |
| 53 | `BeamIncidentIntensity` | 4954-4980 | 27 | 24 / 3 / 0 (+5) | BEAM | beam on | incident intensity Pi(x,y,z) |
| 54 | `TagCellsForRefinement` | 4982-5010 | 29 | 24 / 5 / 0 (+0) | COLUMN (+PROD 6) | AMR; band only if removal && surface.refine_depth>0 | base tagging + surface band |
| 55 | `ParseDamageSettings` | 5013-5110 | 98 | 92 / 0 / 6 (+0) | DAMAGE | always (keys default off) | damage.* Drucker-Prager parse and checks |
| 56 | `RefreshKappaEff` | 5112-5128 | 17 | 17 / 0 / 0 (+0) | DAMAGE | init/regrid; also called by both removal paths and damage | kappa_eff = kappa exp(-a_D D) (plotted) |
| 57 | `DamageGrowthRHS` | 5130-5142 | 13 | 13 / 0 / 0 (+0) | DAMAGE | damage.enabled | damage-rate law |
| 58 | `IntegrateDamageRK4` | 5144-5160 | 17 | 17 / 0 / 0 (+0) | DAMAGE | damage.enabled | RK4 of the damage law |
| 59 | `UpdateDamageAfterMechanics` | 5162-5271 | 110 | 99 / 0 / 11 (+0) | DAMAGE | damage.enabled | Drucker-Prager damage update from stress_mf or prescribed stress |
| 60 | `ParseCohesiveSettings` | 5273-5324 | 52 | 45 / 4 / 3 (+0) | DAMAGE | always (cohesive.enabled default 0) | cohesive.* parse |
| 61 | `UpdateCohesiveAfterMechanics` | 5326-5389 | 64 | 56 / 0 / 8 (+0) | DAMAGE | cohesive.enabled | bilinear GB cohesive law (no feedback to the solve) |
| 62 | `ResetCohesiveDiagnostics` | 5391-5398 | 8 | 8 / 0 / 0 (+0) | DAMAGE | cohesive.enabled | zero gb_* fields |
| 63 | `ParseSpallSettings` | 5400-5811 | 412 | 217 / 172 / 23 (+0) | PROD (+COLUMN 77, DAMAGE 65, DIAG 58, BEAM 25, DEAD 11); CLOSED 102 | always; most keys after `if (!spall_enabled) return` (:5681) | spall.* keys: stress source/law, surface_normal, melt_aware, CSVs, surface_T, flake keys, removal_energy, damage-law keys, h_spall_max, coherence, seed_void |
| 64 | `ParseVaporSettings` | 5813-5847 | 35 | 33 / 0 / 2 (+0) | BEAM | always (vapor.enabled default 0) | vapor.* parse |
| 65 | `ParseBitSettings` | 5854-5877 | 24 | 22 / 0 / 2 (+5) | COLUMN | always (bit.enabled default 0) | bit.* parse |
| 66 | `PrincipalStressRatio` | 5883-5897 | 15 | 15 / 0 / 0 (+4) | DAMAGE | damage-law spall (column path) | sigma_1/sigma_3 ratio |
| 67 | `GatherColumnSurface` | 5907-5940 | 34 | 33 / 1 / 0 (+8) | COLUMN | BuildSurfaceTColumns, Sp_field site, UpdateSpClusters | top solid cell per column (second copy of ColumnTopSolid; domain top if no removal); binds phi unused |
| 68 | `ColumnSurfaceCos` | 5958-5990 | 33 | 33 / 0 / 0 (+16) | COLUMN; CLOSED 49 | spall.surface_normal=1 | per-column wall-slope cos(theta); preceded by UpdateSpAfterMechanics's stranded header (5942-5949) |
| 69 | `UpdateSpAfterMechanics` | 5992-6343 | 352 | 225 / 105 / 22 (+0) | COLUMN (+MECH 11); CLOSED 18 | spallation.model=sp_weibull && owner (runs on the production deck, output plotted only) | criterion site A: Sp = K_I/K_Ic at each column top into Sp_field; recon witnesses |
|  |   lambda T_recon_at | 6181-6188 | 8 |  | COLUMN; CLOSED 1 | recon | lerp in the top half-cell |
|  |   lambda T0_at | 6214-6218 | 5 |  | COLUMN | - | T_init with T_ref fallback |
|  |   lambda sigma_at_depth | 6219-6274 | 56 |  | COLUMN (+MECH 5) | - | prescribed / local thermoelastic (law or scalar) / stress_mf stress at depth |
| 70 | `UpdateSpClusters` | 6353-6557 | 205 | 163 / 31 / 11 (+8) | COLUMN | sp_weibull && weibull.enabled && owner (runs on the production deck) | top-plane union-find of Sp>=1 into Sp_cluster_id; A_crit gate; clusters CSV |
| 71 | `LawSigma` | 6591-6595 | 5 | 5 / 0 / 0 (+32) | PROD | spall.stress_law=grimsel | the one Grimsel thermoelastic stress entry |
| 72 | `RockLaw` | 6604-6616 | 13 | 5 / 8 / 0 (+7) | PROD | spall.stress_law=grimsel | phase-0 law pointer or nullptr |
| 73 | `BuildEnthalpyTable` | 6626-6630 | 5 | 5 / 0 / 0 (+2) | ZHANG | material block | global-material H<->T table |
| 74 | `SetPhaseFractionsConstAlpha` | 6635-6642 | 8 | 8 / 0 / 0 (+3) | ZHANG | no material, no microstructure | closed-form phase fractions |
| 75 | `ConvertTtoH` | 6644-6696 | 53 | 49 / 4 / 0 (+0) | PROD (+ZHANG 15) | Initialize | T->H by phase table / global table / closed form |
| 76 | `InvertHtoT` | 6698-6760 | 63 | 60 / 3 / 0 (+0) | ZHANG (+PROD 29) | every step; regrid; removal | H->T and phase fractions, same three-way dispatch; production is only the per-phase branch 6716-6723, the rest is the global-material/const-alpha path |
| 77 | `InitSurfaceMask` | 6762-6778 | 17 | 17 / 0 / 0 (+0) | PROD | init/regrid/refresh | surface_mf: domain top, or phi window under removal |
| 78 | `UpdateSurfaceMaskFromPhi` | 6793-6806 | 14 | 14 / 0 / 0 (+13) | COLUMN | removal (init + both removal paths) | surface_mf = phi window; read only when follow_mask=0 and by the ledger's missing-column count |
| 79 | `ApplySeedVoid` | 6816-6844 | 29 | 29 / 0 / 0 (+8) | DIAG | spall.seed_void given (debug) | void seeded index boxes once |
| 80 | `InitPhi` | 6849-6882 | 34 | 34 / 0 / 0 (+3) | PROD | removal | flat level set, removed=0, reset per-mode fields |
| 81 | `RefreshMicrostructureDerivedFields` | 6884-6891 | 8 | 8 / 0 / 0 (+0) | MECH | microstructure && max_level>0 && el.type!=disable | mask + InitializeMicrostructure on every level |
| 82 | `RepairMechanicsAfterRegrid` | 6893-6915 | 23 | 19 / 0 / 4 (+0) | MECH | el.type != disable | refresh model and zero stress/strain after regrid |
| 83 | `FiniteMatrix` | 6917-6923 | 7 | 7 / 0 / 0 (+0) | DAMAGE | damage.enabled | finite-entries check |
| 84 | `AlphaBetaFraction` | 6926-6936 | 11 | 11 / 0 / 0 (+0) | MECH | mechanics / free strain | smoothstep alpha-beta fraction |
| 85 | `ThermalFreeStrain` | 6938-6949 | 12 | 12 / 0 / 0 (+0) | MECH | mechanics / free strain | beta(T-Tref) + alpha-beta strain |
| 86 | `CellToNodeAverageInDomain` | 6951-6981 | 31 | 31 / 0 / 0 (+0) | MECH | UpdateModel | in-domain cell-to-node average |
| 87 | `UpdateMicrostructureFreeStrain` | 6983-7006 | 24 | 24 / 0 / 0 (+0) | MECH | UpdateModel; InitializeMicrostructure | fills thermal_free_strain_mf (plotted) |
| 88 | `InitializeMicrostructure` | 7013-7404 | 392 | 289 / 92 / 11 (+5) | MICRO (+PROD 150, DAMAGE 28) | microstructure (Initialize, Regrid, AMR mechanics refresh) | phase tables, Voronoi/expression phase+grain ids, property fill, ghost fill, GB flags, damage/cohesive reset, k_eff, flaws |
| 89 | `InitializeFlaws` | 7414-7546 | 133 | 91 / 36 / 6 (+8) | MICRO | weibull.enabled | hashed per-cell Weibull flaw length |
| 90 | `enum class SolverKind` | 7575-7575 | 1 | 1 / 0 / 0 (+0) | ZHANG | type | explicit \| implicit |
| 91 | `WeibullSizeFactor` | 7752-7759 | 8 | 6 / 2 / 0 (+2) | COLUMN | weibull.size_volume=flux (face list refuses it) | (V/V_ref)^(1/m) size factor |
| 92 | `struct SurfacePatch` | 7786-8357 | 572 | 426 / 131 / 15 (+12) | FLAME (+PROD 154); CLOSED 95 | surface_patch.* (queryclass) | the patch settings (~70 members) and their 460-line Parse: mode, robin_form, face_list/face_heat, pad shield, jet closure + D2b-D2v options, feet, h_expr |
|  |   CheckPadRing | 7887-7895 | 9 |  | FLAME; CLOSED 9 | pad_shield_open>=0 | pad ring must exist (called on both Parse exits) |
|  |   SurfacePatch::Parse | 7897-8356 | 460 |  | FLAME (+PROD 109); CLOSED 53 | surface_patch.enabled | every surface_patch.* key, ~50 aborts |
| 93 | `struct MineralPhase` | 8365-8402 | 38 | 27 / 11 / 0 (+5) | PROD (+DEAD 2) | type | per-phase scalars, Grimsel law, alpha-beta; alpha_atten parsed, never read |
| 94 | `enum class MicroMode` | 8403-8403 | 1 | 1 / 0 / 0 (+0) | MICRO | type | voronoi \| expression |
| 95 | `enum class SpallSampleMode` | 8525-8525 | 1 | 1 / 0 / 0 (+7) | DAMAGE | type | top_cell \| column_max \| shallow_max (damage law) |
| 96 | `MeltSkeleton` | 8611-8618 | 8 | 8 / 0 / 0 (+5) | BEAM; CLOSED 13 | spall.melt_aware=1 (column path only) | f(Lambda_L) = clamp(1 - LL/lambda_crit) |
| 97 | `removal_enabled` | 8686-8686 | 1 | 1 / 0 / 0 (+0) | PROD | - | spall \|\| vapor |
| 98 | `enum class SpallationModel` | 8692-8692 | 1 | 1 / 0 / 0 (+4) | PROD | type | damage_law \| sp_weibull |
| 99 | `enum class KIcLaw` | 8715-8715 | 1 | 1 / 0 / 0 (+9) | PROD | type | nasseri \| linear \| grimsel_table |
| 100 | `enum class DetachMode` | 8888-8888 | 1 | 1 / 0 / 0 (+7) | COLUMN | type | per_face \| connected_cluster |
| 101 | `enum class SurfaceTMode` | 8964-8964 | 1 | 1 / 0 / 0 (+41) | PROD | type | cell_average \| diagnostic \| reconstructed (the stored mode is never read) |
| 102 | `RawSurfaceX` | 9000-9010 | 11 | 11 / 0 / 0 (+21) | PROD | BuildSurfaceTColumns; FrAssignTs; RecedeAlongLines | two-cell x = ln((T_top-T0)/(T_1-T0)) |
| 103 | `ReconFromX` | 9020-9109 | 90 | 31 / 57 / 2 (+8) | PROD | same callers | the reconstruction algebra: T_surf, ell, h_f (c_f, h_spall_max cap, flake_h_max bound), T_bar |
| 104 | `SurfaceGasTempCap` | 9117-9124 | 8 | 8 / 0 / 0 (+6) | PROD | reconstruction callers | gas-temperature ceiling max(T_flame, jet_T_nozzle) |
| 105 | `SplitMix64` | 9129-9135 | 7 | 7 / 0 / 0 (+3) | MICRO | weibull.enabled | hash for the flaw RNG |
| 106 | `ParseSpallationModel` | 9137-9372 | 236 | 157 / 62 / 17 (+0) | PROD (+COLUMN 73, MICRO 56, DAMAGE 10) | always; sp_* after `if model != sp_weibull return` (:9164) | model, confining.p, buckle_span, sp.*, kic.*, prescribed sigma/T, weibull.* |

## Table 2. Non-function blocks (members, includes), runs of 3 lines or more

These are the lines between functions: data members, nested enums, the two `#include`s of
`Removal.H`/`FaceRemoval.H` (`:6621-6622`), and the file header. With Table 1 they cover
all 9,375 lines (`work/lineage_by_line.tsv` has the per-line assignment). Members are still
interleaved with functions (B10): face-list members sit at `:3286-3309`, between
`FaceEntry` and `FrPathOffsets`.

| lines | count | lineage | code/comment |
|---|---|---|---|
| 1-78 | 78 | PROD | 41/30 |
| 83-85 | 3 | PROD | 1/0 |
| 3286-3302 | 17 | PROD | 13/4 |
| 3303-3309 | 7 | FLAME | 3/4 |
| 4888-4907 | 20 | BEAM | 1/18 |
| 6617-6623 | 7 | PROD | 2/3 |
| 7547-7566 | 20 | PROD | 15/2 |
| 7567-7580 | 14 | ZHANG | 8/3 |
| 7583-7598 | 16 | BEAM | 10/4 |
| 7599-7607 | 9 | PROD | 3/5 |
| 7608-7613 | 6 | DIAG | 1/4 |
| 7614-7681 | 68 | FLAME | 40/28 |
| 7682-7714 | 33 | PROD | 8/25 |
| 7715-7749 | 35 | COLUMN | 2/32 |
| 7762-7773 | 12 | DIAG | 10/0 |
| 8403-8414 | 12 | MICRO | 11/1 |
| 8415-8418 | 4 | PROD | 2/2 |
| 8419-8434 | 16 | MICRO | 8/7 |
| 8435-8442 | 8 | PROD | 3/5 |
| 8443-8449 | 7 | MECH | 4/2 |
| 8450-8501 | 52 | DAMAGE | 36/12 |
| 8502-8506 | 5 | PROD | 0/5 |
| 8507-8528 | 22 | DAMAGE | 14/7 |
| 8529-8537 | 9 | COLUMN CLOSED | 1/7 |
| 8538-8550 | 13 | COLUMN | 0/13 |
| 8551-8557 | 7 | COLUMN CLOSED | 1/6 |
| 8558-8562 | 5 | DIAG | 1/3 |
| 8563-8592 | 30 | PROD | 4/25 |
| 8593-8598 | 6 | COLUMN CLOSED | 1/5 |
| 8599-8605 | 7 | BEAM CLOSED | 2/5 |
| 8619-8622 | 4 | COLUMN CLOSED | 1/3 |
| 8623-8638 | 16 | COLUMN | 3/12 |
| 8639-8643 | 5 | PROD | 4/0 |
| 8644-8659 | 16 | BEAM | 8/7 |
| 8660-8675 | 16 | COLUMN | 9/7 |
| 8676-8678 | 3 | BEAM | 3/0 |
| 8679-8685 | 7 | PROD | 1/5 |
| 8687-8721 | 35 | PROD | 11/19 |
| 8722-8758 | 37 | COLUMN | 17/18 |
| 8759-8763 | 5 | PROD | 2/3 |
| 8764-8774 | 11 | DEAD | 0/11 |
| 8775-8795 | 21 | COLUMN CLOSED | 5/16 |
| 8797-8827 | 31 | DIAG | 10/20 |
| 8828-8858 | 31 | COLUMN | 9/21 |
| 8859-8880 | 22 | MICRO | 7/13 |
| 8881-8914 | 34 | COLUMN | 9/23 |
| 8915-8922 | 8 | DIAG | 2/5 |
| 8923-8978 | 56 | PROD | 5/49 |
| 9373-9375 | 3 | PROD | 3/0 |

---

## Report 1. The four conduction kernels (B6)

`Advance` dispatches at `:1457-1460`: if `microstructure_enabled`, `AdvanceMicrostructure`;
else if `material` and `solver = implicit`, `AdvanceMaterialImplicit`; else if `material`,
`AdvanceMaterial`; else `AdvanceConstAlpha`. A run never mixes kernels. `Initialize:1136`
aborts on `material` together with `microstructure`.

| | `AdvanceConstAlpha` 1472-1511 | `AdvanceMaterial` 1513-1633 | `AdvanceMaterialImplicit` 1642-1771 | `AdvanceMicrostructure` 1788-2085 |
|---|---|---|---|---|
| runs when | no `material`, no `microstructure` | `material.type`, `solver = explicit` | `material.type`, `solver = implicit` | `microstructure.nphases` (**production**) |
| conductivity | constant `hc.heat.alpha` x `Numeric::Laplacian` (`:1496`) | `mat->kappa(T)` per neighbour, harmonic face (`:1582-1598`) | `mat->kappa(T_old)` harmonic B-coefficients (`:1709-1723`) | per-cell `kappa_mf` (Grimsel-refreshed by `UpdatePropertyFields`) x `exp(-a_D D)`, harmonic, GB contact resistance in series (`k_face` `:1930-1966`) |
| enthalpy units / inversion | H in kelvin (`H = T - T_ref`, `ConvertTtoH:6680-6690`). Adds `Q` [W/m^3] and `q_loss/dz` to `a_th lap T` [K/s] with no 1/(rho c_p) (`:1508`), so it is consistent only if the input takes rho c_p = 1 | global `material_table` | BE on T, then `H = material_table.H_of_T(T_new)` (`:1768`) | per-phase `phase_tables` (Grimsel `Table::Build`, `:7051-7066`) |
| voids | none | void-aware (`removed` ghosts, `:1537-1538`, `nvoid` `:1585-1591`) | ignores voids | void-aware |
| surface patch | none | prescribed_T or constant h on the **domain** top `k = bigEnd` (`:1611-1619`) | none (Neumann hard-coded on every face, `:1739-1741`) | `SurfaceCellFlux` at the mask column top (follow_mask) in four forms, plus face-list `FrFaceFlux` (`:2003-2033`) |
| losses | inline on `surface_mf` (`:1504-1507`) | inline on `surface_mf` (`:1625-1628`) | inline, explicit (`:1700-1704`) | via `SurfaceCellFlux` (cell form `:4035-4044` for out-of-patch tops) |
| beam path | `BeamSource(..., DX[2])`, geometric path z0 - z | same, geometric: **solid cells below a cavity are attenuated through the void** (the E1 leak survives on this path) | same, geometric | `SolidFacePath` solid count (`:2048-2061`) |
| AMR / level | any | any | asserts `finest_level == 0` (`:1645`) | owner-gated surface terms (`:1825`) |

**Shared copy-paste.** The beam-and-loss local block appears in all four kernels: `:1478-1485`,
`:1519-1526`, `:1651-1658` and `:1845-1852`.
The inline loss expression appears three times: `:1506`, `:1627` and `:1702`. The fourth
kernel uses `SurfaceCellFlux` instead, and `SurfaceCellFlux` repeats the expression in its own
three forms (`:3979, :4027, :4041`). `kc = kappa exp(-a_D D)` appears 8 times: `:1923`,
`:1943`, `:2500`, `:3152`, `:3757`, `:4149`, `:5125` and `FaceRemoval.H:454`. The
`a_D == 0` exp-skip exists only at `:1923/:1943/:3757`. The march (`:2500`), PatchDiagnostics
(`:3152`) and the ledger (`:4149`) call `exp` unconditionally.

**Production** runs `AdvanceMicrostructure` and nothing else. `solver = implicit` combined with
a microstructure is **silently ignored**: the dispatch checks `microstructure_enabled` first,
and nothing aborts on it (see report 4). `AdvanceConstAlpha`, `AdvanceMaterial` and
`AdvanceMaterialImplicit` are reached only without a microstructure block, which means the
Step 2-4 regressions (const-alpha) and the Zhang validation (`material.type`, explicit or
implicit). Hu runs on the microstructure kernel through its `prescribed_T` branch (`:1996`).

**After a ZHANG merge onto the phase-law path** (Zhang as `microstructure.phase0.law = zhang`,
a single-phase expression microstructure, its `Material` going through the same polymorphic
`Table::Build` the Grimsel law uses at `:7051-7066`, and `UpdatePropertyFields` generalised
from "grimsel" to "any T-dependent law"), these become redundant:
- `AdvanceMaterial` (121 lines).
- `AdvanceConstAlpha` (40). The Step 2-4 regressions move to a constant one-phase
  microstructure. Note that their H is in kelvin today, so the references re-pin.
- `BuildEnthalpyTable` and `SetPhaseFractionsConstAlpha`.
- The global-material and closed-form branches of `ConvertTtoH` (`:6676-6690`) and
  `InvertHtoT` (`:6724-6757`).
- The `material` + `microstructure` abort (`:1136`).
- Members `material`, `material_table`, `SolverKind`, `implicit_*`, and `T_ref/T_m/T_v/lh_*`
  (the const-alpha `phase.*` keys, `:155-161`).
- Parse `:174-181` and `:252-267`.
- `UpdateModel`'s homogeneous `alpha[]/eta` branch (`:1361-1398`, MECH, reachable only
  without a microstructure).

`AdvanceMaterialImplicit` has no microstructure equivalent: it is single-level, Neumann-only and
void-blind. Either it is ported (an implicit microstructure kernel) or the Zhang implicit leg is
re-run explicit. That is a decision, not a mechanical merge.

---

## Report 2. The surface-temperature machinery (B3)

Every place a column-top or face surface temperature is formed, and whether production
(face list) reaches it:

| site | lines | what | production (face list) |
|---|---|---|---|
| `RawSurfaceX` | 9000-9010 | two-cell x = ln((T_top-T0)/(T_1-T0)) | **yes**: called by `FrAssignTs` (`FaceRemoval.H:235`) and `RecedeAlongLines` (`FaceRemoval.H:501`) per surface cell along its line |
| `ReconFromX` | 9020-9109 | T_surf, ell_eff, h_f (c_f, `h_spall_max` cap, `flake_h_max` bound), T_bar | **yes**: same two callers (`FaceRemoval.H:240, :504`) |
| `SurfaceGasTempCap` | 9117-9124 | ceiling max(T_flame, jet_T_nozzle) | **yes** (`FaceRemoval.H:229, :469`) |
| `FrAssignTs` (FaceRemoval.H) | 162-253 | per listed cell T_s (carry from removal, or fresh recon, or fallback flag) | **yes**, called from `BuildFlameColumns:2232` |
| `FrFaceFlux` (FaceRemoval.H) | 138-156 | every listed face's flux at its cell's T_s; fallback = face Newton on e* | **yes**: kernel `:2030`, march `:2624/:2666`, ledger `:4187`, `FaceListDiagnostics:3919` |
| `SurfaceCellFlux` surface form | 3974-3986 | flux + losses at a supplied T_surf | **yes**, through `FrFaceFlux` (`T_surf_in = Ts`) |
| `SurfaceCellFlux` face form (Newton) | 3987-4022 | infer T_f on a linear half-cell | **yes**, only as `FrFaceFlux`'s fallback (`FaceRemoval.H:148`) and `RecedeAlongLines:523` |
| `SurfaceCellFlux` pinned branch | 4023-4033 | T_s = T_pin | no (face list requires `robin_form = surface`, `:999`) |
| `SurfaceCellFlux` cell form | 4035-4044 | flux on the cell average + losses on `surf_c` | **yes, for out-of-patch column tops only** (in-patch tops are all listed and their result is zeroed at `:2020`) |
| `FaceH` / `FaceHr` / `FaceListHeat` | 3663-3691, 3848-3901 | per-face h (the heat rule; no temperature) | **yes** |
| `BuildSurfaceTColumns` | 2866-3033 | per-column recon into `surf_T_col` etc. via `GatherColumnSurface` + `ColumnSolidRunBelowTop` | **runs** (`:1434`, keyed on `surface_T`), but under the face list its output only feeds diagnostics: the Sp_field site, the ledger's `surf_q_mean` (`:4287-4292`), and `SurfaceRobinT`, whose consumers discard it |
| `SharedSurfaceRecon` | 3046-3062 | rebuild a SurfaceRecon from the column arrays | Sp_field site (`:6159`) and the column removal path (`Removal.H:1316`, after the dispatch). Diagnostic only on production |
| `SurfaceRobinT` | 3069-3074 | column T_surf for `robin_form = surface` | called by kernel `:2004-2009` (result zeroed for listed cells, `:2020`), march `q_col :2607-2610` (discarded when `fl_on_m`, `:2617-2629`), ledger `:4167-4174` (zeroed `:4179`), `PatchDiagnostics:3155-3159` (live, but a column-form diagnostic) |
| `BuildPinEffective` / `PinRule` | 3076-3093, 3098-3115 | the pinned T_pin | no (refused by the face list) |
| ell_flux q_net (flake_length = flux) | 2974-3014 | flux-based heated depth | no (`fr_abort :1052`); CLOSED |
| Sp_field `T_recon_at` | 6181-6187 | lerp of T_surf to the top-cell average over dz/2 | runs (plotted only) |
| `Removal.H` column scan `ScanTAtZv` / `T_at_didx` | Removal.H 68-103, ~1390-1425 | the scan profile | no (after the `face_list` return at `Removal.H:543-547`) |

**Column-top sweeps per production step** (B3 said "about eight"):
- `GatherColumnSurface` three times: `BuildSurfaceTColumns:2882`, `UpdateSpAfterMechanics:6074`
  and `UpdateSpClusters:6378`. All three feed diagnostics only.
- `ColumnSolidRunBelowTop` once (`:2883`).
- `Removal.H`'s inline `k_top_global` once (`Removal.H:512-537`, before the dispatch). Production
  needs it.
- `ColumnTopSolid` in the kernel once (`:1828-1830`). Production needs it.

That is **5 top sweeps + 1 run sweep, each with an MPI reduce**. Two are needed (before and
after removal). The three `GatherColumnSurface` calls and the run sweep exist only for the
COLUMN diagnostics.

---

## Report 3. The criterion sites (B4)

| site | where | reached on production? |
|---|---|---|
| A: Sp_field | `UpdateSpAfterMechanics` `:5992-6343`. Builds `sigma_at_depth` (`:6219-6276`: prescribed / law / scalar / `stress_mf`), `K_I` at `:6296`, `K_Ic` at the tip `T_recon_at(a*icth)` `:6303-6316`, writes `Sp_field` `:6319`. Called at `Advance:1446`. | **Runs every step** (guard: `sp_weibull && owner`). On the face list its output is **plotted and hashed only**: `Sp_field`, `T_surf_field`, `ell_eff_field` and `recon_*`. `FaceRemoval.H` never reads `Sp_field_mf`. The only other readers are `UpdateSpClusters` (`:6386`, diagnostic) and the column path's detach gate (`Removal.H:1100-1104`). It does **not** apply the P2-15 buckling condition or `melt_aware`. Without removal (spall.enabled = 0) it is the only criterion, but the face list requires `spall.enabled`, and the Kant calibration now runs the face list (`kant/grimsel_calibration/input` `FILE =`s `model_form_flame.inputs`). |
| B: column depth scan | `Removal.H:1225-1570` (`UpdateRemovalAfterCohesive`) through the shared pieces `ScanTAtZv` (`Removal.H:68`), `ScanKIcPa` (108), `ScanSigmaT` (132) and `DepthScanH` (158) | **No**. `Removal.H:543-547` returns into `RecedeAlongLines` when `face_list = 1` |
| C: line rule | `FaceRemoval.H:331-935` (`RecedeAlongLines`). Criterion at `FaceRemoval.H:560-591` (same `ScanSigmaT`/`DepthScanH`, with `f_skel = 1.0` and `sn_cth = sn_icth = 1.0` passed as literals). Buckling withdrawal at `:601-636` (`BuckleRatio :315`) | **Yes, this is the production criterion** |

Sites B and C share one copy of the arithmetic, the P2-9 pieces at the top of `Removal.H`. Site
A remains a second hand-written copy: `sigma_at_depth` `:6219-6276` repeats `ScanSigmaT`'s two
branches verbatim, and its `K_Ic` switch `:6307-6316` repeats `ScanKIcPa`. The "two-firing-site
invariant" comments at `:6035-6039` and `:6202-6207` now describe A against B, and B is not the
production path. A port or deletion of site A should either call `ScanSigmaT`/`ScanKIcPa` or go
entirely.

---

## Report 4. Feature interlocks (B14): every abort that refuses a key combination

Value-range checks (finite, > 0, enum spelling) are left out. These are the 185 `Util::Abort`s
in `work/p2_aborts.tsv` filtered to combinations.

**(a) `MMWSpalling::Parse` (`:88-1125`)**

| line | refused combination |
|---|---|
| 120 | any of 21 P2-10 retired keys present (`:98-118`: lateral_*, per_cell_removal, side_face_*, step_diag_radius, jet_wall_*, jet_body*) |
| 130 | `spallation.sp.constraint_depth` present (P2-16) |
| 627 / 630 / 633 / 647 / 651 | damage / cohesive / spall / vapor / sp_weibull without `microstructure.nphases` |
| 641 | `kic.law = grimsel_table` with no phase |
| 666 | `surface.follow_mask` without spall or vapor |
| 670 | `follow_mask` with `surface_patch.mode = prescribed_T` |
| 675 | `follow_mask` or ledger without microstructure |
| 678 | ledger with `amr.max_level > 0` |
| 692 / 697 | AMR + removal without `surface.refine_depth > 0` / with `nsubsteps != 1` |
| 704 / 713 | `robin_form = face` without microstructure / with `epsilon_high < epsilon` |
| 723 / 726 / 730 / 733 | `robin_form = pinned` without microstructure / follow_mask / (spall && sp_weibull) / with `eps_hi < eps` |
| 749 / 752 / 756 / 759 / 763 | `robin_form = surface`: the same four, plus without `spall.surface_T = reconstructed` |
| 784 | `removal_events_csv` without (spall && sp_weibull) |
| 792 | `surface_T` != cell_average, or `removal_energy = flake`, with `spallation.model != sp_weibull` |
| 802 | `removal_energy = flake` without `spall.enabled` |
| 842 | `h_expr` without `follow_mask` |
| 896 | `foot_body_clearance` without removal |
| 950 / 955 | `spall.buckle_span > 0` without `face_list` / without `surface_T = reconstructed` |
| 986-991 | `face_heat = path` with `pad_shield_open >= 0`, `jet_free_surface`, `nozzle_descent = feet`, `jet_above_h`, `robin_form != surface`, or `jet_closure != enthalpy` |
| 995 / 999 / 1002 / 1006 | `face_list` without (convective_flame && h_expr) / robin_form surface / microstructure / follow_mask |
| 1033-1061 (`fr_abort`, defined 1028) | `face_list` with vapor, bit, feet + foot_body_clearance, `melt_aware`, `surface_normal`, `stress_source != local_thermoelastic`, `max_level > 0`, a model other than sp_weibull spall, `surface_T != reconstructed`, `removal_energy != flake`, `flake_length = flux` / `flake_fraction = weibull` / `size_volume = flux`, `flake_coherence_length > 0`, or connected_cluster with `A_crit > 0` |

**(b) `SurfacePatch::Parse` / `CheckPadRing` (`:7887-8356`)**

| line | refused combination |
|---|---|
| 7890 (`CheckPadRing`, called on both exits `:8338`, `:8355`) | `pad_shield_open >= 0` without `foot_r_outer > foot_r_inner` |
| 7915 / 7927 | `face_heat` without `face_list` / `barrier_expr` without `face_heat = path` |
| 7954 | `mode = prescribed_T` with legacy `type != prescribed_temperature` |
| 7970 / 7979 | `robin_form = face` / `pinned` / `surface` without `mode = convective_flame` |
| 8064 / 8067 | `jet_closure = enthalpy` without `h_expr` / with `T_flame_expr` |
| 8105 / 8111 | `jet_T_ent_mode` exhaust or rec_fixed without `jet_entrained_mass` / rec_fixed without `jet_T_rec_fixed` |
| 8134 / 8137 / 8145 | `jet_decay_diameter = momentum` without `jet_stagnation = decay` / without an explicit `jet_core_length`; `jet_De_ref` without momentum |
| 8164 | `jet_fs_aspect`/`exponent` without `jet_free_surface` |
| 8190 / 8194 / 8198 | `jet_above_h` without `jet_T_above` / the converse / `jet_T_above <= jet_T_ent` |
| 8230 | any of 9 D2c/D2v jet keys without `jet_closure = enthalpy` |
| 8243 / 8246 / 8249 | feet without h_expr / `jet_free_surface` without feet / `foot_body_clearance` without feet |
| 8254 / 8257 / 8260 | `h_expr` without convective_flame / without (T_flame_expr or enthalpy) / with `h_conv` or `T_flame` |
| 8277 / 8313 | feet with `nozzle_z0`/`nozzle_feed` / `foot_body_clearance` with `foot_rule != pads` |

**(c) Elsewhere: the `Parse*Settings` helpers and run time**

| line | refused combination |
|---|---|
| 1136 (`Initialize`) | `material` together with `microstructure` (run time, not Parse) |
| 2305 / 2473 (FeetDescent, JetEnthalpyMarch) | feet / jet enthalpy on a non-owner level (run time) |
| 3706 (`BuildFaceList`) | face_list with non-cubic cells (run time, every step) |
| 5069 (`ParseDamageSettings`) | `damage.enabled` with `el.type = disable` and no prescribed stress |
| 5612 (`ParseSpallSettings`) | `flake_h_max > 0 && flake_h_max >= h_spall_max`. **Parse-order defect**: `spall.h_spall_max` is not read until `:5752`, after the early return at `:5681`, so this compares against the member default 0.010 (`:8570`), not the deck's value. A deck with `h_spall_max = 0.015, flake_h_max = 0.012` aborts wrongly, and `h_spall_max = 0.001, flake_h_max = 0.005` passes. Production (0.0015 vs 0.015) passes either way. |
| 5652 / 5673 | `removal_energy = flake` without `surface_T = reconstructed` / (`surface_T` recon or flake) without `local_thermoelastic` |
| 5862 (`ParseBitSettings`) | `bit.enabled` without spall or vapor |
| 9170 / 9173 / 9333 (`ParseSpallationModel`) | sp_weibull with cohesive / with damage; `size_volume = flux` without `V_ref` |
| 6024 / 6031 (`UpdateSpAfterMechanics`, run time) | `local_thermoelastic` with phase 0 lacking E, mu, or physical nu |
| `Removal.H:735` (run time) | sp_weibull spall, neither prescribed_T nor local_thermoelastic, with mechanics disabled |

**Combinations nothing refuses** (silently inert or silently ignored on production):
- `solver = implicit` with a microstructure is ignored (dispatch `:1457`).
- `spallation.sp.prescribed_stress.*` and `prescribed_T_surface.*` only change `Sp_field`
  (`use_presc_s`/`use_presc_T` `:6006-6007`, read at `:6147` and `:6221-6225`); the face list never reads them.
- `spall.prescribed_ratio`, `spall.sample != top_cell`, `spall.h_col_events_csv` and
  `spall.clusters_csv` act on the column path only. The deck comment
  `model_form_flame.inputs:6` says so for `sample`.
- `pinned_idle_cycles` under `surface` is only a `Util::Message` (`:770-774`).
- A beam (`beam.P0 > 0`) with the face list is not refused, but the beam deposits through
  column tops only (report 5).

---

## Report 5. BeamSource and the beam energy balance

- `BeamSource` is at `:4908-4947`: its 20-line header sits in the gap at `:4888-4907` and its
  `AMREX_GPU_HOST_DEVICE` at `:4907`. It has a 10-line lead comment and a 30-line body
  (`:4918-4947`). `BeamIncidentIntensity` is `:4949-4980`, `SolidFacePath` `:4315-4331`,
  `ColumnSolidAboveCount` `:4613-4665`, and `BeamEnergyBalance` `:4716-4887` (32-line header
  + 140 body).
- The beam is parsed by `Numeric::MMWBeam` (`MMWBeam.H:72-77`: `beam.alpha_expr` as a parser
  of T only, default "10"). The power schedule is `MMWBeam::power_at` (`MMWBeam.H:93-100`).
- The beam's lineage in this file: the beam counters (Parse `:140-149`, members `:7583-7598`),
  the kernel block (`:1833-1841`, `:2035-2061`) and the vapour parse.

**Is alpha(T) energy-conserving? No, only for constant alpha.** `BeamSource` uses
`a = alpha_f(T_cell)` (`:4932`) both for the cell's own absorption and for the attenuation over
the whole material path above it, `path = SolidFacePath(n_above)` (`:4933-4943`). The analytic
per-cell integral `(e^{-a s_top} - e^{-a s_bot})/dz` telescopes to `1 - e^{-a s_N}` down a
column only if every cell has the same `a`. With alpha(T) each cell attenuates with its own
`a_k` instead of the cumulative optical depth `sum_j a_j dz`, so the column can over- or
under-deposit. The code says so itself (`BeamEnergyBalance` header `:4723-4727`, "needs the
cumulative optical depth ... required before real alpha(T) MMW work"). The closure check that
would catch it, `P_dep <= P_inc (1 + 1e-9)`, is applied to the **uniform** profile only
(`:4859-4861`). The Gaussian profile carries a further O((dz/z_R)^2) error, because `Pi` uses
the geometric depth z0 - z (`:4974-4975`).
- On the ZHANG kernels the path is geometric (`dz_abs = -1`), so void cavities attenuate the
  beam (report 1).
- On the face list the beam enters only through column tops (`ktc`, `sab`). Listed side, up
  and down faces receive nothing, and `FrFaceFlux` has no beam term. A beam port must decide
  this.

**Per-phase `alpha_attenuation`.** It is parsed with `query_required` at `:347-348` into
`MineralPhase::alpha_atten` (`:8375`), so every microstructure deck must carry it.
`grep -rn alpha_atten src/` finds only those two lines, so it is **never read**. The production
rock file says the same (`rock_grimsel.inputs`: "beam only; never read (MODEL_CARD §12 item
11)"). The beam's attenuation comes from `beam.alpha_expr` alone, with no phase dependence.

---

## Report 6. Everything DEAD (CODE_CENSUS C.7.4 updated)

**Computed, then discarded.**

On the production step:
- `AdvanceMicrostructure`: `SurfaceCellFlux` is evaluated for every in-patch top cell
  (`:2003-2009`), then its `q_in/q_out` are zeroed for every listed cell (`:2020`). Every
  in-patch top cell is listed, because `BuildFaceList` uses the same radius test.
- `JetEnthalpyMarch`, under the face list:
  - the `jTc`/`jkc` gather with its two `ReduceRealSum` (`:2485-2504`) feeds only the
    `SurfaceCellFlux` + `PinRule` call in `q_col` (`:2603-2611`);
  - that call's result is then replaced by the listed faces' sum (`:2617-2629`).
- `EnergyLedger`: `SurfaceCellFlux` per listed cell is zeroed at `:4179`.
- `RefreshKappaEff` called by `FaceRemoval.H:929` (and `Removal.H:2576`) writes `kappa_eff_mf`.
  The kernel then overwrites every cell of it in the same step at `:2072`. Only the
  init/regrid calls matter.
- The Sp_field site, `UpdateSpClusters`, `BuildSurfaceTColumns` and three `GatherColumnSurface`
  sweeps run every step, and their only consumers are plotted fields and thermo columns
  (report 3). They are diagnostic, not dead, but nothing in the production trajectory reads
  them.
- `PatchDiagnostics` (`:3122-3196`) runs under the face list (`:2083-2084`, `sp_on &&
  use_expr`) and reports `patch_P_robin` from the column-form `SurfaceCellFlux`. That is a
  flux the face-list kernel does not apply.

Locals:
- `pc` (`:1926`, `(void)` `:1928`), `p_n` (`:1945-1947`), `q_rob`/`q_los` (`:2063`) in the
  kernel.
- The `phi` patch bound and unread at `:1916-1917` (kernel), `:4779-4780`
  (`BeamEnergyBalance`) and `:5920` (`GatherColumnSurface`).
- `face_top` at `:3912`/`:3920` (`FaceListDiagnostics`).
- `z_lo` (`:6855-6856`, `InitPhi`), `DX` (`:5175`, `UpdateDamageAfterMechanics`).

Parameters:
- `dt` of `BuildSurfaceTColumns` (`(void) dt` `:2950`) and `z_n` of `FaceListPath`
  (`(void) z_n` `:3518`).
- `dt` of `UpdateSpAfterMechanics` / `UpdateSpClusters` / `UpdateCohesiveAfterMechanics`.
- The redundant `(void)a_time; (void)a_step;` in `TimeStepComplete` (`:1192-1193`), which
  uses both above.

Members written or parsed and never read:
- `SurfacePatch::T_ext` (`:7822`, parsed `:8041`).
- `fl_wtop` (`:3290`, written `:3808`, `:3816`).
- `sp_surface_T_mode` (`:8965`, written `:5555-5567`; only the two bools are read).
- `MineralPhase::alpha_atten` (report 5).
- The legacy `surface_patch.type` key (`:7938-7940`), checked only for equality under
  prescribed_T (`:7953-7955`).

Unreachable or inert:
- Invariant (b) of `BeamEnergyBalance` exists only as comment (`:4740-4742`, `:4822-4828`).
- The k_eff bulk diagnostic counts phases rank-locally with no MPI reduction
  (`:7368-7378`, "fine for Step 6 (no MPI tests)"). It is inert for one phase but
  decomposition-dependent for several. It is a defect, not dead code.

**Stranded or stale text.**
- `:622` "P1b: the acceptance test for the support rule is that" (unfinished; still there).
- `:1401-1404`: the Advance header ("no material -> AdvanceConstAlpha") sits above
  `SurfaceRecon` and omits the microstructure branch.
- `:2846-2853`: the C1 pin-rule header sits above `BuildSurfaceTColumns`; it belongs to
  `BuildPinEffective` at `:3076`.
- `:3035-3039`: `SurfaceRobinT`'s header sits above `SharedSurfaceRecon`.
- `:3064-3068`: the orphan "P1 -- THE REFERENCE TEMPERATURE" fragment (still there).
- `:3198-3233`: `SurfaceCellFlux`'s 36-line header sits **694 lines** above the function
  (`:3928`); `:3234` is the `RiserWeight` tombstone.
- `:3939-3954`: the "P1f ... the lateral criterion needs the face temperature" rationale. The
  lateral criterion was deleted in P2-10; `T_face_out` is now read by `FrFaceFlux`.
- `:4479-4484`: `ColumnTopSolid`'s header sits above `UpdatePropertyFields`, 105 lines before
  the body (`:4589`). The census said 500 lines; the distance shrank but the comment is still
  stranded.
- `:5588-5598`: the "P1 -- THE GRADIENT-AWARE DRIVING STRESS, `absolute` (default)" block
  describes a key with no parse (deleted at both sites, `:6209`).
- `:5899-5905`: `GatherColumnSurface`'s header says "phi >= 0". The rule has been
  removed-only since D2q-1 (`:5930`).
- `:5942-5949`: `UpdateSpAfterMechanics`'s header sits above `ColumnSurfaceCos`.
- `:6574-6578`: `LawSigma` says the Sp_field site sits in "`amrex::ParallelFor` under
  `[=, this]`". It is `amrex::LoopOnCpu` with `[&]` (`:6136`), and the comment points at the
  deleted lateral criterion.
- `:8360-8361`: the duplicated "Step 5: Voronoi mineral microstructure" line.
- `:8764-8774` and `:8796`: P1/P1b support-rule text with no members (`spall.support_rule`
  is parsed nowhere).
- Stale line references:
  - `:637-638` (":206", ":293")
  - `:787`
  - `:2888` (`kc_cell`)
  - `:2967` (`Removal.H:1503-1517`)
  - `:2991` (":7920-7925")
  - `:4500` (":1436")
  - `:5478-5481` (`Removal.H:691`, ":4668")
  - `:5669`
  - `:6584` (`Removal.H:1833`)
  - `:7215` (":5380")
  - `:7703` (`Removal.H:589-594`)
  - `:8932` (`Removal.H:846-858`)

**C.7.4, item by item**

| C.7.4 item | now |
|---|---|
| kernel `pc`, `p_n`, `q_rob`, `q_los` discarded | **still** (`:1926-1928`, `:1945-1947`, `:2063`) |
| `BuildSurfaceTColumns` `dt` unused | **still** (`:2950`) |
| Removal.H PASS 3 temperature array | Removal.H; not re-checked here (part 3) |
| `SurfacePatch::Parse` legacy `type` | **still** (`:7938-7955`) |
| `BeamEnergyBalance` invariant (b) comment only | **still** (`:4740-4742`, `:4822-4828`) |
| Removal.H P1b support-rule comment blocks | Removal.H; grep still finds P1b text at `Removal.H:640, 1448, 1611, 2058, 2213, 2321, 2490` |
| `SurfaceRobinT` orphaned comment | **still** (`:3064-3068`) |
| `ColumnTopSolid` comment 500 lines above | **still stranded**, now 105 lines above (`:4479-4484` vs `:4589`) |
| `SurfaceCellFlux` `AMREX_GPU_HOST_DEVICE` attached to nothing | **gone** (`SurfaceCellFlux :3928` has no macro); its header comment is now stranded instead |
| Parse "P1b: the acceptance test ..." | **still** (`:622`) |
| `BuildFlameColumns` pad-ring block pasted twice | **still** (`:2166-2184` and `:2205-2222`) plus a third code copy in `FaceHr :3686-3689` |
| `UpdateSpAfterMechanics` "READ 1 of 2" pasted twice | **gone** (once, `:6278`) |

---

## Report 7. Line budget

From `work/p2_ranges.py`. Every line is counted once. CLOSED lines are counted inside their
lineage and also flagged.

| lineage | lines | code | comment | blank | of which CLOSED |
|---|---|---|---|---|---|
| PROD-SHARED | 2,969 | 1,827 | 1,003 | 139 | 0 |
| FLAME | 1,686 | 1,264 | 396 | 26 | 170 |
| BEAM | 526 | 307 | 200 | 19 | 45 |
| ZHANG | 414 | 313 | 72 | 29 | 0 |
| MICRO | 543 | 349 | 164 | 30 | 0 |
| MECH | 354 | 299 | 35 | 20 | 0 |
| DAMAGE | 640 | 523 | 63 | 54 | 0 |
| COLUMN | 1,634 | 992 | 573 | 69 | 223 |
| DIAG | 564 | 418 | 126 | 20 | 0 |
| DEAD (stranded comment) | 45 | 0 | 42 | 3 | 0 |
| **total** | **9,375** | 6,292 | 2,674 | 409 | 438 |

**Removal scenario** (DAMAGE + COLUMN-only + CLOSED switches + ZHANG + stranded text):

| removed | lines |
|---|---|
| DAMAGE | 640 |
| COLUMN (incl. its 223 CLOSED) | 1,634 |
| ZHANG | 414 |
| CLOSED outside COLUMN: FLAME 170 (pad shield, above-plane pair, rec_fixed, robin_form = face aborts) + BEAM 45 (`melt_aware`) | 215 |
| DEAD stranded text | 45 |
| **total removed** | **2,948** |
| **remaining** | **6,427** (6,472 if the stranded text is kept) |

Adjustments to that 6,427:
- The ZHANG merge also makes `UpdateModel`'s homogeneous branch redundant (`:1361-1398`, about
  39 MECH lines) and adds roughly 40-80 lines: a generic T-dependent phase law in
  `UpdatePropertyFields` and the `law = zhang` parse. The net is about **6,430 ± 20** after
  the merge.
- Further FLAME cuts that are column-only on production and are not in the brief's list:
  - the pinned closure (refused by the face list): about 135 lines (`:717-740`,
    `:828-838`, `:1891-1897`, `BuildPinEffective`, `PinRule`, `:4023-4033`,
    `:7614-7632`, `:7976-7998`);
  - `face_heat = path`, which the user OK'd deleting: about 385 lines (`FaceListPath`,
    `FrPathOffsets`, `FaceCentre`, `:3303-3309`, `:3855-3863`, `:975-992`, `:1096-1111`,
    `:7810-7815`, `:7911-7934`, the march's path branch `:2631-2672`);
  - `PatchDiagnostics` (80).
  - Together these bring the file to about **5,830**.
- Deleting COLUMN also deletes plotted fields: `Sp_field`, `Sp_cluster_id`, `T_surf_field`,
  `ell_eff_field`, `h_spall_field`, `crack_event_count` and the `recon_*` thermo columns. That
  re-pins the gate header. It also requires editing the production callers
  `FaceRemoval.H:928` (`UpdateSurfaceMaskFromPhi`) and `Advance:1434/:1446/:1450`.
- MECH (354) stays by the brief's rules but is unreachable on production (`el.type = disable`).

---

## Report 8. The two closed-study branches at the criterion/removal sites a beam port will edit

The face list refuses both (`fr_abort`, `:1037` and `:1038`). The production criterion
therefore never takes either branch, but both still thread through the shared criterion
pieces the face list calls, as literal 1.0 arguments. Every site follows.

**`spall.melt_aware` (f_skel / MeltSkeleton)**
- `MMWSpalling.H:5502-5526`: parse of `spall.melt_aware` and `melt.lambda_crit`, with the
  range abort at `:5525`.
- `MMWSpalling.H:8599-8605`: members `melt_aware` and `melt_lambda_crit`, and their comment.
- `MMWSpalling.H:8606-8618`: `MeltSkeleton` and its lead comment.
- `MMWSpalling.H:1037`: the face-list refusal.
- `Removal.H:128-146`: `ScanSigmaT`'s **`f_skel` parameter**, which multiplies both branches
  (`:138`, `:143`). The face list calls it with the literal `1.0` at `FaceRemoval.H:577`. A
  port deletes the parameter and that argument; multiplying by exactly 1.0 is an identity, so
  the deletion is bit-safe.
- `Removal.H:1017-1057`: the `LL_col_prof` gather and its `ReduceRealSum`, and `LL_col_p`.
- `Removal.H:1366-1386`: `LLcol_i` and the `LL_at_didx` lambda.
- `Removal.H:1464-1479`: `f_skel = MeltSkeleton(LL_at_didx(...))`.
- `Removal.H:1501-1502`: `f_skel` in the `stress_mf` branch.
- `Removal.H:1527`: `f_skel` passed to `ScanSigmaT`.
- Site A (`UpdateSpAfterMechanics`) **never applied** the skeleton factor, so the two sites
  disagree under the key.
- Inputs: 9 set it per C.7.3, and `model_form_flame.inputs` sets `spall.melt_aware = 0`. By the
  P2-10 pattern the key moves to the retired-key list at `:98-118`.

**`spall.surface_normal` (the criterion half; the recession half is already superseded)**
- `MMWSpalling.H:5448-5501`: parse plus the `Util::Message` at `:5493-5500`.
- `MMWSpalling.H:8593-8598`: member `spall_surface_normal`.
- `MMWSpalling.H:8619-8622`: `SP_SN_COS_MIN`.
- `MMWSpalling.H:8551-8557`: `sn_recession_superseded`.
- `MMWSpalling.H:938-940`: its thermo registration.
- `MMWSpalling.H:1038`: the face-list refusal.
- `MMWSpalling.H:5942-5990`: `ColumnSurfaceCos`, including the stranded Sp header.
- Site A:
  - `:6081-6090` (`sn_cos_col`, `sn_cos_p`)
  - `:6189-6194` (`icth`)
  - `:6220` (`dv = depth * icth`)
  - `:6300-6304` (`T_kic` at `a_loc * icth`)
- `Removal.H:548-558`: `sn_cos_col`/`sn_cos_p`, after the face-list return.
- `Removal.H:158-216`: `DepthScanH`'s **`sn_cth`/`sn_icth` parameters**, used at `:174-175`,
  `:179`, `:201-204`, `:209-211` and `:215`. The face list passes `1.0, 1.0` at
  `FaceRemoval.H:591`. A port deletes both parameters and the two literals; this is bit-safe.
- `Removal.H:1224-1231`: per-column `sn_cth`/`sn_icth`.
- `Removal.H:1460-1463`: `z_v = z_abs * sn_icth`.
- `Removal.H:1544`: the `DepthScanH` call.
- `Removal.H:1558-1569`: `h_vert = h_col * sn_icth`, the superseded recession.
- `Removal.H:1807-1810`: `sn_active` switches off the coherence cap.
- `Removal.H:2000-2009`: the comment.
- `Removal.H:2403-2424`: `sn_recession_superseded`.
- Inputs: 4 set it per C.7.3, two of them inside the 2,639-file sweep (`unit/beam_void_closure`
  and the Meier `input_2d_dev` leg). Removing it re-pins those.

---

## Appendix. What from CODE_CENSUS C.2 no longer exists

| C.2 entry | status now |
|---|---|
| C.2.4 `UpdateColumnPhase`, `AssertStepSplit`, `WallStepDiagnostic` | **deleted** (0 grep hits) |
| C.2.5 the P2-1 wall law: `WJ_*`, `WallMu`, `WallZFill`, `WallJetBranch`, `WallJetH`, `SideFaceH`, `ColumnWallCorner`, `PatchAzimuthFraction`, `UpdateWallLawState` | **deleted** (0 hits; `side_face_h_law` is a retired-key abort, `:107`) |
| C.2.5 the D2z helpers `ColPhi`, `PairDelta`, `StepCellWeight`, `StaircaseFactor`, `AreaOfDelta`, `StepConductance` | **deleted** (0 hits) |
| C.2.3 `RiserWeight` | **deleted** (tombstone comment `:3234`) |
| C.2.9 `PerCellRemoval` | **deleted** (`spall.per_cell_removal` retired, `:102`) |
| C.2.1 `ParseSpallSettings` P1f parse-time `SurfaceCellFlux` self-test (B14) | **deleted** |
| C.2.2 kernel `col_in_patch` side-face block; C.2.4 march `SideHit` gather; C.2.10 the ledger's side-face/area-law/D2x/D2z/P1e statistics | **deleted**; the face list (`FrFaceFlux`) replaces all four side-face copies (B2 is closed) |
| C.2.1 `SurfacePatch` side-face block, `side_face_closure` (`SideTForm`) | **deleted** (tombstone `:7797-7798`) |
| P2-2e constraint depth (`spallation.sp.constraint_depth`) | **deleted** in P2-16 (retired-key abort `:130`) |
| carried state `side_T_eff`, `lat_accum_col`, `col_t_last`/`col_dt`/`col_nvoid`, `wall_corner`, `side_cols_above`, `wall_step_*` | **gone** (0 hits). `T0_col_prof` survives only in `Removal.H` |
| C.2.3 claim "SurfaceCellFlux: 20 positional params, 8 call sites" | now 18 positional + 4 defaulted (`:3928-3937`), **7** call sites (kernel, march, PatchDiagnostics, ledger, `FrFaceFlux` x2, `RecedeAlongLines`) |
| sizes | `Parse` 1,100 → 1,039; `AdvanceMicrostructure` 583 → 298; `JetEnthalpyMarch` 569 → 376; `EnergyLedger` 699 → 250; `ParseSpallSettings` 508 → 412; `BuildFlameColumns` 170 → 137 (body); file 11,033 → 9,375 |
| new since C.2 (none of these are in CODE_CENSUS) | face list: `FaceEntry`, `BuildFaceList`, `FaceListHeat`, `FaceListTg`, `FaceListDiagnostics`, `FaceH`/`FaceHr`, `FaceWeight`, `FrGauss`/`FrSmoothAt`/`FrSmoothedGradient`, `FrInitFraction`, `FillFractionGhosts`; path: `FrPathOffsets`, `FaceListPath`, `FaceCentre`; P2-15 buckle members + Parse block `:943-965`; retired-key guards `:90-130` |
