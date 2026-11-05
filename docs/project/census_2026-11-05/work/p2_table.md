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
