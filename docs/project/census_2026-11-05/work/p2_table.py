#!/usr/bin/env python3
"""Emit the function table (markdown) for 2_regions_main.md: split_raw.tsv counts +
lineage composition from p2_lineage_by_line.tsv + hand-written guard/purpose (ANN)."""
import collections
W = 'docs/project/census_2026-11-05/work/'
L = {}
for l in open(W + 'p2_lineage_by_line.tsv'):
    a = l.rstrip('\n').split('\t'); L[int(a[0])] = (a[1], a[2] == '1')
rows = [l.rstrip('\n').split('\t') for l in open(W + 'p2_split_raw.tsv')]
P = 'production (Meier/Kant face-list deck)'
ANN = {
'MMWSpalling': ('always', 'constructors; the ParmParse one calls Parse'),
'Parse': ('always (queryclass)', 'every key: retired-key aborts, base parses, the Parse*Settings calls, microstructure/phase parse, 60 RegisterNewFab and 125 RegisterIntegratedVariable calls, 42 aborts, the cross-key abort matrix, face-list validation'),
'Initialize': ('always', 'enthalpy table (material only), base inits, T_init copy, phi/seed void (removal), mask, microstructure, T->H'),
'StoreInitialTemperature': ('microstructure (T_init_mf exists)', 'the one write of T_init_mf, the stress-free reference'),
'TimeStepBegin': ('always; refresh only if microstructure && max_level>0 && el.type!=disable', 'base hooks + microstructure refresh for AMR mechanics'),
'TimeStepComplete': ('always; same refresh guard', 'base hooks + the same refresh again (B15 double refresh)'),
'Regrid': ('AMR regrid', 'removal repair (lev>0), mask, microstructure re-init, H->T, mechanics repair'),
'UpdateModel': ('el.type != disable (Mechanics base skips it otherwise)', 'nodal elastic model from cell E/mu/free strain; Hu prescribed_T top-node override; homogeneous alpha[] branch'),
'struct SurfaceRecon': ('type', 'result of ReconFromX (T_surf, ell, h_f, T_bar, flags); preceded by a stranded, stale Advance header (1401-1404)'),
'Advance': ('always', 'the step: surface T columns, property refresh, mechanics, damage, cohesive, Sp_field, clusters, removal, swap, one conduction kernel, H->T'),
'AdvanceConstAlpha': ('no material block, no microstructure', 'constant-alpha Laplacian + beam + cell-average losses (Step 2-4 regressions)'),
'AdvanceMaterial': ('material.type set, solver=explicit, no microstructure', 'explicit kappa(T) harmonic-face conduction, void-aware, prescribed_T/const flame patch, beam, losses'),
'AdvanceMaterialImplicit': ('material.type set, solver=implicit', 'backward-Euler MLABecLaplacian, single level, Neumann-only, ignores voids'),
'AdvanceMicrostructure': ('microstructure.nphases set (' + P + ')', 'THE production kernel: face conductivities, top-face SurfaceCellFlux, face-list FrFaceFlux, beam, H update; calls BuildFaceList/BuildFlameColumns/BuildPinEffective and the diagnostics'),
'BuildFlameColumns': ('surface_patch.enabled && owner (microstructure)', 'per-column flame_h/Tg/s (constants or h_expr at stand-off from FeetDescent or prescribed feed); pad shield x2; calls FaceListHeat, FrAssignTs, JetEnthalpyMarch, FaceListTg'),
'FaceListTg': ('face_list=1 (fl_ready)', 'copies column flame_Tg to each face (position); path: per-face Tg and the fl_s/fl_Tg fields'),
'FeetDescent': ('surface_patch.nozzle_descent=feet', 'burner-on-feet nozzle height from pad quantiles; arms foot_clear_k (foot_body_clearance)'),
'JetEnthalpyMarch': ('surface_patch.jet_closure=enthalpy', 'serial radial wall-jet enthalpy march: stagnation state, per-bin draw, dilution, profile CSV, 1e-9 invariant abort'),
'BuildSurfaceTColumns': ('spall.surface_T in {diagnostic, reconstructed} && owner', 'per-column two-cell exponential reconstruction (surf_T/ell/h_f/T_bar/flag) + flat-hold + flake_length=flux q_net length'),
'SharedSurfaceRecon': ('callers: Sp_field site (recon_on), column removal path', 'rebuilds a SurfaceRecon from the per-column arrays'),
'SurfaceRobinT': ('returns 0 unless robin_form=surface', 'column-top reconstructed surface T for SurfaceCellFlux'),
'BuildPinEffective': ('robin_form=pinned && sp_on', 'pin_T_eff/pin_reason per column via PinRule'),
'PinRule': ('robin_form=pinned', 'the C1 idle-clock pin rule for one column'),
'PatchDiagnostics': ('owner && (pinned || (sp_on && use_expr))', 'pinned branch counters + patch_P_robin by the COLUMN form (not the flux the face list applies)'),
'struct FaceEntry': ('type', 'one listed solid/void face; preceded by SurfaceCellFlux\'s stranded 37-line header (3198-3234) and the P2-8 banner'),
'FrPathOffsets': ('face_heat=path', 'Dijkstra neighbour offsets |o|<=4'),
'FaceListPath': ('surface_patch.face_heat=path', 'P2-11b: surface-path distance from stagnation by Dijkstra, barrier, h at r:=s_path'),
'FaceCentre': ('callers only under face_heat=path', 'a face centre (FaceListHeat inlines the same arithmetic)'),
'FrGauss': ('face_list=1', 'Gaussian weight for the smoothed fraction'),
'FrSmoothAt': ('face_list=1', '5x5x5 smoothed solid fraction'),
'FrSmoothedGradient': ('face_list=1', 'central-difference gradient of the smoothed fraction (the normal)'),
'FrInitFraction': ('face_list=1 (once)', 'solid fraction 1/0 from removed_mf'),
'FillFractionGhosts': ('face_list=1', 'fraction ghost rule: gas above, mirror at symmetry planes, clamp elsewhere'),
'FaceWeight': ('face_list=1', 'F3 area weight 1/(|nx|+|ny|+|nz|)'),
'FaceH': ('face_list=1, face_heat=position', 'F2 h at a face centre (r from the patch centre)'),
'FaceHr': ('face_list=1', 'h_expr at (x,y,r,s) with jet zeroing above the plane and the pad shield'),
'BuildFaceList': ('face_list=1 && sp_on && removal && !prescribed_T', 'F1+F3: every solid/void face in the patch, gathered/sorted on every rank, per-cell index'),
'FaceListHeat': ('face_list=1 (inside BuildFlameColumns, use_expr)', 'F2 over the list (position) or FaceListPath (path); fl_h field; top-face mismatch witness'),
'FaceListDiagnostics': ('face_list=1 (fl_on)', 'per-class robin power fl_P_*, fl_id_march'),
'SurfaceCellFlux': ('every non-void owner cell (microstructure); FrFaceFlux on the face list', 'the one Robin/loss routine: surface, face-Newton, pinned and cell forms'),
'EnergyLedger': ('energy_ledger.enabled && lev==0 (on the production deck)', 'compensated H sums, Robin/loss re-evaluation (face-list mirror), side-power CSV, ledger_err/_solid'),
'SolidFacePath': ('beam on && void-aware', 'Beer-Lambert path (n_above+1/2)dz'),
'SurfaceOwner': ('always', 'finest level owns the surface'),
'RepairRemovalStateAfterRegrid': ('Regrid, lev>0 && removal (AMR; face list refuses max_level>0)', 'rebuild column-linear phi / clean removed after regrid'),
'CheckSurfaceCovered': ('column removal path, max_level>0', 'abort if the surface left the refined band'),
'UpdatePropertyFields': ('some phase law=grimsel (' + P + ')', 'per-step Grimsel refresh of kappa, E, mu, rho, Cp, secant beta + prop_* witnesses; preceded by ColumnTopSolid\'s stranded header (4479-4484)'),
'ColumnTopSolid': ('void-aware && (beam || follow_mask || ledger); BuildFaceList; mech clearance; tagging', 'top solid cell per column, ReduceIntMax'),
'ColumnSolidAboveCount': ('beam on && void-aware', 'solid cells above each (col,k), ncols*nz reduction'),
'ColumnSolidRunBelowTop': ('BuildSurfaceTColumns; column removal path', 'contiguous solid run below each column top'),
'BeamEnergyBalance': ('beam power > 0 (microstructure kernel epilogue)', 'E1 beam closure P_inc/P_dep/P_trans, invariant (a), uniform double-count guard'),
'BeamSource': ('every kernel; returns 0 at P=0', 'Beer-Lambert x Gaussian/uniform source, analytic per-cell integral'),
'BeamIncidentIntensity': ('beam on', 'incident intensity Pi(x,y,z)'),
'TagCellsForRefinement': ('AMR; band only if removal && surface.refine_depth>0', 'base tagging + surface band'),
'ParseDamageSettings': ('always (keys default off)', 'damage.* Drucker-Prager parse and checks'),
'RefreshKappaEff': ('init/regrid; also called by both removal paths and damage', 'kappa_eff = kappa exp(-a_D D) (plotted)'),
'DamageGrowthRHS': ('damage.enabled', 'damage-rate law'),
'IntegrateDamageRK4': ('damage.enabled', 'RK4 of the damage law'),
'UpdateDamageAfterMechanics': ('damage.enabled', 'Drucker-Prager damage update from stress_mf or prescribed stress'),
'ParseCohesiveSettings': ('always (cohesive.enabled default 0)', 'cohesive.* parse'),
'UpdateCohesiveAfterMechanics': ('cohesive.enabled', 'bilinear GB cohesive law (no feedback to the solve)'),
'ResetCohesiveDiagnostics': ('cohesive.enabled', 'zero gb_* fields'),
'ParseSpallSettings': ('always; most keys after `if (!spall_enabled) return` (:5681)', 'spall.* keys: stress source/law, surface_normal, melt_aware, CSVs, surface_T, flake keys, removal_energy, damage-law keys, h_spall_max, coherence, seed_void'),
'ParseVaporSettings': ('always (vapor.enabled default 0)', 'vapor.* parse'),
'ParseBitSettings': ('always (bit.enabled default 0)', 'bit.* parse'),
'PrincipalStressRatio': ('damage-law spall (column path)', 'sigma_1/sigma_3 ratio'),
'GatherColumnSurface': ('BuildSurfaceTColumns, Sp_field site, UpdateSpClusters', 'top solid cell per column (second copy of ColumnTopSolid; domain top if no removal); binds phi unused'),
'ColumnSurfaceCos': ('spall.surface_normal=1', 'per-column wall-slope cos(theta); preceded by UpdateSpAfterMechanics\'s stranded header (5942-5949)'),
'UpdateSpAfterMechanics': ('spallation.model=sp_weibull && owner (runs on the production deck, output plotted only)', 'criterion site A: Sp = K_I/K_Ic at each column top into Sp_field; recon witnesses'),
'UpdateSpClusters': ('sp_weibull && weibull.enabled && owner (runs on the production deck)', 'top-plane union-find of Sp>=1 into Sp_cluster_id; A_crit gate; clusters CSV'),
'LawSigma': ('spall.stress_law=grimsel', 'the one Grimsel thermoelastic stress entry'),
'RockLaw': ('spall.stress_law=grimsel', 'phase-0 law pointer or nullptr'),
'BuildEnthalpyTable': ('material block', 'global-material H<->T table'),
'SetPhaseFractionsConstAlpha': ('no material, no microstructure', 'closed-form phase fractions'),
'ConvertTtoH': ('Initialize', 'T->H by phase table / global table / closed form'),
'InvertHtoT': ('every step; regrid; removal', 'H->T and phase fractions, same three-way dispatch; production is only the per-phase branch 6716-6723, the rest is the global-material/const-alpha path'),
'InitSurfaceMask': ('init/regrid/refresh', 'surface_mf: domain top, or phi window under removal'),
'UpdateSurfaceMaskFromPhi': ('removal (init + both removal paths)', 'surface_mf = phi window; read only when follow_mask=0 and by the ledger\'s missing-column count'),
'ApplySeedVoid': ('spall.seed_void given (debug)', 'void seeded index boxes once'),
'InitPhi': ('removal', 'flat level set, removed=0, reset per-mode fields'),
'RefreshMicrostructureDerivedFields': ('microstructure && max_level>0 && el.type!=disable', 'mask + InitializeMicrostructure on every level'),
'RepairMechanicsAfterRegrid': ('el.type != disable', 'refresh model and zero stress/strain after regrid'),
'FiniteMatrix': ('damage.enabled', 'finite-entries check'),
'AlphaBetaFraction': ('mechanics / free strain', 'smoothstep alpha-beta fraction'),
'ThermalFreeStrain': ('mechanics / free strain', 'beta(T-Tref) + alpha-beta strain'),
'CellToNodeAverageInDomain': ('UpdateModel', 'in-domain cell-to-node average'),
'UpdateMicrostructureFreeStrain': ('UpdateModel; InitializeMicrostructure', 'fills thermal_free_strain_mf (plotted)'),
'InitializeMicrostructure': ('microstructure (Initialize, Regrid, AMR mechanics refresh)', 'phase tables, Voronoi/expression phase+grain ids, property fill, ghost fill, GB flags, damage/cohesive reset, k_eff, flaws'),
'InitializeFlaws': ('weibull.enabled', 'hashed per-cell Weibull flaw length'),
'enum class SolverKind': ('type', 'explicit | implicit'),
'WeibullSizeFactor': ('weibull.size_volume=flux (face list refuses it)', '(V/V_ref)^(1/m) size factor'),
'struct SurfacePatch': ('surface_patch.* (queryclass)', 'the patch settings (~70 members) and their 460-line Parse: mode, robin_form, face_list/face_heat, pad shield, jet closure + D2b-D2v options, feet, h_expr'),
'struct MineralPhase': ('type', 'per-phase scalars, Grimsel law, alpha-beta; alpha_atten parsed, never read'),
'enum class MicroMode': ('type', 'voronoi | expression'),
'enum class SpallSampleMode': ('type', 'top_cell | column_max | shallow_max (damage law)'),
'MeltSkeleton': ('spall.melt_aware=1 (column path only)', 'f(Lambda_L) = clamp(1 - LL/lambda_crit)'),
'removal_enabled': ('-', 'spall || vapor'),
'enum class SpallationModel': ('type', 'damage_law | sp_weibull'),
'enum class KIcLaw': ('type', 'nasseri | linear | grimsel_table'),
'enum class DetachMode': ('type', 'per_face | connected_cluster'),
'enum class SurfaceTMode': ('type', 'cell_average | diagnostic | reconstructed (the stored mode is never read)'),
'RawSurfaceX': ('BuildSurfaceTColumns; FrAssignTs; RecedeAlongLines', 'two-cell x = ln((T_top-T0)/(T_1-T0))'),
'ReconFromX': ('same callers', 'the reconstruction algebra: T_surf, ell, h_f (c_f, h_spall_max cap, flake_h_max bound), T_bar'),
'SurfaceGasTempCap': ('reconstruction callers', 'gas-temperature ceiling max(T_flame, jet_T_nozzle)'),
'SplitMix64': ('weibull.enabled', 'hash for the flaw RNG'),
'ParseSpallationModel': ('always; sp_* after `if model != sp_weibull return` (:9164)', 'model, confining.p, buckle_span, sp.*, kic.*, prescribed sigma/T, weibull.*'),
}
SUB = {  # (parent, start, end, name, guard, purpose)
'Parse': [(88, 130, '  retired-key guards', 'pp.contains(key)', 'abort on any P2-10/P2-16 retired key'),
          (292, 624, '  microstructure block', 'microstructure.nphases', 'phase parse + field registrations (order = plotfile layout)'),
          (626, 941, '  cross-key matrix + thermo registrations', 'various', 'the late abort matrix and keyed thermo columns'),
          (943, 965, '  buckle_span block', 'spall.buckle_span>0', 'P2-15 aborts + bk_* witnesses'),
          (967, 1112, '  face-list block', 'surface_patch.face_list=1', 'path aborts, F-requirements, fr_abort refusals, fl_/fr_ registrations')],
'AdvanceMicrostructure': [(1930, 1966, '  lambda k_face', 'per cell', 'harmonic face conductivity, void zero, GB contact resistance')],
'JetEnthalpyMarch': [(2601, 2629, '  lambda q_col', 'per column', 'column Robin draw; under the face list the SurfaceCellFlux/PinRule part is discarded and the listed faces are billed'),
                     (2662, 2668, '  lambda q_face', 'face_heat=path', 'one face billed alone')],
'EnergyLedger': [(4092, 4097, '  lambda neumaier', '-', 'compensated sum')],
'UpdateSpAfterMechanics': [(6181, 6188, '  lambda T_recon_at', 'recon', 'lerp in the top half-cell'),
                           (6214, 6218, '  lambda T0_at', '-', 'T_init with T_ref fallback'),
                           (6219, 6274, '  lambda sigma_at_depth', '-', 'prescribed / local thermoelastic (law or scalar) / stress_mf stress at depth')],
'struct SurfacePatch': [(7887, 7895, '  CheckPadRing', 'pad_shield_open>=0', 'pad ring must exist (called on both Parse exits)'),
                        (7897, 8356, '  SurfacePatch::Parse', 'surface_patch.enabled', 'every surface_patch.* key, ~50 aborts')],
}
order = ['PROD','FLAME','BEAM','ZHANG','MICRO','MECH','DAMAGE','COLUMN','DIAG','DEAD']
def comp(a, b):
    c = collections.Counter(L[x][0] for x in range(a, b + 1))
    cl = sum(1 for x in range(a, b + 1) if L[x][1])
    main = c.most_common(1)[0][0]
    rest = ', '.join(f'{k} {v}' for k, v in sorted(c.items(), key=lambda kv: -kv[1]) if k != main)
    s = main + (f' (+{rest})' if rest else '')
    if cl: s += f'; CLOSED {cl}'
    return s
out = ['| # | name | lines | count | code / comment / blank (+lead comment) | lineage | guard (runs only if) | one-line purpose |',
       '|---|---|---|---|---|---|---|---|']
n = 0
for r in rows:
    if r[0] != '2': continue
    nm, a, b, cnt, code, com, bl, lc = r[1], int(r[2]), int(r[3]), r[4], r[5], r[6], r[7], r[8]
    ls = int(r[9])
    g, p = ANN.get(nm, ('?', '?'))
    if nm == 'MMWSpalling' and a == 86: p = 'ParmParse constructor -> Parse'
    n += 1
    g = g.replace('|', '\\|'); p = p.replace('|', '\\|')
    out.append(f'| {n} | `{nm}` | {a}-{b} | {cnt} | {code} / {com} / {bl} (+{lc}) | {comp(ls, b)} | {g} | {p} |')
    for (sa, sb, sn, sg, sp) in SUB.get(nm, []):
        sg = sg.replace('|', '\\|'); sp = sp.replace('|', '\\|')
        out.append(f'|  | {sn} | {sa}-{sb} | {sb-sa+1} |  | {comp(sa, sb)} | {sg} | {sp} |')
open(W + 'p2_table.md', 'w').write('\n'.join(out) + '\n')
print(n, 'rows')
