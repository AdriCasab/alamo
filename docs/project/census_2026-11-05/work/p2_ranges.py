#!/usr/bin/env python3
"""Line-by-line lineage assignment of src/Integrator/MMWSpalling.H and the line budget.
Every line 1..N gets exactly one lineage and a CLOSED flag (a default-off key whose study
closed, per CODE_CENSUS C.7.3 / the face-list refusals). Base assignment per function
(lead comment + body, from split_raw.tsv) or per explicit member range; later entries in
OVERRIDES win. Prints coverage gaps (must be none) and the budget."""
import sys, collections
SRC = 'src/Integrator/MMWSpalling.H'
lines = open(SRC).read().split('\n')
if lines and lines[-1] == '': lines = lines[:-1]
N = len(lines)
# comment/blank classification (same rules as split.py)
import re
def kinds():
    out = []; inb = False
    for l in lines:
        s = l.strip(); code = False; com = inb; i = 0
        while i < len(l):
            if inb:
                j = l.find('*/', i)
                if j < 0: i = len(l); break
                inb = False; i = j + 2; continue
            if l.startswith('//', i): com = True; break
            if l.startswith('/*', i): inb = True; com = True; i += 2; continue
            if not l[i].isspace(): code = True
            i += 1
        out.append('blank' if (s == '' and not com) else ('code' if code else 'comment'))
    return out
K = kinds()
lin = [None] * (N + 1); closed = [False] * (N + 1)
def put(a, b, L, c=False):
    for x in range(a, b + 1): lin[x] = L; closed[x] = c
# --- base: functions from the splitter (lead comment start .. end)
FUNC = {
 'MMWSpalling':'PROD','Parse':'PROD','Initialize':'PROD','StoreInitialTemperature':'PROD',
 'TimeStepBegin':'PROD','TimeStepComplete':'PROD','Regrid':'PROD','UpdateModel':'MECH',
 'struct SurfaceRecon':'PROD','Advance':'PROD','AdvanceConstAlpha':'ZHANG','AdvanceMaterial':'ZHANG',
 'AdvanceMaterialImplicit':'ZHANG','AdvanceMicrostructure':'PROD','BuildFlameColumns':'FLAME',
 'FaceListTg':'PROD','FeetDescent':'FLAME','JetEnthalpyMarch':'FLAME','BuildSurfaceTColumns':'COLUMN',
 'SharedSurfaceRecon':'COLUMN','SurfaceRobinT':'COLUMN','BuildPinEffective':'FLAME','PinRule':'FLAME',
 'PatchDiagnostics':'DIAG','struct FaceEntry':'PROD','FrPathOffsets':'FLAME','FaceListPath':'FLAME',
 'FaceCentre':'FLAME','FrGauss':'PROD','FrSmoothAt':'PROD','FrSmoothedGradient':'PROD',
 'FrInitFraction':'PROD','FillFractionGhosts':'PROD','FaceWeight':'PROD','FaceH':'PROD','FaceHr':'PROD',
 'BuildFaceList':'PROD','FaceListHeat':'PROD','FaceListDiagnostics':'DIAG','SurfaceCellFlux':'PROD',
 'EnergyLedger':'DIAG','SolidFacePath':'BEAM','SurfaceOwner':'PROD','RepairRemovalStateAfterRegrid':'COLUMN',
 'CheckSurfaceCovered':'COLUMN','UpdatePropertyFields':'PROD','ColumnTopSolid':'PROD',
 'ColumnSolidAboveCount':'BEAM','ColumnSolidRunBelowTop':'COLUMN','BeamEnergyBalance':'BEAM',
 'BeamSource':'BEAM','BeamIncidentIntensity':'BEAM','TagCellsForRefinement':'COLUMN',
 'ParseDamageSettings':'DAMAGE','RefreshKappaEff':'DAMAGE','DamageGrowthRHS':'DAMAGE',
 'IntegrateDamageRK4':'DAMAGE','UpdateDamageAfterMechanics':'DAMAGE','ParseCohesiveSettings':'DAMAGE',
 'UpdateCohesiveAfterMechanics':'DAMAGE','ResetCohesiveDiagnostics':'DAMAGE','ParseSpallSettings':'PROD',
 'ParseVaporSettings':'BEAM','ParseBitSettings':'COLUMN','PrincipalStressRatio':'DAMAGE',
 'GatherColumnSurface':'COLUMN','ColumnSurfaceCos':'COLUMN','UpdateSpAfterMechanics':'COLUMN',
 'UpdateSpClusters':'COLUMN','LawSigma':'PROD','RockLaw':'PROD','BuildEnthalpyTable':'ZHANG',
 'SetPhaseFractionsConstAlpha':'ZHANG','ConvertTtoH':'PROD','InvertHtoT':'PROD','InitSurfaceMask':'PROD',
 'UpdateSurfaceMaskFromPhi':'COLUMN','ApplySeedVoid':'DIAG','InitPhi':'PROD',
 'RefreshMicrostructureDerivedFields':'MECH','RepairMechanicsAfterRegrid':'MECH','FiniteMatrix':'DAMAGE',
 'AlphaBetaFraction':'MECH','ThermalFreeStrain':'MECH','CellToNodeAverageInDomain':'MECH',
 'UpdateMicrostructureFreeStrain':'MECH','InitializeMicrostructure':'MICRO','InitializeFlaws':'MICRO',
 'enum class SolverKind':'ZHANG','WeibullSizeFactor':'COLUMN','struct SurfacePatch':'FLAME',
 'struct MineralPhase':'PROD','enum class MicroMode':'MICRO','enum class SpallSampleMode':'DAMAGE',
 'MeltSkeleton':'BEAM','removal_enabled':'PROD','enum class SpallationModel':'PROD',
 'enum class KIcLaw':'PROD','enum class DetachMode':'COLUMN','enum class SurfaceTMode':'PROD',
 'RawSurfaceX':'PROD','ReconFromX':'PROD','SurfaceGasTempCap':'PROD','SplitMix64':'MICRO',
 'ParseSpallationModel':'PROD'}
CLOSEDF = {'ColumnSurfaceCos', 'MeltSkeleton'}
rows = [l.rstrip('\n').split('\t') for l in open('docs/project/census_2026-11-05/work/p2_split_raw.tsv')]
for r in rows:
    if r[0] != '2': continue
    nm = r[1]
    if nm not in FUNC: print('UNMAPPED', nm, r[2]); continue
    put(int(r[9]), int(r[3]), FUNC[nm], nm in CLOSEDF)
# --- member / gap ranges and sub-block overrides (start, end, lineage, closed)
OV = [
 (1, 78, 'PROD'),
 # Parse sub-blocks
 (140, 149, 'BEAM'), (155, 161, 'ZHANG'), (174, 181, 'ZHANG'), (191, 202, 'DAMAGE'),
 (210, 231, 'COLUMN', True), (232, 237, 'BEAM'), (248, 251, 'COLUMN'), (252, 267, 'ZHANG'),
 (268, 270, 'MECH'), (284, 327, 'MICRO'), (393, 414, 'MECH'), (424, 437, 'MICRO'),
 (469, 479, 'MECH'), (483, 517, 'DAMAGE'), (543, 548, 'BEAM'), (561, 566, 'COLUMN'),
 (567, 572, 'MICRO'), (573, 577, 'COLUMN'), (579, 621, 'COLUMN'), (622, 622, 'DEAD'),
 (680, 700, 'COLUMN'), (701, 716, 'FLAME', True), (717, 740, 'FLAME'), (779, 786, 'DIAG'),
 (804, 827, 'DIAG'), (828, 904, 'FLAME'), (931, 940, 'COLUMN', True),
 (975, 992, 'FLAME'), (1096, 1111, 'FLAME'), (1112, 1125, 'BEAM'),
 # stranded comments
 (1401, 1404, 'DEAD'), (2846, 2853, 'DEAD'), (3064, 3068, 'DEAD'),
 # AdvanceMicrostructure: the beam path (solid-above gather, BeamSource call)
 (1833, 1841, 'BEAM'), (2035, 2061, 'BEAM'), (1891, 1897, 'FLAME'),
 # face-list members, then the path members
 (3286, 3302, 'PROD'), (3303, 3309, 'FLAME'), (3310, 3310, 'PROD'),
 # FaceListHeat's path branch
 (3858, 3864, 'FLAME'),
 # SurfaceCellFlux's pinned branch
 (4023, 4033, 'FLAME'),
 # BeamSource's header gap
 (4888, 4907, 'BEAM'),
 # TagCellsForRefinement: the base calls
 (4982, 4987, 'PROD'),
 # ParseSpallSettings sub-blocks
 (5447, 5501, 'COLUMN', True), (5502, 5526, 'BEAM', True), (5527, 5541, 'DIAG'),
 (5588, 5598, 'DEAD'), (5620, 5637, 'COLUMN', True), (5683, 5747, 'DAMAGE'),
 (5763, 5766, 'COLUMN', True), (5768, 5810, 'DIAG'),
 # UpdateSpAfterMechanics: the stress_mf (mechanics) branch
 (6061, 6062, 'MECH'), (6125, 6126, 'MECH'), (6270, 6276, 'MECH'),
 (6080, 6090, 'COLUMN', True), (6188, 6194, 'COLUMN', True),
 # the include of Removal.H / FaceRemoval.H
 (6617, 6623, 'PROD'),
 # ConvertTtoH / InvertHtoT global-material + const-alpha branches
 (6676, 6690, 'ZHANG'), (6724, 6757, 'ZHANG'),
 # InitializeMicrostructure: phase tables + property fill are production
 (7020, 7072, 'PROD'), (7190, 7286, 'PROD'), (7334, 7361, 'DAMAGE'),
 # members after InitializeFlaws
 (7547, 7566, 'PROD'), (7567, 7573, 'ZHANG'), (7574, 7574, 'ZHANG'), (7576, 7580, 'ZHANG'),
 (7581, 7582, 'MECH'), (7583, 7598, 'BEAM'), (7599, 7607, 'PROD'), (7608, 7613, 'DIAG'),
 (7614, 7681, 'FLAME'), (7682, 7714, 'PROD'), (7715, 7749, 'COLUMN'), (7760, 7761, 'FLAME'),
 (7762, 7773, 'DIAG'),
 # SurfacePatch: the shared parts are production
 (7774, 7809, 'PROD'), (7821, 7829, 'PROD'), (7897, 7911, 'PROD'), (7935, 7975, 'PROD'),
 (8013, 8048, 'PROD'), (8252, 8268, 'PROD'), (8358, 8359, 'PROD'),
 (8360, 8361, 'DEAD'),
 # members 8404..8978
 (8404, 8414, 'MICRO'), (8415, 8418, 'PROD'), (8419, 8434, 'MICRO'), (8435, 8442, 'PROD'),
 (8443, 8449, 'MECH'), (8450, 8501, 'DAMAGE'), (8502, 8506, 'PROD'), (8507, 8528, 'DAMAGE'),
 (8529, 8537, 'COLUMN', True), (8538, 8550, 'COLUMN'), (8551, 8557, 'COLUMN', True),
 (8558, 8562, 'DIAG'), (8563, 8592, 'PROD'), (8593, 8598, 'COLUMN', True),
 (8599, 8605, 'BEAM', True), (8619, 8622, 'COLUMN', True), (8623, 8638, 'COLUMN'),
 (8639, 8643, 'PROD'), (8644, 8659, 'BEAM'), (8660, 8675, 'COLUMN'), (8676, 8678, 'BEAM'),
 (8679, 8685, 'PROD'), (8687, 8687, 'PROD'), (8693, 8705, 'PROD'), (8716, 8721, 'PROD'),
 (8722, 8758, 'COLUMN'), (8759, 8763, 'PROD'), (8764, 8774, 'DEAD'), (8775, 8795, 'COLUMN', True),
 (8796, 8796, 'DEAD'), (8797, 8827, 'DIAG'), (8828, 8858, 'COLUMN'), (8859, 8880, 'MICRO'),
 (8889, 8914, 'COLUMN'), (8915, 8922, 'DIAG'), (8923, 8963, 'PROD'), (8965, 8978, 'PROD'),
 # ParseSpallationModel sub-blocks
 (9166, 9175, 'DAMAGE'), (9242, 9268, 'COLUMN'), (9269, 9317, 'MICRO'), (9318, 9341, 'COLUMN'),
 (9342, 9348, 'MICRO'), (9349, 9370, 'COLUMN'),
 (9373, N, 'PROD'),
 # closed FLAME switches (C.7.3): pad_shield_open (D2aa), jet_above_h/jet_T_above (D2v),
 # jet_T_ent_mode = rec_fixed (D2e)
 (986, 986, 'FLAME', True), (989, 989, 'FLAME', True), (2116, 2116, 'FLAME', True),
 (2151, 2155, 'FLAME', True), (2166, 2191, 'FLAME', True), (2205, 2222, 'FLAME', True),
 (2568, 2568, 'FLAME', True), (3671, 3672, 'FLAME', True), (3686, 3689, 'FLAME', True),
 (7816, 7820, 'FLAME', True), (7840, 7842, 'FLAME', True), (7850, 7860, 'FLAME', True),
 (7874, 7896, 'FLAME', True), (7999, 8012, 'FLAME', True), (8095, 8098, 'FLAME', True),
 (8107, 8118, 'FLAME', True), (8178, 8200, 'FLAME', True),
]
for o in OV:
    put(o[0], o[1], o[2], o[3] if len(o) > 3 else False)
# single blank/gap lines between functions take the lineage of the next assigned line
for x in range(N, 0, -1):
    if lin[x] is None:
        nx = x + 1
        while nx <= N and lin[nx] is None: nx += 1
        lin[x] = lin[nx] if nx <= N else 'PROD'
gaps = [x for x in range(1, N + 1) if lin[x] is None]
if gaps:
    # compress
    seg = []; s = gaps[0]; p = s
    for g in gaps[1:]:
        if g != p + 1: seg.append((s, p)); s = g
        p = g
    seg.append((s, p)); print('GAPS', seg)
tot = collections.Counter(); kd = collections.defaultdict(collections.Counter)
cl = collections.Counter()
for x in range(1, N + 1):
    tot[lin[x]] += 1; kd[lin[x]][K[x - 1]] += 1
    if closed[x]: cl[lin[x]] += 1
order = ['PROD','FLAME','BEAM','ZHANG','MICRO','MECH','DAMAGE','COLUMN','DIAG','DEAD']
print(f"{'lineage':8} {'lines':>6} {'code':>6} {'comment':>8} {'blank':>6} {'closed':>7}")
for L in order:
    print(f"{L:8} {tot[L]:6d} {kd[L]['code']:6d} {kd[L]['comment']:8d} {kd[L]['blank']:6d} {cl[L]:7d}")
print('TOTAL', sum(tot.values()), 'closed total', sum(cl.values()))
allK = collections.Counter(K); print('file code/comment/blank', allK)
# removal scenario
rem = 0
for x in range(1, N + 1):
    if lin[x] in ('DAMAGE', 'COLUMN', 'ZHANG', 'DEAD') or closed[x]: rem += 1
print('removed (DAMAGE+COLUMN+ZHANG+DEAD+closed)', rem, 'remaining', N - rem)
rem2 = sum(1 for x in range(1, N + 1) if lin[x] in ('DAMAGE','COLUMN','ZHANG') or closed[x])
print('removed without DEAD', rem2, 'remaining', N - rem2)
# dump per-line lineage for audit
with open('docs/project/census_2026-11-05/work/p2_lineage_by_line.tsv', 'w') as f:
    for x in range(1, N + 1):
        f.write(f"{x}\t{lin[x]}\t{int(closed[x])}\t{K[x-1]}\n")
