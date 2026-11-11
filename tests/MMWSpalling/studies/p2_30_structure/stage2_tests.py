#!/usr/bin/env python3
"""P2-30 stage 2 (deliverable 5): the lock trim of unit/lateral_flux_lock (S/LOCK_TRIM.txt carries the reasons), and the one
pointer line in the parent's head. Every replacement must match exactly once or nothing is written."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
LFL = os.path.join(ROOT, "tests", "MMWSpalling", "unit", "lateral_flux_lock", "test")
PARENT = os.path.join(ROOT, "src", "Integrator", "MMWSpalling.H")

CUTS = [
    ('''    # P2-28 (R11's twin): no column-only name survives anywhere in src/ (comments included)
    src_all = "".join(t for r, t in _src.tree().items() if r.endswith((".H", ".cpp", ".cc")))
    COL = ["UpdateRemovalAfterCohesive", "BuildSurfaceTColumns", "UpdateSpAfterMechanics", "UpdateSpClusters",
           "PatchDiagnostics", "GatherColumnSurface", "ColumnSolidRunBelowTop", "PinRule", "BuildPinEffective",
           "SurfaceRobinT", "SharedSurfaceRecon", "RepairRemovalStateAfterRegrid", "CheckSurfaceCovered",
           "rc_ell_prev", "DetachMode", "E_sup_local"]
    left = [(n, len(re.findall(r"(?<![A-Za-z_])%s" % n, src_all))) for n in COL
            if re.search(r"(?<![A-Za-z_])%s" % n, src_all)]
    need(not left and not _src.exists(REM),
         "the column path is GONE from src/: Removal.H deleted, %d column-only names left %s" % (len(left), left))
''', ''),
    ('''        # P2-27 (census C1b): the three fields deregistered
        "C1b": r"\\bD_mf\\b|\\bkappa_eff_mf\\b|\\bk_eff_mf\\b|\\bic_damage\\b|\\bRefreshKappaEff\\b|\\bkeff_out\\b",
''', ''),
    ('''    nlist = sum(_src.header().count(k) for k in ("\\"spall.lateral_drive\\"", "\\"surface_patch.side_face_flux\\"",
                                    "\\"surface_patch.jet_wall_march\\""))
    need(nlist == 3, "the retired-key list itself exists and names them (%d of 3 probes)" % nlist)
''', ''),
    ('''    # P2-10: P2-8's face Newton (F4) had no caller after P2-9 and is deleted.
    need(ndecl(r"void\\s+FaceListFlux\\s*\\(") == 0,
         "`FaceListFlux` is DELETED (P2-10; FrFaceFlux heats every listed face)")
''', ''),
]


def main():
    t = open(LFL).read()
    for old, new in CUTS:
        n = t.count(old)
        if n != 1:
            sys.exit("expected 1 match, found %d of:\n%s" % (n, old[:160]))
        t = t.replace(old, new, 1)
    p = open(PARENT).read()
    head = "// See in-main-tex-you-will-quizzical-treasure.md.\n"
    if p.count(head) != 1:
        sys.exit("parent head line not unique")
    p = p.replace(head, head + "// Dated history moved out of these sources (P2-30): docs/project/CODE_HISTORY.md.\n", 1)
    open(LFL, "w").write(t)
    open(PARENT, "w").write(p)
    print("lateral_flux_lock: %d clauses trimmed (S/LOCK_TRIM.txt); parent head: one pointer line" % len(CUTS))


if __name__ == "__main__":
    main()
