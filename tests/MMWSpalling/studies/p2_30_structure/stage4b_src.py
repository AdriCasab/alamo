#!/usr/bin/env python3
"""P2-30 stage 4b (deliverable 11): the three true duplicates of the surface patch set-up (inv_dz2, sp_presc, sp_r2) come
from one function, read by the kernel and by EnergyLedger; the two "was I handed column tops" locals are renamed have_tops;
L7's lateral_flux_lock clause. `sp_on` and the list_state() sites are NOT touched; no flux expression is rewritten.
Every replacement must match exactly once or nothing is written."""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
D = os.path.join(ROOT, "src", "Integrator", "MMWSpalling")
LFL = os.path.join(ROOT, "tests", "MMWSpalling", "unit", "lateral_flux_lock", "test")

SHARED_ANCHOR = """    // =====================================================================
    // P2-18 A3: THE STEP PROLOGUE both paths share (formerly the first lines
"""
SHARED_NEW = """    // P2-30 (SR 14): the surface patch's three constants, ONE definition, read by the kernel and by EnergyLedger.
    // ⚠ `sp_on` is NOT here: the kernel's carries `owner` and the ledger's does not (equal on one level only).
    struct PatchConst { Set::Scalar inv_dz2; bool presc; Set::Scalar r2; };
    PatchConst PatchConstants(const Set::Scalar* DX) const
    {
        PatchConst p;
        p.inv_dz2 = 1.0 / (DX[2] * DX[2]);
        p.presc   = surface_patch.mode == SurfacePatch::Mode::PrescribedTemperature;
        p.r2      = surface_patch.radius * surface_patch.radius;
        return p;
    }

"""

EDITS = {
    "Shared.H": [(SHARED_ANCHOR, SHARED_NEW + SHARED_ANCHOR)],
    "Heat.H": [
        ("        const Set::Scalar inv_dz2 = 1.0 / (DX[2] * DX[2]);\n",
         "        const PatchConst patch_c = PatchConstants(DX);\n"
         "        const Set::Scalar inv_dz2 = patch_c.inv_dz2;\n"),
        ("""        const bool        sp_presc = surface_patch.mode ==
                                     SurfacePatch::Mode::PrescribedTemperature;
""", "        const bool        sp_presc = patch_c.presc;\n"),
        ("        const Set::Scalar sp_r2    = surface_patch.radius * surface_patch.radius;\n",
         "        const Set::Scalar sp_r2    = patch_c.r2;\n"),
    ],
    "Ledger.H": [
        ("        const Set::Scalar inv_dz2 = 1.0 / (dz * dz);\n",
         "        const PatchConst patch_c = PatchConstants(DX);\n"
         "        const Set::Scalar inv_dz2 = patch_c.inv_dz2;\n"),
        ("""        const bool        sp_presc = surface_patch.mode ==
                                     SurfacePatch::Mode::PrescribedTemperature;
""", "        const bool        sp_presc = patch_c.presc;\n"),
        ("        const Set::Scalar sp_r2    = surface_patch.radius * surface_patch.radius;\n",
         "        const Set::Scalar sp_r2    = patch_c.r2;\n"),
    ],
}

L7_ANCHOR = """    # ---- P2-10: NOTHING DELETED SURVIVES ------------------------------------
"""
L7_CLAUSE = '''    # ---- P2-30 (L7, SR 14): THE SURFACE PATCH'S THREE CONSTANTS, ONE DEFINITION ---------------------------
    print("\\nP2-30 (L7) THE PATCH CONSTANTS -- one definition, read by the kernel and the ledger")
    code_l7 = "\\n".join(l.split("//")[0] for l in _src.text().split("\\n"))
    if "--plant-patch-def" in sys.argv:          # the control: a second definition must turn the clause RED
        code_l7 += "\\n    PatchConst PatchConstants(const Set::Scalar* DX) const { return PatchConst(); }\\n"
    n_pc = len(re.findall(r"PatchConst\\s+PatchConstants\\s*\\(", code_l7))
    adv_l7 = re.search(r"void\\s+AdvanceMicrostructure\\s*\\(.*?\\n    \\}\\n", code_l7, re.S)
    led_l7 = re.search(r"void\\s+EnergyLedger\\s*\\(.*?\\n    \\}\\n", code_l7, re.S)
    reads_l7 = lambda b: (b is not None and "PatchConstants(DX)" in b.group(0) and "= patch_c.inv_dz2" in b.group(0)
                          and "= patch_c.presc" in b.group(0) and "= patch_c.r2" in b.group(0))
    n_dz2 = len(re.findall(r"\\binv_dz2\\s*=\\s*1\\.0\\s*/", code_l7))
    n_r2 = len(re.findall(r"\\bsp_r2\\s*=\\s*surface_patch\\.radius", code_l7))
    n_prs = len(re.findall(r"\\bsp_presc\\s*=\\s*surface_patch\\.mode", code_l7))
    need(n_pc == 1 and reads_l7(adv_l7) and reads_l7(led_l7) and n_dz2 == 1 and n_r2 == 0 and n_prs == 0,
         "L7 `PatchConstants` defined once (%d) and read by the kernel and EnergyLedger; inv_dz2 formed once (%d), "
         "sp_r2 and sp_presc formed nowhere else (%d, %d)" % (n_pc, n_dz2, n_r2, n_prs))

'''


def main():
    new = {}
    for fn, reps in EDITS.items():
        p = os.path.join(D, fn)
        t = open(p).read()
        for old, rep in reps:
            n = t.count(old)
            if n != 1:
                sys.exit("%s: expected 1 match, found %d of:\n%s" % (fn, n, old[:200]))
            t = t.replace(old, rep, 1)
        new[p] = t
    # the two "was I handed column tops" locals -> have_tops, each through the end of its own function; the list_state()
    # void_aware (BeamEnergyBalance) is a different meaning and stays
    p = os.path.join(D, "BeamProps.H")
    L = open(p).read().split("\n")
    defs = [i for i, s in enumerate(L) if "const bool void_aware = !k_top_col.empty();" in s]
    if len(defs) != 2:
        sys.exit("BeamProps.H: expected 2 definitions of the handed-tops void_aware, found %d" % len(defs))
    total = 0
    for i0 in defs:
        i1 = next(j for j in range(i0, len(L)) if L[j] == "    }")
        for j in range(i0, i1 + 1):
            L[j], n = re.subn(r"\bvoid_aware\b", "have_tops", L[j])
            total += n
    new[p] = "\n".join(L)
    left = len(re.findall(r"\bvoid_aware\b", new[p]))
    # the lock
    t = open(LFL).read()
    if t.count(L7_ANCHOR) != 1:
        sys.exit("lateral_flux_lock anchor")
    new[LFL] = t.replace(L7_ANCHOR, L7_CLAUSE + L7_ANCHOR, 1)
    for k, v in new.items():
        open(k, "w").write(v)
    print("4b applied: %s; have_tops renames %d (the two definitions and their uses); void_aware left in BeamProps.H (the list_state one): %d" % (
        ", ".join(sorted(EDITS)), total, left))


if __name__ == "__main__":
    main()
