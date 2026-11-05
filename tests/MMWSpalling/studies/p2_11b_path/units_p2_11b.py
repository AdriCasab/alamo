#!/usr/bin/env python3
"""P2-11b: every kept unit test on the Stage B binary (P2-10's `units` harness, output pointed HERE), compared with
P2-10's UNIT_c2 (binary C2): same rc, identical verdict lines. Expected differences: lateral_flux_lock (its FrFaceFlux
count restated 4 -> 5) and the new face_heat (no C2 row). Writes UNITS_p2_11b_vs_c2.txt."""
import importlib.machinery
import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)
ld = importlib.machinery.SourceFileLoader("r10u", os.path.join(STUD, "p2_10_retire", "run_p2_10.py"))
R10 = importlib.util.module_from_spec(importlib.util.spec_from_loader("r10u", ld))
ld.exec_module(R10)
C2 = os.path.join(STUD, "p2_10_retire", "UNIT_c2")
R10.HERE = HERE
R10.units("p2_11b")
B = os.path.join(HERE, "UNIT_p2_11b")
rc = lambda d: {l.split()[0]: l.split()[1] for l in open(os.path.join(d, "rc.txt")).read().splitlines()[1:]}
ra, rb = rc(C2), rc(B)
L = ["P2-11b -- kept unit tests on the Stage B binary vs P2-10's UNIT_c2: rc and verdict lines", "=" * 100]
bad = []
for n in sorted(set(ra) | set(rb)):
    va = R10.verdict_lines(os.path.join(C2, n + ".log")) if n in ra else None
    vb = R10.verdict_lines(os.path.join(B, n + ".log")) if n in rb else None
    same = ra.get(n) == rb.get(n) and va == vb
    L.append(f"  {n:28s} rc {ra.get(n, '-'):>5s} -> {rb.get(n, '-'):>5s}   verdict lines {'identical' if same else 'DIFFER'}")
    if not same and n not in ("lateral_flux_lock", "face_heat"):
        bad.append(n)
L.append(f"kept units (lateral_flux_lock, face_heat expected to differ): {'PASS' if not bad else 'FAIL: ' + ', '.join(bad)}")
txt = "\n".join(L)
print(txt)
open(os.path.join(HERE, "UNITS_p2_11b_vs_c2.txt"), "w").write(txt + "\n")
