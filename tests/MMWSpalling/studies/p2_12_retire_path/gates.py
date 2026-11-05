#!/usr/bin/env python3
"""P2-12 -- retire P2-11b's surface-path heat mode (surface_patch.face_heat = path, barrier_expr).

The deletion restores src/Integrator/MMWSpalling.H and MMWSpalling/FaceRemoval.H to P2-11b's own
pre-edit snapshot (studies/p2_11b_path/src_pre/, the P2-10 C2 source), restores the two lines of
unit/lateral_flux_lock/test that P2-11b restated (FrFaceFlux count 5 -> 4), and moves unit/face_heat/
into retired_units/. The pre-deletion files are kept in src_p2_11b/.

Gates (every one is a no-change claim; P2-10's comparators and exclusions, loaded as a module with
the output folder pointed HERE so no P2-10 or P2-11b record is touched):
  d1  source: the two headers byte-identical to studies/p2_11b_path/src_pre/          (chain.sh)
  d2  tree sweep 0 of 2639 vs p2_10_c2 AND vs p2_11b; controls d2aa 26, d2t 863       (chain.sh)
  d3  gates.py W        the witness set W on this binary vs P2-10's W_c2               -> W_p2_12_vs_c2.txt
  d4  gates.py UNITS    every kept unit test vs P2-10's UNIT_c2: same rc, identical
                        verdict lines, for EVERY test (face_heat absent on both)        -> UNITS_p2_12_vs_c2.txt
  d5  gates.py RETIRED  X1's own argv (studies/p2_11b_path/legs/X1) with stop_time 1e-9 must ABORT
                        and name surface_patch.face_heat (the strict surface_patch parse) -> RETIRED.txt
"""
import os
import shlex
import subprocess
import sys
import importlib.machinery
import importlib.util

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)
ROOT = os.path.abspath(os.path.join(STUD, "..", "..", ".."))
LABEL = "p2_12"


def _load(path, name):
    ld = importlib.machinery.SourceFileLoader(name, path)
    m = importlib.util.module_from_spec(importlib.util.spec_from_loader(name, ld))
    ld.exec_module(m)
    return m


R10 = _load(os.path.join(STUD, "p2_10_retire", "run_p2_10.py"), "r10_p212")
R10.HERE = HERE                     # every output of the harness lands in this study


def W():
    R10.refs(LABEL)
    A = os.path.join(STUD, "p2_10_retire", "W_c2")
    B = os.path.join(HERE, "W_" + LABEL)
    lines = [f"d3 (W) -- W_{LABEL} (this binary) vs P2-10's W_c2 (binary C2), P2-10's comparator and exclusions "
             "(metadata, diff.patch, run.log, argv/BINARY/DROPPED/.done; thermo timers)", "=" * 100]
    ok = True
    for sub in ["M0T"] + ["fr/" + x for x in R10.FR_LEGS] + ["fl/" + x for x in R10.FL_LEGS]:
        o, th = R10.compare_dirs(os.path.join(A, sub), os.path.join(B, sub), lines, sub)
        ok = ok and o
        for rel in th:
            c = R10.thermo_cmp(os.path.join(A, sub, rel), os.path.join(B, sub, rel))
            ok = ok and not c["only_a"]
    lines.append(f"W_{LABEL} vs W_c2: {'PASS' if ok else 'FAIL'}")
    txt = "\n".join(lines)
    print(txt)
    open(os.path.join(HERE, f"W_{LABEL}_vs_c2.txt"), "w").write(txt + "\n")
    return ok


def UNITS():
    R10.units(LABEL)
    A = os.path.join(STUD, "p2_10_retire", "UNIT_c2")
    B = os.path.join(HERE, "UNIT_" + LABEL)
    rc = lambda d: {l.split()[0]: l.split()[1] for l in open(os.path.join(d, "rc.txt")).read().splitlines()[1:]}
    ra, rb = rc(A), rc(B)
    L = [f"d4 -- every kept unit test on this binary vs P2-10's UNIT_c2 (binary C2): rc and verdict lines", "=" * 100]
    bad = []
    for n in sorted(set(ra) | set(rb)):
        va = R10.verdict_lines(os.path.join(A, n + ".log")) if n in ra else None
        vb = R10.verdict_lines(os.path.join(B, n + ".log")) if n in rb else None
        same = ra.get(n) == rb.get(n) and va == vb
        L.append(f"  {n:28s} rc {ra.get(n, '-'):>5s} -> {rb.get(n, '-'):>5s}   verdict lines "
                 f"{'identical (%d)' % len(vb) if same else 'DIFFER'}")
        if not same:
            bad.append(n)
    L.append(f"  face_heat present on this binary: {'face_heat' in rb} (must be False)")
    if "face_heat" in rb:
        bad.append("face_heat present")
    L.append(f"d4 (every test, no exception): {'PASS' if not bad else 'FAIL: ' + ', '.join(bad)}")
    txt = "\n".join(L)
    print(txt)
    open(os.path.join(HERE, f"UNITS_{LABEL}_vs_c2.txt"), "w").write(txt + "\n")
    return not bad


def RETIRED():
    src = os.path.join(STUD, "p2_11b_path", "legs", "X1", "argv.sh")
    toks = shlex.split(open(src).read())
    d = os.path.join(HERE, "smoke", "X1_retired")
    os.makedirs(d, exist_ok=True)
    out = []
    for t in toks:
        k = t.split("=", 1)[0]
        if k == "plot_file":
            out.append("plot_file=" + os.path.join(d, "run"))
        elif k == "stop_time":
            out.append("stop_time=1e-9")
        else:
            out.append(t)
    open(os.path.join(d, "argv.sh"), "w").write(" ".join(shlex.quote(c) for c in out) + "\n")
    with open(os.path.join(d, "run.log"), "w") as f:
        rc = subprocess.call(out, cwd=ROOT, stdout=f, stderr=subprocess.STDOUT)
    log = open(os.path.join(d, "run.log"), errors="replace").read()
    named = [l.strip() for l in log.splitlines() if "face_heat" in l or "barrier_expr" in l]
    ok = rc != 0 and any("face_heat" in l for l in named)
    L = ["d5 -- X1's argv (face_heat=path, barrier_expr) on this binary, stop_time 1e-9", "=" * 100,
         f"  source argv: {os.path.relpath(src, ROOT)}", f"  rc {rc}",
         "  log lines naming the retired keys:"] + ["    " + l[:200] for l in named[:10]] + \
        [f"d5: {'PASS' if ok else 'FAIL'} (must abort and name surface_patch.face_heat)"]
    txt = "\n".join(L)
    print(txt)
    open(os.path.join(HERE, "RETIRED.txt"), "w").write(txt + "\n")
    return ok


if __name__ == "__main__":
    ok = {"W": W, "UNITS": UNITS, "RETIRED": RETIRED}[sys.argv[1]]()
    sys.exit(0 if ok else 1)
