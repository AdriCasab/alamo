#!/usr/bin/env python3
"""P2-12 REVERSED (user, 2026-10-30): X1's surface-path mode is restored, default-off, as the physically preferred
structure for wall heating (it follows the gas along the surface; G1 maps heat by position in space). The sources
are restored byte-for-byte from ../src_p2_11b/ (the post-P2-11b working tree saved before the deletion), the lock
test likewise, and unit/face_heat is moved back.

Gates (no-change claims against P2-11b's own records; P2-10's comparators, output pointed HERE):
  r1  source == ../src_p2_11b (3 files)                                                  (chain_restore.sh)
  r1b the binary carries the X1 strings again                                           (chain_restore.sh)
  r2  tree 0 of 2639 vs p2_11b and p2_10_c2; controls d2aa 26, d2t 863                    (chain_restore.sh)
  r3  W            the witness set vs P2-10's W_c2                                       -> W_p2_12r_vs_c2.txt
  r4  UNITS        every kept unit test vs P2-11b's UNIT_p2_11b (face_heat included):
                   same rc, identical verdict lines, every test                          -> UNITS_p2_12r_vs_p2_11b.txt
  r5  X1REPRO      X1's own argv (studies/p2_11b_path/legs/X1) to 15 s: plotfiles 00000 and 01875 byte-identical
                   to the original X1 leg (AMReX metadata excluded); thermo rows to 15 s identical on every
                   non-timer column; event-log and jet-profile rows with time < 15 s identical  -> X1REPRO.txt
"""
import csv
import hashlib
import os
import shlex
import subprocess
import sys
import importlib.machinery
import importlib.util

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(os.path.dirname(HERE))
ROOT = os.path.abspath(os.path.join(STUD, "..", "..", ".."))
LABEL = "p2_12r"
TIMERS = {"fr_t_removal", "fr_t_step", "fl_t_path"}


def _load(path, name):
    ld = importlib.machinery.SourceFileLoader(name, path)
    m = importlib.util.module_from_spec(importlib.util.spec_from_loader(name, ld))
    ld.exec_module(m)
    return m


R10 = _load(os.path.join(STUD, "p2_10_retire", "run_p2_10.py"), "r10_p212r")
R10.HERE = HERE


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def W():
    R10.refs(LABEL)
    A = os.path.join(STUD, "p2_10_retire", "W_c2")
    B = os.path.join(HERE, "W_" + LABEL)
    lines = [f"r3 (W) -- W_{LABEL} (restored binary) vs P2-10's W_c2, P2-10's comparator and exclusions", "=" * 100]
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
    A = os.path.join(STUD, "p2_11b_path", "UNIT_p2_11b")
    B = os.path.join(HERE, "UNIT_" + LABEL)
    rc = lambda d: {l.split()[0]: l.split()[1] for l in open(os.path.join(d, "rc.txt")).read().splitlines()[1:]}
    ra, rb = rc(A), rc(B)
    L = ["r4 -- every kept unit test on the restored binary vs P2-11b's UNIT_p2_11b: rc and verdict lines", "=" * 100]
    bad = []
    for n in sorted(set(ra) | set(rb)):
        va = R10.verdict_lines(os.path.join(A, n + ".log")) if n in ra else None
        vb = R10.verdict_lines(os.path.join(B, n + ".log")) if n in rb else None
        same = ra.get(n) == rb.get(n) and va == vb
        L.append(f"  {n:28s} rc {ra.get(n, '-'):>5s} -> {rb.get(n, '-'):>5s}   verdict lines "
                 f"{'identical (%d)' % len(vb) if same else 'DIFFER'}")
        if not same:
            bad.append(n)
    L.append(f"r4 (every test, face_heat included): {'PASS' if not bad else 'FAIL: ' + ', '.join(bad)}")
    txt = "\n".join(L)
    print(txt)
    open(os.path.join(HERE, f"UNITS_{LABEL}_vs_p2_11b.txt"), "w").write(txt + "\n")
    return not bad


def _rows_before(path, tmax, tcol="time"):
    rows = list(csv.reader(open(path)))
    h = rows[0]
    i = h.index(tcol)
    return h, [r for r in rows[1:] if r and float(r[i]) < tmax]


def X1REPRO(run=True):
    O = os.path.join(STUD, "p2_11b_path", "legs", "X1")
    d = os.path.join(HERE, "X1_15")
    os.makedirs(d, exist_ok=True)
    toks = shlex.split(open(os.path.join(O, "argv.sh")).read())
    out = []
    for t in toks:
        k = t.split("=", 1)[0]
        if k == "plot_file":
            out.append("plot_file=" + os.path.join(d, "run"))
        elif k == "stop_time":
            out.append("stop_time=15.0")
        else:
            out.append(t)
    open(os.path.join(d, "argv.sh"), "w").write(" ".join(shlex.quote(c) for c in out) + "\n")
    if run:
        with open(os.path.join(d, "run.log"), "w") as f:
            rc = subprocess.call(out, cwd=ROOT, stdout=f, stderr=subprocess.STDOUT)
    else:
        rc = 0 if "STEP" in open(os.path.join(d, "run.log")).read() else 1   # compare-only on an existing rerun
    L = ["r5 -- X1's own argv on the restored binary to 15 s vs the original X1 leg (P2-11b binary f787140d)",
         "=" * 100, f"  rc {rc}"]
    ok = rc == 0
    # plotfiles
    for pf in ("00000cell", "00000node", "01875cell", "01875node"):
        a, b = os.path.join(O, "run", pf), os.path.join(d, "run", pf)
        if not os.path.isdir(b):
            L.append(f"  {pf}: MISSING in the rerun")
            ok = False
            continue
        fa = {os.path.relpath(os.path.join(dp, f), a): os.path.join(dp, f) for dp, _, fs in os.walk(a) for f in fs}
        fb = {os.path.relpath(os.path.join(dp, f), b): os.path.join(dp, f) for dp, _, fs in os.walk(b) for f in fs}
        diff = [k for k in set(fa) & set(fb) if md5(fa[k]) != md5(fb[k])]
        only = sorted(set(fa) ^ set(fb))
        L.append(f"  {pf}: {len(set(fa) & set(fb))} files, {len(diff)} differ, {len(only)} only on one side")
        ok = ok and not diff and not only
    # thermo
    ta = open(os.path.join(O, "run", "thermo.dat")).read().splitlines()
    tb = open(os.path.join(d, "run", "thermo.dat")).read().splitlines()
    ha = [c.strip() for c in ta[0].split("\t")]
    hb = [c.strip() for c in tb[0].split("\t")]
    ra = [r.split("\t") for r in ta[1:] if r.strip() and float(r.split("\t")[0]) <= 15.0 + 1e-9]
    rb = [r.split("\t") for r in tb[1:] if r.strip() and float(r.split("\t")[0]) <= 15.0 + 1e-9]
    cols = [c for c in ha if c not in TIMERS]
    nbad = 0
    if ha != hb or len(ra) != len(rb):
        L.append(f"  thermo: header same {ha == hb}, rows {len(ra)} vs {len(rb)}")
        ok = False
    else:
        for x, y in zip(ra, rb):
            for c in cols:
                i = ha.index(c)
                if x[i] != y[i]:
                    nbad += 1
        L.append(f"  thermo: {len(ra)} rows to 15 s, {len(cols)} non-timer columns, {nbad} differing cells")
        ok = ok and nbad == 0
    # event log, jet profile. ⚠ DECLARED FIX (first run, CHAIN.log 19:11): the window "time < 15 s" counted
    # 2 original rows stamped 14.999999999999234 -- the step that STARTS at t = 15 s (the log stamps a step's
    # start, SR 7(iii)), which a run stopped at 15.0 never executes. The claim is therefore PREFIX identity:
    # the rerun's file equals the original's first len(rerun) rows, and the rerun ends at the last step
    # before the 15 s plotfile.
    for fn in ("run_removal_events.csv", "run_jet_profile.csv"):
        pa, pb = os.path.join(O, fn), os.path.join(d, fn)
        xa = list(csv.reader(open(pa)))
        xb = list(csv.reader(open(pb)))
        prefix = xa[:len(xb)] == xb
        t_last = float(xb[-1][0]) if len(xb) > 1 else -1.0
        # the event log writes every step (last = 14.992 s); the jet profile writes once per
        # jet_profile_interval = 1 s (last = the first step of the 14 s interval) -- second declared fix
        lo = 14.99 if fn == "run_removal_events.csv" else 14.0
        ends = lo <= t_last < 15.0 - 1e-6
        L.append(f"  {fn}: rerun {len(xb) - 1} rows, an exact prefix of the original's {len(xa) - 1}: {prefix}; "
                 f"last rerun step stamped {t_last:.6f} s (must be the step before 15 s): {ends}")
        ok = ok and prefix and ends
    L.append(f"r5: {'PASS' if ok else 'FAIL'}")
    txt = "\n".join(L)
    print(txt)
    open(os.path.join(HERE, "X1REPRO.txt"), "w").write(txt + "\n")
    return ok


if __name__ == "__main__":
    ok = {"W": W, "UNITS": UNITS, "X1REPRO": X1REPRO, "X1CMP": lambda: X1REPRO(run=False)}[sys.argv[1]]()
    sys.exit(0 if ok else 1)
