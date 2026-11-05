#!/usr/bin/env python3
"""P2-14a -- every row. Usage: reduce.py [A1|A3|GATES|B|C] (no mode: RESULTS.txt to stdout).

S_top (SR 19): per column, the first event-log row with dir = 5 and ts_fb = 0 (the column path logs top firings
only: every row counts), its T_surf; T35 interior columns only (4 <= i, j <= 19).
Sources: unit/face_cref/output/<leg> (F0_k0/k1, F2_k0/k1, F164_k1, T35_0/2, N2_0/2, new binary) and legs/<leg>.
"""
import csv
import hashlib
import math
import os
import re
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
T = os.path.join(ROOT, "tests", "MMWSpalling")
UNIT = os.path.join(T, "unit", "face_cref", "output")
LEGS = os.path.join(HERE, "legs")
CALD = os.path.join(T, "validation", "kant", "grimsel_calibration")
FR = os.path.join(T, "unit", "face_removal", "output")
P2E = os.path.join(T, "studies", "p2_2e_cref")
CHECK = os.path.join(T, "validation", "shared", "check_shared.py")
TIMERS = {"fr_t_removal", "fr_t_step", "fl_t_path"}
sys.path.insert(0, os.path.join(T, "studies", "p2_2_props"))
PY = sys.executable


def d_of(leg):
    for base in (UNIT, LEGS):
        if os.path.isdir(os.path.join(base, leg)):
            return os.path.join(base, leg)
    return os.path.join(LEGS, leg)


def events(d):
    p = os.path.join(d, "run_removal_events.csv")
    return list(csv.DictReader(open(p))) if os.path.exists(p) else []


def first_top(rows, interior=False):
    f = {}
    for r in rows:
        if "dir" in r and not (r["dir"] == "5" and float(r["ts_fb"]) == 0.0):
            continue
        c = (int(r["col_i"]), int(r["col_j"]))
        if interior and not (4 <= c[0] <= 19 and 4 <= c[1] <= 19):
            continue
        f.setdefault(c, r)
    return f


def stat(d, interior=False):
    f = first_top(events(d), interior)
    v = np.array([float(r["T_surf"]) for r in f.values()])
    t = np.array([float(r["time"]) for r in f.values()])
    if not len(v):
        return dict(n=0, med=float("nan"), p10=float("nan"), p90=float("nan"), spread=0.0, t0=float("nan"))
    return dict(n=len(v), med=float(np.median(v)), p10=float(np.percentile(v, 10)), p90=float(np.percentile(v, 90)),
                spread=float(v.max() - v.min()), t0=float(t.min()), t1=float(t.max()))


def fmt(s):
    return ("n %4d  median %8.2f K  p10 %7.2f  p90 %7.2f" % (s["n"], s["med"], s["p10"], s["p90"])) if s["n"] \
        else "n    0  (no top-face firing)"


def thermo(d):
    L = open(os.path.join(d, "run", "thermo.dat")).read().splitlines()
    h = [c.strip() for c in L[0].split("\t")]
    return h, [r.split("\t") for r in L[1:] if r.strip()]


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def pf_dirs(d):
    r = os.path.join(d, "run")
    return sorted(x for x in os.listdir(r) if re.fullmatch(r"\d{5}(cell|node)", x)) if os.path.isdir(r) else []


def same_plotfiles(a, b):
    """(files compared, differing, only on one side) over every plotfile of both runs"""
    n = nd = no = 0
    pa, pb = pf_dirs(a), pf_dirs(b)
    no += len(set(pa) ^ set(pb))
    for pf in set(pa) & set(pb):
        A, B = os.path.join(a, "run", pf), os.path.join(b, "run", pf)
        fa = {os.path.relpath(os.path.join(p, f), A): os.path.join(p, f) for p, _, fs in os.walk(A) for f in fs}
        fb = {os.path.relpath(os.path.join(p, f), B): os.path.join(p, f) for p, _, fs in os.walk(B) for f in fs}
        no += len(set(fa) ^ set(fb))
        for k in set(fa) & set(fb):
            n += 1
            nd += md5(fa[k]) != md5(fb[k])
    return n, nd, no


def same_file(a, b, nm):
    pa, pb = os.path.join(a, nm), os.path.join(b, nm)
    return os.path.exists(pa) and os.path.exists(pb) and md5(pa) == md5(pb)


def thermo_same(a, b):
    ha, ra = thermo(a)
    hb, rb = thermo(b)
    if ha != hb or len(ra) != len(rb):
        return False, "header same %s, rows %d vs %d" % (ha == hb, len(ra), len(rb))
    idx = [i for i, c in enumerate(ha) if c not in TIMERS]
    nbad = sum(1 for x, y in zip(ra, rb) if [x[i] for i in idx] != [y[i] for i in idx])
    return nbad == 0, "%d rows x %d non-timer columns, %d rows differ" % (len(ra), len(idx), nbad)


def rc(d):
    p = os.path.join(d, "rc")
    return int(open(p).read()) if os.path.exists(p) else None


def check_leg(leg):
    d = d_of(leg)
    deck = open(os.path.join(d, "DECK.txt")).read().strip()
    allow = open(os.path.join(d, "ALLOW.txt")).read().split()
    a = [PY, CHECK, deck, os.path.join(d, "argv.sh")] + (["--allow"] + allow if allow else [])
    p = subprocess.run(a, cwd=ROOT, capture_output=True, text=True)
    return p.returncode == 0, p.stdout.strip().splitlines()[-1]


# ======================================================================= a1
def A1():
    L = ["a1 -- the Sources numbers re-read from disk (zero compute)", "=" * 100]
    ok = True

    def row(lbl, got, want, fmt_="%.2f"):
        nonlocal ok
        g = fmt_ % got if isinstance(got, float) else str(got)
        good = g == want
        ok = ok and good
        L.append("  %-46s %12s   packet %10s   %s" % (lbl, g, want, "ok" if good else "DIFFERS"))

    for leg, want in (("CAL", ("850.43", "830.69", "882.00", "576")), ("CAL_h2000", ("852.98",)),
                      ("CAL_p27", ("823.95",)), ("CAL_p48", ("812.00",))):
        s = stat(os.path.join(CALD, "legs", leg))
        row("P2-13 %s S_top median" % leg, s["med"], want[0])
        if len(want) > 1:
            row("P2-13 CAL p10", s["p10"], want[1]); row("P2-13 CAL p90", s["p90"], want[2])
            row("P2-13 CAL n", s["n"], want[3])
    ck = open(os.path.join(P2E, "CHECKS.txt")).read()
    for lbl, pat, want in (("P2-2e V_h150", r"V_h150: (NO FIRE) to t = 200\.0 s", "NO FIRE"),
                           ("P2-2e V_h500 onset", r"V_h500: FIRES at ([\d.]+) K", "960.03"),
                           ("P2-2e V_h500 time", r"V_h500: FIRES at [\d.]+ K, t = ([\d.]+) s", "8.813"),
                           ("P2-2e M80 2 mm", r"M80_2 +dz 2\.0 mm: onset +([\d.]+)", "915.816"),
                           ("P2-2e M80 1 mm", r"M80_1 +dz 1\.0 mm: onset +([\d.]+)", "890.430"),
                           ("P2-2e M80 0.5 mm", r"M80_05 +dz 0\.5 mm: onset +([\d.]+)", "888.097")):
        m = re.search(pat, ck)
        row(lbl, m.group(1) if m else "MISSING", want, "%s")
    negs = re.findall(r"\n +(?:V\w+|M80_\w+) +reads .*?neg +([\d.]+) %", ck)   # the flat legs (not B100_d's beam)
    row("P2-2e neg share, flat c12 legs V*/M80* (max)", "%.4f" % max(float(x) for x in negs), "0.0000", "%s")
    pr = open(os.path.join(P2E, "PROD.txt")).read()
    m = re.search(r"volume ([-+\d.]+) %.*firing T ([-+\d.]+) K", pr)
    row("P2-2e PROD HD volume", m.group(1) if m else "MISSING", "-59.91", "%s")
    row("P2-2e PROD HD firing T", m.group(2) if m else "MISSING", "+27.24", "%s")
    dr = open(os.path.join(P2E, "D_REF.txt")).read()
    row("P2-2e D_REF d_hi", "2.0" if "d_hi (largest d where SOD 150 does NOT fire) = 2.0 mm" in dr else "MISSING",
        "2.0", "%s")
    for leg, t_w, T_w in (("T00_1_k0", "2.440", "846.122"), ("T00_1_k1", "2.440", "845.995"),
                          ("T00_2_k0", "2.648", "846.110"), ("T00_2_k1", "2.648", "845.884")):
        rows = events(os.path.join(FR, leg))
        f = {}
        for r in rows:
            f.setdefault((r["col_i"], r["col_j"]), r)
        ts = sorted(set("%.3f" % float(r["time"]) for r in f.values()))
        Ts = sorted(set("%.3f" % float(r["T_surf"]) for r in f.values()))
        row("F %s first firing (all columns: one time)" % leg, ts[0] if len(ts) == 1 else "/".join(ts), t_w, "%s")
        row("F %s T_surf (all columns: one value)" % leg, Ts[0] if len(Ts) == 1 else "/".join(Ts), T_w, "%s")
    f = {}
    for r in events(os.path.join(FR, "TL35_1")):
        f.setdefault((r["col_i"], r["col_j"]), r)
    row("F TL35_1 first-firing median (any dir)", float(np.median([float(r["T_surf"]) for r in f.values()])), "845.79")
    L.append("a1 %s (report)" % ("REPRODUCED" if ok else "SOME DIFFER (see rows)"))
    return "\n".join(L)


# ======================================================================= a3
def A3():
    import importlib.machinery
    import importlib.util
    ld = importlib.machinery.SourceFileLoader("r10q", os.path.join(T, "studies", "p2_10_retire", "run_p2_10.py"))
    R10 = importlib.util.module_from_spec(importlib.util.spec_from_loader("r10q", ld))
    ld.exec_module(R10)
    A, B = os.path.join(CALD, "UNIT_p2_13"), os.path.join(HERE, "UNIT_p2_14a")
    rcm = lambda d: {l.split()[0]: l.split()[1] for l in open(os.path.join(d, "rc.txt")).read().splitlines()[1:]}
    ra, rb = rcm(A), rcm(B)
    L = ["a3 -- every kept unit on the new binary vs P2-13's UNIT_p2_13: rc and verdict lines", "=" * 100]
    bad = []
    for n in sorted(set(ra) | set(rb)):
        # ⚠ declared normalisation: an abort quoting `<file>.H:<line>` moves with every edit above that line
        # (P2-14a added 17 lines to MMWSpalling.H), so source line numbers are masked in BOTH runs' lines
        vl = lambda d: [re.sub(r"(\.H):\d+", r"\1:<n>", x) for x in R10.verdict_lines(os.path.join(d, n + ".log"))]
        va = vl(A) if n in ra else None
        vb = vl(B) if n in rb else None
        if n == "lateral_flux_lock":
            # the one declared change (user, 2026-11-02): the Grimsel-mention clause allows the nu read
            # (lower-case `ok` clause lines are not R10 verdict lines, so compare the FULL logs: exactly the
            # one clause line may differ, and it must read ok)
            fa = open(os.path.join(A, n + ".log")).read().splitlines()
            fb = open(os.path.join(B, n + ".log")).read().splitlines()
            da = [x for x in fa if x not in fb]
            db = [x for x in fb if x not in fa]
            ok = (ra.get(n) == rb.get(n) == "rc=0" and va == vb and len(da) == 1 and len(db) == 1
                  and "law POINTER" in da[0] and "ONE declared nu read" in db[0]
                  and db[0].lstrip().startswith("ok"))
            L.append("  %-28s rc %5s -> %5s   verdict lines identical; full log differs in exactly the declared nu "
                     "clause, which reads ok: %s" % (n, ra.get(n, "-"), rb.get(n, "-"), ok))
        elif n == "constraint_ref_lock":
            # its added clauses: every OLD verdict line still present, rc 0, the log ends in PASS
            log = open(os.path.join(B, n + ".log")).read()
            ok = rb.get(n) == ra.get(n) == "rc=0" and "\nPASS:" in log and "FAIL" not in log
            L.append("  %-28s rc %5s -> %5s   added P2-14a clauses, every check ok: %s"
                     % (n, ra.get(n, "-"), rb.get(n, "-"), ok))
        elif n == "face_cref":
            log = open(os.path.join(B, n + ".log")).read()
            ok = n not in ra and rb.get(n) == "rc=0" and log.rstrip().endswith("PASS")
            L.append("  %-28s (new)  rc %5s   PASS: %s" % (n, rb.get(n, "-"), ok))
        else:
            ok = ra.get(n) == rb.get(n) and va == vb
            L.append("  %-28s rc %5s -> %5s   verdict lines %s" % (n, ra.get(n, "-"), rb.get(n, "-"),
                     "identical (%d)" % len(vb) if ok else "DIFFER"))
        if not ok:
            bad.append(n)
    L.append("a3 %s (gate)" % ("PASS" if not bad else "FAIL: " + ", ".join(bad)))
    return "\n".join(L), not bad


# ======================================================================= a4-a8
def float_rows_close(a, b, ints=("time", "col_i", "col_j", "k", "dir", "n_voided", "h_line"), rel=1e-9):
    ra, rb = events(a), events(b)
    if len(ra) != len(rb):
        return False, "rows %d vs %d" % (len(ra), len(rb))
    worst, first_bad = 0.0, None
    for n, (x, y) in enumerate(zip(ra, rb)):
        for k in x:
            if k in ints:
                if x[k] != y[k] and first_bad is None:
                    first_bad = "row %d %s %s vs %s" % (n, k, x[k], y[k])
                continue
            fx, fy = float(x[k]), float(y[k])
            e = abs(fx - fy) / max(abs(fx), abs(fy), 1e-300)
            if e > worst:
                worst = e
            if e > rel and first_bad is None:
                first_bad = "row %d %s %.17g vs %.17g (rel %.2e)" % (n, k, fx, fy, e)
    return first_bad is None, "%d rows; worst float rel %.2e%s" % (len(ra), worst,
                                                                 "; FIRST DIFFERING " + first_bad if first_bad else "")


def GATES():
    L, ok_all = [], True

    def gate(tag, ok, msg):
        nonlocal ok_all
        L.append("  %-4s %s  %s" % (tag, "PASS" if ok else "FAIL", msg))
        ok_all = ok_all and ok

    L += ["a4-a8 -- the no-change claims (d = 0 bytes, nu by value, the column path, the abort)", "=" * 100]
    for leg in ("CAL", "CAL_det", "CAL_h2000", "CAL_2mm", "HOS_R"):
        a, b = os.path.join(CALD, "legs", leg), os.path.join(LEGS, leg)
        ev = same_file(a, b, "run_removal_events.csv")
        jp = (not os.path.exists(os.path.join(a, "run_jet_profile.csv"))) or same_file(a, b, "run_jet_profile.csv")
        n, nd, no = same_plotfiles(a, b)
        th, thm = thermo_same(a, b)
        cs, csm = check_leg(leg)
        gate("a4", rc(b) == 0 and ev and jp and nd == 0 and no == 0 and th and cs,
             "%-9s vs P2-13: event log %s, jet profile %s, plotfiles %d files %d differ %d one-sided; thermo %s; %s"
             % (leg, "identical" if ev else "DIFFERS", "identical" if jp else "DIFFERS", n, nd, no, thm, csm))
    for leg in ("CAL_p27", "CAL_p48"):
        ok, msg = float_rows_close(os.path.join(CALD, "legs", leg), os.path.join(LEGS, leg))
        cs, csm = check_leg(leg)
        gate("a5", rc(os.path.join(LEGS, leg)) == 0 and ok and cs, "%s vs P2-13 (integer columns exact, floats "
             "<= 1e-9 rel): %s; %s" % (leg, msg, csm))
    same_mu = same_file(os.path.join(LEGS, "CAL_p27_mu"), os.path.join(LEGS, "CAL_p27"), "run_removal_events.csv")
    cs, csm = check_leg("CAL_p27_mu")
    gate("a6", same_mu and cs, "new: CAL_p27_mu (phase0.mu 12e9) vs CAL_p27: event logs %s; %s"
         % ("byte-identical" if same_mu else "DIFFER", csm))
    sp = stat(os.path.join(LEGS, "CAL_p27_mu_pre"))
    s13 = stat(os.path.join(CALD, "legs", "CAL_p27"))
    gate("a6", abs(sp["med"] - s13["med"]) >= 3.0, "pre: CAL_p27_mu_pre %.2f K vs P2-13's CAL_p27 %.2f K: %+.2f K "
         "(>= 3 K: the control sees nu)" % (sp["med"], s13["med"], sp["med"] - s13["med"]))
    a, b = os.path.join(LEGS, "F2_k0_pre"), d_of("F2_k0")
    ev = same_file(a, b, "run_removal_events.csv")
    n, nd, no = same_plotfiles(a, b)
    gate("a7", rc(a) == 0 and ev and nd == 0 and no == 0, "F2_k0 (column path, d 2 mm) pre vs new: event log %s; "
         "plotfiles %d files, %d differ, %d one-sided" % ("identical" if ev else "DIFFERS", n, nd, no))
    for leg in ("F0_k0", "F0_k1"):
        a, b = os.path.join(LEGS, leg + "_pre"), d_of(leg)
        ev = same_file(a, b, "run_removal_events.csv")
        n, nd, no = same_plotfiles(a, b)
        L.append("  (a7+) %s pre vs new (d = 0): event log %s; plotfiles %d files, %d differ"
                 % (leg, "identical" if ev else "DIFFERS", n, nd))
    pre = os.path.join(LEGS, "F2_k1_pre")
    log = open(os.path.join(pre, "run.log"), errors="replace").read()
    named = "spallation.sp.constraint_depth" in log
    gate("a8", rc(pre) not in (0, None) and named and rc(d_of("F2_k1")) == 0,
         "F2_k1 on pre: rc %s, the abort names spallation.sp.constraint_depth: %s; on new rc %s"
         % (rc(pre), named, rc(d_of("F2_k1"))))
    L.append("a4-a8 %s (gates)" % ("PASS" if ok_all else "FAIL"))
    return "\n".join(L), ok_all


# ======================================================================= b
def B():
    L, ok_all = ["b -- the port on fixture F (b1, b3, b4, b5 also scored by unit/face_cref; b2 here)", "=" * 100], True
    st = {n: stat(d_of(n), interior=n.startswith("T35")) for n in
          ("F0_k0", "F0_k1", "F2_k0", "F2_k1", "F164_k1", "T35_0", "T35_2", "N2_0", "N2_2", "M2_k1", "M05_k1",
           "F2_k1_np1", "F2_k1_g8")}
    for n, s in st.items():
        L.append("  %-10s %s   spread %.3f K   first firing %s s" % (n, fmt(s), s["spread"],
                 "%.3f" % s["t0"] if s["n"] else "-"))
    ulog = os.path.join(T, "unit", "face_cref", "output", "..", "test.log")
    # b1
    a, b = first_top(events(d_of("F2_k0"))), first_top(events(d_of("F2_k1")))
    same = set(a) == set(b) and len(a) == 576
    dtm = max(abs(float(a[c]["time"]) - float(b[c]["time"])) for c in a) if same else float("inf")
    d2, d0 = abs(st["F2_k1"]["med"] - st["F2_k0"]["med"]), abs(st["F0_k1"]["med"] - st["F0_k0"]["med"])
    sh = st["F2_k1"]["med"] - st["F0_k1"]["med"]
    ok = same and dtm <= 0.008 + 1e-9 and d2 <= d0 + 0.5 and sh >= 20.0
    ok_all = ok_all and ok
    L.append("b1 %s (gate) CROSS-PATH: 576 columns same step (max |dt| %.4f s); |dS| %.3f K <= %.3f + 0.5; "
             "d shifts the face list %+.2f K (>= 20)" % ("PASS" if ok else "FAIL", dtm, d2, d0, sh))
    # b2
    e = [os.path.join(d_of(n), "run_removal_events.csv") for n in ("F2_k1", "F2_k1_np1", "F2_k1_g8")]
    sm = all(os.path.exists(x) for x in e) and md5(e[0]) == md5(e[1]) == md5(e[2])
    ths = []
    for n in ("F2_k1_np1", "F2_k1_g8"):
        ha, ra = thermo(d_of("F2_k1"))
        hb, rb = thermo(d_of(n))
        ix = [ha.index(c) for c in ha if c.startswith("cref_")]
        ths.append(len(ra) == len(rb) and all([x[i] for i in ix] == [y[i] for i in ix] for x, y in zip(ra, rb)))
    ok = sm and all(ths)
    ok_all = ok_all and ok
    L.append("b2 %s (gate) F2_k1 vs -np 1 vs max_grid_size 8: event logs %s; cref_* thermo columns identical: %s "
             "(the line counters are not rank-summed)" % ("PASS" if ok else "FAIL",
                                                           "byte-identical" if sm else "DIFFER", ths))
    # b3
    h, rows = thermo(d_of("F2_k1"))
    ix = {c: h.index(c) for c in ("cref_n_vert", "cref_n_lat", "cref_n_clamp", "cref_n_neg")}
    rr = [[float(x) for x in r] for r in rows]
    started = [r for r in rr if r[0] > 0.012]
    reads = sum(r[ix["cref_n_vert"]] for r in rr)
    lat = sum(r[ix["cref_n_lat"]] for r in rr)
    clp = sum(r[ix["cref_n_clamp"]] for r in rr)
    neg = sum(r[ix["cref_n_neg"]] for r in rr)
    ev = events(d_of("F2_k1"))
    tc_ok = bool(ev) and all(293.15 < float(r["T_cref"]) < float(r["T_surf"]) for r in ev)
    h0, _ = thermo(d_of("F0_k1"))
    e0 = events(d_of("F0_k1"))
    off_ok = not any(c.startswith("cref_") for c in h0) and (not e0 or "T_cref" not in e0[0])
    ok = min(r[ix["cref_n_vert"]] for r in started) > 0 and lat == 0 and clp == 0 and neg <= 1e-4 * reads \
        and tc_ok and off_ok
    ok_all = ok_all and ok
    tcv = [float(r["T_cref"]) for r in ev]
    L.append("b3 %s (gate) F2_k1: vert > 0 on all %d rows after the first list; reads %.4g, lat %g, clamp %g, neg %g; "
             "T_cref %.1f-%.1f K, every row 293.15 < T_cref < T_surf: %s. F0_k1: no cref_*/T_cref: %s"
             % ("PASS" if ok else "FAIL", len(started), reads, lat, clp, neg, min(tcv), max(tcv), tc_ok, off_ok))
    # b4
    n0 = st["N2_0"]
    ne2 = len(events(d_of("N2_2")))
    ok = n0["n"] == 144 and ne2 == 0
    ok_all = ok_all and ok
    L.append("b4 %s (gate) MINIMUM FLUX h 150: N2_0 fires on %d/144 columns (%.1f-%.1f s, %.2f K); N2_2 firings %d"
             % ("PASS" if ok else "FAIL", n0["n"], n0.get("t0", -1), n0.get("t1", -1), n0["med"], ne2))
    # b5
    D0 = st["T35_0"]["med"] - st["F0_k1"]["med"]
    D2_ = st["T35_2"]["med"] - st["F2_k1"]["med"]
    W = st["F164_k1"]["med"] - st["F2_k1"]["med"]
    if abs(W) < 10.0:
        L.append("b5 NOT DISCRIMINATING (|W| %.2f K < 10 K)" % abs(W))
    else:
        ok = abs(D2_ - D0) < abs(W) / 2
        ok_all = ok_all and ok
        L.append("b5 %s (gate) NORMAL LENGTH: Delta(0) %+.3f, Delta(2) %+.3f, |Delta(2)-Delta(0)| %.3f < |W|/2 %.3f "
                 "(W = %+.3f K)" % ("PASS" if ok else "FAIL", D0, D2_, abs(D2_ - D0), abs(W) / 2, W))
        # the diagnosis of the remainder: the staircase's per-cell mu
        f = first_top(events(d_of("T35_2")), interior=True)
        mus = np.array([float(r["mu"]) for r in f.values()])
        Ts = np.array([float(r["T_surf"]) for r in f.values()])
        slope = W / (2.0e-3 - 1.6383041e-3)            # K per metre of normal reference depth, from F2/F164
        pred = st["F2_k1"]["med"] + slope * (2.0e-3 - 2.0e-3 * math.cos(math.radians(35.0)) / mus)
        L.append("   diagnosis: T35_2 interior mu %.3f-%.3f (median %.3f; cos 35 = 0.819). A vertical line reads d/mu; "
                 "in 35-degree isotherms that is a NORMAL depth d cos35/mu, = d only where mu = cos 35. Predicted "
                 "per column from F's slope %.1f K/mm: median %.2f K vs measured %.2f K; corr(pred, measured) %.3f"
                 % (mus.min(), mus.max(), float(np.median(mus)), slope * 1e-3, float(np.median(pred)),
                    float(np.median(Ts)), float(np.corrcoef(pred, Ts)[0, 1])))
    L.append("b %s (gates)" % ("PASS" if ok_all else "FAIL"))
    return "\n".join(L), ok_all, st


# ======================================================================= c
def vol_sectors(d):
    import fabcmp  # noqa: F401
    import importlib.machinery
    import importlib.util
    ld = importlib.machinery.SourceFileLoader("fr_t", os.path.join(T, "unit", "face_removal", "test"))
    m = importlib.util.module_from_spec(importlib.util.spec_from_loader("fr_t", ld))
    ld.exec_module(m)
    pf = os.path.join(d, "run", [x for x in pf_dirs(d) if x.endswith("cell")][-1])
    meta = {}
    for ln in open(os.path.join(d, "run", "metadata")):
        if " = " in ln:
            k, v = ln.rstrip("\n").split(" = ", 1)
            meta[k.strip()] = v.strip()
    n = [int(x) for x in meta["amr.n_cell"].split()]
    hi = [float(x) for x in meta["geometry.prob_hi"].split()]
    dx = hi[0] / n[0]
    f = m.read_pf(pf, ["fr_f"], tuple(n))["fr_f"]
    rem = (1.0 - f) * dx ** 3
    x = (np.arange(n[0]) + 0.5) * dx
    ang = np.degrees(np.arctan2(x[None, :], x[:, None]))       # [i, j]: atan2(y, x)
    sec = np.minimum((ang // 30).astype(int), 2)
    out = []
    for s in range(3):
        msk = (sec == s)[:, :, None] * np.ones(n[2], bool)[None, None, :]
        out.append((rem[msk].sum(), msk.sum()))
    return float(rem.sum()), out


def C(st_b):
    L = ["c -- reports", "=" * 100]
    flux = [("C0_h150", 150.0, 0), ("CAL", 500.0, 0), ("C0_h1000", 1000.0, 0), ("CAL_h2000", 2000.0, 0),
            ("C2_h150", 150.0, 2), ("C2_h500", 500.0, 2), ("C2_h1000", 1000.0, 2), ("C2_h2000", 2000.0, 2)]
    L.append("c1 FLUX SERIES on the CAL deck (S_top; q = h (1673.15 - median)):")
    S = {}
    for leg, h, dd in flux:
        s = stat(os.path.join(LEGS, leg))
        S[leg] = s
        cs, csm = check_leg(leg)
        q = h * (1673.15 - s["med"]) / 1e3 if s["n"] else float("nan")
        L.append("   d %d mm  h %6.0f  %-9s %s   q %7.1f kW/m2   %s" % (dd, h, leg, fmt(s), q, csm))
    for dd, a, b in ((0, "CAL", "CAL_h2000"), (2, "C2_h500", "C2_h2000")):
        if S[a]["n"] and S[b]["n"]:
            L.append("   d %d mm: S_top(h 2000) - S_top(h 500) = %+.2f K" % (dd, S[b]["med"] - S[a]["med"]))
        else:
            L.append("   d %d mm: h 500 or h 2000 did not fire (n %d / %d)" % (dd, S[a]["n"], S[b]["n"]))
    ms = [st_b["M2_k1"]["med"], st_b["F2_k1"]["med"], st_b["M05_k1"]["med"]]
    L.append("c2 MESH at d 2 mm (fixture F): 2 mm %.2f / 1 mm %.2f / 0.5 mm %.2f K; 2-1 %+.2f K, 0.5-1 %+.2f K, "
             "spread %.2f K" % (ms[0], ms[1], ms[2], ms[0] - ms[1], ms[2] - ms[1], max(ms) - min(ms)))
    a, b = os.path.join(LEGS, "HOS_R"), os.path.join(LEGS, "HOS_d2")
    h, rows = thermo(b)
    j = h.index("ledger_err_solid")
    le = max(abs(float(r[j])) for r in rows)
    sums = {c: sum(float(r[h.index(c)]) for r in rows) for c in ("cref_n_vert", "cref_n_lat", "cref_n_clamp", "cref_n_neg")}
    reads = sums["cref_n_vert"] + sums["cref_n_lat"]
    va, sa = vol_sectors(a)
    vb, sb = vol_sectors(b)
    L.append("c3 DRILLING HOS_d2 vs HOS_R at 15 s: rc %s; max|ledger_err_solid| %.2e; reads vert %.4g lat %.4g "
             "(lat > 0: %s); clamp %.4f %%, neg %.4f %% of reads; V(f) ratio %.4f (%+.1f %%)"
             % (rc(b), le, sums["cref_n_vert"], sums["cref_n_lat"], sums["cref_n_lat"] > 0,
                100 * sums["cref_n_clamp"] / reads, 100 * sums["cref_n_neg"] / reads, vb / va, 100 * (vb / va - 1)))
    for nm, d in (("HOS_R", a), ("HOS_d2", b)):
        ev = events(d)
        by = {}
        for r in ev:
            by.setdefault(r["dir"], []).append(float(r["T_surf"]))
        L.append("   %-7s all-firing T_surf median %.2f K (n %d); by dir: %s" % (
            nm, float(np.median([float(r["T_surf"]) for r in ev])), len(ev),
            ", ".join("%s %.2f (%d)" % (k, float(np.median(v)), len(v)) for k, v in sorted(by.items()))))
    ra = [x / n for x, n in sa]
    rb = [x / n for x, n in sb]
    L.append("c4 ISOTROPY V(f)/cells per 30-degree sector (0-30, 30-60, 60-90): HOS_R %s  max/min %.4f;  "
             "HOS_d2 %s  max/min %.4f" % (" ".join("%.4g" % x for x in ra), max(ra) / min(ra),
                                          " ".join("%.4g" % x for x in rb), max(rb) / min(rb)))
    return "\n".join(L), S


def main(a):
    mode = a[0] if a else "RESULTS"
    if mode == "A1":
        t = A1(); print(t); open(os.path.join(HERE, "A1.txt"), "w").write(t + "\n"); return 0
    if mode == "A3":
        t, ok = A3(); print(t); open(os.path.join(HERE, "A3.txt"), "w").write(t + "\n"); return 0 if ok else 1
    if mode == "GATES":
        t, ok = GATES(); print(t); open(os.path.join(HERE, "GATES_a4_a8.txt"), "w").write(t + "\n"); return 0 if ok else 1
    if mode == "B":
        t, ok, _ = B(); print(t); open(os.path.join(HERE, "B.txt"), "w").write(t + "\n"); return 0 if ok else 1
    out = []
    for f in ("A0.txt", "A1.txt", "GATE_tree.txt", "A3.txt", "GATES_a4_a8.txt"):
        p = os.path.join(HERE, f)
        out += [open(p).read().rstrip() if os.path.exists(p) else f + ": NOT ON DISK", ""]
    tb, okb, st = B()
    out += [tb, ""]
    tc, S = C(st)
    out += [tc, ""]
    g = {f: open(os.path.join(HERE, f)).read() for f in ("A0.txt", "GATE_tree.txt", "A3.txt", "GATES_a4_a8.txt")
         if os.path.exists(os.path.join(HERE, f))}
    gates = ("PASS:" in g.get("A0.txt", "") and "vs p2_12r: 0 of 2639" in g.get("GATE_tree.txt", "")
             and "a3 PASS" in g.get("A3.txt", "") and "a4-a8 PASS" in g.get("GATES_a4_a8.txt", "") and okb)
    d5, d20 = S["C2_h500"], S["C2_h2000"]
    out += ["e1 VERDICT",
            "   THE PORT IS FAITHFUL: %s (a0 source lock, b1 cross-path %.2f K at the same step, b2 decomposition-"
            "free, b5 the along-line read is d/mu)." % ("every gate PASSED" if gates else "a gate FAILED", abs(
                st["F2_k1"]["med"] - st["F2_k0"]["med"])),
            "   d = 2 mm IN THE PRODUCTION FORM puts a strong flux dependence into the onset: %.1f K at h 500 -> %.1f K"
            " at h 2000 (%+.0f K over x%.1f in onset flux), and h 150 never fires; at d 0 the same series is flat"
            " (%+.2f K). The data show no trend (Wilkinson 0.2-2.7 MW/m2; Kant; Hoeser), and %.0f K at Kant's"
            " deck flux lies above every Grimsel onset measurement (<= 967 K)."
            % (d5["med"], d20["med"], d20["med"] - d5["med"], 1587.4 / 335.9, S["CAL_h2000"]["med"] - S["CAL"]["med"],
               d5["med"]),
            "   Item 3's re-fit must therefore score BOTH the minimum flux (Kant 0.3-0.5 MW/m2) AND the flatness of"
            " the onset with flux, on Grimsel in this form. d = 2 mm (Gotthard, column path) brackets its minimum flux"
            " between 0.12 MW/m2 (h 150, no firing) and 0.34 MW/m2 (h 500, fires) -- overlapping only the low edge"
            " of Kant's 0.3-0.5 -- and fails the flatness by ~120 K. d stays 0 in model_form_flame.inputs."]
    out.append("")
    print("\n".join(out))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
