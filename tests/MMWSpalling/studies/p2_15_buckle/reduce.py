#!/usr/bin/env python3
"""P2-15 -- every row. Usage: reduce.py [A1|A3|A4|A6|B] (no mode: RESULTS.txt to stdout).

S_top (SR 19): per column the first event-log row with dir = 5 and ts_fb = 0, its T_surf. Fired fraction: columns
with >= 1 top-face row / 576. Blocked fraction: sum(bk_n_block) / sum(bk_n_eval) over thermo rows (per-step
columns, SR 7(i); a RATIO of two sums over the same rows, so plot_int cancels).
"""
import csv
import hashlib
import importlib.machinery
import importlib.util
import math
import os
import re
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
T = os.path.join(ROOT, "tests", "MMWSpalling")
LEGS = os.path.join(HERE, "legs")
P14 = os.path.join(T, "studies", "p2_14a_cref_face")
CALD = os.path.join(T, "validation", "kant", "grimsel_calibration")
FB = os.path.join(T, "unit", "face_buckle")
CHECK = os.path.join(T, "validation", "shared", "check_shared.py")
PY = sys.executable
CF = 0.4


def load(name, path):
    ld = importlib.machinery.SourceFileLoader(name, path)
    m = importlib.util.module_from_spec(importlib.util.spec_from_loader(name, ld))
    ld.exec_module(m)
    return m


DS = load("ds15", os.path.join(HERE, "derive_span.py"))
R14 = load("r14", os.path.join(P14, "reduce.py"))


def events(d):
    p = os.path.join(d, "run_removal_events.csv")
    return list(csv.DictReader(open(p))) if os.path.exists(p) else []


def first_top(rows):
    f = {}
    for r in rows:
        if r["dir"] == "5" and float(r["ts_fb"]) == 0.0:
            f.setdefault((r["col_i"], r["col_j"]), r)
    return f


def stat(rows, key="T_surf"):
    v = np.array([float(r[key]) for r in first_top(rows).values()])
    if not len(v):
        return dict(n=0, med=float("nan"), p10=float("nan"), p90=float("nan"))
    return dict(n=len(v), med=float(np.median(v)), p10=float(np.percentile(v, 10)), p90=float(np.percentile(v, 90)))


def fmt(s, u="K", k=1.0):
    return ("n %4d  median %8.2f %s  p10 %7.2f  p90 %7.2f" % (s["n"], k * s["med"], u, k * s["p10"], k * s["p90"])) \
        if s["n"] else "n    0  (no top-face firing)"


def thermo(d):
    L = open(os.path.join(d, "run", "thermo.dat")).read().splitlines()
    h = [c.strip() for c in L[0].split("\t")]
    return h, [[float(x) for x in r.split("\t")] for r in L[1:] if r.strip()]


def tsum(d, c):
    h, rows = thermo(d)
    return sum(r[h.index(c)] for r in rows) if c in h else float("nan")


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def check_leg(leg):
    d = os.path.join(LEGS, leg)
    deck = open(os.path.join(d, "DECK.txt")).read().strip()
    allow = open(os.path.join(d, "ALLOW.txt")).read().split()
    a = [PY, CHECK, deck, os.path.join(d, "argv.sh")] + (["--allow"] + allow if allow else [])
    p = subprocess.run(a, cwd=ROOT, capture_output=True, text=True)
    return p.returncode == 0, p.stdout.strip().splitlines()[-1]


def rc(d):
    p = os.path.join(d, "rc")
    return int(open(p).read()) if os.path.exists(p) else None


# ---------------------------------------------------------------- a1
def A1():
    L = ["a1 -- ell_eff at first top firing on P2-14a's span-off CAL legs (zero compute) and the implied bk", "=" * 100]
    for leg, h in (("C0_h150", 150), ("CAL", 500), ("C0_h1000", 1000), ("CAL_h2000", 2000)):
        f = first_top(events(os.path.join(P14, "legs", leg)))
        ell = np.array([float(r["ell_eff"]) for r in f.values()])
        Ts = np.array([float(r["T_surf"]) for r in f.values()])
        cl = DS.C_L(float(np.median(Ts)))
        x = CF * float(np.median(ell))
        L.append("  h %4d %-9s ell_eff %.3f mm (p10 %.3f, p90 %.3f)  T_s %.1f K  C_L %.2f  C_L c_f ell %.2f cm  "
                 "bk @12 mm %.3f, @15 mm %.3f" % (h, leg, 1e3 * np.median(ell), 1e3 * np.percentile(ell, 10),
                                                  1e3 * np.percentile(ell, 90), np.median(Ts), cl, 100 * cl * x,
                                                  cl * x / 0.012, cl * x / 0.015))
    L.append("  packet's offline estimates: C_L c_f ell = 5.1 / 1.5 / 0.76 / 0.38 cm at h 150 / 500 / 1000 / 2000")
    L.append("a1 (report)")
    return "\n".join(L)


# ---------------------------------------------------------------- a3
def mask(x):
    """source line numbers out of abort text (P2-14a DR-6), incl. the fixed-width tail a test may print, where only
    the last digits of `<file>.H:<line> (<function>)` survive"""
    x = re.sub(r"(\.H):\d+", r"\1:<n>", x)
    return re.sub(r"\d+ \((operator\(\)|[A-Za-z_]+)\)", r"<n> (\1)", x)


def A3():
    R10 = load("r10q15", os.path.join(T, "studies", "p2_10_retire", "run_p2_10.py"))
    A, B = os.path.join(P14, "UNIT_p2_14a"), os.path.join(HERE, "UNIT_p2_15")
    rcm = lambda d: {l.split()[0]: l.split()[1] for l in open(os.path.join(d, "rc.txt")).read().splitlines()[1:]}
    ra, rb = rcm(A), rcm(B)
    L = ["a3 -- every kept unit on the P2-15 binary vs P2-14a's UNIT_p2_14a (the last accepted state; source line "
         "numbers masked, P2-14a DR-6): rc and verdict lines; the two locks also by full log", "=" * 100]
    bad = []
    for n in sorted(set(ra) | set(rb)):
        vl = lambda d: [mask(x) for x in R10.verdict_lines(os.path.join(d, n + ".log"))]
        if n in ("buckle_lock", "face_buckle"):
            log = open(os.path.join(B, n + ".log")).read()
            ok = n not in ra and rb.get(n) == "rc=0" and re.search(r"^PASS", log, re.M) is not None
            L.append("  %-28s (new)  rc %5s   PASS: %s" % (n, rb.get(n, "-"), ok))
        else:
            va = vl(A) if n in ra else None
            vb = vl(B) if n in rb else None
            ok = ra.get(n) == rb.get(n) and va == vb
            extra = ""
            if n in ("constraint_ref_lock", "lateral_flux_lock"):
                fa = [mask(x) for x in open(os.path.join(A, n + ".log")).read().splitlines()]
                fb = [mask(x) for x in open(os.path.join(B, n + ".log")).read().splitlines()]
                ok = ok and fa == fb
                extra = "; full log identical: %s" % (fa == fb)
            L.append("  %-28s rc %5s -> %5s   verdict lines %s%s" % (n, ra.get(n, "-"), rb.get(n, "-"),
                     "identical (%d)" % len(vb) if va == vb else "DIFFER", extra))
        if not ok:
            bad.append(n)
    L.append("a3 %s (gate)" % ("PASS" if not bad else "FAIL: " + ", ".join(bad)))
    return "\n".join(L), not bad


# ---------------------------------------------------------------- a4
def A4():
    L, ok = ["a4 -- span = 0 is inert: the P2-15 binary vs the frozen runs", "=" * 100], True
    for leg in ("CAL", "HOS_R"):
        a, b = os.path.join(CALD, "legs", leg), os.path.join(LEGS, leg)
        ev = R14.same_file(a, b, "run_removal_events.csv")
        n, nd, no = R14.same_plotfiles(a, b)
        th, thm = R14.thermo_same(a, b)
        cs, csm = check_leg(leg)
        o = rc(b) == 0 and ev and nd == 0 and no == 0 and th and cs
        ok = ok and o
        L.append("  %s  %-6s vs P2-13: event log %s; plotfiles %d files, %d differ, %d one-sided; thermo %s; %s"
                 % ("PASS" if o else "FAIL", leg, "identical" if ev else "DIFFERS", n, nd, no, thm, csm))
    a, b = os.path.join(HERE, "ref_p2_14a", "F0_k1"), os.path.join(T, "unit", "face_cref", "output", "F0_k1")
    ev = R14.same_file(a, b, "run_removal_events.csv")
    n, nd, no = R14.same_plotfiles(a, b)
    o = ev and nd == 0 and no == 0
    ok = ok and o
    L.append("  %s  F0_k1 (unit/face_cref, re-run by a3 on this binary) vs P2-14a's (archived ref_p2_14a/): event log "
             "%s; plotfiles %d files, %d differ, %d one-sided" % ("PASS" if o else "FAIL",
                                                                   "identical" if ev else "DIFFERS", n, nd, no))
    L.append("a4 %s (gate)" % ("PASS" if ok else "FAIL"))
    return "\n".join(L), ok


# ---------------------------------------------------------------- a6
def A6():
    sp = float(open(os.path.join(HERE, "SPAN.txt")).read().split()[0])
    lo, hi = sp / 4e-3, sp / 1e-3
    return ("a6 -- L_span %.4g m = %.1f-%.1f grains of 1-4 mm (Kant: feldspar up to ~4 mm); inside the lab torch "
            "spalls' 10-20 grains (Rauenzahn & Tester): %s; field spalls are wider (several cm, Walsh & Lomov) (report)"
            % (sp, lo, hi, "overlaps" if (lo <= 20 and hi >= 10) else "NO"))


# ---------------------------------------------------------------- b
def B():
    lp = os.path.join(HERE, "UNIT_p2_15", "face_buckle.log")      # written by a3's run of the unit
    log = open(lp).read() if os.path.exists(lp) else ""
    L = ["b1-b4 -- unit/face_buckle (its a3 run's log, UNIT_p2_15/face_buckle.log) plus the cross-binary check "
         "against P2-14a", "=" * 100]
    L += ["  " + x.strip() for x in log.splitlines() if x.strip().startswith(("PASS", "FAIL", "REPORT"))]
    a = events(os.path.join(HERE, "ref_p2_14a", "F0_k1"))
    b = events(os.path.join(FB, "output", "F_S"))
    sh = [k for k in (a[0].keys() if a else []) if k in (b[0] if b else {})]
    same = len(a) == len(b) > 0 and all(x[k] == y[k] for x, y in zip(a, b) for k in sh)
    L.append("  %s  F_S vs P2-14a's F0_k1 (other binary): every shared event-log column of every row identical "
             "(%d rows, %d columns)" % ("PASS" if same else "FAIL", len(b), len(sh)))
    ok = same and log.rstrip().endswith("PASS")
    L.append("b %s (gates)" % ("PASS" if ok else "FAIL"))
    return "\n".join(L), ok


# ---------------------------------------------------------------- c
def leg_row(leg):
    d = os.path.join(LEGS, leg)
    rows = events(d)
    st = stat(rows)
    ev, bl = tsum(d, "bk_n_eval"), tsum(d, "bk_n_block")
    bk = [float(r["bk"]) for r in rows if "bk" in r]
    cs, csm = check_leg(leg)
    return st, ev, bl, bk, cs, csm


def C():
    L = ["c -- reports (span %s m from SPAN.txt)" % open(os.path.join(HERE, "SPAN.txt")).read().split()[0], "=" * 100]
    L.append("c1 GOTTHARD (by construction: the span is derived on 130/140) and c1b BETHEL (the independent check):")
    for leg in ("KG50_S", "KG80_S", "KG110_S", "KG130_S", "KG140_S", "KG150_S", "KG160_S", "KB150_S", "KB160_S"):
        st, ev, bl, bk, cs, csm = leg_row(leg)
        L.append("   %-8s fired %3d/576 = %5.1f %%   %s   blocked %s of %d evaluated   bk median %s   %s" % (
            leg, st["n"], 100.0 * st["n"] / 576, fmt(st), "%.4f" % (bl / ev) if ev > 0 else "-", ev,
            "%.3f" % float(np.median(bk)) if bk else "-", csm))
    sp = float(open(os.path.join(HERE, "SPAN.txt")).read().split()[0])
    # ⚠ THE FLAT DECK IS ALL-OR-NOTHING. A blocked firing leaves the field untouched, so until the deck's EARLIEST
    # firing passes, the floor stays flat and every skin keeps thickening (bk only grows): the earliest firing of the
    # span-off twin decides the span-on deck. Its bk (flat floor, mu = 1) is printed for every pair; the medians of
    # later first firings come from a history the blocked deck never has (neighbours' staircase, mu 0.82-0.86).
    L.append("   THE DECK'S DISCRIMINANT -- the span-off twin's EARLIEST firing, its bk at the span (mu = 1 on the still-flat "
             "floor) vs the span-on outcome:")
    need = {}
    for off, on, want in (("KG130_0", "KG130_S", "fire"), ("KG140_0", "KG140_S", "silent"), ("KG150_0", "KG150_S", "silent"),
                          ("KG160_0", "KG160_S", "silent"), ("KB150_0", "KB150_S", "fire"), ("KB160_0", "KB160_S", "silent")):
        r = events(os.path.join(LEGS, off))[0]
        bk0 = DS.C_L(float(r["T_surf"])) * CF * float(r["ell_eff"]) * float(r["mu"]) / sp
        n_on = len(first_top(events(os.path.join(LEGS, on))))
        need[off] = (bk0 * sp, want)
        L.append("     %-8s earliest firing t %7.3f s  ell %.3f mm  bk %.3f  -> %-8s %3d/576   (Kant: %s)" % (
            off, float(r["time"]), 1e3 * float(r["ell_eff"]), bk0, on, n_on, want))
    lo = max(v for v, w in need.values() if w == "fire")
    hi = min(v for v, w in need.values() if w == "silent")
    lo_g = need["KG130_0"][0]
    hi_g = min(need[k][0] for k in ("KG140_0", "KG150_0", "KG160_0"))
    L.append("   spans that reproduce Gotthard alone: [%.2f, %.2f) mm; BOTH of Kant's rocks: %s" % (
        1e3 * lo_g, 1e3 * hi_g, "[%.3f, %.3f) mm (width %.3f mm)" % (1e3 * lo, 1e3 * hi, 1e3 * (hi - lo)) if lo < hi
        else "EMPTY (fire needs >= %.3f mm, silence < %.3f mm)" % (1e3 * lo, 1e3 * hi)))
    L.append("c2 FLUX SERIES on Grimsel (span on vs P2-14a's span-off twins):")
    for leg, twin, h in (("S_h150", "C0_h150", 150), ("S_h500", "CAL", 500), ("S_h1000", "C0_h1000", 1000),
                         ("S_h2000", "CAL_h2000", 2000)):
        st, ev, bl, bk, cs, csm = leg_row(leg)
        fa = first_top(events(os.path.join(P14, "legs", twin)))
        fb = first_top(events(os.path.join(LEGS, leg)))
        both = set(fa) & set(fb)
        same = sum(1 for c in both if fa[c]["T_surf"] == fb[c]["T_surf"] and fa[c]["time"] == fb[c]["time"])
        L.append("   h %4d %-8s fired %3d/576  %s  blocked %s of %d  bk median %s;  columns firing in both: %d, "
                 "identical first firing (time and T_surf) on %d   %s" % (
                     h, leg, st["n"], fmt(st), "%.4f" % (bl / ev) if ev > 0 else "-", ev,
                     "%.3f" % float(np.median(bk)) if bk else "-", len(both), same, csm))
    a, b = os.path.join(LEGS, "HOS_R"), os.path.join(LEGS, "HOS_S")
    hb, rows = thermo(b)
    le = max(abs(r[hb.index("ledger_err_solid")]) for r in rows)
    ev, bl = tsum(b, "bk_n_eval"), tsum(b, "bk_n_block")
    cls = {c: tsum(b, "bk_n_block_" + c) for c in ("floor", "wall", "roof")}
    va, sa = R14.vol_sectors(a)
    vb, sb = R14.vol_sectors(b)
    ea, eb = events(a), events(b)
    thr = [float(r["bk"]) for r in eb if r["through"] == "1"]
    ra = [x / n for x, n in sa]
    rb = [x / n for x, n in sb]
    sa_, sb_ = stat(ea), stat(eb)
    L.append("c3 DRILLING HOS_S vs HOS_R at 15 s (report only; uptake caveat): rc %s; max|ledger_err_solid| %.2e; "
             "V(f) ratio %.4f; blocked %.4f of %.4g evaluated (sampled rows, plot_int 6) -- by class floor %.4g, wall "
             "%.4g, roof %.4g; S_top (dir 5) %.2f vs %.2f K; sector max/min %.4f vs %.4f; through rows %d, all bk "
             "finite: %s" % (rc(b), le, vb / va, bl / ev if ev else float("nan"), ev, cls["floor"], cls["wall"],
                             cls["roof"], sb_["med"], sa_["med"], max(rb) / min(rb), max(ra) / min(ra), len(thr),
                             all(math.isfinite(x) for x in thr)))
    return "\n".join(L)


def main(a):
    mode = a[0] if a else "RESULTS"
    for m, fn, f in (("A1", A1, "A1.txt"), ("A6", A6, "A6.txt")):
        if mode == m:
            t = fn(); print(t); open(os.path.join(HERE, f), "w").write(t + "\n"); return 0
    for m, fn, f in (("A3", A3, "A3.txt"), ("A4", A4, "A4.txt"), ("B", B, "B.txt")):
        if mode == m:
            t, ok = fn(); print(t); open(os.path.join(HERE, f), "w").write(t + "\n"); return 0 if ok else 1
    out = []
    for f in ("A0.txt", "A1.txt", "GATE_tree.txt", "A3.txt", "A4.txt", "A5.txt", "A6.txt", "B.txt"):
        p = os.path.join(HERE, f)
        out += [open(p).read().rstrip() if os.path.exists(p) else f + ": NOT ON DISK", ""]
    out += [C(), ""]
    out += ["e1 VERDICT",
            "   ⚠ The flat deck is ALL-OR-NOTHING: its earliest firing's bk predicts all ten span-on outcomes (fire iff "
            "< 1). P1's partial fired fractions (60-95 %, 5-40 %) cannot occur on it.",
            "   B IS FAITHFUL: a0 (lock), b1-b4 (A unmoved where B passes, byte for byte; N blocked; witnesses; bk / ell "
            "mesh-free to 0.05 %), a2-a4 (span = 0 inert).",
            "   THE SPAN, derived on Gotthard 130/140 (by construction): 11.0 mm (~3-11 grains). Gotthard 50-130 fire, "
            "140-160 are silent.",
            "   BETHEL: KB150_S fires on 0 of 576 columns where Kant's Bethel spalled at 150 mm (its span-off control "
            "fires on all 576) -- P1b's REFUTED reading at the derived span. The deck's own discriminant (above) says "
            "why: under the conductivity-only transfer Bethel at 150 mm and Gotthard at 140 mm sit at the SAME bk "
            "(1.090 vs 1.092), so one span can satisfy both rocks only in a ~0.02 mm window at ~12.0 mm. The model does "
            "not separate the two rocks the way Kant's tests do; the likeliest cause is the strength side kept at "
            "Grimsel's (Bethel's onset is ~30 K lower) -- a finding for item 3, not a fit to make here.",
            "   FLUX SERIES on Grimsel: where the span passes, the onset is untouched (h 1000 / 2000: every column's "
            "first firing identical in time and T_surf); h 150 AND h 500 are fully blocked, so with this span the "
            "production rock's minimum flux lies between ~0.41 and ~0.82 MW/m2 -- above Kant's 0.3-0.5 for his rocks.",
            "   DRILLING (report only): HOS_S removes the same volume (ratio 1.0003); 0.12 % of evaluations blocked, "
            "all on walls; under the uptake-scaled heat the flake rarely approaches the span.",
            "   d CANNOT be retired on this evidence: the thickness condition reproduces Gotthard's bracket only by "
            "construction and fails the independent rock. spall.buckle_span stays 0 in production."]
    out.append("")
    print("\n".join(out))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
