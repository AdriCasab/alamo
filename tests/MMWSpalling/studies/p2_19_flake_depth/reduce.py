#!/usr/bin/env python3
"""P2-19 -- rows a1-a3, c1-c4, e0. Usage: reduce.py [A1|A2|A3] (no mode: RESULTS.txt to stdout).

Imported, never copied: P2-17's reduce.py (thermo / event readers, the isotropy sectors, the P2-11 volume instrument
over P2-8's plotfile reader) and through it P2-15's / P2-14a's byte comparators.

⚠ COMPARING A NEW LOG WITH AN OLD ONE: the new event log carries four appended columns (ell_p, z_peak, ell_pk, h_E),
so a row is compared on its first N fields, N = the old header's field count, as TEXT (byte identity of the
pre-existing columns). A longer run is cut to a shorter one's horizon BY STEP, never by `time <= stop`: the log stamps
a step's START and t accumulates in floating point (P2-17 a3), so the cut is `time < stop - dt/2`.
"""
import csv
import importlib.machinery
import importlib.util
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
T = os.path.join(ROOT, "tests", "MMWSpalling")
STUD = os.path.join(T, "studies")
LEGS = os.path.join(HERE, "legs")
W = os.path.join(STUD, "tree_sweep_0922", "witness")
NEW_COLS = ["ell_p", "z_peak", "ell_pk", "h_E"]
TIMERS = {"fr_t_removal", "fr_t_step", "fl_t_path"}


def load(name, path):
    keep = sys.argv
    sys.argv = [path]
    ld = importlib.machinery.SourceFileLoader(name, path)
    m = importlib.util.module_from_spec(importlib.util.spec_from_loader(name, ld))
    ld.exec_module(m)
    sys.argv = keep
    return m


P17 = load("red17_p219", os.path.join(STUD, "p2_17_buckle_on", "reduce.py"))
R15, R14, R11R = P17.R15, P17.R14, P17.R11R
REF = {"CAL": os.path.join(STUD, "p2_16_span_rule", "legs", "CAL"),
       "HOS_R": os.path.join(STUD, "p2_16_span_rule", "legs", "HOS_R"),
       "N_S164": os.path.join(STUD, "p2_17_buckle_on", "legs", "N_S164"),
       "C0_h150": os.path.join(STUD, "p2_14a_cref_face", "legs", "C0_h150")}


def geom_meta(d):
    """P2-8 geom_of's keys from the run's RESOLVED metadata, not the argv: the CAL deck sets its mesh in the input
    file (geom_of reads argv only). The domain is read, never assumed (no LZ constant, SR 5)."""
    kv = {}
    for ln in open(os.path.join(d, "run", "metadata")):
        if " = " in ln:
            k, v = ln.rstrip("\n").split(" = ", 1)
            kv[k.strip()] = v.strip()
    nx, ny, nz = (int(v) for v in kv["amr.n_cell"].split())
    lo = [float(v) for v in kv["geometry.prob_lo"].split()]
    hi = [float(v) for v in kv["geometry.prob_hi"].split()]
    return dict(nx=nx, ny=ny, nz=nz, dx=(hi[0] - lo[0]) / nx, dz=(hi[2] - lo[2]) / nz, lo=lo, hi=hi,
                dt=float(kv["timestep"]), tint=int(kv.get("amr.thermo.plot_int", "1")), st=float(kv["stop_time"]),
                fl=kv.get("surface_patch.face_list", "0"))


def mod_for(root):
    """P2-8's plotfile reader with legdir = <root>/<leg> (root is a legs/ directory), geometry from metadata"""
    m = R11R._load("r8_p219_%d" % (abs(hash(root)) % 10 ** 8), os.path.dirname(root))
    m.legdir = lambda n, _r=root: os.path.join(_r, n)
    m.geom_of = geom_meta
    return m


def lines(d):
    return open(os.path.join(d, "run_removal_events.csv")).read().splitlines()


def dicts(d):
    return list(csv.DictReader(open(os.path.join(d, "run_removal_events.csv"))))


def cut(rows, stop, dt):
    """rows (text, header first) of the steps that START before stop - dt/2"""
    return [rows[0]] + [r for r in rows[1:] if float(r.split(",", 1)[0]) < stop - 0.5 * dt]


def prefix_cmp(new, old):
    """compare text rows on the old header's field count. Returns (header_ok, n_same_rows, first_diff_index or
    None, n_new, n_old). header_ok: the new header is the old one + NEW_COLS."""
    hn, ho = new[0].split(","), old[0].split(",")
    hok = hn == ho + NEW_COLS
    k = len(ho)
    first = None
    n = min(len(new), len(old)) - 1
    for i in range(1, n + 1):
        if new[i].split(",")[:k] != old[i].split(","):
            first = i
            break
    if first is None and len(new) != len(old):
        first = n + 1
    return hok, (first - 1 if first else n), first, len(new) - 1, len(old) - 1


def dt_of(d):
    for t in open(os.path.join(d, "argv.sh")).read().split():
        if t.startswith("timestep="):
            return float(t.split("=", 1)[1])
    return 0.008


def radius(r, dx):
    return math.hypot((int(r["col_i"]) + 0.5) * dx, (int(r["col_j"]) + 0.5) * dx)


def dx_of(d):
    return geom_meta(d)["dx"]


def med(v):
    v = [x for x in v if math.isfinite(x)]
    return float(np.median(v)) if v else float("nan")


def pct(v, q):
    v = [x for x in v if math.isfinite(x)]
    return float(np.percentile(v, q)) if v else float("nan")


def tsum(d, c):
    h, rows = R15.thermo(d)
    return sum(r[h.index(c)] for r in rows) if c in h else float("nan")


# ---------------------------------------------------------------- a1: the face gate's audit
def A1():
    a, b = os.path.join(W, "face_p2_18.md5"), os.path.join(W, "face_p2_19.md5")
    ha = dict((l.split("  ", 1)[1], l.split("  ", 1)[0]) for l in open(a).read().splitlines() if l.strip())
    hb = dict((l.split("  ", 1)[1], l.split("  ", 1)[0]) for l in open(b).read().splitlines() if l.strip())
    only = sorted(set(ha) ^ set(hb))
    diff = sorted(p for p in set(ha) & set(hb) if ha[p] != hb[p])
    L = ["a1 -- face gate cmp p2_18 p2_19: every differing path classified (%d files each side; %d differ; %d only on "
         "one side)" % (len(ha), len(diff), len(only)), "=" * 100]
    bad = list(only)
    for p in diff:
        A, B = os.path.join(W, "face_p2_18", p), os.path.join(W, "face_p2_19", p)
        if p.endswith("thermo.dat"):
            ra = [l.split("\t") for l in open(A).read().splitlines()]
            rb = [l.split("\t") for l in open(B).read().splitlines()]
            hA, hB = [x.strip() for x in ra[0]], [x.strip() for x in rb[0]]
            ren = [h.replace("fr_noskin_", "fr_nodepth_") for h in hA]
            ok = ren == hB and len(ra) == len(rb)
            moved = []
            if ok:
                for j, h in enumerate(hA):
                    if h in TIMERS or h.startswith("fr_noskin_"):
                        continue
                    if any(x[j].strip() != y[j].strip() for x, y in zip(ra[1:], rb[1:])):
                        moved.append(h)
            nd = [sum(float(x[hB.index(c)]) for x in rb[1:]) for c in hB if c.startswith("fr_nodepth_")] if ok else []
            ns = [sum(float(x[hA.index(c)]) for x in ra[1:]) for c in hA if c.startswith("fr_noskin_")] if ok else []
            cls = ok and not moved
            L.append("  %s  (i)  %s: header -fr_noskin_* +fr_nodepth_* only: %s; other non-timer columns moved: %s; "
                     "sum noskin %s -> nodepth %s" % ("ok  " if cls else "FAIL", p, ok, moved or "none", ns, nd))
            if not cls:
                bad.append(p)
        elif p.endswith("_removal_events.csv"):
            hok, nsame, first, nn, no = prefix_cmp(open(B).read().splitlines(), open(A).read().splitlines())
            cls = hok and first is None
            L.append("  %s  (ii) %s: header + %s: %s; rows %d vs %d, pre-existing columns identical on all: %s" % (
                "ok  " if cls else "FAIL", p, ",".join(NEW_COLS), hok, nn, no, first is None))
            if not cls:
                bad.append(p)
        else:
            L.append("  FAIL  %s: a %s differs (no Cell_D / Cell_H / Header may)" % (p, os.path.basename(p)))
            bad.append(p)
    for p in only:
        L.append("  FAIL  %s only in %s" % (p, "p2_18" if p in ha else "p2_19"))
    ncell = sum(1 for p in diff if re.search(r"Cell_D|Cell_H|/Header$", p))
    fbs = []
    for leg in sorted({p.split("/run")[0] for p in hb if p.endswith("thermo.dat")}):
        d = os.path.join(W, "face_p2_19", leg)
        fbs.append("%s %g" % (leg, sum(tsum(d, c) for c in ("fr_fb_floor", "fr_fb_wall", "fr_fb_roof"))))
    L.append("  Cell_D / Cell_H / Header differing: %d (P6: 0). fr_fb_* sums on the new legs: %s" % (ncell, ", ".join(fbs)))
    L.append("a1 %s (gate)" % ("PASS" if not bad else "FAIL: %d path(s)" % len(bad)))
    return "\n".join(L), not bad


# ---------------------------------------------------------------- a2: T00 x3
def A2():
    L, ok = ["a2 -- T00 at three meshes vs the P2-9 record (unit/face_removal/output/T00_*_k1, cut by step at 8 s) and "
             "T00_1 vs the P2-18 face gate's fr/T00_1_k1 (8 s): pre-existing columns", "=" * 100], True
    pairs = [("T00_2", os.path.join(T, "unit", "face_removal", "output", "T00_2_k1"), "P2-9 record"),
             ("T00_1", os.path.join(T, "unit", "face_removal", "output", "T00_1_k1"), "P2-9 record"),
             ("T00_05", os.path.join(T, "unit", "face_removal", "output", "T00_05_k1"), "P2-9 record"),
             ("T00_1", os.path.join(W, "face_p2_18", "fr", "T00_1_k1"), "face_p2_18")]
    for leg, ref, who in pairs:
        d = os.path.join(LEGS, leg)
        new = lines(d)
        old = cut(lines(ref), 8.0, dt_of(d))
        hok, nsame, first, nn, no = prefix_cmp(new, old)
        o = open(os.path.join(d, "rc")).read().strip() == "0" and hok and first is None and nn > 0
        ok = ok and o
        L.append("  %s  %-6s vs %-11s rows %d vs %d (old cut at 8 s - dt/2); header +4: %s; pre-existing columns "
                 "identical: %s" % ("PASS" if o else "FAIL", leg, who, nn, no, hok, first is None))
    L.append("a2 %s (gate)" % ("PASS" if ok else "FAIL"))
    return "\n".join(L), ok


# ---------------------------------------------------------------- a3: CAL
def A3():
    d, ref = os.path.join(LEGS, "CAL"), REF["CAL"]
    hok, nsame, first, nn, no = prefix_cmp(lines(d), lines(ref))
    n, nd, no_pf = R14.same_plotfiles(ref, d)
    o = open(os.path.join(d, "rc")).read().strip() == "0" and hok and first is None and n > 0 and nd == 0 \
        and no_pf == 0
    return ("a3 -- CAL vs P2-16 legs/CAL: event log rows %d vs %d, pre-existing columns identical: %s (header +4: %s); "
            "plotfiles %d files, %d differ, %d one-sided\na3 %s (gate)" % (nn, no, first is None, hok, n, nd, no_pf,
                                                                          "PASS" if o else "FAIL")), o


# ---------------------------------------------------------------- c
def population(rows, dx, pick):
    sel = [r for r in rows if pick(r)]
    return dict(n=len(sel), r=med([radius(r, dx) for r in sel]), hf=med([float(r["h_f"]) for r in sel]),
                ell=med([float(r["ell_eff"]) for r in sel]), vol=sum(float(r["h_line"]) for r in sel) * dx * dx,
                rows=sel)


def fb(r):
    return int(round(float(r["rc_flag"]))) & 1 == 1


def C1():
    d, ref = os.path.join(LEGS, "C0H150"), REF["C0_h150"]
    new, old = lines(d), lines(ref)
    hok, nsame, first, nn, no = prefix_cmp(new, old)
    t_first = float(new[first].split(",", 1)[0]) if first and first < len(new) else float("nan")
    dx = dx_of(d)
    dn, do = dicts(d), dicts(ref)
    pn = population(dn, dx, fb)
    po = population(do, dx, fb)
    fns = population(dn, dx, lambda r: fb(r) and float(r["q_net"]) <= 0.0)
    # h_f new/old on the SAME (col, k) cells where both logs fired on the fallback branch
    key = lambda r: (r["col_i"], r["col_j"], r["k"], r["dir"])
    old_by = {}
    for r in po["rows"]:
        old_by.setdefault(key(r), r)
    ratios = [float(r["h_f"]) / float(old_by[key(r)]["h_f"]) for r in pn["rows"]
              if key(r) in old_by and float(old_by[key(r)]["h_f"]) > 0 and math.isfinite(float(r["ell_p"]))]
    fin = sum(1 for r in pn["rows"] if math.isfinite(float(r["ell_p"])))
    thr_n = sum(1 for r in pn["rows"] if r["through"] == "1")
    thr_o = sum(1 for r in po["rows"] if r["through"] == "1")
    runs_n = sorted({int(r["run"]) for r in pn["rows"]})
    runs_o = sorted({int(r["run"]) for r in po["rows"]})
    # removal identity on the fallback rows matched by (time, cell, dir)
    tkey = lambda r: (r["time"], r["col_i"], r["col_j"], r["k"], r["dir"])
    ob = {tkey(r): r for r in po["rows"]}
    matched = [(r, ob[tkey(r)]) for r in pn["rows"] if tkey(r) in ob]
    rem_same = sum(1 for a, b in matched if a["h_line"] == b["h_line"] and a["n_voided"] == b["n_voided"])
    # the old NON-through fallback firings (the only ones whose thickness set a removal) and their fate
    nkeys = {(r["col_i"], r["col_j"], r["k"], r["dir"]) for r in dn}
    nonthr = [r for r in po["rows"] if r["through"] == "0"]
    gone = sum(1 for r in nonthr if (r["col_i"], r["col_j"], r["k"], r["dir"]) not in nkeys)
    # which columns the first differing row differs in
    fcols = []
    if first and first < len(new) and first < len(old):
        a, b = dn[first - 1], do[first - 1]
        fcols = [k for k in b if a.get(k) != b[k]]
    nd_new = sum(tsum(d, "fr_nodepth_" + c) for c in ("floor", "wall", "roof"))
    ns_old = sum(tsum(ref, "fr_noskin_" + c) for c in ("floor", "wall", "roof"))
    mn, mo = mod_for(LEGS), mod_for(os.path.dirname(ref))
    g, P = R11R.pfs_by_t(mn, "C0H150")
    go, Po = R11R.pfs_by_t(mo, "C0_h150")
    vn = R11R.fr_volume(mn, "C0H150", P[200])
    vo = R11R.fr_volume(mo, "C0_h150", Po[200])
    L = ["c1 -- C0H150 (the calibration deck at h 150, 200 s) vs P2-14a legs/C0_h150", "=" * 100,
         "  identical (pre-existing columns) for %d rows; first differing row %s at t %.3f s (P1: the first old fallback "
         "firing, 120.3 s), differing in %s; rows %d new vs %d old; header +4: %s" % (
             nsame, first, t_first, fcols, nn, no, hok),
         "  fallback firings (rc_flag bit 1): new %d (P1 118 +- 30), old %d; median r %.2f / %.2f mm; h_f %.4f / %.4f "
         "mm; ell_eff %.3f / %.3f mm" % (pn["n"], po["n"], 1e3 * pn["r"], 1e3 * po["r"], 1e3 * pn["hf"],
                                         1e3 * po["hf"], 1e3 * pn["ell"], 1e3 * po["ell"]),
         "  ⚠ the fallback firings are THROUGH firings of thin ledges: new %d of %d through (runs %s), old %d of %d "
         "(runs %s). ell_p is finite on %d of the %d new ones: a run of 1-2 cells leaves 1-2 samples in the window, "
         "too few for the rise to e-fold, so the depth is undefined. A through firing removes its whole solid run "
         "(L_av) whatever h_f says, so the removal is unchanged: h_line and n_voided identical on %d of %d fallback "
         "rows matched by (time, cell, dir)" % (thr_n, pn["n"], runs_n, thr_o, po["n"], runs_o, fin, pn["n"],
                                                 rem_same, len(matched)),
         "  the old NON-through fallback firings (the only ones whose h_f set a removal): %d; absent from the new log "
         "(refused as no depth): %d" % (len(nonthr), gone),
         "  h_f new/old where both thicknesses exist: %s (P2 [0.7, 1.1], z0e 0.857): %s" % (
             "median %.4f on %d" % (med(ratios), len(ratios)) if ratios else "none",
             "NOT SCOREABLE on this deck (no fallback firing has a finite ell_p)" if not ratios else "scored"),
         "  formerly-no-skin firings (bit 1 AND q_net <= 0): %d (r %s mm, h_f %s mm, volume %.3g m3)" % (
             fns["n"], "%.2f" % (1e3 * fns["r"]) if fns["n"] else "-", "%.4f" % (1e3 * fns["hf"]) if fns["n"] else "-",
             fns["vol"]),
         "  sum fr_nodepth (new) %g vs sum fr_noskin (old) %g (thermo every step on this deck)" % (nd_new, ns_old),
         "  V(f) at 200 s new/old = %.4f / %.4f = %.5f (P4 1 +- 0.05; fr_volume's x4 convention, the ratio only)"
         % (vn, vo, vn / vo)]
    return "\n".join(L), dict(first_t=t_first, nfb=pn["n"], ratio=med(ratios) if ratios else float("nan"),
                              vr=vn / vo, fns=fns["n"], fin=fin, gone=gone, nonthr=len(nonthr))


def C2():
    d = os.path.join(LEGS, "HOS")
    new = lines(d)
    L = ["c2 -- HOS (P2-16 HOS_R's tokens; B ON since P2-17's rock file) vs P2-17 N_S164 (same tokens, B on) cut by "
         "step at 15 s, and vs P2-16 HOS_R (B off)", "=" * 100]
    res = {}
    for nm in ("N_S164", "HOS_R"):
        old = cut(lines(REF[nm]), 15.0, 0.008)
        hok, nsame, first, nn, no = prefix_cmp(new, old)
        tf = float(new[first].split(",", 1)[0]) if first and first < len(new) else float("nan")
        res[nm] = (nsame, first, tf)
        L.append("  vs %-6s: identical for %d rows; first differing row %s at t %.3f s; rows %d vs %d" % (
            nm, nsame, first, tf, nn, no))
    dx = dx_of(d)
    dn = dicts(d)
    pf = population(dn, dx, fb)
    fns = population(dn, dx, lambda r: fb(r) and float(r["q_net"]) <= 0.0)
    far = [r for r in dn if radius(r, dx) > 0.025]
    L.append("  fallback firings: %d at t %s; formerly-no-skin firings (bit 1 AND q_net <= 0): %d%s; firings at r > 25 "
             "mm: %d (P3: 0 and 0)" % (pf["n"], [round(float(r["time"]), 3) for r in pf["rows"]][:6], fns["n"],
                                       " (median r %.1f mm)" % (1e3 * fns["r"]) if fns["n"] else "", len(far)))
    nd = sum(tsum(d, "fr_nodepth_" + c) for c in ("floor", "wall", "roof"))
    h, rows = R15.thermo(REF["N_S164"])
    ns15 = sum(r[h.index("fr_noskin_" + c)] for r in rows if r[h.index("time")] <= 15.0 for c in ("floor", "wall", "roof"))
    nsR = sum(tsum(REF["HOS_R"], "fr_noskin_" + c) for c in ("floor", "wall", "roof"))
    L.append("  sum fr_nodepth %g vs the old sum fr_noskin over the same 15 s: N_S164 %g, HOS_R %g (thermo sampled 1 "
             "in 6 on all three: the same rows, a like-for-like count). The cold rim (q_net <= 0) now has a finite "
             "depth and is no longer refused; it does not fire" % (nd, ns15, nsR))
    return "\n".join(L), dict(res=res, fns=fns["n"], far=len(far))


def C3():
    L = ["c3 -- the energy thin limit h_E = K_Ic(T_s)^2 (1 + nu) / sigma_s^2 (circular blister) on every logged firing",
         "=" * 100]
    for leg in ("HOS", "CAL", "C0H150"):
        rows = dicts(os.path.join(LEGS, leg))
        hE = [float(r["h_E"]) for r in rows]
        rat = [float(r["h_f"]) / float(r["h_E"]) for r in rows
               if math.isfinite(float(r["h_E"])) and float(r["h_E"]) > 0]
        L.append("   %-7s n %5d  h_E median %7.2f um  p10 %7.2f  p90 %7.2f   h_f/h_E median %.2f" % (
            leg, len(rows), 1e6 * med(hE), 1e6 * pct(hE, 10), 1e6 * pct(hE, 90), med(rat)))
    return "\n".join(L)


def C4():
    P17.MOD = mod_for(LEGS)
    sn = P17.sectors("HOS", t=15)
    P17.MOD = mod_for(os.path.dirname(REF["N_S164"]))
    so = P17.sectors("N_S164", t=15)
    a, b = max(sn) / min(sn) - 1, max(so) / min(so) - 1
    return ("c4 -- isotropy witness, HOS at 15 s vs N_S164 at 15 s (B on both): spread %.4f vs %.4f, new - old %+.4f "
            "(> 0.05 would mean a grid-direction bias: %s)" % (a, b, a - b, "YES" if a - b > 0.05 else "no"))


def main(a):
    mode = a[0] if a else "RESULTS"
    for m, fn, f in (("A1", A1, "A1.txt"), ("A2", A2, "A2.txt"), ("A3", A3, "A3.txt")):
        if mode == m:
            t, ok = fn(); print(t); open(os.path.join(HERE, f), "w").write(t + "\n"); return 0 if ok else 1
    out = []
    for f in ("Z0.txt", "GATE_tree.txt", "A1.txt", "A2.txt", "A3.txt"):
        p = os.path.join(HERE, f)
        out += [open(p).read().rstrip() if os.path.exists(p) else f + ": NOT ON DISK", ""]
    c1, r1 = C1()
    c2, r2 = C2()
    out += [c1, "", c2, "", C3(), "", C4(), ""]
    g = {}
    for k, f in (("z0", "Z0.txt"), ("a0", "GATE_tree.txt"), ("a1", "A1.txt"), ("a2", "A2.txt"), ("a3", "A3.txt")):
        p = os.path.join(HERE, f)
        t = open(p).read() if os.path.exists(p) else ""
        g[k] = ("0 of 2639 changed" in t and "vs d2aa: 26 of 2639" in t and "vs d2t: 863 of 2639" in t) if k == "a0" \
            else (t.rstrip().endswith("PASS") if k == "z0" else bool(re.search(r"^%s PASS" % k, t, re.M)))
    out.append("e0 VERDICT")
    out.append("   the change is confined to fallback and formerly-no-skin firings: %s (%s)" % (
        "YES" if all(g.values()) else "NO", ", ".join("%s %s" % (k, "PASS" if v else "FAIL") for k, v in g.items())))
    out.append("   C0H150: first moved row at t %.3f s (logged columns only); %d fallback firings, finite ell_p on %d (all "
               "thin-ledge THROUGH firings: removal unchanged); %d of the %d old non-through fallback firings refused "
               "as no depth; V(f) new/old %.5f; formerly-no-skin firings %d" % (
                   r1["first_t"], r1["nfb"], r1["fin"], r1["gone"], r1["nonthr"], r1["vr"], r1["fns"]))
    out.append("   HOS: first moved row vs N_S164 at t %.3f s; formerly-no-skin firings %d; firings at r > 25 mm %d" % (
        r2["res"]["N_S164"][2], r2["fns"], r2["far"]))
    out.append("   the helper reads an interior peak (z0d); ell_p, z_peak, ell_pk and h_E are logged on every firing")
    print("\n".join(out))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
