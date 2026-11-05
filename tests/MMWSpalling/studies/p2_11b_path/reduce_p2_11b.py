#!/usr/bin/env python3
"""P2-11b -- reductions. Identity work uses P2-10's comparator (`studies/p2_10_retire/run_p2_10.py`:
tree_files / compare_dirs / thermo_cmp, its exclusions unchanged), loaded as a module with its output folder
pointed HERE, so no P2-10 record is touched. Shape, volume and plotfile reads are P2-11's reducer
(`studies/p2_11_heat/reduce_p2_11.py`), loaded twice, as P2-11 did with P2-8's.

  reduce_p2_11b.py W          a3 (second half): the witness set W on this binary into W_p2_11b/, vs P2-10's W_c2
  reduce_p2_11b.py IDENTITY   a4 (G1_C vs G1), a4b (P47n_C vs P47n, binary B)        -> IDENTITY.txt
  reduce_p2_11b.py A6         a6: the X legs' first plotfile                          -> A6.txt
  reduce_p2_11b.py SSTAR      b1                                                       -> SSTAR.txt
  reduce_p2_11b.py            c1-c5, e1 (stdout; the packet redirects it into RESULTS.txt)
"""
import csv
import importlib.machinery
import importlib.util
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)
ROOT = os.path.abspath(os.path.join(STUD, "..", "..", ".."))


def _load(path, name):
    ld = importlib.machinery.SourceFileLoader(name, path)
    m = importlib.util.module_from_spec(importlib.util.spec_from_loader(name, ld))
    ld.exec_module(m)
    return m


R10 = _load(os.path.join(STUD, "p2_10_retire", "run_p2_10.py"), "r10")
R11 = _load(os.path.join(STUD, "p2_11_heat", "reduce_p2_11.py"), "r11")
P11 = R11.P11                       # P2-11's legs (G1, binary C2)
P9 = R11.P9                         # P2-9's legs (P47n, binary B)
PB = R11._load("r8_p211b", HERE)    # this study's legs
H0G = 7272.3958
TARGET_V = 33.0
# P2-10's retired thermo set: the 67 columns M0 carried and M0T did not (P2-10 GATE_a1.txt line 7) and the
# columns P2-10 dropped from face-list runs (run_p2_10.DECLARED["c2"]); the 11 lateral event-log columns
RETIRED_COLS = set(open(os.path.join(STUD, "p2_10_retire", "GATE_a1.txt")).read().split(
    "columns only in M0 (67): ", 1)[1].split("\n", 1)[0].split()) | R10.DECLARED["c2"]
assert len(RETIRED_COLS) == 67 + 8, len(RETIRED_COLS)


# ------------------------------------------------------------------ a3 (W) ----
def W():
    label = "p2_11b"
    R10.HERE = HERE                 # outputs into this study
    R10.refs(label)
    A = os.path.join(STUD, "p2_10_retire", "W_c2")
    B = os.path.join(HERE, "W_" + label)
    lines = [f"a3 (W) -- W_{label} (this binary, face_heat unset) vs P2-10's W_c2 (binary C2), P2-10's comparator and "
             "exclusions (metadata, diff.patch, run.log, argv/BINARY/DROPPED/.done; thermo timers)", "=" * 100]
    ok = True
    for sub in ["M0T"] + ["fr/" + x for x in R10.FR_LEGS] + ["fl/" + x for x in R10.FL_LEGS]:
        o, th = R10.compare_dirs(os.path.join(A, sub), os.path.join(B, sub), lines, sub)
        ok = ok and o
        for rel in th:
            c = R10.thermo_cmp(os.path.join(A, sub, rel), os.path.join(B, sub, rel))
            ok = ok and not c["only_a"]
    lines.append(f"W_{label} vs W_c2: {'PASS' if ok else 'FAIL'}")
    txt = "\n".join(lines)
    print(txt)
    open(os.path.join(HERE, f"W_{label}_vs_c2.txt"), "w").write(txt + "\n")


# ------------------------------------------------------------------ a4, a4b ----
def events_cmp(pa, pb):
    """the event logs on their common columns: (rows a, rows b, common, only a, only b, rows differing)"""
    ra = list(csv.reader(open(pa)))
    rb = list(csv.reader(open(pb)))
    ha, hb = ra[0], rb[0]
    common = [c for c in ha if c in hb]
    ia = [ha.index(c) for c in common]
    ib = [hb.index(c) for c in common]
    nd = 0
    if len(ra) == len(rb):
        for x, y in zip(ra[1:], rb[1:]):
            if [x[i] for i in ia] != [y[i] for i in ib]:
                nd += 1
    return len(ra) - 1, len(rb) - 1, len(common), [c for c in ha if c not in hb], [c for c in hb if c not in ha], nd


def identity():
    L = ["P2-11b (a4, a4b) -- the controls on the Stage B binary", "=" * 100]
    # a4: G1_C vs G1, every file (P2-10's exclusions: metadata, diff.patch, run.log, argv/DROPPED/.done)
    lines = []
    ok4, _ = R10.compare_dirs(P11.legdir("G1"), PB.legdir("G1_C"), lines, "G1 vs G1_C")
    L += ["a4  " + lines[0]] + lines[1:]
    L.append(f"a4  {'PASS' if ok4 else 'FAIL'}")
    # a4b: P47n_C vs P47n (binary B)
    a, b = P9.legdir("P47n"), PB.legdir("P47n_C")
    fa, fb = R10.tree_files(a), R10.tree_files(b)
    pfs = sorted(r for r in set(fa) | set(fb) if "/Level_" in r or r.endswith("/Header"))
    pdiff = [r for r in pfs if r not in fa or r not in fb or R10.md5(fa[r]) != R10.md5(fb[r])]
    L.append(f"a4b plotfiles: {len(pfs)} Cell_D_* / Header files, {len(pdiff)} differ or missing"
             + (": " + ", ".join(pdiff[:8]) if pdiff else ""))
    na, nb, nc, oa, ob, nd = events_cmp(os.path.join(a, "run_removal_events.csv"),
                                        os.path.join(b, "run_removal_events.csv"))
    L.append(f"a4b event log: rows {na} / {nb}, {nc} common columns, {nd} rows differ on them; only in P47n {oa}; "
             f"only in P47n_C {ob or 'none'}")
    c = R10.thermo_cmp(os.path.join(a, "run", "thermo.dat"), os.path.join(b, "run", "thermo.dat"))
    extra = sorted(set(c["only_a"]) - RETIRED_COLS)
    L.append(f"a4b thermo: rows {c['rows']}, common non-timer columns differing {c['diff_cols'] or 'none'}, order "
             f"same {c.get('order_same')}, only in P47n_C {c['only_b'] or 'none'}; only in P47n {len(c['only_a'])}, "
             f"outside P2-10's retired set {extra or 'none'}")
    lat = set(oa) <= set(R10.LAT_CSV)
    ok4b = (not pdiff and na == nb and nd == 0 and lat and not ob and c["rows"][0] == c["rows"][1]
            and not c["diff_cols"] and not c["only_b"] and not extra and c.get("order_same", False))
    L.append(f"a4b {'PASS' if ok4b else 'FAIL'} (report)")
    txt = "\n".join(L)
    print(txt)
    open(os.path.join(HERE, "IDENTITY.txt"), "w").write(txt + "\n")


# ------------------------------------------------------------------ the X1 formula ----
D = 0.0079


def x1_h(S):
    """the X legs' h_expr at (r := s_path, s := H), in the parser's left-to-right order"""
    C = float("%.6f" % (S * H0G))

    def h(s, H):
        p = math.pow(s / 0.0079, 1.285)
        return C * 533 * (533 + 0.3575 * 44 * p) / math.pow(533 + 44 * p, 2) * math.pow(0.0632 / max(0.0632, H), 1.3185305455)
    return h


def leg_S(n):
    e = R11.argv_kv(PB, n)["surface_patch.h_expr"]
    return float(e.split("*", 1)[0]) / H0G


def thermo_rows(mod, n):
    hdr, rows = mod.thermo(mod.legdir(n))
    return hdr, [[float(x) for x in r] for r in rows]


def H_of_pf(mod, n, pf):
    """the stagnation standoff the step used: z_n - top face of the axis column (the innermost bin's only
    column at dx = dr = 1 mm), from the plotfile's mask"""
    d = mod.legdir(n)
    g = mod.geom_of(d)
    rem = mod.read_pf(d, pf, g, comps=("removed",))["removed"]
    ks = np.nonzero(~(rem[0, 0, :] > 0.5))[0]
    zn = float(R11.argv_kv(mod, n)["surface_patch.nozzle_z0"])
    return zn - (ks.max() + 1) * g["dz"]


# ------------------------------------------------------------------ a6 ----
def a6():
    L = ["P2-11b (a6) -- the X legs' first plotfile after t = 0: fl_h vs the X1 formula at the logged fl_s and the "
         "step's H (= z_n - the axis column's top face, read from the same plotfile; thermo's fl_H_stag beside it)",
         "=" * 110]
    ok = True
    for n in ("X1_S14", "X1_S20", "X1_S10", "X1_S28", "X1"):
        if not PB.done(n):
            continue
        g, P = R11.pfs_by_t(PB, n)
        t = min(k for k in P if k > 0)
        d = PB.legdir(n)
        f = PB.read_pf(d, P[t], g, comps=("fl_h", "fl_w", "fl_s"))
        H = H_of_pf(PB, n, P[t])
        h = x1_h(leg_S(n))
        lst = f["fl_w"] > 0
        rel = []
        for (i, j, k) in np.argwhere(lst):
            s = f["fl_s"][i, j, k]
            want = h(s, H) if s >= 0 else 0.0
            got = f["fl_h"][i, j, k]
            rel.append(0.0 if want == got else (abs(got / want - 1) if want > 0 else float("inf")))
        rel = np.array(rel)
        share = float((rel <= 1e-12).mean())
        hdr, rows = thermo_rows(PB, n)
        it = hdr.index("time")
        near = min(rows, key=lambda r: abs(r[it] - t))
        src, blk = near[hdr.index("fl_n_src")], near[hdr.index("fl_n_blocked")]
        Hth = near[hdr.index("fl_H_stag")]
        x = (np.arange(g["nx"]) + 0.5) * g["dx"]
        R = np.hypot(x[:, None], x[None, :])[:, :, None] * np.ones((1, 1, g["nz"]))
        depth = (0.15 - (np.arange(g["nz"]) + 0.5) * g["dz"])[None, None, :] * np.ones(lst.shape)
        deep = lst & (depth > 0.005) & (f["fl_s"] >= 0)
        sr = float(np.median(f["fl_s"][deep] / R[deep])) if deep.any() else float("nan")
        o = share >= 0.99 and src > 0 and blk > 0 and sr > 1
        ok &= o
        L.append(f"  {n:7s} t {t:3d} s: listed {lst.sum()}, |rel| <= 1e-12 on {share:.4f} (max {rel.max():.2e}); "
                 f"H {1e3 * H:.1f} mm (thermo fl_H_stag at t {near[it]:g}: {1e3 * Hth:.1f}); fl_n_src {src:g}, "
                 f"fl_n_blocked {blk:g}; median fl_s/r over listed cells > 5 mm deep {sr:.3f} (n {deep.sum()})  "
                 f"{'PASS' if o else 'FAIL'}")
    L.append(f"a6 {'PASS' if ok else 'FAIL'}")
    txt = "\n".join(L)
    print(txt)
    open(os.path.join(HERE, "A6.txt"), "w").write(txt + "\n")


# ------------------------------------------------------------------ b1 ----
def sstar():
    L = ["P2-11b (b1) -- S*_X1 by the rule: log V(f)(30 s, x4) vs log S between X1_S14 and X1_S20 to 33.0 cm3",
         "=" * 100]
    pts = []
    for n in ("X1_S10", "X1_S14", "X1_S20", "X1_S28"):
        if PB.done(n):
            g, P = R11.pfs_by_t(PB, n)
            v = R11.fr_volume(PB, n, P[30])
            pts.append((round(leg_S(n), 4), v, n))
            L.append(f"  {n:7s} S {leg_S(n):.4f}  V(f) 30 s {v:.3f} cm3")
    pts.sort()
    star = None
    for (s1, v1, n1), (s2, v2, n2) in zip(pts[:-1], pts[1:]):
        if (v1 - TARGET_V) * (v2 - TARGET_V) <= 0:
            x = math.log(s1) + (math.log(TARGET_V) - math.log(v1)) * \
                (math.log(s2) - math.log(s1)) / (math.log(v2) - math.log(v1))
            star = round(math.exp(x), 4)
            L.append(f"  bracket {n1}-{n2}: S* = {math.exp(x):.6f} -> {star:.4f}   (P1 predicted 0.17-0.21)")
            break
    if star is None:
        if any(PB.done(n) for n in ("X1_S10", "X1_S28")):
            L.append("  NOT BRACKETED after the one added point: STOP and report (b1 FAIL)")
        else:
            lo = max(p[1] for p in pts) if pts else None
            L.append("  NOT BRACKETED -- the rule adds ONE point: %s" %
                     ("X1_S28 (S 0.28)" if lo is not None and lo < TARGET_V else "X1_S10 (S 0.10)"))
        L.append("b1 not resolved")
    else:
        L.append(f"  C* = {star * H0G:.6f}")
        L.append(f"b1 PASS (interior, {len(pts) - 2} added point)")
        L.append(f"S_STAR={star:.4f}")
    txt = "\n".join(L)
    print(txt)
    open(os.path.join(HERE, "SSTAR.txt"), "w").write(txt + "\n")


# ------------------------------------------------------------------ c1-c5 ----
HER_V = {0: 0.0, 15: 16.0, 30: 33.0, 60: 81.0, 90: 135.0}
WIN = [(0, 15), (15, 30), (30, 60), (60, 90), (90, 120)]


def shares(mod, n, t):
    """realised conductance shares fl_h x fl_w: within r < 10 mm; the top-10-mm WALL band (listed, smoothed
    n_z <= cos 45, cell centre <= 10 mm below the original top) -- predict_x1.py's definitions"""
    g, P = R11.pfs_by_t(mod, n)
    f = mod.read_pf(mod.legdir(n), P[t], g, comps=("fr_f", "fl_h", "fl_w"))
    nz, _, zero = R11.smoothed_nz(f["fr_f"])
    lst = f["fl_w"] > 0
    G = f["fl_h"] * f["fl_w"]
    x = (np.arange(g["nx"]) + 0.5) * g["dx"]
    R = np.hypot(x[:, None], x[None, :])[:, :, None] * np.ones((1, 1, g["nz"]))
    depth = (0.15 - (np.arange(g["nz"]) + 0.5) * g["dz"])[None, None, :] * np.ones(lst.shape)
    floor = lst & (nz > math.cos(math.radians(45))) & ~zero
    band = lst & ~floor & (depth <= 0.010)
    tot = G[lst].sum()
    return G[lst & (R < 0.010)].sum() / tot, G[band].sum() / tot


def offline_ratios():
    """a1/a2 primary ratios from OFFLINE_X1.txt"""
    out = {}
    on = False
    for ln in open(os.path.join(HERE, "OFFLINE_X1.txt")):
        if ln.strip().startswith("PRIMARY"):
            on = True
            continue
        if ln.strip().startswith("SECONDARY"):
            break
        if on and "|" in ln and ln.split("|")[0].strip().isdigit():
            p = ln.split("|")
            t = int(p[0])
            out[t] = (float(p[1].split()[2]), float(p[2].split()[2]))
    return out


def growth(Vt):
    rates = [((Vt[b] - Vt[a]) / ((b - a) / 60.0)) if (a in Vt and b in Vt) else float("nan") for a, b in WIN]
    return rates, rates[3] / rates[0]


def c4_rows(n, S):
    """lateral firings (dir != +z) per 15 s window with h_cal = fl_h at the firing cell / S < 500. The event log
    carries no h, so fl_h is read at the cell from the plotfile at the window's start, or its end if the cell is
    not listed at the start; firings at cells listed in neither are counted apart."""
    g, P = R11.pfs_by_t(PB, n)
    ev = PB.events(PB.legdir(n))
    d = PB.legdir(n)
    cache = {}

    def F(t):
        if t not in cache:
            cache[t] = PB.read_pf(d, P[t], g, comps=("fl_h", "fl_w")) if t in P else None
        return cache[t]
    rows = []
    for a in range(0, int(round(g["st"])), 15):
        hc, miss = [], 0
        for e in ev:
            te = float(e["time"])
            if not (a < te <= a + 15) or e.get("dir") == "5":
                continue
            i, j, k = int(e["col_i"]), int(e["col_j"]), int(e["k"])
            v = None
            for tt in (a, a + 15):
                f = F(tt)
                if f is not None and f["fl_w"][i, j, k] > 0:
                    v = f["fl_h"][i, j, k]
                    break
            if v is None:
                miss += 1
                continue
            hc.append(v / S)
        hc = np.array(hc)
        rows.append((a, len(hc) + miss, miss, int((hc < 500).sum()) if len(hc) else 0,
                     float(np.median(hc)) if len(hc) else float("nan"), float(hc.min()) if len(hc) else float("nan")))
    return rows


def results():
    P = print
    P("P2-11b -- X1 (the surface-path heat rule) on Hoeser: c1-c5 (binary Stage B; G1 is P2-11's, binary C2 -- a4 "
      "shows G1_C reproduces it; P47n is P2-9's, binary B)")
    P("=" * 116)
    if not PB.done("X1"):
        P("  X1 not done")
        return
    S = leg_S("X1")
    A = R11.series(PB, "X1")
    G = R11.series(P11, "G1")
    B = R11.series(P9, "P47n")
    HER = dict(zip([15, 30, 60, 90, 120], zip(P9.HER_DEPTH, P9.HER_DMAX)))
    P(f"\n  c1  X1 (S* {S:.4f}) vs G1 (S* 0.1636) vs Hoeser 47 MPa. depth = fr_f axis depth (4-line spread); max Ø = mask")
    P("        t |  depth X1 (spr)    G1 (spr)  her  X1/G1  X1/her |  maxØ X1   G1  her  X1-G1 |  V(f) X1     G1   her")
    agree = True
    for t in (15, 30, 60, 90, 120):
        a, b = A[t], G[t]
        hd, hm = HER[t]
        if abs(a["depth"] / hd - 1) > 0.15 or (t >= 60 and abs(a["dmax"] - hm) > 6):
            agree = False
        P(f"     {t:4d} | {a['depth']:7.2f} ({a['spread']:.2f}) {b['depth']:7.2f} ({b['spread']:.2f}) {hd:4.0f} "
          f"{a['depth'] / b['depth']:6.3f} {100 * (a['depth'] / hd - 1):+6.1f}% | {a['dmax']:6.1f} {b['dmax']:5.1f} "
          f"{hm:4.0f} {a['dmax'] - b['dmax']:+6.1f} | {a['V']:7.2f} {b['V']:7.2f} {HER_V.get(t, float('nan')):5.1f}")
    for nm, X in (("X1", A), ("G1", G), ("P47n", B)):
        gr = (X[120]["depth"] - X[15]["depth"]) / max(X[120]["dmax"] - X[15]["dmax"], 1e-9)
        P(f"    {nm:5s} depth/diameter growth 15->120 s {gr:.2f} (hers 1.8): depth {X[15]['depth']:.1f} -> "
          f"{X[120]['depth']:.1f} mm, max Ø {X[15]['dmax']:.1f} -> {X[120]['dmax']:.1f} mm")
    P(f"    X1 SHAPE AGREES with Hoeser (max Ø within 6 mm at 60/90/120 s AND depth within 15 % every time): "
      f"{'YES' if agree else 'NO'}")

    P("\n  c2  incremental removal rate per window [cm3/min] (V(f), x4) and the 60-90 / 0-15 ratio; absorbed power "
      "ledger_P_robin x 4 [kW]")
    P("        leg | " + " ".join(f"{a:3d}-{b:3d}" for a, b in WIN) + " | 60-90/0-15 | P_abs 15 / 30 / 60 / 90 / 120 s")
    for nm, X, mod, leg in (("X1", A, PB, "X1"), ("G1", G, P11, "G1"), ("P47n", B, P9, "P47n"), ("Hoeser", None, None, None)):
        Vt = {0: 0.0, **{t: X[t]["V"] for t in X}} if X else HER_V
        rates, ratio = growth(Vt)
        pw = ""
        if mod is not None:
            hdr, rows = thermo_rows(mod, leg)
            it, ip = hdr.index("time"), hdr.index("ledger_P_robin")
            pw = " / ".join("%.2f" % (4e-3 * min(rows, key=lambda r: abs(r[it] - t))[ip]) for t in (15, 30, 60, 90, 120))
        P(f"     {nm:6s} | " + " ".join(f"{r:7.1f}" for r in rates) + f" | {ratio:10.2f} | {pw}")

    P("\n  c3  realised share ratios X1/G1 (each on its own shape) vs the offline a1/a2 (on G1's shape)")
    P("        t |  within 10 mm: X1     G1   X1/G1  a1 |  wall band: X1     G1   X1/G1  a2")
    off = offline_ratios()
    for t in (60, 90, 120):
        x, gg = shares(PB, "X1", t), shares(P11, "G1", t)
        o = off.get(t, (float("nan"), float("nan")))
        P(f"     {t:4d} |          {x[0]:.3f}  {gg[0]:.3f}  {x[0] / gg[0]:5.2f} {o[0]:5.2f} |          {x[1]:.3f}  "
          f"{gg[1]:.3f}  {x[1] / gg[1]:5.2f} {o[1]:5.2f}")

    P(f"\n  c4  X1 lateral firings (dir != +z) per 15 s window: h_cal = fl_h at the firing cell / S* ({S:.4f}); "
      f"Kant's minimum 500 W/m2K. (fl_h from the window's start plotfile, else its end; 'no h' = listed in neither)")
    P("     window  | firings  no h  below 500 | median h_cal  min h_cal [W/m2K]")
    for a, nl, miss, nb, med, mn in c4_rows("X1", S):
        P(f"     {a:3d}-{a + 15:3d} | {nl:7d} {miss:5d} {nb:10d} | {med:12.0f} {mn:10.0f}")
    P("    the FIELD bound (covers every firing): h_cal over every HEATED listed cell (fl_s >= 0, fl_h > 0) of each "
      "plotfile -- if its minimum is > 500, no face was below Kant's minimum at that time")
    g, PF = R11.pfs_by_t(PB, "X1")
    for tt in sorted(k for k in PF if k > 0):
        f = PB.read_pf(PB.legdir("X1"), PF[tt], g, comps=("fl_h", "fl_w", "fl_s"))
        hc = f["fl_h"][(f["fl_w"] > 0) & (f["fl_s"] >= 0) & (f["fl_h"] > 0)] / S
        P(f"     t {tt:3d} s: {hc.size} heated cells, h_cal min {hc.min():.0f}, p05 {np.percentile(hc, 5):.0f}, "
          f"median {np.median(hc):.0f}; share below 500 {(hc < 500).mean():.3f}, below 1000 {(hc < 1000).mean():.3f}")

    P("\n  c5  X1 to 120 s: bookkeeping, isotropy, cost")
    hdr, rows = thermo_rows(PB, "X1")
    mx = lambda c: max(abs(r[hdr.index(c)]) for r in rows if r[hdr.index("time")] > 0)
    mean = lambda h_, rr, c: float(np.mean([r[h_.index(c)] for r in rr if r[h_.index("time")] > 0]))
    s = PB.leg_summary("X1")
    secs, spr = s["sh"]["at"](30)
    hg, rg = thermo_rows(PB, "G1_C") if PB.done("G1_C") else (None, None)
    P(f"    ledger_err_solid max {mx('ledger_err_solid'):.2e} (<= 5e-3); fl_id_march max {mx('fl_id_march'):.2e}, "
      f"fl_id_ledger max {mx('fl_id_ledger'):.2e} (<= 1e-9); fr_multi max {mx('fr_multi'):g}; fl_n_unreached max "
      f"{mx('fl_n_unreached'):g}")
    P(f"    Ø by 30 deg sector at 30 mm depth: {R11.P9.fmt_sec(secs)}  spread {spr if spr is None else '%.2f' % spr} mm")
    med = lambda h_, rr, c: float(np.median([r[h_.index(c)] for r in rr if r[h_.index("time")] > 0]))
    # ⚠ medians: X1's run was suspended ~2 h mid-run (host sleep); the timed step that spanned it carries the
    # whole gap, which a mean would spread over every step
    P(f"    fl_t_path median {1e3 * med(hdr, rows, 'fl_t_path'):.2f} ms/step (max {1e3 * mx('fl_t_path'):.2f}); "
      f"fr_t_step median X1 {1e3 * med(hdr, rows, 'fr_t_step'):.1f} ms vs G1_C "
      + (f"{1e3 * med(hg, rg, 'fr_t_step'):.1f} ms" if hg else "n/a")
      + f"; fl_s_max at 120 s {1e3 * rows[-1][hdr.index('fl_s_max')]:.1f} mm, fl_H_stag {1e3 * rows[-1][hdr.index('fl_H_stag')]:.1f} mm")
    R30, R60 = A[30]["depth"] / G[30]["depth"], A[60]["depth"] / G[60]["depth"]
    P(f"\n  e1 inputs: matched 30 s volume (X1 {A[30]['V']:.2f}, G1 {G[30]['V']:.2f} cm3): depth X1/G1 at 30 / 60 / 90 / "
      f"120 s {R30:.3f} / {R60:.3f} / {A[90]['depth'] / G[90]['depth']:.3f} / {A[120]['depth'] / G[120]['depth']:.3f}; "
      f"max Ø X1 - G1 at 120 s {A[120]['dmax'] - G[120]['dmax']:+.1f} mm")


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "RESULTS"
    {"W": W, "IDENTITY": identity, "A6": a6, "SSTAR": sstar, "RESULTS": results}[mode]()


if __name__ == "__main__":
    main()
