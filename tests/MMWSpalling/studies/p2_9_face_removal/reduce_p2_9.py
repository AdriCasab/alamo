#!/usr/bin/env python3
"""P2-9 rev 2 -- reductions of the Hoeser legs. Every quantity is computed by
P2-8's reducer (`studies/p2_8_face_list/reduce_p2_8.py`, loaded twice: once on P2-8's
legs, once re-pointed at this study's), so a P2-8 / P2-9 comparison is the SAME code on
both sides; P2-9's own columns (fr_*) are read here.

  reduce_p2_9.py GATE    b0 (K0g vs P2-8 K0), b1 (M0 witnesses)       -> GATE_b.txt
  reduce_p2_9.py MESH    b2, c1 (M0 vs L0), c2, c3 (M0c, M0f vs M0)    -> MESH.txt
  reduce_p2_9.py SSTAR   d1 (the packet's rule)                        -> SSTAR.txt
  reduce_p2_9.py         d2, d3, d4, the e1 table, the P2-10 deletion list (stdout;
                         the packet's command redirects it into RESULTS.txt)
VOLUME. P2-8's V counts voided cells from the event log. Under the line rule a cell
also recedes PARTLY (solid fraction f), so `V(f)` = sum over cells of (1 - f) dx^3 x 4
is reported beside it, and d1 uses V(f) -- the volume that left. P2-8 legs have no
f; their V is the void count.
"""
import importlib.util
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)
P8DIR = os.path.join(STUD, "p2_8_face_list")


def _load(name, legroot):
    spec = importlib.util.spec_from_file_location(name, os.path.join(P8DIR, "reduce_p2_8.py"))
    m = importlib.util.module_from_spec(spec)
    sys.path.insert(0, P8DIR)
    spec.loader.exec_module(m)
    m.legdir = lambda n, _r=legroot: os.path.join(_r, "legs", n)
    return m


P8 = _load("r8_p28", P8DIR)
P9 = _load("r8_p29", HERE)
H0G = 7272.3958
TARGET_V = 33.0


def fr_volume(n, pf=None):
    """V(f) [cm3, x4]: sum of (1 - f) dx^3 over the initially solid cells."""
    d = P9.legdir(n)
    g = P9.geom_of(d)
    pfs = P9.plotfiles(d)
    pf = pf or pfs[-1]
    f = P9.read_pf(d, pf, g, comps=("fr_f",))["fr_f"]
    vc = g["dx"] * g["dx"] * g["dz"]
    return (f.size - f.sum()) * vc * 4e6


def rows_of(n, mod=P9):
    hdr, rows = mod.thermo(mod.legdir(n))
    it = hdr.index("time")
    return hdr, [r for r in rows if float(r[it]) > 0]


def colmax(hdr, rows, c, absval=False):
    if c not in hdr:
        return None
    i = hdr.index(c)
    v = [float(r[i]) for r in rows]
    return max(abs(x) for x in v) if absval else max(v)


def colsum(hdr, rows, c):
    i = hdr.index(c)
    return sum(float(r[i]) for r in rows)


def gate():
    out = ["P2-9 rev 2 (b0, b1) -- the edge-cell fix on the key-off leg, and M0's witnesses",
           "=" * 96]
    P = out.append
    if P9.done("K0g"):
        a, b = P8.leg_summary("K0"), P9.leg_summary("K0g")
        P(f"  b0  K0 (P2-8, f4c27b27) vs K0g (binary B, key off): V {a['V']:.3f} -> {b['V']:.3f} cm3 "
          f"({100 * (b['V'] / a['V'] - 1):+.2f} %); axis {a['sh']['axis'] * 1e3:.1f} -> "
          f"{b['sh']['axis'] * 1e3:.1f} mm; max Ø {a['sh']['dmax']:.1f} -> {b['sh']['dmax']:.1f} mm; "
          f"ledger_err_solid max {a['led']:.3e} -> {b['led']:.3e}")
        P(f"      (P0b predicted < 2 % in volume and < 1 mm in axis depth)")
    else:
        P("  b0  K0g not run")
    if P9.done("M0"):
        s = P9.leg_summary("M0")
        hdr, rows = rows_of("M0")
        ev = P9.events(P9.legdir("M0"))
        g = s["g"]
        n6 = sum(1 for e in ev if int(float(e["regime"])) == 6)
        nwall = sum(1 for e in ev if e.get("dir") in ("0", "1", "2", "3"))
        nroof = sum(1 for e in ev if e.get("dir") == "4")
        lid = 0
        for e in ev:
            if e.get("dir") != "5" or int(e["k"]) != g["nz"] - 1:
                continue
            x = (int(e["col_i"]) + 0.5) * g["dx"]
            y = (int(e["col_j"]) + 0.5) * g["dx"]
            if math.hypot(x, y) > 0.02:
                lid += 1
        mm = colmax(hdr, rows, "fr_multi")
        idr = colmax(hdr, rows, "fr_id_remap")
        ok = (s["led"] <= 5e-3 and s["idm"] <= 1e-9 and s["idl"] <= 1e-9 and mm == 0
              and idr <= 1e-12 and n6 == 0 and lid == 0)
        P(f"  b1  M0: ledger_err_solid max {s['led']:.3e} (<= 5e-3); fl_id_march {s['idm']:.2e}, "
          f"fl_id_ledger {s['idl']:.2e} (<= 1e-9); fr_multi max {mm:g} (0); fr_id_remap max {idr:.2e} "
          f"(<= 1e-12); regime-6 rows {n6} (0); firings on the lidded original surface "
          f"(+z line, top cell, r > 20 mm) {lid} (0)")
        P(f"      witnesses the line path runs: wall firings {nwall}, roof firings {nroof}, "
          f"floor {sum(1 for e in ev if e.get('dir') == '5')}")
        P(f"  b1  {'PASS' if ok else 'FAIL'}")
    else:
        P("  b1  M0 not run")
    txt = "\n".join(out)
    print(txt)
    open(os.path.join(HERE, "GATE_b.txt"), "w").write(txt + "\n")


def classes(ev):
    c = {"floor": [0, 0], "wall": [0, 0], "roof": [0, 0]}
    for e in ev:
        k = {"5": "floor", "4": "roof"}.get(e.get("dir"), "wall")
        c[k][0] += 1
        c[k][1] += int(float(e["n_voided"]))
    return c


def share(hdr, rows, num, den):
    a = sum(colsum(hdr, rows, x) for x in num)
    b = sum(colsum(hdr, rows, x) for x in den)
    return a / b if b else float("nan")


def mesh():
    out = ["P2-9 rev 2 (b2, c1-c3) -- the mesh legs, M0 vs P2-8's L0, isotropy", "=" * 96]
    P = out.append
    for n in ("M0c", "M0f"):
        P(f"  b2  {n} reaches 30 s: {'rc 0' if P9.done(n) else 'NOT DONE'}")
    if P9.done("M0"):
        a, b = P8.leg_summary("L0"), P9.leg_summary("M0")
        hdr, rows = rows_of("M0")
        h8, r8 = rows_of("L0", P8)
        ev = P9.events(P9.legdir("M0"))
        cl = classes(ev)
        g = b["g"]
        rem = P9.read_pf(P9.legdir("M0"), b["pf"], g)["removed"]
        pc = P8.plate_census(rem, g)
        last = rows[-1]
        P(f"  c1  V (voids) L0 {a['V']:.3f} -> M0 {b['V']:.3f} cm3, ratio {b['V'] / a['V']:.3f}; "
          f"V(f) M0 {fr_volume('M0'):.3f} cm3")
        P(f"      mouth Ø {a['sh']['mouth']:.1f} -> {b['sh']['mouth']:.1f} mm; max Ø {a['sh']['dmax']:.1f} -> "
          f"{b['sh']['dmax']:.1f} mm; axis {a['sh']['axis'] * 1e3:.1f} -> {b['sh']['axis'] * 1e3:.1f} mm; "
          f"r17.5/axis {a['sh']['r175']:.3f} -> {b['sh']['r175']:.3f}")
        P(f"      firings / voids by line class: " + ", ".join(
            f"{k} {v[0]} / {v[1]}" for k, v in cl.items()))
        P(f"      sheets at 30 s (the code: last thermo row) fr_sheet {float(last[hdr.index('fr_sheet')]):g}, "
          f"fr_misaligned {float(last[hdr.index('fr_misaligned')]):g} of fr_ncell "
          f"{float(last[hdr.index('fr_ncell')]):g}; P2-8's PLATE rule on M0's mask: "
          f"{pc['pl_cells']} cells / {pc['pl_faces']} faces ({100 * pc['pl_faces'] / max(pc['faces'], 1):.2f} %)"
          f"  (L0: 116 / 276)")
        ps = float(last[hdr.index("fl_P_side")])
        ps8 = float(r8[-1][h8.index("fl_P_side")])
        P(f"      fl_P_side M0 / L0 at 30 s {ps:.4g} / {ps8:.4g} W = {ps / ps8:.3f}")
        nf = ["fr_n_floor", "fr_n_wall", "fr_n_roof"]
        P(f"      shares over the sampled thermo rows (ratios, SR 7): fr_clip / firings "
          f"{share(hdr, rows, ['fr_clip'], nf):.4f}; fr_cap / firings {share(hdr, rows, ['fr_cap'], nf):.4f}; "
          f"fallback T_s / surface cells {share(hdr, rows, ['fr_fb_floor', 'fr_fb_wall', 'fr_fb_roof'], ['fr_ncell']):.4f} "
          f"(floor {share(hdr, rows, ['fr_fb_floor'], ['fr_ncell']):.4f}, wall "
          f"{share(hdr, rows, ['fr_fb_wall'], ['fr_ncell']):.4f}, roof {share(hdr, rows, ['fr_fb_roof'], ['fr_ncell']):.4f}); "
          f"noskin / surface cells {share(hdr, rows, ['fr_noskin_floor', 'fr_noskin_wall', 'fr_noskin_roof'], ['fr_ncell']):.4f}; "
          f"misaligned / surface cells {share(hdr, rows, ['fr_misaligned'], ['fr_ncell']):.4f}; "
          f"overlap / firings {share(hdr, rows, ['fr_overlap'], nf):.4f}")
        tint = b["g"]["tint"]
        Ed = colsum(hdr, rows, "fr_E_drop") * tint
        Ef = float(last[hdr.index("ledger_E_flake")])
        P(f"      fr_E_drop (sampled sum x thermo.plot_int {tint}, an ESTIMATE) / ledger_E_flake = {Ed / Ef:.2e}")
        ts = [float(r[hdr.index("fr_t_step")]) for r in rows]
        tr = [float(r[hdr.index("fr_t_removal")]) for r in rows]
        w8 = eval(open(os.path.join(P8.legdir("L0"), ".done")).read())["wall"]
        steps = int(round(g["st"] / g["dt"]))
        P(f"      step wall time: M0 mean {np.mean(ts):.4f} s (removal {np.mean(tr):.4f} s); L0 "
          f"{w8 / steps:.4f} s (its run wall / steps) -> ratio {np.mean(ts) / (w8 / steps):.2f}")
        P(f"  c2  M0 Ø by 30 deg sector at 5 mm:  {P8.fmt_sec(b['sh']['at'](5)[0])}  spread "
          f"{b['sh']['at'](5)[1]:.2f} mm (<= 3)")
        P(f"  c2  M0 Ø by 30 deg sector at 10 mm: {P8.fmt_sec(b['sh']['at'](10)[0])}  spread "
          f"{b['sh']['at'](10)[1]:.2f} mm (<= 3)")
        for n, vb, db, ab in (("M0c", 15, 4, None), ("M0f", 5, 2, 1.5)):
            if not P9.done(n):
                P(f"  c3  {n} not done")
                continue
            c = P9.leg_summary(n)
            dv = 100 * abs(fr_volume(n) / fr_volume("M0") - 1)
            dd = abs(c["sh"]["dmax"] - b["sh"]["dmax"])
            da = abs(c["sh"]["axis"] - b["sh"]["axis"]) * 1e3
            okv = dv <= vb and dd <= db and (ab is None or da <= ab)
            P(f"  c3  {n} vs M0: V(f) {fr_volume(n):.3f} vs {fr_volume('M0'):.3f} cm3 ({dv:.1f} %, <= {vb}); "
              f"max Ø {c['sh']['dmax']:.1f} vs {b['sh']['dmax']:.1f} ({dd:.1f} mm, <= {db}); axis "
              f"{c['sh']['axis'] * 1e3:.1f} vs {b['sh']['axis'] * 1e3:.1f} ({da:.1f} mm"
              f"{', <= %.1f' % ab if ab else ''})  {'within' if okv else 'OUTSIDE'}")
    txt = "\n".join(out)
    print(txt)
    open(os.path.join(HERE, "MESH.txt"), "w").write(txt + "\n")


def sstar():
    out = ["P2-9 rev 2 (d1) -- S*new by the rule: log V(f)(30 s, x4) vs log S between M0 (0.10) and "
           "C15n (0.15) to 33.0 cm3", "=" * 96]
    P = out.append
    pts = []
    for n, S in (("M0", 0.10), ("C15n", 0.15), ("C07n", 0.07), ("C20n", 0.20)):
        if P9.done(n):
            v = fr_volume(n)
            pts.append((S, v, n))
            P(f"  {n:5s} S {S:.2f}  V(f) {v:.3f} cm3 (voids {P9.leg_summary(n)['V']:.3f})")
    pts.sort()
    star = None
    for (s1, v1, n1), (s2, v2, n2) in zip(pts[:-1], pts[1:]):
        if (v1 - TARGET_V) * (v2 - TARGET_V) <= 0:
            x = math.log(s1) + (math.log(TARGET_V) - math.log(v1)) * \
                (math.log(s2) - math.log(s1)) / (math.log(v2) - math.log(v1))
            star = round(math.exp(x), 4)
            P(f"  bracket {n1}-{n2}: S* = {math.exp(x):.6f} -> {star:.4f}")
            break
    if star is None:
        lo = min(p[1] for p in pts) if pts else None
        P("  NOT BRACKETED -- the rule adds ONE point: %s" %
          ("C20n (S 0.20)" if lo is not None and lo < TARGET_V else "C07n (S 0.07)"))
    else:
        P(f"  C* = {star * H0G:.6f}; |S* - 0.1314| = {abs(star - 0.1314):.4f} -> "
          f"{'P47n RUNS' if abs(star - 0.1314) > 0.0039 else 'P47n skipped (P47r stands in)'}"
          f"   (P3 predicted 0.10-0.14)")
        P(f"S_STAR={star:.4f}")
    txt = "\n".join(out)
    print(txt)
    open(os.path.join(HERE, "SSTAR.txt"), "w").write(txt + "\n")


def series(mod, n):
    d = mod.legdir(n)
    g = mod.geom_of(d)
    rows = []
    for pf in mod.plotfiles(d)[1:]:
        s = mod.leg_summary(n, mod.plotfiles(d).index(pf))
        vf = fr_volume(n, pf) if mod is P9 else None
        rows.append((s["t"], s, vf))
    return g, rows


def results():
    P = print
    P("P2-9 rev 2 -- the line rule on Hoeser: d2, d3, d4, e1")
    P("=" * 110)
    if P9.done("P47r") and P8.done("P47"):
        P("\n  d2  P47r (line rule) vs P2-8 P47 (same C* 955.592808)")
        P("        t |   V P2-8  V P47r(voids)  V(f) |  axis 8  axis 9 | maxØ 8  maxØ 9 | mouth 8 mouth 9 | sheets/faces 9")
        _, r8 = series(P8, "P47")
        _, r9 = series(P9, "P47r")
        h, rows = rows_of("P47r")
        it = h.index("time")
        for (t8, a, _), (t9, b, vf) in zip(r8, r9):
            near = min(rows, key=lambda r: abs(float(r[it]) - t9))
            shr = float(near[h.index("fr_sheet")]) / max(float(near[h.index("fr_ncell")]), 1)
            P(f"     {t9:5.0f} | {a['V']:8.2f} {b['V']:12.2f} {vf:6.2f} | {a['sh']['axis'] * 1e3:6.1f} "
              f"{b['sh']['axis'] * 1e3:7.1f} | {a['sh']['dmax']:6.1f} {b['sh']['dmax']:7.1f} | "
              f"{a['sh']['mouth']:7.1f} {b['sh']['mouth']:7.1f} | {shr:.4f}")
    leg = "P47n" if P9.done("P47n") else "P47r"
    if P9.done(leg):
        P(f"\n  d3  {leg} vs Hoeser 47 MPa (Fig. 5.10)")
        P("        t |  V(f)  V her |  axis  her     % |  maxØ  her     Δ | mouth")
        _, r9 = series(P9, leg)
        byt = {round(t): (s, vf) for t, s, vf in r9}
        Vt = {}
        agree = True
        for t, hd, hm in zip(P8.HER_T, P8.HER_DEPTH, P8.HER_DMAX):
            if round(t) not in byt:
                continue
            s, vf = byt[round(t)]
            Vt[t] = vf
            ax = s["sh"]["axis"] * 1e3
            dpc = 100 * (ax / hd - 1)
            dd = s["sh"]["dmax"] - hm
            if abs(dpc) > 15 or (t >= 60 and abs(dd) > 6):
                agree = False
            P(f"     {t:5.0f} | {vf:6.2f} {P8.HER_V.get(t, float('nan')):6.1f} | {ax:5.1f} {hd:4.0f} {dpc:+6.1f} | "
              f"{s['sh']['dmax']:5.1f} {hm:4.0f} {dd:+6.1f} | {s['sh']['mouth']:5.1f}")
        if 30.0 in Vt and 60.0 in Vt and 90.0 in Vt:
            P(f"    V60/V30 {Vt[60.0] / Vt[30.0]:.2f} (hers 2.45); V90/V30 {Vt[90.0] / Vt[30.0]:.2f} (hers 4.09)")
        if 15.0 in byt and 120.0 in byt:
            s15, s120 = byt[15][0], byt[120][0]
            gr = (s120["sh"]["axis"] - s15["sh"]["axis"]) * 1e3 / max(s120["sh"]["dmax"] - s15["sh"]["dmax"], 1e-9)
            P(f"    depth growth / diameter growth 15->120 s: {gr:.2f} (hers 1.8)")
        P(f"    SHAPE AGREES (max Ø within 6 mm at 60/90/120 s AND axis within 15 % every time): "
          f"{'YES' if agree else 'NO'}")
    for n in ("P47r", "P47n"):
        if not P9.done(n):
            continue
        s = P9.leg_summary(n)
        h, rows = rows_of(n)
        ev = P9.events(P9.legdir(n))
        g = s["g"]
        first50 = None
        for e in ev:
            x = (int(e["col_i"]) + 0.5) * g["dx"]
            y = (int(e["col_j"]) + 0.5) * g["dx"]
            if math.hypot(x, y) >= 0.05 - 0.5 * g["dx"]:
                first50 = float(e["time"])
                break
        P(f"\n  d4  {n} at 120 s: ledger_err_solid max {s['led']:.3e} (5e-3); identities march {s['idm']:.2e} "
          f"ledger {s['idl']:.2e} (1e-9); fr_multi max {colmax(h, rows, 'fr_multi'):g}; "
          f"first firing at r >= 50 mm: {first50 if first50 is not None else 'none'} s")
        for dmm in (30, 50):
            secs, spr = s["sh"]["at"](dmm)
            if secs:
                P(f"      Ø by sector at {dmm} mm: {P8.fmt_sec(secs)}  spread "
                  f"{spr if spr is None else '%.2f' % spr} (bar 3)")
    P("\n  P2-10 DELETION LIST (deliverable 9; rev 1 takeaway 7, from the code as built) -- under the key")
    P("  these are unreachable:")
    for x in ("Removal.H's lateral criterion (the `lat_on` block) and its transitional skip `if (k == kT) continue;`",
              "regime 6, `lat_accum_col` and the lateral thermal shift block (P1g-1/P1g-2)",
              "the four spall.lateral_* keys",
              "the per-column vertical scan, winner loop and PASS 3 remap (they survive as the column path's",
              "  callers of the shared ScanTAtZv / ScanKIcPa / ScanSigmaT / DepthScanH / RemapProfile)",
              "BuildSurfaceTColumns and surf_x_hold (diagnostics only under the key) and SurfaceRobinT's physics readers",
              "FaceListFlux (F4's Newton BC branch: 0 calls)",
              "group B (the nine side_face_* settings) and the old side-face blocks of kernel, ledger, march",
              "P2-8's +-1 indicator gradient for F3 (replaced by FrSmoothedGradient)"):
        P("    - " + x)


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "RESULTS"
    {"GATE": gate, "MESH": mesh, "SSTAR": sstar, "RESULTS": results}[mode]()


if __name__ == "__main__":
    main()
