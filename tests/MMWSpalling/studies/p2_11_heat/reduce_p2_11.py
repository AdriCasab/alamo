#!/usr/bin/env python3
"""P2-11 -- reductions of the G1 probe. Shape, volume and plotfile reads are P2-8's reducer
(`studies/p2_8_face_list/reduce_p2_8.py`), loaded once re-pointed at P2-9's legs (the control P47n) and
once at this study's, as P2-9 did -- so G1 and P47n are reduced by the SAME code.

  reduce_p2_11.py A1       a1: fl_h of the first plotfile after t = 0 vs the deck formula      -> A1.txt
  reduce_p2_11.py SSTAR    b1: S*_G1 by the packet's rule                                      -> SSTAR.txt
  reduce_p2_11.py          c1, c2, c3, e1 (stdout; the packet's command redirects it into RESULTS.txt)

AXIS DEPTH (packet, Sources): from `fr_f`, 150 mm minus the mean solid length (sum over k of f dz) of the
four innermost lines i, j in {0, 1}. MAX Ø: the mask's (P2-9 d3, P2-8 `shape`). V(f): sum of (1 - f) dx^3, x4.
"""
import importlib.util
import math
import os
import shlex
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


P9 = _load("r8_p29", os.path.join(STUD, "p2_9_face_removal"))     # P47n, binary B
P11 = _load("r8_p211", HERE)                                        # the G1 legs, binary C2
H0G = 7272.3958
TARGET_V = 33.0
D = 0.0079
SPREAD = {15: 0.18, 30: 0.50, 60: 0.65, 90: 1.41, 120: 0.65}       # P47n's own 4-line spread (packet)


# ------------------------------------------------------------------ helpers ----
def argv_kv(mod, n):
    kv = {}
    for t in shlex.split(open(os.path.join(mod.legdir(n), "argv.sh")).read()):
        if "=" in t:
            k, v = t.split("=", 1)
            kv[k] = v
    return kv


def h_of(mod, n):
    """the leg's h(r, s) from its own argv: the G1 (local) or the P47n (average) form, constant C read
    from the deck string. Evaluated in the parser's left-to-right order."""
    e = argv_kv(mod, n)["surface_patch.h_expr"]
    pre = "if((r>0.02)and(s<0.0633),0.001,"
    assert e.startswith(pre), e
    C = float(e[len(pre):].split("*", 1)[0])
    local = "0.3575" in e

    def h(r, s):
        if r > 0.02 and s < 0.0633:
            return 0.001
        stand = math.pow(0.0632 / max(0.0632, s), 1.3185305455)
        p = math.pow(r / 0.0079, 1.285)
        if local:
            return C * 533 * (533 + 0.3575 * 44 * p) / math.pow(533 + 44 * p, 2) * stand
        return C * 533 / (533 + 44 * p) * stand
    return h, C / H0G, local


def pfs_by_t(mod, n):
    d = mod.legdir(n)
    g = mod.geom_of(d)
    return g, {round(mod.t_of(pf, g)): pf for pf in mod.plotfiles(d)}


def fr_volume(mod, n, pf):
    d = mod.legdir(n)
    g = mod.geom_of(d)
    f = mod.read_pf(d, pf, g, comps=("fr_f",))["fr_f"]
    return (f.size - f.sum()) * g["dx"] * g["dx"] * g["dz"] * 4e6


def axis_depth(mod, n, pf):
    """[mm] 150 - mean solid length of lines (0,0), (0,1), (1,0), (1,1); and the spread (max - min)."""
    d = mod.legdir(n)
    g = mod.geom_of(d)
    f = mod.read_pf(d, pf, g, comps=("fr_f",))["fr_f"]
    Lz = (g["hi"][2] - g["lo"][2]) * 1e3
    deps = [Lz - f[i, j, :].sum() * g["dz"] * 1e3 for i in (0, 1) for j in (0, 1)]
    return float(np.mean(deps)), float(max(deps) - min(deps))


def faces(rem):
    """open faces per cell: in-domain void neighbours + the domain top (the face list's rule; P2-10 b4w).
    Returns {dir: bool array} for dir 0 -x .. 5 +z, restricted to solid cells."""
    V = np.pad(rem > 0.5, 1, mode="constant", constant_values=False)
    solid = ~(rem > 0.5)
    out = {}
    for dr in range(6):
        ax, sg = dr // 2, (1 if dr % 2 else -1)
        sl = [slice(1, -1)] * 3
        sl[ax] = slice(1 + sg, V.shape[ax] - 1 + sg)
        o = V[tuple(sl)].copy()
        if dr == 5:
            o[:, :, -1] = True
        out[dr] = o & solid
    return out


def face_pos(i, j, k, dr, g):
    dx, dz = g["dx"], g["dz"]
    lo = g["lo"]
    ox = 0.0 if dr == 0 else (1.0 if dr == 1 else 0.5)
    oy = 0.0 if dr == 2 else (1.0 if dr == 3 else 0.5)
    x = lo[0] + (i + ox) * dx
    y = lo[1] + (j + oy) * dx
    zf = lo[2] + (k + 1) * dz if dr == 5 else (lo[2] + k * dz if dr == 4 else lo[2] + (k + 0.5) * dz)
    return x, y, zf


def smoothed_nz(F):
    """the code's smoothed normal's z component (P2-10 b4w's offline stencil, b4w-verified): ghost rule
    x0 = y0 = 0 -> low x, y mirrored; above the top gas (0); other sides the nearest cell; normalised
    Gaussian sigma 1 cell, +-2; central differences. Returns n_z and the zero-gradient mask."""
    P = np.pad(F, ((3, 0), (0, 0), (0, 0)), mode="symmetric")
    P = np.pad(P, ((0, 3), (0, 0), (0, 0)), mode="edge")
    P = np.pad(P, ((0, 0), (3, 0), (0, 0)), mode="symmetric")
    P = np.pad(P, ((0, 0), (0, 3), (0, 0)), mode="edge")
    P = np.pad(P, ((0, 0), (0, 0), (3, 0)), mode="edge")
    P = np.pad(P, ((0, 0), (0, 0), (0, 3)), mode="constant", constant_values=0.0)
    e1, e2 = math.exp(-0.5), math.exp(-2.0)
    gs = 1.0 + 2.0 * e1 + 2.0 * e2
    w = {0: 1.0 / gs, 1: e1 / gs, 2: e2 / gs}
    S = P
    for ax in range(3):
        L = S.shape[ax] - 4
        acc = 0.0
        for a in range(-2, 3):
            sl = [slice(None)] * 3
            sl[ax] = slice(2 + a, 2 + a + L)
            acc = acc + w[abs(a)] * S[tuple(sl)]
        S = acc
    gx = (S[2:, 1:-1, 1:-1] - S[:-2, 1:-1, 1:-1]) / 2.0
    gy = (S[1:-1, 2:, 1:-1] - S[1:-1, :-2, 1:-1]) / 2.0
    gz = (S[1:-1, 1:-1, 2:] - S[1:-1, 1:-1, :-2]) / 2.0
    gn = np.sqrt(gx * gx + gy * gy + gz * gz)
    zero = ~(gn > 0)
    gn = np.where(zero, 1.0, gn)
    return -gz / gn, (np.abs(gx) + np.abs(gy) + np.abs(gz)) / gn, zero


# ------------------------------------------------------------------ a1 ----
def a1_leg(mod, n, L):
    g, P = pfs_by_t(mod, n)
    t = min(k for k in P if k > 0)
    pf = P[t]
    d = mod.legdir(n)
    f = mod.read_pf(d, pf, g, comps=("removed", "fl_h", "fl_w"))
    h, S, local = h_of(mod, n)
    zn = float(argv_kv(mod, n)["surface_patch.nozzle_z0"])
    fo = faces(f["removed"])
    listed = f["fl_w"] > 0
    rels = []
    nmiss = 0
    for (i, j, k) in np.argwhere(listed):
        best = None
        for dr in range(6):
            if fo[dr][i, j, k]:
                x, y, zf = face_pos(i, j, k, dr, g)
                v = h(math.sqrt(x * x + y * y), zn - zf)
                best = v if best is None else max(best, v)
        if best is None:
            nmiss += 1
            rels.append(float("inf"))
            continue
        fh = f["fl_h"][i, j, k]
        rels.append(abs(fh / best - 1) if best > 0 else (0.0 if fh == 0 else float("inf")))
    rels = np.array(rels)
    share = float((rels <= 1e-12).mean())
    ok = share >= 0.99
    L.append(f"  {n:7s} t {t:3d} s ({'G1 local' if local else 'P47n average'} form, S {S:.4f}): listed {listed.sum()}, "
             f"|rel| <= 1e-12 on {share:.4f} (median {np.median(rels):.2e}, max {np.max(rels):.2e}); "
             f"listed cells with no open face {nmiss}  {'PASS' if ok else 'FAIL'}")
    return ok


def a1():
    L = ["P2-11 (a1) -- the first plotfile after t = 0: plotfile fl_h (max over a cell's open faces) vs the deck "
         "formula at each open face's centre", "=" * 110]
    ok = True
    for n in ("G1_S14", "G1_S20", "G1_S10", "G1_S28", "G1"):
        if P11.done(n):
            ok &= a1_leg(P11, n, L)
    L.append("  control (the reducer itself, on binary B's P47n with ITS average form):")
    a1_leg(P9, "P47n", L)
    L.append(f"a1 {'PASS' if ok else 'FAIL'}")
    txt = "\n".join(L)
    print(txt)
    open(os.path.join(HERE, "A1.txt"), "w").write(txt + "\n")


# ------------------------------------------------------------------ b1 ----
def sstar():
    L = ["P2-11 (b1) -- S*_G1 by the rule: log V(f)(30 s, x4) vs log S between G1_S14 and G1_S20 to 33.0 cm3",
         "=" * 100]
    pts = []
    for n in ("G1_S10", "G1_S14", "G1_S20", "G1_S28"):
        if P11.done(n):
            g, P = pfs_by_t(P11, n)
            S = h_of(P11, n)[1]
            v = fr_volume(P11, n, P[30])
            pts.append((round(S, 4), v, n))
            L.append(f"  {n:7s} S {S:.4f}  V(f) 30 s {v:.3f} cm3")
    pts.sort()
    star = None
    for (s1, v1, n1), (s2, v2, n2) in zip(pts[:-1], pts[1:]):
        if (v1 - TARGET_V) * (v2 - TARGET_V) <= 0:
            x = math.log(s1) + (math.log(TARGET_V) - math.log(v1)) * \
                (math.log(s2) - math.log(s1)) / (math.log(v2) - math.log(v1))
            star = round(math.exp(x), 4)
            L.append(f"  bracket {n1}-{n2}: S* = {math.exp(x):.6f} -> {star:.4f}   (P1 predicted 0.150-0.170)")
            break
    if star is None:
        added = [n for n in ("G1_S10", "G1_S28") if P11.done(n)]
        if added:
            L.append(f"  NOT BRACKETED after the one added point ({added[0]}): STOP and report (b1 FAIL)")
        else:
            lo = max(p[1] for p in pts) if pts else None
            L.append("  NOT BRACKETED -- the rule adds ONE point: %s" %
                     ("G1_S28 (S 0.28)" if lo is not None and lo < TARGET_V else "G1_S10 (S 0.10)"))
        L.append("b1 not resolved")
    else:
        L.append(f"  C* = {star * H0G:.6f}")
        L.append(f"b1 PASS (interior, {len(pts) - 2} added point)")
        L.append(f"S_STAR={star:.4f}")
    txt = "\n".join(L)
    print(txt)
    open(os.path.join(HERE, "SSTAR.txt"), "w").write(txt + "\n")


# ------------------------------------------------------------------ c1-c3, e1 ----
def series(mod, n):
    g, P = pfs_by_t(mod, n)
    out = {}
    for t in sorted(P):
        if t <= 0:
            continue
        s = mod.leg_summary(n, mod.plotfiles(mod.legdir(n)).index(P[t]))
        dep, spr = axis_depth(mod, n, P[t])
        out[t] = dict(depth=dep, spread=spr, dmax=s["sh"]["dmax"], V=fr_volume(mod, n, P[t]),
                      mask_axis=s["sh"]["axis"] * 1e3)
    return out


def conductance(mod, n, t):
    g, P = pfs_by_t(mod, n)
    d = mod.legdir(n)
    f = mod.read_pf(d, P[t], g, comps=("fr_f", "fl_h", "fl_w", "removed"))
    nz, l1, zero = smoothed_nz(f["fr_f"])
    G = f["fl_h"] * f["fl_w"]
    listed = f["fl_w"] > 0
    x = (np.arange(g["nx"]) + 0.5) * g["dx"] + g["lo"][0]
    y = (np.arange(g["ny"]) + 0.5) * g["dx"] + g["lo"][1]
    R = np.hypot(x[:, None], y[None, :])[:, :, None] * np.ones((1, 1, g["nz"]))
    floor = listed & (nz > math.cos(math.radians(45))) & ~zero
    tot = G[listed].sum()
    fo = faces(f["removed"])
    nf = sum(fo[dr].astype(int) for dr in range(6))
    w_off = np.where(zero, 1.0, 1.0 / np.where(zero, 1.0, l1))
    agree = float((np.abs(f["fl_w"] - nf * w_off) <= 1e-9)[listed].mean())
    return dict(axis=G[listed & (R < 0.010)].sum() / tot, floor=G[floor].sum() / tot, agree=agree)


def c3_rows(mod, n, S):
    h, _, _ = h_of(mod, n)
    zn = float(argv_kv(mod, n)["surface_patch.nozzle_z0"])
    g = mod.geom_of(mod.legdir(n))
    ev = mod.events(mod.legdir(n))
    rows = []
    for a in range(0, int(round(g["st"])), 15):
        hc = []
        for e in ev:
            te = float(e["time"])
            if not (a < te <= a + 15) or e.get("dir") == "5":
                continue
            x, y, zf = face_pos(int(e["col_i"]), int(e["col_j"]), int(e["k"]), int(e["dir"]), g)
            hc.append(h(math.sqrt(x * x + y * y), zn - zf) / S)
        hc = np.array(hc)
        rows.append((a, len(hc), int((hc < 500).sum()) if len(hc) else 0,
                     float(np.median(hc)) if len(hc) else float("nan"),
                     float(hc.min()) if len(hc) else float("nan")))
    return rows


def results():
    P = print
    P("P2-11 -- the G1 probe: c1, c2, c3, e1 (binary C2 for G1; the control P47n is P2-9's, binary B)")
    P("=" * 110)
    if not P11.done("G1"):
        P("  G1 not done")
        return
    S = h_of(P11, "G1")[1]
    A = series(P11, "G1")
    B = series(P9, "P47n")
    HER = dict(zip([15, 30, 60, 90, 120], zip(P9.HER_DEPTH, P9.HER_DMAX)))
    P(f"\n  c1  G1 (S* {S:.4f}) vs P47n (S* 0.1431) vs Hoeser 47 MPa (Fig. 5.10). depth = fr_f axis depth "
      f"(4 innermost lines; spread in brackets); max Ø = mask")
    P("        t |  depth G1 (spr)   P47n (spr)  her   G1/P47n  G1/her |  maxØ G1  P47n  her   ΔØ |  V(f) G1   P47n   her")
    agree = True
    for t in (15, 30, 60, 90, 120):
        a, b = A[t], B[t]
        hd, hm = HER[t]
        if abs(a["depth"] / hd - 1) > 0.15 or (t >= 60 and abs(a["dmax"] - hm) > 6):
            agree = False
        P(f"     {t:4d} | {a['depth']:7.2f} ({a['spread']:.2f}) {b['depth']:7.2f} ({b['spread']:.2f}) {hd:4.0f}  "
          f"{a['depth'] / b['depth']:7.3f} {100 * (a['depth'] / hd - 1):+6.1f}% | {a['dmax']:6.1f} {b['dmax']:5.1f} "
          f"{hm:4.0f} {a['dmax'] - b['dmax']:+5.1f} | {a['V']:7.2f} {b['V']:7.2f} {P9.HER_V.get(float(t), float('nan')):5.1f}")
    for nm, X in (("G1", A), ("P47n", B)):
        gr = (X[120]["depth"] - X[15]["depth"]) / max(X[120]["dmax"] - X[15]["dmax"], 1e-9)
        P(f"    {nm:5s} V60/V30 {X[60]['V'] / X[30]['V']:.2f}, V90/V30 {X[90]['V'] / X[30]['V']:.2f} (hers 2.45, 4.09); "
          f"depth growth / diameter growth 15->120 s {gr:.2f} (hers 1.8); depth 15->120 {X[15]['depth']:.1f} -> "
          f"{X[120]['depth']:.1f} mm, max Ø {X[15]['dmax']:.1f} -> {X[120]['dmax']:.1f} mm")
    P(f"    G1 SHAPE AGREES with Hoeser (max Ø within 6 mm at 60/90/120 s AND depth within 15 % every time): "
      f"{'YES' if agree else 'NO'}")
    P("\n  c2  conductance (fl_h x fl_w) share within r < 10 mm of the axis, and on floor faces (smoothed n_z > cos 45)")
    P("        t |  axis G1   P47n  ratio |  floor G1   P47n  ratio | fl_w offline agreement G1 / P47n")
    for t in (15, 30, 60, 90, 120):
        a, b = conductance(P11, "G1", t), conductance(P9, "P47n", t)
        P(f"     {t:4d} | {a['axis']:7.3f} {b['axis']:6.3f} {a['axis'] / b['axis']:6.2f} | {a['floor']:8.3f} "
          f"{b['floor']:6.3f} {a['floor'] / b['floor']:6.2f} | {a['agree']:.4f} / {b['agree']:.4f}")
    P(f"\n  c3  G1 lateral firings (dir != +z) per 15 s window: calorimeter-equivalent h = model h at the firing face / S* "
      f"({S:.4f}); Kant's minimum 500 W/m2K")
    P("     window  | firings  below 500 | median h_cal  min h_cal [W/m2K]")
    for a, nl, nb, med, mn in c3_rows(P11, "G1", S):
        P(f"     {a:3d}-{a + 15:3d} | {nl:7d} {nb:9d} | {med:12.0f} {mn:10.0f}")
    R30 = A[30]["depth"] / B[30]["depth"]
    R60 = A[60]["depth"] / B[60]["depth"]
    dO = A[120]["dmax"] - B[120]["dmax"]
    if R30 < 0.97 or R60 < 0.97 or dO > 3:
        v = "STOP (wrong direction)"
    elif R30 > 1.25 or R60 > 1.25:
        v = "OVERSHOOT => GO, X1's predictions re-derived offline first"
    elif 1.04 <= R30 <= 1.25 and 1.04 <= R60 <= 1.25 and dO <= -3:
        v = "GO"
    elif R30 < 1.04 and R60 < 1.04 and dO > -3:
        v = "STOP (insensitive)"
    else:
        v = "none of the bars -- report; the user decides"
    P(f"\n  e1  matched 30 s volume (G1 {A[30]['V']:.2f}, P47n {B[30]['V']:.2f} cm3): R30 = {R30:.3f}, R60 = {R60:.3f}; "
      f"ΔØ(120 s) = {dO:+.1f} mm")
    P(f"      GO: R in [1.04, 1.25] at both AND ΔØ <= -3; STOP insensitive: R < 1.04 at both AND ΔØ > -3; "
      f"STOP wrong direction: R < 0.97 at either OR ΔØ > +3; R > 1.25 overshoot")
    P(f"  e1  VERDICT: {v}")


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "RESULTS"
    {"A1": a1, "SSTAR": sstar, "RESULTS": results}[mode]()


if __name__ == "__main__":
    main()
