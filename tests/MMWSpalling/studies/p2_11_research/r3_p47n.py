#!/usr/bin/env python3
"""P2-11 research R3 (zero compute): where does P2-9's Höser run P47n put its heat, and what would an
'exposure' rule (Goldstein's flat-plate law evaluated at the distance ALONG THE SURFACE from the stagnation point,
standoff to the stagnation point, the flange as a gas barrier) do to the same shapes?

Reads P47n plotfiles (binary B, quarter domain 50x50x150 at 1 mm, axis at x = y = 0, rock top z = 0.15 m,
nozzle z_n = 0.2132 m fixed). Conductance per surface cell G = h * w * dx^2 (w = fl_w, the face list's area
weight; h = today's map at the cell's top face, recomputed from the deck formula and checked against fl_h).
Gas temperature is NOT included (both rules share the march) -- shares are CONDUCTANCE shares.
Writes R3.txt beside this script."""
import os, sys, math, heapq, csv
import importlib.machinery, importlib.util
import numpy as np
from scipy.ndimage import gaussian_filter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
ld = importlib.machinery.SourceFileLoader("frtest", os.path.join(ROOT, "tests/MMWSpalling/unit/face_removal/test"))
m = importlib.util.module_from_spec(importlib.util.spec_from_loader("frtest", ld)); ld.exec_module(m)
LEG = os.path.join(ROOT, "tests/MMWSpalling/studies/p2_9_face_removal/legs/P47n")
SHAPE = (50, 50, 150); DX = 1e-3; ZTOP = 0.15; ZN = 0.2132; D = 0.0079
CSTAR = 1040.679839; S = 0.1431          # P2-9 d1: C* = S* x 7272.3958
HOESER = {15: (18, 40), 30: (31, 47), 60: (48, 56), 90: (61, 64), 120: (77, 73)}   # axis depth, max diam [mm]

def gold(rho):            # Goldstein 1986 radial factor, = 1 at the stagnation point
    return 533.0 / (533.0 + 44.0 * (rho / D) ** 1.285)
def gold_loc(rho):        # the LOCAL value implied by Goldstein's AREA-AVERAGE: (1/2r) d(r^2 Nu_avg)/dr, = 1 at r = 0
    u = 44.0 * (rho / D) ** 1.285
    return 533.0 * (533.0 + (1.0 - 1.285 / 2.0) * u) / (533.0 + u) ** 2
RNOZ = D / 2.0            # Rauenzahn Fig 6.20 (low p, Z/R = 10), Q/Q_stag vs surface distance / R_noz; beyond 15: e-fold 10 R
RZ = [(0, 1.0), (0.6, 1.0), (1.0, .90), (1.7, .75), (2.0, .62), (2.4, .62), (3.0, .53), (3.4, .51), (4.0, .45), (5.0, .36),
      (6.0, .34), (7.0, .30), (8.0, .28), (9.0, .25), (10, .24), (12, .20), (15, .14)]
def rau(rho):
    x = rho / RNOZ
    if x >= 15: return 0.14 * math.exp(-(x - 15) / 10.0)
    return float(np.interp(x, [a for a, b in RZ], [b for a, b in RZ]))
def rau_scaled(rho, H):   # Rauenzahn's cavity SHAPE with its half-width (3.4 R_noz on Fig 6.20) rescaled to 0.23 H
    return rau(rho * (3.4 * RNOZ) / (0.23 * H))
def stand(s):
    return (0.0632 / max(0.0632, s)) ** 1.3185318
def blocked_lid(R, s):
    return (R > 0.02) & (s < 0.0633)
def h_now(r, s):          # the deck's h_expr, token for token
    if r > 0.02 and s < 0.0633: return 0.001
    return CSTAR * gold(r) * stand(s)

def main():
    pfs = {round(m.pf_time(p)): p for p in m.plotfiles("../../" + os.path.relpath(LEG, os.path.join(ROOT, "tests/MMWSpalling")) ) } \
        if False else {round(m.pf_time(p)): p for p in sorted(__import__("glob").glob(os.path.join(LEG, "run", "[0-9]" * 5 + "cell")))}
    ev = list(csv.DictReader(open(os.path.join(LEG, "run_removal_events.csv"))))
    L = [__doc__.strip().split("\n")[0], "=" * 110]
    ii = (np.arange(SHAPE[0]) + 0.5) * DX
    X, Y, Z = np.meshgrid(ii, ii, (np.arange(SHAPE[2]) + 0.5) * DX, indexing="ij")
    R = np.hypot(X, Y)
    for t in (15, 30, 60, 90, 120):
        f = m.read_pf(pfs[t], ["fr_f", "fl_w", "fl_h", "T_surf_field"], SHAPE)
        F, W, FH = f["fr_f"], f["fl_w"], f["fl_h"]
        # ---- shape: void radius per depth
        void = F < 0.5
        prof = []
        for k in range(SHAPE[2] - 1, -1, -1):
            v = void[:, :, k]
            if not v.any(): break
            prof.append(((ZTOP - (k + 0.5) * DX) * 1e3, 2e3 * (R[:, :, k][v].max() + 0.5 * DX)))
        depth = (ZTOP - Z[void].min() + 0.5 * DX) * 1e3 if void.any() else 0.0
        # ---- smoothed normal (mirror at the symmetry planes x=0, y=0; void above the top)
        Fp = np.concatenate([F[::-1], F], 0); Fp = np.concatenate([Fp[:, ::-1], Fp], 1)
        Fp = np.concatenate([Fp, np.zeros((Fp.shape[0], Fp.shape[1], 3))], 2)
        Fs = gaussian_filter(Fp, 1.0, truncate=2.0)[SHAPE[0]:, SHAPE[1]:, :SHAPE[2]]
        gx, gy, gz = np.gradient(Fs)
        nrm = np.sqrt(gx**2 + gy**2 + gz**2) + 1e-30
        nz = -gz / nrm
        surf = W > 0
        # ---- today's h at each listed cell's top-face position, vs the plotfile's fl_h (max over faces)
        s_face = ZN - (Z + 0.5 * DX)
        Hn = np.vectorize(h_now)(R, s_face)
        G_now = np.where(surf, FH * W * DX * DX, 0.0)
        # ---- exposure: Dijkstra along surface cells (26-neighbour), from the stagnation cell
        idx = np.argwhere(surf)
        ax = idx[(idx[:, 0] == 0) & (idx[:, 1] == 0)]
        st = tuple(ax[np.argmin(ax[:, 2])])                 # deepest listed cell on the axis
        z_st = (st[2] + 1) * DX
        dist = np.full(SHAPE, np.inf); dist[st] = math.hypot(0.5 * DX, 0.5 * DX)
        # the flange: original-top cells at r > 20 mm are covered by steel -> not part of the gas path
        blocked = surf & (np.arange(SHAPE[2])[None, None, :] == SHAPE[2] - 1) & (R > 0.02)
        pq = [(dist[st], st)]
        nb = [(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1) if (a, b, c) != (0, 0, 0)]
        while pq:
            d0, (i, j, k) = heapq.heappop(pq)
            if d0 > dist[i, j, k]: continue
            for a, b, c in nb:
                p, q, r_ = i + a, j + b, k + c
                if p < 0 or q < 0 or r_ < 0 or p >= SHAPE[0] or q >= SHAPE[1] or r_ >= SHAPE[2]: continue
                if not surf[p, q, r_] or blocked[p, q, r_]: continue
                nd = d0 + DX * math.sqrt(a * a + b * b + c * c)
                if nd < dist[p, q, r_]:
                    dist[p, q, r_] = nd; heapq.heappush(pq, (nd, (p, q, r_)))
        s_st = ZN - z_st
        H_exp = np.where(surf & np.isfinite(dist), CSTAR * np.vectorize(gold)(np.where(np.isfinite(dist), dist, 1.0)) * stand(s_st), 0.0)
        G_exp = H_exp * W * DX * DX
        # same total conductance (the level S is re-matched in any real leg)
        G_exp *= G_now.sum() / max(G_exp.sum(), 1e-30)
        dd = np.where(np.isfinite(dist), dist, 1.0)
        def norm(G):
            return G * G_now.sum() / max(G.sum(), 1e-30)
        V = {"today (avg-as-local, r, own s)": G_now, "exposure (avg form, path, s_stag)": G_exp,
             "LOCAL Goldstein at r, own s": norm(np.where(surf, CSTAR * np.vectorize(gold_loc)(R) * np.vectorize(stand)(s_face) * W * DX * DX * (~blocked_lid(R, s_face)), 0.0)),
             "LOCAL Goldstein on path, s_stag": norm(np.where(surf & np.isfinite(dist), CSTAR * np.vectorize(gold_loc)(dd) * stand(s_st) * W * DX * DX, 0.0)),
             "cavity shape, half-width 0.23 H, path": norm(np.where(surf & np.isfinite(dist), CSTAR * np.vectorize(lambda q: rau_scaled(q, s_st))(dd) * stand(s_st) * W * DX * DX, 0.0)),
             "Rauenzahn cavity CFD on path, s_stag": norm(np.where(surf & np.isfinite(dist), CSTAR * np.vectorize(rau)(dd) * stand(s_st) * W * DX * DX, 0.0))}
        floor = surf & (nz > math.cos(math.radians(45)))
        depth_c = (ZTOP - Z) * 1e3
        L.append("t %3d s: model axis depth %.0f mm (Höser %d), max diam %.1f mm (Höser %d); stagnation cell %s, standoff to it %.1f mm"
                 % (t, depth, HOESER[t][0], max(p[1] for p in prof), HOESER[t][1], st, s_st * 1e3))
        L.append("   fl_h vs the recomputed map at the top-face position: median |rel| %.2e on cells with an exposed +z face"
                 % np.median(np.abs(FH[surf & floor] / np.maximum(Hn[surf & floor], 1e-30) - 1)))
        tot = G_now.sum()
        L.append("   conductance share  FLOOR (n_z > cos45): today %.3f  exposure %.3f   | WALL: today %.3f  exposure %.3f"
                 % (G_now[floor].sum() / tot, G_exp[floor].sum() / tot, G_now[surf & ~floor].sum() / tot, G_exp[surf & ~floor].sum() / tot))
        for nm, G in V.items():
            sel10 = surf & ~floor & (depth_c <= 10)
            L.append("   %-38s floor share %.3f | mouth-band (top 10 mm) wall share %.3f | share within r<10 mm %.3f"
                     % (nm, G[floor].sum() / tot, G[sel10].sum() / tot, G[surf & (R < 0.010)].sum() / tot))
        bands = [(0, 10), (10, 20), (20, 40), (40, 60), (60, 100)]
        row = []
        for a, b in bands:
            sel = surf & (depth_c > a) & (depth_c <= b)
            row.append("%d-%d mm %.3f->%.3f" % (a, b, G_now[sel].sum() / tot, G_exp[sel].sum() / tot))
        L.append("   share by depth below the original top (today->exposure): " + "; ".join(row))
        # ratio exposure/today per wall cell, at the mouth band
        sel = surf & ~floor & (depth_c <= 10) & (G_now > 0)
        if sel.any():
            L.append("   mouth-band wall cells (top 10 mm): exposure/today conductance median %.2f; path length median %.0f mm vs radius %.0f mm"
                     % (np.median(G_exp[sel] / G_now[sel]), np.median(dist[sel]) * 1e3, np.median(R[sel]) * 1e3))
        # ---- calorimeter-equivalent h at wall firings near t (t-7.5, t+7.5], vs Kant's h_min 500 W/m2K
        win = [e for e in ev if t - 7.5 < float(e["time"]) <= t + 7.5 and e["dir"] in ("0", "1", "2", "3", "4")]
        hc = []
        for e in win:
            i, j, k, d = int(e["col_i"]), int(e["col_j"]), int(e["k"]), int(e["dir"])
            x = (i + 0.5 + (0.5 if d == 1 else -0.5 if d == 0 else 0)) * DX
            y = (j + 0.5 + (0.5 if d == 3 else -0.5 if d == 2 else 0)) * DX
            zf = (k + 0.5) * DX if d < 4 else k * DX
            hc.append(h_now(math.hypot(x, y), ZN - zf) / S)
        hc = np.array(hc)
        if len(hc):
            L.append("   lateral firings in (t-7.5, t+7.5]: %d; calorimeter-equivalent h (= model h / S) median %.0f W/m2K, "
                     "share below Kant's h_min 500: %.3f, below 1000: %.3f" % (len(hc), np.median(hc), (hc < 500).mean(), (hc < 1000).mean()))
        L.append("   profile (depth mm: diam mm): " + ", ".join("%.0f:%.0f" % p for p in prof[::max(1, len(prof) // 10)]))
    txt = "\n".join(L); print(txt)
    open(os.path.join(HERE, "R3.txt"), "w").write(txt + "\n")

if __name__ == "__main__":
    main()
