#!/usr/bin/env python3
"""P2-32 planning, zero compute: the alternative -- the flake is the depth where the TEMPERATURE RISE has fallen to
e^{-c_f} (0.6703) of the surface rise (= c_f ell on an exponential skin at every surface temperature, no property
table in the rule). Same 1-D columns as stress_depth.py / kant_bracket.py; L_req = C_L(T_s) x; bk = L_req / 12 mm.
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stress_depth as S

F = math.exp(-0.4)


def depth_rise(z, T, Ts, T0, f=F):
    r = np.concatenate([[Ts], T]) - T0
    zz = np.concatenate([[0.0], z])
    m = r <= f * r[0]
    if not m.any():
        return float("nan")
    i = int(np.argmax(m))
    return zz[i - 1] + (zz[i] - zz[i - 1]) * (r[i - 1] - f * r[0]) / (r[i - 1] - r[i])


def kant(leg, t_max=200.0, dz=0.1e-3, depth=0.064):
    a, rho = {"G": (1923.5, 2594.0), "B": (1540.6, 2566.0)}[leg[1]]
    h, th = {130: (747, 1250), 140: (690, 1175), 150: (633, 1100), 160: (576, 1025)}[int(leg[2:])]
    Tg, T0 = th + 273.15, 293.15
    lam = lambda T: a / (T - 273.15 + 563.0)
    n = int(round(depth / dz))
    z = (np.arange(n) + 0.5) * dz
    T = np.full(n, T0)
    dt = 0.4 * dz * dz / (lam(T0) / (rho * S.cp(T0)))
    t, nxt, rows = 0.0, 820.0, []
    while t < t_max:
        k = lam(T[0])
        Ts = T[0]
        for _ in range(30):
            g = h * (Tg - Ts) - 10.0 * (Ts - T0) - 0.8 * S.SIG * (Ts ** 4 - T0 ** 4) - k * (Ts - T[0]) / (0.5 * dz)
            Ts -= g / (-h - 10.0 - 4 * 0.8 * S.SIG * Ts ** 3 - k / (0.5 * dz))
        if Ts >= nxt:
            q = h * (Tg - Ts) - 10.0 * (Ts - T0) - 0.8 * S.SIG * (Ts ** 4 - T0 ** 4)
            rows.append((t, Ts, S.C_L(Ts) * 0.4 * lam(Ts) * (Ts - T0) / q, S.C_L(Ts) * depth_rise(z, T, Ts, T0)))
            nxt += 10.0
        L = lam(T)
        Lf = 2.0 * L[:-1] * L[1:] / (L[:-1] + L[1:])
        Fx = Lf * (T[1:] - T[:-1]) / dz
        div = np.zeros(n)
        div[:-1] += Fx
        div[1:] -= Fx
        div[0] += k * (Ts - T[0]) / (0.5 * dz)
        T = T + dt * (div / dz) / (rho * S.cp(T))
        t += dt
    return rows


p = print
p("KANT'S BRACKET under the rise-fraction depth: L_req at the first firing (820 K) and its minimum over 200 s")
res = {}
LOG = {"KG130": 9.915, "KG140": 12.010, "KG150": 14.956, "KG160": 19.295, "KB150": 11.992, "KB160": 15.430}
for leg in ("KG130", "KG140", "KG150", "KG160", "KB150", "KB160"):
    rows = kant(leg)
    i = int(np.argmin([r[3] for r in rows]))
    res[leg] = (rows[0][3], rows[i][3])
    p("   %-6s today %7.3f mm (logged %6.3f) | rise-fraction at onset %7.3f mm (x %.3f of today's) | minimum %7.3f mm at T_s %.0f K; at 900 / 1000 K: %s" % (
        leg, 1e3 * rows[0][2], LOG[leg], 1e3 * rows[0][3], rows[0][3] / rows[0][2], 1e3 * rows[i][3], rows[i][1],
        " / ".join("%.1f" % (1e3 * r[3]) for r in rows if abs(r[1] - 900) < 5 or abs(r[1] - 1000) < 5)))
for lab, k in (("at the first firing", 0), ("minimum over the history", 1)):
    g = (res["KG130"][k], min(res["KG140"][k], res["KG150"][k], res["KG160"][k]))
    b = (res["KB150"][k], res["KB160"][k])
    lo, hi = max(g[0], b[0]), min(g[1], b[1])
    p("   %s: Gotthard [%.3f, %.3f) mm, Bethel [%.3f, %.3f) mm -> one span [%.3f, %.3f): %s" % (
        lab, 1e3 * g[0], 1e3 * g[1], 1e3 * b[0], 1e3 * b[1], 1e3 * lo, 1e3 * hi,
        "EXISTS, width %.2f %%" % (100 * (hi / lo - 1)) if lo < hi else "EMPTY by %.2f %%" % (100 * (lo / hi - 1))))

p("\nTHE FLAT FLAME FIXTURE at its first firing (846 K): the rise-fraction depth, exact and read from cell samples")
T0 = 293.15
t, Ts, z, T, _ = S.column(T0, 20.0, 846.0, h=700.0, Tg=1866.0, hc=10.0, Tamb=293.15, snaps=(846.0,))[0]
ex = depth_rise(z, T, Ts, T0)
p("   exact %.4f mm (0.4 x the tangent length = %.4f mm)" % (1e3 * ex, 0.4e3 * S.lam(Ts) * (Ts - T0) / (700.0 * (1866.0 - Ts) - 10.0 * (Ts - T0) - 0.8 * S.SIG * (Ts ** 4 - T0 ** 4))))
for dzc in (2e-3, 1e-3, 0.5e-3, 0.25e-3):
    zc, Tc = S.samples(z, T, dzc, 15)
    xx = math.log((Tc[0] - T0) / (Tc[1] - T0))
    two = dzc / xx
    walk = S.fall_depth(np.concatenate([[0.0], zc]), np.concatenate([[Ts - T0], Tc - T0]), F * (Ts - T0))
    # the hybrid: the reconstruction's exponential inside the top cell (heat enters), the samples below it
    if 0.4 * two <= dzc:
        hyb = 0.4 * two
    else:
        zs = np.concatenate([[dzc], zc[1:]]); rs = np.concatenate([[(Ts - T0) * math.exp(-xx)], Tc[1:] - T0])
        hyb = S.fall_depth(zs, rs, F * (Ts - T0))
    p("   dz %5.3f mm: 0.4 x two-cell %.4f mm (%+.1f %%) | sample walk from the face %.4f mm (%+.1f %%) | exponential in the top cell, samples below %.4f mm (%+.1f %%)" % (
        1e3 * dzc, 0.4e3 * two, 100 * (0.4 * two / ex - 1), 1e3 * walk, 100 * (walk / ex - 1), 1e3 * hyb, 100 * (hyb / ex - 1)))

p("\nTHE BEAM COLUMN under the rise-fraction depth: bk from the criterion (846 K) to the solidus")
T0b = 297.5
for d, qq in ((0.5, 2.0), (1.0, 2.0), (2.0, 2.0), (4.0, 2.0), (8.0, 2.0), (2.0, 4.0)):
    snaps = (846.0, 900.0, 1000.0, 1100.0, 1200.0, 1273.0, 1350.0, 1410.0)
    res = S.column(T0b, 60.0, 1410.5, q=qq * 1e6, alpha=1.0 / (d * 1e-3), hc=56.0, Tamb=297.5, snaps=snaps)
    line, first = [], None
    for t, Ts, z, T, lv in res:
        x = depth_rise(z, T, Ts, T0b)
        bk = S.C_L(Ts, T0b) * x / 0.012
        if first is None and bk <= 1.0:
            first = (lv, t)
        line.append("%.0f K: x %.2f mm bk %.2f" % (lv, 1e3 * x, bk))
    t, Ts, z, T, _ = res[0]
    ex = depth_rise(z, T, Ts, T0b)
    row = []
    for n in (16, 32, 64, 128, 256):
        dzc = 0.05 / n
        zc, Tc = S.samples(z, T, dzc, int(0.015 / dzc) + 1)
        k = S.lam(Tc[0]); Tsf = Tc[0]
        for _ in range(30):
            g = -56.0 * (Tsf - 297.5) - 0.8 * S.SIG * (Tsf ** 4 - 297.5 ** 4) - k * (Tsf - Tc[0]) / (0.5 * dzc)
            Tsf -= g / (-56.0 - 4 * 0.8 * S.SIG * Tsf ** 3 - k / (0.5 * dzc))
        row.append("%.3f" % (1e3 * S.fall_depth(np.concatenate([[0.0], zc]), np.concatenate([[Tsf - T0b], Tc - T0b]), F * (Tsf - T0b))))
    p("   D %.1f mm, %g MW/m2: %s" % (d, qq, "; ".join(line)))
    p("      => %s | at the criterion exact x %.3f mm; from cell samples at 3.125 ... 0.195 mm: %s" % (
        "spall at T_s %.0f K, t %.2f s" % first if first else "no spall before the solidus: MELT", 1e3 * ex, " ".join(row)))
