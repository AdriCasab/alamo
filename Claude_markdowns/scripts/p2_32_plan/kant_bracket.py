#!/usr/bin/env python3
"""P2-32 planning, zero compute: Kant's bracket legs under the stress-fraction depth (f = 0.75), by a 1-D column.
For each leg (Kant 2016 Table 4 h and flame temperature, P2-15 run.py KT; Gotthard / Bethel conduction constants):
the span the leg REQUIRES, L_req = C_L(T_s) x, (a) at the first firing (T_s 820 K, the logged level on these legs) under
today's rule (0.4 x the tangent length) and under the stress depth, and (b) its MINIMUM over the heating history after
that (to 200 s, the legs' stop), with the surface temperature where the minimum sits. A leg spalls at a span L iff
L_req <= L at some time: (a) is P2-16's discriminant (earliest firing), (b) is what a run with the span ON would do.
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stress_depth as S

KT = {130: (747, 1250), 140: (690, 1175), 150: (633, 1100), 160: (576, 1025)}
ROCK = {"G": (1923.5, 2594.0), "B": (1540.6, 2566.0)}
LOGGED = {"KG130": 9.915, "KG140": 12.010, "KG150": 14.956, "KG160": 19.295, "KB150": 11.992, "KB160": 15.430}


def run(leg, f=0.75, t_max=200.0, dz=0.1e-3, depth=0.064):
    a, rho = ROCK[leg[1]]
    h, th = KT[int(leg[2:])]
    Tg, T0 = th + 273.15, 293.15
    lam = lambda T: a / (T - 273.15 + 563.0)
    n = int(round(depth / dz))
    z = (np.arange(n) + 0.5) * dz
    T = np.full(n, T0)
    dt = 0.4 * dz * dz / (lam(T0) / (rho * S.cp(T0)))
    t, Ts, nxt, rows = 0.0, T0, 820.0, []
    while t < t_max:
        k = lam(T[0])
        Ts = T[0]
        for _ in range(30):
            g = h * (Tg - Ts) - 10.0 * (Ts - T0) - 0.8 * S.SIG * (Ts ** 4 - T0 ** 4) - k * (Ts - T[0]) / (0.5 * dz)
            Ts -= g / (-h - 10.0 - 4 * 0.8 * S.SIG * Ts ** 3 - k / (0.5 * dz))
        if Ts >= nxt:
            q = h * (Tg - Ts) - 10.0 * (Ts - T0) - 0.8 * S.SIG * (Ts ** 4 - T0 ** 4)
            x = S.depth_true(z, T, Ts, T0, f)
            rows.append((t, Ts, S.C_L(Ts) * 0.4 * lam(Ts) * (Ts - T0) / q, S.C_L(Ts) * x, x))
            nxt += 10.0
        L = lam(T)
        Lf = 2.0 * L[:-1] * L[1:] / (L[:-1] + L[1:])
        F = Lf * (T[1:] - T[:-1]) / dz
        div = np.zeros(n)
        div[:-1] += F
        div[1:] -= F
        div[0] += k * (Ts - T[0]) / (0.5 * dz)
        T = T + dt * (div / dz + 0.0) / (rho * S.cp(T))
        t += dt
    return rows, Ts


print("leg     | at 820 K: t s, today's L_req mm (logged), stress-depth L_req mm | minimum of the stress-depth L_req after onset: mm at T_s K, t s | T_s at 200 s")
res = {}
for leg in ("KG130", "KG140", "KG150", "KG160", "KB150", "KB160"):
    rows, Tend = run(leg)
    t0, Ts0, Lt0, Ls0, x0 = rows[0]
    i = int(np.argmin([r[3] for r in rows]))
    res[leg] = (Ls0, rows[i][3])
    print("%-7s | %6.2f  %7.3f (%6.3f)  %7.3f (x %.3f mm) | %7.3f at %6.0f K, %6.1f s (x %.3f mm) | %6.0f K" % (
        leg, t0, 1e3 * Lt0, LOGGED[leg], 1e3 * Ls0, 1e3 * x0, 1e3 * rows[i][3], rows[i][1], rows[i][0], 1e3 * rows[i][4], Tend))
    if leg == "KG140":
        print("          KG140 history (T_s K: today's L_req, stress L_req mm): " + "  ".join(
            "%.0f: %.1f, %.1f" % (r[1], 1e3 * r[2], 1e3 * r[3]) for r in rows[::3]))
for lab, k in (("at the first firing (P2-16's discriminant)", 0), ("at the minimum over the history (the span ON for 200 s)", 1)):
    g = (res["KG130"][k], min(res["KG140"][k], res["KG150"][k], res["KG160"][k]))
    b = (res["KB150"][k], res["KB160"][k])
    lo, hi = max(g[0], b[0]), min(g[1], b[1])
    print("%s: Gotthard window [%.3f, %.3f) mm, Bethel window [%.3f, %.3f) mm -> one span in [%.3f, %.3f) mm: %s" % (
        lab, 1e3 * g[0], 1e3 * g[1], 1e3 * b[0], 1e3 * b[1], 1e3 * lo, 1e3 * hi, "EXISTS" if lo < hi else "EMPTY by %.1f %%" % (100 * (lo / hi - 1))))
