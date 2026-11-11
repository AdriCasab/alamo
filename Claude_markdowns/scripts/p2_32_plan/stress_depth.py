#!/usr/bin/env python3
"""P2-32 planning, zero compute (2026-11-11): the stress-fraction flake depth, measured before the packet states a bar.

The Grimsel law is transcribed from src/Numeric/Material/Grimsel.H (Table 3 read from the file itself; strain = the
integral of alpha from 20 C, held above 1000 C; sigma_thermal = E(T_hot) [eps(T_hot) - eps(T_cold)] / (1 - nu);
C_L = pi / sqrt(6 (1 + nu) eps_th)). Checked against two logged values: C_L(851.21 K) = 12.727 (P2-15 CHAIN.log:23)
and C_L(758.15 K) = 15.243 (P2-16 RULES.txt c0b).
Thermal properties: rho 2650, c_p = 1370 - 178000/(theta + 271), lambda = 2000/(theta + 563), theta in C.
A 1-D explicit column at dz 0.05 mm supplies "true" profiles: (a) the flat flame fixture (unit/face_removal: h 700,
T_g 1866 K, eps 0.8, h_conv 10, T0 293.15 K); (b) the beam column (studies/beam_column_1110: Beer-Lambert q (1-R)
already absorbed, eps 0.8, h_conv 56, T0 297.5 K).
"""
import math
import os
import re
import sys

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
s = open(os.path.join(ROOT, "src", "Numeric", "Material", "Grimsel.H")).read()


def arr(n):
    m = re.search(n + r"\[NT3\]\s*=\s*\{([^}]*)\}", s, re.S)
    return np.array([float(x) for x in m.group(1).split(",") if x.strip()])


TH, AL, EE = arr("T3_th"), arr("T3_alpha"), arr("T3_EoE0")
E0, NU, RHO, SIG, R = 53.3e9, 0.30, 2650.0, 5.670374e-8, 0.531294
TG = np.arange(20.0, 1000.0 + 1e-9, 0.25)
EPS = np.concatenate([[0.0], np.cumsum(0.5 * (np.interp(TG[1:], TH, AL) + np.interp(TG[:-1], TH, AL)) * 1e-6 * 0.25)])


def strain(T):
    return np.interp(np.clip(np.asarray(T, float) - 273.15, 20.0, 1000.0), TG, EPS)


def sig(T, T0):
    return R * E0 * np.interp(np.asarray(T, float) - 273.15, TH, EE) * (strain(T) - strain(T0)) / (1 - NU)


def C_L(Ts, T0=293.15):
    return math.pi / math.sqrt(6.0 * (1 + NU) * float(strain(Ts) - strain(T0)))


def cp(T):
    return 1370.0 - 178000.0 / (T - 273.15 + 271.0)


def lam(T):
    return 2000.0 / (T - 273.15 + 563.0)


def fall_depth(zs, vs, target):
    """ProfileFallDepth's rule on samples (zs[0] = 0 the face): ln-linear between the bracketing samples"""
    for a in range(len(zs) - 1):
        if vs[a + 1] <= target:
            if vs[a + 1] > 0:
                return zs[a] + (zs[a + 1] - zs[a]) * math.log(vs[a] / target) / math.log(vs[a] / vs[a + 1])
            return zs[a] + (zs[a + 1] - zs[a]) * (vs[a] - target) / (vs[a] - vs[a + 1])
    return float("nan")


def depth_true(z, T, Ts, T0, f):
    sg = sig(np.concatenate([[Ts], T]), T0)
    zz = np.concatenate([[0.0], z])
    i = np.argmax(sg <= f * sg[0])
    if not (sg <= f * sg[0]).any():
        return float("nan")
    return zz[i - 1] + (zz[i] - zz[i - 1]) * (sg[i - 1] - f * sg[0]) / (sg[i - 1] - sg[i])


def samples(z, T, dzc, n):
    """cell AVERAGES of the fine profile on cells of size dzc, at their centres"""
    dzf = z[1] - z[0]
    m = int(round(dzc / dzf))
    k = min(n, len(T) // m)
    return (np.arange(k) + 0.5) * dzc, T[:k * m].reshape(k, m).mean(axis=1)


def est_sigma(zc, Tc, Ts, T0, f):
    """A: ProfileFallDepth on the STRESS samples (face value from T_s, then the cells)"""
    return fall_depth(np.concatenate([[0.0], zc]), np.concatenate([[float(sig(Ts, T0))], sig(Tc, T0)]), f * float(sig(Ts, T0)))


def est_rise(zc, Tc, Ts, T0, f):
    """B: ln-linear in the temperature RISE between samples, the stress inverted exactly inside the bracket"""
    zz = np.concatenate([[0.0], zc]); rr = np.concatenate([[Ts - T0], Tc - T0])
    tgt = f * float(sig(Ts, T0))
    for a in range(len(zz) - 1):
        if float(sig(T0 + rr[a + 1], T0)) <= tgt:
            lo, hi = 0.0, 1.0
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                r = rr[a] * (rr[a + 1] / rr[a]) ** mid if rr[a + 1] > 0 else rr[a] + mid * (rr[a + 1] - rr[a])
                if float(sig(T0 + r, T0)) > tgt:
                    lo = mid
                else:
                    hi = mid
            return zz[a] + 0.5 * (lo + hi) * (zz[a + 1] - zz[a])
    return float("nan")


def column(T0, t_max, stop_at, dz=0.05e-3, depth=0.05, h=0.0, Tg=0.0, hc=0.0, Tamb=293.15, q=0.0, alpha=0.0, eps=0.8,
           snaps=()):
    """explicit 1-D column; the face temperature from the face balance on the top half cell (Newton);
    returns [(t, Ts, z, T)] at the first time Ts >= each value of `snaps` (K), and stops at stop_at or t_max"""
    n = int(round(depth / dz))
    z = (np.arange(n) + 0.5) * dz
    ze = dz * np.arange(n + 1)
    S = q * (np.exp(-alpha * ze[:-1]) - np.exp(-alpha * ze[1:])) / dz if alpha > 0 else np.zeros(n)
    T = np.full(n, T0)
    dt = 0.4 * dz * dz / (lam(T0) / (RHO * cp(T0)))
    out, todo, t, Ts = [], list(snaps), 0.0, T0
    while t < t_max:
        k = lam(T[0])
        Ts = T[0]
        for _ in range(30):
            g = h * (Tg - Ts) - hc * (Ts - Tamb) - eps * SIG * (Ts ** 4 - Tamb ** 4) - k * (Ts - T[0]) / (0.5 * dz)
            dg = -h - hc - 4 * eps * SIG * Ts ** 3 - k / (0.5 * dz)
            Ts -= g / dg
        while todo and Ts >= todo[0]:
            out.append((t, Ts, z.copy(), T.copy(), todo.pop(0)))
        if Ts >= stop_at:
            break
        L = lam(T)
        Lf = 2.0 * L[:-1] * L[1:] / (L[:-1] + L[1:])
        F = Lf * (T[1:] - T[:-1]) / dz
        div = np.zeros(n)
        div[:-1] += F
        div[1:] -= F
        div[0] += k * (Ts - T[0]) / (0.5 * dz)
        T = T + dt * (div / dz + S) / (RHO * cp(T))
        t += dt
    return out


def main():
    p = print
    p("law check: C_L(851.21 K) = %.3f (log 12.727); C_L(758.15 K) = %.3f (log 15.243); R sigma(850 K) = %.1f MPa (UCS 169)" % (
        C_L(851.21), C_L(758.15), float(sig(850.0, 293.15)) / 1e6))
    T0 = 293.15
    p("\n1. the fraction that reproduces 0.4 ell on an exponential skin, by surface temperature")
    for Ts in (800.0, 820.0, 843.757, 846.0, 850.43, 860.0, 880.0, 900.0, 1000.0, 1200.0):
        p("   T_s %8.3f K: f = %.4f;  with f = 0.75 the depth is %.4f ell (x %.3f of 0.4 ell)" % (
            Ts, float(sig(T0 + (Ts - T0) * math.exp(-0.4), T0) / sig(Ts, T0)),
            depth_true(np.linspace(1e-4, 4, 40000), T0 + (Ts - T0) * np.exp(-np.linspace(1e-4, 4, 40000)), Ts, T0, 0.75),
            depth_true(np.linspace(1e-4, 4, 40000), T0 + (Ts - T0) * np.exp(-np.linspace(1e-4, 4, 40000)), Ts, T0, 0.75) / 0.4))
    p("\n2. exponential skin ell 1.5 mm, T_s 846 K, cell-average samples: the two estimators against the exact depth")
    zf = (np.arange(1200) + 0.5) * 0.0125e-3
    Tf = T0 + (846.0 - T0) * np.exp(-zf / 1.5e-3)
    ex = depth_true(zf, Tf, 846.0, T0, 0.75)
    for dzc in (2e-3, 1e-3, 0.5e-3, 0.25e-3):
        zc, Tc = samples(zf, Tf, dzc, 7)
        a, b = est_sigma(zc, Tc, 846.0, T0, 0.75), est_rise(zc, Tc, 846.0, T0, 0.75)
        p("   dz %5.3f mm: exact %.4f mm; A (ln-linear in stress) %.4f (%+.1f %%); B (ln-linear in rise, stress inverted) %.4f (%+.1f %%)" % (
            1e3 * dzc, 1e3 * ex, 1e3 * a, 100 * (a / ex - 1), 1e3 * b, 100 * (b / ex - 1)))
    p("\n3. the flat flame fixture at its first firing (h 700, T_g 1866 K; the face reaches 846 K)")
    t, Ts, z, T, _ = column(T0, 20.0, 846.0, h=700.0, Tg=1866.0, hc=10.0, Tamb=293.15, snaps=(846.0,))[0]
    ex = depth_true(z, T, Ts, T0, 0.75)
    q_net = 700.0 * (1866.0 - Ts) - 10.0 * (Ts - 293.15) - 0.8 * SIG * (Ts ** 4 - 293.15 ** 4)
    p("   t %.3f s (P2-31 m1: 2.38 s at 0.5 mm); tangent depth k dT/q_net = %.3f mm; exact stress depth at f 0.75 = %.4f mm; bk = C_L x / 12 mm = %.3f" % (
        t, 1e3 * lam(Ts) * (Ts - T0) / q_net, 1e3 * ex, C_L(Ts) * ex / 0.012))
    p("   today's rule on the same samples: 0.4 x the two-cell length (P2-31 m1 measured ell 1.623 / 1.480 / 1.437 mm => 0.649 / 0.592 / 0.575 mm)")
    for dzc in (2e-3, 1e-3, 0.5e-3, 0.25e-3):
        zc, Tc = samples(z, T, dzc, 15)
        two = dzc / math.log((Tc[0] - T0) / (Tc[1] - T0))
        a, b = est_sigma(zc, Tc, Ts, T0, 0.75), est_rise(zc, Tc, Ts, T0, 0.75)
        p("   dz %5.3f mm: two-cell ell %.3f mm -> 0.4 ell = %.4f mm | A %.4f mm (%+.1f %% of exact) | B %.4f mm (%+.1f %%)" % (
            1e3 * dzc, 1e3 * two, 0.4e3 * two, 1e3 * a, 100 * (a / ex - 1), 1e3 * b, 100 * (b / ex - 1)))
    p("\n4. the beam column (T0 297.5 K, eps 0.8, h_conv 56): the stress depth and bk from the criterion (846 K) to the solidus (1410.5 K)")
    T0b = 297.5
    for d, qq in ((0.5, 2.0), (1.0, 2.0), (2.0, 2.0), (4.0, 2.0), (8.0, 2.0), (2.0, 4.0)):
        snaps = (846.0, 900.0, 1000.0, 1100.0, 1200.0, 1273.0, 1350.0, 1410.0)
        res = column(T0b, 60.0, 1410.5, q=qq * 1e6, alpha=1.0 / (d * 1e-3), hc=56.0, Tamb=297.5, snaps=snaps)
        p("   D %.1f mm, %g MW/m2" % (d, qq))
        first_pass = None
        for t, Ts, z, T, lv in res:
            x = depth_true(z, T, Ts, T0b, 0.75)
            bk = C_L(Ts, T0b) * x / 0.012
            if first_pass is None and bk <= 1.0:
                first_pass = lv
            ipk = int(np.argmax(T))
            p("      T_s %6.0f K at t %6.3f s: stress depth %.3f mm, C_L %.2f, bk %.3f %s; peak %.1f K at %.3f mm" % (
                lv, t, 1e3 * x, C_L(Ts, T0b), bk, "pass" if bk <= 1 else "WITHDRAWN", T[ipk], 1e3 * z[ipk]))
        p("      => %s" % ("spall at T_s %.0f K" % first_pass if first_pass else "no spall before the solidus: MELT"))
        t, Ts, z, T, _ = res[0]
        ex = depth_true(z, T, Ts, T0b, 0.75)
        row = []
        for n in (16, 32, 64, 128, 256, 512):
            dzc = 0.05 / n
            if dzc < 2 * (z[1] - z[0]) - 1e-12:
                row.append("%8s" % "-")
                continue
            zc, Tc = samples(z, T, dzc, int(0.015 / dzc) + 1)
            k = lam(Tc[0])
            Tsf = Tc[0]                                   # the face balance on the coarse cell (losses only)
            for _ in range(30):
                g = -56.0 * (Tsf - 297.5) - 0.8 * SIG * (Tsf ** 4 - 297.5 ** 4) - k * (Tsf - Tc[0]) / (0.5 * dzc)
                Tsf -= g / (-56.0 - 4 * 0.8 * SIG * Tsf ** 3 - k / (0.5 * dzc))
            row.append("%8.3f" % (1e3 * est_rise(zc, Tc, Tsf, T0b, 0.75)))
        p("      at the criterion: exact %.3f mm; estimator B from cell samples with the face-balance T_s, cells 3.125 ... 0.098 mm: %s" % (1e3 * ex, " ".join(row)))


if __name__ == "__main__":
    main()
