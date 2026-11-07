#!/usr/bin/env python3
"""P2-23 z0 -- the offline numbers, each with its formula. Usage: z0.py > Z0.txt

The Grimsel law is transcribed from src/Numeric/Material/Grimsel.H (rho 2650; c_p = 1370 - 178000/(theta + 271);
lambda = 2000/(theta + 563); Table 3 alpha / E/E0 / K_Ic, clamped at both ends; the strain = the trapezoid integral of
alpha on a 1 K grid from 20 C to 1000 C, held above; nu 0.30). theta in C.
"""
import importlib.machinery
import importlib.util
import math
import os
import sys

import numpy as np
from scipy.linalg import solve_banded

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
T3_th = [20, 100, 200, 300, 400, 450, 500, 550, 600, 650]
T3_alpha = [8.4, 8.7, 9.1, 9.5, 16.4, 19.9, 23.4, 26.9, 28.3, 28.4]
NU, RHO = 0.30, 2650.0
LAM_HE11, W0_AP, STANDOFF = 1.070687e-2, 0.0127, 0.0079
SIG = 5.670374419e-8


def cp(T):
    return 1370.0 - 178000.0 / ((T - 273.15) + 271.0)


def lam(T):
    return 2000.0 / ((T - 273.15) + 563.0)


def kappa(T):
    return lam(T) / (RHO * cp(T))


def alpha_th(th):
    return 1e-6 * np.interp(th, T3_th, T3_alpha)          # clamped at both ends, as Interp


_GRID = np.arange(20.0, 1000.0 + 0.5, 1.0)
_CUM = np.concatenate([[0.0], np.cumsum(0.5 * (alpha_th(_GRID[:-1]) + alpha_th(_GRID[1:])))])


def strain(T):
    th = T - 273.15
    if th <= 20.0:
        return 0.0
    i = int(th - 20.0)
    if i >= len(_GRID) - 1:
        return float(_CUM[-1])
    thi = 20.0 + i
    return float(_CUM[i] + 0.5 * (alpha_th(thi) + alpha_th(th)) * (th - thi))


def bk(Ts, T0, x, L):
    e = strain(Ts) - strain(T0)
    CL = math.pi / math.sqrt(6.0 * (1.0 + NU) * e)
    return CL, e, CL * x / L


def load(name, path):
    keep = sys.argv
    sys.argv = [path]
    ld = importlib.machinery.SourceFileLoader(name, path)
    m = importlib.util.module_from_spec(importlib.util.spec_from_loader(name, ld))
    ld.exec_module(m)
    sys.argv = keep
    return m


def schedule():
    import shlex
    a = shlex.split(open(os.path.join(ROOT, "tests", "MMWSpalling", "studies", "p2_22_inrock_spread", "legs", "A",
                                      "argv.sh")).read())
    te = [float(x) for x in next(t for t in a if t.startswith("beam.schedule.t_end=")).split("=", 1)[1].split()]
    P = [float(x) for x in next(t for t in a if t.startswith("beam.schedule.P=")).split("=", 1)[1].split()]
    return lambda t: P[min(len(P) - 1, int(np.searchsorted(te, t, side="right")))]


def column(P_of_t, w0, pairs, a0, R, n_idx, L=0.05, dz=2e-4, t_end=400.0, dt=0.25, T0=297.5, Tm=1510.5,
           eps=0.8, h=56.0):
    """the 1-D column on the axis with the GRIMSEL law (no latent heat: the Grimsel phase has none), the table on the
    cold alpha through the optical-depth sum, the in-rock spread (z_R x n, P2-22), backward Euler"""
    nc = int(round(L / dz))
    zc = (np.arange(nc) + 0.5) * dz
    zR = n_idx * math.pi * w0 * w0 / LAM_HE11
    geo = 1.0 / (1.0 + (zc / zR) ** 2)
    Ts_, rs_ = np.array([p[0] for p in pairs]), np.array([p[1] for p in pairs])
    T = np.full(nc, T0)
    t = 0.0
    while t < t_end - 1e-9:
        a = a0 * np.interp(T, Ts_, rs_)
        tb = np.cumsum(a * dz)
        tt = tb - a * dz
        Q = (1.0 - R) * 2.0 * P_of_t(t) / (math.pi * w0 * w0) * geo * (np.exp(-tt) - np.exp(-tb)) / dz
        ql = eps * SIG * (T[0] ** 4 - T0 ** 4) + h * (T[0] - T0)
        ad = RHO * cp(T) / dt
        k = lam(T)
        kf = 2 * k[:-1] * k[1:] / (k[:-1] + k[1:]) / dz ** 2
        diag = ad.copy()
        diag[:-1] += kf
        diag[1:] += kf
        ab = np.zeros((3, nc))
        ab[0, 1:] = -kf
        ab[1] = diag
        ab[2, :-1] = -kf
        rhs = ad * T + Q
        rhs[0] -= ql / dz
        T = solve_banded((1, 1), ab, rhs)
        t += dt
        if T.max() >= Tm:
            rise = T - T0
            ip = int(np.argmax(T))
            ef = next((zc[i] for i in range(ip, nc) if rise[i] <= rise[ip] / math.e), float("nan"))
            return dict(t=t, z_melt=zc[int(np.argmax(T >= Tm))], z_pk=zc[ip], ell_pk_abs=ef, T_top=T[0])
    return None


def main():
    print(__doc__.split("\n\n")[0])
    print("=" * 110)
    zR = math.pi * W0_AP ** 2 / LAM_HE11
    w = W0_AP * math.sqrt(1.0 + (STANDOFF / zR) ** 2)
    print("(1) the spot: HE11 1/e^2 diameter 25.4 mm => w0 12.7 mm at the aperture, z_R = pi w0^2/lambda = %.1f mm; "
          "standoff 7.9 mm (DOE Appendix Table 2) => w = 12.7 sqrt(1 + (7.9/%.1f)^2) = %.2f mm (the deck: 12.88)" % (
              1e3 * zR, 1e3 * zR, 1e3 * w))
    Ts = np.linspace(297.0, 1600.0, 2607)
    ks = kappa(Ts)
    print("(2) kappa = lambda/(rho c_p) (Grimsel law) over 297-1600 K: max %.4e m2/s at %.1f K (room T), %.3e at 900 K, "
          "%.3e at 1500 K" % (ks.max(), Ts[int(np.argmax(ks))], kappa(900.0), kappa(1500.0)))
    for dxm, dtd, lab in ((3.125e-3, 0.5, "U (unit/face_beam)"), (1.5625e-3, 0.2, "B0, B0F (the deck)"),
                          (1.5625e-3, 0.016, "P"), (0.78125e-3, 0.05, "B1")):
        lim = dxm ** 2 / (6.0 * ks.max())
        print("    dx %.4f mm: explicit limit dx^2/(6 kappa_max) = %.4f s; %-20s dt %.3f s (%.0f %% of the limit)" % (
            1e3 * dxm, lim, lab, dtd, 100 * dtd / lim))
    td = 0.6970986507425017
    a = 2.0 * math.pi * math.sqrt(1.0) * (td / 2.19e-3)
    print("(3) P's alpha = 2 pi sqrt(1) (0.6970986507425017 / 2.19e-3) = %r 1/m (exactly 2000: %s); R = 0; absorption "
          "depth 1/alpha = %.2f mm" % (a, a == 2000.0, 1e3 / a))
    I0 = 2.0 * 150.0 / (math.pi * 0.01288 ** 2)
    f_top = 1.0 - math.exp(-2000.0 * 1.5625e-3)
    rate = I0 * f_top / (1.5625e-3 * RHO * cp(297.5))
    print("(4) P at 150 W: peak incident (= absorbed, R 0) flux 2P/(pi w0^2) = %.3f MW/m2; the top 1.5625 mm cell takes "
          "1 - e^(-alpha dz) = %.4f of it => %.0f K/s at room T (rho c_p %.3e J/m3K). Grimsel's minimum net flux for B at "
          "the 12 mm span: 0.432 MW/m2 (P2-17 Z2)" % (I0 / 1e6, f_top, rate, RHO * cp(297.5)))
    print("    the conduction length after t s: sqrt(kappa t) = %.2f / %.2f / %.2f mm at 3 / 7 / 15 s -- the layer B reads "
          "grows past the 0.5 mm absorption depth within a second" % tuple(1e3 * math.sqrt(ks.max() * t)
                                                                           for t in (3, 7, 15)))
    for dxm in (3.125e-3, 1.5625e-3, 0.78125e-3):
        print("(5) the window 15 mm at dx %.4f mm: %.1f cells (L.Dmax = floor = %d)" % (
            1e3 * dxm, 0.015 / dxm, int(0.015 / dxm)))
    print("(6) bk_lb by hand: C_L = pi/sqrt(6 (1 + nu) eps_th), eps_th = strain(T_s) - strain(297.5 K), x = mu c_f "
          "(window) = 0.4 x 0.015 m on a top line (mu 1), L 0.012 m")
    for Tsv in (850.0, 900.0, 1000.0, 1200.0, 1400.0):
        CL, e, b = bk(Tsv, 297.5, 0.4 * 0.015, 0.012)
        print("    T_s %6.1f K: eps_th %.5f, C_L %.3f, bk_lb %.3f" % (Tsv, e, CL, b))
    print("    ⚠ THE CLAMP, read from Grimsel.H (the packet said 'the table clamps at 873 K, so 900 K and above read the "
          "SAME strain'): Table 3's alpha, E/E0 and K_Ic are flat above 650 C (923 K; 600 and 650 C nearly equal), but "
          "the STRAIN integrates the held alpha (28.4e-6/K) up to 1000 C (1273 K) and is held only above that -- so "
          "bk_lb keeps falling to 1273 K, as the rows above show")
    print("    the layer limit: B blocks when C_L mu c_f ell > L, i.e. ell > L/(C_L c_f) = %.2f mm at T_s 850 K, %.2f mm "
          "at 1200 K (mu 1)" % (0.012 / (bk(850.0, 297.5, 1, 1)[0] * 0.4) * 1e3,
                                0.012 / (bk(1200.0, 297.5, 1, 1)[0] * 0.4) * 1e3))
    # (7) the 1-D column
    R21 = load("run221_for_p223", os.path.join(ROOT, "tests", "MMWSpalling", "studies", "p2_21_loss_table", "run.py"))
    pairs = R21.table_pairs()
    n = math.sqrt(5.02)
    a0 = 2.0 * math.pi * n * (0.01065 / LAM_HE11)
    R = ((n - 1) / (n + 1)) ** 2
    P = schedule()
    print("(7) the 1-D column on the axis, GRIMSEL law (no latent heat), the measured ramp, w0 12.88 mm, the table on "
          "alpha %.4f (R %.4f), z_R x n = %.3f (in-rock spread), eps 0.8, h 56, 50 mm, 0.2 mm cells, dt 0.25 s" % (
              a0, R, n))
    for lab, pr in (("table", pairs), ("flat", [(295.15, 1.0)])):
        o = column(P, 0.01288, pr, a0, R, n)
        if o is None:
            print("    %-6s no melt (T_m 1510.5) by 400 s" % lab)
        else:
            print("    %-6s first melt %.1f s at %.2f mm (peak %.2f mm; top %.0f K); the peak's absolute e-fold depth "
                  "%.1f mm" % (lab, o["t"], 1e3 * o["z_melt"], 1e3 * o["z_pk"], o["T_top"], 1e3 * o["ell_pk_abs"]))
    print("    beside the Zhang material at w0 12.88 mm (P2-22 leg A, 3-D): first melt 50 s at the peak 1.95 mm, "
          "bc_ell_p 22.0 mm")
    return 0


if __name__ == "__main__":
    sys.exit(main())
