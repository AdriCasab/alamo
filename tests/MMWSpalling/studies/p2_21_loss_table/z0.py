#!/usr/bin/env python3
"""P2-21 z0 -- the offline numbers, each with its formula or its figure reading. Usage: z0.py > Z0.txt

  (1) the 17 Hartlieb pairs as the deck parses them (K), and alpha(T) at 28 GHz at 22 / 500 / 700 / 995 C
  (2) U1's fractions (the two-layer column) and the OLD per-cell kernel's split, by hand
  (3) the spot at the rock from the HE11 aperture and the report's standoff
  (4) Fig. 34 re-read (doe/FIG34.txt, digitised from the re-fetched PDF): incident power and the 28 GHz pick-off
      at beam-on and just before the collapse => the observed transmitted-FRACTION ratio
  (5) the 1-D column (Zhang's laws, 0.2 mm cells, adiabatic bottom, eps 0.7 / 0.3 above T_m, h 56, T_amb 297.5,
      R from eps_r, the on-axis Gaussian with its waist growth, alpha(T_old) through the optical-depth sum,
      backward Euler dt 0.25 s): onset (the first step whose peak reaches T_m), the shallowest melting cell then,
      the transmitted fraction e^-tau_N then and its ratio to the cold value -- for T1a-c, T1d, the flat control,
      the m-scan at omega 20 mm (m = 1, 1.5, 2, 3.5, 8.88), and the same at Fig. 34's MEASURED power ramp; T2's
      surface T at 50 / 100 / 150 s (#28c, 800 W, 25 mm)
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
SIG = 5.670374419e-8
LAM = 1.070687e-2                   # the Zhang deck's c / 28 GHz
ER, TD = 5.02, 0.01065              # the deck's (= the rock file's) Barre pair
TM = 1510.5
T0 = 297.5


def load(name, path):
    keep = sys.argv
    sys.argv = [path]
    ld = importlib.machinery.SourceFileLoader(name, path)
    m = importlib.util.module_from_spec(importlib.util.spec_from_loader(name, ld))
    ld.exec_module(m)
    sys.argv = keep
    return m


RUN = load("run221", os.path.join(HERE, "run.py"))
PAIRS = RUN.table_pairs()
ALPHA = 2.0 * math.pi * math.sqrt(ER) * (TD / LAM)          # the code's association (MMWBeam.H, R6)
N = math.sqrt(ER)
R = ((N - 1.0) / (N + 1.0)) ** 2


def r_of(T, pairs):
    Ts = np.array([p[0] for p in pairs])
    rs = np.array([p[1] for p in pairs])
    return np.interp(T, Ts, rs)            # linear, clamped at both ends (= MMWBeam::alpha_at)


# ------------------------------------------------------------------ Zhang's laws (Zhang.H defaults, P2-20 z0)
def rho(T):
    return np.where(T >= 1273.0, 2300.0, -0.2877 * (T - 273.15) + 2884.0)


def kappa(T):
    return 1.0 / (0.3154 + 3.8e-4 * (np.minimum(T, 1273.0) - 273.15))


def cp(T):
    th = T - 273.15
    return np.where(T <= 1510.5, -2.73e-4 * th * th + 0.894 * th + 788.0, 0.184 * th + 1249.0)


def rhocp_eff(T, fw=200.0, Lm=3.4e5):
    c = cp(T) + np.where((T >= TM - 0.5 * fw) & (T <= TM + 0.5 * fw), Lm / fw, 0.0)
    c = c + np.where((T >= 3233.0) & (T <= 3503.0), 4.8e6 / (3503.0 - 3233.0), 0.0)
    return rho(T) * c


def column(P_of_t, w0, pairs, L=0.05, dz=2e-4, t_end=400.0, dt=0.25, report=()):
    n = int(round(L / dz))
    zc = (np.arange(n) + 0.5) * dz                       # cell centres below the top face (index 0 = top)
    zR = math.pi * w0 * w0 / LAM
    geo = (w0 / (w0 * np.sqrt(1.0 + (zc / zR) ** 2))) ** 2
    T = np.full(n, T0)
    t, out, onset, tau0 = 0.0, {}, None, None
    rep = sorted(report)
    while t < t_end - 1e-9:
        a = ALPHA * r_of(T, pairs)
        tb = np.cumsum(a * dz)
        tt = tb - a * dz
        if tau0 is None:
            tau0 = tb[-1]
        I0 = 2.0 * P_of_t(t) / (math.pi * w0 * w0)
        Q = (1.0 - R) * I0 * geo * (np.exp(-tt) - np.exp(-tb)) / dz
        eps = 0.7 if T[0] <= TM else 0.3
        ql = eps * SIG * (T[0] ** 4 - T0 ** 4) + 56.0 * (T[0] - T0)
        ad = rhocp_eff(T) / dt
        k = kappa(T)
        kf = 2 * k[:-1] * k[1:] / (k[:-1] + k[1:]) / dz ** 2
        diag = ad.copy()
        diag[:-1] += kf
        diag[1:] += kf
        ab = np.zeros((3, n))
        ab[0, 1:] = -kf
        ab[1] = diag
        ab[2, :-1] = -kf
        rhs = ad * T + Q
        rhs[0] -= ql / dz
        T = solve_banded((1, 1), ab, rhs)
        t += dt
        if onset is None and T.max() >= TM:
            im = int(np.argmax(T >= TM))
            onset = dict(t=t, z_melt=zc[im], z_pk=zc[int(np.argmax(T))], tr=math.exp(-tb[-1]),
                         ratio=math.exp(-tb[-1]) / math.exp(-tau0), a_max=a.max(), tau=tb[-1])
        while rep and t >= rep[0] - 1e-9:
            out[rep.pop(0)] = T[0]
    return onset, out, tau0


def fig34():
    p = os.path.join(HERE, "doe", "FIG34.txt")
    rows = []
    for ln in open(p):
        v = ln.split()
        if len(v) == 4 and v[0].isdigit():
            try:
                rows.append((float(v[0]), float(v[1]), float(v[2])))
            except ValueError:
                pass
    return rows


def ramp():
    """the legs' schedule (run.py ramp_schedule, the user's decision 2026-11-06) as P(t) for the 1-D column"""
    te, P = RUN.ramp_schedule(480.0)          # the digitised trace covers 480 s after beam-on; held after
    return lambda t: P[min(len(P) - 1, int(np.searchsorted(te, t, side="right")))]


def main():
    print(__doc__.split("\n\n")[0])
    print("=" * 110)
    print("(1) the table (validation/shared/loss_table_hartlieb.inputs) as parsed: %d pairs (T [K], r)" % len(PAIRS))
    print("    " + "; ".join("%.2f %.2f" % p for p in PAIRS))
    print("    alpha = 2 pi sqrt(eps_r) (tan_delta / lambda) = 2 pi sqrt(%.2f) (%.5f / %.6e) = %.4f 1/m (cold); R = "
          "((n-1)/(n+1))^2 = %.4f" % (ER, TD, LAM, ALPHA, R))
    for tc in (22, 500, 700, 995):
        print("    alpha(%4d C = %7.2f K) = %.4f x r %.2f = %6.2f 1/m" % (tc, tc + 273.15, ALPHA,
                                                                       float(r_of(tc + 273.15, PAIRS)),
                                                                       ALPHA * float(r_of(tc + 273.15, PAIRS))))
    T1 = 474.15 + (1.0 - 0.98) / (1.15 - 0.98) * (573.15 - 474.15)
    print("    r returns to 1 at %.1f K (the 0.98 dip: b1's 486 K)" % T1)
    print("    b2: at the IC 297.5 K r = %.5f => 1 - exp(-alpha r 0.05) = %.5f (bar 0.5034 +- 0.002)" % (
        float(r_of(297.5, PAIRS)), 1 - math.exp(-ALPHA * float(r_of(297.5, PAIRS)) * 0.05)))

    print("\n(2) U1, the two-layer column (alpha 4 1/m over the top 0.5 m, 1 1/m below; the beam enters at z0 = 1):")
    print("    top 1 - e^-(4 x 0.5) = %.6f; bottom e^-2 (1 - e^-(1 x 0.5)) = %.6f; transmitted e^-2.5 = %.6f" % (
        1 - math.exp(-2), math.exp(-2) * (1 - math.exp(-0.5)), math.exp(-2.5)))
    ob = math.exp(-0.5) * (1 - math.exp(-0.5))
    print("    the OLD per-cell kernel (each bottom cell's own alpha 1 over the WHOLE path above it, s from 0.5 to 1): "
          "bottom e^-0.5 (1 - e^-0.5) = %.4f, total %.4f (> 1: the block absorbs more than arrives)" % (
              ob, 1 - math.exp(-2) + ob))

    print("\n(3) the spot at the rock (DOE p. 976 / p. 1134 of the text dump): HE11 1/e^2 diameter 25.4 mm at the "
          "32.5 mm aperture => w0 12.7 mm, z_R = pi w0^2 / lambda")
    zR = math.pi * 0.0127 ** 2 / LAM
    for s in (0.5 * 0.0325, 2.0 * 0.0325):
        print("    standoff %.2f mm (%.1f aperture diameters): w = w0 sqrt(1 + (s/z_R)^2) = %.1f mm (z_R %.1f mm)" % (
            1e3 * s, s / 0.0325, 1e3 * 0.0127 * math.sqrt(1 + (s / zR) ** 2), 1e3 * zR))

    rows = fig34()
    print("\n(4) Fig. 34 re-read (doe/fig34_digitise.py -> doe/FIG34.txt, the lower panel; 10 s windows):")
    for ln in open(os.path.join(HERE, "doe", "FIG34.txt")).read().splitlines()[-4:]:
        print("    " + ln)
    print("    ⚠ the packet's Sources read 'beam on ~200 s at ~2 kW rising to ~2.3 kW': the digitised trace is 1.46 kW at "
          "beam-on, 1.96 kW just before the collapse (mean 1.62 kW); ~2.0-2.3 kW is the level AFTER onset (510-1000 s)")
    print("    the transmitted-FRACTION ratio (pick-off / incident) does not fall monotonically: it dips to 0.79 ~90 s "
          "after beam-on, then RISES above its beam-on value (1.17) before the collapse (1.05) -- the planner's "
          "'>= 0.75' holds as an envelope")

    print("\n(5) the 1-D column (on axis; 3-D runs ~2x later, P2-20):")
    print("    %-34s %8s %8s %8s %8s %8s %8s" % ("case", "onset s", "z_melt", "z_peak", "e^-tau", "/ cold",
                                                 "a_max"))

    def line(lab, o):
        if o is None:
            print("    %-34s %8s" % (lab, "none"))
        else:
            print("    %-34s %8.1f %6.2fmm %6.2fmm %8.3f %8.3f %8.1f" % (lab, o["t"], 1e3 * o["z_melt"], 1e3 * o["z_pk"],
                                                                     o["tr"], o["ratio"], o["a_max"]))
    P2 = lambda t: 2000.0
    for lab, w0 in (("T1a  w0 13 mm, table, 2 kW", 0.013), ("T1b  w0 16.5 mm, table, 2 kW", 0.0165),
                    ("T1c  w0 20 mm, table, 2 kW", 0.020)):
        line(lab, column(P2, w0, PAIRS)[0])
    line("T1d  w0 16.5 mm, 100 mm block", column(P2, 0.0165, PAIRS, L=0.10)[0])
    flat = [(295.15, 1.0)]
    for w0 in (0.0165, 0.020):
        line("flat w0 %.1f mm, 2 kW" % (1e3 * w0), column(P2, w0, flat)[0])
    print("    the m-scan at w0 20 mm, 2 kW (r_m = 1 + (r - 1)(m - 1)/(r_max - 1)):")
    for m in (1.0, 1.5, 2.0, 3.5, 8.88):
        line("  m = %.2f" % m, column(P2, 0.020, RUN.scaled(m) if m < 8.879 else PAIRS)[0])
    print("    (second agent's 1-D, spot and read time unstated: first melt 8.1 / 6.6 / 5.6 / 3.9 / 1.9 mm, "
          "transmission vs cold 1.00 / 0.77 / 0.63 / 0.43 / 0.22)")
    Pr = ramp()
    print("    the same at Fig. 34's MEASURED incident power (run.py ramp_schedule: the legs' schedule, the user's decision "
          "2026-11-06):")
    for lab, w0, pairs in (("  w0 13 mm, table", 0.013, PAIRS), ("  w0 16.5 mm, table", 0.0165, PAIRS),
                           ("  w0 20 mm, table", 0.020, PAIRS), ("  w0 16.5 mm, flat", 0.0165, flat),
                           ("  w0 20 mm, flat", 0.020, flat)):
        line(lab, column(Pr, w0, pairs, t_end=600.0)[0])
    print("    the m-scan at w0 20 mm at the measured power (the M legs' schedule):")
    for m in (1.0, 1.5, 2.0, 3.5, 8.88):
        line("  m = %.2f" % m, column(Pr, 0.020, RUN.scaled(m) if m < 8.879 else PAIRS, t_end=600.0)[0])
    print("\n    T2 (#28c, 800 W, the 25 mm disc, w0 14.1 mm): surface T [K] at 50 / 100 / 150 s")
    for lab, pairs in (("table", PAIRS), ("flat", flat)):
        o, out, _ = column(lambda t: 800.0, 0.0141, pairs, L=0.025, t_end=150.0, report=(50, 100, 150))
        print("    %-6s %s   (#28c uncorrected 987 / 1141 / 1201 K)" % (lab, " / ".join("%.0f" % out[t]
                                                                                   for t in (50, 100, 150))))
    return 0


if __name__ == "__main__":
    sys.exit(main())
