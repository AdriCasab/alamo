#!/usr/bin/env python3
"""P2-13 a0 -- the offline onset twin (the advisor's, 2026-10-30), ported.

Uniform stress over the flaw, Tada weight, the code's own Grimsel evaluators:

    W * (R * sigma(T) + p*nu/(1-nu)) * sqrt(a/pi) = K_Ic(T)

  sigma(T) = E(T) * strain(T) / (1 - nu)      E, strain, K_Ic: columns of
             studies/p2_2_props/KIC_CPP.txt, dumped by kic_dump.cc straight out
             of src/ (Grimsel.H Table 3; E0 53.3 GPa, K_Ic0 1.5; no caps)
  W        = 2 INT_0^{pi/2} g(sin th) dth, g = 1.3 - 0.3 xi^1.25 = 3.525546
  R        = spallation.sp.stress_relief = 0.531294 (P2-2d: UCS / sigma(850 K))
  p nu/(1-nu): the confinement term, Removal.H:41 (nu 0.30)

It reproduces the five Sources numbers of the P2-13 packet to 0.05 K (a0) and prints
the onset at every leg's flaw and pressure, i.e. the predictions P1, P3 and P5 rest on.
Usage: offline_onset.py [> A0.txt]   (exit 1 on a failed row)
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", "..", ".."))
KIC = os.path.join(ROOT, "tests", "MMWSpalling", "studies", "p2_2_props", "KIC_CPP.txt")
NU, R, M = 0.30, 0.531294, 13.0
A0 = 8.0e-6


def weight(n=20000):
    g = lambda x: 1.3 - 0.3 * (x ** 1.25 if x > 0.0 else 0.0)
    h = (math.pi / 2) / n
    s = sum((1 if i in (0, n) else (4 if i % 2 else 2)) * g(math.sin(i * h)) for i in range(n + 1))
    return 2 * s * h / 3


W = weight()
T_, KIC_, E_, EPS_ = [], [], [], []
with open(KIC) as f:
    next(f)
    for ln in f:
        t = ln.split()
        T_.append(float(t[0])); KIC_.append(float(t[3])); E_.append(float(t[4])); EPS_.append(float(t[5]))


def at(T, ys):
    """linear read of the 1 K dump, clamped as the evaluators clamp"""
    if T <= T_[0]:
        return ys[0]
    if T >= T_[-1]:
        return ys[-1]
    i = max(0, min(len(T_) - 2, int(T - T_[0])))
    f = (T - T_[i]) / (T_[i + 1] - T_[i])
    return ys[i] + f * (ys[i + 1] - ys[i])


def sigma(T):
    return at(T, E_) * at(T, EPS_) / (1.0 - NU)


def onset(a=A0, p=0.0, r=R):
    """first T in 300-1200 K where the drive reaches the toughness (bisected); None if never"""
    f = lambda T: W * (r * sigma(T) + p * NU / (1.0 - NU)) * math.sqrt(a / math.pi) - at(T, KIC_) * 1e6
    prev, T = f(300.0), 301.0
    while T <= 1200.0:
        cur = f(T)
        if prev < 0.0 <= cur:
            lo, hi = T - 1.0, T
            for _ in range(80):
                m = 0.5 * (lo + hi)
                lo, hi = (m, hi) if f(m) < 0.0 else (lo, m)
            return 0.5 * (lo + hi)
        prev, T = cur, T + 1.0
    return None


def main():
    a_med = A0 * math.log(2.0) ** (1.0 / M)
    rows = [("p 0, a 8 um", onset(), 842.39), ("p 1 MPa", onset(p=1.0e6), 841.29),
            ("p 27 MPa", onset(p=27.0e6), 818.47), ("p 48 MPa", onset(p=48.0e6), 806.92),
            ("p 0, Weibull median a %.3f um" % (a_med * 1e6), onset(a=a_med), 848.66)]
    print("P2-13 a0 -- the offline onset twin (uniform stress, Tada weight, Grimsel.H via KIC_CPP.txt)")
    print("=" * 96)
    print("  W = %.6f   R = %.6f   nu = %.2f   m = %g   (K_Ic, E, strain: %s)"
          % (W, R, NU, M, os.path.relpath(KIC, ROOT)))
    bad = 0
    for nm, got, want in rows:
        ok = got is not None and abs(got - want) <= 0.05
        bad += not ok
        print("  %-34s %9.2f K   advisor %7.2f K   %+6.3f K   %s"
              % (nm, got, want, got - want, "PASS" if ok else "FAIL"))
    print("\n  for the legs (not bars): onset at the Weibull median flaw a0*(ln 2)^(1/m)")
    for a0 in (6.0e-6, 10.0e-6):
        am = a0 * math.log(2.0) ** (1.0 / M)
        print("    a0 %4.1f um: median flaw %.3f um -> %.2f K   (deterministic a = a0: %.2f K)"
              % (a0 * 1e6, am * 1e6, onset(a=am), onset(a=a0)))
    for p in (27.0e6, 48.0e6):
        print("    Weibull median, p %2.0f MPa: %.2f K (shift %+.2f K)"
              % (p / 1e6, onset(a=a_med, p=p), onset(a=a_med, p=p) - onset(a=a_med)))
    print("\na0 %s (%d of %d rows within 0.05 K)" % ("PASS" if not bad else "FAIL", len(rows) - bad, len(rows)))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
