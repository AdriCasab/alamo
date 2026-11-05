#!/usr/bin/env python3
"""P2-17 ord 1 (zero compute) -- z1 and z2. Writes Z1.txt, Z2.txt.

z1  THE CAP-DEAD CROSSOVER. Under B a fired, non-through line passes iff C_L(T_s) mu c_f ell <= L; the flake cap
    acts iff mu c_f ell > flake_h_max (it compares the along-line h_f = c_f ell with flake_h_max / mu,
    FaceRemoval.H:536-538 and ReconFromX). So the cap can act on a PASSING line only where L / C_L(T_s) > 1.5 mm,
    i.e. above T_s*(L) with C_L(T_s*) = L / 1.5 mm. C_L is P2-15's (studies/p2_15_buckle/derive_span.py, Grimsel.H
    Table 3, nu 0.30 = rock_grimsel.inputs' grimsel.nu), imported. Beside it: the on-disk max firing T_surf (every
    event-log row, any direction) on P2-16's CAL and HOS_R.

z2  THE MINIMUM FLUX AT A SPAN (P2-16 DR-4). L_req(h) at each span-off CAL flux leg's EARLIEST firing is P2-15's
    discriminant as P2-16's rules.py imports it (rules.l_req); h_min(L) is log-log interpolated between the two
    legs bracketing L. The net flux at h_min is the code's own face balance (SurfaceCellFlux, MMWSpalling.H:3928,
    the surface branch at the reconstructed T_s):
        q_net = h (T_g - T_s) - eps sigma (T_s^4 - T_amb^4) - h_conv (T_s - T_amb)          [W/m2]
    with T_g 1673.15 K, T_amb 293.15 K, h_conv 10 W/m2K (the CAL deck) and eps 0.8 (the rock file); T_s is the
    earliest firing's T_surf, linear in ln h between the two bracketing legs. The packet's form (no h_conv term)
    is printed beside it. The CAL deck's own q_net (h 500 at its own earliest firing) is the reference.
"""
import importlib.machinery
import importlib.util
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
STUD = os.path.join(ROOT, "tests", "MMWSpalling", "studies")


def load(name, path):
    ld = importlib.machinery.SourceFileLoader(name, path)
    m = importlib.util.module_from_spec(importlib.util.spec_from_loader(name, ld))
    ld.exec_module(m)
    return m


RU = load("rules17", os.path.join(STUD, "p2_16_span_rule", "rules.py"))
DS = RU.R15.DS
CAP = 1.5e-3
T_G, T_AMB, H_CONV, EPS, SIG = 1673.15, 293.15, 10.0, 0.8, 5.670374419e-8
T_TAB_END = 273.15 + 1000.0                  # Grimsel.H strain table ends at 1000 degC; C_L is flat above it
P14L = os.path.join(STUD, "p2_14a_cref_face", "legs")
FLUX_LEGS = (("C0_h150", 150.0), ("CAL", 500.0), ("C0_h1000", 1000.0), ("CAL_h2000", 2000.0))
P16L = os.path.join(STUD, "p2_16_span_rule", "legs")


def t_star(L):
    """T_s* with C_L(T_s*) = L / CAP (C_L falls with T); None if C_L never reaches it below the table's end"""
    want = L / CAP
    if DS.C_L(T_TAB_END) > want:
        return None
    lo, hi = 400.0, T_TAB_END
    for _ in range(100):
        m = 0.5 * (lo + hi)
        lo, hi = (m, hi) if DS.C_L(m) > want else (lo, m)
    return 0.5 * (lo + hi)


def max_ts(leg):
    rows = RU.R15.events(os.path.join(P16L, leg))
    v = sorted(float(r["T_surf"]) for r in rows)
    return v[-1], v[int(0.99 * (len(v) - 1))], len(v)


def q_net(h, Ts, conv=True):
    return h * (T_G - Ts) - EPS * SIG * (Ts ** 4 - T_AMB ** 4) - (H_CONV * (Ts - T_AMB) if conv else 0.0)


def z1():
    L = ["z1 -- the cap-dead crossover T_s*(L): C_L(T_s*) = L / 1.5 mm (P2-15's C_L; table ends %.2f K, "
         "C_L(end) = %.3f)" % (T_TAB_END, DS.C_L(T_TAB_END)), "=" * 100]
    for span in (0.012, 0.0132, 0.01777):
        ts = t_star(span)
        L.append("   L %6.2f mm  L/1.5 mm = %.3f  T_s* = %s" % (1e3 * span, span / CAP,
                 "%.1f K" % ts if ts else "NONE below the table's end (the cap cannot act on a passing line at "
                                           "any T_s the law represents)"))
    for leg in ("CAL", "HOS_R"):
        mx, p99, n = max_ts(leg)
        L.append("   on disk (P2-16 legs/%s, %d firings): max T_surf %.1f K, p99 %.1f K" % (leg, n, mx, p99))
    L.append("z1 (report)")
    return "\n".join(L)


def z2():
    pts = []
    for leg, h in FLUX_LEGS:
        v, r = RU.l_req(P14L, leg)
        pts.append((h, v, float(r["T_surf"]), leg))
    L = ["z2 -- the minimum flux at a span (P2-16 DR-4): earliest-firing L_req on P2-14a's span-off CAL flux legs",
         "=" * 100]
    for h, v, Ts, leg in pts:
        L.append("   %-9s h %6.0f  T_s %.2f K  L_req %7.3f mm  q_net %.4f MW/m2" % (leg, h, Ts, 1e3 * v,
                                                                              1e-6 * q_net(h, Ts)))
    cal = next(p for p in pts if p[3] == "CAL")
    q_cal = q_net(cal[0], cal[2])
    for span in (0.0103, 0.012, 0.0132, 0.01777):
        br = [(a, b) for a, b in zip(pts[:-1], pts[1:]) if (a[1] - span) * (b[1] - span) <= 0]
        if not br:
            L.append("   L %6.2f mm: not bracketed by the flux legs" % (1e3 * span))
            continue
        (h1, v1, T1, n1), (h2, v2, T2, n2) = br[0]
        f = (math.log(span) - math.log(v1)) / (math.log(v2) - math.log(v1))
        hm = math.exp(math.log(h1) + f * (math.log(h2) - math.log(h1)))
        Ts = T1 + (math.log(hm) - math.log(h1)) / (math.log(h2) - math.log(h1)) * (T2 - T1)
        q, q0 = q_net(hm, Ts), q_net(hm, Ts, conv=False)
        L.append("   L %6.2f mm: between %s and %s  h_min %.1f W/m2K  T_s %.2f K  q_net,min %.4f MW/m2 (without "
                 "h_conv %.4f)  -> %s the CAL deck's %.4f MW/m2 (h 500): the h 500 deck %s with B on" % (
                     1e3 * span, n1, n2, hm, Ts, 1e-6 * q, 1e-6 * q0, "ABOVE" if q > q_cal else "BELOW",
                     1e-6 * q_cal, "would NOT spall" if q > q_cal else "would spall"))
    L.append("z2 (report)")
    return "\n".join(L)


def main():
    for fn, f in ((z1, "Z1.txt"), (z2, "Z2.txt")):
        t = fn()
        print(t)
        open(os.path.join(HERE, f), "w").write(t + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
