#!/usr/bin/env python3
"""P2-11 (a0) -- the LOCAL Goldstein profile, checked against the AREA-AVERAGE it is derived from.

Goldstein, Behbahani & Heppelmann 1986 give the average over a disc of radius r:
    g_avg(r) = 533 / (533 + 44 (r/D)^1.285)            (= 1 at r = 0)
The local value is (1/2r) d(r^2 g_avg)/dr, which the packet writes in closed form as
    g_loc(r) = 533 (533 + 0.3575 u) / (533 + u)^2,     u = 44 (r/D)^1.285,  0.3575 = 1 - 1.285/2.
Checks: (i) the disc mean (2/R^2) int_0^R r g_loc(r) dr equals g_avg(R) at R/D = 0.5 ... 32 (rel <= 1e-6),
by quadrature, independent of the closed form's derivation; (ii) the closed form equals a numerical
derivative of r^2 g_avg; (iii) g_loc(17.5 mm) = 0.716 +- 0.001. Writes A0.txt beside the study.
"""
import math
import os

from scipy.integrate import quad

HERE = os.path.dirname(os.path.abspath(__file__))
D = 0.0079


def g_avg(r):
    return 533.0 / (533.0 + 44.0 * (r / D) ** 1.285)


def g_loc(r):
    u = 44.0 * (r / D) ** 1.285
    return 533.0 * (533.0 + 0.3575 * u) / (533.0 + u) ** 2


def main():
    L = ["P2-11 (a0) -- g_loc = (1/2r) d(r^2 g_avg)/dr; its disc mean must return Goldstein's average", "=" * 96]
    ok = True
    L.append("   R/D   disc mean of g_loc   g_avg(R)        rel       | closed form vs d(r^2 g_avg)/dr / 2r")
    for x in (0.5, 1, 2, 4, 8, 16, 32):
        R = x * D
        # substitution r = R t^(1/0.285)-free: integrate in t = r/R; the integrand is smooth except the
        # r^0.285 cusp at 0, which quad handles with points near 0 and a tight tolerance
        val, err = quad(lambda t: t * g_loc(t * R), 0.0, 1.0, epsabs=0, epsrel=1e-13, limit=500,
                        points=[1e-6, 1e-3])
        mean = 2.0 * val
        rel = abs(mean / g_avg(R) - 1)
        ok &= rel <= 1e-6
        h = R * 1e-5
        num = ((R + h) ** 2 * g_avg(R + h) - (R - h) ** 2 * g_avg(R - h)) / (2 * h) / (2 * R)
        rel2 = abs(g_loc(R) / num - 1)
        ok &= rel2 <= 1e-6
        L.append(f"  {x:5.1f}   {mean:.12f}       {g_avg(R):.12f}  {rel:.2e}   | {g_loc(R):.10f} vs {num:.10f}"
                 f"  rel {rel2:.2e}")
    gl = g_loc(0.0175)
    ok17 = abs(gl - 0.716) <= 0.001
    L.append(f"  g_loc(0) = {g_loc(0.0):.6f} (the stagnation value unchanged)")
    L.append(f"  g_loc(17.5 mm) = {gl:.5f} (bar 0.716 +- 0.001: {'ok' if ok17 else 'OUT'}); g_avg(17.5 mm) "
             f"= {g_avg(0.0175):.5f}; her crater 0.7095-0.7632")
    L.append(f"a0 {'PASS' if ok and ok17 else 'FAIL'}")
    txt = "\n".join(L)
    print(txt)
    open(os.path.join(os.path.dirname(HERE), "A0.txt"), "w").write(txt + "\n")


if __name__ == "__main__":
    main()
