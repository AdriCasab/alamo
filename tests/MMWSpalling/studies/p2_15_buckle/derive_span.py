#!/usr/bin/env python3
"""P2-15 a5 -- derive the buckling span on Kant 2016's Gotthard bracket (spalls at SOD 130, silent at 140).

For each span-off leg KG130_0 / KG140_0 (Gotthard thermal, Grimsel strength): per column the first top-face
non-fallback firing; ell_eff and T_surf medians. C_L = pi / sqrt(6 (1 + nu) eps_th), eps_th = strain(T_s) -
strain(293.15 K) from Grimsel.H's Table 3 alpha (re-typed below, trapezoid on the code's 1 K grid; self-checked to
C_L = 12.9 +- 0.1 at 843.76 K). L_lo = C_L c_f ell(130), L_hi = C_L c_f ell(140), L_span = sqrt(L_lo L_hi) (mu = 1:
flat floor). Writes SPAN.txt (first token = the span in metres, 6 significant figures, no hand rounding).
The four span-off controls KG150_0, KG160_0, KB150_0, KB160_0 must fire on >= 90 % of the 576 columns.
Usage: derive_span.py   (exit 1 on a failed gate)
"""
import csv
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
NU, CF, T0 = 0.30, 0.4, 293.15
# Grimsel.H:77-85, Table 3: theta [degC], alpha [1e-6/K] (all ten rows)
TH = [20.0, 100.0, 200.0, 300.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0]
AL = [8.4, 8.7, 9.1, 9.5, 16.4, 19.9, 23.4, 26.9, 28.3, 28.4]


def alpha(th):
    if th <= TH[0]:
        return 1e-6 * AL[0]
    if th >= TH[-1]:
        return 1e-6 * AL[-1]
    for i in range(len(TH) - 1):
        if th <= TH[i + 1]:
            f = (th - TH[i]) / (TH[i + 1] - TH[i])
            return 1e-6 * (AL[i] + f * (AL[i + 1] - AL[i]))


TAB = [0.0]
for i in range(1, 981):                      # 20 .. 1000 degC at 1 K, as Grimsel.H builds strain_tab
    TAB.append(TAB[-1] + 0.5 * (alpha(20.0 + i - 1) + alpha(20.0 + i)))


def strain(T):
    th = T - 273.15
    if th <= 20.0:
        return 0.0
    i = int(th - 20.0)
    if i >= 980:
        return TAB[980]
    thi = 20.0 + i
    return TAB[i] + 0.5 * (alpha(thi) + alpha(th)) * (th - thi)


def C_L(T_s):
    return math.pi / math.sqrt(6.0 * (1.0 + NU) * (strain(T_s) - strain(T0)))


def first_top(leg):
    f = {}
    for r in csv.DictReader(open(os.path.join(HERE, "legs", leg, "run_removal_events.csv"))):
        if r["dir"] == "5" and float(r["ts_fb"]) == 0.0:
            f.setdefault((r["col_i"], r["col_j"]), r)
    return f


def main():
    bad = []
    L = ["a5 -- the buckling span, derived on Kant 2016's Gotthard bracket (130 spalls, 140 silent)", "=" * 100]
    c0 = C_L(843.76)
    ok = abs(c0 - 12.9) <= 0.1
    L.append("  self-check: eps_th(843.76 K) = %.5e, C_L = %.3f (bar 12.9 +- 0.1): %s" % (
        strain(843.76) - strain(T0), c0, "ok" if ok else "FAIL"))
    if not ok:
        bad.append("C_L self-check")
    Ls = {}
    for leg in ("KG130_0", "KG140_0"):
        f = first_top(leg)
        ell = np.array([float(r["ell_eff"]) for r in f.values()])
        Ts = np.array([float(r["T_surf"]) for r in f.values()])
        cl = C_L(float(np.median(Ts)))
        Ls[leg] = cl * CF * float(np.median(ell))
        L.append("  %-8s n %3d  ell_eff median %.4f mm (p10 %.4f, p90 %.4f)  T_surf median %.2f K  C_L %.3f  "
                 "C_L c_f ell = %.5f m" % (leg, len(f), 1e3 * np.median(ell), 1e3 * np.percentile(ell, 10),
                                           1e3 * np.percentile(ell, 90), np.median(Ts), cl, Ls[leg]))
    lo, hi = Ls["KG130_0"], Ls["KG140_0"]
    sep = lo < hi
    span = math.sqrt(lo * hi)
    L.append("  L_lo %.5f m < L_hi %.5f m: %s;  L_span = sqrt(L_lo L_hi) = %.6g m" % (lo, hi, sep, span))
    if not sep:
        bad.append("no span separates the bracket")
    for leg in ("KG150_0", "KG160_0", "KB150_0", "KB160_0"):
        n = len(first_top(leg))
        okc = n >= 0.9 * 576
        L.append("  control %-8s fires on %3d/576 columns (%.1f %%): %s" % (leg, n, 100.0 * n / 576, "ok" if okc else "SILENT"))
        if not okc:
            bad.append("control " + leg)
    L.append("a5 %s (gate)%s" % ("PASS" if not bad else "FAIL", "" if not bad else ": " + "; ".join(bad)))
    txt = "\n".join(L)
    print(txt)
    open(os.path.join(HERE, "A5.txt"), "w").write(txt + "\n")
    if not bad:
        open(os.path.join(HERE, "SPAN.txt"), "w").write("%.6g  m  (L_lo %.6g, L_hi %.6g; derive_span.py, A5.txt)\n"
                                                        % (span, lo, hi))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
