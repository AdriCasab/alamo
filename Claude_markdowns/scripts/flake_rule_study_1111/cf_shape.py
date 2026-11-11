#!/usr/bin/env python3
"""Does the flake constant c matter beyond roughness? RISE and ISO with c = 0.4 / 0.2 / 0.1 (the flake base at
e^-c of the surface rise), each with the span Kant's bracket gives it (span/c is what Kant fixes: 12.02 mm at 0.4).
Flame: thickness at the fixture's first firing. Beam column: bk at the criterion and its minimum before a cut-off
temperature, law V0 and V1 (study.py). Kant 2016: theta_SP ~ 600 C, theta_FP ~ 550 C => c ~ 0.09 (0.19 at 100 K)."""
import math, os, pickle, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import study as S
P = S.profiles()
def iso_depth(s, T0, Tstar):
    zz = np.concatenate([[0.0], s["z"]]); TT = np.concatenate([[s["Ts"]], s["T"]])
    return S.first_cross(zz, TT, Tstar) if s["Ts"] > Tstar else float("nan")
print("flake base temperature at the onset (846 K surface, T0 293 K): c 0.4 -> %.0f K, 0.2 -> %.0f K, 0.1 -> %.0f K; Kant measured ~823 K under ~873 K" % tuple(293.15 + 552.85 * math.exp(-c) for c in (0.4, 0.2, 0.1)))
s = P["PRE293"][0]
print("fixture, first firing: thickness %s mm for c 0.4 / 0.2 / 0.1" % " / ".join("%.3f" % (1e3 * iso_depth(s, 293.15, 293.15 + math.exp(-c) * (s["Ts"] - 293.15))) for c in (0.4, 0.2, 0.1)))
for law in (S.LAWS[0], S.LAWS[1]):
    for rule in ("RISE", "ISO"):
        for c in (0.4, 0.2, 0.1):
            span = 12.02e-3 * c / 0.4
            row = []
            for name, d, qq in S.BEAM:
                b0, x0, best, first = None, None, 9e9, None
                for s in P[name]:
                    Tstar = 297.5 + math.exp(-c) * ((s["Ts"] if rule == "RISE" else 846.0) - 297.5)
                    x = iso_depth(s, 297.5, Tstar)
                    bk = law.C_L(s["Ts"], 297.5) * x / span
                    if b0 is None: b0, x0 = bk, x
                    if math.isfinite(bk) and bk < best: best = bk
                    if first is None and bk <= 1.0 and float(law.score(s["Ts"], 297.5)) >= 1.0: first = s["level"]
                row.append("%s x %.2f bk %.2f min %.2f %s" % (name, 1e3 * x0, b0, best, ("S@%.0f" % first) if first else "M"))
            print("%s %-4s c %.1f span %5.2f mm | %s" % (law.name, rule, c, 1e3 * span, " | ".join(row)))
