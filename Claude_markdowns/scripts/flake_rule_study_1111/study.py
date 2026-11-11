#!/usr/bin/env python3
"""Flake-rule study (2026-11-11, zero compute): what each candidate thickness rule reads, and where they part.

One thermal tool (an explicit 1-D column, dz 0.1 mm flame / 0.05 mm beam; rho c_p(T), lambda(T) as coded) supplies
profiles; every rule and every law variant is evaluated on the SAME stored profiles.

LAW VARIANTS (the stress law sigma = R E(T) [eps(T) - eps(T0)] / (1 - nu), eps = integral of alpha, held above 1000 C):
  V0  as coded (Grimsel.H Table 3: alpha 28.4e-6 and E/E0 0.50 held flat above 600-650 C)
  V1  alpha falls back after the quartz inversion: 28.3e-6 at 575 C -> 14e-6 at 600 C -> 10e-6 from 650 C
      (ILLUSTRATIVE: dilatometer curves spike near 570 C and fall to a lower beta-quartz value; no Grimsel curve exists)
  V2  E/E0 keeps falling: 0.50 at 600 C -> 0.25 at 900 C, held (ILLUSTRATIVE: cooled-specimen trend)
  V3  V1 and V2 together
V1-V3 differ from V0 only above 575 C (848 K), so the onset level and the R calibration are unchanged.

RULES (the flake's thickness x; every rule is given its OWN span from Kant's bracket, as a re-derivation would):
  TAN   0.4 x the tangent length k (T_s - T0) / q_net              (the layer rule; flame only)
  RISE  depth where the temperature rise falls to e^-0.4 of the surface rise
  ISO   depth of a FIXED isotherm: the temperature at which RISE sits at the onset, T* = T0 + e^-0.4 (846 K - T0)
  SIG   depth where the stress falls to 0.75 of its surface value
  WB    Rauenzahn & Tester's Weibull median: z with int_0^z z'^2 sigma^m dz' = 1/2 of the total, m = 13
  WBs   the slab Weibull median: int_0^z sigma^m dz' = 1/2 of the total, m = 13
Buckling: bk = C_L(T_s) x / L_span, C_L = pi / sqrt(6 (1 + nu) [eps(T_s) - eps(T0)]); a line fires iff bk <= 1.
"""
import math
import os
import pickle
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
E0, NU, SIGSB, R, M = 53.3e9, 0.30, 5.670374e-8, 0.531294, 13.0
TH0 = [20.0, 100.0, 200.0, 300.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0]
AL0 = [8.4, 8.7, 9.1, 9.5, 16.4, 19.9, 23.4, 26.9, 28.3, 28.4]
EE0 = [1.00, 1.00, 1.00, 0.93, 0.85, 0.75, 0.65, 0.58, 0.50, 0.50]
KI0 = [1.00, 1.00, 0.90, 0.82, 0.78, 0.77, 0.75, 0.65, 0.60, 0.60]


class Law:
    def __init__(self, name, al=None, ee=None):
        self.name = name
        self.tha, self.al = (TH0, AL0) if al is None else al
        self.the, self.ee = (TH0, EE0) if ee is None else ee
        self.tg = np.arange(20.0, 1000.0 + 1e-9, 0.25)
        a = np.interp(self.tg, self.tha, self.al) * 1e-6
        self.eps = np.concatenate([[0.0], np.cumsum(0.5 * (a[1:] + a[:-1]) * 0.25)])

    def strain(self, T):
        return np.interp(np.clip(np.asarray(T, float) - 273.15, 20.0, 1000.0), self.tg, self.eps)

    def E(self, T):
        return E0 * np.interp(np.asarray(T, float) - 273.15, self.the, self.ee)

    def sig(self, T, T0):
        return R * self.E(T) * (self.strain(T) - self.strain(T0)) / (1 - NU)

    def C_L(self, Ts, T0):
        return math.pi / math.sqrt(6.0 * (1 + NU) * float(self.strain(Ts) - self.strain(T0)))

    def score(self, T, T0, a=8.0e-6):
        kic = 1.5e6 * np.interp(np.asarray(T, float) - 273.15, TH0, KI0)
        return 1.1215 * self.sig(T, T0) * math.sqrt(math.pi * a) / kic


AL1 = ([20, 100, 200, 300, 400, 450, 500, 550, 575, 600, 650], [8.4, 8.7, 9.1, 9.5, 16.4, 19.9, 23.4, 26.9, 28.3, 14.0, 10.0])
EE2 = ([20, 100, 200, 300, 400, 450, 500, 550, 600, 900], [1.00, 1.00, 1.00, 0.93, 0.85, 0.75, 0.65, 0.58, 0.50, 0.25])
LAWS = [Law("V0"), Law("V1", al=AL1), Law("V2", ee=EE2), Law("V3", al=AL1, ee=EE2)]
RULES = ("TAN", "RISE", "ISO", "SIG", "WB", "WBs")


def cp(T):
    return 1370.0 - 178000.0 / (T - 273.15 + 271.0)


def first_cross(zz, v, target):
    m = v <= target
    if not m.any():
        return float("nan")
    i = int(np.argmax(m))
    if i == 0:
        return 0.0
    return zz[i - 1] + (zz[i] - zz[i - 1]) * (v[i - 1] - target) / (v[i - 1] - v[i])


def rules(law, z, T, Ts, T0, ell_tan):
    """the five thicknesses [m] on one profile (z cell centres, T cell values, Ts the face value)"""
    zz = np.concatenate([[0.0], z])
    TT = np.concatenate([[Ts], T])
    out = {"TAN": 0.4 * ell_tan if ell_tan and ell_tan > 0 else float("nan")}
    out["RISE"] = first_cross(zz, TT - T0, math.exp(-0.4) * (Ts - T0))
    out["ISO"] = first_cross(zz, TT - T0, math.exp(-0.4) * (846.0 - T0)) if Ts - T0 > math.exp(-0.4) * (846.0 - T0) else float("nan")
    sg = np.maximum(law.sig(TT, T0), 0.0)
    out["SIG"] = first_cross(zz, sg, 0.75 * sg[0]) if sg[0] > 0 else float("nan")
    w = (sg / sg.max()) ** M
    for key, wt in (("WB", zz ** 2 * w), ("WBs", w)):
        c = np.concatenate([[0.0], np.cumsum(0.5 * (wt[1:] + wt[:-1]) * np.diff(zz))])
        out[key] = float(np.interp(0.5 * c[-1], c, zz)) if c[-1] > 0 else float("nan")
    return out


def column(T0, t_max, lam_a, rho, levels, dz, depth, h=0.0, Tg=0.0, hc=0.0, Tamb=293.15, q=0.0, alpha=0.0, eps=0.8,
           T_init=None):
    lam = lambda T: lam_a / (T - 273.15 + 563.0)
    n = int(round(depth / dz))
    z = (np.arange(n) + 0.5) * dz
    ze = dz * np.arange(n + 1)
    S = q * (np.exp(-alpha * ze[:-1]) - np.exp(-alpha * ze[1:])) / dz if alpha > 0 else np.zeros(n)
    T = np.full(n, T0) if T_init is None else T_init(z)
    dt = 0.4 * dz * dz / (lam(293.15) / (rho * cp(293.15)))
    out, todo, t = [], list(levels), 0.0
    while t < t_max and todo:
        k = lam(T[0])
        Ts = T[0]
        for _ in range(30):
            g = h * (Tg - Ts) - hc * (Ts - Tamb) - eps * SIGSB * (Ts ** 4 - Tamb ** 4) - k * (Ts - T[0]) / (0.5 * dz)
            Ts -= g / (-h - hc - 4 * eps * SIGSB * Ts ** 3 - k / (0.5 * dz))
        while todo and Ts >= todo[0]:
            qn = h * (Tg - Ts) - hc * (Ts - Tamb) - eps * SIGSB * (Ts ** 4 - Tamb ** 4)
            out.append(dict(t=t, Ts=Ts, z=z.copy(), T=T.copy(), level=todo.pop(0), ell_tan=(lam(Ts) * (Ts - T0) / qn) if qn > 0 else None))
        L = lam(T)
        Lf = 2.0 * L[:-1] * L[1:] / (L[:-1] + L[1:])
        F = Lf * (T[1:] - T[:-1]) / dz
        div = np.zeros(n)
        div[:-1] += F
        div[1:] -= F
        div[0] += k * (Ts - T[0]) / (0.5 * dz)
        T = T + dt * (div / dz + S) / (rho * cp(T))
        t += dt
    return out


KT = {130: (747, 1250), 140: (690, 1175), 150: (633, 1100), 160: (576, 1025)}
ROCK = {"G": (1923.5, 2594.0), "B": (1540.6, 2566.0)}
KLEGS = ("KG130", "KG140", "KG150", "KG160", "KB150", "KB160")
BEAM = (("D0.5 Q2", 0.5, 2.0), ("D1 Q2", 1.0, 2.0), ("D2 Q2", 2.0, 2.0), ("D4 Q2", 4.0, 2.0), ("D8 Q2", 8.0, 2.0), ("D2 Q4", 2.0, 4.0))


def profiles():
    pk = os.path.join(HERE, "profiles.pkl")
    if os.path.exists(pk):
        return pickle.load(open(pk, "rb"))
    P = {}
    for leg in KLEGS:
        a, rho = ROCK[leg[1]]
        h, th = KT[int(leg[2:])]
        P[leg] = column(293.15, 200.0, a, rho, [820.0 + 10.0 * i for i in range(60)], 0.1e-3, 0.064, h=h, Tg=th + 273.15, hc=10.0)
    for name, d, qq in BEAM:
        P[name] = column(297.5, 80.0, 2000.0, 2650.0, [846.0 + 20.0 * i for i in range(29)] + [1410.0], 0.05e-3, 0.05,
                         q=qq * 1e6, alpha=1.0 / (d * 1e-3), hc=56.0, Tamb=297.5)
    for Tb in (293.15, 450.0, 550.0):          # the flat flame fixture on a base preheated uniformly to T_b
        P["PRE%d" % Tb] = column(293.15, 60.0, 2000.0, 2650.0, [846.0], 0.05e-3, 0.064, h=700.0, Tg=1866.0, hc=10.0,
                                 T_init=lambda z, Tb=Tb: np.full(z.shape, Tb))
    pickle.dump(P, open(pk, "wb"))
    return P


def main():
    p = print
    P = profiles()
    V0 = LAWS[0]
    p("A. THE LAW AS CODED, and the three variants above 575 C")
    p("   T K   theta C | E/E0   eps 1e-3  R sigma MPa  K_Ic   score  C_L   n = dln(sigma)/dln(dT) |  R sigma MPa under V1, V2, V3 | n under V1, V2, V3")
    for T in (500, 600, 700, 750, 800, 820, 846, 873, 900, 950, 1000, 1100, 1200, 1273, 1350, 1410):
        def n_of(law):
            return float((law.sig(T + 1.0, 293.15) - law.sig(T - 1.0, 293.15)) / 2.0 * (T - 293.15) / law.sig(T, 293.15))
        p("   %5d %6.0f | %5.2f %8.3f %10.1f %6.2f %7.3f %6.2f %8.2f                | %7.1f %7.1f %7.1f      | %5.2f %5.2f %5.2f" % (
            T, T - 273.15, float(V0.E(T)) / E0, 1e3 * float(V0.strain(T) - V0.strain(293.15)), float(V0.sig(T, 293.15)) / 1e6,
            float(np.interp(T - 273.15, TH0, KI0)), float(V0.score(T, 293.15)), V0.C_L(T, 293.15), n_of(V0),
            *[float(l.sig(T, 293.15)) / 1e6 for l in LAWS[1:]], *[n_of(l) for l in LAWS[1:]]))
    Ton = next(T for T in np.arange(700.0, 1000.0, 0.1) if float(V0.score(T, 293.15)) >= 1.0)
    p("   onset: score = 1 at %.1f K (the code fires at 845.5-846.0 K on the fixture)" % Ton)

    for law in LAWS:
        Tw = [T for T in np.arange(700.0, 1410.5, 1.0) if float(law.score(T, 293.15)) >= 1.0]
        segs, a = [], None
        for T in Tw:
            if a is None:
                a = b = T
            elif T - b > 1.5:
                segs.append((a, b)); a = b = T
            else:
                b = T
        if a is not None:
            segs.append((a, b))
        p("   %s: the onset criterion (score >= 1) holds for T_s in %s K; R sigma peaks at %.0f MPa" % (
            law.name, ", ".join("[%.0f, %.0f]" % x for x in segs), max(float(law.sig(T, 293.15)) for T in np.arange(700.0, 1411.0, 1.0)) / 1e6))
    p("\nB. EVERY RULE ON AN EXPONENTIAL SKIN: thickness in units of ell, by surface temperature (law V0 | V1 | V2 | V3)")
    zeta = np.linspace(1e-4, 8.0, 80000)
    for rule in RULES[1:]:
        row = []
        for Ts in (800.0, 846.0, 900.0, 1000.0, 1100.0, 1273.0, 1400.0):
            vals = [rules(l, zeta, 293.15 + (Ts - 293.15) * np.exp(-zeta), Ts, 293.15, None)[rule] for l in LAWS]
            row.append("%4.0f K: %s" % (Ts, " ".join("%.3f" % v for v in vals)))
        p("   %-5s %s" % (rule, " | ".join(row)))

    p("\nC. KANT'S BRACKET (flame; first firing at 820 K; history to 200 s). L_req = C_L(T_s) x [mm]")
    span = {}
    for law in LAWS:
        for rule in RULES:
            on, mn, where = {}, {}, {}
            for leg in KLEGS:
                Ls = []
                for s in P[leg]:
                    x = rules(law, s["z"], s["T"], s["Ts"], 293.15, s["ell_tan"])[rule]
                    Ls.append((law.C_L(s["Ts"], 293.15) * x, s["level"]))
                on[leg] = Ls[0][0]
                good = [(a, b) for a, b in Ls if math.isfinite(a)]
                mn[leg], where[leg] = min(good)
            g = (on["KG130"], min(on["KG140"], on["KG150"], on["KG160"]))
            b = (on["KB150"], on["KB160"])
            lo, hi = max(g[0], b[0]), min(g[1], b[1])
            span[(law.name, rule)] = math.sqrt(g[1] * b[0])               # the knife-edge centre: this rule's own span
            sp = span[(law.name, rule)]
            later = [leg for leg in ("KG140", "KG150", "KG160", "KB160") if mn[leg] < 0.98 * on[leg] and mn[leg] <= sp
                     and float(law.score(where[leg], 293.15)) >= 1.0]
            p("   %s %-5s first firing: KG130 %.2f | KG140 %.2f | KB150 %.2f | KB160 %.2f; one span %s; span %.2f mm; "
              "history minimum / first firing: %s; no-spall legs that would fire later: %s" % (
                  law.name, rule, 1e3 * on["KG130"], 1e3 * on["KG140"], 1e3 * on["KB150"], 1e3 * on["KB160"],
                  ("exists (%.2f %% wide)" % (100 * (hi / lo - 1))) if lo < hi else ("misses by %.2f %%" % (100 * (lo / hi - 1))),
                  1e3 * sp, " ".join("%s %.2f@%.0fK" % (leg[1:], mn[leg] / on[leg], where[leg]) for leg in ("KG140", "KG150", "KG160")),
                  ", ".join(later) if later else "none"))
        p("")

    p("D. THE BEAM COLUMN with each rule's own span: bk at the criterion (846 K), its minimum before the solidus, and the outcome")
    for law in LAWS:
        for rule in RULES[1:]:
            sp = span[(law.name, rule)]
            row = []
            for name, d, qq in BEAM:
                best, first, b0, x0 = 9e9, None, None, None
                for s in P[name]:
                    x = rules(law, s["z"], s["T"], s["Ts"], 297.5, None)[rule]
                    bk = law.C_L(s["Ts"], 297.5) * x / sp
                    if b0 is None:
                        b0, x0 = bk, x
                    if math.isfinite(bk) and bk < best:
                        best = bk
                    if first is None and bk <= 1.0 and float(law.score(s["Ts"], 297.5)) >= 1.0:
                        first = s["level"]
                row.append("%s: x %.2f mm bk %.2f, min %.2f -> %s" % (name, 1e3 * x0, b0, best, ("S@%.0fK" % first) if first else "MELT"))
            p("   %s %-5s (span %.2f mm) %s" % (law.name, rule, 1e3 * sp, " | ".join(row)))
        p("")

    p("E. A SKIN ON A PREHEATED BASE (the fixture's flame on rock already at T_b; the face at 846 K): thickness [mm], law V0")
    for Tb in (293.15, 450.0, 550.0):
        s = P["PRE%d" % Tb][0]
        r = rules(V0, s["z"], s["T"], s["Ts"], 293.15, s["ell_tan"])
        i2 = int(round(2e-3 / (s["z"][1] - s["z"][0])))
        zc = np.array([1e-3, 3e-3])
        Tc = np.array([s["T"][:2 * i2 // 2].mean(), s["T"][i2:2 * i2].mean()]) if True else None
        Tc = np.array([s["T"][:i2].mean(), s["T"][i2:2 * i2].mean()])
        two = 2e-3 / math.log((Tc[0] - 293.15) / (Tc[1] - 293.15))
        p("   T_b %5.0f K, t %5.2f s: TAN %.3f | two-cell at 2 mm (today's) %.3f | RISE %.3f | SIG %.3f | WB %.3f | WBs %.3f   [C_L %.2f: 12 mm span blocks above %.3f mm]" % (
            Tb, s["t"], 1e3 * r["TAN"], 0.4e3 * two, 1e3 * r["RISE"], 1e3 * r["SIG"], 1e3 * r["WB"], 1e3 * r["WBs"],
            V0.C_L(s["Ts"], 293.15), 12.0 / V0.C_L(s["Ts"], 293.15)))


if __name__ == "__main__":
    main()
