#!/usr/bin/env python3
"""Offline regime map, step 2: outcomes on the stored temperature histories (thermal.py).

DECISIONS APPLIED (the user, 2026-11-11): the flake's base is a FIXED temperature, the one it has at the line's onset,
T* = T0 + e^-c (T_on - T0); c re-anchored on Kant & von Rohr 2016 (surface ~600 C before, failure plane ~550 C after a
spall: c = ln(580/530) = 0.090; the single-spall 100 K difference gives 0.19, the bracket); the law above 600 C gets a
tail built from other granites; the span is re-derived on Kant's bracket for each c.

TAILS (all leave the law below 575 C, hence the onset and R, untouched):
  coded  Table 3 held flat above 600-650 C; strain frozen above 1000 C; fires up to the solidus (1410.5 K)
  T1     expansion falls to 7e-6 /K above 600 C (Alzahrani et al. 2022, eight granites, dilatometer: strain +0.14 to
         +0.29 % from 700 to 1000 C => 5-10e-6 /K); E/E0 held at 0.50; brittle limit 1000 C (Violay et al. 2017:
         Westerly brittle to 900 C, ductile at 1000 C, tested hot under confinement)
  T2     T1's expansion; E/E0 falls to 0.30 by 650 C (Zhang et al. 2022, tested hot, unconfined: 0.29 at 600 and 700 C)
  T3     T1's expansion; E/E0 falls linearly from 0.50 at 600 C to 0.10 at 1000 C
Outcome per cell: S spall at the onset; S@T spall later, at surface temperature T; M the line is still withdrawn at the
brittle limit / solidus; C the onset is not reached in 600 s; U onset reached, neither fired nor at the limit by 600 s.
"""
import math, os, pickle, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "flake_rule_study_1111"))
import study as S
import thermal as TH

AL_T = ([20, 100, 200, 300, 400, 450, 500, 550, 600, 650, 1200], [8.4, 8.7, 9.1, 9.5, 16.4, 19.9, 23.4, 26.9, 28.3, 7.0, 7.0])
EE_2 = ([20, 100, 200, 300, 400, 450, 500, 550, 600, 650], [1.00, 1.00, 1.00, 0.93, 0.85, 0.75, 0.65, 0.58, 0.50, 0.30])
EE_3 = ([20, 100, 200, 300, 400, 450, 500, 550, 600, 1000], [1.00, 1.00, 1.00, 0.93, 0.85, 0.75, 0.65, 0.58, 0.50, 0.10])


class LawT(S.Law):
    """a tail law: the strain is NOT frozen at 1000 C (it runs on at the tail's expansion)"""
    def __init__(self, name, al, ee=None):
        S.Law.__init__(self, name, al=al, ee=ee)
        self.tg = np.arange(20.0, 1200.0 + 1e-9, 0.25)
        a = np.interp(self.tg, self.tha, self.al) * 1e-6
        self.eps = np.concatenate([[0.0], np.cumsum(0.5 * (a[1:] + a[:-1]) * 0.25)])
    def strain(self, T):
        return np.interp(np.clip(np.asarray(T, float) - 273.15, 20.0, 1200.0), self.tg, self.eps)


LAWS = {"coded": (S.Law("coded"), 1410.5), "T1": (LawT("T1", AL_T), 1273.15), "T2": (LawT("T2", AL_T, EE_2), 1273.15),
        "T3": (LawT("T3", AL_T, EE_3), 1273.15)}
T0B, TON = 297.5, 843.0


def iso(s, Tstar):
    zz = np.concatenate([[0.0], s["z"]]); TT = np.concatenate([[s["Ts"]], s["T"]])
    return S.first_cross(zz, TT, Tstar) if s["Ts"] > Tstar else float("nan")


def kant_span(c):
    P = S.profiles()
    L = {}
    for leg in ("KG140", "KB150"):
        s = P[leg][0]
        L[leg] = LAWS["coded"][0].C_L(s["Ts"], 293.15) * iso(s, 293.15 + math.exp(-c) * (s["Ts"] - 293.15))
    return math.sqrt(L["KG140"] * L["KB150"])


def outcome(hist, law, cut, c, span, rule):
    if not hist:
        return "C", float("nan"), float("nan")
    b0 = x0 = None
    for s in hist:
        if s["level"] > cut:
            break
        Tstar = T0B + math.exp(-c) * ((TON if rule == "ISO" else s["Ts"]) - T0B)
        x = iso(s, Tstar)
        bk = law.C_L(s["Ts"], T0B) * x / span
        if b0 is None:
            b0, x0 = bk, x
        # the first snapshot IS the onset (843 K; the criterion reads 1 at 842.7-844 K depending on T0): tolerance 0.5 %
        if bk <= 1.0 and float(law.score(s["Ts"], T0B)) >= 0.995:
            return ("S" if s["level"] <= TON + 1 else "S@%d" % s["level"]), b0, x0
    reached = hist[-1]["level"]
    return ("M" if reached >= min(cut, 1410.0) - 25.0 else "U"), b0, x0


def main():
    p = print
    H = pickle.load(open(os.path.join(HERE, "thermal.pkl"), "rb"))
    p("LAW CHECK: R sigma [MPa] and strain [1e-3] at 600 / 700 / 800 / 1000 C")
    for k, (law, cut) in LAWS.items():
        p("   %-6s sigma %s | strain %s | score >= 1 from 843 K up to %s; fires below %.0f K" % (
            k, " / ".join("%.0f" % (float(law.sig(T, 293.15)) / 1e6) for T in (873.15, 973.15, 1073.15, 1273.15)),
            " / ".join("%.2f" % (1e3 * float(law.strain(T) - law.strain(293.15))) for T in (873.15, 973.15, 1073.15, 1273.15)),
            next(("%.0f K" % T for T in np.arange(850.0, 1411.0, 1.0) if float(law.score(T, 293.15)) < 1.0), "the limit"), cut))
    p("   measured total strain, eight Egyptian granites (Alzahrani 2022, six without the outliers): 0.87-1.08 %% at 500 C, 1.32-1.68 %% at 700 C, 1.58-1.91 %% at 1000 C (unconfined, includes crack opening)")
    spans = {c: kant_span(c) for c in (0.4, 0.19, 0.09)}
    p("\nKANT: the span each constant needs (first firings, the knife-edge centre): " + "; ".join("c %.2f -> %.2f mm" % (c, 1e3 * v) for c, v in spans.items()))
    p("flake base temperature under an 843 K surface: " + "; ".join("c %.2f -> %.0f K" % (c, 293.15 + math.exp(-c) * (843 - 293.15)) for c in spans))
    cases = [("TODAY'S FORM, for reference: c 0.40, base at a fraction of the rise, the coded law (span %.2f mm)" % (1e3 * spans[0.4]), "coded", 0.4, "RISE")]
    for c in (0.09, 0.19):
        for k in ("T1", "T2", "T3", "coded"):
            cases.append(("DECIDED FORM: c %.2f, fixed base temperature, tail %s (span %.2f mm)" % (c, k, 1e3 * spans[c]), k, c, "ISO"))
    grids = {}
    for title, k, c, rule in cases:
        law, cut = LAWS[k]
        g = {}
        for dk, D in TH.DS:
            for qk, q in TH.QS:
                g[(dk, qk)] = outcome(H["D%s_Q%s" % (dk, qk)], law, cut, c, spans[c], rule)
        grids[(k, c, rule)] = g
    def show(title, g, what):
        p("\n" + title)
        p("   absorbed flux MW/m2:  " + "".join("%9s" % q for _, q in TH.QS))
        for dk, D in TH.DS:
            p("   D %5.1f mm            " % (1e3 * D) + "".join("%9s" % (g[(dk, qk)][0] if what == "o" else ("%.2f" % g[(dk, qk)][1] if math.isfinite(g[(dk, qk)][1]) else "-")) for qk, _ in TH.QS))
    show(cases[0][0], grids[("coded", 0.4, "RISE")], "o")
    for c in (0.09, 0.19):
        same = all(grids[("T1", c, "ISO")][k][0] == grids[(t, c, "ISO")][k][0] for t in ("T2", "T3", "coded") for k in grids[("T1", c, "ISO")])
        show("DECIDED FORM, c %.2f, tail T1 (span %.2f mm) -- outcome" % (c, 1e3 * spans[c]), grids[("T1", c, "ISO")], "o")
        p("   the same outcome in every cell under T2, T3 and the coded law: %s" % ("YES" if same else "NO"))
        if not same:
            for t in ("T2", "T3", "coded"):
                d = [("D%s_Q%s" % k, grids[("T1", c, "ISO")][k][0], grids[(t, c, "ISO")][k][0]) for k in grids[("T1", c, "ISO")] if grids[("T1", c, "ISO")][k][0] != grids[(t, c, "ISO")][k][0]]
                p("      %s differs in: %s" % (t, d))
        show("   bk at the onset (a line spalls iff <= 1), c %.2f" % c, grids[("T1", c, "ISO")], "b")
    p("\nflake thickness at the onset [mm], c 0.09 (the depth of the fixed-temperature base)")
    p("   absorbed flux MW/m2:  " + "".join("%9s" % q for _, q in TH.QS))
    for dk, D in TH.DS:
        p("   D %5.1f mm            " % (1e3 * D) + "".join("%9s" % ("%.3f" % (1e3 * grids[("T1", 0.09, "ISO")][(dk, qk)][2]) if math.isfinite(grids[("T1", 0.09, "ISO")][(dk, qk)][2]) else "-") for qk, _ in TH.QS))


if __name__ == "__main__":
    main()
