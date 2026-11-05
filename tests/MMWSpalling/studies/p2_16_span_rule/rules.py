#!/usr/bin/env python3
"""P2-16 c0 + c0b (zero compute) -- does the span follow a measured rock property?

The windows are P2-15's deck discriminant, IMPORTED (studies/p2_15_buckle/reduce.py: events, DS.C_L, CF): the
span-off twin's EARLIEST firing decides the flat deck (all-or-nothing, P2-15 DR-5), so the span a leg requires is
    L_req [m] = C_L(T_s) [-] * mu [-] * c_f [-] * ell_eff [m]          (mu = 1 on the still-flat floor; ledger L4)
Gotthard window [L(KG130_0), min L(KG140/150/160_0)); Bethel window [L(KB150_0), L(KB160_0)).

Rules (pre-registered, Delta = 0 only). Predicted Bethel/Gotthard span ratio rho:
  R0 one span for all rocks (grain size; "feldspar up to ~4 mm" in all three)   rho = 1
  R1 L ~ 1/UCS                                                                   rho = 214/152
  R2 L ~ 1/UCS^2 (Griffith length at equal toughness)                            rho = (214/152)^2
A rule's feasible Gotthard window is [G_lo, G_hi) n [B_lo/rho, B_hi/rho). Relative width = (hi - lo) / L_G*, with
L_G* = sqrt(lo hi) the window's geometric mean. SURVIVES iff non-empty, relative width >= 2 %, and L_G* >= c0b's
bound. Grimsel's span: L_Gr = L_G* f(Grimsel)/f(Gotthard), UCS 169 +- 37 MPa carried as a band.

c0b (one-sided; no silent case): L >= C_L(T_s) c_f k_r dT / Q at R&T 1989's lowest Q_max 0.68 MW/m2, dT 465 K,
k_r 2.14 W/mK (Wilkinson 1993), T_s = 293.15 + 465 K; and with the deck's own k_eff 1.6-2.1 (P2-15 DR-12).
⚠ The rules' Barre-scaled level would need Barre's UCS, which no source on disk gives: the bound is compared with
each rule's level at the Kant/Grimsel rocks, as P1b did, and is non-discriminating by construction (one-sided).

Usage: rules.py [legs-dir]   (default: P2-15's legs; c1 passes its own) -> RULES.txt
"""
import importlib.machinery
import importlib.util
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
P15 = os.path.join(ROOT, "tests", "MMWSpalling", "studies", "p2_15_buckle")
ld = importlib.machinery.SourceFileLoader("red15", os.path.join(P15, "reduce.py"))
R15 = importlib.util.module_from_spec(importlib.util.spec_from_loader("red15", ld))
ld.exec_module(R15)

UCS = {"Bethel": 152.0, "Gotthard": 214.0, "Grimsel": 169.0}      # MPa; Kant 2016 Table 2; rock_grimsel.inputs:59
UCS_GR_ERR = 37.0
RULES = (("R0", "one span for all rocks", 0), ("R1", "L ~ 1/UCS", 1), ("R2", "L ~ 1/UCS^2", 2))
MIN_REL = 0.02


def l_req(legdir, leg):
    """L_req [m] at the leg's EARLIEST firing (the first event-log row), P2-15's discriminant verbatim"""
    rows = R15.events(os.path.join(legdir, leg))
    if not rows:
        return None, None
    r = rows[0]
    return R15.DS.C_L(float(r["T_surf"])) * R15.CF * float(r["ell_eff"]) * float(r["mu"]), r


def windows(legdir=os.path.join(P15, "legs"), bethel=("KB150_0", "KB160_0")):
    L = {}
    for leg in ("KG130_0", "KG140_0", "KG150_0", "KG160_0") + tuple(bethel):
        L[leg] = l_req(legdir if leg in bethel else os.path.join(P15, "legs"), leg)[0]
    g = (L["KG130_0"], min(L["KG140_0"], L["KG150_0"], L["KG160_0"]))
    b = (L[bethel[0]], L[bethel[1]])
    return g, b, L


def score(g, b):
    out = []
    for tag, name, n in RULES:
        rho = (UCS["Gotthard"] / UCS["Bethel"]) ** n
        lo, hi = max(g[0], b[0] / rho), min(g[1], b[1] / rho)
        ok = lo < hi
        gm = math.sqrt(lo * hi) if ok else float("nan")
        rel = (hi - lo) / gm if ok else float("nan")
        out.append(dict(tag=tag, name=name, n=n, rho=rho, lo=lo, hi=hi, ok=ok, gm=gm, rel=rel))
    return out


def c0b():
    T_s = 293.15 + 465.0
    cl = R15.DS.C_L(T_s)
    L = lambda k: cl * R15.CF * k * 465.0 / 0.68e6
    return T_s, cl, L(2.14), L(1.6), L(2.1)


def grimsel(gm, n):
    f = lambda u: (UCS["Gotthard"] / u) ** n
    return gm * f(UCS["Grimsel"]), gm * f(UCS["Grimsel"] + UCS_GR_ERR), gm * f(UCS["Grimsel"] - UCS_GR_ERR)


def report(g, b, L, label):
    T_s, cl, b214, b16, b21 = c0b()
    out = ["%s -- the span rules on the windows (mm; earliest-firing L_req, P2-15's discriminant imported)" % label,
           "=" * 100]
    for k, v in L.items():
        out.append("   L_req %-9s %8.3f mm" % (k, 1e3 * v))
    out.append("   Gotthard window [%.3f, %.3f) mm   Bethel window [%.3f, %.3f) mm" % tuple(1e3 * x for x in g + b))
    out.append("   feasible rho = (B_lo/G_hi, B_hi/G_lo) = (%.4f, %.4f)" % (b[0] / g[1], b[1] / g[0]))
    rows = score(g, b)
    for s in rows:
        surv = s["ok"] and s["rel"] >= MIN_REL and s["gm"] >= b214
        if s["ok"]:
            gr, gr_hi_ucs, gr_lo_ucs = grimsel(s["gm"], s["n"])
            e_lo = min(grimsel(s["lo"], s["n"])[1:]); e_hi = max(grimsel(s["hi"], s["n"])[1:])
            out.append("   %s %-22s rho %.4f  L_G in [%.3f, %.3f) mm, width %.2f %% of L_G* %.3f mm  -> Bethel %.3f mm; "
                       "Grimsel %.2f mm (UCS 169+-37 at L_G*: %.2f-%.2f mm; with the window's edges too: %.2f-%.2f mm)  %s" % (
                           s["tag"], s["name"], s["rho"], 1e3 * s["lo"], 1e3 * s["hi"], 100 * s["rel"], 1e3 * s["gm"],
                           1e3 * s["gm"] * s["rho"], 1e3 * gr, 1e3 * min(gr_hi_ucs, gr_lo_ucs),
                           1e3 * max(gr_hi_ucs, gr_lo_ucs), 1e3 * e_lo, 1e3 * e_hi, "SURVIVES" if surv else
                           ("FAILS (width < 2 %: a coincidence, not support)" if s["rel"] < MIN_REL else "FAILS (c0b)")))
        else:
            out.append("   %s %-22s rho %.4f  L_G window EMPTY (needs [%.3f, %.3f) mm)  -> EXCLUDED" % (
                s["tag"], s["name"], s["rho"], 1e3 * s["lo"], 1e3 * s["hi"]))
        s["surv"] = surv
    out.append("c0b Barre/Westerly one-sided bound: C_L(%.2f K) = %.3f; L >= %.2f mm at k_r 2.14 (%.2f-%.2f mm at the "
               "deck's k_eff 1.6-2.1); every rule's L_G* is %s it (report; one-sided, non-discriminating)" % (
                   T_s, cl, 1e3 * b214, 1e3 * b16, 1e3 * b21,
                   "above" if all(s["gm"] >= b214 for s in rows if s["ok"]) else "NOT all above"))
    return "\n".join(out), rows


def main(a):
    g, b, L = windows()
    txt, _ = report(g, b, L, "c0 (Delta = 0, P2-15's accepted legs)")
    print(txt)
    open(os.path.join(HERE, "RULES.txt"), "w").write(txt + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
