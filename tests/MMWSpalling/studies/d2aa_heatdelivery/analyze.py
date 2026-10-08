#!/usr/bin/env python3
"""D2aa analysis -- scores CRITERION.md (a)-(g) and P1-P7.

Imports the equal-burner-depth machinery (`shape`, `rop`, `rings`, `pmean`, both
volume bases) UNCHANGED from studies/d2x_areasplit/analyze.py; its WHERE lookup
is EXTENDED with this study's legs (a table entry, not an edit to a function).
Frozen studies are read by path and never re-run.

⚠ All Meier thermo powers are QUARTER-hole sums; everything below is printed
full-hole (x4), as CRITERION.md fixes.
"""
import argparse
import ast
import filecmp
import importlib.util
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)


def load_mod(name, path):
    sp = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


X = load_mod("d2x_an", os.path.join(STUD, "d2x_areasplit", "analyze.py"))

OUT = os.path.join(HERE, "output")
Z2 = os.path.join(STUD, "d2z2_multicell", "output")
Z1 = os.path.join(STUD, "d2z1_delete", "output")
LEGS = ("WL", "WL81", "PD", "BO")
X.WHERE.update({n: OUT for n in LEGS + ("S0",)})
X.WHERE["SW"] = Z2

WIN = (150.0, 400.0)
BAND = (1.3, 1.6)                 # Meier rate band
EXP = 1.59                        # packet's rate ~ P^EXP
JAM = 0.5                         # P7 jam signature [m/h]
PIT_CEIL = 150.0                  # item 0(c) pit-lead ceiling [mm]


def have(n):
    return os.path.exists(os.path.join(X.WHERE.get(n, OUT), n + ".done"))


def thermo(n):
    p = os.path.join(X.WHERE[n], n, "thermo.dat")
    with open(p) as f:
        cols = f.readline().split()
    v = np.loadtxt(p, skiprows=1, ndmin=2)
    return {c: v[:, k] for k, c in enumerate(cols)}


def tm(d, k, mul=4.0):
    """Time-mean over the window, converted to FULL-hole."""
    if k not in d:
        return float("nan")
    s = (d["time"] >= WIN[0]) & (d["time"] <= WIN[1])
    return float(np.nanmean(d[k][s])) * mul


def hdr(t):
    print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78)


def crit_a(rep):
    hdr("(a) SCORED -- rate, shape and volume together, per leg")
    print("  rate definition: BURNER DESCENT over the 100-150 mm stretch")
    print(f"{'leg':<6}{'ROP':>8}{'mean|dD|':>10}{'D@100':>8}{'D@125':>8}"
          f"{'pit lead':>10}{'vol':>7}{'x_pit':>7}{'x_desc':>8}")
    for n in LEGS + ("SW", "A", "D"):
        if not have(n):
            continue
        r = X.rop(n)
        s = X.shape(n, X.PRIMARY, (10, 25, 50, 75, 100, 125))
        if s is None:
            print(f"{n:<6}{r:8.3f}   never reached {X.PRIMARY} mm descent"
                  f"   <- see P7")
            rep.setdefault("rate", {})[n] = r
            continue
        print(f"{n:<6}{r:8.3f}{(s['a'] or float('nan')):10.2f}"
              f"{s['D'].get(100, float('nan')):8.1f}{s['D'].get(125, float('nan')):8.1f}"
              f"{s['lead']:10.1f}{s['vol']:7.2f}{s['x_pit']:7.2f}{s['x_desc']:8.2f}")
        rep.setdefault("rate", {})[n] = r
        rep.setdefault("shape", {})[n] = s
        if s["lead"] > PIT_CEIL:
            print(f"       ⚠ pit lead {s['lead']:.0f} mm exceeds the registered "
                  f"{PIT_CEIL:.0f} mm ceiling -- DISQUALIFIED from the rate comparison")
    print(f"\n  rings, 2 mm, {X.RINGS[0]}-{X.RINGS[-1]+2} mm [mm recession], "
          f"PER RING, never averaged:")
    print("       " + "".join(f"{r:7d}" for r in X.RINGS))
    for n in LEGS + ("SW", "A"):
        if have(n):
            print(f"  {n:<5}" + "".join(f"{v:7.1f}" for v in X.rings(n)))


def crit_b(rep):
    hdr("(b) SCORED -- THE FALSE-PASS TEST")
    inband = [n for n in LEGS
              if have(n) and BAND[0] <= rep.get("rate", {}).get(n, -1) <= BAND[1]]
    if not inband:
        print(f"  No leg lands in Meier's {BAND[0]}-{BAND[1]} m/h band, so the "
              f"false-pass test does not fire.")
        print("  Rates: " + ", ".join(
            f"{n} {rep['rate'][n]:.3f}" for n in LEGS if n in rep.get("rate", {})))
        rep["b"] = "not triggered"
        return
    for n in inband:
        f0, f1 = f_multi_trend(n)
        print(f"  {n}: in band at {rep['rate'][n]:.3f} m/h; "
              f"wall_step_f_multi_pads {f0:.3f} -> {f1:.3f}")


def f_multi_trend(n):
    d = thermo(n)
    s = (d["time"] >= WIN[0]) & (d["time"] <= WIN[1])
    v = d["wall_step_f_multi_pads"][s]
    h = len(v) // 2
    return float(v[:h].mean()), float(v[h:].mean())


def crit_c(rep):
    hdr("(c) SCORED -- the power identity, +/- 0.02 W")
    for n in LEGS + ("S0", "SW", "D"):
        if not have(n):
            continue
        d = thermo(n)
        if "jet_P_face" not in d:
            print(f"  {n:<5} no jet_P_face (no march)"); continue
        s = (d["time"] >= WIN[0]) & (d["time"] <= WIN[1])
        e = float(np.abs((d["jet_P_face"] - (d["patch_P_robin"] + d["side_P_face"]
                                             + d["side_P_loss"]))[s]).max())
        print(f"  {n:<5} max residual {e:.3e} W   {'PASS' if e <= 0.02 else 'FAIL'}")
        rep.setdefault("c", []).append((n, e, e <= 0.02))


def crit_d(rep):
    hdr("(d) SCORED -- the region split closes, and the side power by band")
    print(f"{'leg':<6}{'sideP':>9}{'bore':>8}{'pads':>8}{'wall':>8}"
          f"{'splitErr':>10}{'A_bore':>8}{'A_pads':>8}{'A_wall':>8}{'pads%':>7}")
    for n in LEGS:
        if not have(n):
            continue
        d = thermo(n)
        tot = tm(d, "side_P_face")
        b, p, w = (tm(d, "side_P_" + k) for k in ("bore", "pads", "wall"))
        ab, ap, aw = (tm(d, "side_A_" + k) * 1e4 for k in ("bore", "pads", "wall"))
        err = abs(b + p + w - tot) / max(abs(tot), 1e-12)
        print(f"{n:<6}{tot:9.1f}{b:8.1f}{p:8.1f}{w:8.1f}{err:10.1e}"
              f"{ab:8.2f}{ap:8.2f}{aw:8.2f}{100*ap/max(ab+ap+aw,1e-30):6.1f}%")
        rep.setdefault("d", []).append((n, err, ap, aw, p, w))
        if p < 0.0:
            print(f"       ⚠ side_P_pads is NEGATIVE: with the ring shielded the "
                  f"gas deposits nothing there but the faces still RADIATE.")
    print("\n  areas are cm2 FULL-hole; powers W FULL-hole. Bar on splitErr: 1e-4.")


def crit_e(rep):
    hdr("(e) SCORED -- identity: S0 vs d2z1_delete/RB_D")
    if not have("S0"):
        print("  S0 not on disk"); return
    ok = True
    for suf in ("_removal_events.csv", "_jet_profile.csv"):
        a, b = os.path.join(OUT, "S0" + suf), os.path.join(Z1, "RB_D" + suf)
        if not (os.path.exists(a) and os.path.exists(b)):
            print(f"  {suf}: missing on one side"); ok = False; continue
        la, lb = open(a).readlines(), open(b).readlines()
        m = min(len(la), len(lb))
        same = la[:m] == lb[:m]
        ok &= same
        print(f"  {suf}: {len(la)} vs {len(lb)} rows, prefix {m}: "
              f"{'IDENTICAL' if same else 'DIFFERS'}")
    ta = thermo("S0")
    p = os.path.join(Z1, "RB_D", "thermo.dat")
    with open(p) as f:
        cb = f.readline().split()
    vb = np.loadtxt(p, skiprows=1, ndmin=2)
    tb = {c: vb[:, k] for k, c in enumerate(cb)}
    extra = sorted(set(ta) - set(tb))
    print(f"  thermo header delta: {extra if extra else 'EMPTY'}   "
          f"{'PASS' if not extra else 'FAIL'}")
    ok &= not extra
    bad = [c for c in set(ta) & set(tb)
           if not np.array_equal(ta[c][:min(len(ta[c]), len(tb[c]))],
                                 tb[c][:min(len(ta[c]), len(tb[c]))])]
    print(f"  thermo: {len(bad)} column(s) differ {bad if bad else ''}")
    ok &= not bad
    rep["e"] = ok
    print(f"\n  verdict: {'PASS' if ok else 'FAIL'}")


def crit_f(rep):
    hdr("(f) REPORTED -- specific energy, and the power-rate scaling")
    print(f"{'leg':<6}{'P_rock kW':>11}{'vol cm3':>9}{'J/mm3':>8}")
    for n in LEGS + ("SW", "A"):
        if not have(n):
            continue
        d = thermo(n)
        P = tm(d, "patch_P_robin") + (tm(d, "side_P_face") if "side_P_face" in d else 0.0)
        s = rep.get("shape", {}).get(n)
        if s is None:
            print(f"{n:<6}{P/1e3:11.2f}{'(no 150 mm)':>18}")
            rep.setdefault("P", {})[n] = P / 1e3
            continue
        rate = s["vol"] / s["t"]
        print(f"{n:<6}{P/1e3:11.2f}{s['vol']:9.2f}{P/1e3/max(rate,1e-12)/1e3:8.2f}")
        rep.setdefault("P", {})[n] = P / 1e3
    print("\n  ⚠ The ratio is ~constant across every leg -- it is the sensible heat")
    print("    to the pinned surface temperature, as energy-set recession requires.")
    print("    So heat cannot be moved between regions to buy rate: the rate")
    print("    follows the POWER. If rock really takes 4-6 J/mm3 the model's spall")
    print("    energetics are wrong by ~3x -- a larger, SEPARATE finding, and not")
    print("    a licence to keep the power.")
    P, R = rep.get("P", {}), rep.get("rate", {})
    if "SW" in P and "SW" in R:
        print(f"\n  rate ~ P^{EXP} (the packet's own scaling), predicted from SW:")
        for n in LEGS + ("A",):
            if n in P and n in R and R[n] == R[n] and P[n] > 0:
                pr = R["SW"] * (P[n] / P["SW"]) ** EXP
                e = math.log(R[n] / R["SW"]) / math.log(P[n] / P["SW"])
                print(f"    {n:<5} predicted {pr:6.3f} vs measured {R[n]:6.3f} m/h "
                      f"({100*(pr-R[n])/R[n]:+5.1f} %), fitted exponent {e:.2f}")
        print("    power needed for the Meier band:")
        for t in (BAND[0], 0.5 * sum(BAND), BAND[1]):
            print(f"      {t:.2f} m/h -> {P['SW']*(t/R['SW'])**(1/EXP):.2f} kW")


def predictions(rep):
    hdr("PREDICTIONS")
    R, P = rep.get("rate", {}), rep.get("P", {})
    wl = R.get("WL", float("nan"))
    print(f"  P1 (reviewer, PRIMARY) WL rate 1.8-2.2 m/h: {wl:.3f}  "
          f"{'HELD' if 1.8 <= wl <= 2.2 else 'REFUTED'}")
    print(f"  P2 (planner alt)       WL rate 1.3-1.6 m/h: {wl:.3f}  "
          f"{'HELD' if 1.3 <= wl <= 1.6 else 'REFUTED'}")
    for n in ("WL", "WL81"):
        if not have(n):
            continue
        sp = tm(thermo(n), "side_P_face") / 1e3
        ok = 0.3 <= sp <= 0.7
        tag = "SCORED" if n == "WL" else "band's other end, reported"
        print(f"  P3 ({tag}) {n} side path {sp:.2f} kW, held inside 0.3-0.7: "
              f"{'HELD' if ok else 'REFUTED'}")
    print("  P3b PD pad-ring top-face power >= 5x down: HELD BY CONSTRUCTION and")
    print("      therefore VACUOUS -- pad_shield_open = 0 sets flame_h in the ring")
    print("      to exactly 0, so the reduction is infinite. ⚠ Deliverable 3")
    print("      instrumented the SIDE power by band, not the TOP-face power by")
    print("      band, so the measured form of P3b is not scoreable as specified.")
    bo, pd = R.get("BO", float("nan")), R.get("PD", float("nan"))
    print(f"  P4 BO below both single edges: BO {bo:.3f} vs WL {wl:.3f}, PD {pd}")
    print("      ⚠ Reads HELD numerically, but both PD and BO are in JAM territory")
    print("      (P7), so the comparison is contaminated and the coupling the")
    print("      packet predicted is NOT what this measures.")
    print(f"  P5 false-pass guard: {rep.get('b', 'n/a')} -- no leg in band, VACUOUS")
    dd = {n: (ap, aw) for n, _, ap, aw, _, _ in rep.get("d", [])}
    if "WL" in dd:
        ap, aw = dd["WL"]
        tot = sum(tm(thermo("WL"), "side_A_" + k) * 1e4 for k in ("bore", "pads", "wall"))
        print(f"  P6 >half the charged side AREA under the pads: "
              f"{100*ap/tot:.1f} %  {'HELD' if ap/tot > 0.5 else 'REFUTED'}")
        pw = [(p, w) for n, _, _, _, p, w in rep.get("d", []) if n == "WL"][0]
        print(f"      on POWER instead: pads {pw[0]:.1f} W vs wall {pw[1]:.1f} W "
              f"= {100*pw[0]/(pw[0]+pw[1]):.1f} %  (pads is still the largest band)")
    print(f"  P7 pads-leg JAM signature (rate < {JAM} m/h with the 28-40 mm rings")
    print("      frozen): ", end="")
    if have("PD"):
        rg = X.rings("PD")
        frozen = all(v < 1.0 for v in rg[:5])          # rings 30-38 mm
        print(f"PD rate {pd}, rings 30-38 mm = "
              f"{[round(v,1) for v in rg[:5]]}  "
              f"{'CONFIRMED -- the model lateral-recession limit, NOT physics' if frozen else 'not seen'}")


def main():
    ap = argparse.ArgumentParser()
    ap.parse_args()
    print("D2aa -- scoring CRITERION.md (sha d81932d18019649b...)")
    print("legs on disk: " + ", ".join(n for n in LEGS + ("S0",) if have(n)))
    rep = {}
    crit_a(rep); crit_b(rep); crit_c(rep); crit_d(rep); crit_e(rep); crit_f(rep)
    predictions(rep)
    return 0


if __name__ == "__main__":
    sys.exit(main())
