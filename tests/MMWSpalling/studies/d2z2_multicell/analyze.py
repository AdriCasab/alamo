#!/usr/bin/env python3
"""D2z-2 analysis -- scores CRITERION.md (a)-(j) and the seven predictions.

  analyze.py                 every finished leg
  analyze.py --quick         skip the offline population passes

Imported UNCHANGED, as CRITERION.md fixes them:
  * `col_rate` and the .done reader from studies/d2j0b_plane/analyze.py
  * `shape`, `rop`, `rings`, `pmean`, the volume bases from
    studies/d2x_areasplit/analyze.py -- which already caps depths 50 mm above
    the burner. Its WHERE lookup is EXTENDED with this study's legs (a table
    entry, not an edit to a function).
Frozen studies are read by path and never re-run.
"""
import argparse
import ast
import filecmp
import importlib.util
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


J = load_mod("d2j0b_an", os.path.join(STUD, "d2j0b_plane", "analyze.py"))
X = load_mod("d2x_an",   os.path.join(STUD, "d2x_areasplit", "analyze.py"))
DO = load_mod("d2z2_off", os.path.join(HERE, "delta_offline.py"))

col_rate = J.col_rate
meta_of = J.R.J.meta_of

OUT = os.path.join(HERE, "output")
SH = os.path.join(STUD, "d2z_stepheight", "output")     # frozen D2z legs
Z1 = os.path.join(STUD, "d2z1_delete", "output")        # frozen whole-cell

WHERE = {"R2": OUT, "R2w": OUT, "R2d": OUT, "R1d": OUT,
         "SW": OUT, "SWs": OUT, "SW0": OUT, "S0": OUT,
         "S2": SH, "S2w": SH, "S1": SH, "S1w": SH, "S2f": SH,
         "W2": Z1, "W1": Z1, "RB_D": Z1}
# Extend D2x's lookup so its shape/rop/rings work on this study's Meier legs.
X.WHERE.update({n: OUT for n in ("SW", "SWs", "SW0", "S0")})

WIN = (50.0, 200.0)
R_GATE, R_WIDE = 0.025, 0.030
FROZEN_GAP_25, FROZEN_GAP_30 = 1.39, -0.71     # D2z step-law, r<25 / r<30 [%]


def have(n):
    return os.path.exists(os.path.join(WHERE[n], n + ".done"))


def meta(n):
    return ast.literal_eval(open(os.path.join(WHERE[n], n + ".done")).read())


def events(n):
    m = meta(n)
    c = (np.arange(m["nxy"]) + 0.5) * m["dz"]
    r = np.hypot(*np.meshgrid(c, c, indexing="ij"))
    a = np.loadtxt(os.path.join(WHERE[n], n + "_removal_events.csv"),
                   delimiter=",", skiprows=1, ndmin=2)
    return m, r, a


def thermo(n):
    p = os.path.join(WHERE[n], n, "thermo.dat")
    with open(p) as f:
        cols = f.readline().split()
    v = np.loadtxt(p, skiprows=1, ndmin=2)
    return {c: v[:, k] for k, c in enumerate(cols)}


def tmean(th, key, win=WIN):
    if key not in th:
        return float("nan")
    s = (th["time"] >= win[0]) & (th["time"] <= win[1])
    return float(np.nanmean(th[key][s])) if s.any() else float("nan")


def disk_rate(n, rmax=R_GATE, win=WIN):
    m, r, a = events(n)
    return float(col_rate(r, a, win)[r < rmax].mean())


def pct(x):
    return f"{100*x:+.2f} %"


def hdr(t):
    print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78)


# ---------------------------------------------------------------- (a) repairs
def crit_a(rep):
    hdr("(a) SCORED -- the diagnostic repairs are diagnostics-only")
    if not (have("R2") and os.path.exists(os.path.join(SH, "S2_removal_events.csv"))):
        print("  R2 not on disk"); return
    same = filecmp.cmp(os.path.join(OUT, "R2_removal_events.csv"),
                       os.path.join(SH, "S2_removal_events.csv"), shallow=False)
    rep["a"] = same
    print(f"  R2_removal_events.csv vs frozen S2: "
          f"{'BYTE-IDENTICAL' if same else 'DIFFERS'}   "
          f"{'PASS' if same else 'FAIL'}")
    print("  (this is what licenses scoring R2d/R1d against the frozen S2/S1)")
    if not same:
        a = open(os.path.join(OUT, "R2_removal_events.csv")).readlines()
        b = open(os.path.join(SH, "S2_removal_events.csv")).readlines()
        d = [i for i, (x, y) in enumerate(zip(a, b)) if x != y]
        print(f"  ⚠ {len(a)} vs {len(b)} rows; first differing row "
              f"{d[0] if d else 'n/a'}")


# ------------------------------------------------------------- (b) estimator
def crit_b(rep, quick):
    hdr("(b) SCORED -- the estimator, both populations specified identically")
    if quick:
        print("  skipped (--quick)"); return
    jobs = [("R2w", R_GATE, 0.0, 0.0), ("R2w", R_WIDE, 0.0, 0.0)]
    if have("SW"):
        jobs.append(("SW", 0.050, None, None))
    for n, rmax, x0, y0 in jobs:
        if not have(n):
            continue
        m = meta(n)
        if x0 is None:                     # Meier: patch centre from the argv
            x0 = y0 = m["lxy"] / 2.0
        th = thermo(n)
        live = tmean(th, "delta_mean")
        both = {}
        for req, tag in ((False, "PRE-REGISTERED (all columns, phi=1)"),
                         (True,  "MEASURED live population (taller has fired)")):
            o = DO.run(os.path.join(OUT, n + "_removal_events.csv"), m["nxy"],
                       m["dz"], rmax, [WIN],
                       f"{n} offline r<{rmax*1e3:.0f} mm -- {tag}",
                       x0=x0, y0=y0, require_tall_fired=req)
            both[req] = o[f"{int(WIN[0])}-{int(WIN[1])}"]["tw_mean"]
        rel = (live - both[False]) / both[False]
        ok = abs(rel) <= 0.10
        rel2 = (live - both[True]) / both[True]
        print(f"  {n} r<{rmax*1e3:.0f}: live {live:.4f}")
        print(f"      vs PRE-REGISTERED offline {both[False]:.4f}  rel {pct(rel)}"
              f"   {'PASS' if ok else 'FAIL'}   <- SCORED")
        print(f"      vs live-population offline {both[True]:.4f}  rel {pct(rel2)}"
              f"   {'(within 10 %)' if abs(rel2) <= 0.10 else '(outside 10 %)'}"
              f"   <- reported")
        rep.setdefault("b", []).append((n, rmax, live, both[False], rel, ok, rel2))


# ------------------------------------------------------------- (c) the rim test
def crit_c(rep):
    hdr("(c) SCORED -- THE RIM TEST: does the derived form reopen D2z's cure?")
    if not (have("R1d") and have("R2d")):
        print("  needs R1d and R2d"); return
    for rmax, tag, frozen in ((R_GATE, "r<25", FROZEN_GAP_25),
                              (R_WIDE, "r<30", FROZEN_GAP_30)):
        vf, vc = disk_rate("R1d", rmax), disk_rate("R2d", rmax)
        g = 100.0 * (vf - vc) / vc
        # the same pair on the frozen D2z legs, recomputed not quoted
        fz = 100.0 * (disk_rate("S1", rmax) - disk_rate("S2", rmax)) / disk_rate("S2", rmax)
        print(f"  {tag}: R1d {vf:.4f} vs R2d {vc:.4f}  gap {g:+.2f} %   "
              f"(frozen S1/S2 recomputed {fz:+.2f} %, packet {frozen:+.2f} %)")
        if tag == "r<25":
            d = abs(g - fz)
            print(f"        interior is 0.0 % multi-cell, so n = 1 and the derived")
            print(f"        form is EXACT there: |shift| = {d:.2f} pts "
                  f"(bar 0.5)   {'PASS' if d <= 0.5 else 'FAIL'}")
            print(f"        |gap| <= 5 %: {'PASS' if abs(g) <= 5 else 'FAIL'}")
            rep["c25"] = (g, fz, d)
        else:
            print(f"        SCORED vs the frozen {fz:+.2f} %: "
                  f"|gap| = {abs(g):.2f} %  "
                  f"{'PASS -- the rim cure was population + drain, and is real'
                     if abs(g) <= 5 else
                     'FAIL -- part of the cure was the double count'}")
            rep["c30"] = (g, fz)


# ------------------------------------------------------------- (d) drain meter
def crit_d(rep):
    hdr("(d) SCORED -- the drain meter agrees with the model it reports")
    print(f"{'leg':<5}{'q_in':>11}{'qeff_in':>11}{'flux g':>9}{'step_g_mean':>13}"
          f"{'|diff|':>9}{'n_in':>8}{'n':>8}{'g_max':>8}{'f_multi':>9}")
    for n in ("R2", "R2w", "R2d", "R1d", "SW", "SWs", "SW0"):
        if not have(n):
            continue
        th = thermo(n)
        qi, qe = tmean(th, "wall_step_q_in"), tmean(th, "wall_step_q_eff_in")
        gm = tmean(th, "step_g_mean")
        fl = qe / qi if qi else float("nan")
        print(f"{n:<5}{qi:11.1f}{qe:11.1f}{fl:9.4f}{gm:13.4f}{abs(fl-gm):9.4f}"
              f"{tmean(th,'wall_step_n_in'):8.1f}{tmean(th,'wall_step_n'):8.1f}"
              f"{tmean(th,'step_g_max'):8.4f}{tmean(th,'wall_step_f_multi'):9.4f}")
        rep.setdefault("d", []).append((n, fl, gm, abs(fl - gm)))
    print("\n  bar: |flux-weighted - face-weighted| within 5 % -- and see")
    print("  CRITERION.md (d): these are two DIFFERENT means of the same g,")
    print("  and q_raw changes sign, so the bar measures their CORRELATION.")
    for n, fl, gm, d in rep.get("d", []):
        rel = d / gm if gm else float("nan")
        print(f"    {n:<5} rel {100*rel:6.2f} %   {'PASS' if rel <= 0.05 else 'FAIL'}")
    print("\n  the check with teeth -- ungated vs in-gate, and the contamination:")
    for n in ("R2", "R2d", "SW"):
        if not have(n):
            continue
        th = thermo(n)
        u = tmean(th, "wall_step_q_eff") / tmean(th, "wall_step_q")
        g = tmean(th, "wall_step_q_eff_in") / tmean(th, "wall_step_q_in")
        print(f"    {n:<5} ungated ratio {u:.4f}  vs in-gate {g:.4f}   "
              f"(D2z read 0.843 ungated where in-patch was 0.335)")


# ----------------------------------------------------------- (e) power identity
def crit_e(rep):
    hdr("(e) SCORED -- the power identity, +/- 0.02 W, on the Meier legs")
    for n in ("SW", "SWs", "SW0", "S0"):
        if not have(n):
            continue
        th = thermo(n)
        if "jet_P_face" not in th:
            print(f"  {n:<5} NOT SCOREABLE: no jet_P_face (no march)"); continue
        s = (th["time"] >= WIN[0]) & (th["time"] <= WIN[1])
        e = float(np.abs((th["jet_P_face"] - (th["patch_P_robin"] + th["side_P_face"]
                                              + th["side_P_loss"]))[s]).max())
        print(f"  {n:<5} max residual {e:.4e} W   {'PASS' if e <= 0.02 else 'FAIL'}"
              f"   (6-sig-fig floor 0.015 W)")
        rep.setdefault("e", []).append((n, e, e <= 0.02))


# ----------------------------------------------------------------- (f) identity
def crit_f(rep):
    hdr("(f) SCORED -- identity: S0 vs d2z1_delete/RB_D")
    if not have("S0"):
        print("  S0 not on disk"); return
    ok = True
    for suf in ("_removal_events.csv", "_jet_profile.csv"):
        a, b = os.path.join(OUT, "S0" + suf), os.path.join(Z1, "RB_D" + suf)
        if not (os.path.exists(a) and os.path.exists(b)):
            print(f"  {suf}: missing on one side"); continue
        la, lb = open(a).readlines(), open(b).readlines()
        n = min(len(la), len(lb))
        s = la[:n] == lb[:n]
        ok &= s
        print(f"  {suf}: {len(la)} vs {len(lb)} rows, prefix {n}: "
              f"{'IDENTICAL' if s else 'DIFFERS'}")
    ta, tb = thermo("S0"), thermo("RB_D")
    extra = sorted(set(ta) - set(tb))
    print(f"  thermo header delta: {extra if extra else 'EMPTY'}   "
          f"{'PASS' if not extra else 'FAIL -- a new column appeared with the keys off'}")
    ok &= not extra
    bad = [c for c in set(ta) & set(tb)
           if not np.array_equal(ta[c][:min(len(ta[c]), len(tb[c]))],
                                 tb[c][:min(len(ta[c]), len(tb[c]))])]
    print(f"  thermo: {len(bad)} column(s) differ {bad if bad else ''}")
    ok &= not bad
    rep["f"] = ok
    print(f"\n  verdict: {'PASS' if ok else 'FAIL'}")
    print("  ⚠ `metadata` and `diff.patch` are run provenance (compile date,")
    print("    plot_file path, the ParmParse table -- query_default ADDS a key")
    print("    to it -- and a start stamp). The 2639-file tree list does not")
    print("    hash them either; excluding them is the same rule.")


# ------------------------------------------------------- (g)(h)(i) Meier shape
def crit_ghi(rep):
    hdr("(g)(h)(i) Meier at equal burner depth")
    legs = [n for n in ("SW", "SWs", "SW0") if have(n)]
    if not legs:
        print("  no Meier leg on disk yet -- (g), (h), (i), P5-P7 unscored")
        return
    ref = [n for n in ("A", "D") if X.have(n)]
    print(f"  legs {legs}, reference {ref}, primary descent {X.PRIMARY} mm")
    print(f"\n{'leg':<5}{'mean|dD|':>10}{'D@100':>8}{'D@125':>8}{'ROP':>8}"
          f"{'pit lead':>10}{'vol':>9}{'x_pit':>8}{'x_desc':>8}")
    S = {}
    for n in legs + ref:
        d = X.shape(n, X.PRIMARY, (10, 25, 50, 75, 100, 125))
        if d is None:
            print(f"{n:<5}  (no t at {X.PRIMARY} mm descent)"); continue
        S[n] = d
        print(f"{n:<5}{(d['a'] if d['a'] else float('nan')):10.2f}"
              f"{d['D'].get(100, float('nan')):8.1f}{d['D'].get(125, float('nan')):8.1f}"
              f"{X.rop(n):8.3f}{d['lead']:10.1f}{d['vol']:9.1f}"
              f"{d['x_pit']:8.2f}{d['x_desc']:8.2f}")
    print("\n  ⚠ rate definition: burner DESCENT slope over the 100-150 mm")
    print("    stretch (A reads 1.251 there; 1.296 over the 150-400 s window).")
    print("  ⚠ Meier Fig. 8.8 is a visible width = a LOWER BOUND on diameter;")
    print("    band 85-95 +/- 5 mm. No agreement claim tighter than +/- 5 mm.")
    print(f"\n  ring table, 2 mm rings {X.RINGS[0]}-{X.RINGS[-1]+2} mm "
          f"[mm recession], per ring, never averaged:")
    print("       " + "".join(f"{r:7d}" for r in X.RINGS))
    for n in legs + ref:
        if X.have(n):
            print(f"  {n:<5}" + "".join(f"{v:7.2f}" for v in X.rings(n)))
    if "SW" in S and "SWs" in S:
        d = abs(S["SW"]["a"] - S["SWs"]["a"])
        print(f"\n  (h) law A/B: wetted {S['SW']['a']:.2f} vs step "
              f"{S['SWs']['a']:.2f} mm, |diff| {d:.2f} mm")
        print(f"      {'DISTINGUISHED: ' + ('wetted' if S['SW']['a'] < S['SWs']['a'] else 'step') + ' wins'
                 if d > 0.25 else
                 'UNDECIDED (<= 0.25 mm, a scoring-depth systematic); tie-break: wetted'}")
        rep["h"] = d
    if "SW" in S and "SW0" in S:
        print(f"\n  (i) REPORTED, derived vs D2z g: mean|dD| "
              f"{S['SW']['a']:.2f} vs {S['SW0']['a']:.2f} mm, "
              f"ROP {X.rop('SW'):.3f} vs {X.rop('SW0'):.3f} m/h, "
              f"vol {S['SW']['vol']:.1f} vs {S['SW0']['vol']:.1f} cm3")
    rep["S"] = S


# ------------------------------------------------------- (j) the multi-cell witness
def crit_j(rep, quick):
    hdr("(j) REQUIRED WITNESS -- the multi-cell population, per leg")
    print(f"{'leg':<5}{'f_multi':>10}{'side_P_face':>13}{'P_multi':>10}"
          f"{'multi share':>13}{'g_max':>8}{'delta_mean':>11}{'delta_faces':>12}")
    for n in ("R2", "R2w", "R2d", "R1d", "SW", "SWs", "SW0"):
        if not have(n):
            continue
        th = thermo(n)
        sp, sm = tmean(th, "side_P_face"), tmean(th, "side_P_face_multi")
        print(f"{n:<5}{tmean(th,'wall_step_f_multi'):10.4f}{sp:13.3f}{sm:10.3f}"
              f"{(sm/sp if sp else float('nan')):13.4f}"
              f"{tmean(th,'step_g_max'):8.4f}{tmean(th,'delta_mean'):11.4f}"
              f"{tmean(th,'delta_faces'):12.1f}")
    print("\n  ⚠ step_g_max is bounded by (n+1)/n, NOT by 1: delta spans")
    print("    [n-1, n+1] through the sub-cell phase. The bound FALLS with n")
    print("    (2.00 at n=1 -> 1.14 at n=7); D2z's ROSE (1.99 -> 7.73).")
    if have("SW"):
        fm = tmean(thermo("SW"), "wall_step_f_multi")
        print(f"\n  P5: SW wall_step_f_multi = {fm:.4f}  "
              f"{'HELD (>= 0.20)' if fm >= 0.20 else 'REFUTED (< 0.20) -- (i) is UNINFORMATIVE'}")
        rep["p5"] = fm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    print("D2z-2 -- scoring CRITERION.md (sha 46571c3c884b976a...)")
    print("legs on disk: " + ", ".join(n for n in WHERE if have(n)))
    print("conductance per leg: " + ", ".join(
        f"{n}={meta(n).get('conductance', '?')}" for n in
        ("R2", "R2w", "R2d", "R1d", "SW", "SWs", "SW0", "S0") if have(n)))
    rep = {}
    crit_a(rep); crit_b(rep, a.quick); crit_c(rep); crit_d(rep)
    crit_e(rep); crit_f(rep); crit_ghi(rep); crit_j(rep, a.quick)

    hdr("PREDICTIONS")
    if "a" in rep:
        print(f"  P1 R2 byte-identical to S2: {'HELD' if rep['a'] else 'REFUTED'}")
    for n, rmax, live, off, rel, ok, rel2 in rep.get("b", []):
        print(f"  P2 {n} r<{rmax*1e3:.0f} live-vs-offline {pct(rel)}: "
              f"{'HELD' if ok else 'REFUTED'}   "
              f"(live-population reading {pct(rel2)})")
    if have("R2d"):
        s = filecmp.cmp(os.path.join(OUT, "R2d_removal_events.csv"),
                        os.path.join(SH, "S2_removal_events.csv"), shallow=False)
        print(f"  P3a R2d byte-identical to S2 (derived inert at 2 mm): "
              f"{'HELD' if s else 'REFUTED'}")
        if not s:
            print("      MECHANISM (measured on R2's own log, 50-200 s): 42.9 % of")
            print("      stepped lateral faces at 2 mm are CRATER WALL -- the taller")
            print("      column never fired, so col_dt = 0 and a = 0 -- and the two")
            print("      forms disagree on 71.7 % of those. The packet assumed")
            print("      delta = n there (phi_S = phi_T = 1), where both return 1.")
            print("      But phi_T = 1 only because the TALLER column never fired;")
            print("      the SHORTER column is ACTIVE, so phi_S < 1 and")
            print("      delta = n + phi_S - 1 < n. D2z's a->0 branch returns")
            print("      min(delta,1) = 1; the derived one returns delta/n ~ 0.93.")
            print("      Small per face, but on ~31 % of all stepped faces, so the")
            print("      trajectory diverges (first difference: row 343).")
    if have("R1d"):
        s = filecmp.cmp(os.path.join(OUT, "R1d_removal_events.csv"),
                        os.path.join(SH, "S1_removal_events.csv"), shallow=False)
        print(f"  P3b R1d DIFFERS from S1 (the change reached a multi-cell face): "
              f"{'REFUTED -- identical, so the test is PLUMBING not physics' if s else 'HELD'}")
    if "c30" in rep:
        g = rep["c30"][0]
        print(f"  P4 rim holds at r<30: |gap| {abs(g):.2f} % "
              f"{'HELD' if abs(g) <= 5 else 'REFUTED'}")
    if "p5" in rep:
        print(f"  P5 SW f_multi {rep['p5']:.3f} "
              f"{'HELD' if rep['p5'] >= 0.20 else 'REFUTED'}")
    S = rep.get("S", {})
    if "SW" in S:
        d100 = S["SW"]["D"].get(100, float("nan"))
        print(f"  P6 SW D@100 = {d100:.1f} mm "
              f"{'HELD (< 87)' if d100 < 87 else 'REFUTED (>= 87)'}")
        r = X.rop("SW")
        print(f"  P7 SW descent {r:.3f} m/h "
              f"{'HELD (1.3-1.6)' if 1.3 <= r <= 1.6 else 'REFUTED'}")
    else:
        print("  P5-P7 need the Meier legs")
    return 0


if __name__ == "__main__":
    sys.exit(main())
