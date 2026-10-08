#!/usr/bin/env python3
"""D2z analysis -- scores CRITERION.md (a)-(f2) and reports (g)/(h).

  analyze.py                    every finished leg
  analyze.py --no-populations   skip (f2) (the slow offline pass)

`col_rate` and the metadata reader are IMPORTED UNCHANGED from
`studies/d2j0b_plane/analyze.py`, as CRITERION.md fixes them. Frozen and
re-baselined studies are read by path and never re-run.
"""
import argparse
import importlib.util
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)
ROOT = os.path.abspath(os.path.join(STUD, "..", "..", ".."))

# --- imported unchanged: col_rate + the .done reader -------------------------
_sp = importlib.util.spec_from_file_location(
    "d2j0b_analyze", os.path.join(STUD, "d2j0b_plane", "analyze.py"))
A = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(A)
col_rate = A.col_rate                  # events -> per-column m/h over a window
meta_of = A.R.J.meta_of                # ast.literal_eval of <prefix>.done

_sp2 = importlib.util.spec_from_file_location(
    "delta_offline", os.path.join(HERE, "delta_offline.py"))
DO = importlib.util.module_from_spec(_sp2)
_sp2.loader.exec_module(DO)

OUT = os.path.join(HERE, "output")
REF = os.path.join(STUD, "d2z1_delete", "output")

WIN = (50.0, 200.0)                    # the SINGLE pre-registered window
BLOCKS = [(50.0, 100.0), (100.0, 150.0), (150.0, 200.0)]
R_GATE, R_WIDE = 0.025, 0.030
KAPPA, RHO, CP = 1.5, 2750.0, 790.0
ALPHA = KAPPA / (RHO * CP)

# (a)'s reference: the whole-cell pair, recomputed here rather than quoted.
LEGS = {"S2": OUT, "S1": OUT, "S2w": OUT, "S1w": OUT, "S2f": OUT, "S0": OUT,
        "SW": OUT, "W2": REF, "W1": REF, "RB_D": REF}


def have(n):
    return os.path.exists(os.path.join(LEGS[n], n + ".done"))


def events(n):
    d = LEGS[n]
    m = meta_of(os.path.join(d, n))
    dz = m["dz"]
    c = (np.arange(m["nxy"]) + 0.5) * dz
    X, Y = np.meshgrid(c, c, indexing="ij")
    r = np.hypot(X, Y)                 # disk at the ORIGIN corner
    a = np.loadtxt(os.path.join(d, n + "_removal_events.csv"), delimiter=",",
                   skiprows=1, ndmin=2)
    return m, r, a


def thermo(n):
    d = LEGS[n]
    p = os.path.join(d, n, "thermo.dat")
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


# ---------------------------------------------------------------- (a) the gate
def crit_a(rep):
    print("\n" + "=" * 78)
    print("(a) SCORED -- THE GATE: interior mesh gap, r < 25 mm, SINGLE 50-200 s")
    print("=" * 78)
    rows = []
    for fine, coarse, lab in (("S1", "S2", "step"), ("S1w", "S2w", "wetted"),
                              ("W1", "W2", "whole-cell ref (f = 0.2)")):
        if not (have(fine) and have(coarse)):
            continue
        for rmax, tag in ((R_GATE, "r<25"), (R_WIDE, "r<30")):
            vf, vc = disk_rate(fine, rmax), disk_rate(coarse, rmax)
            g = (vf - vc) / vc
            rows.append((lab, tag, vf, vc, g))
    print(f"{'law':<26}{'disc':<7}{'1 mm':>9}{'2 mm':>9}{'gap':>11}")
    for lab, tag, vf, vc, g in rows:
        print(f"{lab:<26}{tag:<7}{vf:9.4f}{vc:9.4f}{pct(g):>11}")
    gate = {lab: g for lab, tag, _, _, g in rows if tag == "r<25"}
    rep["gate"] = gate
    print("\nper-block (REPORTED -- inherits block noise from BOTH meshes):")
    for fine, coarse, lab in (("S1", "S2", "step"), ("S1w", "S2w", "wetted"),
                              ("W1", "W2", "whole-cell ref")):
        if not (have(fine) and have(coarse)):
            continue
        b = []
        for w in BLOCKS:
            vf = disk_rate(fine, R_GATE, w)
            vc = disk_rate(coarse, R_GATE, w)
            b.append(100.0 * (vf - vc) / vc)
        print(f"  {lab:<22}" + "".join(f"{x:+8.2f}" for x in b) + "   [%]")
    print("\nverdict:")
    for lab, g in gate.items():
        if lab.startswith("whole-cell"):
            print(f"  {lab:<26} {pct(g)}   (the bar to beat)")
            continue
        v = ("PASS <=2 % (as good as the fudge)" if abs(g) <= 0.02 else
             "PASS <=5 % (cured, but worse than f = 0.2)" if abs(g) <= 0.05 else
             "FAIL > 5 %")
        print(f"  {lab:<26} {pct(g)}   {v}")


# --------------------------------------------------------------- (b) the drain
def crit_b(rep):
    print("\n" + "=" * 78)
    print("(b) SCORED -- the flat-floor lateral drain falls (EFFECTIVE pair)")
    print("=" * 78)
    if not (have("S2") and have("S2f")):
        return
    ms = meta_of(os.path.join(OUT, "S2"))
    dx = dz = ms["dz"]
    out = {}
    for n in ("S2", "S2w", "S2f"):
        if not have(n):
            continue
        th = thermo(n)
        q_raw = tmean(th, "wall_step_q")
        q_eff = tmean(th, "wall_step_q_eff")
        nf = tmean(th, "wall_step_n")
        out[n] = dict(q_raw=q_raw, q_eff=q_eff, n=nf,
                      W_raw=q_raw * nf * dx * dz,
                      W_eff=(q_eff if q_eff == q_eff else q_raw) * nf * dx * dz)
    print(f"{'leg':<6}{'q_raw[W/m2]':>14}{'q_eff[W/m2]':>14}{'n_faces':>10}"
          f"{'W_raw':>10}{'W_eff':>10}")
    for n, d in out.items():
        qe = f"{d['q_eff']:14.1f}" if d["q_eff"] == d["q_eff"] else f"{'(off)':>14}"
        print(f"{n:<6}{d['q_raw']:14.1f}{qe}{d['n']:10.2f}"
              f"{d['W_raw']:10.4f}{d['W_eff']:10.4f}")
    if "S2" in out and "S2f" in out:
        # S2f has the key off, so its q_eff IS its q_raw (g = 1 everywhere).
        red = out["S2f"]["W_raw"] / out["S2"]["W_eff"]
        rep["drain_x"] = red
        print(f"\n  S2f / S2(eff) = {red:.2f}x   bar >= 8x   "
              f"{'PASS' if red >= 8.0 else 'FAIL'}")
        # decomposition: per-face g against population
        gf = out["S2"]["q_eff"] / out["S2"]["q_raw"] if out["S2"]["q_raw"] else float("nan")
        popn = out["S2f"]["n"] / out["S2"]["n"] if out["S2"]["n"] else float("nan")
        print(f"  decomposition: per-face g = {gf:.3f} ({1/gf:.2f}x)  x  "
              f"population n(S2f)/n(S2) = {popn:.2f}x  = {popn/gf:.2f}x")
        print(f"  (expected ~4.0x per-face at delta = 0.157, ~4x population)")


# ------------------------------------------------------ (c) the power identity
def crit_c(rep):
    print("\n" + "=" * 78)
    print("(c) SCORED -- the power identity, +/- 0.01 W")
    print("=" * 78)
    miss = []
    for n in ("S2", "S1", "S2w", "S1w", "S2f", "S0", "SW"):
        if not have(n):
            continue
        th = thermo(n)
        if "jet_P_face" not in th:
            miss.append(n)
            continue
        s = (th["time"] >= WIN[0]) & (th["time"] <= WIN[1])
        err = (th["jet_P_face"] - (th["patch_P_robin"] + th["side_P_face"]
                                   + th["side_P_loss"]))[s]
        e = float(np.abs(err).max())
        print(f"  {n:<5} max |jet_P_face - (robin + side + loss)| = {e:.3e} W   "
              f"{'PASS' if e <= 0.01 else 'FAIL'}")
    if miss:
        print(f"\n  ⚠ NOT SCOREABLE on {', '.join(miss)}: no `jet_P_face` column.")
        print("    The arbiter never runs the march (no jet_closure), so the two")
        print("    MARCH call sites are exercised only by the Meier legs (S0/SW),")
        print("    and on S0 the step key is OFF. CRITERION.md calls (c) 'the only")
        print("    check that catches a kernel<->march split' and asks for it on")
        print("    'every leg'; on this geometry that is not possible.")
    rep["c_missing"] = miss
    print("\n  ledger closure (all legs, independent of the march):")
    for n in ("S2", "S1", "S2w", "S1w", "S2f", "S0"):
        if not have(n):
            continue
        th = thermo(n)
        s = (th["time"] >= WIN[0]) & (th["time"] <= WIN[1])
        print(f"    {n:<5} max |ledger_err| = {float(np.abs(th['ledger_err'][s]).max()):.3e}")


# -------------------------------------------------------------- (d) the estimator
def crit_d(rep):
    print("\n" + "=" * 78)
    print("(d) SCORED -- THE ESTIMATOR: live delta_mean vs the offline both-column")
    print("=" * 78)
    band = {2e-3: (0.137, 0.171, 0.157), 1e-3: (0.197, 0.319, 0.270)}
    for n in ("S2", "S2w", "S1", "S1w"):
        if not have(n):
            continue
        m = meta_of(os.path.join(OUT, n))
        th = thermo(n)
        live = tmean(th, "delta_mean")
        fb = tmean(th, "delta_fallback_cols")
        nsf = tmean(th, "side_faces")
        off = DO.run(os.path.join(OUT, n + "_removal_events.csv"), m["dz"],
                     R_GATE, [WIN], f"{n}  offline both-column estimator")
        o = off[f"{int(WIN[0])}-{int(WIN[1])}"]["tw_mean"]
        rel = (live - o) / o
        lo, hi, agg = band[m["dz"]]
        print(f"\n  {n}: live delta_mean = {live:.4f}   offline t-wtd = {o:.4f}"
              f"   rel = {pct(rel)}   {'PASS' if abs(rel) <= 0.10 else 'FAIL'} (bar ~10 %)")
        print(f"     REPORTED vs item-0a whole-cell band {lo:.3f}-{hi:.3f} "
              f"(agg {agg:.3f}): {'inside' if lo <= live <= hi else 'OUTSIDE'}"
              f"  -- below the band is a RESULT, not a failure")
        frac = fb / nsf if nsf else float("nan")
        print(f"     delta_fallback_cols = {fb:.2f} of {nsf:.1f} side faces "
              f"= {100*frac:.2f} %   bar < 5 %   {'PASS' if frac < 0.05 else 'FAIL'}")
        rep.setdefault("delta", {})[n] = (live, o, rel, frac)


# ------------------------------------------------------------ (e) the area law
def crit_e(rep):
    print("\n" + "=" * 78)
    print("(e) SCORED -- WHICH AREA LAW, with the 1 pt tie rule")
    print("=" * 78)
    g = rep.get("gate", {})
    if "step" not in g or "wetted" not in g:
        return
    gs, gw = 100 * g["step"], 100 * g["wetted"]
    print(f"  step   gap {gs:+.2f} %   {'inside' if abs(gs) <= 5 else 'OUTSIDE'} +/-5 %")
    print(f"  wetted gap {gw:+.2f} %   {'inside' if abs(gw) <= 5 else 'OUTSIDE'} +/-5 %")
    d = abs(gs - gw)
    print(f"  |difference| = {d:.2f} pt   (aggregate floor ~0.5-1 pt)")
    if abs(gs) > 5 and abs(gw) > 5:
        v = "BOTH REFUTED -> the answer is D2j's continuous front"
    elif d <= 1.0:
        v = "UNDECIDED (<= 1 pt); pre-registered preference: wetted"
    else:
        v = f"{'step' if abs(gs) < abs(gw) else 'wetted'} wins"
    print(f"  verdict: {v}")
    rep["law"] = v
    print("\n  riser_area_ratio (time-mean, 50-200 s):")
    for n in ("S2", "S2w", "S1", "S1w"):
        if have(n):
            th = thermo(n)
            print(f"    {n:<5} {tmean(th, 'riser_area_ratio'):10.4f}"
                  f"   side_P_face = {tmean(th, 'side_P_face'):8.3f} W"
                  f"   side_faces = {tmean(th, 'side_faces'):7.1f}")


# ----------------------------------------------------------------- (f) identity
def crit_f(rep):
    print("\n" + "=" * 78)
    print("(f) SCORED -- identity: S0 vs d2z1_delete/RB_D, over the common prefix")
    print("=" * 78)
    if not have("S0"):
        return
    import hashlib
    a = os.path.join(OUT, "S0_removal_events.csv")
    b = os.path.join(REF, "RB_D_removal_events.csv")
    la, lb = open(a).readlines(), open(b).readlines()
    n = min(len(la), len(lb))
    same = la[:n] == lb[:n]
    print(f"  removal log: {len(la)} vs {len(lb)} rows, common prefix {n}: "
          f"{'IDENTICAL' if same else 'DIFFER'}")
    ta, tb = thermo("S0"), thermo("RB_D")
    print(f"  thermo columns: {len(ta)} vs {len(tb)}, "
          f"{'same set' if set(ta) == set(tb) else 'DIFFERENT SET'}")
    bad = 0
    for k in set(ta) & set(tb):
        m = min(len(ta[k]), len(tb[k]))
        if not np.array_equal(ta[k][:m], tb[k][:m]):
            bad += 1
            print(f"    differs: {k}")
    print(f"  thermo: {bad} column(s) differ over the common prefix")
    # plotfiles
    nd = 0
    for root, _, fs in os.walk(os.path.join(OUT, "S0")):
        for f in fs:
            p = os.path.join(root, f)
            q = p.replace(os.path.join(OUT, "S0"), os.path.join(REF, "RB_D"))
            if not os.path.exists(q):
                continue
            # `metadata` and `diff.patch` are RUN PROVENANCE, not simulation
            # state: metadata records the compile date, the plot_file path, the
            # ParmParse table (⚠ `query_default` ADDS a key to that table, so
            # the two new keys appear there at their defaults even when unset)
            # and a start-time stamp; diff.patch is the git working-tree diff.
            # None of them is read back by the solver. The 2639-file tree gate
            # does not hash them either (0 `metadata` lines in d2z1.md5), so
            # excluding them here is the same rule, not a relaxation.
            if f in ("thermo.dat", "metadata", "diff.patch"):
                continue
            ha = hashlib.md5(open(p, "rb").read()).hexdigest()
            hb = hashlib.md5(open(q, "rb").read()).hexdigest()
            if ha != hb:
                nd += 1
                print(f"    differs: {os.path.relpath(p, OUT)}")
    print(f"  plotfile/Cell_D: {nd} file(s) differ")
    ok = same and bad == 0 and nd == 0 and set(ta) == set(tb)
    rep["identity"] = ok
    print(f"\n  verdict: {'PASS' if ok else 'FAIL'}")


# ------------------------------------------------------- (f2) step populations
def crit_f2(rep):
    print("\n" + "=" * 78)
    print("(f2) REQUIRED WITNESS -- step populations, every leg (read-only)")
    print("=" * 78)
    for n in ("S2", "S2w", "S2f", "W2", "S1", "S1w", "W1"):
        if not have(n):
            continue
        m = meta_of(os.path.join(LEGS[n], n))
        DO.run(os.path.join(LEGS[n], n + "_removal_events.csv"), m["dz"],
               R_GATE, [(50, 100), (100, 150), (150, 200), WIN],
               f"{n}  dz = {m['dz']*1e3:.0f} mm")


# --------------------------------------------------------- (g) the front proxy
def crit_g(rep):
    print("\n" + "=" * 78)
    print("(g) REPORTED -- the front-proxy error, recomputed from S2/S1")
    print("=" * 78)
    for n in ("S2", "S1"):
        if not have(n):
            continue
        m = meta_of(os.path.join(OUT, n))
        th = thermo(n)
        d = tmean(th, "delta_mean")
        dz = m["dz"]
        # a discrete riser of height delta*dz against a true continuous front:
        # the proxy pays the full vertical face, the front pays the ramp.
        step = d
        ramp = np.sqrt(1.0 + d * d) - 1.0
        print(f"  {n}: dz = {dz*1e3:.0f} mm, delta = {d:.4f}  ->  "
              f"step/ramp = {step/ramp:.2f}x")
    print("  (the packet carried 1.35-1.75 as a claim; at these delta the")
    print("   proxy-to-ramp factor is far larger because delta << 1 makes the")
    print("   ramp second order -- see RESULTS.md)")


def population_split(rep):
    """The evidence for (b)/(d): the live witnesses have NO radius gate, so they
    average the flat floor together with the 25-30 mm rim. Nothing spalls beyond
    30 mm on this geometry, so `r < 30` IS the full live population."""
    import math
    print("\n" + "=" * 78)
    print("POPULATION SPLIT -- why (b) and (d) read what they read")
    print("=" * 78)
    print("The live `delta_mean` / `wall_step_q_eff` are means over EVERY step")
    print("face. The gate disc is r < 25 mm. Offline, on each leg's own log:\n")
    print(f"{'leg':<6}{'interior r<25':>28}{'full r<30 (= live set)':>28}"
          f"{'live':>10}")
    print(f"{'':<6}{'delta   meanG   g>1%':>28}{'delta   meanG   g>1%':>28}")
    for n in ("S2", "S2w", "S1", "S1w"):
        if not have(n):
            continue
        m = meta_of(os.path.join(OUT, n))
        dx = m["dz"]
        th = thermo(n)
        live = tmean(th, "delta_mean")
        cells = []
        for rmax in (R_GATE, R_WIDE):
            cols = DO.voids(os.path.join(OUT, n + "_removal_events.csv"), dx, rmax)
            keys = set(cols)
            pairs = [(c, (c[0] + d[0], c[1] + d[1])) for c in keys
                     for d in ((1, 0), (0, 1)) if (c[0] + d[0], c[1] + d[1]) in keys]
            ds, gs, namp, tot = [], [], 0, 0
            t = WIN[0]
            while t <= WIN[1]:
                for (ca, cb) in pairs:
                    na, pa, _ = DO.state(cols[ca], t)
                    nb, pb, _ = DO.state(cols[cb], t)
                    nd = na - nb
                    if nd == 0:
                        continue
                    d = abs(nd + (pa - pb)) if nd > 0 else abs(-nd + (pb - pa))
                    tall = cols[ca] if nd > 0 else cols[cb]
                    k = sum(1 for x in tall if x <= t)
                    if k < 2:
                        continue
                    dti = tall[k - 1] - tall[k - 2]
                    if dti <= 0:
                        continue
                    a = dx * (dx / dti) / ALPHA
                    # StepConductance, reimplemented here from the closed form
                    g = (min(d, 1.0) if a <= 1e-12 else
                         max(0.0, (1 - math.exp(-d * a)) / (1 - math.exp(-a))))
                    ds.append(d); gs.append(g); tot += 1
                    if g > 1.0:
                        namp += 1
                t += 0.25
            cells.append((sum(ds) / tot, sum(gs) / tot, 100.0 * namp / tot, max(gs))
                         if tot else (float("nan"),) * 4)
        (di, gi, ai, mi), (df, gf, af, mf) = cells
        print(f"{n:<6}{di:10.3f}{gi:8.3f}{ai:9.1f}  "
              f"{df:12.3f}{gf:8.3f}{af:9.1f}  {live:10.4f}")
        rep.setdefault("pop", {})[n] = (di, gi, ai, mi, df, gf, af, mf, live)
    print("\n⚠ `g > 1 %` is the AMPLIFICATION defect, derived here and not taken")
    print("  on trust: StepConductance normalises by (1 - e^-a), so g(1) = 1 and")
    print("  ANY delta > 1 returns g > 1 -- the base face of a multi-cell step is")
    print("  charged MORE than a full face. Max g reached:")
    for n, v in rep.get("pop", {}).items():
        print(f"    {n:<5} interior max g = {v[3]:.3f}   full max g = {v[7]:.3f}")
    print("  On the GATE DISC it is 0.0 % on every leg, so criterion (a) is")
    print("  untouched. It lives entirely in the 25-30 mm rim, which is exactly")
    print("  where the step law's rim bonus was measured -- so that bonus is NOT")
    print("  attributable until the normalisation is settled. See RESULTS.md.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-populations", action="store_true")
    a = ap.parse_args()
    print("D2z step-height closure -- scoring CRITERION.md")
    print(f"window {WIN[0]:.0f}-{WIN[1]:.0f} s, gate disc r < {R_GATE*1e3:.0f} mm")
    print("legs on disk: " + ", ".join(n for n in LEGS if have(n)))
    rep = {}
    crit_a(rep); crit_b(rep); crit_c(rep); crit_d(rep); crit_e(rep); crit_f(rep)
    if not a.no_populations:
        population_split(rep)
        crit_f2(rep)
    crit_g(rep)
    print("\n" + "=" * 78)
    print("PREDICTIONS")
    print("=" * 78)
    g = rep.get("gate", {})
    for lab in ("step", "wetted"):
        if lab in g:
            print(f"  P1 ({lab:<6}) |gap| = {abs(100*g[lab]):.2f} %  "
                  f"{'HELD' if abs(g[lab]) <= 0.05 else 'REFUTED'}")
    if "drain_x" in rep:
        print(f"  P3 drain {rep['drain_x']:.2f}x  "
              f"{'HELD' if rep['drain_x'] >= 8 else 'REFUTED'}")
    if "step" in g and "wetted" in g:
        sp, wp = abs(g["step"]) <= 0.05, abs(g["wetted"]) <= 0.05
        print(f"  P4 step passes and wetted does not: "
              f"{'HELD' if (sp and not wp) else 'REFUTED'}"
              f"   (step {'pass' if sp else 'fail'}, wetted {'pass' if wp else 'fail'})")
    print("  P5 requires SW (Meier) -- not run at the time of this pass"
          if not have("SW") else "  P5 scored below")
    return 0


if __name__ == "__main__":
    sys.exit(main())
