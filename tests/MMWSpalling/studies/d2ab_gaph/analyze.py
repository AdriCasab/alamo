#!/usr/bin/env python3
"""D2ab analysis -- scores CRITERION.md (a)-(f) and P1-P7.

Imports the equal-burner-depth machinery (`shape`, `rop`, `rings`, `descent`,
`t_at`, both volume bases) UNCHANGED from studies/d2x_areasplit/analyze.py, as
CRITERION.md requires; its WHERE lookup is EXTENDED with this study's legs (a
table entry, not an edit to a function). Frozen studies are read by path and
never re-run.

⚠ All Meier thermo powers are QUARTER-hole sums; everything below is printed
full-hole (x4).
"""
import ast
import filecmp
import importlib.util
import os

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
AA = os.path.join(STUD, "d2aa_heatdelivery", "output")
Z2 = os.path.join(STUD, "d2z2_multicell", "output")
Z1 = os.path.join(STUD, "d2z1_delete", "output")

NEW = ("GL", "GH", "WLr", "WLv", "SWc")
X.WHERE.update({n: OUT for n in NEW + ("S0",)})
X.WHERE.update({"WL": AA, "WL81": AA, "SW": Z2, "RB_D": Z1})

WIN = (150.0, 400.0)
BAND = (1.3, 1.6)                  # Meier rate band
EXP = 1.59                         # rate ~ P^EXP

# ---- PRE-REGISTERED, from PREFLIGHT.md section 4. DO NOT EDIT. -------------
PRED = {"GL": dict(h=38.4, rate=2.87, D100=87.1, D125=84.7, lead=63.5, vol=1.2),
        "GH": dict(h=169.3, rate=3.16, D100=87.3, D125=85.1, lead=54.3, vol=1.1)}
BAR_A = 0.15                       # (a) rate within 15 % of prediction
P2_BAND = (2.3, 2.9)               # the PACKET's prediction; item 0 disagreed
WLV_BAND = (1.85, 2.00)            # (b) bar 1
BAR_C = 5.0                        # (c) percent


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


def rel(a, b):
    return 100.0 * (a - b) / b if b else float("nan")


# --------------------------------------------------------------------------
def crit_a(rep):
    hdr("(a) SCORED -- the derived bracket against its OWN pre-registered "
        "prediction")
    print("  ⚠ This scores THE DERIVATION, not Meier. A leg outside Meier's")
    print("    band that matches its own prediction is a PASS.")
    print(f"  h came from PREFLIGHT.md (sha frozen before the first leg): "
          f"GL {PRED['GL']['h']}, GH {PRED['GH']['h']} W/m^2K\n")
    print(f"{'leg':<5}{'h':>7}{'pred':>8}{'meas':>8}{'rel %':>8}  "
          f"{'verdict':<8}{'D@100':>8}{'D@125':>8}{'lead':>8}{'vol':>7}")
    for n in ("GL", "GH"):
        if not have(n):
            print(f"{n:<5}  (not on disk)")
            continue
        r = X.rop(n)
        s = X.shape(n, X.PRIMARY, (10, 25, 50, 75, 100, 125))
        p = PRED[n]
        d = rel(r, p["rate"])
        ok = abs(d) <= BAR_A * 100
        rep.setdefault("rate", {})[n] = r
        rep.setdefault("shape", {})[n] = s
        rep.setdefault("a", []).append((n, ok))
        dd = (f"{s['D'].get(100, float('nan')):8.1f}{s['D'].get(125, float('nan')):8.1f}"
              f"{s['lead']:8.1f}{s['vol']:7.2f}") if s else "   never reached 150 mm"
        print(f"{n:<5}{p['h']:7.1f}{p['rate']:8.2f}{r:8.3f}{d:+8.1f}  "
              f"{'PASS' if ok else 'FAIL':<8}{dd}")
    print(f"\n  predicted shape (REPORTED, not scored): "
          + ", ".join(f"{n} D@100 {PRED[n]['D100']:.1f} D@125 {PRED[n]['D125']:.1f} "
                      f"lead {PRED[n]['lead']:.0f} vol {PRED[n]['vol']:.1f}"
                      for n in ("GL", "GH")))
    print("  ⚠ PREFLIGHT.md recorded a DISAGREEMENT with the packet's P2 band")
    print(f"    ({P2_BAND[0]}-{P2_BAND[1]} m/h) BEFORE launch. Both are scored.")


def crit_b(rep):
    hdr("(b) SCORED -- the attribution pair, ONE edge each")
    print("  WL81 (frozen: const T 1463.5, ALL faces)")
    print("    -> WLv (live T, ALL faces)        the frozen-T fix alone")
    print("    -> WLr (live T, WALL only)        the face assignment alone\n")
    print(f"{'leg':<6}{'T mode':>9}{'faces':>7}{'ROP':>8}{'sideP':>9}"
          f"{'bore':>8}{'pads':>8}{'wall':>8}{'T_exh':>8}{'T_wall':>8}")
    for n, mode, faces in (("WL81", "1463.5", "all"), ("WLv", "live", "all"),
                           ("WLr", "live", "r>40")):
        if not have(n):
            print(f"{n:<6}  (not on disk)")
            continue
        d = thermo(n)
        r = X.rop(n)
        rep.setdefault("rate", {})[n] = r
        rep.setdefault("side", {})[n] = [tm(d, k) for k in
                                         ("side_P_face", "side_P_bore",
                                          "side_P_pads", "side_P_wall")]
        sp = rep["side"][n]
        tw = tm(d, "side_T_wall", 1.0)
        print(f"{n:<6}{mode:>9}{faces:>7}{r:8.3f}{sp[0]:9.1f}{sp[1]:8.1f}"
              f"{sp[2]:8.1f}{sp[3]:8.1f}{tm(d,'jet_T_exhaust',1.0):8.1f}{tw:8.1f}")
    rv, rr = rep.get("rate", {}).get("WLv"), rep.get("rate", {}).get("WLr")
    if rv is not None:
        ok1 = WLV_BAND[0] <= rv <= WLV_BAND[1]
        print(f"\n  bar 1: WLv in {WLV_BAND[0]}-{WLV_BAND[1]} m/h -> {rv:.3f}  "
              f"{'PASS' if ok1 else 'FAIL'}")
        rep.setdefault("b", []).append(("WLv band", ok1))
    if rv is not None and rr is not None:
        ok2 = rr > rv
        print(f"  bar 2: WLr strictly above WLv -> {rr:.3f} vs {rv:.3f}  "
              f"{'PASS' if ok2 else 'FAIL'}")
        rep.setdefault("b", []).append(("WLr > WLv", ok2))


def crit_c(rep):
    hdr("(c) SCORED -- the face assignment is demonstrated, not asserted")
    print("  ⚠ The packet's literal form is UNSCOREABLE: the only all-inherited")
    print("    leg is SW, and SW was frozen BEFORE the region columns existed")
    print("    (0 of 3 in its header). `SWc` was run to supply that baseline.")
    print("    The substitute was fixed in CRITERION.md BEFORE launch.\n")
    print(f"{'leg':<6}{'h':>8}{'r_min':>8}{'bore':>9}{'pads':>9}{'wall':>9}"
          f"{'sideP':>9}")
    got = {}
    for n in ("SWc", "GL", "GH", "WLr", "WLv"):
        if not have(n):
            continue
        d = thermo(n)
        m = ast.literal_eval(open(os.path.join(X.WHERE[n], n + ".done")).read())
        v = [tm(d, k) for k in ("side_P_bore", "side_P_pads", "side_P_wall",
                                "side_P_face")]
        got[n] = v
        print(f"{n:<6}{m['side_face_h']:>8}{m['side_face_h_r_min']:>8}"
              f"{v[0]:9.1f}{v[1]:9.1f}{v[2]:9.1f}{v[3]:9.1f}")

    # (c1) against SWc, the packet's own test.
    print(f"\n  (c1) bore and pads within {BAR_C:.0f} % of SWc (the inherited map):")
    if "SWc" in got:
        for n in ("GL", "GH", "WLr"):
            if n not in got:
                continue
            db, dp = rel(got[n][0], got["SWc"][0]), rel(got[n][1], got["SWc"][1])
            ok = abs(db) <= BAR_C and abs(dp) <= BAR_C
            print(f"       {n:<5} bore {db:+7.2f} %   pads {dp:+7.2f} %   "
                  f"{'PASS' if ok else 'FAIL'}")
            rep.setdefault("c", []).append((f"c1 {n}", ok))
    else:
        print("       SWc not on disk -- (c1) NOT SCOREABLE")

    # (c2) the gate is insensitive to h, across a 4.4x span.
    print(f"\n  (c2) the SAME two bands agree within {BAR_C:.0f} % across "
          f"GL/GH/WLr, although h spans 38.4 -> 169.3 (4.4x):")
    ns = [n for n in ("GL", "GH", "WLr") if n in got]
    if len(ns) >= 2:
        for i, k in enumerate(("bore", "pads")):
            vs = [got[n][i] for n in ns]
            spread = 100.0 * (max(vs) - min(vs)) / abs(np.mean(vs)) if np.mean(vs) else float("nan")
            ok = abs(spread) <= BAR_C
            print(f"       {k:<5} " + "  ".join(f"{n} {v:8.1f}" for n, v in zip(ns, vs))
                  + f"   spread {spread:6.2f} %   {'PASS' if ok else 'FAIL'}")
            rep.setdefault("c", []).append((f"c2 {k}", ok))
        vs = [got[n][2] for n in ns]
        print(f"       wall  " + "  ".join(f"{n} {v:8.1f}" for n, v in zip(ns, vs))
              + "   <- this is the band that MAY move")

    # (c3) the positive control: the gate is not a no-op.
    print(f"\n  (c3) POSITIVE CONTROL -- WLr vs WLv at the SAME h = 81 and the")
    print(f"       same live T: bore and pads MUST differ by more than {BAR_C:.0f} %,")
    print(f"       or (c1)/(c2) would pass vacuously:")
    if "WLr" in got and "WLv" in got:
        db, dp = rel(got["WLr"][0], got["WLv"][0]), rel(got["WLr"][1], got["WLv"][1])
        ok = abs(db) > BAR_C and abs(dp) > BAR_C
        print(f"       bore {db:+8.1f} %   pads {dp:+8.1f} %   "
              f"{'PASS' if ok else 'FAIL'}")
        rep.setdefault("c", []).append(("c3 control", ok))


def crit_d(rep):
    hdr("(d) SCORED -- the power identity, +/- 0.02 W")
    print("  the one check that catches a kernel/march/ledger split, which is")
    print("  why the live temperature is resolved ONCE per step and all three")
    print("  sites read that member.")
    for n in NEW + ("S0",):
        if not have(n):
            continue
        d = thermo(n)
        if "jet_P_face" not in d:
            print(f"  {n:<5} no jet_P_face (no march)")
            continue
        s = (d["time"] >= WIN[0]) & (d["time"] <= WIN[1])
        e = float(np.abs((d["jet_P_face"] - (d["patch_P_robin"] + d["side_P_face"]
                                             + d["side_P_loss"]))[s]).max())
        ok = e <= 0.02
        print(f"  {n:<5} max residual {e:.3e} W   {'PASS' if ok else 'FAIL'}"
              f"    ledger_err max {np.abs(d['ledger_err']).max():.2e}")
        rep.setdefault("d", []).append((n, ok))


def _same(a, pa, b, pb, label, rep, key):
    ea, eb = os.path.exists(pa), os.path.exists(pb)
    if not (ea and eb):
        print(f"  {label}: missing ({a} {ea}, {b} {eb})")
        return
    ok = filecmp.cmp(pa, pb, shallow=False)
    if not ok:                      # report the prefix that does match
        la, lb = open(pa).readlines(), open(pb).readlines()
        n = sum(1 for x, y in zip(la, lb) if x == y)
        print(f"  {label}: {len(la)} vs {len(lb)} rows, matching prefix {n}  DIFFER")
    else:
        print(f"  {label}: IDENTICAL")
    rep.setdefault(key, []).append((label, ok))


def crit_e(rep):
    hdr("(e) SCORED -- identity")
    print("  S0 vs d2z1_delete/RB_D  (no new key at all):")
    for suf in ("_removal_events.csv", "_jet_profile.csv"):
        _same("S0", os.path.join(OUT, "S0" + suf),
              "RB_D", os.path.join(Z1, "RB_D" + suf), "    " + suf, rep, "e")
    hdr2 = lambda p: open(os.path.join(p, "thermo.dat")).readline().split()
    for a, pa, b, pb in (("S0", os.path.join(OUT, "S0"), "RB_D", os.path.join(Z1, "RB_D")),
                         ("SWc", os.path.join(OUT, "SWc"), "SW", os.path.join(Z2, "SW"))):
        if not (os.path.exists(pa) and os.path.exists(pb)):
            continue
        ha, hb = hdr2(pa), hdr2(pb)
        delta = set(ha) ^ set(hb)
        va = np.loadtxt(os.path.join(pa, "thermo.dat"), skiprows=1, ndmin=2)
        vb = np.loadtxt(os.path.join(pb, "thermo.dat"), skiprows=1, ndmin=2)
        nr = min(len(va), len(vb))
        diff = [c for c in set(ha) & set(hb)
                if not np.array_equal(va[:nr, ha.index(c)], vb[:nr, hb.index(c)])]
        ok = not delta and not diff
        print(f"\n  {a} vs {b}: thermo header delta "
              f"{'EMPTY' if not delta else sorted(delta)}; "
              f"{len(diff)} shared column(s) differ  {'PASS' if ok else 'FAIL'}")
        if diff:
            print(f"    differing: {sorted(diff)[:8]}")
        rep.setdefault("e", []).append((f"{a} vs {b} thermo", ok))
    _same("SWc", os.path.join(OUT, "SWc_removal_events.csv"),
          "SW", os.path.join(Z2, "SW_removal_events.csv"),
          "  SWc vs SW _removal_events.csv", rep, "e")
    print("\n  ⚠ CRITERION DEFECT, recorded not reinterpreted. CRITERION.md (e)")
    print("    says SWc's header delta 'must also be EMPTY'. That sentence is")
    print("    WRONG and contradicts the check made minutes earlier that put")
    print("    SWc in the campaign: frozen SW predates D2aa, so SWc necessarily")
    print("    GAINS D2aa's seven region columns. Scored as written (FAIL); the")
    print("    substantive reading is the one that matters -- 0 SHARED columns")
    print("    differ, so SWc reproduces frozen SW exactly on every column both")
    print("    carry, which is what re-validating the baseline required.")
    print("\n  ⚠ `metadata` and `diff.patch` are run provenance (compile date,")
    print("    plot_file path, the ParmParse table, a start stamp) and are")
    print("    excluded -- the same rule the 2639-file tree list applies.")
    print("  ⚠ SWc legitimately gains NO column: side_T_wall registers only")
    print("    when side_face_h is given, and SWc gives none.")
    print("\n  tree gate vs witness/d2z1.md5: run guarded_post_all.sh separately.")


def depth_profile(n):
    """Deliverable 3: rate against DEPTH in 50 mm stretches."""
    out = []
    for lo in (50, 100, 150, 200, 250, 300):
        t0, t1 = X.t_at(n, lo), X.t_at(n, lo + 50)
        out.append((lo, (0.050 / (t1 - t0) * 3600.0) if (t0 and t1) else None))
    return out


def crit_f(rep):
    hdr("(f) REPORTED -- Meier, with residuals stated not hidden")
    print("  ⚠ RATE IS NOT STEADY. A single number hides that, so it is")
    print("    reported against DEPTH in 50 mm stretches (deliverable 3).")
    print("    The SCORED metric stays the 100-150 mm descent stretch.\n")
    legs = [n for n in NEW + ("WL", "WL81", "SW", "A", "D") if have(n)]
    print(f"  {'leg':<6}{'ROP':>8}  " + "".join(f"{lo:>4}-{lo+50:<4}" for lo in
                                                (50, 100, 150, 200, 250, 300)))
    print(f"  {'':6}{'[m/h]':>8}  " + "".join(f"{'[m/h]':>9}" for _ in range(6)))
    for n in legs:
        pr = depth_profile(n)
        print(f"  {n:<6}{X.rop(n):8.3f}  "
              + "".join(f"{(v if v else float('nan')):9.2f}" for _, v in pr))
    print(f"\n  shape at equal burner descent ({X.PRIMARY} mm):")
    print(f"  {'leg':<6}{'mean|dD|':>10}{'D@100':>8}{'D@125':>8}{'pit lead':>10}"
          f"{'vol':>7}{'x_pit':>7}{'x_desc':>8}")
    for n in legs:
        s = X.shape(n, X.PRIMARY, (10, 25, 50, 75, 100, 125))
        if s is None:
            print(f"  {n:<6}   never reached {X.PRIMARY} mm descent")
            continue
        print(f"  {n:<6}{(s['a'] or float('nan')):10.2f}"
              f"{s['D'].get(100, float('nan')):8.1f}"
              f"{s['D'].get(125, float('nan')):8.1f}{s['lead']:10.1f}"
              f"{s['vol']:7.2f}{s['x_pit']:7.2f}{s['x_desc']:8.2f}")
    print(f"\n  rings, 2 mm, {X.RINGS[0]}-{X.RINGS[-1]+2} mm [mm recession], "
          f"PER RING, never averaged:")
    print("        " + "".join(f"{r:7d}" for r in X.RINGS))
    for n in legs:
        print(f"  {n:<5} " + "".join(f"{v:7.1f}" for v in X.rings(n)))
    print("\n  ⚠ Meier Fig. 8.8 is a VISIBLE WIDTH = a LOWER BOUND on diameter,")
    print("    band 85-95 +/- 5 mm. No agreement claim tighter than +/- 5 mm.")
    print("  ⚠ Meier's rate is ONE-SIDED (operator advance).")


def specific_energy(rep):
    hdr("SPECIFIC ENERGY and the power scaling")
    print(f"  {'leg':<6}{'P_rock kW':>11}{'P_top kW':>10}{'P_side kW':>11}"
          f"{'vol cm3':>9}{'J/mm3':>8}{'A_side cm2':>12}")
    for n in [x for x in NEW + ("WL", "WL81", "SW", "A") if have(x)]:
        d = thermo(n)
        pf = tm(d, "jet_P_face") / 1e3
        ps = tm(d, "side_P_face") / 1e3
        s = X.shape(n, X.PRIMARY, (100,))
        asd = sum(tm(d, k) / 4.0 * 1e4 for k in ("side_A_bore", "side_A_pads",
                                                 "side_A_wall")) if \
            "side_A_bore" in d else float("nan")
        v = s["vol"] if s else float("nan")
        # ⚠ `vol` is in LITRES (the d2x basis: x_desc = vol/(V_L*d/V_M)), and
        # specific energy needs a volumetric RATE, so it is
        #     J/mm3 = P[W] * t_150[s] / (vol[L] * 1e6 mm3/L).
        # An earlier revision printed P/vol, which is W/mm3 -- off by the time
        # basis and by 1e3. Fixed here; the corrected numbers are the ones in
        # RESULTS.md.
        t150 = s["t"] if s else float("nan")
        je = (pf * 1e3 * t150 / (v * 1e6)) if (v == v and v) else float("nan")
        print(f"  {n:<6}{pf:11.2f}{pf-ps:10.2f}{ps:11.2f}{v:9.2f}"
              f"{je:8.2f}{asd:12.2f}")
    print("\n  ⚠ P_rock is the FULL-hole jet_P_face; J/mm3 is P * t_150 / V at")
    print("    equal burner descent, so it is a rate-independent ratio.")


def predictions(rep):
    hdr("PREDICTIONS")
    R = rep.get("rate", {})

    def say(k, txt, ok):
        print(f"  {k} {txt}: {'HELD' if ok else 'REFUTED'}")

    for n in ("GL", "GH"):
        if n in R:
            d = rel(R[n], PRED[n]["rate"])
            say("P1", f"{n} within 15 % of {PRED[n]['rate']:.2f} "
                      f"(got {R[n]:.3f}, {d:+.1f} %)", abs(d) <= 15)
    for n in ("GL", "GH"):
        if n in R:
            say("P2", f"(the PACKET's) {n} in {P2_BAND[0]}-{P2_BAND[1]} m/h "
                      f"(got {R[n]:.3f})", P2_BAND[0] <= R[n] <= P2_BAND[1])
    if "WLv" in R:
        say("P3", f"WLv in {WLV_BAND[0]}-{WLV_BAND[1]} (got {R['WLv']:.3f})",
            WLV_BAND[0] <= R["WLv"] <= WLV_BAND[1])
    if "WLv" in R and "WLr" in R:
        say("P3", f"WLr above WLv ({R['WLr']:.3f} vs {R['WLv']:.3f})",
            R["WLr"] > R["WLv"])
    if all(k in R for k in ("GL", "GH", "WLr")):
        lo, hi = min(R["GL"], R["GH"]), max(R["GL"], R["GH"])
        say("P4", f"the derived bracket contains 81, so GL/GH bracket WLr "
                  f"({lo:.3f} <= {R['WLr']:.3f} <= {hi:.3f})",
            lo <= R["WLr"] <= hi)
        print("       (on paper at pre-flight: 38.4 < 81 < 169.3, already true)")
    # P5: the charged riser area keeps growing as h falls.
    print("\n  P5 charged riser area vs h (D2aa: 68 -> 133 -> 143 cm2):")
    pts = []
    for n in ("GH", "WLr", "GL"):
        if not have(n):
            continue
        d = thermo(n)
        m = ast.literal_eval(open(os.path.join(OUT, n + ".done")).read())
        a = sum(tm(d, k) / 4.0 * 1e4 for k in ("side_A_bore", "side_A_pads",
                                               "side_A_wall"))
        pts.append((float(m["side_face_h"]), a, n))
        print(f"       {n:<5} h = {m['side_face_h']:>7}  charged area {a:8.2f} cm2")
    if len(pts) >= 2:
        pts.sort(reverse=True)         # descending h
        grow = all(b[1] >= a[1] for a, b in zip(pts, pts[1:]))
        say("P5", "the area keeps GROWING as h falls (refuted = it saturates)",
            grow)
    # P6: the live temperature settles.
    print("\n  P6 the live exhaust settles (std < 2 % of mean, no growing "
          "oscillation):")
    for n in ("GL", "GH", "WLr", "WLv"):
        if not have(n):
            continue
        d = thermo(n)
        s = (d["time"] >= WIN[0]) & (d["time"] <= WIN[1])
        v = d["jet_T_exhaust"][s]
        cv = 100.0 * v.std() / v.mean()
        h = len(v) // 2
        drift = rel(float(v[h:].mean()), float(v[:h].mean()))
        ok = cv < 2.0
        print(f"       {n:<5} T_exh {v.mean():7.1f} K  std {cv:5.2f} %  "
              f"2nd-half drift {drift:+5.2f} %   {'HELD' if ok else 'REFUTED'}")
        rep.setdefault("p6", []).append((n, ok))
    print("       ⚠ Refuted ⇒ the feedback loop needs damping before any leg")
    print("         is scored (less heat to the wall ⇒ hotter gas ⇒ more heat).")
    # P7: no leg reaches Meier's band.
    ns = [n for n in NEW if n in R]
    inband = [n for n in ns if BAND[0] <= R[n] <= BAND[1]]
    print()
    say("P7", f"NO leg reaches Meier's {BAND[0]}-{BAND[1]} m/h band "
              f"(rates: {', '.join(f'{n} {R[n]:.3f}' for n in ns)})", not inband)
    if not inband:
        print("       ⇒ the residual is the BORE and PAD-RING lateral faces, in")
        print("         impinging flow, which this model cannot treat (PD jammed")
        print("         on lateral recession). The successor is LATERAL")
        print("         RECESSION, not another coefficient.")


def main():
    print(f"D2ab -- scoring CRITERION.md (sha from CRITERION.sha256)")
    print(f"legs on disk: "
          + ", ".join(n for n in NEW + ("S0",) if have(n)))
    rep = {}
    crit_a(rep)
    crit_b(rep)
    crit_c(rep)
    crit_d(rep)
    crit_e(rep)
    crit_f(rep)
    specific_energy(rep)
    predictions(rep)
    hdr("SCORECARD")
    for k in ("a", "b", "c", "d", "e"):
        rows = rep.get(k, [])
        if not rows:
            continue
        npass = sum(1 for _, ok in rows if ok)
        print(f"  ({k}) {npass}/{len(rows)} pass"
              + ("" if npass == len(rows) else
                 "   FAILED: " + ", ".join(str(n) for n, ok in rows if not ok)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
