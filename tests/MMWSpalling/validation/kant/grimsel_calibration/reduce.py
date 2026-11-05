#!/usr/bin/env python3
"""P2-13 -- reduce the Grimsel calibration legs.

Usage: reduce.py [A1 | A2 | A3 | B1 | D1]     (no mode: every row -> stdout, i.e. RESULTS.txt)

THE STATISTIC (packet): per column (col_i, col_j), the first event-log row with dir = 5 (top face) and
ts_fb = 0; its T_surf (the reconstructed surface temperature the criterion fired on). Reported: n, median,
p10, p90; the count of columns whose FIRST firing was not a top face or was a fallback; the all-firing median.
⚠ The advisor's G1 number (842.74 K) is a DIFFERENT statistic: the first non-fallback firing per column in ANY
direction. On G1 most columns are first reached from the side as the bore widens, so the two differ there;
on this flat deck they coincide when every first firing is a top face (b1). Both are printed.
"""
import csv
import math
import os
import re
import shlex
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", "..", ".."))
LEGS = os.path.join(HERE, "legs")
PY = sys.executable
CHECK = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "shared", "check_shared.py")
G1 = os.path.join(ROOT, "tests", "MMWSpalling", "studies", "p2_11_heat", "legs", "G1")
LAB = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "hoeser", "input_lab")
BAND = (826.15, 873.15)
BANKED, G1_ALL, P22F = 843.757, 851.55, 852.61
TIMERS = {"fr_t_removal", "fr_t_step", "fl_t_path"}
CAL_LEGS = ["CAL", "CAL_det", "CAL_h2000", "CAL_p27", "CAL_p48", "CAL_2mm", "CAL_05mm"]


def events(path):
    return list(csv.DictReader(open(path)))


def stat(rows):
    """the packet's statistic plus the advisor's variant; rows in file (time) order"""
    first, first_q, first_any_nofb = {}, {}, {}
    for r in rows:
        c = (int(r["col_i"]), int(r["col_j"]))
        first.setdefault(c, r)
        if float(r["ts_fb"]) == 0.0:
            first_any_nofb.setdefault(c, r)
            if r["dir"] == "5":
                first_q.setdefault(c, r)
    T = lambda d: np.array([float(r["T_surf"]) for r in d.values()])
    q, adv, al = T(first_q), T(first_any_nofb), np.array([float(r["T_surf"]) for r in rows])
    bad_first = sum(1 for r in first.values() if not (r["dir"] == "5" and float(r["ts_fb"]) == 0.0))
    s = dict(n=len(q), cols_fired=len(first), not_top_first=bad_first, n_rows=len(rows))
    for nm, x in (("", q), ("adv_", adv), ("all_", al)):
        s[nm + "n"] = len(x)
        s[nm + "med"] = float(np.median(x)) if len(x) else float("nan")
        s[nm + "p10"] = float(np.percentile(x, 10)) if len(x) else float("nan")
        s[nm + "p90"] = float(np.percentile(x, 90)) if len(x) else float("nan")
    s["above_clamp"] = float(np.mean(al > 873.15)) if len(al) else float("nan")
    s["first"] = first
    return s


def fmt(s, pre=""):
    return "n %4d  median %8.2f K  p10 %7.2f  p90 %7.2f" % (s[pre + "n"], s[pre + "med"], s[pre + "p10"], s[pre + "p90"])


def meta(d):
    m = {}
    p = os.path.join(d, "run", "metadata")
    for ln in open(p):
        if " = " in ln:
            k, v = ln.rstrip("\n").split(" = ", 1)
            m[k.strip()] = v.strip()
    return m


def thermo(d):
    L = open(os.path.join(d, "run", "thermo.dat")).read().splitlines()
    h = [c.strip() for c in L[0].split("\t")]
    return h, [r.split("\t") for r in L[1:] if r.strip()]


def leg(name):
    return os.path.join(LEGS, name)


def leg_stat(name):
    return stat(events(os.path.join(leg(name), "run_removal_events.csv")))


# ------------------------------------------------------------------------- A1
def A1():
    s = stat(events(os.path.join(G1, "run_removal_events.csv")))
    L = ["a1 -- G1's event log by the statistic (zero compute): %s" % os.path.relpath(G1, ROOT), "=" * 100,
         "  packet statistic (first dir=5, ts_fb=0 row per column):  " + fmt(s),
         "  advisor's (first non-fallback firing, ANY direction):     " + fmt(s, "adv_"),
         "  all firings:                                             " + fmt(s, "all_"),
         "  columns fired %d; first firing NOT a non-fallback top face: %d (%.1f %%)"
         % (s["cols_fired"], s["not_top_first"], 100.0 * s["not_top_first"] / s["cols_fired"])]
    ok_adv = abs(s["adv_med"] - 842.74) < 0.005 and s["adv_n"] == 2076
    ok_all = abs(s["all_med"] - 851.55) < 0.005
    L.append("  Sources: first-firing 842.74 K over 2,076 columns -> reproduced by the ADVISOR's statistic: %s (%.3f)"
             % ("yes" if ok_adv else "NO", s["adv_med"]))
    L.append("           all-firing 851.55 K -> %s (%.3f)" % ("yes" if ok_all else "NO", s["all_med"]))
    L.append("  ⚠ the packet's own statistic gives %.2f K on G1 (%d columns): G1's columns are mostly first reached from"
             % (s["med"], s["n"]))
    L.append("    the side as the bore widens. On the flat calibration deck the two coincide iff b1's top-face count holds.")
    L.append("a1 %s (report): the reducer reads the log as the advisor did; the statistic named in the packet is a"
             " different column subset on G1" % ("REPRODUCED" if ok_adv and ok_all else "NOT REPRODUCED"))
    return "\n".join(L)


# ------------------------------------------------------------------------- A2
def check(deck, argv, allow, run=None):
    a = [PY, CHECK, deck, argv] + (["--allow"] + allow if allow else []) + (["--run", run] if run else [])
    p = subprocess.run(a, cwd=ROOT, capture_output=True, text=True)
    return p.returncode, p.stdout.rstrip()


def allow_of(d):
    p = os.path.join(d, "ALLOW.txt")
    return open(p).read().split() if os.path.exists(p) else []


def A2():
    L = ["a2 -- check_shared on CAL (a 1-step run, legs/SMOKE) and on HOS_R (input_lab re-pointed)", "=" * 100]
    ok = True
    for nm, deck in (("SMOKE", os.path.join(HERE, "input")), ("HOS_R", LAB)):
        rc, out = check(deck, os.path.join(leg(nm), "argv.sh"), allow_of(leg(nm)))
        L += ["  " + x for x in out.splitlines()]
        ok = ok and rc == 0
    # non-vacuity: an undeclared override MUST fail (i); a declared one must pass
    fake = os.path.join(leg("SMOKE"), "argv_control.sh")
    open(fake, "w").write(open(os.path.join(leg("SMOKE"), "argv.sh")).read().rstrip() + " weibull.m=20\n")
    rc_c, out_c = check(os.path.join(HERE, "input"), fake, [], os.path.join(leg("SMOKE"), "run"))
    ctrl = rc_c != 0 and "(i)   FAIL" in out_c
    L.append("  control: SMOKE's argv + an undeclared weibull.m=20 -> (i) %s" % ("FAILS, as it must" if ctrl else "PASSED: the check is vacuous"))
    ok = ok and ctrl
    L.append("a2 %s (gate)" % ("PASS" if ok else "FAIL"))
    return "\n".join(L), ok


# ------------------------------------------------------------------------- A3
def md5(p):
    import hashlib
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def A3():
    O, d = G1, leg("HOS_R")
    L = ["a3 -- HOS_R (input_lab re-pointed + G1's argv minus the shared keys) vs G1, to 15 s", "=" * 100]
    rc = int(open(os.path.join(d, "rc")).read())
    L.append("  rc %d; binary %s (G1 ran on P2-11's binary; same default-path code, P2-12 r2/r5)"
             % (rc, open(os.path.join(d, "BINARY.txt")).read().split()[0][:16]))
    ok = rc == 0
    for pf in ("00000cell", "00000node", "01875cell", "01875node"):
        a, b = os.path.join(O, "run", pf), os.path.join(d, "run", pf)
        if not os.path.isdir(b):
            L.append("  %s: MISSING" % pf); ok = False; continue
        fa = {os.path.relpath(os.path.join(p, f), a): os.path.join(p, f) for p, _, fs in os.walk(a) for f in fs}
        fb = {os.path.relpath(os.path.join(p, f), b): os.path.join(p, f) for p, _, fs in os.walk(b) for f in fs}
        diff = [k for k in set(fa) & set(fb) if md5(fa[k]) != md5(fb[k])]
        only = sorted(set(fa) ^ set(fb))
        L.append("  plotfile %s: %d files, %d differ, %d only on one side" % (pf, len(set(fa) & set(fb)), len(diff), len(only)))
        ok = ok and not diff and not only
    ha, ra = thermo(os.path.join(O))
    hb, rb = thermo(d)
    ra = [r for r in ra if float(r[0]) <= 15.0 + 1e-9]
    rb = [r for r in rb if float(r[0]) <= 15.0 + 1e-9]
    if ha != hb or len(ra) != len(rb):
        L.append("  thermo: header same %s, rows %d vs %d" % (ha == hb, len(ra), len(rb))); ok = False
    else:
        idx = [i for i, c in enumerate(ha) if c not in TIMERS]
        nbad = sum(1 for x, y in zip(ra, rb) if [x[i] for i in idx] != [y[i] for i in idx])
        L.append("  thermo: %d rows to 15 s, %d non-timer columns, %d rows differ" % (len(ra), len(idx), nbad))
        ok = ok and nbad == 0
    for nm in ("run_removal_events.csv", "run_jet_profile.csv"):
        A = open(os.path.join(O, nm)).read().splitlines()
        B = open(os.path.join(d, nm)).read().splitlines()
        pre = A[:len(B)] == B
        tcol = A[0].split(",").index("time")
        n15 = 1 + sum(1 for r in A[1:] if float(r.split(",")[tcol]) < 15.0 - 1e-9)
        L.append("  %s: %d lines, exact prefix of G1's %s; G1 has %d lines with time < 15 s"
                 % (nm, len(B), "YES" if pre else "NO", n15))
        ok = ok and pre and len(B) >= n15 - 0
    L.append("a3 %s (gate)" % ("PASS" if ok else "FAIL"))
    return "\n".join(L), ok


# ------------------------------------------------------------------------- B1
def b1_leg(name):
    d = leg(name)
    rc = int(open(os.path.join(d, "rc")).read())
    m = meta(d)
    n = [int(x) for x in m["amr.n_cell"].split()]
    hi = [float(x) for x in m["geometry.prob_hi"].split()]
    dz = hi[2] / n[2]
    s = leg_stat(name)
    ncol = n[0] * n[1]
    share = s["n"] / ncol
    # the floor at stop: per column, the cells voided by all its removals (summed n_voided). Conservative
    # against SR 4's k_top - n_voided reading, which gave 1 mm more on CAL_05mm (P2-13 notes).
    vox = {}
    for r in events(os.path.join(d, "run_removal_events.csv")):
        c = (int(r["col_i"]), int(r["col_j"]))
        vox[c] = vox.get(c, 0) + int(r["n_voided"])
    floor = (n[2] - max(vox.values())) * dz if vox else hi[2]
    h, rows = thermo(d)
    j = h.index("ledger_err_solid")
    le = max(abs(float(r[j])) for r in rows)
    rcs, out = check(os.path.join(HERE, "input"), os.path.join(d, "argv.sh"), allow_of(d))
    ok = rc == 0 and share >= 0.95 and le <= 5e-3 and floor >= 0.010 and rcs == 0
    return ok, s, ("  %-10s rc %d  top-first columns %4d/%4d (%.1f %%)  not-top-first %d  max|ledger_err_solid| %.2e"
                   "  lowest surface %.1f mm above base  check_shared %s  -> %s"
                   % (name, rc, s["n"], ncol, 100 * share, s["not_top_first"], le, floor * 1e3,
                      "PASS" if rcs == 0 else "FAIL", "PASS" if ok else "FAIL"))


def B1(names=None):
    names = names or [n for n in CAL_LEGS if os.path.exists(os.path.join(leg(n), ".done"))]
    L = ["b1 -- each leg: rc 0; >= 95 % columns with a first top-face non-fallback firing; ledger_err_solid <= 5e-3;"
         " floor >= 10 mm above the base at stop; plus check_shared on the leg's own argv", "=" * 100]
    S, ok = {}, True
    for n in names:
        o, s, ln = b1_leg(n)
        S[n] = s; ok = ok and o; L.append(ln)
    missing = [n for n in CAL_LEGS if n not in names]
    if missing:
        L.append("  NOT RUN: " + ", ".join(missing)); ok = False
    L.append("b1 %s (gate)" % ("PASS" if ok else "FAIL"))
    return "\n".join(L), ok, S


# ------------------------------------------------------------------------- D1
def interp_a(pts, target=850.0):
    """log-log interpolation of a0 vs onset median to the target; pts = [(a0, T)]"""
    (a1, t1), (a2, t2) = pts
    x = (math.log(target) - math.log(t1)) / (math.log(t2) - math.log(t1))
    return math.exp(math.log(a1) + x * (math.log(a2) - math.log(a1)))


def D1():
    L = ["d1 -- recalibrate the flaw scale once (only because c1 failed)", "=" * 100]
    s6, s10 = leg_stat("RECAL_6"), leg_stat("RECAL_10")
    L.append("  RECAL_6  a0 6.0 um : " + fmt(s6))
    L.append("  RECAL_10 a0 10.0 um: " + fmt(s10))
    br = min(s6["med"], s10["med"]) < 850.0 < max(s6["med"], s10["med"])
    L.append("  bracket 850 K: %s" % ("yes" if br else "NO"))
    if not br:
        L.append("d1 FAIL (gate): the level cannot be reached with one flaw scale")
        return "\n".join(L), False
    a = interp_a([(6.0e-6, s6["med"]), (10.0e-6, s10["med"])])
    a4 = float("%.4g" % a)
    open(os.path.join(HERE, "RECAL_X.txt"), "w").write("%.4g\n" % a4)
    L.append("  log-log interpolated a0 = %.4g m -> RECAL_X.txt" % a4)
    if os.path.exists(os.path.join(leg("RECAL_x"), ".done")):
        sx = leg_stat("RECAL_x")
        ok = abs(sx["med"] - 850.0) <= 3.0
        L.append("  RECAL_x: " + fmt(sx) + "  |median - 850| = %.2f K -> %s" % (abs(sx["med"] - 850.0), "PASS" if ok else "FAIL"))
        L.append("d1 %s (gate)" % ("PASS" if ok else "FAIL"))
        return "\n".join(L), ok
    L.append("d1 PENDING: run RECAL_x")
    return "\n".join(L), None


# ------------------------------------------------------------------------- RESULTS
def RESULTS():
    out = []
    for f in ("A0.txt", "A1.txt", "A2.txt", "A3.txt", "UNITS_p2_13_vs_p2_12r.txt"):
        p = os.path.join(HERE, f)
        out.append(open(p).read().rstrip() if os.path.exists(p) else "%s: NOT ON DISK" % f)
        out.append("")
    txt, ok_b1, S = B1()
    out += [txt, ""]
    L = ["c -- the onset (the packet statistic, first top-face non-fallback firing per column)", "=" * 100]
    for n in CAL_LEGS:
        if n in S:
            L.append("  %-10s %s   all-firing median %.2f K" % (n, fmt(S[n]), S[n]["all_med"]))
    c = S.get("CAL")
    verdict_c1 = None
    if c:
        inside = BAND[0] <= c["med"] <= BAND[1]
        verdict_c1 = inside
        L.append("")
        L.append("c1 THE LEVEL: CAL first-firing median %.2f K vs the Grimsel onset band %.2f-%.2f K -> %s (gate)"
                 % (c["med"], BAND[0], BAND[1], "INSIDE (PASS)" if inside else "OUTSIDE (FAIL: d1 runs)"))
        L.append("   P1 predicted 845-856 K: %s" % ("held" if 845 <= c["med"] <= 856 else "MISSED"))
        L.append("   (advisor's any-direction statistic on CAL: %.2f K; identical iff every first firing is a top face)" % c["adv_med"])
        if "CAL_det" in S:
            dd = S["CAL_det"]["med"] - BANKED
            L.append("c2 CAL_det median %.2f K - banked %.3f = %+.2f K   (P2: within +-5 K: %s)"
                     % (S["CAL_det"]["med"], BANKED, dd, "held" if abs(dd) <= 5 else "MISSED"))
            d3 = c["med"] - S["CAL_det"]["med"]
            L.append("c3 CAL - CAL_det = %+.2f K   (what the Weibull draw moves; P3 +4..+9 K, twin +6.3: %s)"
                     % (d3, "held" if 4 <= d3 <= 9 else "MISSED"))
        if "CAL_h2000" in S:
            th = S["CAL_h2000"]["med"]
            d4 = th - c["med"]
            q1, q2 = 500.0 * (1673.15 - c["med"]), 2000.0 * (1673.15 - th)
            e = math.log((th - 293.15) / (c["med"] - 293.15)) / math.log(q2 / q1)
            L.append("c4 CAL_h2000 - CAL = %+.2f K; onset flux q = h(T_flame - T_onset): %.0f -> %.0f kW/m2 (x%.2f);"
                     " implied exponent dln(dT)/dln(q) = %+.4f   (P4 +10..+40 K: %s; above +25 K is a finding against"
                     " data that show no trend)" % (d4, q1 / 1e3, q2 / 1e3, q2 / q1, e, "held" if 10 <= d4 <= 40 else "MISSED"))
        for nm, kant, kmod, lo, hi in (("CAL_p27", -55, -25, -28, -20), ("CAL_p48", -70, -45, -42, -30)):
            if nm in S:
                d = S[nm]["med"] - c["med"]
                L.append("c5 %s - CAL = %+.2f K   (Kant measured ~%d K [fig], his model %d K; P5 %d..%d K: %s)%s"
                         % (nm, d, kant, kmod, lo, hi, "held" if lo <= d <= hi else "MISSED",
                            "  ⚠ ZERO: L5 FAILED" if abs(d) < 0.01 else ""))
        ms = [(n, S[n]["med"]) for n in ("CAL_2mm", "CAL", "CAL_05mm") if n in S]
        if len(ms) == 3:
            sp = max(v for _, v in ms) - min(v for _, v in ms)
            L.append("c6 mesh: 2 mm %.2f / 1 mm %.2f / 0.5 mm %.2f K; spread %.2f K = %.2f %%   (P6 within +-8 K: %s)"
                     % (ms[0][1], ms[1][1], ms[2][1], sp, 100 * sp / c["med"], "held" if all(abs(v - c["med"]) <= 8 for _, v in ms) else "MISSED"))
        d7a, d7b = c["all_med"] - G1_ALL, c["all_med"] - P22F
        L.append("c7 CAL all-firing median %.2f K: vs G1 %.2f %+.2f K, vs P2-2f Meier %.2f %+.2f K   (P7 within +-10 K: %s)"
                 % (c["all_med"], G1_ALL, d7a, P22F, d7b, "held" if abs(d7a) <= 10 and abs(d7b) <= 10 else "MISSED"))
        L.append("c8 share of CAL firings above 873.15 K (the table clamp): %.1f %% of %d" % (100 * c["above_clamp"], c["n_rows"]))
    out += L
    out.append("")
    if verdict_c1 is False:
        txt, ok = D1()
        out += [txt, ""]
    out.append("e1 VERDICT inputs: c1 %s; b1 %s" % ({None: "not run", True: "INSIDE", False: "OUTSIDE"}[verdict_c1],
                                                    "PASS" if ok_b1 else "FAIL"))
    if verdict_c1 is True and ok_b1:
        out += ["",
                "e1 VERDICT: THE ONE GRIMSEL CALIBRATION STANDS UNCHANGED.",
                "   R 0.531294, flaw scale 8 um (weibull.a0_ig = spallation.sp.a0), m 13, Grimsel laws and Table 3:",
                "   in the production form (surface Robin, losses, Weibull, flake rule, face list) the Kant onset is",
                "   %.2f K, inside the Grimsel band; d1 was not run, so G1's S* needs no re-derivation." % c["med"],
                "   Files: validation/shared/rock_grimsel.inputs + model_form_flame.inputs, included by",
                "   validation/kant/grimsel_calibration/input and validation/hoeser/input_lab (hence G1's deck and",
                "   unit/face_heat). NOT yet included: the legacy Kant inputs (sp_kant_onset*), the Meier decks and",
                "   every tree input (untouched by the guardrails; their retirement needs the user's permission).",
                "   The merge step carries: d_ref (refused by the face list; its 2.0 mm was calibrated on the column",
                "   path, Kant 2016 Gotthard properties) and the Meier deck, which must be rebuilt to include both files."]
    return "\n".join(out)


def main(a):
    mode = a[0] if a else "RESULTS"
    if mode == "A1":
        t = A1(); print(t); open(os.path.join(HERE, "A1.txt"), "w").write(t + "\n")
    elif mode == "A2":
        t, ok = A2(); print(t); open(os.path.join(HERE, "A2.txt"), "w").write(t + "\n"); return 0 if ok else 1
    elif mode == "A3":
        t, ok = A3(); print(t); open(os.path.join(HERE, "A3.txt"), "w").write(t + "\n"); return 0 if ok else 1
    elif mode == "B1":
        t, ok, _ = B1(); print(t); return 0 if ok else 1
    elif mode == "D1":
        t, ok = D1(); print(t); open(os.path.join(HERE, "D1.txt"), "w").write(t + "\n"); return 0 if ok else 1
    else:
        print(RESULTS())
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
