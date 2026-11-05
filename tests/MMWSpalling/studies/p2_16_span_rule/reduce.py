#!/usr/bin/env python3
"""P2-16 -- every row. Usage: reduce.py [A0|A2|A3|A4|C1] (no mode: RESULTS.txt to stdout).

Helpers are IMPORTED: P2-15's reduce.py (events, first_top, stat, fmt, check_leg logic, mask, the discriminant through
rules.py), P2-14a's (same_file, thermo_same), P2-10's run_p2_10.py (verdict_lines). S_top (SR 19): per column the
first event-log row with dir = 5 and ts_fb = 0, its T_surf; the level is the median over columns.
"""
import importlib.machinery
import importlib.util
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
T = os.path.join(ROOT, "tests", "MMWSpalling")
LEGS = os.path.join(HERE, "legs")
P15 = os.path.join(T, "studies", "p2_15_buckle")
CALD = os.path.join(T, "validation", "kant", "grimsel_calibration")
CHECK = os.path.join(T, "validation", "shared", "check_shared.py")
PY = sys.executable


def load(name, path):
    ld = importlib.machinery.SourceFileLoader(name, path)
    m = importlib.util.module_from_spec(importlib.util.spec_from_loader(name, ld))
    ld.exec_module(m)
    return m


RU = load("rules16", os.path.join(HERE, "rules.py"))
R15 = RU.R15
R14 = R15.R14
REFS = {"CAL": (os.path.join(CALD, "legs", "CAL"), "P2-13"),
        "HOS_R": (os.path.join(CALD, "legs", "HOS_R"), "P2-13"),
        "F0_k1": (os.path.join(T, "retired_units", "face_cref", "output", "F0_k1"), "P2-15 binary (face_cref's run)"),
        "KB150_0": (os.path.join(P15, "legs", "KB150_0"), "P2-15")}


def rc(d):
    p = os.path.join(d, "rc")
    return int(open(p).read()) if os.path.exists(p) else None


def check_leg(leg):
    d = os.path.join(LEGS, leg)
    deck = open(os.path.join(d, "DECK.txt")).read().strip()
    allow = open(os.path.join(d, "ALLOW.txt")).read().split()
    a = [PY, CHECK, deck, os.path.join(d, "argv.sh")] + (["--allow"] + allow if allow else [])
    p = subprocess.run(a, cwd=ROOT, capture_output=True, text=True)
    return p.returncode == 0, (p.stdout.strip().splitlines() or ["(no output)"])[-1]


def kept_plotfiles(ref, new):
    """every plotfile the NEW leg kept (the guardrail prunes all but the last except CAL/HOS_R) must exist in the
    reference and be byte-identical: (files compared, differing, missing on either side)"""
    n = nd = no = 0
    for pf in R14.pf_dirs(new):
        A, B = os.path.join(ref, "run", pf), os.path.join(new, "run", pf)
        if not os.path.isdir(A):
            no += 1
            continue
        fa = {os.path.relpath(os.path.join(p, f), A): os.path.join(p, f) for p, _, fs in os.walk(A) for f in fs}
        fb = {os.path.relpath(os.path.join(p, f), B): os.path.join(p, f) for p, _, fs in os.walk(B) for f in fs}
        no += len(set(fa) ^ set(fb))
        for k in set(fa) & set(fb):
            n += 1
            nd += R15.md5(fa[k]) != R15.md5(fb[k])
    return n, nd, no


def logtail(p, k=1):
    return open(p).read().rstrip().splitlines()[-k:] if os.path.exists(p) else ["%s: NOT ON DISK" % p]


# ---------------------------------------------------------------- a0
def A0():
    L = ["a0 -- source half of unit/retired_keys (+ mutations, --pre src_pre/), lateral_flux_lock's restated clause and "
         "its mutation control", "=" * 100]
    ok = True
    p = os.path.join(HERE, "A0_retired_keys.log")
    t = open(p).read() if os.path.exists(p) else ""
    L += ["  " + x for x in t.splitlines() if x.startswith(("  ok", "  FAIL", "  RED", "  GRN", "PASS", "FAIL"))]
    ok = ok and t.rstrip().splitlines()[-1].startswith("PASS") if t else False
    p = os.path.join(HERE, "A0_lfl.log")
    t = open(p).read() if os.path.exists(p) else ""
    keep = ("STORED field", "law-path twins", "T0_at_didx, unconditionally", "no reference resolver")
    L += ["  lateral_flux_lock: " + x.strip() for x in t.splitlines() if any(k in x for k in keep)]
    lfl = bool(t) and t.rstrip().splitlines()[-1].startswith(("PASS", "      in src/")) and "  FAIL" not in t
    L.append("  lateral_flux_lock overall: %s" % ("PASS" if lfl else "FAIL"))
    ok = ok and lfl
    c = os.path.join(HERE, "lfl_control", "CONTROL.txt")
    ct = open(c).read() if os.path.exists(c) else ""
    hold = "control HOLDS" in ct
    L.append("  lateral_flux_lock mutation control (lfl_control/CONTROL.txt): %s" % ("HOLDS (the restated clause goes "
             "RED, nothing else)" if hold else "BROKEN"))
    ok = ok and hold
    L.append("a0 %s (gate)" % ("PASS" if ok else "FAIL"))
    return "\n".join(L), ok


# ---------------------------------------------------------------- a2
RESTATED = ("subtract their resolved reference", "references come from the ONE resolver", "`Tref_at` returns T0_at",
            "subtract the STORED field T0_at", "law-path twins pass the same", "T0_at_didx, unconditionally",
            "no reference resolver or constraint-reference name")
RETIRED = ("constraint_ref_lock", "face_cref")


def A2():
    R10 = load("r10q16", os.path.join(T, "studies", "p2_10_retire", "run_p2_10.py"))
    A, B = os.path.join(P15, "UNIT_p2_15"), os.path.join(HERE, "UNIT_p2_16")
    rcm = lambda d: {l.split()[0]: l.split()[1] for l in open(os.path.join(d, "rc.txt")).read().splitlines()[1:]}
    ra, rb = rcm(A), rcm(B)
    L = ["a2 -- every kept unit on the P2-16 binary vs P2-15's UNIT_p2_15 (source line numbers masked): rc and verdict "
         "lines; lateral_flux_lock by full log", "=" * 100]
    bad = []
    for n in sorted(set(ra) | set(rb)):
        if n in RETIRED:
            ok = n not in rb and os.path.isdir(os.path.join(T, "retired_units", n))
            L.append("  %-28s rc %5s -> (retired: moved to tests/MMWSpalling/retired_units/)  %s" % (
                n, ra.get(n, "-"), "ok" if ok else "STILL RUN"))
        elif n == "retired_keys":
            log = open(os.path.join(B, n + ".log")).read()
            ok = n not in ra and rb.get(n) == "rc=0" and re.search(r"^PASS", log, re.M) is not None
            L.append("  %-28s (new)  rc %5s   PASS: %s" % (n, rb.get(n, "-"), ok))
        else:
            vl = lambda d: [R15.mask(x) for x in R10.verdict_lines(os.path.join(d, n + ".log"))]
            va = vl(A) if n in ra else None
            vb = vl(B) if n in rb else None
            ok = ra.get(n) == rb.get(n) and va == vb
            extra = ""
            if n == "lateral_flux_lock":
                fa = [R15.mask(x) for x in open(os.path.join(A, n + ".log")).read().splitlines()]
                fb = [R15.mask(x) for x in open(os.path.join(B, n + ".log")).read().splitlines()]
                da = [x for x in fa if x not in fb]
                db = [x for x in fb if x not in fa]
                only = all(any(k in x for k in RESTATED) for x in da + db)
                ok = ok and only and len(fa) - len(da) == len(fb) - len(db)
                extra = "; full log: %d line(s) out, %d in, all in the restated clause: %s" % (len(da), len(db), only)
            L.append("  %-28s rc %5s -> %5s   verdict lines %s%s" % (n, ra.get(n, "-"), rb.get(n, "-"),
                     "identical (%d)" % len(vb) if va == vb else "DIFFER", extra))
        if not ok:
            bad.append(n)
    L.append("a2 %s (gate)" % ("PASS" if not bad else "FAIL: " + ", ".join(bad)))
    return "\n".join(L), not bad


# ---------------------------------------------------------------- a3
def A3():
    L, ok = ["a3 -- d = 0 was inert: the P2-16 binary vs the frozen runs (argv verbatim)", "=" * 100], True
    for leg, (ref, who) in REFS.items():
        b = os.path.join(LEGS, leg)
        ev = R14.same_file(ref, b, "run_removal_events.csv")
        n, nd, no = kept_plotfiles(ref, b)
        th, thm = R14.thermo_same(ref, b)
        cs, csm = check_leg(leg) if leg != "F0_k1" else (True, "(fixture deck: no shared files)")
        o = rc(b) == 0 and ev and n > 0 and nd == 0 and no == 0 and th and cs
        ok = ok and o
        L.append("  %s  %-8s vs %s: rc %s; event log %s; kept plotfiles %d files, %d differ, %d missing; thermo %s; %s"
                 % ("PASS" if o else "FAIL", leg, who, rc(b), "identical" if ev else "DIFFERS", n, nd, no, thm, csm))
    L.append("a3 %s (gate)" % ("PASS" if ok else "FAIL"))
    return "\n".join(L), ok


# ---------------------------------------------------------------- a4
def A4():
    p = os.path.join(HERE, "A4_retired_keys.log")
    t = open(p).read() if os.path.exists(p) else ""
    ok = bool(t) and t.rstrip().splitlines()[-1].startswith("PASS")
    L = ["a4 -- unit/retired_keys run half (one step of the gradient_blind fixture, -np 1)", "=" * 100]
    L += ["  " + x for x in t.splitlines() if x.startswith(("  ok", "  FAIL"))]
    L.append("a4 %s (gate)" % ("PASS" if ok else "FAIL"))
    return "\n".join(L), ok


# ---------------------------------------------------------------- c1
def C1():
    L = ["c1 -- SENSITIVITY, not evidence: the onset offset one span would need (R_30 = %s, R_BETHEL.txt)"
         % open(os.path.join(HERE, "R_BETHEL.txt")).read().strip(), "=" * 100]
    s0 = R15.stat(R15.events(os.path.join(LEGS, "CAL")))
    s1 = R15.stat(R15.events(os.path.join(LEGS, "CAL_m30")))
    d = s1["med"] - s0["med"]
    L.append("   CAL     %s" % R15.fmt(s0))
    L.append("   CAL_m30 %s   %s" % (R15.fmt(s1), check_leg("CAL_m30")[1]))
    L.append("   dS_top = %+.2f K (P2 predicted -30 +- 3 K: %s)" % (d, "inside" if abs(d + 30.0) <= 3.0 else "OUTSIDE"))
    g0, b0, _ = RU.windows()
    rows = []
    for leg in ("KB150_m30", "KB160_m30"):
        v, r = RU.l_req(LEGS, leg)
        rows.append(v)
        L.append("   %-9s earliest firing %s   %s" % (leg, ("t %.3f s  T_surf %.2f K  ell %.3f mm  mu %s  L_req %.3f mm"
                 % (float(r["time"]), float(r["T_surf"]), 1e3 * float(r["ell_eff"]), r["mu"], 1e3 * v)) if r
                 else "NONE by stop_time", check_leg(leg)[1]))
    if None in rows:
        L.append("   Bethel window at Delta -30 K: not scoreable (a leg did not fire)")
        return "\n".join(L), None
    g, b, Lw = RU.windows(LEGS, ("KB150_m30", "KB160_m30"))
    mv = [(b[i] - b0[i]) / b0[i] for i in (0, 1)]
    L.append("   Bethel window: Delta 0 [%.3f, %.3f) mm -> Delta -30 K [%.3f, %.3f) mm; edges move %+.3f %% / %+.3f %% "
             "(P2 predicted < 1 %%: %s)" % (1e3 * b0[0], 1e3 * b0[1], 1e3 * b[0], 1e3 * b[1], 100 * mv[0], 100 * mv[1],
                                           "held" if max(abs(x) for x in mv) < 0.01 else
                                           ("MISSED, > 2 %: the cancellation fails in the binary" if max(abs(x) for x in mv) > 0.02
                                            else "MISSED (1-2 %)")))
    # why it moved: L_req = C_L(T_s) c_f ell mu, so split the ratio (Delta -30 / Delta 0) into its two factors.
    # The offline cancellation (P2) assumed ell ~ dT_s at fixed flux; a transient skin grows ~ sqrt(t) instead.
    for l0, l1 in (("KB150_0", "KB150_m30"), ("KB160_0", "KB160_m30")):
        v0, r0 = RU.l_req(os.path.join(P15, "legs"), l0)
        v1, r1 = RU.l_req(LEGS, l1)
        T0s, T1s = float(r0["T_surf"]), float(r1["T_surf"])
        L.append("   %s/%s: dT_s x %.4f, C_L x %.4f, ell x %.4f, sqrt(t) x %.4f -> L_req x %.4f" % (
            l1, l0, (T1s - 293.15) / (T0s - 293.15), R15.DS.C_L(T1s) / R15.DS.C_L(T0s),
            float(r1["ell_eff"]) / float(r0["ell_eff"]), (float(r1["time"]) / float(r0["time"])) ** 0.5, v1 / v0))
    # the onset offset one span (R0) would need for a >= 2 % window, linear in Delta between the two Bethel states
    gh = g0[1]
    blo_need = gh * (1 - RU.MIN_REL / 2) / (1 + RU.MIN_REL / 2)     # (gh - b)/sqrt(gh b) = 2 % (to 1e-4)
    slope = (b[0] - b0[0]) / -30.0
    L.append("   one span (R0) needs Bethel's lower edge <= %.3f mm for a 2 %% window: an onset offset of %.1f K "
             "(linear between Delta 0 and -30; Kant 2016 Fig. 11/13 reads ~30 K [fig])" % (
                 1e3 * blo_need, (blo_need - b0[0]) / slope))
    txt, sc = RU.report(g, b, Lw, "   at Delta -30 K")
    L += ["   " + x for x in txt.splitlines()[2:] if x.strip().startswith(("Gotthard", "feasible", "R0", "R1", "R2"))]
    return "\n".join(L), dict(d=d, mv=mv, sc=sc)


def main(a):
    mode = a[0] if a else "RESULTS"
    for m, fn, f in (("A0", A0, "A0.txt"), ("A2", A2, "A2.txt"), ("A3", A3, "A3.txt"), ("A4", A4, "A4.txt")):
        if mode == m:
            t, ok = fn(); print(t); open(os.path.join(HERE, f), "w").write(t + "\n"); return 0 if ok else 1
    if mode == "C1":
        t, _ = C1(); print(t); open(os.path.join(HERE, "C1.txt"), "w").write(t + "\n"); return 0
    out = []
    for f in ("A0.txt", "RULES.txt", "GATE_tree.txt", "A2.txt", "A3.txt", "A4.txt"):
        p = os.path.join(HERE, f)
        out += [open(p).read().rstrip() if os.path.exists(p) else f + ": NOT ON DISK", ""]
    c1, cd = C1()
    out += [c1, ""]
    gates = {}
    for f, k in (("A0.txt", "a0"), ("GATE_tree.txt", "a1"), ("A2.txt", "a2"), ("A3.txt", "a3"), ("A4.txt", "a4")):
        p = os.path.join(HERE, f)
        t = open(p).read() if os.path.exists(p) else ""
        gates[k] = ("0 of 2639 changed" in t and "vs d2aa: 26 of 2639" in t and "vs d2t: 863 of 2639" in t) \
            if k == "a1" else bool(re.search(r"^%s PASS" % k, t, re.M))
    g, b, Lw = RU.windows()
    _, sc = RU.report(g, b, Lw, "")
    surv = [s for s in sc if s.get("surv")]
    out.append("e1 VERDICT")
    out.append("   d retired: %s (%s)" % ("YES" if all(gates.values()) else "NO",
                                          ", ".join("%s %s" % (k, "PASS" if v else "FAIL") for k, v in gates.items())))
    out.append("   rules at Delta 0 (c0): " + "; ".join("%s %s" % (s["tag"], "SURVIVES" if s.get("surv") else
                                                                  ("excluded" if not s["ok"] else "fails")) for s in sc))
    if cd:
        out.append("   why (c1): dS_top %+.2f K for R x %.3f; Bethel's window edges move %+.3f %% / %+.3f %% for a 30 K "
                   "onset offset" % (cd["d"], float(open(os.path.join(HERE, "R_BETHEL.txt")).read().split()[0])
                                     / 0.531294, 100 * cd["mv"][0], 100 * cd["mv"][1]))
    out.append("   Grimsel's span, two routes (for the step that switches B on; the span stays 0 here):")
    out.append("     (a) Gotthard carried over (conductivity-only transfer, P2-15): [%.3f, %.3f) mm" % (1e3 * g[0], 1e3 * g[1]))
    for s in surv:
        gr = RU.grimsel(s["gm"], s["n"])
        e_lo = min(RU.grimsel(s["lo"], s["n"])[1:]); e_hi = max(RU.grimsel(s["hi"], s["n"])[1:])
        out.append("     (b) %s %s: %.2f mm (UCS 169 +- 37 MPa at L_G*: %.2f-%.2f mm; with the window's edges: %.2f-%.2f mm)"
                   % (s["tag"], s["name"], 1e3 * gr[0], 1e3 * min(gr[1:]), 1e3 * max(gr[1:]), 1e3 * e_lo, 1e3 * e_hi))
    if not surv:
        out.append("     (b) no rule survives")
    print("\n".join(out))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
