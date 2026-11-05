#!/usr/bin/env python3
"""P2-17 -- every row. Usage: reduce.py [A1|A2|A3|A4|A5] (no mode: RESULTS.txt to stdout).

Imported, never copied: the volume instrument fr_volume / pfs_by_t (studies/p2_11_heat/reduce_p2_11.py, over
P2-8's plotfile reader: geometry from each leg's OWN argv, no LZ constant), the S* rule (log V vs log S between the
bracketing legs, P2-11's sstar), the byte comparators (P2-14a same_file / thermo_same via P2-15's reduce), T_s*
(offline.py).

⚠ SR 7(i): every bk_* and fr_cap thermo column is PER-STEP (reset each step: FaceRemoval.H:341 and the P2-15 bk
reset) and the Hoeser deck writes thermo every 6th step (amr.thermo.plot_int = 6), so a thermo SUM is a 1-in-6
sample. A blocked FRACTION is a ratio of two sums over the same rows (the sampling cancels); a COUNT is printed as
sampled. e0's "the cap never acted" therefore also reads the event log, which has one row per removal: rc_flag bit 16
is h_f_bounded (the flake cap, spall.flake_h_max), bit 8 h_f_capped (the 15 mm search depth) -- complete, not sampled.
"""
import importlib.machinery
import importlib.util
import math
import os
import re
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
T = os.path.join(ROOT, "tests", "MMWSpalling")
STUD = os.path.join(T, "studies")
LEGS = os.path.join(HERE, "legs")
P16L = os.path.join(STUD, "p2_16_span_rule", "legs")
CHECK = os.path.join(T, "validation", "shared", "check_shared.py")
PY = sys.executable
TARGET_V = 13.0                      # Hoeser 0 MPa, 26 cm3/min x 30 s (hoeser_2016_drilling_series.csv)
SPAN = {"N": 0.012, "A": 0.0132}


def load(name, path):
    keep = sys.argv
    sys.argv = [path]
    ld = importlib.machinery.SourceFileLoader(name, path)
    m = importlib.util.module_from_spec(importlib.util.spec_from_loader(name, ld))
    ld.exec_module(m)
    sys.argv = keep
    return m


RUN = load("run17r", os.path.join(HERE, "run.py"))
OFF = load("off17r", os.path.join(HERE, "offline.py"))
R15 = OFF.RU.R15
R14 = R15.R14
R11R = load("red11_p217", os.path.join(STUD, "p2_11_heat", "reduce_p2_11.py"))
MOD = R11R._load("r8_p217", HERE)                 # P2-8's reader, legdir = HERE/legs/<n>


def rc(leg):
    p = os.path.join(LEGS, leg, "rc")
    return int(open(p).read()) if os.path.exists(p) else None


def done(leg):
    return os.path.exists(os.path.join(LEGS, leg, ".done"))


def check_leg(leg, extra=()):
    d = os.path.join(LEGS, leg)
    deck = open(os.path.join(d, "DECK.txt")).read().strip()
    allow = open(os.path.join(d, "ALLOW.txt")).read().split()
    a = [PY, CHECK, deck, os.path.join(d, "argv.sh")] + (["--allow"] + allow if allow else []) + list(extra)
    p = subprocess.run(a, cwd=ROOT, capture_output=True, text=True)
    return p.returncode == 0, p.stdout


def thermo(leg):
    return R15.thermo(os.path.join(LEGS, leg))


def tcol(leg, c):
    h, rows = thermo(leg)
    return [r[h.index(c)] for r in rows] if c in h else None


def volume(leg, t=30):
    g, P = R11R.pfs_by_t(MOD, leg)
    return R11R.fr_volume(MOD, leg, P[t])


def E_src(leg):
    """ledger_E_src (a cumulative STATE column) at the last thermo row (t = 30 s), x 4 (quarter domain) [kJ]"""
    h, rows = thermo(leg)
    return 4.0 * rows[-1][h.index("ledger_E_src")] / 1e3, rows[-1][h.index("time")]


def events(leg):
    return R15.events(os.path.join(LEGS, leg))


def flags(rows):
    n8 = n16 = 0
    for r in rows:
        f = int(round(float(r["rc_flag"])))
        n8 += bool(f & 8)
        n16 += bool(f & 16)
    return n8, n16


# ---------------------------------------------------------------- a1
def A1():
    L, ok = ["a1 -- CAL1 (P2-16 CAL argv on the EDITED deck) vs P2-16 legs/CAL (= P2-13)", "=" * 100], True
    a, b = os.path.join(P16L, "CAL"), os.path.join(LEGS, "CAL1")
    ev = R14.same_file(a, b, "run_removal_events.csv")
    n, nd, no = R14.same_plotfiles(a, b)
    th, thm = R14.thermo_same(a, b)
    cs, cso = check_leg("CAL1")
    i_line = [l for l in cso.splitlines() if l.strip().startswith("(i)")]
    one = bool(i_line) and "(1 declared" in i_line[0]
    mv = None
    for ln in open(os.path.join(b, "run", "metadata")):
        if ln.startswith("spall.buckle_span = "):
            mv = ln.split(" = ", 1)[1].strip()
    o = rc("CAL1") == 0 and ev and n > 0 and nd == 0 and no == 0 and th and cs and one and mv is not None \
        and float(mv) == 0.0
    L.append("  rc %s; event log %s; plotfiles %d files, %d differ, %d one-sided; thermo %s" % (
        rc("CAL1"), "identical" if ev else "DIFFERS", n, nd, no, thm))
    L += ["  " + x for x in cso.strip().splitlines()]
    L.append("  metadata spall.buckle_span = %s" % mv)
    L.append("a1 %s (gate)" % ("PASS" if o else "FAIL"))
    return "\n".join(L), o


# ---------------------------------------------------------------- a2
def A2():
    log = os.path.join(HERE, "face_heat.log")
    t = open(log).read() if os.path.exists(log) else ""
    rcp = os.path.join(HERE, "face_heat.rc")
    r = int(open(rcp).read()) if os.path.exists(rcp) else None
    out = os.path.join(T, "unit", "face_heat", "output")
    nb = ne = 0.0
    per = []
    for leg in sorted(os.listdir(out)) if os.path.isdir(out) else []:
        tp = os.path.join(out, leg, "run", "thermo.dat")
        if not os.path.exists(tp):
            continue
        h, rows = R15.thermo(os.path.join(out, leg))
        if "bk_n_eval" in h:
            e = sum(x[h.index("bk_n_eval")] for x in rows)
            bl = sum(x[h.index("bk_n_block")] for x in rows)
            ne += e
            nb += bl
            per.append("%s %g/%g" % (leg, bl, e))
    ref = os.path.join(STUD, "p2_16_span_rule", "UNIT_p2_16", "face_heat.log")
    rl = [R15.mask(x) for x in open(ref).read().splitlines()] if os.path.exists(ref) else []
    nl = [R15.mask(x) for x in t.splitlines()]
    moved = [x for x in nl if x not in rl]
    o = r == 0 and t.rstrip().endswith("PASS")
    L = ["a2 -- unit/face_heat on the edited shared files (input_lab now includes the span 0.012)", "=" * 100,
         "  rc %s; verdict %s" % (r, (t.rstrip().splitlines() or ["(no log)"])[-1]),
         "  B on its legs: sum bk_n_block / bk_n_eval = %g / %g (every thermo row; %s)" % (nb, ne, ", ".join(per)),
         "  log lines differing from P2-16's UNIT_p2_16/face_heat.log (line numbers masked): %d%s" % (
             len(moved), "".join("\n      " + x.strip()[:160] for x in moved[:12]))]
    L.append("a2 %s (gate)" % ("PASS" if o else "FAIL"))
    return "\n".join(L), o


# ---------------------------------------------------------------- a3
def A3():
    a, b = os.path.join(P16L, "HOS_R"), os.path.join(LEGS, "O_S164")
    ra = open(os.path.join(a, "run_removal_events.csv")).read().splitlines()
    rb = open(os.path.join(b, "run_removal_events.csv")).read().splitlines()
    head_ok = ra[0] == rb[0]
    # ⚠ the event log stamps a step's START (SR 7(iii)) and t accumulates in floating point, so the step that
    # STARTS at t = 15 s prints 14.999999999999234 and passes "<= 15.0" although HOS_R (stop 15) never ran it.
    # The steps HOS_R ran are those starting before 15.0 - dt/2 (dt 0.008). First scoring with "<= 15.0":
    # A3.first.txt (4 rows of that one step over).
    rb15 = [rb[0]] + [x for x in rb[1:] if float(x.split(",", 1)[0]) < 15.0 - 0.004]
    ev = rb15 == ra
    extra = [x for x in rb[1:] if 15.0 - 0.004 <= float(x.split(",", 1)[0]) <= 15.0]
    pf = "01875"
    nd = no = n = 0
    for suf in ("cell", "node"):
        A, B = os.path.join(a, "run", pf + suf), os.path.join(b, "run", pf + suf)
        fa = {os.path.relpath(os.path.join(p, f), A): os.path.join(p, f) for p, _, fs in os.walk(A) for f in fs}
        fb = {os.path.relpath(os.path.join(p, f), B): os.path.join(p, f) for p, _, fs in os.walk(B) for f in fs}
        no += len(set(fa) ^ set(fb))
        for k in set(fa) & set(fb):
            n += 1
            nd += R15.md5(fa[k]) != R15.md5(fb[k])
    o = rc("O_S164") == 0 and head_ok and ev and n > 0 and nd == 0 and no == 0
    L = ["a3 -- O_S164 (O tokens over the edited files) vs P2-16 HOS_R to 15 s", "=" * 100,
         "  rc %s; event log: header %s; rows of the steps HOS_R ran (start < 15.0 - dt/2): %d vs HOS_R's %d -> %s"
         "; %d row(s) of the step starting at t = 15 s (stamped %s) excluded" % (
             rc("O_S164"), "same" if head_ok else "DIFFERS", len(rb15) - 1, len(ra) - 1,
             "byte-identical" if ev else "DIFFER", len(extra), extra[0].split(",", 1)[0] if extra else "-"),
         "  15 s plotfile (%s cell + node): %d files, %d differ, %d one-sided" % (pf, n, nd, no)]
    L.append("a3 %s (gate)" % ("PASS" if o else "FAIL"))
    return "\n".join(L), o


# ---------------------------------------------------------------- a4
def A4():
    v = volume("O_S140")
    o = rc("O_S140") == 0 and abs(v - 24.922) <= 0.001
    return ("a4 -- O_S140 V(f) at 30 s = %.4f cm3 vs P2-11 G1_S14 24.922 cm3 (+-0.001): %s\na4 %s (gate)" % (
        v, "agree" if o else "DIFFER", "PASS" if o else "FAIL")), o


# ---------------------------------------------------------------- a5
def A5():
    ok_cs, out = check_leg("N_S070")
    meta = None
    for ln in open(os.path.join(LEGS, "N_S070", "run", "metadata")):
        if ln.startswith("spall.buckle_span = "):
            meta = ln.split(" = ", 1)[1].strip()
    ev = tcol("N_S070", "bk_n_eval")
    se = sum(ev) if ev else 0.0
    ii = [l for l in out.splitlines() if l.strip().startswith("(ii)")]
    o = rc("N_S070") == 0 and ok_cs and meta is not None and float(meta) == 0.012 and se > 0
    L = ["a5 -- N_S070 reaches B", "=" * 100,
         "  rc %s; check_shared %s; %s" % (rc("N_S070"), "PASS" if ok_cs else "FAIL", ii[0].strip() if ii else "-"),
         "  metadata spall.buckle_span = %s; sum bk_n_eval over thermo rows = %g (sampled 1 in 6)" % (meta, se)]
    L.append("a5 %s (gate)" % ("PASS" if o else "FAIL"))
    return "\n".join(L), o


# ---------------------------------------------------------------- c, e
def sstar(form):
    pts = []
    for n, (f, S) in RUN.ladder().items():
        if f == form and done(n) and rc(n) == 0:
            pts.append((S, volume(n), E_src(n)[0], n))
    pts.sort()
    for (s1, v1, e1, n1), (s2, v2, e2, n2) in zip(pts[:-1], pts[1:]):
        if (v1 - TARGET_V) * (v2 - TARGET_V) <= 0:
            fr = (math.log(TARGET_V) - math.log(v1)) / (math.log(v2) - math.log(v1))
            s = math.exp(math.log(s1) + fr * (math.log(s2) - math.log(s1)))
            e = math.exp(math.log(e1) + fr * (math.log(e2) - math.log(e1)))
            nexp = (math.log(v2) - math.log(v1)) / (math.log(s2) - math.log(s1))
            return dict(pts=pts, s=s, e=e, n=nexp, br=(n1, n2))
    return dict(pts=pts, s=None, e=None, n=None, br=None)


def sectors(leg, t=30, nsec=8, rmax=0.025):
    g, P = R11R.pfs_by_t(MOD, leg)
    f = MOD.read_pf(MOD.legdir(leg), P[t], g, comps=("fr_f",))["fr_f"]
    dx = g["dx"]
    x = g["lo"][0] + (np.arange(g["nx"]) + 0.5) * dx
    y = g["lo"][1] + (np.arange(g["ny"]) + 0.5) * dx
    X, Y = np.meshgrid(x, y, indexing="ij")
    r = np.hypot(X, Y)
    ang = np.degrees(np.arctan2(Y, X))
    sec = np.minimum((ang // (90.0 / nsec)).astype(int), nsec - 1)
    rem = (1.0 - f).sum(axis=2)                     # removed cells per column
    out = []
    for s in range(nsec):
        m = (sec == s) & (r <= rmax)
        out.append(rem[m].sum() / (m.sum() * g["nz"]))
    return out


def C1():
    L = ["c1 -- blocked fraction sum bk_n_block / sum bk_n_eval (thermo rows, sampled 1 in 6: the ratio cancels it); "
         "n_eval sampled", "=" * 100]
    for n, (f, S) in RUN.ladder().items():
        if f == "O" or not done(n):
            continue
        ev, bl = tcol(n, "bk_n_eval"), tcol(n, "bk_n_block")
        cls = {c: sum(tcol(n, "bk_n_block_" + c)) for c in ("floor", "wall", "roof")}
        se, sb = sum(ev), sum(bl)
        L.append("   %-7s S %.4f  blocked %s of n_eval %g (sampled)%s   floor %g wall %g roof %g" % (
            n, S, "%.4f" % (sb / se) if se > 0 else "-", se, "  (n_eval < 100: counts only)" if se < 100 else "",
            cls["floor"], cls["wall"], cls["roof"]))
    # ⚠ bk_n_eval counts a fired line once PER STEP: a line whose onset keeps passing while B blocks it is
    # re-evaluated (and re-blocked) every step, so the fraction above is EVALUATION-weighted. The per-firing
    # measure is the removal count itself, N (or A) vs O at the same S (event-log rows = removals, complete).
    L.append("   per FIRING (event-log removals, complete): N/O and A/O at equal S")
    for S in RUN.S_GRID["O"]:
        o = events(RUN.name_of("O", S)) if done(RUN.name_of("O", S)) else None
        if not o:
            continue
        parts = []
        for f in ("N", "A"):
            n = RUN.name_of(f, S)
            if done(n):
                parts.append("%s/O %.3f (%d / %d)" % (f, len(events(n)) / len(o), len(events(n)), len(o)))
        L.append("     S %.4f  %s" % (S, "   ".join(parts)))
    return "\n".join(L)


def C2C3():
    L = ["c2/c3 -- S* (V(f) 30 s vs %.1f cm3, P2-11's rule), the gap 1/S*, ledger_E_src x4 at S*, the exponent n on "
         "the bracket" % TARGET_V, "=" * 100]
    st = {}
    for form in ("O", "N", "A"):
        r = sstar(form)
        st[form] = r
        for S, v, e, n in r["pts"]:
            L.append("   %-7s S %.4f  V(f) 30 s %8.3f cm3  E_src x4 %9.2f kJ   %s" % (n, S, v, e, check_leg(n)[0]
                                                                                      and "check_shared PASS"
                                                                                      or "check_shared FAIL"))
        if r["s"]:
            L.append("   S*(%s) = %.5f  (bracket %s-%s)  gap 1/S* = %.3f  E_src at S* %.2f kJ  n = %.3f" % (
                form, r["s"], r["br"][0], r["br"][1], 1.0 / r["s"], r["e"], r["n"]))
        else:
            L.append("   S*(%s): NOT BRACKETED by the grid" % form)
    if st["O"]["s"] and st["N"]["s"]:
        L.append("   S*(N)/S*(O) = %.4f;  gap %.2fx -> %.2fx;  E_src at S* N/O = %.4f;  n(N) - n(O) = %+.3f" % (
            st["N"]["s"] / st["O"]["s"], 1 / st["O"]["s"], 1 / st["N"]["s"], st["N"]["e"] / st["O"]["e"],
            st["N"]["n"] - st["O"]["n"]))
    if st["A"]["s"] and st["N"]["s"]:
        L.append("   S*(A)/S*(N) = %.4f" % (st["A"]["s"] / st["N"]["s"]))
    return "\n".join(L), st


def C4():
    L = ["c4 -- isotropy witness: removed fraction per 11.25 deg sector (r <= 25 mm, 30 s) / the sector's cell census; "
         "spread = max/min - 1", "=" * 100]
    sp = {}
    for S in (0.10, 0.07):
        for f in ("O", "N"):
            n = RUN.name_of(f, S)
            if not done(n):
                continue
            s = sectors(n)
            sp[(f, S)] = max(s) / min(s) - 1.0
            L.append("   %-7s spread %.4f   sectors %s" % (n, sp[(f, S)], " ".join("%.4f" % x for x in s)))
        if ("O", S) in sp and ("N", S) in sp:
            d = sp[("N", S)] - sp[("O", S)]
            L.append("   S %.3f: N - O = %+.4f  (> 0.05 would mean B adds a grid-direction bias: %s)" % (
                S, d, "YES" if d > 0.05 else "no"))
    return "\n".join(L)


def E0():
    ts_star = {f: OFF.t_star(SPAN[f]) for f in SPAN}
    L = ["e0 -- the cap is dead as a bound: sum fr_cap (thermo, sampled 1 in 6) and the event log's rc_flag bit 16 "
         "(h_f_bounded, every removal) on every N/A leg; positive control on O_S070 / O_S087; max firing T_s vs T_s*",
         "=" * 100]
    ok_na, ok_ctrl, ok_ts = True, True, True
    for n, (f, S) in RUN.ladder().items():
        if not done(n):
            continue
        fc = tcol(n, "fr_cap")
        ev = events(n)
        n8, n16 = flags(ev)
        mx = max(float(r["T_surf"]) for r in ev) if ev else float("nan")
        tag = ""
        if f in SPAN:
            ts = ts_star[f]
            dead = sum(fc) == 0 and n16 == 0
            below = ts is None or mx < ts
            ok_na = ok_na and dead
            ok_ts = ok_ts and below
            tag = "T_s* %s: max T_s %s" % ("%.1f K" % ts if ts else "none", "below" if below else "ABOVE")
        elif S in (0.07, 0.087):
            ctrl = sum(fc) > 0 or n16 > 0
            ok_ctrl = ok_ctrl and ctrl
            tag = "positive control: %s" % ("the cap acted" if ctrl else "NEVER ACTED")
        L.append("   %-7s sum fr_cap %4g (sampled)  events %5d, bit16 (cap) %4d, bit8 (15 mm) %d  max T_s %.1f K  %s"
                 % (n, sum(fc), len(ev), n16, n8, mx, tag))
    L.append("   cap count zero on every N/A leg: %s; positive control > 0 on O_S070 and O_S087: %s; every N/A max "
             "T_s below its T_s*: %s" % (ok_na, ok_ctrl, ok_ts))
    v = ok_na and ok_ctrl and ok_ts
    L.append("e0 %s (verdict)" % ("the cap is DEAD as a bound under B" if v else
                                  ("NOT shown dead: " + ", ".join(k for k, x in (("a cap count on N/A", not ok_na),
                                   ("zero positive control (deadness rests on z1 alone)", not ok_ctrl),
                                   ("a firing above T_s*", not ok_ts)) if x))))
    return "\n".join(L), dict(na=ok_na, ctrl=ok_ctrl, ts=ok_ts)


def main(a):
    mode = a[0] if a else "RESULTS"
    for m, fn, f in (("A1", A1, "A1.txt"), ("A2", A2, "A2.txt"), ("A3", A3, "A3.txt"), ("A4", A4, "A4.txt"),
                     ("A5", A5, "A5.txt")):
        if mode == m:
            t, ok = fn(); print(t); open(os.path.join(HERE, f), "w").write(t + "\n"); return 0 if ok else 1
    out = []
    for f in ("Z0.txt", "Z1.txt", "Z2.txt", "A1.txt", "A2.txt", "A3.txt", "A4.txt", "A5.txt"):
        p = os.path.join(HERE, f)
        out += [open(p).read().rstrip() if os.path.exists(p) else f + ": NOT ON DISK", ""]
    out += [C1(), ""]
    c2, st = C2C3()
    out += [c2, "", C4(), ""]
    e0, ev = E0()
    out += [e0, ""]
    gates = {k: bool(re.search(r"^%s PASS" % k, open(os.path.join(HERE, k.upper() + ".txt")).read(), re.M))
             if os.path.exists(os.path.join(HERE, k.upper() + ".txt")) else False
             for k in ("z0", "a1", "a2", "a3", "a4", "a5")}
    out.append("e1 VERDICT")
    out.append("   B on is the production form: %s (%s)" % ("YES" if all(gates.values()) else "NO", ", ".join(
        "%s %s" % (k, "PASS" if v else "FAIL") for k, v in gates.items())))
    out.append("   the cap: %s" % e0.splitlines()[-1].split(" ", 1)[1])
    if st["O"]["s"] and st["N"]["s"]:
        r = st["N"]["s"] / st["O"]["s"]
        br = "(c) > 1.5" if r > 1.5 else ("(a) 1.02-1.5" if r >= 1.02 else ("(b) < 1.02" if r >= 0.99 else
                                                                          "BELOW 0.99: B's withdrawals raised removal"))
        out.append("   P4: S*(N)/S*(O) = %.4f -> branch %s; the gap 1/S* %.2fx (O) -> %.2fx (N)" % (
            r, br, 1 / st["O"]["s"], 1 / st["N"]["s"]))
    if st["A"]["s"] and st["N"]["s"]:
        out.append("   P5: S*(A)/S*(N) = %.4f" % (st["A"]["s"] / st["N"]["s"]))
    z2 = open(os.path.join(HERE, "Z2.txt")).read() if os.path.exists(os.path.join(HERE, "Z2.txt")) else ""
    m12 = [l.strip() for l in z2.splitlines() if l.strip().startswith("L  12.00 mm")]
    out.append("   the declared span's minimum flux (z2): %s" % (m12[0] if m12 else "-"))
    print("\n".join(out))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
