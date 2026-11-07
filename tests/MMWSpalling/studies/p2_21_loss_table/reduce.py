#!/usr/bin/env python3
"""P2-21 -- RESULTS.txt: the gate logs, then c5-c10 and e0 from the legs' per-step thermo columns (no plotfile is
read, so pruning cannot move a number) and c9 from T2's test log. Usage: reduce.py > RESULTS.txt

Every T-leg ran at Fig. 34's MEASURED incident power (run.py ramp_schedule, the user's decision 2026-11-06), not
the packet's constant 2 kW. The observed #19 record: melt starts 2-4 mm below the surface (DOE p. 53); onset
~300 +- 50 s after beam-on (Fig. 34: the collapse at ~505 s, beam-on ~210 s => ~295 s); the transmitted FRACTION
stays >= 0.79 of its beam-on value until the collapse (doe/FIG34.txt; 1.05 just before it).
"""
import importlib.machinery
import importlib.util
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
LEGS = os.path.join(HERE, "legs")
ZT = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "zhang", "zhang_oglesby", "test")
OBS_T, OBS_DT = 295.0, 50.0                # Fig. 34: collapse ~505 s - beam-on ~210 s; the packet's +- 50 s
OBS_Z = (2.0e-3, 4.0e-3)
OBS_TR = 0.75                              # the planner's envelope (digitised minimum 0.79)
BUCKLE = 2.4e-3


def load(name, path):
    keep = sys.argv
    sys.argv = [path]
    ld = importlib.machinery.SourceFileLoader(name, path)
    m = importlib.util.module_from_spec(importlib.util.spec_from_loader(name, ld))
    ld.exec_module(m)
    sys.argv = keep
    return m


Z = load("zt221", ZT)


def run_dir(leg):
    return os.path.join(LEGS, leg, "run")


def have(leg):
    return os.path.exists(os.path.join(run_dir(leg), "thermo.dat"))


def mm(x):
    return "%.2f mm" % (1e3 * x) if x is not None and math.isfinite(x) else "undefined"


def closure_stats(leg):
    h, rows, c = Z.table_rows(run_dir(leg))
    cl = c("beam_closure_err")
    o = Z.onset(run_dir(leg))
    at = o["beam_closure_err"] if o else cl[-1]
    cold = [x for x, t in zip(cl, c("bc_T_peak")) if t < 310.0]
    return dict(first=cl[0], cold=max(abs(x) for x in cold) if cold else float("nan"), at=at, lo=min(cl), hi=max(cl),
                colerr=max(c("beam_col_err_max")))


def legs_block():
    L = ["c5-c7, c10 -- at ONSET (the first step row with a finite bc_z_melt), every T-leg at the MEASURED power",
         "=" * 110,
         "   %-5s %-28s %8s %9s %9s %9s %9s %8s %8s %9s %9s" % (
             "leg", "setup", "onset s", "z_melt", "z_peak", "tr TOTAL", "tr CENT", "tau_max", "a_max", "ell_p",
             "ell_pk")]
    setup = {"T1a": "w0 13 mm, table", "T1b": "w0 16.5 mm, table", "T1c": "w0 20 mm, table (m 8.88)",
             "T0": "w0 20 mm, FLAT (m 1)", "M2": "w0 20 mm, m 2", "M3.5": "w0 20 mm, m 3.5",
             "T1d": "w0 16.5 mm, table, 100 mm", "T1x": "T1b explicit dt 0.03 s"}
    rec = {}
    for leg in ("T1a", "T1b", "T1c", "T0", "M2", "M3.5", "T1d", "T1x"):
        if not have(leg):
            L.append("   %-5s (not run)" % leg)
            continue
        r = Z.onset_record(run_dir(leg))
        rec[leg] = r
        if r is None:
            h, rows, c = Z.table_rows(run_dir(leg))
            L.append("   %-5s %-28s   no melt by %.0f s (peak %.1f K at %s)" % (
                leg, setup[leg], rows[-1][0], c("bc_T_peak")[-1], mm(c("bc_z_peak")[-1])))
            continue
        L.append("   %-5s %-28s %8.1f %9s %9s %9.3f %9.3f %8.3f %8.1f %9s %9s" % (
            leg, setup[leg], r["t"], mm(r["z_melt"]), mm(r["z_peak"]), r["tr_tot"], r["tr_c"], r["tau_max"],
            r["alpha_max"], mm(r["ell_p"]), mm(r["ell_pk"])))
    L.append("   observed #19: onset ~%.0f +- %.0f s; melt starts 2-4 mm deep; transmitted fraction >= %.2f of beam-on "
             "(digitised min 0.79, 1.05 just before the collapse). Control: p2_20 Z3b (flat, w0 20 mm, 2 kW) 318 s, "
             "9.0 mm" % (OBS_T, OBS_DT, OBS_TR))
    L.append("   c10: buckling compares the E-FOLD depth (ell_p face-referenced, ell_pk peak-referenced) with "
             "L/(C_L mu c_f) = %.1f mm (P2-19), not the peak depth" % (1e3 * BUCKLE))
    return "\n".join(L), rec


def b1r_block():
    L = ["b1r -- beam_closure_err per step (the lateral residual of a column sum under a spreading Gaussian; "
         "report) and beam_col_err_max (gate on T1b)", "=" * 110]
    for leg in ("T1b", "T1c", "T0", "T1a", "M2", "M3.5", "T1d", "T1x"):
        if have(leg):
            s = closure_stats(leg)
            L.append("   %-5s first %+.2e, cold rows within %.2e, at onset %+.2e, range %+.2e .. %+.2e;  max col_err "
                     "%.1e" % (leg, s["first"], s["cold"], s["at"], s["lo"], s["hi"], s["colerr"]))
    return "\n".join(L)


def c8_block(rec):
    L = ["c8 -- the scan at w0 20 mm read together: depth (c5) and TOTAL transmission ratio (c7) vs DOE 2-4 mm AND "
         ">= %.2f" % OBS_TR, "=" * 110]
    window = []
    for lab, leg in (("m 1", "T0"), ("m 2", "M2"), ("m 3.5", "M3.5"), ("m 8.88", "T1c")):
        r = rec.get(leg)
        if r is None:
            L.append("   %-7s %-5s no melt / not run" % (lab, leg))
            continue
        dz = OBS_Z[0] <= r["z_melt"] <= OBS_Z[1]
        tr = r["tr_tot"] >= OBS_TR
        window.append(dz and tr)
        L.append("   %-7s %-5s onset %6.1f s, z_melt %s (%s), total %.3f / central %.3f (%s)  => %s" % (
            lab, leg, r["t"], mm(r["z_melt"]), "in" if dz else "out", r["tr_tot"], r["tr_c"], "in" if tr else "out",
            "BOTH" if dz and tr else "not both"))
    L.append("   => the window that lands depth AND transmission is %s across the declared scan" % (
        "NON-EMPTY" if any(window) else "EMPTY"))
    return "\n".join(L), any(window)


def c9_block():
    p = os.path.join(HERE, "TEST_T2.log")
    if not os.path.exists(p):
        return "c9 -- T2 not run"
    L = ["c9 -- T2 (#28c at Oglesby's spot w0 14.1 mm, with the table; record: a glassed disc)", "=" * 110]
    for ln in open(p).read().splitlines():
        t = ln.split()
        if len(t) >= 3 and t[0].isdigit() and int(t[0]) in (50, 100, 150):
            L.append("   t %3s s: raw %s K, brightness %s K vs #28c uncorrected %s K" % (t[0], t[1], t[2],
                                                                                    t[3] if len(t) > 3 else "-"))
    return "\n".join(L)


def gate(path, pat):
    p = os.path.join(HERE, path)
    return os.path.exists(p) and re.search(pat, open(p).read(), re.M) is not None


def main():
    out = []
    for f in ("Z1_beam_tau.log", "Z4_retired_keys.log", "GATE_tree.txt", "GATE_face_cmp.txt",
              "A2_beam_void_closure.txt", "A2_beam.log", "A3.txt", "TEST_T1b.log", "TEST_T1x.log", "TEST_scan.log"):
        p = os.path.join(HERE, f)
        txt = open(p).read().rstrip() if os.path.exists(p) else "NOT ON DISK"
        if f == "GATE_face_cmp.txt":
            txt = "\n".join(txt.splitlines()[-6:])
        out += ["== %s" % f, txt, ""]
    out += [b1r_block(), ""]
    lb, rec = legs_block()
    out += [lb, ""]
    c8, win = c8_block(rec)
    out += [c8, "", c9_block(), ""]
    g = {
        "z1-z3": gate("Z1_beam_tau.log", r"^PASS$"),
        "z4": gate("Z4_retired_keys.log", r"^PASS"),
        "a0": gate("GATE_tree.txt", r"vs p2_18: 0 of 2639 changed, 0 only in new, 0 only in ref"),
        "a1": gate("GATE_face_cmp.txt", r"VERDICT: IDENTICAL"),
        "a2": gate("A2_beam_void_closure.txt", r"^a2 beam_void_closure PASS") and gate("A2_beam.log", r"PASS"),
        "a3": gate("A3.txt", r"^a3 PASS"),
        "b1-b2": gate("TEST_T1b.log", r"^PASS$"),
        "b3": gate("TEST_T1x.log", r"PASS  b3"),
    }
    out.append("e0 VERDICT")
    out.append("   the deposition integrates a varying coefficient exactly and every flat deck is untouched: %s (%s)" % (
        "YES" if all(g.values()) else "NO", ", ".join("%s %s" % (k, "PASS" if v else "FAIL") for k, v in g.items())))
    t1 = [rec.get(l) for l in ("T1a", "T1b", "T1c")]
    t1 = [r for r in t1 if r]
    if t1:
        dep = all(OBS_Z[0] <= r["z_melt"] <= OBS_Z[1] for r in t1)
        tim = any(abs(r["t"] - OBS_T) <= OBS_DT for r in t1)
        tra = all(r["tr_tot"] >= OBS_TR for r in t1)
        out.append("   the table on #19 (T1a-c): depth in 2-4 mm on all three: %s (%s); onset within 300 +- 50 s at any "
                   "spot: %s (%s); total transmission >= %.2f on all three: %s (%s)" % (
                       dep, ", ".join(mm(r["z_melt"]) for r in t1), tim, ", ".join("%.0f s" % r["t"] for r in t1),
                       OBS_TR, tra, ", ".join("%.2f" % r["tr_tot"] for r in t1)))
        reading = ("B (all three land)" if dep and tim and tra else
                   "C (the depth misses)" if not dep else
                   "A (the depth lands; " + ("onset" if not tim else "") + (" and " if not tim and not tra else "")
                   + ("transmission" if not tra else "") + " miss)")
        out.append("   P8 mechanical reading on T1a-c: %s; c8's window %s => %s" % (
            reading, "non-empty" if win else "empty", "D also holds (no m lands depth AND transmission)" if not win
            else "D does not hold"))
        # read with: the 0.78 mm cell, and the 2 s step. At the onset row the top ~2 mm is a near-isothermal plateau
        # that crosses T_m within ONE step, so c5 (the shallowest cell >= T_m) catches the plateau's upper edge; the
        # melt INITIATES at the peak (bc_z_peak at onset), which T1x (dt 0.03 s) confirms on T1b's case
        ini = [r["z_peak"] for r in t1]
        dep_i = all(OBS_Z[0] <= z <= OBS_Z[1] for z in ini)
        out.append("   read with: cell 0.78 mm (centres 0.39, 1.17, 1.95, 2.73 mm ...), step 2 s. The INITIATION depth (the "
                   "peak at onset): %s -- in 2-4 mm on all three: %s; T1x (dt 0.03 s) on T1b's case: first melt %s = its "
                   "peak %s, against T1b's c5 %s at dt 2 s" % (
                       ", ".join(mm(z) for z in ini), dep_i, mm(rec["T1x"]["z_melt"]) if rec.get("T1x") else "-",
                       mm(rec["T1x"]["z_peak"]) if rec.get("T1x") else "-", mm(rec["T1b"]["z_melt"]) if rec.get("T1b")
                       else "-"))
        out.append("   ⚠ the TOTAL transmission carries the closure residual (b1r: +%.1f / +%.1f / +%.1f %% of P_inc at "
                   "onset on T1a / T1b / T1c), of the same size as its distance from the 0.75 envelope" % tuple(
                       100 * closure_stats(l)["at"] for l in ("T1a", "T1b", "T1c")))
    print("\n".join(out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
