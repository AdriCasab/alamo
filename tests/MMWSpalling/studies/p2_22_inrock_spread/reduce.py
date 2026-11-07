#!/usr/bin/env python3
"""P2-22 -- RESULTS.txt. Usage: reduce.py > RESULTS.txt
Gates (unchanged fixtures, eps_r 1): z1-z4, a0 tree, a1 face, a2. Reports: C1 moved by how much; the legs A / X / F at
the TRUE spot (omega0 0.01288, Appendix Table 2 standoff 0.31") against P2-21's T1a (omega0 0.013, vacuum spread) and
the #19 record: first melt 2-4 mm below the surface, a pool at ~295 +- 50 s after beam-on, whole-block transmission
within the etalon envelope (0.79-1.17 of beam-on; a BOUND, not a curve: the DOE documents the disc Fabry-Perot).
"""
import importlib.machinery, importlib.util, math, os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__))
P21 = os.path.join(HERE, "..", "p2_21_loss_table")
def load(name, path):
    keep = sys.argv; sys.argv = [path]
    ld = importlib.machinery.SourceFileLoader(name, path)
    m = importlib.util.module_from_spec(importlib.util.spec_from_loader(name, ld)); ld.exec_module(m)
    sys.argv = keep; return m
R21 = load("r21", os.path.join(P21, "reduce.py"))
Z = R21.Z
mm = R21.mm
def rd(base, leg): return os.path.join(base, "legs", leg, "run")
def have(base, leg): return os.path.exists(os.path.join(rd(base, leg), "thermo.dat"))
def cstats(run):
    h, rows, c = Z.table_rows(run); cl = c("beam_closure_err"); o = Z.onset(run)
    at = o["beam_closure_err"] if o else cl[-1]
    cold = [x for x, t in zip(cl, c("bc_T_peak")) if t < 310.0]
    return dict(first=cl[0], cold=max(abs(x) for x in cold) if cold else float("nan"), at=at, lo=min(cl), hi=max(cl),
                colerr=max(c("beam_col_err_max")))
def gate(path, pat):
    p = os.path.join(HERE, path); return os.path.exists(p) and re.search(pat, open(p).read(), re.M) is not None

out = []
for f in ("Z1_beam_tau.log", "Z4_retired_keys.log", "GATE_tree.txt", "GATE_face_cmp.txt", "A2_beam_void_closure.txt",
          "A2_beam.log", "C1_vs_Z1.txt", "TEST_A.log", "TEST_X.log"):
    p = os.path.join(HERE, f); txt = open(p).read().rstrip() if os.path.exists(p) else "NOT ON DISK"
    if f == "GATE_face_cmp.txt": txt = "\n".join(txt.splitlines()[-6:])
    out += ["== %s" % f, txt, ""]

legs = [("T1a (P2-21)", P21, "T1a", "w0 13 mm, table, VACUUM spread"),
        ("A", HERE, "A", "w0 12.88 mm, table, in-rock spread"),
        ("X", HERE, "X", "A explicit dt 0.03 s"),
        ("F", HERE, "F", "w0 12.88 mm, FLAT, in-rock spread"),
        ("T0 (P2-21)", P21, "T0", "w0 20 mm, FLAT, vacuum spread")]
out += ["legs -- at ONSET (the first step row with a finite bc_z_melt); closure = beam_closure_err (the lateral residual)",
        "=" * 118,
        "   %-12s %-34s %8s %8s %8s %8s %8s %7s %8s %8s %9s %9s" % ("leg", "setup", "onset s", "z_melt", "z_peak",
        "tr TOT", "tr CENT", "tau_mx", "clos@on", "clos_max", "ell_p", "ell_pk")]
rec = {}
for lab, base, leg, setup in legs:
    if not have(base, leg):
        out.append("   %-12s (not run)" % lab); continue
    run = rd(base, leg); r = Z.onset_record(run); s = cstats(run); rec[lab] = (r, s)
    if r is None:
        h, rows, c = Z.table_rows(run)
        out.append("   %-12s %-34s no melt by %.0f s (peak %.1f K at %s); closure max %+.2e" % (
            lab, setup, rows[-1][0], c("bc_T_peak")[-1], mm(c("bc_z_peak")[-1]), s["hi"])); continue
    out.append("   %-12s %-34s %8.1f %8s %8s %8.3f %8.3f %7.3f %+8.1e %+8.1e %9s %9s" % (
        lab, setup, r["t"], mm(r["z_melt"]), mm(r["z_peak"]), r["tr_tot"], r["tr_c"], r["tau_max"], s["at"], s["hi"],
        mm(r["ell_p"]), mm(r["ell_pk"])))
out.append("   observed #19: first melt 2-4 mm below the surface (one sentence, DOE p. 53); pool at ~295 +- 50 s; "
           "whole-block transmission 0.79-1.17 of beam-on (etalon envelope, a bound)")
out.append("")
A = rec.get("A", (None, None)); T = rec.get("T1a (P2-21)", (None, None)); X = rec.get("X", (None, None)); F = rec.get("F", (None, None))
out.append("verdict")
g = {"z1-z3": gate("Z1_beam_tau.log", r"^PASS$"), "z4": gate("Z4_retired_keys.log", r"^PASS"),
     "a0": gate("GATE_tree.txt", r"vs p2_18: 0 of 2639 changed, 0 only in new, 0 only in ref"),
     "a1": gate("GATE_face_cmp.txt", r"VERDICT: IDENTICAL"),
     "a2": gate("A2_beam_void_closure.txt", r"^a2 beam_void_closure PASS") and gate("A2_beam.log", r"PASS")}
out.append("   every eps_r 1 fixture untouched: %s (%s)" % ("YES" if all(g.values()) else "NO",
           ", ".join("%s %s" % (k, "PASS" if v else "FAIL") for k, v in g.items())))
if A[0] and T[0]:
    out.append("   closure residual at onset, w0 13 mm: vacuum spread %+.2e -> in-rock %+.2e (x%.2f); run max %+.2e -> %+.2e "
               "(pre-registered: ~1/5)" % (T[1]["at"], A[1]["at"], A[1]["at"] / T[1]["at"], T[1]["hi"], A[1]["hi"]))
    out.append("   onset %.1f -> %.1f s; initiation depth (peak at onset) %s -> %s; whole-block transmission ratio %.3f -> %.3f; "
               "central %.3f -> %.3f" % (T[0]["t"], A[0]["t"], mm(T[0]["z_peak"]), mm(A[0]["z_peak"]), T[0]["tr_tot"],
                                         A[0]["tr_tot"], T[0]["tr_c"], A[0]["tr_c"]))
if A[0]:
    r = A[0]; dep = 2e-3 <= r["z_peak"] <= 4e-3; tim = abs(r["t"] - 295) <= 50; tra = 0.79 <= r["tr_tot"]
    out.append("   #19 at the true spot, table: initiation depth %s (%s 2-4 mm); first melt at %.0f s (%s 295 +- 50; a pool "
               "forms LATER than first melt, so the true gap is smaller than %.1fx); whole-block transmission %.3f (%s the "
               "0.79 envelope, +-%.1e closure)" % (mm(r["z_peak"]), "IN" if dep else "OUT of", r["t"], "IN" if tim else
               "OUT of", 295.0 / r["t"], r["tr_tot"], "inside" if tra else "below", abs(A[1]["at"])))
if X[0]:
    out.append("   X (dt 0.03 s): first melt %s = its peak %s at %.1f s (A: %.1f s, %+.1f %%)" % (
        mm(X[0]["z_melt"]), mm(X[0]["z_peak"]), X[0]["t"], A[0]["t"] if A[0] else float("nan"),
        100 * (X[0]["t"] / A[0]["t"] - 1) if A[0] else float("nan")))
if F[0]:
    out.append("   F (flat alpha 14 at the true spot): first melt %s (peak %s) at %.0f s, transmission %.3f, closure %+.1e" % (
        mm(F[0]["z_melt"]), mm(F[0]["z_peak"]), F[0]["t"], F[0]["tr_tot"], F[1]["at"]))
elif "F" in rec:
    out.append("   F (flat alpha 14 at the true spot): no melt within 400 s")
print("\n".join(out))
