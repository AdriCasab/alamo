#!/usr/bin/env python3
"""P2-11b -- the legs (ACTIVE_STEP.md `## Legs`). Binary: the Stage B build. Each leg writes argv.sh +
DROPPED.txt (binary sha256, what differs from its source argv) and a .done sentinel (SR 8).

  G1_C      studies/p2_11_heat/legs/G1/argv.sh token for token, plot_file re-pointed           120 s
  P47n_C    studies/p2_9_face_removal/legs/P47n/argv.sh minus P2-10's 13 retired tokens        120 s
  X1_S14, X1_S20   G1's argv with h_expr -> the LOCAL Goldstein form at C = S x 7272.3958 (no lid) and two
            tokens added: surface_patch.face_heat=path, surface_patch.barrier_expr=(r>0.02)and(z>0.15)   30 s
  X1        the same at S*_X1 (SSTAR.txt)                                                         120 s
  X1_S10 / X1_S28  the S* rule's ONE added point, only if SSTAR.txt asks for it

  run_p2_11b.py LEG ...            run (serial, -np 4; a leg with .done is skipped)
  run_p2_11b.py --smoke LEG ...    the leg's argv with stop_time 0.1 into smoke/<leg> (no .done)
  run_p2_11b.py --argv LEG         print the argv
"""
import argparse
import hashlib
import os
import shlex
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)
ROOT = os.path.abspath(os.path.join(STUD, "..", "..", ".."))
G1 = os.path.join(STUD, "p2_11_heat", "legs", "G1", "argv.sh")
P47N = os.path.join(STUD, "p2_9_face_removal", "legs", "P47n", "argv.sh")
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
SSTAR = os.path.join(HERE, "SSTAR.txt")
H0G = 7272.3958
G1_EXPR = ("if((r>0.02)and(s<0.0633),0.001,1189.763953*533*(533+0.3575*44*pow(r/0.0079,1.285))"
           "/pow(533+44*pow(r/0.0079,1.285),2)*pow(0.0632/max(0.0632,s),1.3185305455))")
X1_EXPR = ("{C}*533*(533+0.3575*44*pow(r/0.0079,1.285))/pow(533+44*pow(r/0.0079,1.285),2)"
           "*pow(0.0632/max(0.0632,s),1.3185305455)")
BARRIER = "(r>0.02)and(z>0.15)"
RETIRED = [   # P2-10's M0T_DELETE (studies/p2_10_retire/run_p2_10.py)
    "surface_patch.side_face_flux", "surface_patch.side_face_step_height",
    "surface_patch.side_face_area_law", "surface_patch.side_face_h", "surface_patch.side_face_T",
    "surface_patch.side_face_staircase", "spall.lateral_drive", "spall.lateral_removal",
    "spall.lateral_subcell", "spall.lateral_thermal", "surface_patch.side_face_h_datum",
    "surface_patch.side_face_h_r_min", "surface_patch.side_face_h_law"]


def s_star():
    if os.path.exists(SSTAR):
        for ln in open(SSTAR):
            if ln.startswith("S_STAR="):
                return float(ln.split("=", 1)[1])
    return None


def legs():
    L = {"G1_C": ("G1", None, "120.0"), "P47n_C": ("P47n", None, "120.0")}
    for nm, S in (("X1_S14", 0.14), ("X1_S20", 0.20), ("X1_S10", 0.10), ("X1_S28", 0.28)):
        L[nm] = ("G1", S, "30.0")
    st = s_star()
    if st is not None:
        L["X1"] = ("G1", st, "120.0")
    return L


def sha(path, n=24):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()[:n]


def argv_of(name, outdir, stop_override=None):
    src, S, stop = legs()[name]
    toks = shlex.split(open(G1 if src == "G1" else P47N).read())
    want, add = {}, []
    if S is not None:
        want["surface_patch.h_expr"] = X1_EXPR.format(C="%.6f" % (S * H0G))
        want["stop_time"] = stop
        add = [("surface_patch.face_heat", "path"), ("surface_patch.barrier_expr", BARRIER)]
    if stop_override:
        want["stop_time"] = stop_override
    out, rep, gone = [], [], []
    for t in toks:
        if "=" not in t or t.startswith("-"):
            out.append(t)
            continue
        k = t.split("=", 1)[0]
        if src == "P47n" and k in RETIRED:
            gone.append(t)
            continue
        if k == "plot_file":
            out.append("plot_file=" + os.path.join(outdir, "run"))
        elif k in want:
            if k == "surface_patch.h_expr" and src == "G1":
                assert t == "surface_patch.h_expr=" + G1_EXPR, t
            out.append(f"{k}={want[k]}")
            if out[-1] != t:
                rep.append(k)
        else:
            out.append(t)
    if src == "P47n":
        assert sorted(t.split("=", 1)[0] for t in gone) == sorted(RETIRED), gone
    out += [f"{k}={v}" for k, v in add]
    assert BIN in out
    return out, rep, [k for k, _ in add], gone, src


def write_prov(d, name, cmd, rep, add, gone, src):
    with open(os.path.join(d, "argv.sh"), "w") as f:
        f.write(" ".join(shlex.quote(c) for c in cmd) + "\n")
    srcp = G1 if src == "G1" else P47N
    with open(os.path.join(d, "DROPPED.txt"), "w") as f:
        f.write(f"P2-11b leg {name}: from {os.path.relpath(srcp, ROOT)} (sha256 {sha(srcp)}), token for token\n"
                f"binary bin/mmwspalling-3d-g++ sha256 {sha(BIN, 64)}\n\n-- keys REPLACED in place\n")
        for k in rep or ["(none)"]:
            f.write(f"   {k}\n")
        f.write("\n-- keys ADDED\n")
        for k in add or ["(none)"]:
            f.write(f"   {k}\n")
        f.write("\n-- tokens DROPPED" + (" (retired in P2-10):\n" if gone else ": none\n"))
        for t in gone:
            f.write(f"   {t}\n")
        f.write("\nplot_file re-pointed.\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("legs", nargs="+")
    ap.add_argument("--argv", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    a = ap.parse_args()
    L = legs()
    for n in a.legs:
        if n not in L:
            sys.exit(f"unknown or not-yet-available leg {n}; have {sorted(L)}")
    if not a.argv:
        subprocess.Popen(["caffeinate", "-dis", "-w", str(os.getpid())])   # SR 11
    for n in a.legs:
        d = os.path.join(HERE, "smoke" if a.smoke else "legs", n)
        cmd, rep, add, gone, src = argv_of(n, d, "0.1" if a.smoke else None)
        if a.argv:
            print(" \\\n  ".join(shlex.quote(c) for c in cmd))
            print(f"# replaced {rep}; added {add}; dropped {len(gone)}")
            continue
        if not a.smoke and os.path.exists(os.path.join(d, ".done")):
            print(f"  {n:7s} already done (.done) -- skipped", flush=True)
            continue
        os.makedirs(d, exist_ok=True)
        write_prov(d, n, cmd, rep, add, gone, src)
        p = os.path.join(d, ".done")
        if os.path.exists(p):
            os.remove(p)
        t0 = time.time()
        with open(os.path.join(d, "run.log"), "w") as f:
            rc = subprocess.call(cmd, stdout=f, stderr=subprocess.STDOUT, cwd=ROOT)
        wall = time.time() - t0
        if rc == 0 and not a.smoke:
            open(p, "w").write(f"{{'rc': 0, 'wall': {wall:.1f}}}\n")
        print(f"  {n:7s} rc={rc} wall {wall:7.1f} s  replaced {rep} added {add} dropped {len(gone)}", flush=True)
        if rc != 0:
            sys.exit(f"{n} failed (rc {rc}); see {d}/run.log")


if __name__ == "__main__":
    main()
