#!/usr/bin/env python3
"""P2-11 -- the G1 probe: Hoeser 47 MPa with Goldstein's profile used as a LOCAL value.

Every leg is P2-9's P47n argv (`studies/p2_9_face_removal/legs/P47n/argv.sh`, binary B) token for token,
with `plot_file` re-pointed, `stop_time` set per leg, and ONLY `surface_patch.h_expr` replaced by the G1
formula at C = S x 7272.3958 (6 decimals).
⚠ DECLARED (not in the packet's text): the 13 tokens P2-10 retired (M0T_DELETE, the exact list P2-10's
M0T dropped from M0's argv) are DROPPED -- binary C2 aborts on each of them. P2-10 a1 (on substance) and
a6 showed the face-list path does not read them (M0 -> M0T -> C2 identical), so P47n stays the control.
Each leg writes argv.sh + DROPPED.txt + the binary hash (SR 8) and a .done sentinel.

LEGS (ACTIVE_STEP.md `## Legs`)
  G1_S14  S 0.14  30 s     G1_S20  S 0.20  30 s     G1  S*_G1 (SSTAR.txt)  120 s
  G1_S10 / G1_S28: the S* rule's ONE added point, only if SSTAR.txt asks for it

  run_p2_11.py G1_S14 G1_S20 ...   the legs (serial, -np 4)
  run_p2_11.py --argv <leg>        print the argv and exit
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
P47N = os.path.join(STUD, "p2_9_face_removal", "legs", "P47n", "argv.sh")
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
BIN_SHA = "01c49ab554cdeab0"
SSTAR = os.path.join(HERE, "SSTAR.txt")
H0G = 7272.3958
G1_EXPR = ("if((r>0.02)and(s<0.0633),0.001,{C}*533*(533+0.3575*44*pow(r/0.0079,1.285))"
           "/pow(533+44*pow(r/0.0079,1.285),2)*pow(0.0632/max(0.0632,s),1.3185305455))")
P47N_EXPR = ("if((r>0.02)and(s<0.0633),0.001,1040.679839*533/(533+44*pow(r/0.0079,1.285))"
             "*pow(0.0632/max(0.0632,s),1.3185305455))")
# P2-10's M0T_DELETE (studies/p2_10_retire/run_p2_10.py), by key
RETIRED = [
    "surface_patch.side_face_flux", "surface_patch.side_face_step_height",
    "surface_patch.side_face_area_law", "surface_patch.side_face_h", "surface_patch.side_face_T",
    "surface_patch.side_face_staircase", "spall.lateral_drive", "spall.lateral_removal",
    "spall.lateral_subcell", "spall.lateral_thermal", "surface_patch.side_face_h_datum",
    "surface_patch.side_face_h_r_min", "surface_patch.side_face_h_law"]
ORDER = ["G1_S14", "G1_S20", "G1_S10", "G1_S28", "G1"]


def s_star():
    if not os.path.exists(SSTAR):
        return None
    for ln in open(SSTAR):
        if ln.startswith("S_STAR="):
            return float(ln.split("=", 1)[1])
    return None


def legs():
    L = {"G1_S14": (0.14, "30.0"), "G1_S20": (0.20, "30.0"), "G1_S10": (0.10, "30.0"),
         "G1_S28": (0.28, "30.0")}
    st = s_star()
    if st is not None:
        L["G1"] = (st, "120.0")
    return L


def h_expr(S):
    return G1_EXPR.format(C="%.6f" % (S * H0G))


def argv_of(name, outdir):
    S, stop = legs()[name]
    toks = shlex.split(open(P47N).read())
    want = {"surface_patch.h_expr": h_expr(S), "stop_time": stop}
    out, rep, gone = [], [], []
    for t in toks:
        if "=" not in t or t.startswith("-"):
            out.append(t)
            continue
        k = t.split("=", 1)[0]
        if k in RETIRED:
            gone.append(t)
            continue
        if k == "plot_file":
            out.append("plot_file=" + os.path.join(outdir, "run"))
        elif k in want:
            if k == "surface_patch.h_expr":
                assert t == "surface_patch.h_expr=" + P47N_EXPR, t     # the control's h, as the packet quotes it
            out.append(f"{k}={want[k]}")
            if out[-1] != t:
                rep.append(k)
        else:
            out.append(t)
    assert sorted(t.split("=", 1)[0] for t in gone) == sorted(RETIRED), gone
    assert BIN in out
    return out, rep, gone


def sha(path, n=24):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()[:n]


def write_prov(d, name, cmd, rep, gone):
    with open(os.path.join(d, "argv.sh"), "w") as f:
        f.write(" ".join(shlex.quote(c) for c in cmd) + "\n")
    with open(os.path.join(d, "DROPPED.txt"), "w") as f:
        f.write(f"P2-11 leg {name}: from studies/p2_9_face_removal/legs/P47n/argv.sh (sha256 {sha(P47N)}), "
                f"token for token\nbinary bin/mmwspalling-3d-g++ sha256 {sha(BIN, 64)}\n\n-- keys REPLACED in place\n")
        for k in rep:
            f.write(f"   {k}\n")
        f.write("\n-- keys ADDED\n   (none)\n\n-- keys DROPPED (retired in P2-10; C2 aborts on each):\n")
        for t in gone:
            f.write(f"   {t}\n")
        f.write("\nplot_file re-pointed.\n")


def run_one(d, cmd):
    os.makedirs(d, exist_ok=True)
    p = os.path.join(d, ".done")
    if os.path.exists(p):
        os.remove(p)
    t0 = time.time()
    log = os.path.join(d, "run.log")
    with open(log, "w") as f:
        rc = subprocess.call(cmd, stdout=f, stderr=subprocess.STDOUT, cwd=ROOT)
    wall = time.time() - t0
    if rc == 0:
        open(p, "w").write(f"{{'rc': 0, 'wall': {wall:.1f}}}\n")
    return rc, wall


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("legs", nargs="+")
    ap.add_argument("--argv", action="store_true")
    a = ap.parse_args()
    L = legs()
    for n in a.legs:
        if n not in L:
            sys.exit(f"unknown or not-yet-available leg {n}; have {sorted(L)}")
    if not a.argv:
        h = sha(BIN, 64)
        if not h.startswith(BIN_SHA):
            sys.exit(f"binary {h[:16]} is not P2-10's C2 ({BIN_SHA})")
        subprocess.Popen(["caffeinate", "-dis", "-w", str(os.getpid())])   # SR 11
    for n in a.legs:
        d = os.path.join(HERE, "legs", n)
        cmd, rep, gone = argv_of(n, d)
        if a.argv:
            print(" \\\n  ".join(shlex.quote(c) for c in cmd))
            print(f"# replaced {rep}; dropped {len(gone)} retired tokens")
            continue
        if os.path.exists(os.path.join(d, ".done")):
            print(f"  {n:7s} already done (.done) -- skipped", flush=True)
            continue
        os.makedirs(d, exist_ok=True)
        write_prov(d, n, cmd, rep, gone)
        rc, wall = run_one(d, cmd)
        print(f"  {n:7s} S {L[n][0]:.4f} stop {L[n][1]} rc={rc} wall {wall:7.1f} s  replaced {rep}", flush=True)
        if rc != 0:
            sys.exit(f"{n} failed (rc {rc}); see {d}/run.log")


if __name__ == "__main__":
    main()
