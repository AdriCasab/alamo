#!/usr/bin/env python3
"""P2-8 -- the face-list legs on Hoeser (47 MPa growth, Fig. 5.10).

Every leg is the base leg's argv (`studies/whatif_uptake/legs/s100/argv.sh`,
P2-7 HG at S = 0.10, token for token) with ONLY the keys of the packet's leg
table replaced in place, plus `plot_file`. Each leg writes argv.sh + DROPPED.txt
(standing rule 8). h_expr on EVERY leg is the base string with the lid edited to
`s<0.0633` (a1) and its leading constant replaced by C = S x 7272.3958.

LEGS (ACTIVE_STEP.md `## Legs`)
  K0      C 727.239580                         30 s, plot every 3750 steps
  L0      K0 + face_list = 1
  L0c     L0 + amr.n_cell = 25 25 75           plot 1875
          ⚠ + amr.blocking_factor = 1, DECLARED: 25 is not a multiple of the
            deck's blocking factor 2 and AMReX aborts at init (P2-5 hit the same
            defect and measured this fix decomposition-neutral); at max_level 0
            it constrains only the box decomposition.
  L0f     L0 + amr.n_cell = 100 100 300, timestep 0.004, max_grid_size 50, plot 7500
  C15/C20/C25  S 0.15 / 0.20 / 0.25 with the key
  P47     S* (SSTAR.txt, the d1 rule) with the key, 120 s, plot 1875
  P47off  S*, key OFF, 120 s, plot 1875

MODES
  run_p2_8.py K0 L0 ...       the legs (serial, -np 4, sentinel .done)
  run_p2_8.py --probe ALL     ONE step each into probe/<leg> (a5). Before SSTAR.txt
                              exists, P47/P47off are probed with the S = 0.20
                              constant standing in for C* (same string shape).
  run_p2_8.py --argv <leg>    print the argv and exit
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
BASE = os.path.join(STUD, "whatif_uptake", "legs", "s100", "argv.sh")
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
SSTAR = os.path.join(HERE, "SSTAR.txt")
H0G = 7272.3958
BASE_C = "727.239580"
BASE_LID = "s<0.0637"
LID = "s<0.0633"
ORDER = ["K0", "L0", "L0c", "L0f", "C15", "C20", "C25", "P47", "P47off"]


def hexpr(C):
    return ("if((r>0.02)and(%s),0.001,%s*533/(533+44*pow(r/0.0079,1.285))"
            "*pow(0.0632/max(0.0632,s),1.3185305455))" % (LID, C))


def c_of(S):
    return "%.6f" % (S * H0G)


def s_star(probe=False):
    if os.path.exists(SSTAR):
        for ln in open(SSTAR):
            if ln.startswith("S_STAR="):
                return float(ln.split("=", 1)[1])
    if probe:
        return 0.20
    sys.exit(f"{SSTAR} missing -- run reduce_p2_8.py SSTAR first (d1)")


def legs(probe=False):
    FL = {"surface_patch.face_list": "1"}
    L = {
        "K0":  {"surface_patch.h_expr": hexpr(BASE_C)},
        "L0":  dict(FL, **{"surface_patch.h_expr": hexpr(BASE_C)}),
        "L0c": dict(FL, **{"surface_patch.h_expr": hexpr(BASE_C),
                           "amr.n_cell": "25 25 75", "amr.plot_int": "1875",
                           "amr.blocking_factor": "1"}),
        "L0f": dict(FL, **{"surface_patch.h_expr": hexpr(BASE_C),
                           "amr.n_cell": "100 100 300", "timestep": "0.004",
                           "amr.max_grid_size": "50", "amr.plot_int": "7500"}),
        "C15": dict(FL, **{"surface_patch.h_expr": hexpr("1090.859370")}),
        "C20": dict(FL, **{"surface_patch.h_expr": hexpr("1454.479160")}),
        "C25": dict(FL, **{"surface_patch.h_expr": hexpr("1818.098950")}),
    }
    need_star = probe or os.path.exists(SSTAR)
    if need_star:
        Cs = c_of(s_star(probe))
        L["P47"] = dict(FL, **{"surface_patch.h_expr": hexpr(Cs),
                               "stop_time": "120.0", "amr.plot_int": "1875"})
        L["P47off"] = {"surface_patch.h_expr": hexpr(Cs),
                       "stop_time": "120.0", "amr.plot_int": "1875"}
    # the C values are the packet's literals; assert S x H0G reproduces them
    for nm, S in (("C15", 0.15), ("C20", 0.20), ("C25", 0.25)):
        assert c_of(S) in L[nm]["surface_patch.h_expr"], nm
    assert c_of(0.10) == BASE_C
    return L


def base_tokens():
    toks = shlex.split(open(BASE).read())
    for n, t in enumerate(toks):
        if t.endswith("input_lab"):
            return toks[:n + 1], toks[n + 1:]
    sys.exit("base argv: no input_lab deck path")


def argv_of(name, outdir, probe=False, extra=None):
    pre, kv = base_tokens()
    want = dict(legs(probe)[name])
    if extra:
        want.update(extra)
    keys = [t.split("=", 1)[0] for t in kv]
    assert len(keys) == len(set(keys)), "base argv has a duplicate key"
    # the K0 string is the base string with ONLY the lid edited
    base_h = next(t for t in kv if t.startswith("surface_patch.h_expr=")).split("=", 1)[1]
    assert hexpr(BASE_C) == base_h.replace(BASE_LID, LID), "K0's h_expr is not the base's with the lid edit"
    out, rep = [], []
    for t in kv:
        k = t.split("=", 1)[0]
        if k == "plot_file":
            out.append("plot_file=" + os.path.join(outdir, "run"))
        elif k in want:
            out.append(f"{k}={want[k]}")
            if t != out[-1]:
                rep.append(k)
        else:
            out.append(t)
    add = [k for k in want if k not in keys]
    out += [f"{k}={want[k]}" for k in add]
    return pre + out, rep, add


def sha(path, n=24):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()[:n]


def write_prov(d, name, cmd, rep, add):
    with open(os.path.join(d, "argv.sh"), "w") as f:
        f.write(" ".join(shlex.quote(c) for c in cmd) + "\n")
    with open(os.path.join(d, "DROPPED.txt"), "w") as f:
        f.write(f"P2-8 leg {name}: base = studies/whatif_uptake/legs/s100/argv.sh "
                f"(sha256 {sha(BASE)}), token for token\n")
        f.write(f"binary bin/mmwspalling-3d-g++ sha256 {sha(BIN)}\n\n")
        f.write("-- keys REPLACED in place\n")
        for k in rep or ["(none)"]:
            f.write(f"   {k}\n")
        f.write("\n-- keys ADDED\n")
        for k in add or ["(none)"]:
            f.write(f"   {k}\n")
        if name == "L0c":
            f.write("\n⚠ amr.blocking_factor=1 is NOT in the packet's leg table: DECLARED\n"
                    "  (25 is not a multiple of the deck's blocking factor 2; AMReX aborts).\n")
        f.write("\n-- keys DROPPED: none. plot_file re-pointed.\n")


def unused_block(log):
    L = open(log, errors="replace").read().splitlines()
    out, on = [], False
    for l in L:
        if l.startswith("Unused ParmParse Variables"):
            on = True
            continue
        if on:
            if not l.strip():
                break
            out.append(l.strip())
    return out


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
    return rc, wall, log


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("legs", nargs="+")
    ap.add_argument("--probe", action="store_true")
    ap.add_argument("--argv", action="store_true")
    a = ap.parse_args()
    if not a.argv:
        subprocess.Popen(["caffeinate", "-dis", "-w", str(os.getpid())])   # SR 11
    L = legs(a.probe)
    names = [n for n in ORDER if n in L] if a.legs == ["ALL"] else a.legs
    for n in names:
        if n not in L:
            sys.exit(f"unknown or not-yet-available leg {n}; have {sorted(L)}")
    for n in names:
        if a.probe:
            d = os.path.join(HERE, "probe", n)
            cmd, rep, add = argv_of(n, d, True, {"stop_time": "1e-9"})
        else:
            d = os.path.join(HERE, "legs", n)
            cmd, rep, add = argv_of(n, d)
        if a.argv:
            print(" \\\n  ".join(shlex.quote(c) for c in cmd))
            print(f"# replaced {rep}  added {add}")
            continue
        os.makedirs(d, exist_ok=True)
        write_prov(d, n, cmd, rep, add)
        rc, wall, log = run_one(d, cmd)
        print(f"  {n:6s} rc={rc} np={cmd[cmd.index('-np') + 1]} wall {wall:8.1f} s"
              f"  replaced {rep} added {add}", flush=True)
        if a.probe:
            print(f"         unused ParmParse: {unused_block(log)}", flush=True)


if __name__ == "__main__":
    main()
