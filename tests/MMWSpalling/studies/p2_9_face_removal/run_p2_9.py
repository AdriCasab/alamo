#!/usr/bin/env python3
"""P2-9 rev 2 -- the line rule on Hoeser (47 MPa, Fig. 5.10).

Every leg is a P2-8 leg's argv (`studies/p2_8_face_list/legs/<from>/argv.sh`, binary
f4c27b27, token for token) with ONLY `plot_file` and the keys of the packet's leg
table replaced in place; binary = this step's (B). Each leg writes argv.sh +
DROPPED.txt (SR 8).

LEGS (ACTIVE_STEP.md `## Legs`)
  K0g    from K0  (key off)   --                                  30 s / 3750
  M0     from L0              --                                  30 s / 3750
  M0c    from L0c             -- (L0c already carries n_cell 25 25 75 and
                                  blocking_factor 1)              30 s / 1875
  M0f    from L0f             -- (n_cell 100 100 300, dt 0.004, mgs 50)  30 s / 7500
  C15n   from C15             -- (h_expr constant 1090.859370)    30 s / 3750
  P47r   from P47             -- (h_expr constant 955.592808, 120 s)  / 1875
  P47n   from P47             h_expr constant C*new = S*new x 7272.3958, ONLY if
                              |S*new - 0.1314| > 0.0039 (SSTAR.txt)
  C07n / C20n  the d1 rule's ONE added point (S 0.07 / 0.20), from C15 with the
              constant replaced -- only if SSTAR.txt asks for it

MODES
  run_p2_9.py M0 K0g ...      the legs (serial, -np 4, sentinel .done)
  run_p2_9.py --probe ALL     ONE step each into probe/<leg> (a7)
  run_p2_9.py --argv <leg>    print the argv and exit
"""
import argparse
import hashlib
import os
import re
import shlex
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)
ROOT = os.path.abspath(os.path.join(STUD, "..", "..", ".."))
P8 = os.path.join(STUD, "p2_8_face_list", "legs")
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
SSTAR = os.path.join(HERE, "SSTAR.txt")
H0G = 7272.3958
S_P28 = 0.1314
ORDER = ["K0g", "M0", "M0c", "M0f", "C15n", "C07n", "C20n", "P47r", "P47n"]
FROM = {"K0g": "K0", "M0": "L0", "M0c": "L0c", "M0f": "L0f", "C15n": "C15",
        "C07n": "C15", "C20n": "C15", "P47r": "P47", "P47n": "P47"}


def s_star():
    if not os.path.exists(SSTAR):
        return None
    for ln in open(SSTAR):
        if ln.startswith("S_STAR="):
            return float(ln.split("=", 1)[1])
    return None


def with_C(hexpr, C):
    """the P2-8 h_expr string with its leading constant (the `,C*533/` token)
    replaced -- the lid and every other token unchanged."""
    out, n = re.subn(r",(\d+\.\d+)\*533/", ",%s*533/" % C, hexpr, count=1)
    assert n == 1, hexpr
    return out


def legs():
    L = {"K0g": {}, "M0": {}, "M0c": {}, "M0f": {}, "C15n": {}, "P47r": {}}
    base_h = dict(t.split("=", 1) for t in shlex.split(open(os.path.join(P8, "C15", "argv.sh")).read())
                  if t.startswith("surface_patch.h_expr="))["surface_patch.h_expr"]
    L["C07n"] = {"surface_patch.h_expr": with_C(base_h, "%.6f" % (0.07 * H0G))}
    L["C20n"] = {"surface_patch.h_expr": with_C(base_h, "%.6f" % (0.20 * H0G))}
    st = s_star()
    if st is not None and abs(st - S_P28) > 0.0039:
        p47h = dict(t.split("=", 1) for t in shlex.split(open(os.path.join(P8, "P47", "argv.sh")).read())
                    if t.startswith("surface_patch.h_expr="))["surface_patch.h_expr"]
        L["P47n"] = {"surface_patch.h_expr": with_C(p47h, "%.6f" % (st * H0G))}
    return L


def argv_of(name, outdir, extra=None):
    toks = shlex.split(open(os.path.join(P8, FROM[name], "argv.sh")).read())
    want = dict(legs()[name])
    if extra:
        want.update(extra)
    out, rep = [], []
    keys = set()
    for t in toks:
        if "=" not in t or t.startswith("-"):
            out.append(t)
            continue
        k = t.split("=", 1)[0]
        keys.add(k)
        if k == "plot_file":
            out.append("plot_file=" + os.path.join(outdir, "run"))
        elif k in want:
            out.append(f"{k}={want[k]}")
            if out[-1] != t:
                rep.append(k)
        else:
            out.append(t)
    add = [k for k in want if k not in keys]
    out += [f"{k}={want[k]}" for k in add]
    assert BIN in out, "the P2-8 argv does not name bin/mmwspalling-3d-g++"
    return out, rep, add


def sha(path, n=24):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()[:n]


def write_prov(d, name, cmd, rep, add):
    src = os.path.join(P8, FROM[name], "argv.sh")
    with open(os.path.join(d, "argv.sh"), "w") as f:
        f.write(" ".join(shlex.quote(c) for c in cmd) + "\n")
    with open(os.path.join(d, "DROPPED.txt"), "w") as f:
        f.write(f"P2-9 leg {name}: from studies/p2_8_face_list/legs/{FROM[name]}/argv.sh "
                f"(sha256 {sha(src)}), token for token\n")
        f.write(f"binary bin/mmwspalling-3d-g++ sha256 {sha(BIN)}\n\n-- keys REPLACED in place\n")
        for k in rep or ["(none)"]:
            f.write(f"   {k}\n")
        f.write("\n-- keys ADDED\n")
        for k in add or ["(none)"]:
            f.write(f"   {k}\n")
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
    L = legs()
    names = [n for n in ORDER if n in L and n not in ("C07n", "C20n")] if a.legs == ["ALL"] else a.legs
    for n in names:
        if n not in L:
            sys.exit(f"unknown or not-yet-available leg {n}; have {sorted(L)}")
    for n in names:
        if a.probe:
            d = os.path.join(HERE, "probe", n)
            cmd, rep, add = argv_of(n, d, {"stop_time": "1e-9"})
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
            txt = open(log, errors="replace").read()
            notread = "NOT READ under the line rule" in txt
            print(f"         unused ParmParse: {unused_block(log)}; unread-keys line: {notread}",
                  flush=True)


if __name__ == "__main__":
    main()
