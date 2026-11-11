#!/usr/bin/env python3
"""P2-29 stage 1 -- a7 + m2: THE FIVE FIXTURES, old deck on bin_pre vs the production-kernel deck on bin_1.

Runs (into S/legs/fix_old/<name>/run and S/legs/fix_new/<name>/run, argv + binary recorded beside each) the decks of
unit/beam, unit/equilibrium, unit/stefan, unit/convective_patch, then each test on both runs; unit/beam_tau's test runs
its own decks (--deck input --bin bin_pre: the old, into its output/, archived first by chain.sh; --deck input_phase:
the new, into output_phase/).
  a7 [gate]   every test PASSES on the NEW run, and its printed Temp range lies inside [200, 5000] K (the table).
  m2 [report] each test's scored number, old vs new; convective_patch: max |T_new - T_old| over every cell of every
              plotfile (prediction: the top cell within 1e-6 K -- the same problem on two kernels).
Usage: fixtures.py run | score
"""
import glob
import os
import re
import subprocess
import sys

import numpy as np
import yt

yt.set_log_level(50)
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
PY = "/Users/tzetze20/Desktop/code/.venv/bin/python"
U = os.path.join(ROOT, "tests", "MMWSpalling", "unit")
BIN = {"old": os.path.join(HERE, "bin_pre", "mmwspalling-3d-g++"), "new": os.path.join(HERE, "bin_1", "mmwspalling-3d-g++")}
DECK = {"old": "input", "new": "input_phase"}
NAMES = ["beam", "equilibrium", "stefan", "convective_patch"]
MPI = ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4"]
SCORE = re.compile(r"relative error|Temp range|Robustness|transient|profile|U1 |U2 |PASS|FAIL|lambda|front|"
                   r"T_sim|integral|Lambda_L", re.I)


def leg(side, name):
    return os.path.join(HERE, "legs", "fix_" + side, name)


def run():
    rc_all = 0
    for side in ("old", "new"):
        for n in NAMES:
            d = leg(side, n)
            os.makedirs(d, exist_ok=True)
            argv = MPI + [BIN[side], os.path.join("tests/MMWSpalling/unit", n, DECK[side]),
                          "plot_file=" + os.path.join(d, "run")]
            open(os.path.join(d, "argv.sh"), "w").write(" ".join(argv) + "\n")
            subprocess.run(["shasum", "-a", "256", BIN[side]], stdout=open(os.path.join(d, "BINARY.txt"), "w"))
            r = subprocess.run(argv, cwd=ROOT, capture_output=True, text=True)
            open(os.path.join(d, "run.log"), "w").write((r.stdout or "") + (r.stderr or ""))
            print("  run %-4s %-17s rc %d" % (side, n, r.returncode))
            rc_all |= r.returncode
    return 1 if rc_all else 0


def tests():
    out = {}
    for side in ("old", "new"):
        for n in NAMES:
            r = subprocess.run([PY, os.path.join(U, n, "test"), os.path.join(leg(side, n), "run")], cwd=ROOT,
                               capture_output=True, text=True)
            out[(side, n)] = (r.returncode, (r.stdout or "") + (r.stderr or ""))
        argv = [PY, os.path.join(U, "beam_tau", "test"), "--deck", DECK[side]] + (
            ["--bin", BIN["old"]] if side == "old" else ["--bin", BIN["new"]])
        r = subprocess.run(argv, cwd=ROOT, capture_output=True, text=True)
        out[(side, "beam_tau")] = (r.returncode, (r.stdout or "") + (r.stderr or ""))
    for (side, n), (rc, txt) in out.items():
        open(os.path.join(HERE, "K1_fix_%s_%s.log" % (n, side)), "w").write(txt)
    return out


def trange(txt):
    lo, hi = [], []
    for m in re.finditer(r"Temp range[^:]*:\s*([-\d.eE+]+)\s*\.\.\s*([-\d.eE+]+)", txt):
        lo.append(float(m.group(1)))
        hi.append(float(m.group(2)))
    for m in re.finditer(r"U1 Temp range: ([-\d.]+) \.\. ([-\d.]+) K \(t = 0\), ([-\d.]+) \.\. ([-\d.]+) K", txt):
        lo += [float(m.group(1)), float(m.group(3))]
        hi += [float(m.group(2)), float(m.group(4))]
    for m in re.finditer(r"T_top in \[([-\d.]+), ([-\d.]+)\]", txt):
        lo.append(float(m.group(1)))
        hi.append(float(m.group(2)))
    return (min(lo), max(hi)) if lo else None


def conv_diff():
    A = sorted(glob.glob(os.path.join(leg("old", "convective_patch"), "run", "*cell")))
    B = sorted(glob.glob(os.path.join(leg("new", "convective_patch"), "run", "*cell")))
    worst = top = 0.0
    for a, b in zip(A, B):
        ga, gb = yt.load(a), yt.load(b)
        ca = ga.covering_grid(0, ga.domain_left_edge, ga.domain_dimensions)["Temp"].to_value()
        cb = gb.covering_grid(0, gb.domain_left_edge, gb.domain_dimensions)["Temp"].to_value()
        worst = max(worst, float(np.abs(ca - cb).max()))
        top = max(top, float(np.abs(ca[:, :, -1] - cb[:, :, -1]).max()))
    return len(A), len(B), worst, top


def score():
    out = tests()
    bad = []
    for n in NAMES + ["beam_tau"]:
        print("=" * 110)
        for side in ("old", "new"):
            rc, txt = out[(side, n)]
            tr = trange(txt)
            print("  %-17s %-3s rc %d  T range %s" % (n, side, rc, ("%.3f .. %.3f K" % tr) if tr else "-"))
            for ln in txt.splitlines():
                if SCORE.search(ln) and "Loading" not in ln:
                    print("      " + ln.strip()[:170])
            if side == "new":
                if rc != 0:
                    bad.append("%s: the test fails on the new run (rc %d)" % (n, rc))
                if tr is None or tr[0] < 200.0 or tr[1] > 5000.0:
                    bad.append("%s: T range %s not inside [200, 5000] K" % (n, tr))
    na, nb, worst, top = conv_diff()
    print("=" * 110)
    print("m2 convective_patch, the same problem on two kernels: %d / %d plotfiles; max |T_new - T_old| every cell "
          "%.3e K, the top cell %.3e K (prediction: within 1e-6 K)" % (na, nb, worst, top))
    print("RESULT a7: %s" % ("FAIL: " + "; ".join(bad) if bad else "PASS (every fixture test passes on its "
                                                                   "production-kernel run, T inside the table)"))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit({"run": run, "score": score}[sys.argv[1]]())
