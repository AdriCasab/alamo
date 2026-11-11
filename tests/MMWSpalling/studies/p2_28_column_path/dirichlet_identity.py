#!/usr/bin/env python3
"""P2-28 a5 -- THE DIRICHLET IDENTITY (the plumbing): Hu's deck on the column path vs on the face list.

  legs/hu_column/input   = the P2-27 deck as it was (site A, no face list): the COLUMN kernel's prescribed term
  validation/hu/hu_spall_onset/input_granite2_lefm = the P2-28 deck: the FACE LIST's T_presc top faces
Both on the stage-1 binary, whose two call sites carry the SAME corrected expression 2 kc (T_s - T_c)/dz^2 with the same
operands: nothing is removed on either, the losses are zero, so `Temp` must agree at every cell of every plotfile.
Bar: max |dT| <= 1e-9 K. Prediction: bit-identical (0). The factor itself is unit/face_dirichlet's (the exact solution).
Usage: dirichlet_identity.py [--check]     (writes legs/{hu_column,hu_face}/{run, argv.sh, BINARY.txt, run.log})
"""
import hashlib
import os
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
sys.path.insert(0, HERE)
import pf  # noqa: E402

BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
LEGS = {"hu_column": os.path.join(HERE, "legs", "hu_column", "input"),
        "hu_face": os.path.join(ROOT, "tests", "MMWSpalling", "validation", "hu", "hu_spall_onset",
                                "input_granite2_lefm")}
BAR = 1.0e-9


def run(tag):
    d = os.path.join(HERE, "legs", tag)
    os.makedirs(d, exist_ok=True)
    for p in (os.path.join(d, "run"),):
        if os.path.isdir(p):
            shutil.rmtree(p)
    for f in os.listdir(d):
        if f.startswith("run_") and f.endswith(".csv"):
            os.remove(os.path.join(d, f))
    argv = ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4", BIN, LEGS[tag],
            "plot_file=" + os.path.join(d, "run")]
    open(os.path.join(d, "argv.sh"), "w").write(" ".join(argv) + "\n")
    open(os.path.join(d, "BINARY.txt"), "w").write("sha256 %s\n" % hashlib.sha256(open(BIN, "rb").read()).hexdigest())
    with open(os.path.join(d, "run.log"), "w") as log:
        rc = subprocess.call(argv, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    print("  %-10s rc %d" % (tag, rc), flush=True)
    return rc


def main():
    if "--check" not in sys.argv:
        for tag in LEGS:
            if run(tag) != 0:
                print("FAIL: %s exited non-zero" % tag)
                return 1
    C = pf.plotfiles(os.path.join(HERE, "legs", "hu_column", "run") + "/")
    F = pf.plotfiles(os.path.join(HERE, "legs", "hu_face", "run") + "/")
    if not C or len(C) != len(F):
        print("FAIL: %d column vs %d face plotfiles" % (len(C), len(F)))
        return 1
    worst, n_nonzero = 0.0, 0
    for (sc, pc), (sf, pff) in zip(C, F):
        tc, a, _ = pf.load(pc, ["Temp"])
        tf, b, _ = pf.load(pff, ["Temp"])
        if sc != sf or tc != tf:
            print("FAIL: plotfile mismatch %s / %s" % (pc, pff))
            return 1
        d = float(np.nanmax(np.abs(a["Temp"] - b["Temp"])))
        n_nonzero += d > 0.0
        worst = max(worst, d)
    tl, tT, _ = pf.load(F[-1][1], ["Temp"])
    print("a5 -- Hu column vs Hu face list (stage-1 binary): %d plotfiles to t = %.1f s; max |dT| = %.3e K over all "
          "cells and times (%d plotfiles not bit-identical); the face run's top cell max %.2f K at the end"
          % (len(C), tl, worst, n_nonzero, float(np.nanmax(tT["Temp"][-1]))))
    ok = worst <= BAR
    print("%s: bar <= %.0e K; prediction bit-identical: %s" % ("PASS" if ok else "FAIL", BAR,
                                                              "MET" if worst == 0.0 else "NOT met"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
