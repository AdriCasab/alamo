#!/usr/bin/env python3
"""P2-9 finding: the removal step leaves temp_mf's ghosts STALE for the kernel.

After a firing the remap writes H, InvertHtoT recovers T on VALID cells, and
Util::RealFillBoundary refills ghosts -- but RealFillBoundary is a bare
FillBoundary(): box-to-box only, NOT the periodic copy and NOT the physical
boundary condition. The integrator refills those at the start of the NEXT step,
so the conduction kernel of the firing step reads the PRE-remap temperature
across every domain side (and across a periodic side).

Repro: the unit/face_removal fixture (uniform gas, uniform h, a flat floor --
every column identical), one box, x/y periodic, to the first firing (t = 2.44 s,
step 306), a plotfile every step. A correct step keeps every column identical.
  key 0 (the column path, unchanged)        -> row/column 0 differ after the firing
  key 1 (the line path, FrFillTempGhosts)   -> every column identical
Writes GHOST_REPRO.txt.
"""
import os
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "tests", "MMWSpalling", "unit", "face_removal"))
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
INPUT = os.path.join(ROOT, "tests", "MMWSpalling", "unit", "face_removal", "input")
OUT = os.path.join(HERE, "ghost")


def load_reader():
    import importlib.machinery
    import importlib.util
    p = os.path.join(ROOT, "tests", "MMWSpalling", "unit", "face_removal", "test")
    ld = importlib.machinery.SourceFileLoader("frtest", p)
    m = importlib.util.module_from_spec(importlib.util.spec_from_loader("frtest", ld))
    ld.exec_module(m)
    return m.read_pf


def main():
    read_pf = load_reader()
    lines = [__doc__.strip().split("\n")[0], "=" * 88]
    for key in (0, 1):
        d = os.path.join(OUT, "k%d" % key)
        if os.path.isdir(d):
            shutil.rmtree(d)
        os.makedirs(d)
        argv = [BIN, INPUT, "surface_patch.face_list=%d" % key, "amr.max_grid_size=64",
                "stop_time=2.456", "amr.plot_int=1", "geometry.is_periodic=1 1 0",
                "plot_file=" + os.path.join(d, "run")]
        open(os.path.join(d, "argv.sh"), "w").write(
            " ".join("'%s'" % a if " " in a else a for a in argv) + "\n")
        with open(os.path.join(d, "run.log"), "w") as log:
            rc = subprocess.call(argv, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
        lines.append("key %d  rc %d  (%s)" % (key, rc, "column path" if key == 0 else "line path"))
        for st in (305, 306, 307):
            T = read_pf(os.path.join(d, "run", "%05dcell" % st), ["Temp"], (24, 24, 64))["Temp"]
            sp = [T[:, :, k].max() - T[:, :, k].min() for k in (63, 62, 61, 50)]
            lines.append("   step %d  max - min over columns at k 63/62/61/50: %s K"
                         % (st, "  ".join("%.6f" % v for v in sp)))
        T = read_pf(os.path.join(d, "run", "00306cell"), ["Temp"], (24, 24, 64))["Temp"]
        kk = 62 if key == 0 else 63
        lines.append("   step 306, k %d: T - min at (i, j) = (0, 0) %.4f  (0, 12) %.4f  (12, 0) %.4f  "
                     "(12, 12) %.4f  (23, 23) %.4f K" % (
                         kk, *(T[i, j, kk] - T[:, :, kk].min()
                               for i, j in ((0, 0), (0, 12), (12, 0), (12, 12), (23, 23)))))
    # keep the three plotfiles the table reads; the other ~600 are 6.7 GB
    for key in (0, 1):
        run = os.path.join(OUT, "k%d" % key, "run")
        for nm in os.listdir(run):
            if nm[:5].isdigit() and nm[:5] not in ("00305", "00306", "00307"):
                shutil.rmtree(os.path.join(run, nm))
    out = "\n".join(lines)
    print(out)
    open(os.path.join(HERE, "GHOST_REPRO.txt"), "w").write(out + "\n")


if __name__ == "__main__":
    main()
