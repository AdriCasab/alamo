#!/usr/bin/env python3
"""P2-28 a3 (ii) -- SITE A vs THE FACE LIST on Kant's verification deck, in kelvin (stage 1 only: site A dies at stage 2).

The SAME statistic on every run (the old sp_kant_onset/test's onset_homogeneous): the first plotfile time at which the
top slab's MAXIMUM of the criterion field reaches 1, interpolated linearly between the bracketing plotfiles (0.5 s
apart), and the top cells' MEAN temperature there, dT = mean - T_I (293.15 K). The criterion field is Sp_field for site
A and fr_score for the face list.
  A_cell   legs/KV_A_cell:  the old input_verification exactly as recorded (site A at spall.surface_T = cell_average,
           the default it never set; the cell Robin form)
  F        the face deck's run (validation/kant/sp_kant_onset/output_face/verification, from test_face): fr_score
  F_siteA  the SAME run's Sp_field: site A evaluated every step on the same temperature field, at the same
           (reconstructed) surface temperature -- the face deck sets spall.surface_T = reconstructed
The GATE (declared, P2-28 stage 1): |dT(F) - dT(F_siteA)| <= 10 K (prediction <= 3 K) -- two implementations of one
rule on one input. ⚠ A separate site-A run at the reconstructed temperature (the old deck + spall.surface_T =
reconstructed) SEGFAULTS at step 1 on the P2-27 binary too (site A's reconstruction without surface.follow_mask; the
column path, pre-existing, dies at stage 2): F_siteA is that comparison without a second run.
Reported beside it: A_cell (the recorded deck, site A at the CELL-AVERAGE surface temperature and the cell Robin form:
its gap is the surface-temperature form, not the site) and the face list's event-row statistic (P2-13's).
Each is also scored against Kant Eq. 15 (445.6 K).
Usage: siteA_vs_face.py [--check]
"""
import csv
import hashlib
import math
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
KD = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "kant", "sp_kant_onset")
OLD = os.path.join(KD, "input_verification")
FACE = os.path.join(KD, "output_face", "verification")
T_I = 293.15
TH15 = 1.5e6 * (1.0 - 0.26) / (2.0 * 1.763 * math.sqrt(20.0e-6 / math.pi) * 35.0e9 * 8.0e-6)
GATE, PRED = 10.0, 3.0
LEGS = {"KV_A_cell": []}


def run(tag, extra):
    d = os.path.join(HERE, "legs", tag)
    if os.path.isdir(d):
        shutil.rmtree(d)
    os.makedirs(d)
    argv = ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4", BIN, OLD] + extra + \
           ["plot_file=" + os.path.join(d, "run")]
    open(os.path.join(d, "argv.sh"), "w").write(" ".join(argv) + "\n")
    open(os.path.join(d, "BINARY.txt"), "w").write("sha256 %s\n" % hashlib.sha256(open(BIN, "rb").read()).hexdigest())
    with open(os.path.join(d, "run.log"), "w") as log:
        rc = subprocess.call(argv, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    print("  %-10s rc %d" % (tag, rc), flush=True)
    return rc


def onset(run_dir, field):
    prev = None
    for _, p in pf.plotfiles(run_dir + "/"):
        t, a, _ = pf.load(p, [field, "Temp"])
        top_f, top_T = a[field][-1], a["Temp"][-1]
        mx, dT = float(np.nanmax(top_f)), float(np.nanmean(top_T)) - T_I
        if mx >= 1.0:
            if prev is None:
                return t, dT
            t0, mx0, dT0 = prev
            f = (1.0 - mx0) / (mx - mx0) if mx > mx0 else 1.0
            return t0 + f * (t - t0), dT0 + f * (dT - dT0)
        prev = (t, mx, dT)
    return None, None


def event_stat(path):
    first = {}
    for r in csv.DictReader(open(path)):
        if float(r["ts_fb"]) == 0.0 and r["dir"] == "5":
            first.setdefault((r["col_i"], r["col_j"]), float(r["T_surf"]))
    return float(np.median(list(first.values()))) - T_I if first else None


def main():
    if "--check" not in sys.argv:
        for tag, extra in LEGS.items():
            if run(tag, extra) != 0:
                print("FAIL: %s exited non-zero" % tag)
                return 1
    rows = {}
    rows["A_cell"] = onset(os.path.join(HERE, "legs", "KV_A_cell", "run"), "Sp_field")
    rows["F"] = onset(FACE, "fr_score")
    rows["F_siteA"] = onset(FACE, "Sp_field")
    ev = event_stat(FACE + "_removal_events.csv")
    print("a3 (ii) -- site A vs the face list, Kant verification (h 150): onset dT by ONE statistic "
          "(top-slab max criterion >= 1, interpolated; dT = top-cell mean - 293.15 K); Eq. 15 = %.1f K" % TH15)
    for k in ("A_cell", "F", "F_siteA"):
        t, d = rows[k]
        print("  %-8s t %s  dT %s  (vs Eq. 15 %s)" % (k, "%.2f s" % t if t is not None else "never",
                                                     "%.2f K" % d if d is not None else "-",
                                                     "%+.2f %%" % (100 * (d - TH15) / TH15) if d is not None else "-"))
    print("  the face list's event-row statistic (P2-13: median first top-face T_surf - T_I): %s"
          % ("%.2f K (%+.2f %% vs Eq. 15)" % (ev, 100 * (ev - TH15) / TH15) if ev is not None else "none"))
    if rows["F"][1] is None or rows["F_siteA"][1] is None:
        print("FAIL: an onset is missing")
        return 1
    gap = rows["F"][1] - rows["F_siteA"][1]
    print("  GATE |dT(F) - dT(F_siteA)| = %.2f K (bar %.0f K, prediction %.0f K): %s" % (
        abs(gap), GATE, PRED, "PASS" if abs(gap) <= GATE else "FAIL"))
    if rows["A_cell"][1] is not None:
        print("  reported: dT(F) - dT(A_cell) = %+.2f K (the recorded deck: site A at the CELL-AVERAGE surface "
              "temperature, the cell Robin form)" % (rows["F"][1] - rows["A_cell"][1]))
    return 0 if abs(gap) <= GATE else 1


if __name__ == "__main__":
    sys.exit(main())
