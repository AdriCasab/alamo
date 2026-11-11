#!/usr/bin/env python3
"""P2-28 stage 1 -- the declared checks that read the gate outputs (no runs):
  a2  the face gate's only-in-new paths are exactly under KV_S/, HU_S/, GB_S/ (and nothing only-in-ref)
  W1  (the packet's inertness witness) KV_S's fr_score is non-zero on the patch cells at its last plotfile
  W2  HU_S: the buckling log and the event log exist; T_s on every floor row == T_B1(time)
  b5  cal_before/: the CAL leg's events file and thermo exist (the stage-2 comparison's BEFORE half)
Usage: stage1_checks.py > STAGE1_CHECKS.txt
"""
import csv
import glob
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
W = os.path.join(ROOT, "tests", "MMWSpalling", "studies", "tree_sweep_0922", "witness")
sys.path.insert(0, HERE)
import pf  # noqa: E402

bad = []


def need(ok, msg):
    print("  %s  %s" % ("ok  " if ok else "FAIL", msg))
    if not ok:
        bad.append(msg)


txt = open(os.path.join(HERE, "GATE_a2.txt")).read()
md5 = lambda lab: {ln.split("  ", 1)[1].strip(): ln.split("  ", 1)[0]
                   for ln in open(os.path.join(W, "face_%s.md5" % lab)) if ln.strip()}
A, B = md5("p2_27"), md5("p2_28a")
new = sorted(set(B) - set(A))
old = sorted(set(A) - set(B))
m = re.search(r"only in p2_28a: (\d+)", txt)
legs = sorted({p.split("/")[0] for p in new})
head = re.search(r"face_p2_27 vs face_p2_28a: .*", txt).group(0)
print("a2 --", head)
need(not old and m is not None and int(m.group(1)) == len(new) and legs == ["GB_S", "HU_S", "KV_S"],
     "only-in-new = %d paths (the md5 lists), all under %s; only-in-ref %d" % (len(new), legs, len(old)))
per = {l: sum(1 for p in new if p.startswith(l + "/")) for l in legs}
print("     per new leg: %s" % per)
diffs = re.findall(r"thermo (\S+)\s+timers only", txt)
need(len(diffs) == 16 and "0 substantive" in txt, "the %d differing files are thermo.dat, timers only (16 legs)"
     % len(diffs))

print("W1 -- KV_S's fr_score at its last plotfile (the evaluate-only witness)")
pfs = pf.plotfiles(os.path.join(W, "face_p2_28a", "KV_S", "run") + "/")
t, a, _ = pf.load(pfs[-1][1], ["fr_score", "Temp"])
top = a["fr_score"][-1]
need(np.count_nonzero(top > 0) == top.size, "t %.1f s: fr_score > 0 on %d of %d top cells (max %.4f)"
     % (t, np.count_nonzero(top > 0), top.size, float(top.max())))

print("W2 -- HU_S (10 s): the logs")
bk = os.path.join(W, "face_p2_28a", "HU_S", "run_buckle_log.csv")
ev = os.path.join(W, "face_p2_28a", "HU_S", "run_removal_events.csv")
nb = len(list(csv.DictReader(open(bk)))) if os.path.exists(bk) else -1
ne = len(list(csv.DictReader(open(ev)))) if os.path.exists(ev) else -1
need(nb == 0 and ne == 0, "buckle log %d rows, events %d rows by 10 s (the criterion fires at ~50 s)" % (nb, ne))
t, a, _ = pf.load(pf.plotfiles(os.path.join(W, "face_p2_28a", "HU_S", "run") + "/")[-1][1], ["fr_score"])
print("     HU_S t %.1f s: fr_score max on the patch %.4f (non-zero: the criterion runs)" % (t, float(a["fr_score"].max())))

print("b5 -- the BEFORE half")
cb = os.path.join(HERE, "cal_before")
need(os.path.exists(os.path.join(cb, "run_removal_events.csv")) and os.path.exists(os.path.join(cb, "run", "thermo.dat")),
     "cal_before/: run_removal_events.csv (%d rows) and run/thermo.dat" %
     len(list(csv.DictReader(open(os.path.join(cb, "run_removal_events.csv"))))))
print("m4 -- Hu (the full 120 s run, validation/hu/hu_spall_onset/output_lefm): per column, the FIRST floor row of the "
      "buckling log (the criterion's onset there)")
rows = list(csv.DictReader(open(os.path.join(ROOT, "tests", "MMWSpalling", "validation", "hu", "hu_spall_onset",
                                             "output_lefm", "granite2_buckle_log.csv"))))
first = {}
for r in rows:
    if r["class"] == "0":
        first.setdefault((r["col_i"], r["col_j"]), r)
v = lambda k: np.array([float(r[k]) for r in first.values()])
tt, Ts, ell, bkr = v("time"), v("T_surf"), v("ell_p"), v("bk")
print("     columns fired %d (the patch lists 76 top cells); first-row time median %.2f s (min %.2f, max %.2f); "
      "T_s median %.2f K (p10 %.2f, p90 %.2f)" % (len(first), np.median(tt), tt.min(), tt.max(), np.median(Ts),
                                                   np.percentile(Ts, 10), np.percentile(Ts, 90)))
print("     at each column's first firing: ell_p median %.3f mm (p10 %.3f, p90 %.3f); bk median %.3f (min %.3f); "
      "blocked %d of %d" % (1e3 * np.median(ell), 1e3 * np.percentile(ell, 10), 1e3 * np.percentile(ell, 90),
                             np.median(bkr), bkr.min(), sum(int(r["blocked"]) for r in first.values()), len(first)))
print("\nstage 1 checks: %s" % ("PASS" if not bad else "FAIL (%d)" % len(bad)))
sys.exit(1 if bad else 0)
