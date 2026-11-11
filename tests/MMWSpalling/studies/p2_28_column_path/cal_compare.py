#!/usr/bin/env python3
"""P2-28 b5 -- THE CALIBRATION LEG ACROSS THE DELETION: CAL (60 s) on the stage-1 binary (cal_before/) vs the stage-2
binary (validation/kant/grimsel_calibration/legs/CAL).
  gate   the event log byte-identical; thermo.dat: the declared columns gone (recon_*, surf_cols_robin, surf_q_mean,
         patch_P_robin), ledger_err renamed ledger_err_raw with identical values, every other column identical row by
         row (the wall-clock timers excluded)
  report the c1 level -- grimsel_calibration/reduce.py's statistic (per column the first dir-5, ts_fb-0 row's T_surf,
         the median) -- beside P2-13's banked 850.43 K (it predates P2-17 ... P2-27 and may have drifted for reasons
         outside this packet)
Usage: cal_compare.py
"""
import csv
import hashlib
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
T = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
import audit_components as AC  # noqa: E402

A = os.path.join(HERE, "cal_before")
B = os.path.join(T, "validation", "kant", "grimsel_calibration", "legs", "CAL")
BANKED = 850.43


def md5(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest()


def level(path):
    first = {}
    for r in csv.DictReader(open(path)):
        if float(r["ts_fb"]) == 0.0 and r["dir"] == "5":
            first.setdefault((r["col_i"], r["col_j"]), float(r["T_surf"]))
    v = np.array(list(first.values()))
    return float(np.median(v)), len(v), float(np.mean(v > 873.15))


def main():
    bad = []
    ea, eb = os.path.join(A, "run_removal_events.csv"), os.path.join(B, "run_removal_events.csv")
    same_ev = md5(ea) == md5(eb)
    if not same_ev:
        bad.append("the event log differs")
    st, gone = AC.thermo_audit(os.path.join(A, "run", "thermo.dat"), os.path.join(B, "run", "thermo.dat"))
    if st != "ok":
        bad.append("thermo: " + st)
    la, na, fa = level(ea)
    lb, nb, fb = level(eb)
    print("b5 -- CAL (60 s), stage-1 binary (cal_before/) vs stage-2 binary (legs/CAL)")
    print("  event log: %s (%d rows)" % ("BYTE-IDENTICAL" if same_ev else "DIFFERS",
                                         sum(1 for _ in open(eb)) - 1))
    print("  thermo.dat: %s; columns gone %s; ledger_err -> ledger_err_raw" % (st, gone))
    print("  c1 level (median first top-face T_surf over %d columns): %.2f K before, %.2f K after; P2-13 banked %.2f K "
          "(diff %+.2f K; reported, not a gate); share above the 873.15 K table end %.1f %%"
          % (nb, la, lb, BANKED, lb - BANKED, 100 * fb))
    print("GATE: %s" % ("PASS" if not bad else "FAIL: " + "; ".join(bad)))
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
