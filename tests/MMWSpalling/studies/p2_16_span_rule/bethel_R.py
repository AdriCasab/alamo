#!/usr/bin/env python3
"""P2-16 ord 1 -- R_30: the stress relief that moves the deterministic onset by DELTA = -30 K (Bethel's offset from
Gotthard read off Kant 2016 Fig. 11/13 [fig]; the c1 SENSITIVITY only). The P2-13 offline twin's onset(a, p, r) is
imported (validation/kant/grimsel_calibration/offline_onset.py), never copied; a 8 um, p 0. R enters the drive only
(not C_L, not ell; ledger L5). Writes R_BETHEL.txt, read by run.py at launch of the *_m30 legs."""
import importlib.machinery
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
P = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "kant", "grimsel_calibration", "offline_onset.py")
ld = importlib.machinery.SourceFileLoader("oo16", P)
OO = importlib.util.module_from_spec(importlib.util.spec_from_loader("oo16", ld))
ld.exec_module(OO)
DELTA = -30.0


def main():
    t0 = OO.onset(OO.A0, 0.0, OO.R)
    want = t0 + DELTA
    lo, hi = OO.R, 4.0 * OO.R                      # a larger R fires cooler
    assert OO.onset(OO.A0, 0.0, hi) < want < t0
    for _ in range(80):
        m = 0.5 * (lo + hi)
        lo, hi = (m, hi) if OO.onset(OO.A0, 0.0, m) > want else (lo, m)
    r = 0.5 * (lo + hi)
    got = OO.onset(OO.A0, 0.0, float("%.6g" % r))
    L = ["%.6g  stress_relief for DELTA %+.0f K  (twin: onset %.3f K at R %.6g -> %.3f K at R %.6g, target %.3f K; "
         "ratio R_30/R %.4f)" % (r, DELTA, t0, OO.R, got, r, want, r / OO.R)]
    open(os.path.join(HERE, "R_BETHEL.txt"), "w").write("\n".join(L) + "\n")
    print(L[0])
    return 0 if abs(got - want) < 0.05 else 1


if __name__ == "__main__":
    sys.exit(main())
