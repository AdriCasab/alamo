#!/usr/bin/env python3
"""Offline regime map, step 1: the temperature histories (independent of the stress law and of the flake rule).
P2-25's 7 x 6 grid (absorption depth D 0.5 ... 71 mm; absorbed flux 0.3 ... 8 MW/m2) on the 1-D column of
flake_rule_study_1111/study.py: T0 297.5 K, eps 0.8, h_conv 56, 50 mm deep (the power past it leaves), flaws off.
Snapshots at the onset (843 K), every 20 K above it, and 1410 K; 600 s cap. dz = D/16 clipped to [0.03, 0.2] mm."""
import os, pickle, sys
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "flake_rule_study_1111"))
import study as S
DS = (("0p5", 0.5e-3), ("1", 1e-3), ("2", 2e-3), ("4", 4e-3), ("8", 8e-3), ("16", 16e-3), ("71", 71.4e-3))
QS = (("0p3", 0.3), ("0p5", 0.5), ("1", 1.0), ("2", 2.0), ("4", 4.0), ("8", 8.0))
LEVELS = [843.0] + [860.0 + 20.0 * i for i in range(28)] + [1410.0]

def one(a):
    dk, D, qk, q = a
    dz = min(max(D / 16.0, 0.03e-3), 0.2e-3)
    out = S.column(297.5, 600.0, 2000.0, 2650.0, LEVELS, dz, 0.05, q=q * 1e6, alpha=1.0 / D, hc=56.0, Tamb=297.5)
    return ("D%s_Q%s" % (dk, qk), out)

if __name__ == "__main__":
    jobs = [(dk, D, qk, q) for dk, D in DS for qk, q in QS]
    with Pool(7) as p:
        res = dict(p.map(one, jobs, chunksize=1))
    pickle.dump(res, open(os.path.join(HERE, "thermal.pkl"), "wb"))
    for k in sorted(res):
        r = res[k]
        print("%-10s levels reached %2d of %d; onset %s s; last T_s %.0f K at %.1f s" % (
            k, len(r), len(LEVELS), ("%.2f" % r[0]["t"]) if r else "-", r[-1]["Ts"] if r else float("nan"), r[-1]["t"] if r else float("nan")))
