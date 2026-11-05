#!/usr/bin/env python3
"""P2-9 rev 2 (a5 diagnostic): is the tilt legs' slow recession REAL, or the
plotfile measurement? Two independent routes over 5-20 s, the same interior
columns (>= 4 cells from every side), for T00 key 1 and the z-dominant tilts:
  EVENTS   sum over firings whose cell sits in an interior column of the volume
           the line removed, h_line dx^2 (h_line = solid fraction consumed along
           the line), per column, x cos(theta_true) / (15 s)
  FIELD    the decrease of sum(fr_f) dx^3 over the same columns between the 5 s
           and 20 s plotfiles, per column, x cos(theta_true) / (15 s)
Plus the per-firing normal recession h_line cos(theta_true) and the firing count,
against the floor's. Writes TILT_CROSSCHECK.txt."""
import csv
import math
import os
import importlib.machinery
import importlib.util

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
T = os.path.join(ROOT, "tests", "MMWSpalling", "unit", "face_removal", "test")
ld = importlib.machinery.SourceFileLoader("frtest", T)
m = importlib.util.module_from_spec(importlib.util.spec_from_loader("frtest", ld))
ld.exec_module(m)
OUT = m.OUT


def main():
    lines = [__doc__.strip().split("\n")[0], "=" * 100]
    legs = [("T00_1_k1", 0.0, 0.0)] + [(m.leg_tilt(nm, "1"), th, ph) for nm, (th, ph) in m.TILTS.items()
                                       if nm != "TL60"]
    ref = None
    for lg, th, ph in legs:
        n = m.normal(th, ph)
        ct = abs(n[2])
        ev = [e for e in m.events(lg) if 5.0 <= float(e["time"]) < 20.0]
        inner = lambda e: 4 <= int(e["col_i"]) <= 19 and 4 <= int(e["col_j"]) <= 19
        evi = [e for e in ev if inner(e)]
        ncol = 16 * 16
        r_ev = sum(float(e["h_line"]) for e in evi) * ct / ncol / 15.0
        pfs = m.plotfiles(lg)
        byt = {round(m.pf_time(p)): p for p in pfs}
        f5 = m.read_pf(byt[5], ["fr_f"], (24, 24, 64))["fr_f"][4:20, 4:20, :].sum()
        f20 = m.read_pf(byt[20], ["fr_f"], (24, 24, 64))["fr_f"][4:20, 4:20, :].sum()
        r_f = (f5 - f20) * 0.001 * ct / ncol / 15.0
        per = np.mean([float(e["h_line"]) for e in evi]) * ct
        if ref is None:
            ref = (r_ev, r_f, per, len(evi))
        lines.append("  %-9s theta %4.1f: EVENTS %.4e  FIELD %.4e m/s (vs floor %+.4f / %+.4f);  firings %5d "
                     "(%+.3f);  per-firing normal recession %.4e m (%+.4f)" % (
                         lg, th, r_ev, r_f, r_ev / ref[0] - 1, r_f / ref[1] - 1, len(evi),
                         len(evi) / ref[3] - 1, per, per / ref[2] - 1))
    txt = "\n".join(lines)
    print(txt)
    open(os.path.join(HERE, "TILT_CROSSCHECK.txt"), "w").write(txt + "\n")


if __name__ == "__main__":
    main()
