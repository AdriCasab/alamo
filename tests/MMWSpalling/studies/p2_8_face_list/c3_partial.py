#!/usr/bin/env python3
"""P2-8 c3 partial (DIAGNOSTIC, not the c3 row): L0f aborted at t = 27.296 s (Removal.H:3231),
so no 30 s plotfile exists. Volume from the event logs (sum n_voided x cell volume x 4, SR 3)
for events strictly before T, at T up to 27 s, for L0f / L0 / L0c. Max Ø and axis depth need a
plotfile and are NOT available for L0f."""
import csv, os
HERE = os.path.dirname(os.path.abspath(__file__))
out = ["P2-8 c3 partial -- event-log volume before T (x4); L0f aborted at 27.296 s, so c3's 30 s bars are unscoreable",
       "=" * 96]
legs = (("L0", 1e-3), ("L0f", 5e-4), ("L0c", 2e-3))
V = {}
for leg, dx in legs:
    ev = list(csv.DictReader(open(os.path.join(HERE, "legs", leg, "run_removal_events.csv"))))
    for T in (15.0, 20.0, 25.0, 27.0):
        sel = [e for e in ev if float(e["time"]) < T]
        nv = sum(int(float(e["n_voided"])) for e in sel)
        n6 = sum(int(float(e["n_voided"])) for e in sel if int(float(e["regime"])) == 6)
        V[leg, T] = nv * dx ** 3 * 4e6
        out.append(f"  {leg:4s} T {T:4.1f}  V {V[leg, T]:7.3f} cm3  lateral {100 * n6 / max(nv, 1):5.1f} %")
for T in (15.0, 20.0, 25.0, 27.0):
    out.append(f"  T {T:4.1f}  L0f/L0 - 1 = {100 * (V['L0f', T] / V['L0', T] - 1):+5.1f} %   "
               f"L0c/L0 - 1 = {100 * (V['L0c', T] / V['L0', T] - 1):+5.1f} %")
print("\n".join(out))
open(os.path.join(HERE, "C3_PARTIAL.txt"), "w").write("\n".join(out) + "\n")
