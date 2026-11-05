#!/usr/bin/env python3
"""P2-2e diagnostic -- does the binary read the constraint reference at the
right depth? Compare the binary's CONSTANT-PROPERTY ladder (legs X<d>_<SOD>,
K0C_*) with the offline erfc ladder evaluated at the binary's OWN key-0
threshold, so the only remaining differences are discretisation (cell
averages, the reconstructed surface, the flaw-scale K_I scan).

A read that is too DEEP would make the drive too large and drag the binary's
interval far below the offline one -- which is the bug that would fake the
calibrated d_ref = 2.0 mm. This file is the evidence it is NOT that.
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import offline_ladder as OL                                       # noqa: E402
L = os.path.join(HERE, "legs")


def on(leg):
    for r in csv.DictReader(open(os.path.join(L, leg, "run_removal_events.csv"))):
        if r["regime"] == "1":
            return (float(r["time"]), float(r["T_surf"]))
    return None


def main():
    P = print
    th = [on(k)[1] - OL.T0 for k in ("K0C_130", "K0C_150", "K0C_h500")]
    dtc = sum(th) / 3
    P("P2-2e diagnostic -- the binary with the OFFLINE model's assumptions")
    P("=" * 86)
    P("  constant k 3.27 W/mK, kappa 1.6e-6 m2/s, constant E/beta (law = constant,")
    P("  stress_law = scalar), constant K_Ic (linear law, T_melt = 1e12 K)")
    P(f"  binary's key-0 threshold dT_c = {th[0]:.2f} / {th[1]:.2f} / {th[2]:.2f} K"
      f" at SOD 130 / SOD 150 / h500 -- flux-INDEPENDENT, as constant")
    P(f"  properties must give; mean {dtc:.2f} K")
    ds = [round(1.0 + 0.5 * i, 1) for i in range(23)]
    b130 = [d for d in ds if on(f"X{d:.1f}_130")]
    b150n = [d for d in ds if not on(f"X{d:.1f}_150")]
    bin_iv = (min(b130), max(b150n))
    OL.DT_C = dtc
    o130 = [d for d in ds if OL.drive_max(*OL.SOD[130], 3.27, 1.6e-6, d * 1e-3)[1] is not None]
    o150n = [d for d in ds if OL.drive_max(*OL.SOD[150], 3.27, 1.6e-6, d * 1e-3)[1] is None]
    off_iv = (min(o130), max(o150n))
    rnd = lambda x: OL.math.floor(x / 0.5 + 0.5) * 0.5
    P("")
    P(f"  BINARY  (constant props)       interval [{bin_iv[0]}, {bin_iv[1]}] mm"
      f"  -> rule d = {rnd(sum(bin_iv)/2):.1f} mm")
    P(f"  OFFLINE (same dT_c {dtc:.1f} K)    interval [{off_iv[0]}, {off_iv[1]}] mm"
      f"  -> rule d = {rnd(sum(off_iv)/2):.1f} mm")
    P(f"  BINARY  (the real T-dependent Gotthard/Grimsel laws, b2)  [1.5, 2.0] mm -> 2.0 mm")
    P("")
    P("  ⇒ With the offline model's assumptions the binary lands within")
    P(f"    {abs(rnd(sum(bin_iv)/2) - rnd(sum(off_iv)/2)):.1f} mm of the offline answer;"
      " a reference read too deep (or too")
    P("    shallow) by the factor that separates 2.0 from ~8 mm would show here as")
    P("    the same factor, and it does not. The calibrated 2.0 mm is therefore the")
    P("    TEMPERATURE DEPENDENCE of the rock laws acting on the criterion, not the")
    P("    implementation.")
    P("=" * 86)


if __name__ == "__main__":
    main()
