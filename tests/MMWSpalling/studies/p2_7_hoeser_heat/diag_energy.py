#!/usr/bin/env python3
"""P2-7 diagnostic (falls out of the runs already made; no new leg).

Where each leg's absorbed energy enters, and how the removal rate evolves.
  floor map   int patch_P_robin dt   (the h_expr Robin power on top faces)
  side faces  int side_P_face dt     (riser / wall faces)
  jet face    int jet_P_face dt      (floor + side as the closure books it)
  absorbed    ledger_E_src, CUMULATIVE, final value (authoritative)
⚠ patch_P_robin / side_P_face / jet_P_face are PER-STEP powers printed every
  amr.thermo.plot_int = 6 steps, so an integral is sum x (dt x 6) -- a sampled
  estimate, never a total by summation of the raw column (SR 7 i). It is
  checked against ledger_E_src, which is exact. All energies x4 (SR 3).
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import reduce_p2_7 as R                         # noqa: E402

WIN = 5.0


def main():
    legs = sys.argv[1:] or ["HA", "HG", "HGd"]
    print("P2-7 diagnostic -- the absorbed energy by face kind, and the rate by window")
    print("=" * 100)
    print(f"  {'leg':>4s} | {'floor map':>10s} {'side faces':>10s} {'jet face':>10s}"
          f" {'E_src':>10s} [kJ x4] | {'(floor+side)/E_src':>18s} {'side share':>10s}")
    out = {}
    for n in legs:
        d = R.legdir(n)
        g = R.geom_of(d)
        th = R.thermo(d)
        dts = g["dt"] * g["tint"]
        Ef = sum(th["patch_P_robin"]) * dts * 4.0
        Es = sum(th["side_P_face"]) * dts * 4.0
        Ej = sum(th["jet_P_face"]) * dts * 4.0
        Esrc = th["ledger_E_src"][-1] * 4.0
        out[n] = (th, g)
        print(f"  {n:>4s} | {Ef/1e3:10.3f} {Es/1e3:10.3f} {Ej/1e3:10.3f} {Esrc/1e3:10.3f}"
              f"         | {(Ef+Es)/Esrc:18.4f} {Es/(Ef+Es)*100:9.1f} %")
    print("\n  FLOOR-MAP POWER patch_P_robin [W, one quadrant] at selected times")
    ts = (0.048, 1.0, 5.0, 10.0, 20.0, 29.9)
    print("  " + f"{'leg':>4s} | " + " ".join(f"{t:>8.2f}s" for t in ts))
    for n in legs:
        th, g = out[n]
        row = []
        for t in ts:
            i = min(range(len(th["time"])), key=lambda q: abs(th["time"][q] - t))
            row.append(th["patch_P_robin"][i])
        print(f"  {n:>4s} | " + " ".join(f"{v:9.2f}" for v in row))
    print(f"\n  REMOVAL RATE by {WIN:g} s window [cm3/min, x4], from the event log (n_voided)")
    for n in legs:
        d = R.legdir(n)
        g = R.geom_of(d)
        ev = R.events(d)
        nb = int(round(g["st"] / WIN))
        b = [0] * nb
        for e in ev:
            t = R._f(e, "time")
            b[min(nb - 1, int(t / WIN))] += R._i(e, "n_voided")
        print(f"  {n:>4s} | " + " ".join(f"{x*g['V_cell']*4e6/WIN*60:7.2f}" for x in b))
    print("=" * 100)


if __name__ == "__main__":
    main()
