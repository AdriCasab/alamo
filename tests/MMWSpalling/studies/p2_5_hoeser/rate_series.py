#!/usr/bin/env python3
"""P2-5 rev 2 (c10, c11) -- the transient rate series and the L 13 ratio.

Volume vs time from the EVENT LOG (sum n_voided * V_cell, x4). The rate at a
time t is the volume accrued in the preceding window, expressed as cm3/min so
it sits beside Hoeser's 47 MPa 64/66/81/90 cm3/min at 15/30/60/90 s.
⚠ Her numbers are at 47 MPa; this deck is at 0 MPa, where she measured
26 cm3/min at 30 s. The SHAPE of the series is what c11 scores, not its level.
"""
import csv
import os
import shlex
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def vcell(leg):
    toks = shlex.split(open(os.path.join(HERE, "legs", leg, "argv.sh")).read())
    kv, i = {}, 0
    while i < len(toks):
        t = toks[i]
        if "=" in t:
            k, v = t.split("=", 1)
            i += 1
            while i < len(toks) and "=" not in toks[i]:
                v += " " + toks[i]
                i += 1
            kv[k] = v
        else:
            i += 1
    nx, _, nz = (int(v) for v in kv["amr.n_cell"].split())
    lo = [float(v) for v in kv["geometry.prob_lo"].split()]
    hi = [float(v) for v in kv["geometry.prob_hi"].split()]
    dx, dz = (hi[0] - lo[0]) / nx, (hi[2] - lo[2]) / nz
    return dx * dx * dz, dz


def vt(leg):
    """[(t, V_cm3_x4, deepest_mm)] cumulative, from the event log."""
    V = os.path.join(HERE, "legs", leg, "run_removal_events.csv")
    vc, dz = vcell(leg)
    out, tot, depth = [], 0, {}
    for r in csv.DictReader(open(V)):
        try:
            t = float(r["time"]); nv = int(float(r["n_voided"] or 0))
            i, j = int(float(r["col_i"])), int(float(r["col_j"]))
        except (ValueError, KeyError):
            continue
        tot += nv
        depth[(i, j)] = depth.get((i, j), 0) + nv
        out.append((t, tot * vc * 4.0 * 1e6, max(depth.values()) * dz * 1e3))
    return out


def at(series, t):
    v, d = 0.0, 0.0
    for tt, vv, dd in series:
        if tt <= t:
            v, d = vv, dd
        else:
            break
    return v, d


def main():
    print("P2-5 rev 2 (c10, c11) -- the transient and the standoff ratio")
    print("=" * 86)
    s120 = vt("G05T120")
    print("  c11  G05T120 (shape G05, K30, 120 s, 1 mm)")
    print(f"    {'t [s]':>6s} {'V cum [cm3]':>12s} {'rate [cm3/min]':>15s}"
          f" {'deepest [mm]':>13s}   her 47 MPa")
    her = {15: 64, 30: 66, 60: 81, 90: 90, 120: None}
    prev_t = prev_v = 0.0
    for t in (15, 30, 60, 90, 120):
        v, d = at(s120, t)
        rate = (v - prev_v) / (t - prev_t) * 60.0
        h = her[t]
        print(f"    {t:6d} {v:12.3f} {rate:15.2f} {d:13.1f}   "
              f"{('-' if h is None else str(h)):>5s}")
        prev_t, prev_v = t, v
    v15, _ = at(s120, 15)
    v90, _ = at(s120, 90)
    r15 = v15 / 15 * 60
    r90 = (v90 - at(s120, 60)[0]) / 30 * 60
    print(f"    rise from the 0-15 s rate to the 60-90 s rate:"
          f" {r90/max(r15,1e-9):.3f}x   (c11 prediction: rising by <= 1.5x)")
    print()
    v13, d13 = at(vt("G05L13"), 30)
    v8, d8 = at(vt("G05"), 30)
    print("  c10  G05L13 / G05 at 30 s -- the standoff exponent")
    print(f"    G05   (L 8,  s0 = 63.2 mm)  V = {v8:8.3f} cm3,"
          f" deepest {d8:6.1f} mm")
    print(f"    G05L13 (L 13, s0 = 102.7 mm) V = {v13:8.3f} cm3,"
          f" deepest {d13:6.1f} mm")
    print(f"    ratio {v13/max(v8,1e-9):.4f}   vs her 47 MPa 0.46 and her flux"
          f" ratio 0.58; prediction 0.45-0.60 (K30 law 0.53 at 850 K)")
    print("=" * 86)


if __name__ == "__main__":
    main()
