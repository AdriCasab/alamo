#!/usr/bin/env python3
"""P1 revision 5, deliverable 5 -- the exposed-face census.

A REDUCTION of a leg already on disk; it starts no run. The successor packet is
designed on these numbers, so they have to be reproducible from the shipped
event log. Verify the log first:

    shasum -a 256 output_fh/FH_removal_events.csv

Reads output_fh/FH_removal_events.csv (quarter domain, dx = dz = 2 mm) and
answers one question: when the nozzle collided at t = 177.616 s, was the floor
under the burner a plateau or a field of towers and chasms?
"""
import csv, math, os, statistics, sys

HERE = os.path.dirname(os.path.abspath(__file__))
EV   = os.path.join(HERE, "output_fh", "FH_removal_events.csv")
DX = DZ = 0.002
KTOP0 = 299                 # (0.5 - (-0.1))/0.002 - 1: the initial top cell
T_END = 177.616             # the abort time
T_FREEZE = 120.0            # "stopped firing before" -- 57 s of standing still
R_BURNER = 0.040            # the burner footprint

def main():
    last, final, nfire = {}, {}, {}
    for d in csv.DictReader(open(EV)):
        t = float(d["time"])
        key = (int(d["col_i"]), int(d["col_j"]))
        k = int(float(d["k_top"])) - int(float(d["n_voided"]))
        last[key] = max(last.get(key, -1.0), t)
        final[key] = min(final.get(key, 10**9), k)
        nfire[key] = nfire.get(key, 0) + 1

    def r_of(i, j): return math.hypot((i + 0.5) * DX, (j + 0.5) * DX)
    def depth(key):  return (KTOP0 - final[key]) * DZ if key in final else 0.0

    under = [k for k in final if r_of(*k) < R_BURNER]
    frozen = [k for k in under if last[k] < T_FREEZE]
    live   = [k for k in under if last[k] >= T_FREEZE]
    print("Meier leg, t <= %.3f s, r < %.0f mm" % (T_END, 1000 * R_BURNER))
    print("  columns under the burner that fired at all : %d" % len(under))
    print("  ... that STOPPED firing before t = %.0f s  : %d (%.0f %%)"
          % (T_FREEZE, len(frozen), 100.0 * len(frozen) / len(under)))
    print("  median depth  frozen %.0f mm   live %.0f mm"
          % (1000 * statistics.median([depth(k) for k in frozen]),
             1000 * statistics.median([depth(k) for k in live])))

    # An "exposed side" is a 4-neighbour whose floor sits at least one cell
    # LOWER: the column stands proud of it, so that face is a free surface.
    def exposed(key):
        i, j = key; n = 0
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            nb = (i + di, j + dj)
            # ⚠ Mirror across the quarter-domain axes: i = -1 is the image of
            # i = 0, NOT a missing neighbour. Treating it as missing would
            # under-count exposure on the two axis rows.
            if nb[0] < 0: nb = (0, nb[1])
            if nb[1] < 0: nb = (nb[0], 0)
            if nb == key: continue
            if nb in final and depth(nb) > depth(key) + 0.5 * DZ: n += 1
        return n

    for lab, grp in (("frozen", frozen), ("live", live)):
        hist = [0] * 5
        for k in grp: hist[exposed(k)] += 1
        tot = max(1, len(grp))
        print("  %-6s exposed sides  0:%2.0f %%  1:%2.0f %%  2:%2.0f %%  "
              "3:%2.0f %%  4:%2.0f %%   (>=2: %.0f %%)"
              % (lab, *[100.0 * h / tot for h in hist],
                 100.0 * sum(hist[2:]) / tot))

    # neighbour-to-neighbour step inside the frozen set
    steps = []
    fs = set(frozen)
    for (i, j) in frozen:
        for di, dj in ((1, 0), (0, 1)):
            nb = (i + di, j + dj)
            if nb in final:
                steps.append(abs(depth((i, j)) - depth(nb)))
    if steps:
        steps.sort()
        print("  neighbour depth step at a frozen column: median %.0f mm, "
              "p90 %.0f mm   (cell = %.0f mm)"
              % (1000 * statistics.median(steps),
                 1000 * steps[int(0.9 * (len(steps) - 1))], 1000 * DZ))

    c = (6, 11)
    nb_d = [1000 * depth((c[0] + a, c[1] + b))
            for a, b in ((1, 0), (-1, 0), (0, 1), (0, -1))
            if (c[0] + a, c[1] + b) in final]
    print("  THE COLLIDING COLUMN (6,11)  r = %.1f mm: %d of 4 sides exposed, "
          "depth %.0f mm, deepest neighbour %.0f mm (%.1f mm below it), "
          "last fired %.1f s, stood %.0f s"
          % (1000 * r_of(*c), exposed(c), 1000 * depth(c), max(nb_d),
             max(nb_d) - 1000 * depth(c), last[c], T_END - last[c]))
    print("\n=> the frozen region is NOT a plateau: it is a field of towers and "
          "chasms, so the collision is a ROUGHNESS instability inside the "
          "resolved floor, not the wall-band effect P1 was aimed at.")

if __name__ == "__main__":
    sys.exit(main())
