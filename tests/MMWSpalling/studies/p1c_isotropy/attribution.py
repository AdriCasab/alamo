#!/usr/bin/env python3
"""P1c -- attribute any residual azimuthal anisotropy to the RULE or the
CRITERION. A reduction; starts no run.

The question this answers, and it is the only one that decides how to read (f):

  if P1c's remaining per-sector spread lives in the REGIME-1 (criterion)
  removals, it was there all along and P1b's diagonal over-drilling was MASKING
  it -- which is exactly what P5 pre-registered. If it lives in the REGIME-5
  (support) removals, the rule is still the anisotropic thing and `min_sides = 3`
  is the wrong definition.

Both legs are read at their FINAL state, so the comparison is matched.
"""
import csv, math, os, statistics, sys

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)
DZ = DX = 0.002
KTOP0 = 299
SECTORS = [0, 15, 30, 45, 60, 75]
BANDS = (("floor  r < 25", 0.000, 0.025),
         ("mid    25-35 ", 0.025, 0.035),
         ("ring   35-44 ", 0.035, 0.044),
         ("wall   44-50 ", 0.044, 0.050))
LEGS = (("P1b  min_sides=2", os.path.join(STUD, "p1b_support_rule",
                                          "output_fh")),
        ("P1c  min_sides=3", os.path.join(HERE, "output_fh")))


def r_of(i, j):
    return math.hypot((i + 0.5) * DX, (j + 0.5) * DX)


def sec(i, j):
    return int(math.degrees(math.atan2(j + 0.5, i + 0.5)) // 15) * 15


def load(out):
    p = os.path.join(out, "FH_removal_events.csv")
    if not os.path.exists(p):
        return None
    ev = []
    for d in csv.DictReader(open(p)):
        ev.append((int(d["col_i"]), int(d["col_j"]),
                   int(float(d["regime"])), int(float(d["n_voided"]))))
    return ev


def main():
    for lab, out in LEGS:
        ev = load(out)
        if ev is None:
            print("%s: MISSING %s\n" % (lab, out)); continue
        print("=" * 78)
        print("%s   %s" % (lab, os.path.relpath(out, STUD)))
        kt = {}
        for i, j, reg, nv in ev:
            kt[(i, j)] = kt.get((i, j), KTOP0) - nv
        for bl, lo, hi in BANDS:
            r1 = {s: 0 for s in SECTORS}
            r5 = {s: 0 for s in SECTORS}
            dep = {s: [] for s in SECTORS}
            for i, j, reg, nv in ev:
                r = r_of(i, j)
                if not (lo <= r < hi):
                    continue
                if reg == 1: r1[sec(i, j)] += 1
                elif reg == 5: r5[sec(i, j)] += 1
            for (i, j), k in kt.items():
                r = r_of(i, j)
                if lo <= r < hi:
                    dep[sec(i, j)].append((KTOP0 - k) * DZ)
            print("  %s" % bl)
            print("    sector            " + "".join("%8d" % s for s in SECTORS)
                  + "    max/min")
            for nm, d in (("criterion (reg 1)", r1), ("support   (reg 5)", r5)):
                v = [d[s] for s in SECTORS]
                sp = (max(v) / min(v)) if min(v) > 0 else None
                print("    %-17s " % nm + "".join("%8d" % x for x in v)
                      + ("   %7.2f" % sp if sp else "        --"))
            share = []
            for s in SECTORS:
                t = r1[s] + r5[s]
                share.append(100.0 * r5[s] / t if t else 0.0)
            print("    support share %%   " + "".join("%7.1f%%" % x for x in share))
            meds = [statistics.median(dep[s]) for s in SECTORS if dep[s]]
            if len(meds) == len(SECTORS) and min(meds) > 0:
                print("    median depth mm  " + "".join(
                    "%8.1f" % (1000 * statistics.median(dep[s]))
                    for s in SECTORS)
                    + "   %7.2f" % (max(meds) / min(meds)))
        print()
    print("=" * 78)
    print("READ IT THIS WAY: compare the two legs' `criterion (reg 1)` max/min in\n"
          "the 35-44 and 44-50 mm bands. If P1c's criterion row is as uneven as\n"
          "its depth row, the residual anisotropy is the CRITERION's and P1b was\n"
          "masking it -- P5. If P1c's `support (reg 5)` row is the uneven one,\n"
          "min_sides = 3 is still the wrong definition and the successor is P1d.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
