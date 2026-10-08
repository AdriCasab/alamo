#!/usr/bin/env python3
"""P1c criteria (f), (j), (k), (l) -- roundness, offline. Starts no run.

⚠⚠ (f) IS HOLE RADIUS AT FIXED DEPTH, PER 15 deg SECTOR -- not depth inside a
fixed radial ring. Depth-in-a-ring reads max/min ~ 1 whenever the ring is
uniformly SHALLOW, so it can pass while saying nothing; radius-at-fixed-depth
cannot. The depth-in-a-ring form survives only as reported (j).

⚠ Per-sector statistics throughout, never a global median -- a global median is
how P1b's criterion (f) passed at 2.0 mm with an 18x spread just outside its cut.

Usage:
    roundness.py --selftest          reproduce the frozen reference values
    roundness.py <event-log.csv>     score one leg
"""
import csv, math, os, statistics, sys

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)
DZ = DX = 0.002
KTOP0 = 299                      # the original surface cell index
DEPTHS = (0.150, 0.200, 0.250)   # (f): the calibrated section, never 100 mm
SECTOR = 15.0                    # degrees
RING_J = (0.044, 0.050)          # (j)
RING_K = (0.035, 0.044)          # (k): where P5 says the hole narrows
BORE_L = (0.040, 0.100, 0.150, 0.200, 0.250)   # (l), report only

REF = {
    "rev7": (os.path.join(STUD, "p0b_one_surface_T", "output_fh_rev7",
                          "FH_removal_events.csv"), (1.015, 1.015, 1.027)),
    "P1b":  (os.path.join(STUD, "p1b_support_rule", "output_fh",
                          "FH_removal_events.csv"), (1.095, 1.095, 1.108)),
}


def load(path):
    """-> (k_top_final per column, per-column support/total removal counts).

    ⚠ `k_top` in the event log is the top BEFORE that removal; the top after it
    is k_top - n_voided. The running state is therefore k_top0 minus the
    cumulative n_voided, which is what this returns.
    """
    kt, tot, sup, tlast = {}, {}, {}, {}
    for d in csv.DictReader(open(path)):
        c = (int(d["col_i"]), int(d["col_j"]))
        nv = int(float(d["n_voided"]))
        reg = int(float(d["regime"]))
        t = float(d["time"])
        kt[c] = kt.get(c, KTOP0) - nv
        tot[c] = tot.get(c, 0) + 1
        if reg == 5:
            sup[c] = sup.get(c, 0) + 1
        tlast[c] = t
    return kt, tot, sup, tlast


def r_of(i, j):
    return math.hypot((i + 0.5) * DX, (j + 0.5) * DX)


def sector_of(i, j):
    return int(math.degrees(math.atan2(j + 0.5, i + 0.5)) // SECTOR) * int(SECTOR)


SECTORS = [0, 15, 30, 45, 60, 75]


def roundness(kt, depth):
    """(f): max radius per sector among columns drilled at least `depth` deep."""
    k_cut = KTOP0 - int(round(depth / DZ))
    per = {s: [] for s in SECTORS}
    for (i, j), k in kt.items():
        if k <= k_cut:
            per[sector_of(i, j)].append(r_of(i, j))
    radii = {s: (max(v) if v else None) for s, v in per.items()}
    counts = {s: len(v) for s, v in per.items()}
    vals = [v for v in radii.values() if v]
    scorable = len(vals) == len(SECTORS)
    score = (max(vals) / min(vals)) if scorable and min(vals) > 0 else None
    return radii, counts, score, scorable


def report(path, label=None):
    kt, tot, sup, tlast = load(path)
    print("%s\n  %d columns with removal" % (label or path, len(kt)))

    print("\n(f) ROUNDNESS -- max hole radius per %.0f deg sector [mm], "
          "bar max/min <= 1.05" % SECTOR)
    scores = []
    for d in DEPTHS:
        radii, counts, score, ok = roundness(kt, d)
        row = "  depth %3.0f mm: " % (1000 * d) + " ".join(
            ("%5.1f" % (1000 * radii[s])) if radii[s] else "   --"
            for s in SECTORS)
        if ok:
            row += "   max/min %6.3f  %s" % (score, "PASS" if score <= 1.05
                                             else "FAIL")
        else:
            empty = [s for s in SECTORS if not radii[s]]
            row += "   ⚠ NOT SCORABLE: sectors %s hold no column at this depth" \
                   % empty
        print(row)
        print("      per-sector column counts: " + " ".join(
            "%5d" % counts[s] for s in SECTORS))
        scores.append(score)

    print("\n(j) depth per sector in the %.0f-%.0f mm ring [mm] -- REPORTED only"
          % (1000 * RING_J[0], 1000 * RING_J[1]))
    _ring(kt, tot, sup, RING_J)
    print("\n(k) the %.0f-%.0f mm ring -- depth, removals, support share"
          % (1000 * RING_K[0], 1000 * RING_K[1]))
    _ring(kt, tot, sup, RING_K)

    print("\n(l) bore radius per sector at fixed depth [mm] -- REPORTED only")
    print("    ⚠ never quoted as a single diameter: a max over all azimuths "
          "picks the diagonal.\n    ⚠ 100 mm is the pause underream (Meier "
          "p. 223), reported but never scored.")
    for d in BORE_L:
        radii, counts, score, ok = roundness(kt, d)
        print("    %3.0f mm: " % (1000 * d) + " ".join(
            ("%5.1f" % (1000 * radii[s])) if radii[s] else "   --"
            for s in SECTORS)
            + ("   max/min %6.3f" % score if ok else "   not scorable"))
    return scores


def _ring(kt, tot, sup, ring):
    per = {s: [] for s in SECTORS}
    for (i, j), k in kt.items():
        r = r_of(i, j)
        if ring[0] <= r < ring[1]:
            per[sector_of(i, j)].append(((KTOP0 - k) * DZ, tot.get((i, j), 0),
                                         sup.get((i, j), 0)))
    meds = []
    print("    sector      n   median depth   removals   support   share")
    for s in SECTORS:
        v = per[s]
        if not v:
            print("    %2d-%2d     --" % (s, s + 15)); continue
        md = statistics.median([a for a, _, _ in v])
        nt = sum(b for _, b, _ in v)
        ns = sum(c for _, _, c in v)
        meds.append(md)
        print("    %2d-%2d  %5d   %8.1f mm   %8d  %8d  %5.1f %%"
              % (s, s + 15, len(v), 1000 * md, nt, ns,
                 100.0 * ns / max(1, nt)))
    if meds and min(meds) > 0:
        print("    per-sector median depth spread: %.2fx" % (max(meds) / min(meds)))


def selftest():
    ok = True
    for lab, (path, want) in REF.items():
        if not os.path.exists(path):
            print("MISSING %s: %s" % (lab, path)); ok = False; continue
        print("=" * 72)
        got = report(path, "SELFTEST %s" % lab)
        print("\n  selftest %s: (f) got %s, want %s"
              % (lab, ["%.3f" % g if g else None for g in got],
                 list(want)))
        for g, w in zip(got, want):
            if g is None or abs(g - w) > 0.0015:
                print("  FAIL: %s does not reproduce %.3f (got %s)"
                      % (lab, w, g)); ok = False
        print()
    print("=" * 72)
    print("SELFTEST %s -- a statistic that cannot reproduce P1b's 1.095/1.095/"
          "1.108\nis blind to the defect and must fail here, before any leg runs."
          % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        sys.exit(selftest())
    p = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        HERE, "output_fh", "FH_removal_events.csv")
    if not os.path.exists(p):
        print("missing: %s" % p); sys.exit(1)
    report(p)
