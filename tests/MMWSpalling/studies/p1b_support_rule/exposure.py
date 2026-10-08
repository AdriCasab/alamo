#!/usr/bin/env python3
"""P1b criterion (k) + witness 1 -- the exposure test applied OFFLINE.

Starts no run. Replays the per-column top-cell index k_top(t) from a shipped
removal-event log and applies the packet's exposure rule

    n_exposed = # of 4-neighbours with (k_top[self] - k_top[nbr])*dz >= step_min

to BOTH geometries: rev 7's healthy floor (drilled to 397 mm, never stalled)
and P1's stalled one. The rule is a repair only if it is INERT on the first.

⚠ This must pass before the Meier leg is launched (packet, Commands block).
"""
import csv, math, os, sys

HERE  = os.path.dirname(os.path.abspath(__file__))
STUD  = os.path.dirname(HERE)
DX = DZ = 0.002
KTOP0 = 299
R_BURNER = 0.040

LEGS = {
    "rev7 (healthy)": os.path.join(STUD, "p0b_one_surface_T",
                                   "output_fh_rev7", "FH_removal_events.csv"),
    "P1 (stalled)":   os.path.join(STUD, "p1_gradient_criterion",
                                   "output_fh", "FH_removal_events.csv"),
}
TIMES = (100.0, 150.0, 178.0)


def replay(path):
    """-> list of (time, col_i, col_j, n_voided), sorted by time."""
    ev = []
    for d in csv.DictReader(open(path)):
        nv = int(float(d["n_voided"]))
        if nv <= 0:
            continue
        ev.append((float(d["time"]), int(d["col_i"]), int(d["col_j"]), nv,
                   int(float(d["k_top"]))))
    ev.sort(key=lambda e: e[0])
    return ev


def ktop_at(ev, t):
    """k_top per column at time t. Cross-checked against the log's own k_top
    field: the event's k_top must equal the running state before the cut."""
    k, bad = {}, 0
    for (te, i, j, nv, klog) in ev:
        if te > t:
            break
        cur = k.get((i, j), KTOP0)
        if cur != klog:
            bad += 1
        k[(i, j)] = cur - nv
    return k, bad


def r_of(i, j):
    return math.hypot((i + 0.5) * DX, (j + 0.5) * DX)


def exposure(k, step_min):
    """-> dict col -> n_exposed, over the columns inside the burner radius."""
    out = {}
    for (i, j), kt in k.items():
        if r_of(i, j) >= R_BURNER:
            continue
        n = 0
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            ni, nj = i + di, j + dj
            # Quarter domain: i = -1 is the MIRROR of i = 0, not a missing
            # neighbour. The mirror has the same k_top, so it is never exposed
            # -- identical to the kernel's domain-edge `continue`.
            if ni < 0 or nj < 0:
                continue
            kn = k.get((ni, nj), KTOP0)
            if (kt - kn) * DZ >= step_min - 1e-12:
                n += 1
        out[(i, j)] = n
    return out


def main():
    print("P1b (k): fraction of columns at r < %.0f mm with >= 2 exposed sides"
          % (1000 * R_BURNER))
    print("denominator = columns inside the burner radius that have fired at "
          "least once by t\n")
    for lab, path in LEGS.items():
        if not os.path.exists(path):
            print("%-16s MISSING: %s" % (lab, path)); continue
        ev = replay(path)
        print("%s   %s" % (lab, os.path.relpath(path, STUD)))
        print("    events %d, t_max %.3f s" % (len(ev), ev[-1][0]))
        # P1c deliverable 6: `bad` was COMPUTED and silently thrown away. It is
        # the cross-check that the replayed k_top agrees with the log's own
        # k_top field, i.e. that this reduction is reading the geometry the run
        # actually had. A cross-check that is never printed is not a check.
        _, _bad = ktop_at(ev, ev[-1][0])
        print("    k_top field vs replayed state: %d mismatches%s"
              % (_bad, "" if _bad == 0 else "   ⚠⚠ THE REPLAY DISAGREES WITH "
                 "THE LOG -- every number below is suspect"))
        # (g) BOTH axes: min_sides in {2, 3} x step_min in {4, 6, 8} mm, plus
        # 2 mm for the record. ⚠ The sweep is a sensitivity statement, not a
        # search for a pass: (k) decides which settings are shippable, and it
        # costs no compute, so no setting needs a Meier run to be excluded.
        for ms in (2, 3):
            print("    min_sides = %d" % ms)
            print("    %-18s" % "step threshold" + "".join(
                "  t=%-11.0f" % t for t in TIMES))
            for cells in (1, 2, 3, 4):
                sm = cells * DZ
                row = []
                for t in TIMES:
                    k, bad = ktop_at(ev, t)
                    exp = exposure(k, sm)
                    tot = len(exp)
                    hit = sum(1 for v in exp.values() if v >= ms)
                    row.append("%5.1f %% (%d/%d)"
                               % (100.0 * hit / max(1, tot), hit, tot))
                print("    %-18s" % ("%d cell = %.0f mm" % (cells, 1000 * sm))
                      + "".join("  %-14s" % r for r in row))
        print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
