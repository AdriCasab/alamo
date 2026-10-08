#!/usr/bin/env python3
"""P1c criteria (b) and (c) -- the REAL firing condition, offline. Gates the leg.

Successor to P1b's exposure.py, which had two planner defects this closes:

  1. it tested `n_exposed >= min_sides` and OMITTED the `|| pillar` branch, so
     it could not see that a 2-opposing-side column fires at ANY min_sides
     (Removal.H:863, :1511);
  2. it scoped everything to r < 40 mm (R_BURNER), so the wall band where the
     grid-diagonal artefact lives was outside the cut entirely.

Three bands, per 15 deg sector, both frozen geometries, min_sides in {2, 3}.

  (b) GATE: rev 7's HEALTHY geometry must qualify 0.0 % on EVERY band and EVERY
      sector at min_sides = 3, with 0 pillar-only. If not, 3-of-4 is not
      grid-isotropic either and the long leg would be uninterpretable.
  (c) GATE: P1's STALLED floor must still qualify >= 20 % at min_sides = 3 with
      a sector max/min <= 2.0, or the reach that cured the stall is gone.
"""
import csv, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)
DZ = DX = 0.002
KTOP0 = 299
STEP_MIN = 0.006          # spall.support_step_min, NOT tuned here
SECTORS = [0, 15, 30, 45, 60, 75]
BANDS = (("floor  r < 40", 0.000, 0.040),
         ("wall   40-50 ", 0.040, 0.050),
         ("mouth  50-60 ", 0.050, 0.060))

LEGS = {
    "rev7 (HEALTHY)": os.path.join(STUD, "p0b_one_surface_T",
                                   "output_fh_rev7", "FH_removal_events.csv"),
    "P1 (STALLED)":   os.path.join(STUD, "p1_gradient_criterion",
                                   "output_fh", "FH_removal_events.csv"),
}


def ktop_final(path):
    """k_top after the last removal. ⚠ the log's k_top is the top BEFORE the
    removal, so the running top is k_top0 minus the cumulative n_voided. The
    field is cross-checked against that running state and the mismatch count is
    PRINTED (P1b's exposure.py computed it and silently threw it away)."""
    k, bad = {}, 0
    for d in csv.DictReader(open(path)):
        c = (int(d["col_i"]), int(d["col_j"]))
        cur = k.get(c, KTOP0)
        if cur != int(float(d["k_top"])):
            bad += 1
        k[c] = cur - int(float(d["n_voided"]))
    return k, bad


def r_of(i, j):
    return math.hypot((i + 0.5) * DX, (j + 0.5) * DX)


def sector_of(i, j):
    return int(math.degrees(math.atan2(j + 0.5, i + 0.5)) // 15) * 15


def fires(k, i, j, min_sides):
    """Removal.H:367-399 + :863 exactly: the four face neighbours, the one-cell
    floor, the length threshold, the domain-edge skip, then
    (n_exposed >= min_sides) || pillar."""
    kt = k[(i, j)]
    nb = ((i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1))
    open_ = [False] * 4
    for d, (ni, nj) in enumerate(nb):
        if ni < 0 or nj < 0:        # quarter-domain mirror: step is 0, never open
            continue
        kn = k.get((ni, nj), KTOP0)
        ncell = kt - kn
        if ncell < 1:
            continue
        if ncell * DZ < STEP_MIN - 1e-12:
            continue
        open_[d] = True
    n = sum(open_)
    pillar = (open_[0] and open_[1]) or (open_[2] and open_[3])
    return (n >= min_sides or pillar), n, pillar


def main():
    fails = []
    for lab, path in LEGS.items():
        if not os.path.exists(path):
            print("MISSING %s" % path); return 1
        k, bad = ktop_final(path)
        print("=" * 78)
        print("%s   %s" % (lab, os.path.relpath(path, STUD)))
        print("  columns with removal %d; k_top field vs replayed state: "
              "%d mismatches" % (len(k), bad))
        for ms in (2, 3):
            print("  min_sides = %d" % ms)
            print("    %-14s   overall        " % "band"
                  + "".join("%7d" % s for s in SECTORS) + "   pillar-only")
            for bl, lo, hi in BANDS:
                cols = [c for c in k if lo <= r_of(*c) < hi]
                if not cols:
                    print("    %-14s   (no columns)" % bl); continue
                hit = po = 0
                per = {s: [0, 0] for s in SECTORS}
                for c in cols:
                    f, n, pil = fires(k, c[0], c[1], ms)
                    s = sector_of(*c)
                    per[s][1] += 1
                    if f:
                        hit += 1
                        per[s][0] += 1
                        if n < ms and pil:
                            po += 1
                row = "    %-14s  %5.1f %% (%4d/%4d)" % (
                    bl, 100.0 * hit / len(cols), hit, len(cols))
                row += "".join(
                    "%6.0f%%" % (100.0 * per[s][0] / max(1, per[s][1]))
                    for s in SECTORS)
                row += "   %5d" % po
                print(row)
                # ---- the gates ----
                if lab.startswith("rev7") and ms == 3:
                    if hit != 0 or po != 0:
                        fails.append("(b) rev 7 at min_sides=3, %s: %d qualify "
                                     "(%d pillar-only) -- must be 0" % (bl, hit, po))
                if lab.startswith("P1 (") and ms == 3 and bl.startswith("floor"):
                    frac = 100.0 * hit / len(cols)
                    secs = [100.0 * per[s][0] / max(1, per[s][1])
                            for s in SECTORS]
                    sp = max(secs) / max(1e-9, min(secs))
                    print("      (c) reach %.1f %% (bar >= 20) ; sector max/min "
                          "%.2f (bar <= 2.0)" % (frac, sp))
                    if frac < 20.0:
                        fails.append("(c) P1 stalled floor reach %.1f %% < 20 %%"
                                     % frac)
                    if sp > 2.0:
                        fails.append("(c) P1 stalled floor sector max/min %.2f "
                                     "> 2.0 -- the reach is itself anisotropic" % sp)
        print()
    print("=" * 78)
    if fails:
        for f in fails: print("FAIL: %s" % f)
        print("\n⚠⚠ A GATE FAILED. Launch no long leg: report this row as the "
              "result.")
        return 1
    print("(b) PASS -- rev 7's healthy geometry qualifies 0.0 % on every band "
          "and sector at\n    min_sides = 3, with 0 pillar-only.")
    print("(c) PASS -- P1's stalled floor keeps enough isotropic reach.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
