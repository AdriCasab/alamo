#!/usr/bin/env python3
"""P1e criteria (b) and (c) -- the staircase estimator, OFFLINE. Gates the leg.

  (c) GATE: the factor must be EXACTLY 1.000 on a grid-aligned wall and on a
      thin pillar. Synthetic, exact, no tolerance -- if this fails the
      correction would move axis-aligned physics, and every existing
      axis-aligned unit test with it.
  (b) GATE: on rev 7's 40-50 mm wall, per 15 deg sector, the estimator must sit
      within 0.05 of the angle-averaged analytic 1/(|cos|+|sin|).
      ⚠ The bar is 0.05, not 0.03: BOTH estimators read systematically high near
      the axes (a stencil bias on a staircase, not noise), measured max dev
      0.022, so 0.03 would pass with little margin.

Also reports the floor figures, pooled over nine snapshots, because the D2z
floor-step populations must be RE-READ rather than assumed.

⚠ The open-face test here is the KERNEL's: a face is open when the neighbour's
top is at least ONE cell lower (n = kT - kS >= 1). It is NOT
spall.support_step_min -- that threshold belongs to the support rule and has
nothing to do with how much area the staircase exposes.
"""
import csv, math, os, statistics, sys

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)
DZ = DX = 0.002
KTOP0 = 299
SECTORS = [0, 15, 30, 45, 60, 75]
WALL = (0.040, 0.050)
FLOOR_R = 0.025
SNAPS = (100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0)

LEGS = {
    "rev7 (control, no support rule)":
        os.path.join(STUD, "p0b_one_surface_T", "output_fh_rev7",
                     "FH_removal_events.csv"),
    "P1c (rough)":
        os.path.join(STUD, "p1c_isotropy", "output_fh",
                     "FH_removal_events.csv"),
}


def factor(nxm, nxp, nym, nyp):
    """The shipped formula, mirrored from MMWSpalling.H StaircaseFactor.
    Riser CELL COUNTS -- the ratio is scale-invariant so dz cancels exactly."""
    hx = float(nxm + nxp)
    hy = float(nym + nyp)
    s = hx + hy
    if s <= 0.0:
        return 1.0
    return math.hypot(hx, hy) / s


def analytic_sector(lo_deg, hi_deg, n=2001):
    """Angle-average of 1/(|cos|+|sin|) over the sector."""
    tot = 0.0
    for i in range(n):
        th = math.radians(lo_deg + (hi_deg - lo_deg) * i / (n - 1.0))
        tot += 1.0 / (abs(math.cos(th)) + abs(math.sin(th)))
    return tot / n


def events(path):
    ev = []
    for d in csv.DictReader(open(path)):
        ev.append((float(d["time"]), int(d["col_i"]), int(d["col_j"]),
                   int(float(d["n_voided"]))))
    ev.sort(key=lambda e: e[0])
    return ev


def ktop_at(ev, t):
    k = {}
    for (te, i, j, nv) in ev:
        if te > t:
            break
        k[(i, j)] = k.get((i, j), KTOP0) - nv
    return k


def risers(k, i, j):
    """(n_xm, n_xp, n_ym, n_yp), 0 where the face is not open.
    ⚠ Quarter domain: i = -1 is the MIRROR of cell 0 (prob_lo = 0.0 puts the
    symmetry plane on the FACE of cell 0, so the mirror of -1-i is cell 0, NOT
    cell 1). A mirror has the same k_top, so its riser is 0 and it is never
    open -- which is why the open-face form needs no mirror rule at all."""
    kt = k[(i, j)]
    out = []
    for (ni, nj) in ((i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)):
        if ni < 0 or nj < 0:
            out.append(0); continue
        n = kt - k.get((ni, nj), KTOP0)
        out.append(n if n >= 1 else 0)
    return out


def r_of(i, j):
    return math.hypot((i + 0.5) * DX, (j + 0.5) * DX)


def sec_of(i, j):
    return int(math.degrees(math.atan2(j + 0.5, i + 0.5)) // 15) * 15


def main():
    fails = []

    # ---------------- (c) the exact cases, synthetic ----------------
    print("(c) GATE -- exact canonical cases (no tolerance)")
    cases = (("grid-aligned wall  (one open face, riser h)", (3, 0, 0, 0), 1.0),
             ("thin pillar        (opposing pair, h each)",  (3, 3, 0, 0), 1.0),
             ("thin pillar in y                          ", (0, 0, 3, 3), 1.0),
             ("45 deg corner      (one x, one y, equal)  ", (3, 0, 3, 0),
              1.0 / math.sqrt(2.0)),
             ("no open face                              ", (0, 0, 0, 0), 1.0))
    for lab, rn, want in cases:
        got = factor(*rn)
        ok = (got == want)
        print("  %s  %s  ->  %.8f   want %.8f   %s"
              % ("ok  " if ok else "FAIL", lab, got, want,
                 "" if ok else "<-- EXACT equality required"))
        if not ok:
            fails.append("(c) %s gave %.8f, want %.8f" % (lab.strip(), got, want))
    # range invariant: the formula can never leave [1/sqrt2, 1]
    lo, hi = 1.0, 0.0
    for a in range(0, 9):
        for b in range(0, 9):
            for c in range(0, 9):
                for d in range(0, 9):
                    if a + b + c + d == 0:
                        continue
                    f = factor(a, b, c, d)
                    lo = min(lo, f); hi = max(hi, f)
    print("  ok    range over all riser combinations 0..8: [%.8f, %.8f] "
          "(must be within [%.8f, 1])" % (lo, hi, 1 / math.sqrt(2)))
    if lo < 1.0 / math.sqrt(2.0) - 1e-15 or hi > 1.0 + 1e-15:
        fails.append("(c) the factor left [1/sqrt2, 1] -- no clamp should be needed")

    # ---------------- (b) rev 7's wall, per sector ----------------
    for lab, path in LEGS.items():
        if not os.path.exists(path):
            print("\nMISSING %s" % path); fails.append("missing %s" % lab); continue
        ev = events(path)
        k = ktop_at(ev, ev[-1][0])
        print("\n%s   (final state, t = %.1f s)" % (lab, ev[-1][0]))
        print("  %-22s" % "40-50 mm wall" + "".join("%9d" % s for s in SECTORS)
              + "     mean   max dev")
        per = {s: [] for s in SECTORS}
        for (i, j) in k:
            if not (WALL[0] <= r_of(i, j) < WALL[1]):
                continue
            rn = risers(k, i, j)
            if sum(rn) == 0:
                continue
            per[sec_of(i, j)].append(factor(*rn))
        est = [statistics.mean(per[s]) if per[s] else float("nan")
               for s in SECTORS]
        ana = [analytic_sector(s, s + 15) for s in SECTORS]
        dev = [abs(e - a) for e, a in zip(est, ana) if e == e]
        print("  %-22s" % "open-face estimator" + "".join("%9.3f" % e for e in est)
              + "   %6.3f   %7.3f" % (statistics.mean([e for e in est if e == e]),
                                      max(dev) if dev else float("nan")))
        print("  %-22s" % "analytic 1/(|c|+|s|)" + "".join("%9.3f" % a for a in ana)
              + "   %6.3f" % statistics.mean(ana))
        print("  %-22s" % "n columns" + "".join("%9d" % len(per[s]) for s in SECTORS))
        if lab.startswith("rev7"):
            if not dev or max(dev) > 0.05:
                fails.append("(b) rev 7 wall max deviation from analytic %.3f "
                             "> 0.05" % (max(dev) if dev else float("nan")))
            else:
                print("  ok    (b) max deviation %.3f <= 0.05" % max(dev))

    # ---------------- the floor, pooled over nine snapshots ----------------
    path = LEGS["rev7 (control, no support rule)"]
    if os.path.exists(path):
        ev = events(path)
        print("\nrev 7 FLOOR (r < %.0f mm), pooled over %d snapshots t = %.0f..%.0f s"
              % (1000 * FLOOR_R, len(SNAPS), SNAPS[0], SNAPS[-1]))
        print("  ⚠ The D2z floor-step populations must be RE-READ, not assumed. "
              "Axis-aligned\n    unit tests are unaffected (factor exactly 1).")
        vals, nfaces = [], {1: 0, 2: 0, 3: 0, 4: 0}
        for t in SNAPS:
            k = ktop_at(ev, t)
            for (i, j) in k:
                if r_of(i, j) >= FLOOR_R:
                    continue
                rn = risers(k, i, j)
                nf = sum(1 for x in rn if x > 0)
                if nf == 0:
                    continue
                nfaces[nf] += 1
                vals.append(factor(*rn))
        if vals:
            tot = len(vals)
            print("  mean %.3f   median %.3f   n = %d"
                  % (statistics.mean(vals), statistics.median(vals), tot))
            print("  open-face histogram: " + "  ".join(
                "%d-face %d (%.0f %%)" % (nf, nfaces[nf], 100.0 * nfaces[nf] / tot)
                for nf in (1, 2, 3, 4) if nfaces[nf]))
            print("  ⚠ the median is 1.000 when most floor steps have a SINGLE "
                  "open face;\n    a two-face corner tower loses 1 - 1/sqrt(2) "
                  "= 29.3 %.")

    print()
    if fails:
        for f in fails: print("FAIL: %s" % f)
        print("\n⚠⚠ A GATE FAILED -- the estimator is wrong. Launch no leg: it "
              "would be uninterpretable.")
        return 1
    print("(b) and (c) PASS -- the estimator is exact on the canonical cases and "
          "within 0.05 of\nanalytic on rev 7's wall. The leg is cleared.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
