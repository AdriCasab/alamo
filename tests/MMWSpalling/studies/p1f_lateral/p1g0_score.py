#!/usr/bin/env python3
"""P1g-0 ord-2 scoring: (w1), (t1), (d1), (x-vs-pre) and the (z) verdict.

Scores the FIXED Leg B (`output_legB_fix`) against the PRE-FIX one
(`output_legB`, P1f revision 11) and the control (`output_rev7sc`).

⚠⚠ Why a separate file from `analyze.py`: that file carries P1f's ACCEPTED
scoring and its paths are cited by P1f's results and review findings. The rows
below are comparisons between two legs, which it was never written to do.

⚠⚠ THE REFUSAL RULE, stated once (it belongs to (d1) alone): a band is refused
only when BOTH sides have fewer than 5 firings. P1f's rule -- refuse if EITHER
side is under 5 -- would have suppressed this step's clearest result, because
`0 against 185` IS the finding and not an under-sampled sector. Where one side
is zero the row prints the counts and a ONE-SIDED bound.

⚠ Every ratio is printed with the magnitude it is a ratio of, and every count is
divided by the grid census of the same sector.

⚠ No decay-length read anywhere: `lat_ell_med` is retired and a two-cell log
ratio cannot resolve ell below ~1 cell. Firing RATE and temperature TAIL only.
"""
import collections, csv, math, os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np
import yt
yt.set_log_level(50)

HERE = os.path.dirname(os.path.abspath(__file__))
PRE = os.path.join(HERE, "output_legB")          # P1f rev 11, pre-fix
FIX = os.path.join(HERE, "output_legB_fix")      # P1g-0, the fix
CTRL = os.path.join(HERE, "output_rev7sc")
OUTF = os.path.join(HERE, "P1G0_SCORED.txt")

# the SAME sector definition as azimuth_witness.py -- printed, never implied
SECT6 = [(0, 15), (15, 30), (30, 45), (45, 60), (60, 75), (75, 90)]
AXIS_S = (0, 5)
DIAG_S = (2, 3)
RBANDS = ((40, 60), (60, 90))
DBANDS = ((0, 100), (100, 200), (200, 300), (300, 400))
SQRT2_PCT = 41.42      # the per-column figure this criterion is meant to beat
KTOP0, DZ, DX = 299, 0.002, 0.002


def ev(base):
    p = os.path.join(base, "FH_removal_events.csv")
    return list(csv.DictReader(open(p))) if os.path.exists(p) else []


def thermo(base, key):
    p = os.path.join(base, "FH", "thermo.dat")
    if not os.path.exists(p):
        return [], []
    h = open(p).readline().rstrip("\n").split("\t")
    if key not in h or "time" not in h:
        return [], []
    it, ik = h.index("time"), h.index(key)
    t, v = [], []
    for ln in open(p).read().split("\n")[1:]:
        f = ln.split("\t")
        if len(f) > max(it, ik):
            t.append(float(f[it])); v.append(float(f[ik]))
    return t, v


def plotfiles(base):
    d = os.path.join(base, "FH")
    if not os.path.isdir(d):
        return []
    return sorted(os.path.join(d, x) for x in os.listdir(d)
                  if x.endswith("cell") and os.path.isdir(os.path.join(d, x)))


def sector_of(i, j):
    th = math.degrees(math.atan2(j + 0.5, i + 0.5))
    for n, (a, b) in enumerate(SECT6):
        if a <= th < b or (b == 90 and th >= 90.0 - 1e-12):
            return n
    return None


def band_of(v, bands):
    for n, (a, b) in enumerate(bands):
        if a <= v < b:
            return n
    return None


def census():
    """columns per (radial band, sector) on the 60x60 quarter grid."""
    c = collections.Counter()
    for i in range(60):
        for j in range(60):
            s = sector_of(i, j)
            if s is None:
                continue
            r = DX * math.hypot(i + 0.5, j + 0.5) * 1e3
            rb = band_of(r, RBANDS)
            if rb is not None:
                c[(rb, s)] += 1
    return c


def tally(base, regime):
    """(radial band, depth band, sector) -> list of T_top at firing."""
    acc = collections.defaultdict(list)
    for r in ev(base):
        try:
            if int(float(r["regime"])) != regime:
                continue
            i, j = int(r["col_i"]), int(r["col_j"])
            kt = int(float(r["k_top"]))
            rd = float(r.get("riser_depth") or 0.0)
            T = float(r["T_top"])
        except (KeyError, ValueError):
            continue
        s = sector_of(i, j)
        if s is None:
            continue
        rr = DX * math.hypot(i + 0.5, j + 0.5) * 1e3
        rb = band_of(rr, RBANDS)
        # ⚠ k_top is the top BEFORE that removal (standing rule 4)
        depth = (KTOP0 - kt) * DZ * 1e3 + rd * 1e3
        db = band_of(depth, DBANDS)
        if rb is None or db is None:
            continue
        acc[(rb, db, s)].append(T)
    return acc


def main():
    out = []

    def P(*a):
        s = " ".join(str(x) for x in a)
        print(s)
        out.append(s)

    for tag, b in (("pre-fix", PRE), ("fixed", FIX), ("control", CTRL)):
        if not os.path.isdir(os.path.join(b, "FH")):
            P("FAIL: %s leg missing at %s" % (tag, b))
            return 1

    cen = census()
    P("=" * 78)
    P("P1g-0 ORD-2 SCORING -- the FIX against the PRE-FIX leg and the control")
    P("=" * 78)
    P("  fixed   : %s" % FIX)
    P("  pre-fix : %s" % PRE)
    P("  control : %s" % CTRL)
    P("  sectors : theta = atan2(y,x) deg, quarter domain, six 15 deg bins "
      "%s" % "  ".join("%d-%d" % s for s in SECT6))
    P("            AXIS = 0-15 and 75-90;  DIAGONAL = 30-45 and 45-60")
    P("  census  : r 40-60 mm per sector %s"
      % "  ".join("%d" % cen[(0, s)] for s in range(6)))
    P("            r 60-90 mm per sector %s"
      % "  ".join("%d" % cen[(1, s)] for s in range(6)))

    acc6 = {"pre": tally(PRE, 6), "fix": tally(FIX, 6)}
    acc1 = {"pre": tally(PRE, 1), "fix": tally(FIX, 1)}

    # ---------------- (e) completion and the ledger ------------------------
    P("")
    P("-" * 78)
    P("(e) COMPLETION AND THE ENERGY INVARIANT  [gate: 500.0 s rc 0, <= 5e-3]")
    for tag, b in (("fixed  ", FIX), ("pre-fix", PRE), ("control", CTRL)):
        log = os.path.join(b, "FH.log")
        txt = open(log).read() if os.path.exists(log) else ""
        tl = [l for l in txt.split("\n") if "TIME =" in l]
        t, e = thermo(b, "ledger_err_solid")
        late = [abs(x) for tt, x in zip(t, e) if tt > 50.0]
        P("  %s last step: %s" % (tag, (tl[-1].strip() if tl else "(none)")))
        P("  %s finalized %d  aborts/collisions %d  max|err| t>50s %.3e -> %s"
          % (tag, txt.count("finalized"),
             sum(txt.lower().count(k) for k in ("abort", "collision",
                                                "ran out of rock")),
             max(late) if late else float("nan"),
             "PASS" if late and max(late) <= 5e-3 else "FAIL"))

    # ---------------- (x) volume, fixed vs pre-fix vs control --------------
    P("")
    P("-" * 78)
    P("(x) REMOVED VOLUME, by regime  [report: vs the control AND vs pre-fix]")
    vol = {}
    for tag, b in (("control", CTRL), ("pre-fix", PRE), ("fixed", FIX)):
        rows = ev(b)
        per = collections.Counter()
        cells = 0
        for r in rows:
            try:
                n = int(float(r["n_voided"])); g = int(float(r["regime"]))
            except (KeyError, ValueError):
                continue
            per[g] += n
            cells += n
        vol[tag] = (len(rows), cells, per)
        P("  %-8s %7d events  %7d cells   %s"
          % (tag, len(rows), cells,
             "  ".join("r%d:%d" % (g, per[g]) for g in sorted(per))))
    if vol["control"][1]:
        P("  fixed / control   = %.3f x      pre-fix / control = %.3f x"
          % (vol["fixed"][1] / vol["control"][1],
             vol["pre-fix"][1] / vol["control"][1]))
    if vol["pre-fix"][2][6]:
        P("  regime 6: pre-fix %d -> fixed %d   (%+.1f %%)"
          % (vol["pre-fix"][2][6], vol["fixed"][2][6],
             100.0 * (vol["fixed"][2][6] / vol["pre-fix"][2][6] - 1.0)))

    # ---------------- (w1) axis vs diagonal per-column rate ---------------
    P("")
    P("-" * 78)
    P("(w1) PER-COLUMN FIRING RATE, axis vs diagonal, fixed vs pre-fix")
    P("     ⚠⚠ a REPORT, not a gate: 'axis counts unchanged' CANNOT pass on a")
    P("        coupled run -- the first diagonal firing changes the geometry,")
    P("        hence the azimuthally binned jet march and every later column.")
    P("        The exact no-op claim lives in the unit test's (u1).")
    P("     two-sided expectation: |axis move| < diagonal GAIN")
    for rb, (r0, r1) in enumerate(RBANDS):
        for reg, acc in ((6, acc6), (1, acc1)):
            ca = sum(cen[(rb, s)] for s in AXIS_S)
            cd = sum(cen[(rb, s)] for s in DIAG_S)
            if not (ca and cd):
                continue
            row = []
            for leg in ("pre", "fix"):
                na = sum(len(acc[leg][(rb, db, s)])
                         for db in range(len(DBANDS)) for s in AXIS_S)
                nd = sum(len(acc[leg][(rb, db, s)])
                         for db in range(len(DBANDS)) for s in DIAG_S)
                row.append((na, nd, na / ca, nd / cd))
            (na0, nd0, pa0, pd0), (na1, nd1, pa1, pd1) = row
            dax = (pa1 / pa0 - 1.0) * 100.0 if pa0 else float("nan")
            ddg = (pd1 / pd0 - 1.0) * 100.0 if pd0 else float("nan")
            P("     r %d-%d reg %d:  axis %d->%d (%.4f->%.4f /col, %+.2f %%)  "
              "diag %d->%d (%.4f->%.4f /col, %+.2f %%)"
              % (r0, r1, reg, na0, na1, pa0, pa1, dax,
                 nd0, nd1, pd0, pd1, ddg))
            if pa0 and pd0:
                ok = abs(dax) < ddg
                P("       |axis move| %.2f %%  vs  diagonal gain %+.2f %%  ->  %s"
                  % (abs(dax), ddg, "as expected" if ok else "⚠ NOT as expected"))

    # ---------------- (t1) the MECHANISM: firing temperature spread -------
    P("")
    P("-" * 78)
    P("(t1) REGIME-6 CELL TEMPERATURE AT FIRING, spread across the six sectors")
    P("     [report; PREDICTION <= 5 K. This tests the fix's MECHANISM: if the")
    P("      spread stays wide, the conductance was not the binding term.]")
    for leg, nm in (("pre", "pre-fix"), ("fix", "fixed")):
        for rb, (r0, r1) in enumerate(RBANDS):
            per = []
            for s in range(6):
                v = [x for db in range(len(DBANDS))
                     for x in acc6[leg][(rb, db, s)]]
                per.append((len(v), float(np.median(v)) if v else 0.0))
            got = [m for n, m in per if n >= 5]
            P("     %-8s r %d-%d:  %s"
              % (nm, r0, r1,
                 "  ".join("%d-%d:n=%d,T=%.1f" % (SECT6[s][0], SECT6[s][1],
                                                  per[s][0], per[s][1])
                           for s in range(6))))
            if len(got) >= 2:
                P("              spread over sectors with n>=5: %.1f K  "
                  "[predict <= 5 K]  %s"
                  % (max(got) - min(got),
                     "PASS" if (max(got) - min(got)) <= 5.0 else "MISS"))
            else:
                P("              fewer than two sectors with n>=5 -- no spread")

    # ---------------- (d1) the deficit by depth band ----------------------
    P("")
    P("-" * 78)
    P("(d1) LATERAL DEFICIT BY DEPTH BAND -- deficit AND deficit/%.2f" % SQRT2_PCT)
    P("     ⚠⚠ P1f printed %.2f/deficit, so a deficit 1.54x LARGER than the"
      % SQRT2_PCT)
    P("        per-column figure read '0.6x' -- INVERTED. Here it is")
    P("        deficit/%.2f, so > 1 means WORSE than per-column." % SQRT2_PCT)
    P("     ⚠⚠ a band is refused ONLY when BOTH sides are under 5 firings.")
    deficits = {}
    for rb, (r0, r1) in enumerate(RBANDS):
        for leg, nm in (("pre", "pre-fix"), ("fix", "fixed")):
            P("")
            P("     %-8s r = %d-%d mm" % (nm, r0, r1))
            P("     %-10s %s %9s %9s %9s %9s"
              % ("depth[mm]",
                 "".join("%8s" % ("%d-%d" % s) for s in SECT6),
                 "axis/col", "diag/col", "deficit", "/%.0f" % SQRT2_PCT))
            for db, (d0, d1) in enumerate(DBANDS):
                cnt = [len(acc6[leg][(rb, db, s)]) for s in range(6)]
                na = sum(cnt[s] for s in AXIS_S)
                nd = sum(cnt[s] for s in DIAG_S)
                ca = sum(cen[(rb, s)] for s in AXIS_S)
                cd = sum(cen[(rb, s)] for s in DIAG_S)
                pa = na / ca if ca else 0.0
                pd = nd / cd if cd else 0.0
                if na < 5 and nd < 5:
                    P("     %-10s %s   -- refused (axis %d, diag %d: BOTH "
                      "under 5) --"
                      % ("%d-%d" % (d0, d1),
                         "".join("%8d" % v for v in cnt), na, nd))
                    continue
                if pa > 0:
                    dfc = 100.0 * (1.0 - pd / pa)
                    tag = "%+8.2f" % dfc
                    rat = "%9.3f" % (dfc / SQRT2_PCT)
                    deficits[(leg, rb, db)] = dfc
                else:
                    dfc = float("nan"); tag = "    n/a"; rat = "      n/a"
                note = ""
                if nd == 0:
                    note = "  [one-sided: diag is 0 of %d cols]" % cd
                elif na == 0:
                    note = "  [one-sided: axis is 0 of %d cols]" % ca
                P("     %-10s %s %9.4f %9.4f %s %s%s"
                  % ("%d-%d" % (d0, d1),
                     "".join("%8d" % v for v in cnt), pa, pd, tag, rat, note))

    # ---------------- (z) the verdict -------------------------------------
    P("")
    P("=" * 78)
    P("(z) VERDICT -- does the deep deficit fall materially at 100-400 mm?")
    P("    direction only, against (n2)'s pre-measured residual.")
    P("    ⚠ P3 pre-registered: the deficit FALLS, magnitude NOT predicted; a")
    P("      residual SURVIVES (the clipped vertical tail, (n2)). A drop from")
    P("      100 % to ~40 % is the fix WORKING, not a partial failure.")
    P("    ⚠ NO prediction was registered for 0-100 mm (P4) -- printed, not scored.")
    moved = []
    for rb, (r0, r1) in enumerate(RBANDS):
        for db, (d0, d1) in enumerate(DBANDS):
            a = deficits.get(("pre", rb, db))
            b = deficits.get(("fix", rb, db))
            if a is None or b is None:
                continue
            P("    r %d-%d  depth %d-%d:  deficit %+.2f %% -> %+.2f %%  "
              "(%+.2f points)%s"
              % (r0, r1, d0, d1, a, b, b - a,
                 "   [P4: not scored]" if db == 0 else ""))
            if db > 0:
                moved.append((d0, d1, a, b))
    if moved:
        fell = [m for m in moved if m[3] < m[2] - 1.0]
        P("")
        P("    bands at depth >= 100 mm: %d scored, %d fell by > 1 point"
          % (len(moved), len(fell)))
        P("    ⇒ %s" % ("the deep deficit FALLS -- P3 confirmed in direction"
                        if len(fell) * 2 >= len(moved) else
                        "⚠⚠ the deep deficit does NOT fall materially: the "
                        "conductance mismatch was NOT the binding term, and "
                        "P1g MUST NOT scale this mode up"))
    else:
        P("    no band at depth >= 100 mm has a scorable deficit on both legs")

    with open(OUTF, "w") as f:
        f.write("\n".join(out) + "\n")
    P("")
    P("written: %s" % OUTF)
    return 0


if __name__ == "__main__":
    sys.exit(main())
