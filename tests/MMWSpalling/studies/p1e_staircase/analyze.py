#!/usr/bin/env python3
"""P1e criteria (e)-(o) -- a reduction, starts no run.

⚠⚠ THE CONTROL IS P1c, NOT P1b. P1c is `min_sides = 3` with the staircase key
OFF; P1e is the same leg with the key ON, so P1c differs by exactly ONE key and
is the single-variable control. P1b differs by TWO (`min_sides = 2` as well), so
a P1b comparison cannot separate the correction from the support rule. Both are
printed, and every row says which it is against.

⚠ Firings come from the EVENT LOG, never from summing a thermo counter
(`amr.thermo.plot_int` is not 1).
⚠ Every azimuthal count is divided by the grid's own per-sector column census --
1.19x in the 35-44 mm ring, 1.53x on the floor - or it measures the mesh.
"""
import csv, math, os, statistics, sys

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)
DZ = DX = 0.002
KTOP0 = 299
SECTORS = [0, 15, 30, 45, 60, 75]
RING_F = (0.035, 0.044)      # (f) and (o)
RING_M = (0.035, 0.040)      # (m)
DEPTHS = (0.150, 0.200, 0.250)
P1_T_END, P1_FREEZE = 177.616, 120.0

LEGS = (("P1e  staircase=1, sides=3", os.path.join(HERE, "output_fh")),
        ("P1c  staircase=0, sides=3  [THE CONTROL: one key]",
         os.path.join(STUD, "p1c_isotropy", "output_fh")),
        ("P1b  staircase=0, sides=2  [two keys apart]",
         os.path.join(STUD, "p1b_support_rule", "output_fh")))


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
        ev.append({"t": float(d["time"]), "i": int(d["col_i"]),
                   "j": int(d["col_j"]), "reg": int(float(d["regime"])),
                   "nv": int(float(d["n_voided"])),
                   "hf": float(d.get("h_f", 0) or 0.0)})
    ev.sort(key=lambda e: e["t"])
    return ev


def thermo(out, key):
    f = os.path.join(out, "FH", "thermo.dat")
    if not os.path.exists(f):
        return None, None
    L = open(f).read().strip().split("\n")
    h = L[0].split()
    if key not in h:
        return h, None
    i = h.index(key)
    return h, [(float(r.split()[h.index("time")]), float(r.split()[i]))
               for r in L[1:] if len(r.split()) > i]


def census(ring):
    """Per-sector column count of the GRID inside a radial ring -- the
    denominator every azimuthal count needs."""
    c = {s: 0 for s in SECTORS}
    for i in range(0, 60):
        for j in range(0, 60):
            if ring[0] <= r_of(i, j) < ring[1]:
                c[sec(i, j)] += 1
    return c


def main():
    fails = []

    def gate(ok, msg):
        print(("  %s  " % ("ok  " if ok else "FAIL")) + msg)
        if not ok: fails.append(msg)

    cen_f = census(RING_F)
    print("grid census, %.0f-%.0f mm ring: %s   max/min %.2f"
          % (1000 * RING_F[0], 1000 * RING_F[1],
             [cen_f[s] for s in SECTORS],
             max(cen_f.values()) / min(cen_f.values())))

    data = {}
    for lab, out in LEGS:
        ev = load(out)
        if ev is None:
            print("\n%s: MISSING" % lab); continue
        data[lab] = (ev, out)

    main_lab = LEGS[0][0]
    if main_lab not in data:
        print("the P1e leg has no event log yet"); return 1
    T_HI = min(ev[-1]["t"] for ev, _ in data.values())
    print("(f) uses ALL t -- the definition that reproduces the packet's rev 7 "
          "1.01,\n    P1b 1.81 and P1c 1.84. ⚠ A `t > 100 s` cut gives P1b 2.57 "
          "and P1c infinity;\n    the 1.25 bar belongs to the all-t statistic.")

    for lab, (ev, out) in data.items():
        t_end = ev[-1]["t"]
        kt = {}
        for e in ev:
            kt[(e["i"], e["j"])] = kt.get((e["i"], e["j"]), KTOP0) - e["nv"]
        print("\n" + "=" * 78)
        print("%s\n  %d events, t_max %.3f s" % (lab, len(ev), t_end))

        # ---- (f) PER-COLUMN criterion firing rate, 35-44 mm, by sector ----
        cnt = {s: 0 for s in SECTORS}
        for e in ev:
            # ⚠⚠ ALL t, NO TIME CUT. This is the definition that REPRODUCES
            # the packet's own three reference values -- rev 7 1.01, P1b 1.81,
            # P1c 1.84 -- and the 1.25 bar is calibrated against them, so it is
            # the definition the gate must use. ⚠ My first version carried a
            # `t > 100 s` cut over from the packet's Sources table (where it
            # belongs to a DIFFERENT measurement, the ell_eff/T_top asymmetry);
            # under that cut P1b reads 2.57 and P1c reads infinity, and the bar
            # would have been compared against the wrong statistic.
            if e["reg"] != 1:
                continue
            if RING_F[0] <= r_of(e["i"], e["j"]) < RING_F[1]:
                cnt[sec(e["i"], e["j"])] += 1
        rate = [cnt[s] / float(cen_f[s]) for s in SECTORS]
        sp = max(rate) / min(rate) if min(rate) > 0 else float("inf")
        print("  (f) criterion firings PER COLUMN (grid census), %.0f-%.0f mm, "
              "ALL t" % (1000 * RING_F[0], 1000 * RING_F[1]))
        print("      " + "  ".join("%6.1f" % x for x in rate)
              + ("   max/min %.2f  [rev7 1.01, P1b 1.81, P1c 1.84; bar <= 1.25]"
                 % sp))
        if lab == main_lab:
            gate(sp <= 1.25, "(f) per-column firing asymmetry %.2f (bar <= 1.25)" % sp)

        # ---- (g) roundness: radius per sector at fixed depth ----
        print("  (g) max hole radius per sector [mm], bar max/min <= 1.10")
        for d in DEPTHS:
            kcut = KTOP0 - int(round(d / DZ))
            per = {s: [] for s in SECTORS}
            for (i, j), k in kt.items():
                if k <= kcut:
                    per[sec(i, j)].append(r_of(i, j))
            rad = [max(per[s]) if per[s] else None for s in SECTORS]
            if all(rad):
                r_sp = max(rad) / min(rad)
                print("      %3.0f mm: " % (1000 * d)
                      + " ".join("%5.1f" % (1000 * x) for x in rad)
                      + "   max/min %.3f" % r_sp)
                if lab == main_lab:
                    gate(r_sp <= 1.10, "(g) roundness at %.0f mm = %.3f "
                         "(bar <= 1.10)" % (1000 * d, r_sp))
            else:
                print("      %3.0f mm: ⚠ NOT SCORABLE, sectors %s empty"
                      % (1000 * d, [s for s, x in zip(SECTORS, rad) if not x]))

        # ---- (l) (m) the stall ----
        last, final = {}, {}
        for e in ev:
            if e["t"] > P1_T_END: break
            c = (e["i"], e["j"])
            last[c] = max(last.get(c, -1.0), e["t"])
            final[c] = final.get(c, KTOP0) - e["nv"]
        under = [c for c in final if r_of(*c) < 0.040]
        frozen = [c for c in under if last[c] < P1_FREEZE]
        print("  (l) frozen at P1's cut: %d of %d   [P1b 0, P1c 65] -- REPORTED"
              % (len(frozen), len(under)))
        ax = [c for c in final if RING_M[0] <= r_of(*c) < RING_M[1]
              and sec(*c) in (0, 75)]
        axf = [c for c in ax if last[c] < P1_FREEZE]
        print("  (m) AXIS columns frozen, %.0f-%.0f mm: %d of %d = %.1f %%   "
              "[pre-registered 10-12 %%; P1c 10.6 %%]"
              % (1000 * RING_M[0], 1000 * RING_M[1], len(axf), len(ax),
                 100.0 * len(axf) / max(1, len(ax))))

        # ---- (j) (k) ----
        axis = [c for c in kt if r_of(*c) < 0.010]
        dmax = max((KTOP0 - kt[c]) * DZ for c in axis) if axis else 0.0
        print("  (j) axis depth (r < 10 mm): %.1f mm" % (1000 * dmax))
        print("  (k) cells removed, ONE QUADRANT (x4 for Meier): %d"
              % sum(e["nv"] for e in ev))

        # ---- (o) the flake bound's grip ----
        rb = [e for e in ev if e["reg"] == 1
              and RING_F[0] <= r_of(e["i"], e["j"]) < RING_F[1]]
        at = [e for e in rb if abs(e["hf"] - 0.0015) < 1e-9]
        if rb:
            print("  (o) 35-44 mm firings at h_f = 1.5 mm (the flake bound): "
                  "%.1f %%   [P1b 64.4, P1c 43.7]" % (100.0 * len(at) / len(rb)))

        # ---- (e) (h) (i) from thermo ----
        _, led = thermo(out, "ledger_err_solid")
        if led:
            worst = max(abs(v) for _, v in led)
            if lab == main_lab:
                gate(worst <= 5e-3, "(e) worst |ledger_err_solid| %.4e "
                     "(bar 5e-3; P1b 7.60e-4). ⚠ This is the four-site "
                     "consistency check: three of four blows it up." % worst)
            else:
                print("  (e) worst |ledger_err_solid| %.4e" % worst)
        _, scm = thermo(out, "side_sc_mean")
        if scm:
            v = [x for _, x in scm if x > 0.0]
            print("  (h) side_sc_mean: median %.3f, min %.3f   [expect 0.78-0.84]"
                  % (statistics.median(v), min(v)))
        # ---- (n) the D2z floor-step populations, RE-READ ----
        # ⚠⚠ D2z's conclusions were drawn with no staircase factor. They must be
        # RE-DERIVED, not assumed, and a shift is NOT a regression -- the factor
        # falls on the 46 % of floor steps with two or more open faces (the other
        # 54 % are single-face and read exactly 1).
        print("  (n) D2z floor-step populations (median over the run):")
        for key in ("delta_mean", "delta_faces", "step_g_mean", "side_faces",
                    "riser_area_ratio", "wall_step_f_multi"):
            _, ser = thermo(out, key)
            if ser:
                v = [x for _, x in ser]
                print("      %-18s %12.5g" % (key, statistics.median(v)))
        _, p35 = thermo(out, "side_P_in35")
        if p35:
            print("  (i) side_P_in35 per 100 s:", end="")
            for w in range(5):
                seg = [x for t, x in p35 if w * 100 <= t < (w + 1) * 100]
                if seg:
                    print("  %.0f-%.0f %.0f %%" % (w * 100, (w + 1) * 100,
                                                   100 * statistics.median(seg)),
                          end="")
            print("   [rev7 14-17, P1b 26-44, P1c 41->95]")

    print("\n" + "=" * 78)
    if fails:
        for f in fails: print("FAIL: %s" % f)
        return 1
    print("all gated rows PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
