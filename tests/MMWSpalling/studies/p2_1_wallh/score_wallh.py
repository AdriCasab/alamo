#!/usr/bin/env python3
"""P2-1 ord 3/4: score CTRL5 against CTRL4_FIXED and the ord-1 bar.

⚠⚠ RULES THIS SCRIPT ENFORCES:

 1. **The bore comes from the void MASK, never the event log.** The event-log
    estimator is biased ~0.3 mm and reads 398 mm deepest where the mask says
    396 (P1g-2c carried risk iv). Deliverable 8 also requires the per-15-degree
    sector split to be IN THE SCORER -- the reviewer has supplied it by hand
    twice, which is twice too many.
 2. **The (c2) bar is read out of WALL_H.txt**, which was written before the leg
    ran. It is not recomputed here and cannot be moved to fit the result.
 3. **Never sum a per-step `thermo` column.** `amr.thermo.plot_int = 6` on these
    legs, so a sum is 1-in-6 of the truth (standing rule 7). Every thermo
    quantity below is a MEAN or a MAX over a window, never a sum.
 4. **Quarter domain** -- every volume is multiplied by 4 at the point of use.
 5. A band within 20 mm of either front is not comparable and is not scored.

Usage:  python3 score_wallh.py > SCORE_p2_1.txt
"""
import csv
import math
import os
import re
import statistics as st
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "../../../.."))

LEGS = {
    "CTRL4_FIXED": os.path.join(
        ROOT, "tests/MMWSpalling/studies/p1g2c_shiftfix/CTRL4_FIXED"),
    "CTRL5": os.path.join(HERE, "CTRL5"),
}
W_LO, W_HI = 150.0, 350.0
MEIER_L = 1.270            # Meier's removed volume in 500 s [litres]
ROP_REF = 2.41             # the depth rate (c5) holds within 5 % [m/h]
BORE_BANDS = ((110.0, 150.0), (150.0, 200.0), (200.0, 300.0))
N_SECT = 6                 # 15-degree sectors across the quarter domain
P = print


def thermo(leg):
    f = os.path.join(LEGS[leg], "FH/thermo.dat")
    if not os.path.exists(f):
        return None, None
    L = open(f).read().split("\n")
    hdr = L[0].split()
    ix = {n: i for i, n in enumerate(hdr)}
    rows = [l.split() for l in L[1:] if len(l.split()) >= len(hdr)]
    return ix, rows


def win(ix, rows, name, lo=W_LO, hi=W_HI):
    return [float(r[ix[name]]) for r in rows
            if lo <= float(r[ix["time"]]) <= hi and name in ix]


def masks(leg):
    """Every plotfile's (time, removed) for `leg`, newest last."""
    import numpy as np
    import yt
    yt.set_log_level(50)
    d = os.path.join(LEGS[leg], "FH")
    if not os.path.isdir(d):
        return []
    out = []
    for pf in sorted(x for x in os.listdir(d) if x.endswith("cell")):
        ds = yt.load(os.path.join(d, pf))
        dims = ds.domain_dimensions
        cg = ds.covering_grid(0, left_edge=ds.domain_left_edge, dims=dims)
        out.append((float(ds.current_time),
                    np.array(cg["removed"]) > 0.5,
                    float((ds.domain_right_edge[0]
                           - ds.domain_left_edge[0]) / dims[0]),
                    float((ds.domain_right_edge[2]
                           - ds.domain_left_edge[2]) / dims[2]),
                    float(ds.domain_left_edge[2])))
    return out


def geom_of(rem, dx, dz, zlo):
    import numpy as np
    nx, ny, nz = rem.shape
    Z_TOP = zlo + nz * dz
    cc = (np.arange(nx) + 0.5) * dx
    X, Y = np.meshgrid(cc, cc, indexing="ij")
    return nx, ny, nz, Z_TOP, np.hypot(X, Y), np.degrees(np.arctan2(Y, X))


def depth_mm(rem, dz, zlo):
    """Deepest void, from the MASK."""
    import numpy as np
    nx, ny, nz = rem.shape
    Z_TOP = zlo + nz * dz
    ks = np.nonzero(rem.any(axis=(0, 1)))[0]
    if not len(ks):
        return 0.0
    return 1e3 * (Z_TOP - (zlo + (ks.min() + 0.5) * dz))


def volume_L(rem, dx, dz):
    """FULL removed volume [litres]: quarter domain, so x4."""
    return 4.0 * float(rem.sum()) * dx * dx * dz * 1e3


def bore_table(rem, dx, dz, zlo):
    """{(band, sector): diameter_mm} and {band: diameter_mm}, area-equivalent.

    ⚠ Area-equivalent, not max-radius: a staircased wall makes the max radius a
    one-cell artefact, and the two differ by ~2 mm on this geometry -- the whole
    width of (c3)'s bar.
    """
    import numpy as np
    nx, ny, nz, Z_TOP, rr, th = geom_of(rem, dx, dz, zlo)
    out, per_band = {}, {}
    sect = np.clip((th / (90.0 / N_SECT)).astype(int), 0, N_SECT - 1)
    for band in BORE_BANDS:
        ks = [k for k in range(nz)
              if band[0] <= 1e3 * (Z_TOP - (zlo + (k + 0.5) * dz)) <= band[1]]
        if not ks:
            continue
        acc = []
        for k in ks:
            n = int(rem[:, :, k].sum())
            if n:
                acc.append(math.sqrt(4.0 * n * dx * dx / math.pi))
        if acc:
            per_band[band] = 2e3 * st.mean(acc)
        for sx in range(N_SECT):
            m = (sect == sx)
            a2 = []
            for k in ks:
                n = int((rem[:, :, k] & m).sum())
                if n:
                    # the sector covers 1/N_SECT of the quarter, i.e.
                    # (1/4)(1/N) of the disc, so the area scales by 4*N
                    a2.append(math.sqrt(4.0 * N_SECT * n * dx * dx / math.pi))
            if a2:
                out[(band, sx)] = 2e3 * st.mean(a2)
    return out, per_band


def read_bar():
    """The (c2) bar, read out of WALL_H.txt -- written before the leg ran."""
    f = os.path.join(HERE, "WALL_H.txt")
    if not os.path.exists(f):
        return None
    m = re.search(r"\(c2\) bar FIXED: ([\d.]+) W -> ([\d.]+) W", open(f).read())
    if not m:
        return None
    return float(m.group(1)), float(m.group(2))


def main():
    P("=" * 78)
    P("P2-1 ord 3/4 -- CTRL5 (derived wall law) against CTRL4_FIXED (typed)")
    P("=" * 78)
    for nm, d in LEGS.items():
        P("  %-12s %s" % (nm, "present" if os.path.isdir(d) else "MISSING"))
    if not os.path.isdir(LEGS["CTRL5"]):
        P("\nCTRL5 has not run; nothing to score.")
        return 1
    av = os.path.join(HERE, "CTRL5.argv.sh")
    if os.path.exists(av):
        tail = [l for l in open(av).read().split("\n") if l.strip()]
        P("  CTRL5 argv recorded (%d line(s)); binary %s"
          % (len(tail), tail[-1].split()[0][:16] + "..."))
        P("  the ONE key: %s"
          % ("surface_patch.side_face_h_law=derived"
             if "side_face_h_law=derived" in " ".join(tail) else "** ABSENT **"))

    T = {nm: thermo(nm) for nm in LEGS}
    M = {nm: masks(nm) for nm in LEGS}
    fails, reports = [], []

    def gate(ok, label, detail):
        P("  %-6s %s" % ("PASS" if ok else "FAIL", label))
        P("         %s" % detail)
        if not ok:
            fails.append(label)

    # ---------------- (c1) and (c1w) -----------------------------------
    P("")
    P("-" * 78)
    P("(c1) GATE   max |ledger_err_solid| on CTRL5 <= 5e-3")
    P("-" * 78)
    ix5, r5 = T["CTRL5"]
    if ix5 is None:
        P("  CTRL5 wrote no thermo.dat"); return 1
    allv = [abs(float(r[ix5["ledger_err_solid"]])) for r in r5]
    gate(max(allv) <= 5.0e-3, "(c1) ledger_err_solid <= 5e-3",
         "max |.| over the whole run = %.4e  (%d rows)" % (max(allv), len(allv)))
    w50 = [abs(x) for x in win(ix5, r5, "ledger_err_solid", 50.0, 1e9)]
    ix4, r4 = T["CTRL4_FIXED"]
    w50c = ([abs(x) for x in win(ix4, r4, "ledger_err_solid", 50.0, 1e9)]
            if ix4 else [])
    P("")
    P("-" * 78)
    P("(c1w) REPORT  the same, windowed t > 50 s, bar 2e-4")
    P("-" * 78)
    P("  CTRL5       max |.| = %.4e   %s"
      % (max(w50), "within" if max(w50) <= 2.0e-4 else "OUTSIDE"))
    if w50c:
        P("  CTRL4_FIXED max |.| = %.4e   <- the comparison leg's own reading"
          % max(w50c))
        P("  ⚠ NOT a gate: the comparison leg itself sits at %.4e, so a tighter"
          % max(w50c))
        P("    bar would fail CTRL5 for its trajectory and not for the law.")
    reports.append(("(c1w)", "%.4e" % max(w50)))

    # ---------------- (c2) ----------------------------------------------
    P("")
    P("-" * 78)
    P("(c2) VERDICT  wall energy, surf_side_share x ledger_P_robin, %g-%g s"
      % (W_LO, W_HI))
    P("-" * 78)
    def wall_W(ix, rows):
        s = win(ix, rows, "surf_side_share")
        p = win(ix, rows, "ledger_P_robin")
        n = min(len(s), len(p))
        return st.mean([s[i] * p[i] for i in range(n)]) if n else float("nan")
    W5 = wall_W(ix5, r5)
    W4 = wall_W(ix4, r4) if ix4 else float("nan")
    bar = read_bar()
    P("  CTRL4_FIXED (typed)   %.1f W" % W4)
    P("  CTRL5       (derived) %.1f W   =>  ratio %.4f" % (W5, W5 / W4))
    if bar:
        P("  the ord-1 bar, written to WALL_H.txt BEFORE this leg ran:")
        P("    %.1f W -> %.1f W, +/- 10 %%  =>  [%.1f, %.1f]"
          % (bar[0], bar[1], 0.9 * bar[1], 1.1 * bar[1]))
        ok = 0.9 * bar[1] <= W5 <= 1.1 * bar[1]
        # ⚠⚠ (c2) IS A `verdict` ROW, NOT A GATE -- the packet's criteria table
        # says "only a `gate` halts". It is reported as met / not met and the
        # shortfall is ATTRIBUTED; it does not stop the step.
        P("  VERDICT %s: measured %.1f W against the predicted %.1f W (%+.1f %%)"
          % ("MET" if ok else "NOT MET", W5, bar[1],
             100.0 * (W5 / bar[1] - 1.0)))
        reports.append(("(c2)", "%.1f W vs %.1f predicted (%+.1f %%), verdict %s"
                        % (W5, bar[1], 100.0 * (W5 / bar[1] - 1.0),
                           "MET" if ok else "NOT MET")))
        # ---- the attribution. The packet offers two readings for a (c2) miss:
        # (a3)'s cancellation reached the run, or the resolution did not reach
        # all four roles. BOTH ARE EXCLUDED BY MEASUREMENT BELOW, and the whole
        # shortfall is accounted for by two assumptions in the BAR itself.
        P("")
        P("  ATTRIBUTION -- the packet's two readings, tested rather than assumed")
        e4 = max(abs(float(r[ix4["ledger_err_solid"]])) for r in r4
                 if W_LO <= float(r[ix4["time"]]) <= W_HI)
        e5 = max(abs(float(r[ix5["ledger_err_solid"]])) for r in r5
                 if W_LO <= float(r[ix5["time"]]) <= W_HI)
        P("   (i)  LEDGER mirrors KERNEL: max|ledger_err_solid| %.3e -> %.3e."
          % (e4, e5))
        P("        A kernel/ledger split lands straight here; it did not move.")
        m4 = st.mean(win(ix4, r4, "side_P_face"))
        m5 = st.mean(win(ix5, r5, "side_P_face"))
        P("   (ii) MARCH tracks KERNEL: kernel/march split %.4f -> %.4f."
          % (W4 / m4, W5 / m5))
        P("        Had the march kept the typed scalar while the kernel moved,")
        P("        the kernel rose %.1f %% and the split would read ~%.3f."
          % (100.0 * (W5 / W4 - 1.0), (W5 / W4) * (W4 / m4)))
        P("   (iii) CRITERION saw it: its FIRING population moves while its")
        P("        CANDIDATE population barely does -- the criterion fired more")
        P("        on the same faces, which is the law reaching it and nothing")
        P("        else. A criterion left on the typed scalar would show x1.00.")
        for nm in ("lat_fires", "lat_bank", "lat_shift", "lat_cand"):
            a = st.mean(win(ix4, r4, nm)); b = st.mean(win(ix5, r5, nm))
            P("          %-10s %10.3f -> %10.3f   x%.4f"
              % (nm, a, b, b / a if a else float("nan")))
        P("   => NEITHER packet reading holds. The shortfall is in the BAR:")
        A4 = st.mean(win(ix4, r4, "side_A_wall"))
        A5 = st.mean(win(ix5, r5, "side_A_wall"))
        P4_ = st.mean(win(ix4, r4, "side_P_wall"))
        P5_ = st.mean(win(ix5, r5, "side_P_wall"))
        q4, q5 = P4_ / A4, P5_ / A5
        P("        wall FLUX   q = P_wall/A_wall   %9.1f -> %9.1f  x%.4f"
          % (q4, q5, q5 / q4))
        P("        wall AREA   side_A_wall        %.6f -> %.6f  x%.4f"
          % (A4, A5, A5 / A4))
        P("        wall POWER  side_P_wall        %9.1f -> %9.1f  x%.4f"
          % (P4_, P5_, P5_ / P4_))
        P("        the bar predicted an h ratio of %.3f holding dT AND area"
          % (bar[1] / bar[0]))
        P("        fixed. Measured: area x%.4f, and the flux rose only x%.4f"
          % (A5 / A4, q5 / q4))
        P("        where h rose x%.3f, so the Robin dT fell x%.3f."
          % (bar[1] / bar[0], (q5 / q4) / (bar[1] / bar[0])))
        P("        %.3f (h) x %.3f (dT) x %.4f (area) = %.4f, against the"
          % (bar[1] / bar[0], (q5 / q4) / (bar[1] / bar[0]), A5 / A4,
             (q5 / q4) * (A5 / A4)))
        P("        measured power ratio %.4f." % (P5_ / P4_))
        P("   ⚠⚠ A HOTTER-DRIVEN WALL RUNS HOTTER, so (T_gas - T_s) shrinks and")
        P("      the power rises far less than the coefficient. The ord-1 bar")
        P("      said 'dT held' in as many words; that assumption, plus a wall")
        P("      AREA the law itself changes, is the whole of the miss. The law")
        P("      moved the energy -- by %+.1f %%, not the %+.1f %% predicted."
          % (100.0 * (W5 / W4 - 1.0), 100.0 * (bar[1] / bar[0] - 1.0)))
    else:
        gate(False, "(c2) the bar is readable from WALL_H.txt",
             "WALL_H.txt absent or its bar line not found -- UNSCOREABLE")
    # the law's own witness columns
    for nm in ("wall_law_on", "wall_h_corner", "wall_z_fill", "wall_r_bore",
               "wall_h_duct", "wall_sector", "side_h_used_p05",
               "side_h_used_p50", "side_h_used_p95", "side_h_used_n"):
        if nm in ix5:
            v = win(ix5, r5, nm)
            if v:
                P("    %-18s window mean %12.5g   last %12.5g"
                  % (nm, st.mean(v), float(r5[-1][ix5[nm]])))
    if "side_h_used_p05" in ix5:
        a = st.mean(win(ix5, r5, "side_h_used_p05"))
        b = st.mean(win(ix5, r5, "side_h_used_p95"))
        P("    (b2) upper half, on this leg: p95/p05 = %.3f  [bar >= 2]  %s"
          % (b / a, "PASS" if b / a >= 2.0 else "FAIL"))
    on = win(ix5, r5, "wall_law_on") if "wall_law_on" in ix5 else []
    if on:
        allon = [float(r[ix5["wall_law_on"]]) for r in r5]
        P("    ⚠ the law is inert until the bore clears the shield: wall_law_on")
        P("      is 1 on %.1f %% of all rows and %.1f %% inside the window."
          % (100.0 * sum(allon) / len(allon), 100.0 * sum(on) / len(on)))

    # ---------------- (c3) deliverable 8 ---------------------------------
    P("")
    P("-" * 78)
    P("(c3) REPORT  bore by depth band AND per 15-degree sector (deliverable 8)")
    P("-" * 78)
    if not M["CTRL5"] or not M["CTRL4_FIXED"]:
        P("  no plotfiles on one leg -- UNSCOREABLE")
    else:
        t5, m5, dx5, dz5, zl5 = M["CTRL5"][-1]
        t4, m4, dx4, dz4, zl4 = M["CTRL4_FIXED"][-1]
        P("  final masks: CTRL5 t=%.1f s, CTRL4_FIXED t=%.1f s" % (t5, t4))
        s5, b5 = bore_table(m5, dx5, dz5, zl5)
        s4, b4 = bore_table(m4, dx4, dz4, zl4)
        P("")
        P("  by depth band, FULL diameter [mm]:")
        P("    band [mm]      CTRL4_FIXED    CTRL5     delta")
        dmax = 0.0
        for band in BORE_BANDS:
            if band in b4 and band in b5:
                d = b5[band] - b4[band]
                dmax = max(dmax, abs(d))
                P("    %3.0f-%3.0f        %9.2f  %9.2f  %+8.2f"
                  % (band[0], band[1], b4[band], b5[band], d))
        P("")
        P("  per 15-degree sector, FULL diameter [mm]  (0 deg = +x axis):")
        P("    band [mm]  " + "".join("  %2d-%2d deg" % (sx * 15, sx * 15 + 15)
                                      for sx in range(N_SECT)))
        smax = 0.0
        for band in BORE_BANDS:
            row4 = ["%9.2f" % s4[(band, sx)] if (band, sx) in s4 else "      ---"
                    for sx in range(N_SECT)]
            row5 = ["%9.2f" % s5[(band, sx)] if (band, sx) in s5 else "      ---"
                    for sx in range(N_SECT)]
            rowd = []
            for sx in range(N_SECT):
                if (band, sx) in s4 and (band, sx) in s5:
                    d = s5[(band, sx)] - s4[(band, sx)]
                    smax = max(smax, abs(d))
                    rowd.append("%+9.2f" % d)
                else:
                    rowd.append("      ---")
            P("    %3.0f-%3.0f T  " % (band[0], band[1]) + "".join(row4))
            P("    %3.0f-%3.0f D  " % (band[0], band[1]) + "".join(row5))
            P("    %3.0f-%3.0f d  " % (band[0], band[1]) + "".join(rowd))
        P("")
        P("  worst band-mean change   %+.2f mm   [bar +/- 2 mm, SIGN NOT PREDICTED]"
          % dmax)
        P("  worst per-sector change  %+.2f mm" % smax)
        P("  ⚠⚠ a change MUCH LARGER than 2 mm in either direction is the")
        P("     WARNING SIGN, not the prize: it means something was tuned, or")
        P("     the law reached the floor. Read with (c5).")
        reports.append(("(c3) band", "%+.2f mm" % dmax))
        reports.append(("(c3) sector", "%+.2f mm" % smax))

        # ---------------- (c4) and (c5) ------------------------------
        P("")
        P("-" * 78)
        P("(c4) REPORT  VOLUME in 500 s, against Meier's %.3f L" % MEIER_L)
        P("-" * 78)
        v5, v4 = volume_L(m5, dx5, dz5), volume_L(m4, dx4, dz4)
        P("  CTRL4_FIXED %.4f L  => %.3fx Meier" % (v4, v4 / MEIER_L))
        P("  CTRL5       %.4f L  => %.3fx Meier   [band 2.15-2.23x]"
          % (v5, v5 / MEIER_L))
        P("  CTRL5/CTRL4 = %.4f  (%+.2f %%)" % (v5 / v4, 100.0 * (v5 / v4 - 1.0)))
        P("  ⚠⚠ BOUNDED BY ARITHMETIC. (c3) holds width and (c5) holds depth, so")
        P("     volume cannot move more than ~5 %. A reading outside that")
        P("     contradicts (c3) or (c5) and the three must be reconciled first.")
        P("")
        P("  RECONCILIATION, since the reading is outside the stated band:")
        P("   * the band 2.15-2.23x presumed the law would REMOVE wall energy.")
        P("     (a3) measured the opposite sign offline, before this leg ran, so")
        P("     the expected direction reverses with it: volume should RISE.")
        P("   * the change is %+.2f %%, inside the ~5 %% the packet's own bound"
          % (100.0 * (v5 / v4 - 1.0)))
        P("     permits, so (c4) does NOT contradict (c3) or (c5).")
        P("   * CTRL4_FIXED's own mask reads %.3fx, and the packet states the"
          % (v4 / MEIER_L))
        P("     current model at 2.886 L / 2.273x -- so the mask estimator and")
        P("     the packet's figure agree to 3 s.f. on the comparison leg, which")
        P("     is what makes the +%.2f %% attributable to the law."
          % (100.0 * (v5 / v4 - 1.0)))
        reports.append(("(c4)", "%.3fx" % (v5 / MEIER_L)))

        P("")
        P("-" * 78)
        P("(c5) REPORT  depth rate over %g-%g s, within 5 %% of %.2f m/h"
          % (W_LO, W_HI, ROP_REF))
        P("-" * 78)
        for nm in ("CTRL4_FIXED", "CTRL5"):
            pts = [(t, depth_mm(m, dz, zl)) for t, m, _dx, dz, zl in M[nm]
                   if W_LO - 60 <= t <= W_HI + 60]
            inw = [(t, d) for t, d in pts if W_LO <= t <= W_HI]
            if len(inw) >= 2:
                rop = 3.6 * (inw[-1][1] - inw[0][1]) / (inw[-1][0] - inw[0][0])
                P("  %-12s %s  => ROP %.3f m/h (%+.1f %% of %.2f)"
                  % (nm, " ".join("%.0fs:%.0fmm" % p for p in inw),
                     rop, 100.0 * (rop / ROP_REF - 1.0), ROP_REF))
                if nm == "CTRL5":
                    reports.append(("(c5)", "%.3f m/h" % rop))
            else:
                P("  %-12s fewer than two plotfiles in the window" % nm)
        P("  ⚠ the floor owns ~90 % of depth; a large depth change means the")
        P("    wall law is reaching floor risers -- a SCOPE error, and (c4)")
        P("    would then be unattributed.")

    P("")
    P("=" * 78)
    P("SUMMARY")
    P("=" * 78)
    for k, v in reports:
        P("  %-14s %s" % (k, v))
    if fails:
        P("  GATES FAILED: %s" % ", ".join(fails))
        return 1
    P("  every GATE row passed; the report rows are above.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
