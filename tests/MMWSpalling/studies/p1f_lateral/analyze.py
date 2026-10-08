#!/usr/bin/env python3
"""P1f Leg B scoring, on the CORRECTED binary (deliverable 9's unit fix).

⚠⚠ Two rules this file exists to enforce:
  * (r)'s DRIVE half is read from the FIELD, unconditioned on firing.
    `lat_ell_med` is RETIRED -- a mean over FIRING columns of a per-column MAX
    over faces, measured to jump +54.8 % at unchanged lat_fires.
  * (r)'s FIRING-RATE half stays on the event log, which IS a rate and is meant
    to be conditioned on firing.
⚠ Every azimuthal count is divided by the grid census, no bare ratios, and any
sector with < 5 firings is refused rather than quoted.
"""
import collections, csv, math, os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np
import yt
yt.set_log_level(50)

HERE = os.path.dirname(os.path.abspath(__file__))
# ⚠⚠ P1g-0: THREE LEGS NOW, and output_legB is NOT moved or overwritten.
#   LEGB_PRE  P1f revision 11, the accepted run -- the PRE-FIX reference. P1f's
#             results, manifest and review findings all cite this path.
#   LEGB_FIX  P1g-0, the same argv with ONE token changed (the plot_file), on
#             the binary with Removal.H:600 passing DX[0] instead of d_n.
#   VOID      the 1e6x unit bug's only witness. NEVER delete it.
# LEGB scores the FIX when it is present and falls back to the pre-fix run, so
# every row below is repointed by this one line rather than by 40 edits.
LEGB_PRE = os.path.join(HERE, "output_legB")
LEGB_FIX = os.path.join(HERE, "output_legB_fix")
LEGB = (LEGB_FIX if os.path.isdir(os.path.join(LEGB_FIX, "FH"))
        else LEGB_PRE)
CTRL = os.path.join(HERE, "output_rev7sc")
VOID = os.path.join(HERE, "output_legB.void_unitbug")
DX = DZ = 0.002
T_REF = 293.15
SCOPE_R = 0.040
BANDS = ((40, 60), (60, 90), (90, 120), (120, 170))
CENSUS_RING, CENSUS_FLOOR = 1.19, 1.53     # for reference; we normalise directly


def thermo(base, key):
    p = os.path.join(base, "FH", "thermo.dat")
    if not os.path.exists(p):
        return [], []
    h = open(p).readline().rstrip("\n").split("\t")
    if key not in h:
        return [], []
    it, ik = h.index("time"), h.index(key)
    t, v = [], []
    for ln in open(p).read().split("\n")[1:]:
        f = ln.split("\t")
        if len(f) <= max(it, ik):
            continue
        t.append(float(f[it])); v.append(float(f[ik]))
    return t, v


def ev(base):
    p = os.path.join(base, "FH_removal_events.csv")
    return list(csv.DictReader(open(p))) if os.path.exists(p) else []


def plotfiles(base):
    d = os.path.join(base, "FH")
    if not os.path.isdir(d):
        return []
    return sorted(os.path.join(d, x) for x in os.listdir(d)
                  if x.endswith("cell") and os.path.isdir(os.path.join(d, x)))


def grid(p):
    ds = yt.load(p)
    cg = ds.covering_grid(level=0, left_edge=ds.domain_left_edge,
                          dims=ds.domain_dimensions)
    return (float(ds.current_time), np.asarray(cg["Temp"]),
            np.asarray(cg["removed"]))


def ktop_of(solid):
    nz = solid.shape[2]
    any_s = solid.any(axis=2)
    return np.where(any_s, nz - 1 - np.argmax(solid[:, :, ::-1], axis=2), -1)


def sect(i, j):
    th = math.degrees(math.atan2(j + 0.5, i + 0.5))
    return min(th, 90.0 - th)          # fold to 0..45: 0 = axis, 45 = diagonal


def main():
    out = []
    P = lambda *a: (print(*a), out.append(" ".join(str(x) for x in a)))
    P("=" * 72)
    P("SCORED LEG : %s" % os.path.basename(LEGB))
    P("  pre-fix  : %s" % os.path.basename(LEGB_PRE))
    P("  control  : %s" % os.path.basename(CTRL))
    P("  ⚠ 'Leg B' below means the SCORED leg named above.")

    # ---------------- (d,e) --------------------------------------------------
    P("=" * 72)
    P("(d,e) COMPLETION AND THE ENERGY INVARIANT")
    log = os.path.join(LEGB, "FH.log")
    txt = open(log).read() if os.path.exists(log) else ""
    tl = [l for l in txt.split("\n") if "TIME =" in l]
    P("  last step line : %s" % (tl[-1] if tl else "(none)"))
    P("  finalized      : %d   aborts/collisions: %d"
      % (txt.count("finalized"),
         sum(txt.lower().count(k) for k in ("abort", "collision",
                                            "ran out of rock"))))
    for tag, base in (("Leg B ", LEGB), ("control", CTRL)):
        t, e = thermo(base, "ledger_err_solid")
        if not t:
            P("  %s : no thermo" % tag); continue
        late = [abs(x) for tt, x in zip(t, e) if tt > 50.0]
        P("  %s : last %+.4e   max|err| t>50s %.3e   bar 5e-3 -> %s"
          % (tag, e[-1], max(late), "PASS" if max(late) <= 5e-3 else "FAIL"))

    # ---------------- (x) removed volume ------------------------------------
    P("=" * 72)
    P("(x) REMOVED VOLUME vs THE CONTROL  [the diagnostic that saw the bug]")
    vol = {}
    for tag, base in (("control", CTRL), ("Leg B", LEGB), ("VOID run", VOID)):
        rows = ev(base)
        if not rows:
            P("  %-9s (no events)" % tag); continue
        v = collections.Counter()
        for r in rows:
            v[int(float(r["regime"]))] += int(float(r["n_voided"]))
        vol[tag] = v
        P("  %-9s %7d events  %8d cells   %s"
          % (tag, len(rows), sum(v.values()),
             "  ".join("r%d:%d" % (g, v[g]) for g in sorted(v))))
    if "control" in vol and "Leg B" in vol:
        c, b = sum(vol["control"].values()), sum(vol["Leg B"].values())
        P("  Leg B / control = %.3f x   [P5: near 1.0, NOT 3.26]" % (b / c))
        if "VOID run" in vol:
            P("  VOID  / control = %.3f x" % (sum(vol["VOID run"].values()) / c))
            r6b = vol["Leg B"].get(6, 0); r6v = vol["VOID run"].get(6, 0)
            if r6b:
                P("  regime-6 cells: %d now vs %d with the bug  (%.0f x fewer)"
                  % (r6b, r6v, r6v / r6b))

    # ---------------- (h4) scope -------------------------------------------
    P("=" * 72
      )
    P("(h4) REGIME-6 FIRINGS INSIDE r < 40 mm  [bar: exactly 0]")
    r6 = [r for r in ev(LEGB) if int(float(r["regime"])) == 6]
    if r6:
        rr = np.array([DX * math.hypot(float(r["col_i"]) + 0.5,
                                       float(r["col_j"]) + 0.5) for r in r6]) * 1e3
        P("  n=%d  min %.1f  p05 %.1f  p50 %.1f  p95 %.1f  max %.1f mm"
          % (len(rr), rr.min(), np.percentile(rr, 5), np.median(rr),
             np.percentile(rr, 95), rr.max()))
        P("  fraction at r < 40 mm: %.4f %%  ->  %s"
          % (100 * (rr < 40).mean(), "PASS" if (rr < 40).sum() == 0 else "FAIL"))
        P("  ⚠⚠ a GREEN scope gate is not evidence of confinement -- it went")
        P("     green on the VOID run too, while the mode ran away OUTSIDE it.")
        P("     (x) and (h1)/(h2) are what say whether it is bounded.")
    else:
        P("  no regime-6 events at all")

    # ---------------- (u-log) what the branch decided ----------------------
    P("=" * 72)
    P("THE BRANCH'S OWN QUANTITIES  [deliverable 10; invisible before]")
    if r6:
        for k, unit in (("q_lat", "-"), ("K_I_lat", "Pa*sqrt(m)"),
                        ("K_Ic_lat", "Pa*sqrt(m)"), ("sigma_lat", "Pa"),
                        ("T_face_lat", "K"), ("ell_lat", "m")):
            if k not in r6[0]:
                continue
            a = np.array([float(r[k]) for r in r6])
            P("  %-11s [%-11s] p05 %.4g  p50 %.4g  p95 %.4g  max %.4g"
              % (k, unit, np.percentile(a, 5), np.median(a),
                 np.percentile(a, 95), a.max()))
        q = np.array([float(r["q_lat"]) for r in r6])
        P("  q_lat < 1 on %d of %d firing events  (must be 0 -- it is the "
          "firing condition)" % (int((q < 1.0).sum()), len(q)))
        kic = np.array([float(r["K_Ic_lat"]) for r in r6])
        P("  K_Ic_lat min %.4g -> %s"
          % (kic.min(), "Pa*sqrt(m), units correct"
             if kic.min() > 1e5 else "⚠⚠ MPa -- THE BUG IS BACK"))

    # ---------------- (h1) propagation -------------------------------------
    P("=" * 72)
    P("(h1) REGIME-6 FIRINGS BY RISER DEPTH x 100 s")
    if r6:
        P("  window      0-5  5-10 10-15 15-20 20-25 25-30 30-35 35-40 40-45  >=45"
          "    p50[mm]  n")
        for w in range(5):
            sel = [r for r in r6 if int(float(r["time"]) // 100) == w]
            c = collections.Counter(min(int(float(r["riser_depth"]) * 1e3 // 5), 9)
                                    for r in sel)
            d = [float(r["riser_depth"]) * 1e3 for r in sel]
            P("  %3d-%3d s  " % (w * 100, (w + 1) * 100)
              + " ".join("%5d" % c.get(b, 0) for b in range(10))
              + "   %6.1f  %d" % (np.median(d) if d else -1, len(sel)))
        dd = np.array([float(r["riser_depth"]) * 1e3 for r in r6])
        P("  all: p50 %.1f  p95 %.1f  max %.1f mm  [VOID run: p95 248, max 374]"
          % (np.median(dd), np.percentile(dd, 95), dd.max()))

    # ---------------- (l) regime x band x 100 s ----------------------------
    P("=" * 72)
    P("(l) REMOVED CELLS BY REGIME x RADIAL BAND x 100 s")
    rows = ev(LEGB)
    tab = collections.Counter()
    for r in rows:
        g = int(float(r["regime"])); w = int(float(r["time"]) // 100)
        rad = DX * math.hypot(float(r["col_i"]) + 0.5,
                              float(r["col_j"]) + 0.5) * 1e3
        for b in BANDS + ((0, 40),):
            if b[0] <= rad < b[1]:
                tab[(g, b, w)] += int(float(r["n_voided"]))
    for g in sorted({k[0] for k in tab}):
        P("  regime %d:" % g)
        for b in ((0, 40),) + BANDS:
            row = [tab.get((g, b, w), 0) for w in range(5)]
            if sum(row):
                P("    %3d-%-3d mm  " % b + " ".join("%7d" % x for x in row)
                  + "   tot %8d" % sum(row))

    # ---------------- (m2) interior cavities -------------------------------
    P("=" * 72)
    P("(m2) INTERIOR CAVITIES  [gate: the depth scan must never read one]")
    for tag, base in (("Leg B", LEGB), ("control", CTRL)):
        pf = plotfiles(base)
        if not pf:
            P("  %s: no plotfiles" % tag); continue
        P("  %s:" % tag)
        worst = 0
        for p in pf:
            t, T, rem = grid(p)
            void = rem > 0.5; solid = ~void
            kt = ktop_of(solid)
            kk = np.arange(rem.shape[2])[None, None, :]
            buried = void & (kk < kt[:, :, None])
            nb = int(buried.sum()); worst = max(worst, nb)
            # deepest contiguous run of buried void below a solid top
            P("    t=%7.2f  buried %8d   void %8d   %6.3f %%"
              % (t, nb, int(void.sum()), 100 * nb / max(int(void.sum()), 1)))
        P("    worst buried-cell count: %d  [VOID run: 72,925 = 25.2 %%]" % worst)

    # ---------------- (r) the azimuthal deficit ----------------------------
    P("=" * 72)
    P("(r) IN-SITU AZIMUTHAL DEFICIT -- drive from the FIELD, rate from the log")
    P("    ⚠ lat_ell_med is RETIRED and is not read here.")
    AX, DG = (0.0, 15.0), (30.0, 45.0)
    # --- rate half, from the event log, census-normalised
    nI = nJ = 200
    cen = collections.Counter()
    for i in range(nI):
        for j in range(nJ):
            rad = DX * math.hypot(i + 0.5, j + 0.5) * 1e3
            th = sect(i, j)
            for b in BANDS:
                if b[0] <= rad < b[1]:
                    if AX[0] <= th < AX[1]:   cen[(b, "ax")] += 1
                    elif DG[0] <= th <= DG[1]: cen[(b, "dg")] += 1
    fir = collections.Counter()
    for r in rows:
        g = int(float(r["regime"]))
        i, j = int(float(r["col_i"])), int(float(r["col_j"]))
        rad = DX * math.hypot(i + 0.5, j + 0.5) * 1e3
        th = sect(i, j)
        for b in BANDS:
            if b[0] <= rad < b[1]:
                k = "ax" if AX[0] <= th < AX[1] else (
                    "dg" if DG[0] <= th <= DG[1] else None)
                if k: fir[(b, k, g)] += 1
    P("  FIRING RATE per column (refusing any sector with < 5 firings):")
    P("    band      reg   axis/col   diag/col   deficit     cols ax/dg   factor vs 41.42%")
    for b in BANDS:
        ca, cd = cen[(b, "ax")], cen[(b, "dg")]
        for g in (1, 6):
            fa, fd = fir[(b, "ax", g)], fir[(b, "dg", g)]
            if fa < 5 or fd < 5:
                P("    %3d-%-3d    %d     -- refused (%d / %d firings) --"
                  % (b[0], b[1], g, fa, fd)); continue
            ra, rd = fa / ca, fd / cd
            df = 100.0 * (ra - rd) / ra
            fac = (41.42 / abs(df)) if abs(df) > 1e-9 else float("inf")
            P("    %3d-%-3d    %d  %9.4f  %9.4f  %+8.2f %%   %5d/%-5d  %8.1f x"
              % (b[0], b[1], g, ra, rd, df, ca, cd, fac))
    # --- drive half, from the FIELD, unconditioned
    pf = plotfiles(LEGB)
    if pf:
        t, T, rem = grid(pf[-1])
        solid = rem <= 0.5
        nx, ny, nz = T.shape
        kt = ktop_of(solid)
        rec = {"ax": [], "dg": []}
        for i in range(1, nx - 1):
            for j in range(1, ny - 1):
                rad = DX * math.hypot(i + 0.5, j + 0.5)
                if rad < SCOPE_R:
                    continue
                th = sect(i, j)
                key = "ax" if AX[0] <= th < AX[1] else (
                    "dg" if DG[0] <= th <= DG[1] else None)
                if key is None:
                    continue
                for k in range(nz):
                    if not solid[i, j, k] or k >= kt[i, j]:
                        continue
                    for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                        ii, jj, bi, bj = i + di, j + dj, i - di, j - dj
                        if not (0 <= ii < nx and 0 <= jj < ny
                                and 0 <= bi < nx and 0 <= bj < ny):
                            continue
                        if solid[ii, jj, k] or not solid[bi, bj, k]:
                            continue
                        dc, db = T[i, j, k] - T_REF, T[bi, bj, k] - T_REF
                        if dc > 0 and db > 0 and dc > db:
                            rec[key].append(DX / math.log(dc / db))
        P("  DRIVE half from the FIELD at t = %.1f s, r >= 40 mm, "
          "unconditioned on firing:" % t)
        if rec["ax"] and rec["dg"]:
            a = np.array(rec["ax"]) * 1e3; d = np.array(rec["dg"]) * 1e3
            P("    axis  n=%5d median %.4f mm | diag n=%5d median %.4f mm | "
              "gap %+.3f %%" % (len(a), np.median(a), len(d), np.median(d),
                                100 * (np.median(d) / np.median(a) - 1)))
            P("    ⚠⚠ under the ABSOLUTE form the criterion is ell-BLIND "
              "(h_f never enters sigma),")
            P("       so this is a CONSISTENCY CHECK, not thesis evidence -- "
              "see the packet on (r).")
        else:
            P("    (no qualifying faces in one of the sectors)")

    # ---------------- (h2) and (p): the bore, per azimuth sector -----------
    P("=" * 72)
    P("(h2)/(p) BORE DIAMETER BY DEPTH AND AZIMUTH SECTOR")
    P("    ⚠⚠ Quarter domain: a radius is one quadrant's, the diameter is 2r.")
    P("    ⚠⚠ NEVER scored at 100 mm depth (pause underream, Meier p. 223);")
    P("       the calibrated section is 150-450 mm. The 50-60 mm MOUTH is")
    P("       reported SEPARATELY -- averaging it into the bore hides the")
    P("       exhaust-heated faces, the likeliest runaway.")
    # ⚠⚠ P1g-0 deliverable 7 -- THE BORE CONVENTION, STATED. The two readings of
    # the same table are the WHOLE of P1f's 98.0-vs-98.2 and 92.6-vs-94.0
    # disagreements, so the band is spelled out rather than implied:
    #
    #   a band (d0, d1) in mm is swept over cells k in [299 - d1/dz,
    #   299 - d0/dz), i.e. `range(k0, k1)` with the UPPER index EXCLUSIVE.
    #   Depth of cell k is (299 - k)*dz, so the band covers depths
    #       d0 + dz  ...  d1   inclusive.
    #   ⇒ THE BAND END IS INCLUDED AND THE BAND START IS EXCLUDED, by one
    #     cell (2 mm). "150-200" is therefore depths 152..200 mm, NOT 150..198.
    #
    # Each cell printed is the MAXIMUM diameter over every void cell in that
    # (sector, band) -- a maximum, never a mean, and 2r because the domain is a
    # quadrant. The depth range actually swept is printed beside each row so the
    # convention can be checked rather than trusted.
    SECT6 = [(0, 15), (15, 30), (30, 45), (45, 60), (60, 75), (75, 90)]
    P("    ⚠⚠ BAND CONVENTION: (d0,d1) sweeps depths d0+dz .. d1 INCLUSIVE --")
    P("       band END included, band START excluded by one cell (dz = 2 mm).")
    P("       '150-200' is depths 152..200 mm. This is P1f's 98.0-vs-98.2.")
    bore = {}
    for tag, base in (("Leg B", LEGB), ("control", CTRL)):
        pf = plotfiles(base)
        if not pf:
            P("  %s: no plotfiles" % tag); continue
        t, T, rem = grid(pf[-1])
        void = rem > 0.5
        nx, ny, nz = rem.shape
        prof = {}
        for kname, krange in (("mouth 50-60", (25, 30)),
                              ("150-200", (75, 100)), ("200-300", (100, 150)),
                              ("300-400", (150, 200)), ("400-450", (200, 225))):
            # depth d maps to k = k_top0 - d/dz, k_top0 = 299
            k0 = 299 - krange[1]; k1 = 299 - krange[0]
            per = []
            for s0, s1 in SECT6:
                rmax = 0.0
                for k in range(max(k0, 0), min(k1, nz)):
                    for i in range(nx):
                        for j in range(ny):
                            if not void[i, j, k]:
                                continue
                            th = math.degrees(math.atan2(j + 0.5, i + 0.5))
                            if not (s0 <= th < s1):
                                continue
                            rr = DX * math.hypot(i + 0.5, j + 0.5)
                            if rr > rmax:
                                rmax = rr
                per.append(2.0 * rmax * 1e3)
            prof[kname] = per
        bore[tag] = prof
        P("  %s (t = %.1f s):" % (tag, t))
        P("    depth[mm]      " + "".join("%8s" % ("%d-%d" % s) for s in SECT6)
          + "    min     max   round")
        for kname, per in prof.items():
            lo = min(x for x in per if x > 0) if any(x > 0 for x in per) else 0.0
            hi = max(per)
            P("    %-14s" % kname + "".join("%8.1f" % x for x in per)
              + "  %6.1f  %6.1f  %6.3f" % (lo, hi, (hi / lo) if lo > 0 else 0.0))
        # the depth range each band ACTUALLY swept, per the stated convention
        # band (c0,c1) in CELLS sweeps k in [299-c1, 299-c0), so the depths are
        # (c0+1)*dz .. c1*dz inclusive -- start excluded, end included.
        P("    depths swept [mm]: " + "  ".join(
            "%s=[%g,%g]" % (kn, (kr[0] + 1) * 2.0, kr[1] * 2.0)
            for kn, kr in (("mouth", (25, 30)), ("150-200", (75, 100)),
                           ("200-300", (100, 150)), ("300-400", (150, 200)),
                           ("400-450", (200, 225)))))
    if "Leg B" in bore and "control" in bore:
        P("  (h2) ceiling check at depths >= 150 mm, bar <= 110 mm:")
        worst = 0.0; worstk = None
        for kname, per in bore["Leg B"].items():
            if kname.startswith("mouth"):
                continue
            if max(per) > worst:
                worst, worstk = max(per), kname
        P("        worst Leg B diameter %.1f mm at %s  ->  %s"
          % (worst, worstk, "PASS" if worst <= 110.0 else "FAIL"))
        P("        mouth 50-60 mm, reported separately: Leg B %.1f, control %.1f mm"
          % (max(bore["Leg B"]["mouth 50-60"]),
             max(bore["control"]["mouth 50-60"])))
        P("  (p) Leg B vs control, per depth band (max diameter):")
        for kname in bore["Leg B"]:
            b, c = max(bore["Leg B"][kname]), max(bore["control"][kname])
            P("        %-14s Leg B %6.1f   control %6.1f   %+6.2f %%"
              % (kname, b, c, 100 * (b / c - 1) if c > 0 else 0.0))

    with open(os.path.join(HERE, "LEGB_SCORED.txt"), "w") as f:
        f.write("\n".join(out) + "\n")
    P("=" * 72)
    P("written: LEGB_SCORED.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
