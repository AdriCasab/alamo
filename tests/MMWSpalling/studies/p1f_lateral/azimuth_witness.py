#!/usr/bin/env python3
"""P1g-0 deliverable 4 -- the AZIMUTHAL WITNESS, offline, BEFORE the re-run.

Two rows, both scored on files already on disk (P1f's `output_legB`), so they
are measured against the PRE-FIX binary and cannot be confused afterwards with
the fix succeeding or failing:

  (n1) the OPEN-LATERAL-FACE CENSUS per sector x depth band -- the CANDIDATE
       population, unconditioned on firing. This decides whether the deep
       diagonals have NOTHING TO ACT ON (a geometry fact no conductance fix can
       touch) or FACES THAT STAY TOO COOL (what the fix addresses).

  (n2) the VERTICAL-REMOVAL TEMPERATURE TAIL per sector x depth band: p95 and
       the survival fraction above 725 K. This is the RESIDUAL the fix does not
       remove, pre-measured so it cannot later be mistaken for the fix failing.

⚠⚠ THE MEAN AND THE MEDIAN ARE NOT REPORTED, DELIBERATELY. They disagree in
SIGN with the tail: the diagonal median runs ~3.8 K HIGHER while its p95 runs
~2.4 K LOWER (a clipped hot tail, not earlier removal), so a mean-based witness
points the next investigation backwards. The bar is p95 and survival, and only
those are tabulated.

⚠ Every count is printed for ALL SIX sectors and divided by the grid census of
the same sector, so any ratio in the write-up can be recomputed from this
output. The sector definition is printed with the table.

⚠ No decay-length read anywhere: `lat_ell_med` is retired and a two-cell log
ratio cannot resolve ell below ~1 cell. This file scores GEOMETRY and the
TEMPERATURE TAIL only.
"""
import collections, csv, math, os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np
import yt
yt.set_log_level(50)

HERE = os.path.dirname(os.path.abspath(__file__))
LEGB = os.path.join(HERE, "output_legB")          # the PRE-FIX reference
OUTF = os.path.join(HERE, "AZIMUTH_WITNESS.txt")

# ⚠⚠ THE SECTOR DEFINITION, stated because P1f's two conventions gave 0.81 and
# 0.17 for the SAME band (the packet refuses to pre-register 0-100 mm for this
# reason). Here: theta = atan2(y, x) in DEGREES over the quarter domain's full
# 0-90 range, binned in six 15 deg sectors. AXIS = 0-15 and 75-90 (both touch a
# grid axis); DIAGONAL = 30-45 and 45-60 (both touch the 45 deg plane).
SECT6 = [(0, 15), (15, 30), (30, 45), (45, 60), (60, 75), (75, 90)]
AXIS_S = (0, 5)          # indices into SECT6
DIAG_S = (2, 3)
RBANDS = ((40, 60), (60, 90))
DBANDS = ((0, 100), (100, 200), (200, 300), (300, 400))
T_TAIL = 725.0           # the survival threshold, from the axis firing temperature


def plotfiles(base):
    d = os.path.join(base, "FH")
    if not os.path.isdir(d):
        return []
    return sorted(os.path.join(d, x) for x in os.listdir(d)
                  if x.endswith("cell") and os.path.isdir(os.path.join(d, x)))


def sector_of(xc, yc):
    """index into SECT6, or None on the exact origin."""
    if xc <= 0.0 and yc <= 0.0:
        return None
    th = math.degrees(math.atan2(yc, xc))
    for n, (a, b) in enumerate(SECT6):
        if a <= th < b or (b == 90 and th >= 90.0 - 1e-12):
            return n
    return None


def band_of(v, bands):
    for n, (a, b) in enumerate(bands):
        if a <= v < b:
            return n
    return None


def main():
    out = []

    def P(*a):
        s = " ".join(str(x) for x in a)
        print(s)
        out.append(s)

    P("=" * 78)
    P("P1g-0 (n1)/(n2) AZIMUTHAL WITNESS -- OFFLINE, on the PRE-FIX output_legB")
    P("=" * 78)
    P("  source      : %s" % LEGB)
    P("  sectors     : theta = atan2(y,x) deg, quarter domain, six 15 deg bins")
    P("                %s" % "  ".join("%d-%d" % s for s in SECT6))
    P("                AXIS = %s and %s;  DIAGONAL = %s and %s"
      % tuple(["%d-%d" % SECT6[i] for i in (AXIS_S[0], AXIS_S[1],
                                            DIAG_S[0], DIAG_S[1])]))
    P("  ⚠ a count is meaningless without its census; every row prints both.")

    ps = plotfiles(LEGB)
    if not ps:
        P("  FAIL: no plotfile under %s/FH" % LEGB)
        return 1
    ds = yt.load(ps[-1])
    cg = ds.covering_grid(level=0, left_edge=ds.domain_left_edge,
                          dims=ds.domain_dimensions)
    rem = np.asarray(cg["removed"])
    nx, ny, nz = rem.shape
    # ⚠ THE LZ TRAP (standing rule 5): every length comes from the DATASET, not
    # a module-level constant. An imported analyser hard-coding LZ = 0.4 would
    # mis-reduce this 0.6 m domain silently.
    # ⚠⚠ DO NOT unit-convert: yt reports these as `code_length` and maps that
    # to CENTIMETRES, so a .to("cm")/100 round-trip returns 0.0012 for a 0.12 m
    # domain and every radius then lands outside every band -- the whole table
    # came back zeros on the first run. ALAMO's code length IS metres, so the
    # raw numbers are metres. The assertion below is what makes that checkable
    # instead of assumed.
    lo = np.asarray(ds.domain_left_edge)
    hi = np.asarray(ds.domain_right_edge)
    dx = (hi[0] - lo[0]) / nx
    dz = (hi[2] - lo[2]) / nz
    if not (1.9e-3 < dx < 2.1e-3 and 1.9e-3 < dz < 2.1e-3):
        P("  FAIL: dx/dz = %.6g / %.6g m, expected 2 mm -- the domain was "
          "mis-read (unit or LZ trap, standing rule 5)" % (dx, dz))
        return 1
    ktop0 = nz - 1                        # the original rock top
    P("  grid        : %d x %d x %d   dx %.4f mm  dz %.4f mm   t = %.1f s"
      % (nx, ny, nz, dx * 1e3, dz * 1e3, float(ds.current_time)))
    P("  depth of cell k = (%d - k)*dz;  depth bands [mm] %s"
      % (ktop0, DBANDS))

    solid = rem < 0.5

    # ---------------- (n1) the open-face census ---------------------------
    # three nested populations, each strictly inside the one before:
    #   wall  : solid, >=1 in-plane neighbour at the same k is VOID
    #   elig  : wall AND the cell BEHIND that face is solid  (the kernel's own
    #           predicate, Removal.H:566-578 -- what the criterion can judge)
    #   two   : elig with >= 2 open in-plane faces (a 45 deg corner)
    #   lip   : wall whose own top neighbour is void/absent (a fresh step)
    cen = collections.Counter()
    wall = collections.Counter()
    elig = collections.Counter()
    two = collections.Counter()
    lip = collections.Counter()
    nbs = ((-1, 0), (1, 0), (0, -1), (0, 1))
    for i in range(nx):
        xc = lo[0] + (i + 0.5) * dx
        for j in range(ny):
            yc = lo[1] + (j + 0.5) * dx
            s = sector_of(xc, yc)
            if s is None:
                continue
            r = math.hypot(xc, yc) * 1e3
            rb = band_of(r, RBANDS)
            if rb is None:
                continue
            cen[(rb, s)] += 1                       # column census
            for k in range(nz):
                if not solid[i, j, k]:
                    continue
                depth = (ktop0 - k) * dz * 1e3
                db = band_of(depth, DBANDS)
                if db is None:
                    continue
                nopen = 0
                nelig = 0
                for di, dj in nbs:
                    ni, nj = i + di, j + dj
                    if not (0 <= ni < nx and 0 <= nj < ny):
                        continue
                    if solid[ni, nj, k]:
                        continue                     # face not open
                    nopen += 1
                    bi, bj = i - di, j - dj          # the cell BEHIND the face
                    if 0 <= bi < nx and 0 <= bj < ny and solid[bi, bj, k]:
                        nelig += 1
                if nopen == 0:
                    continue
                key = (rb, db, s)
                wall[key] += 1
                if nelig:
                    elig[key] += 1
                if nelig >= 2:
                    two[key] += 1
                if k + 1 >= nz or not solid[i, j, k + 1]:
                    lip[key] += 1

    for rb, (r0, r1) in enumerate(RBANDS):
        P("")
        P("-" * 78)
        P("(n1) OPEN-LATERAL-FACE CENSUS, r = %d-%d mm   [candidates, NOT "
          "conditioned on firing]" % (r0, r1))
        P("     column census per sector: %s"
          % "  ".join("%d-%d:%d" % (SECT6[s][0], SECT6[s][1], cen[(rb, s)])
                      for s in range(6)))
        P("     %-9s %-6s %s" % ("depth[mm]", "pop",
                                 "".join("%9s" % ("%d-%d" % s) for s in SECT6)
                                 + "%10s%10s%9s" % ("axis/col", "diag/col",
                                                    "diag/axis")))
        for db, (d0, d1) in enumerate(DBANDS):
            for nm, tab in (("wall", wall), ("elig", elig),
                            ("2face", two), ("lip", lip)):
                row = [tab[(rb, db, s)] for s in range(6)]
                ca = sum(cen[(rb, s)] for s in AXIS_S)
                cd = sum(cen[(rb, s)] for s in DIAG_S)
                na = sum(row[s] for s in AXIS_S)
                nd = sum(row[s] for s in DIAG_S)
                pa = na / ca if ca else 0.0
                pd = nd / cd if cd else 0.0
                P("     %-9s %-6s %s%10.4f%10.4f%9s"
                  % ("%d-%d" % (d0, d1), nm,
                     "".join("%9d" % v for v in row), pa, pd,
                     ("%8.3f" % (pd / pa)) if pa > 0 else "      --"))
            P("")

    # ---------------- (n2) the vertical-removal temperature tail ----------
    p = os.path.join(LEGB, "FH_removal_events.csv")
    if not os.path.exists(p):
        P("  FAIL: no event log at %s" % p)
        return 1
    rows = list(csv.DictReader(open(p)))
    P("-" * 78)
    P("(n2) VERTICAL-REMOVAL (regime 1) TEMPERATURE TAIL, per sector x depth")
    P("     the quantity is T_top AT FIRING -- the cell temperature the")
    P("     vertical mode removed at. p95 and the survival fraction above")
    P("     %.0f K. ⚠⚠ MEAN AND MEDIAN ARE NOT TABULATED: they disagree in"
      % T_TAIL)
    P("     SIGN with the tail (diagonal median runs HIGHER while its p95 runs")
    P("     LOWER -- a clipped hot tail, not earlier removal), so a mean-based")
    P("     witness points the next investigation backwards.")
    acc = collections.defaultdict(list)
    for r in rows:
        try:
            if int(float(r["regime"])) != 1:
                continue
            i, j = int(r["col_i"]), int(r["col_j"])
            T = float(r["T_top"])
            kt = int(float(r["k_top"]))
            rd = float(r.get("riser_depth") or 0.0)
        except (KeyError, ValueError):
            continue
        xc = lo[0] + (i + 0.5) * dx
        yc = lo[1] + (j + 0.5) * dx
        s = sector_of(xc, yc)
        if s is None:
            continue
        rr = math.hypot(xc, yc) * 1e3
        rb = band_of(rr, RBANDS)
        if rb is None:
            continue
        # ⚠ k_top in the event log is the top BEFORE that removal (standing
        # rule 4), which is exactly the cell that was removed.
        depth = (ktop0 - kt) * dz * 1e3 + rd * 1e3
        db = band_of(depth, DBANDS)
        if db is None:
            continue
        acc[(rb, db, s)].append(T)

    for rb, (r0, r1) in enumerate(RBANDS):
        P("")
        P("     r = %d-%d mm" % (r0, r1))
        P("     %-9s %-9s %s" % ("depth[mm]", "stat",
                                 "".join("%9s" % ("%d-%d" % s) for s in SECT6)
                                 + "%10s%10s%9s" % ("axis", "diag", "d-a")))
        for db, (d0, d1) in enumerate(DBANDS):
            cnt = [len(acc[(rb, db, s)]) for s in range(6)]
            p95 = []
            sur = []
            for s in range(6):
                v = acc[(rb, db, s)]
                p95.append(float(np.percentile(v, 95)) if v else 0.0)
                sur.append(100.0 * sum(1 for x in v if x > T_TAIL) / len(v)
                           if v else 0.0)
            va = [x for s in AXIS_S for x in acc[(rb, db, s)]]
            vd = [x for s in DIAG_S for x in acc[(rb, db, s)]]
            a95 = float(np.percentile(va, 95)) if va else 0.0
            d95 = float(np.percentile(vd, 95)) if vd else 0.0
            asu = 100.0 * sum(1 for x in va if x > T_TAIL) / len(va) if va else 0.0
            dsu = 100.0 * sum(1 for x in vd if x > T_TAIL) / len(vd) if vd else 0.0
            P("     %-9s %-9s %s%10d%10d%9s"
              % ("%d-%d" % (d0, d1), "n",
                 "".join("%9d" % v for v in cnt), len(va), len(vd), ""))
            P("     %-9s %-9s %s%10.1f%10.1f%+9.1f"
              % ("", "p95[K]", "".join("%9.1f" % v for v in p95),
                 a95, d95, d95 - a95))
            P("     %-9s %-9s %s%10.2f%10.2f%+9.2f"
              % ("", ">%.0fK[%%]" % T_TAIL,
                 "".join("%9.2f" % v for v in sur), asu, dsu, dsu - asu))
            P("")

    with open(OUTF, "w") as f:
        f.write("\n".join(out) + "\n")
    P("written: %s" % OUTF)
    return 0


if __name__ == "__main__":
    sys.exit(main())
