#!/usr/bin/env python3
"""P2-2f (c6) -- bore and Ø(z) against Meier Fig. 8.8. BARE REPORT, not a
prediction: the packet states the radius is substantially taper-set
(ARCHIVE_DONE L12628), so agreement here is not evidence for the property laws.

THREE STANDING RULES GOVERN THIS REDUCTION AND EACH ONE CHANGED THE CODE.

⚠⚠ RULE: SCORE THE BORE FROM THE VOID MASK, NOT THE EVENT LOG. The log's
`k_top` is the top BEFORE that removal, and a column's final state is only
recoverable by replaying every event in order. The mask is the state. This
reads `removed` out of the last cell plotfile.

⚠⚠ RULE: REPORT PER AZIMUTH SECTOR, NEVER AS A SINGLE DIAMETER. The maximum
radius over ALL azimuths picks the grid diagonal and has already produced one
retracted claim. Per-sector maxima are taken over 9 sectors of the quadrant and
the SECTOR MEAN is the reported diameter, with the sector spread beside it.

⚠⚠ RULE: NEVER SCORE DIAMETER AT 100 mm DEPTH -- that is the pause underream
(Meier p. 223). The calibrated section is 150-450 mm, and rows outside it are
printed greyed but excluded from every statistic.

⚠ THE LZ TRAP (standing rule 5) IS WHY THIS IS SELF-CONTAINED. The existing
`d2c_steady`/`d2b_feet_rop` analyzer pair reaches Ø(z) through `run.lz` taken
from a `.done` file that only those harnesses write, and the family has a
history of a module-level `LZ = 0.4`. This deck is 0.6 m tall
(prob_lo -0.1, prob_hi 0.5). Rather than synthesise a `.done` and inherit that
risk, the geometry here is parsed from `feet_argv.txt` -- the SAME file the legs
are launched from, so it cannot drift from what ran.

Usage: bore_profile.py <leg-tag>          # e.g. m13
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "tests", "MMWSpalling", "studies",
                                "p2_2_props"))
import fabcmp  # noqa: E402

FIG = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "meier",
                   "meier_fig8_8_hole_profile.csv")
ARGV_FILE = os.path.join(HERE, "feet_argv.txt")
N_SEC = 9
CAL = (150.0, 450.0)       # mm, the calibrated section (standing rule 2)
UNDERREAM = 100.0          # mm, the pause underream -- never scored
# ⚠ A depth within FLOOR_PAD of the deepest void is the burner's WORKING FLOOR,
# not finished bore: it is still being cut and reads narrow because the taper
# has not opened yet. Measured on rev 1's m=13 log the wall sits at a flat
# 79.3 mm from 125 to 400 mm and then falls to 72.5 at 450 and 49.6 at 480 --
# that fall is the floor, not a bore error, and averaging it in would charge
# the model for rock it has not reached. Rows inside the pad are printed and
# EXCLUDED, with the reason named on the row.
FLOOR_PAD = 25.0           # mm


def geometry():
    """(nx, ny, nz, dx, dz, z_top) from the launch argv. Multi-value keys are
    read as a group: `amr.n_cell=60 60 300` is three values, and rev 1 had a
    parser that kept only the first and silently truncated 100 mm of rock."""
    txt = open(ARGV_FILE).read()

    def vals(key):
        m = re.search(re.escape(key) + r"=((?:\s*-?[0-9][0-9eE.+-]*)+)", txt)
        return [float(x) for x in m.group(1).split()] if m else None
    n = vals("amr.n_cell")
    lo = vals("geometry.prob_lo")
    hi = vals("geometry.prob_hi")
    if not (n and lo and hi) or min(len(n), len(lo), len(hi)) < 3:
        raise SystemExit("could not parse all three geometry keys from %s"
                         % ARGV_FILE)
    dx = (hi[0] - lo[0]) / n[0]
    dz = (hi[2] - lo[2]) / n[2]
    # MMWSpalling.H:8073 initialises phi to a FLAT top surface at prob_hi[2],
    # so depth below the original rock top is z_top - z.
    return int(n[0]), int(n[1]), int(n[2]), dx, dz, hi[2]


def last_plotfile(tag, leg_root=None):
    """⚠ P2-3: `leg_root` lets another study point this at its own leg dirs.
    The reduction is reused UNCHANGED (packet deliverable 11); only the path is
    a parameter, so there is one Ø(z) estimator and not a second copy."""
    d = (os.path.join(leg_root, "run") if leg_root
         else os.path.join(HERE, "ctrl8", tag, "run"))
    if not os.path.isdir(d):
        return None
    cells = sorted(x for x in os.listdir(d) if x.endswith("cell"))
    return os.path.join(d, cells[-1]) if cells else None


def void_mask(pf, nx, ny, nz):
    """Boolean mask[i][j][k] of removed cells, assembled over every rank-file."""
    hdr = os.path.join(pf, "Header")
    lvl = os.path.join(pf, "Level_0")
    files = sorted(x for x in os.listdir(lvl) if x.startswith("Cell_D_"))
    if not files:
        raise SystemExit("no Cell_D_* under %s" % lvl)
    first = fabcmp.read_fabs_boxed(os.path.join(lvl, files[0]))
    names = fabcmp.var_names(hdr, first[0]["ncomp"])
    if "removed" not in names:
        raise SystemExit("no `removed` component; plotfile carries: %s"
                         % ", ".join(names))
    ci = names.index("removed")
    mask = [[[False] * nz for _ in range(ny)] for _ in range(nx)]
    for fn in files:
        for fab in fabcmp.read_fabs_boxed(os.path.join(lvl, fn)):
            lo, hi = fab["lo"], fab["hi"]
            comp = fab["comp"][ci]
            # AMReX FAB ordering is i fastest, then j, then k.
            p = 0
            for k in range(lo[2], hi[2] + 1):
                for j in range(lo[1], hi[1] + 1):
                    for i in range(lo[0], hi[0] + 1):
                        if comp[p] > 0.5 and 0 <= i < nx and 0 <= j < ny \
                                and 0 <= k < nz:
                            mask[i][j][k] = True
                        p += 1
    return mask


def sector_radii(mask, k, nx, ny, dx):
    """Per-sector max void radius at level k. None where a sector is solid."""
    import math
    out = [0.0] * N_SEC
    seen = [False] * N_SEC
    for i in range(nx):
        xc = (i + 0.5) * dx
        for j in range(ny):
            if not mask[i][j][k]:
                continue
            yc = (j + 0.5) * dx
            a = math.atan2(yc, xc)
            s = min(N_SEC - 1, int(a / (0.5 * math.pi) * N_SEC))
            r = (xc * xc + yc * yc) ** 0.5
            seen[s] = True
            if r > out[s]:
                out[s] = r
    return [out[s] if seen[s] else None for s in range(N_SEC)]


def load_fig():
    rows = []
    for ln in open(FIG):
        if ln.startswith("#") or not ln.strip():
            continue
        p = [x.strip() for x in ln.split(",")]
        rows.append((float(p[0]), float(p[1])))
    return rows


def main(tag=None, leg_root=None):
    if tag is None:
        tag = sys.argv[1] if len(sys.argv) > 1 else "m13"
    if leg_root is None and len(sys.argv) > 2:
        leg_root = sys.argv[2]
    nx, ny, nz, dx, dz, z_top = geometry()
    print("P2-2f (c6) -- bore Ø(z) vs Meier Fig. 8.8   [leg %s]" % tag)
    print("=" * 72)
    print("  geometry from %s: %dx%dx%d, dx %.4f m, dz %.4f m, rock top z=%.3f"
          % (os.path.relpath(ARGV_FILE, ROOT), nx, ny, nz, dx, dz, z_top))
    pf = last_plotfile(tag, leg_root)
    if pf is None:
        print("  not computed -- no plotfile for leg %s (root %s)"
              % (tag, leg_root or "ctrl8/<tag>/run"))
        return 0
    print("  mask from %s (`removed` component, NOT the event log)"
          % os.path.relpath(pf, ROOT))
    mask = void_mask(pf, nx, ny, nz)

    # deepest void cell => the working floor.
    d_floor = 0.0
    for k in range(nz):
        z = (z_top - nz * dz) + (k + 0.5) * dz
        if any(mask[i][j][k] for i in range(nx) for j in range(ny)):
            d_floor = max(d_floor, (z_top - z) * 1e3)
    print("  deepest void (working floor) at %.0f mm; depths within %.0f mm of"
          % (d_floor, FLOOR_PAD))
    print("  it are printed but EXCLUDED -- still being cut, not finished bore.")

    fig = load_fig()
    print()
    print("  depth   model Ø   sectors    Meier Ø   diff     [* = outside the")
    print("   mm        mm     min-max      mm       mm       150-450 cal band]")
    diffs = []
    for d_mm, w_mm in fig:
        z = z_top - d_mm * 1e-3
        k = int((z - (z_top - nz * dz)) / dz)
        if not (0 <= k < nz):
            continue
        rs = [r for r in sector_radii(mask, k, nx, ny, dx) if r is not None]
        if not rs:
            print("  %5.0f      solid      --        %5.0f      --" % (d_mm, w_mm))
            continue
        dia = 2.0 * (sum(rs) / len(rs)) * 1e3
        lo, hi = 2e3 * min(rs), 2e3 * max(rs)
        inband = CAL[0] <= d_mm <= CAL[1]
        flag = "" if inband else " *"
        note = ""
        if abs(d_mm - UNDERREAM) < 1.0:
            note = "  <- PAUSE UNDERREAM, never scored"
        elif d_mm > d_floor - FLOOR_PAD:
            note = "  <- WORKING FLOOR, excluded"
        print("  %5.0f    %6.1f   %5.1f-%5.1f   %5.0f   %+6.1f%s%s"
              % (d_mm, dia, lo, hi, w_mm, dia - w_mm, flag, note))
        if inband and not note:
            diffs.append(dia - w_mm)
    print()
    if diffs:
        md = sum(diffs) / len(diffs)
        print("  CALIBRATED SECTION %.0f-%.0f mm, %d depths:" % (CAL[0], CAL[1],
                                                                 len(diffs)))
        print("    mean  Ø error %+.1f mm" % md)
        print("    worst Ø error %+.1f mm" % max(diffs, key=abs))
        print("    ⚠ the deepest in-band depth still sits in the TAPER above")
        print("      the floor (the wall is flat to within 0.1 mm from 125 to")
        print("      400 mm, then narrows), so it is the largest contributor to")
        print("      the mean. It is kept IN: excluding it would IMPROVE the")
        print("      reported error, and a cut that flatters the result needs a")
        print("      better reason than that it flatters the result. FLOOR_PAD")
        print("      is set from the floor position, not from the residual.")
    else:
        print("  not computed -- no calibrated-band depth had void cells")
    print()
    print("  ⚠⚠ BARE REPORT, NOT A PREDICTION. ARCHIVE_DONE L12628: the radius")
    print("     is substantially TAPER-SET, so a small Ø error here is not")
    print("     evidence that the property laws are right, and a large one is")
    print("     not evidence they are wrong. ⚠ Fig. 8.8 itself is a VISIBLE")
    print("     WIDTH on a saw cut that need not pass through the axis, i.e. a")
    print("     LOWER BOUND on diameter, at about +/-5 mm (see the CSV header).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
