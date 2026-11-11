#!/usr/bin/env python3
"""P2-31 stage 3 (goal 11): THE OVERLAP RULE. `Upd` (FaceRemoval.H's per-cell deltas) carried ONE `h_line`, written by
every line through the cell and again by a +z line at its own cell -- last writer wins -- and read by the column's phi
shift AND the plotted thickness / rate. Now: `h_phiz`, written ONLY at the +z line's own cell and read by the phi shift;
`h_max`, the LARGEST h of the lines through the cell, read by the thickness and the two rate fields. A cell one line
passes through reads the same value as before (max(0, h) = h), so a leg with no overlap is unchanged bit for bit.
Every replacement must match exactly once or NOTHING is written."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
F = "src/Integrator/MMWSpalling/FaceRemoval.H"

REPS = [
    ("""            struct Upd { Set::Scalar dH = 0.0, df = 0.0, H0 = 0.0, f0 = 0.0;
                         Set::Scalar h_line = 0.0; int nl = 0; int phiz = 0; };""",
     """            // ⚠ P2-31: THE OVERLAP RULE. Several lines can pass through one cell (fr_overlap: 792 cells in 4 steps
            // on the gate's 45 deg tilt). h_phiz is the +z line's own recession, written ONLY at that line's own cell
            // and read by the column's phi shift; h_max is the LARGEST h of the lines through the cell, read by the
            // plotted thickness and rate. The additive books (df, dH) are unaffected.
            struct Upd { Set::Scalar dH = 0.0, df = 0.0, H0 = 0.0, f0 = 0.0;
                         Set::Scalar h_max = 0.0, h_phiz = 0.0; int nl = 0; int phiz = 0; };"""),
    ("""                    u.f0 = f[static_cast<std::size_t>(d)];
                    u.h_line = h;""",
     """                    u.f0 = f[static_cast<std::size_t>(d)];
                    u.h_max = std::max(u.h_max, h);"""),
    ("""                    upd[key].phiz = 1;
                    upd[key].h_line = h;""",
     """                    upd[key].phiz = 1;
                    upd[key].h_phiz = h;"""),
    ("""                            // a +z line fired on this cell: the whole column's phi drops by h_line
                            if (p[0] >= bx.smallEnd(0) && p[0] <= bx.bigEnd(0) &&
                                p[1] >= bx.smallEnd(1) && p[1] <= bx.bigEnd(1))
                                for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k)
                                    phi(p[0], p[1], k) -= u.h_line;""",
     """                            // a +z line fired on this cell: the whole column's phi drops by that line's own h
                            if (p[0] >= bx.smallEnd(0) && p[0] <= bx.bigEnd(0) &&
                                p[1] >= bx.smallEnd(1) && p[1] <= bx.bigEnd(1))
                                for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k)
                                    phi(p[0], p[1], k) -= u.h_phiz;"""),
    ("""                            ropU(p[0], p[1], p[2]) = (dt > 0.0) ? u.h_line / dt : 0.0;
                            ev_s(p[0], p[1], p[2]) = 1.0;
                            thk_s(p[0], p[1], p[2]) = u.h_line;
                            rop_s(p[0], p[1], p[2]) = (dt > 0.0) ? u.h_line / dt : 0.0;""",
     """                            ropU(p[0], p[1], p[2]) = (dt > 0.0) ? u.h_max / dt : 0.0;
                            ev_s(p[0], p[1], p[2]) = 1.0;
                            thk_s(p[0], p[1], p[2]) = u.h_max;
                            rop_s(p[0], p[1], p[2]) = (dt > 0.0) ? u.h_max / dt : 0.0;"""),
]


def main():
    p = os.path.join(ROOT, F)
    t = open(p).read()
    for old, rep in REPS:
        n = t.count(old)
        if n != 1:
            sys.exit("%s: expected 1 match, found %d of:\n%s" % (F, n, old[:300]))
        t = t.replace(old, rep, 1)
    import re
    left = re.findall(r"\bu\.h_line\b|upd\[key\]\.h_line\b", t)
    if left:
        sys.exit("an Upd h_line access is left: %s" % left)
    open(p, "w").write(t)
    print("3 applied: %s" % F)


if __name__ == "__main__":
    main()
