#!/usr/bin/env python3
"""Per-COMPONENT comparison of two AMReX multi-FAB `Cell_D_*` files.

P2-2a `(b1)` is barred at component granularity: the 50 strip-derived files may
differ in `Sp_field` ONLY, while `Temp`, `H`, `phi`, `removed`, `eta`,
`T_ref_phase`, `thermal_free_strain`, `beta_phase` and `E_phase` must be
byte-identical per component. A whole-file sha256 cannot express that bar, so
this module reads the FABs and reports which components moved.

⚠ Why this exists at all: rev 5 measured `(b1)` by substituting
`microstructure.phase0.T_ref`, which ALSO redefines `beta_mf`'s secant and the
eigenstrain's `T_ref_phase_mf`. That proxy moved more than change A does, so a
whole-file hash could not tell "the criterion's diagnostic moved" from "the
change leaked into the elastic solve". This tool is the discriminator.

FAB wire format (AMReX, `FArrayBox::writeOn`): an ASCII header line

    FAB ((8, (<bits>)),(8, (<order>)))((lo0,lo1,lo2) (hi0,hi1,hi2) (0,0,0)) <ncomp>

followed immediately by `ncomp * npts` little-endian float64, COMPONENT-MAJOR
(all of component 0 over the box, then component 1, ...). A `Cell_D_*` file
concatenates one such FAB per box owned by that rank-file.
"""

import re
import struct
import sys

HDR = re.compile(
    rb"FAB \(\(8, \([^)]*\)\),\((8|32), \([^)]*\)\)\)"
    rb"\(\((-?\d+),(-?\d+),(-?\d+)\) \((-?\d+),(-?\d+),(-?\d+)\) "
    rb"\((-?\d+),(-?\d+),(-?\d+)\)\) (\d+)\n"
)


def read_fabs(path):
    """-> list of (ncomp, npts, [np-free list of float tuples per component]).

    Returns components as a list of `bytes` slices so an exact comparison is a
    bytes comparison -- no float tolerance is introduced where the bar is
    byte-identity.
    """
    blob = open(path, "rb").read()
    out, pos = [], 0
    while True:
        m = HDR.search(blob, pos)
        if not m:
            break
        lo = tuple(int(m.group(i)) for i in (2, 3, 4))
        hi = tuple(int(m.group(i)) for i in (5, 6, 7))
        ncomp = int(m.group(11))
        npts = 1
        for a, b in zip(lo, hi):
            npts *= (b - a + 1)
        start = m.end()
        need = ncomp * npts * 8
        comp = [blob[start + c * npts * 8: start + (c + 1) * npts * 8]
                for c in range(ncomp)]
        out.append((ncomp, npts, comp))
        pos = start + need
    return out


def read_fabs_boxed(path):
    """-> list of dicts {lo, hi, ncomp, npts, comp: [list-of-float per comp]}.

    ⚠ ADDED, not a change to `read_fabs`: that function's (ncomp, npts, bytes)
    shape is depended on by `unit/pin_clock/test`, where the bar is byte
    identity and decoding to float would weaken it. This variant keeps the BOX
    BOUNDS, which a spatial reduction needs to place a cell, and decodes to
    floats. AMReX FAB ordering is i fastest, then j, then k.
    """
    import struct
    blob = open(path, "rb").read()
    out, pos = [], 0
    while True:
        m = HDR.search(blob, pos)
        if not m:
            break
        lo = tuple(int(m.group(i)) for i in (2, 3, 4))
        hi = tuple(int(m.group(i)) for i in (5, 6, 7))
        ncomp = int(m.group(11))
        npts = 1
        for a, b in zip(lo, hi):
            npts *= (b - a + 1)
        start = m.end()
        need = ncomp * npts * 8
        comp = [list(struct.unpack_from("<%dd" % npts, blob,
                                        start + c * npts * 8))
                for c in range(ncomp)]
        out.append({"lo": lo, "hi": hi, "ncomp": ncomp, "npts": npts,
                    "comp": comp})
        pos = start + need
    return out


def var_names(header_path, ncomp):
    """AMReX plotfile Header: line 0 is the version, line 1 the component
    count, then one variable name per line."""
    lines = open(header_path).read().split("\n")
    n = int(lines[1])
    return lines[2:2 + n] if n == ncomp else lines[2:2 + ncomp]


def compare(path_a, path_b, names=None):
    """-> (ok, moved) where `moved` is a list of (index, name) that differ."""
    fa, fb = read_fabs(path_a), read_fabs(path_b)
    if len(fa) != len(fb):
        return False, [(-1, "FAB COUNT %d vs %d" % (len(fa), len(fb)))]
    ncomp = fa[0][0] if fa else 0
    moved = []
    for c in range(ncomp):
        same = True
        for (na, pa, ca), (nb, pb, cb) in zip(fa, fb):
            if na != nb or pa != pb or ca[c] != cb[c]:
                same = False
                break
        if not same:
            moved.append((c, names[c] if names and c < len(names) else "comp%d" % c))
    return (len(moved) == 0), moved


def main(argv):
    if len(argv) < 3:
        print(__doc__)
        print("usage: fabcmp.py A/Cell_D_00000 B/Cell_D_00000 [A/../Header]")
        return 2
    a, b = argv[1], argv[2]
    names = var_names(argv[3], read_fabs(a)[0][0]) if len(argv) > 3 else None
    ok, moved = compare(a, b, names)
    if ok:
        print("IDENTICAL in all components")
        return 0
    print("MOVED components (%d):" % len(moved))
    for i, n in moved:
        print("  %3d  %s" % (i, n))
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
