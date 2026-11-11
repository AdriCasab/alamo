#!/usr/bin/env python3
"""P2-26 b1, the field-level reading of the moved plotfiles: for every cell plotfile of the two C2 legs whose fabs
moved (GATE_b1.txt), which plotfile FIELDS differ between the archive (pre2_out/, = p2_18 bytes) and the stage-2
sweep's output, in how many cells, and by how much.
⚠ Neither leg writes an event log (beam_void_closure's 169 gate files = 42 plotfiles x 4 fabs + thermo.dat; dev2d's
64 = 16 x 4), so "no event log moved" would be vacuous. The bar is on the REMOVAL STATE instead: `removed`, `phi`,
`Temp`, `H` identical in every moved plotfile (the removal is whole-cell, h = dz, whatever the scan reports), and
the moved fields confined to the scan's per-step thickness reports -- `h_spall_field` (the scan ran in the
surface-normal frame, z_crack = s dz cos(theta); now vertical) and, on a step that a plotfile catches mid-firing,
`spall_thickness` / `RoP_spall` (the candidate thickness h_vert = h_col / cos(theta), now h_col).
Usage: field_diff.py > B1_FIELDS.txt
"""
import os
import re

import numpy as np
import yt

yt.set_log_level(50)
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
T = os.path.join(ROOT, "tests", "MMWSpalling")


def main():
    moved = sorted({re.sub(r"/Level_0/.*", "", ln.split()[1]) for ln in open(os.path.join(HERE, "GATE_b1.txt"))
                    if ln.startswith("  PLT ")})
    fields_moved = set()
    print("P2-26 b1 field-level diff: %d moved cell plotfiles" % len(moved))
    for rel in moved:
        a = yt.load(os.path.join(HERE, "pre2_out", rel))
        b = yt.load(os.path.join(T, rel))
        ga = a.covering_grid(0, a.domain_left_edge, a.domain_dimensions)
        gb = b.covering_grid(0, b.domain_left_edge, b.domain_dimensions)
        out = []
        for f in a.field_list:
            if f[0] != "boxlib":
                continue
            x, y = np.asarray(ga[f]), np.asarray(gb[f])
            nd = int(np.sum(~((x == y) | (np.isnan(x) & np.isnan(y)))))
            if nd:
                fields_moved.add(f[1])
                out.append("%s %d cells, max |d| %.4g (max |old| %.4g)" % (f[1], nd, np.nanmax(np.abs(x - y)),
                                                                           np.nanmax(np.abs(x))))
        print("  %-62s %s" % (rel, "; ".join(out) if out else "no field differs"))
    state = {"removed", "phi", "Temp", "H"}
    allowed = {"h_spall_field", "spall_thickness", "RoP_spall"}
    print("fields that moved anywhere: %s" % sorted(fields_moved))
    print("removal-state fields %s moved in any plotfile: %s" % (sorted(state), sorted(fields_moved & state) or "none"))
    print("b1 fields %s" % ("PASS (only the scan's thickness reports moved; the removal state is identical)"
                            if fields_moved <= allowed and not fields_moved & state else "UNEXPECTED"))
    return 0


if __name__ == "__main__":
    main()
