#!/usr/bin/env python3
"""P2-26 m2 -- what moved on each LIVE deck of AB.txt: the last cell plotfile of the BEFORE and the AFTER leg
(ab_legs/<deck>__{before,after}/run), field by field -- which fields differ, in how many cells, by how much -- and
the removal state's summary on both sides: cells removed (sum of `removed` > 0.5), the deepest removed depth below
the initial top (column-wise), and the hottest cell. "Thickness reports only" = the moved fields are confined to
h_spall_field / spall_thickness / RoP_spall with `removed`, `phi`, `Temp`, `H` identical; anything else is a LIVE
change of the trajectory, recorded with its numbers.
Usage: ab_fields.py > AB_FIELDS.txt
"""
import glob
import os
import re

import numpy as np
import yt

yt.set_log_level(50)
HERE = os.path.dirname(os.path.abspath(__file__))
STATE = {"removed", "phi", "Temp", "H"}
REPORT = {"h_spall_field", "spall_thickness", "RoP_spall"}


def last(d):
    p = sorted(glob.glob(os.path.join(d, "run", "*cell")))
    return p[-1] if p else None


def summary(g, dz):
    rem = np.asarray(g[("boxlib", "removed")]) if ("boxlib", "removed") in g.ds.field_list else None
    T = np.asarray(g[("boxlib", "Temp")])
    n_rem = int(np.sum(rem > 0.5)) if rem is not None else -1
    depth = 0.0
    if rem is not None and n_rem > 0:
        depth = float(np.max(np.sum(rem > 0.5, axis=2))) * dz
    return n_rem, depth, float(np.nanmax(T))


def main():
    live = []
    for ln in open(os.path.join(HERE, "AB.txt")):
        m = re.match(r"^(\S+)\s+\(horizon", ln)
        if m:
            cur = m.group(1)
        if "=> LIVE" in ln:
            live.append(cur)
    print("P2-26 m2 -- the LIVE decks of AB.txt, field by field at the last plotfile of each leg")
    for rel in live:
        tag = rel.replace("/", "__")
        pa = last(os.path.join(HERE, "ab_legs", tag + "__before"))
        pb = last(os.path.join(HERE, "ab_legs", tag + "__after"))
        a, b = yt.load(pa), yt.load(pb)
        a.force_periodicity(); b.force_periodicity()   # a domain edge in float round-off trips yt's region check
        dz = float((a.domain_right_edge[2] - a.domain_left_edge[2]) / a.domain_dimensions[2])
        ga = a.covering_grid(0, a.domain_left_edge, a.domain_dimensions)
        gb = b.covering_grid(0, b.domain_left_edge, b.domain_dimensions)
        moved = []
        for f in a.field_list:
            if f[0] != "boxlib" or f not in b.field_list:
                continue
            x, y = np.asarray(ga[f]), np.asarray(gb[f])
            nd = int(np.sum(~((x == y) | (np.isnan(x) & np.isnan(y)))))
            if nd:
                moved.append((f[1], nd, float(np.nanmax(np.abs(x - y)))))
        names = {m[0] for m in moved}
        sa, sb = summary(ga, dz), summary(gb, dz)
        kind = ("thickness reports only (removal state identical)" if names <= REPORT
                else "the TRAJECTORY moved" if names & STATE else "other diagnostics only (removal state identical)")
        print("\n%s  (%s)" % (rel, os.path.basename(pa)))
        print("  fields moved: %s" % ", ".join("%s %d cells max|d| %.3g" % m for m in sorted(moved)))
        print("  removed cells %d -> %d; deepest column %.4g -> %.4g m; hottest cell %.1f -> %.1f K"
              % (sa[0], sb[0], sa[1], sb[1], sa[2], sb[2]))
        print("  => %s" % kind)
    return 0


if __name__ == "__main__":
    main()
