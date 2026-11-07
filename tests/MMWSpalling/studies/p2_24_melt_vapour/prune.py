#!/usr/bin/env python3
"""P2-24 -- prune one leg exactly as chain.sh's prune() did: md5 manifest of every plotfile file first, then keep the
first cell plotfile at/after onset, 600 s, 1,580 s and the last; drop the rest and the node plotfiles.
Usage: prune.py LEG [LEG ...]"""
import hashlib
import math
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def prune(name):
    d = os.path.join(HERE, "legs", name)
    run = os.path.join(d, "run")
    L = open(os.path.join(run, "thermo.dat")).read().splitlines()
    h = [c.strip() for c in L[0].split("\t")]
    rows = [[float(x) for x in r.split("\t")] for r in L[1:] if r.strip()]
    t_on = next((r[0] for r in rows if math.isfinite(r[h.index("bc_z_melt")])), None)
    man = []
    for p in sorted(os.listdir(run)):
        if p.endswith(("cell", "node")) and os.path.isdir(os.path.join(run, p)):
            for dp, _, fs in os.walk(os.path.join(run, p)):
                for fn in sorted(fs):
                    q = os.path.join(dp, fn)
                    man.append("%s  %s" % (hashlib.md5(open(q, "rb").read()).hexdigest(), os.path.relpath(q, run)))
    open(os.path.join(d, "MANIFEST.md5"), "w").write("\n".join(man) + "\n")
    pfs = sorted(p for p in os.listdir(run) if p.endswith("cell"))

    def ptime(pf):
        H = open(os.path.join(run, pf, "Header")).read().split("\n")
        return float(H[2 + int(H[1]) + 1])
    keep = {pfs[-1]}
    for t in [x for x in (t_on, 600.0, 1580.0) if x is not None]:
        after = [pf for pf in pfs if ptime(pf) >= t - 1e-9]
        if after:
            keep.add(after[0])
    gone = [pf for pf in pfs if pf not in keep]
    for pf in gone:
        shutil.rmtree(os.path.join(run, pf))
    for pf in [p for p in os.listdir(run) if p.endswith("node")]:
        shutil.rmtree(os.path.join(run, pf))
    open(os.path.join(d, "PRUNED.txt"), "w").write(
        "pruned after its test: kept %s (the first at/after onset %s s, 600 s, 1580 s; the last); removed %d cell "
        "plotfiles and the node plotfiles; MANIFEST.md5 lists every plotfile file before the prune\n"
        % (" ".join(sorted(keep)), t_on, len(gone)))
    print("  %-5s pruned: kept %s" % (name, " ".join(sorted(keep))))


if __name__ == "__main__":
    for n in sys.argv[1:]:
        prune(n)
