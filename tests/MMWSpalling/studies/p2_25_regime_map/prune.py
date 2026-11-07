#!/usr/bin/env python3
"""P2-25 -- prune legs after reduce (the P2-24 prune.py with the kept set reduced): an md5 manifest of every plotfile
file FIRST (MANIFEST.md5), then keep the first and the last cell plotfile; drop the rest and every node plotfile.
Usage: prune.py LEG [LEG ...] | --all"""
import hashlib
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def prune(name):
    d = os.path.join(HERE, "legs", name)
    run = os.path.join(d, "run")
    if not os.path.isdir(run) or os.path.exists(os.path.join(d, "PRUNED.txt")):
        return
    man = []
    for p in sorted(os.listdir(run)):
        if p.endswith(("cell", "node")) and os.path.isdir(os.path.join(run, p)):
            for dp, _, fs in os.walk(os.path.join(run, p)):
                for fn in sorted(fs):
                    q = os.path.join(dp, fn)
                    man.append("%s  %s" % (hashlib.md5(open(q, "rb").read()).hexdigest(), os.path.relpath(q, run)))
    open(os.path.join(d, "MANIFEST.md5"), "w").write("\n".join(man) + "\n")
    pfs = sorted(p for p in os.listdir(run) if p.endswith("cell"))
    keep = {pfs[0], pfs[-1]} if pfs else set()
    gone = [pf for pf in pfs if pf not in keep]
    for pf in gone:
        shutil.rmtree(os.path.join(run, pf))
    for pf in [p for p in os.listdir(run) if p.endswith("node")]:
        shutil.rmtree(os.path.join(run, pf))
    open(os.path.join(d, "PRUNED.txt"), "w").write(
        "pruned after reduce: kept %s (the first and the last); removed %d cell plotfiles and the node plotfiles; "
        "MANIFEST.md5 lists every plotfile file before the prune\n" % (" ".join(sorted(keep)), len(gone)))


def main(args):
    if args == ["--all"]:
        args = sorted(n for n in os.listdir(os.path.join(HERE, "legs")) if ".old" not in n)
    for n in args:
        prune(n)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
