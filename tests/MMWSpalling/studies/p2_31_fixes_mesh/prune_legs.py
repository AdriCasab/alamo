#!/usr/bin/env python3
"""P2-31 DISK: delete the plotfile folders (NNNNNcell / NNNNNnode) of legs whose results are recorded elsewhere -- the z0
smokes (SMOKE_Z0.txt) and list N's legs (byte-identical to P2-30's baselines, which stay). MANIFEST.md5 (md5 of every
file before the prune) is written into the leg first; logs, argv, rc, thermo.dat, metadata and CSVs stay.
Usage: prune_legs.py <leg dir>...   (regenerable: re-run the leg on the binary its .BINARY.txt names, S/bins/)"""
import glob
import hashlib
import os
import re
import shutil
import sys


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


tot = 0
for d in sys.argv[1:]:
    if not os.path.isdir(d):
        continue
    man = []
    for dp, _, fs in os.walk(d):
        for fn in sorted(fs):
            p = os.path.join(dp, fn)
            man.append("%s  %s" % (md5(p), os.path.relpath(p, d)))
    pfs = [p for p in glob.glob(os.path.join(d, "**", "*"), recursive=True)
           if os.path.isdir(p) and re.fullmatch(r"\d{5}(cell|node)", os.path.basename(p))]
    if not pfs:
        continue
    open(os.path.join(d, "MANIFEST.md5"), "w").write("\n".join(man) + "\n")
    n = sum(os.path.getsize(os.path.join(dp, f)) for p in pfs for dp, _, fs in os.walk(p) for f in fs)
    for p in pfs:
        shutil.rmtree(p)
    open(os.path.join(d, "PRUNED.txt"), "w").write("plotfile folders removed after the result was recorded: %d (%.2f GB); "
                                                    "md5 of every file before: MANIFEST.md5\n" % (len(pfs), n / 1e9))
    tot += n
print("pruned %.2f GB" % (tot / 1e9))
