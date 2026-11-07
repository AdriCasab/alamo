#!/usr/bin/env python3
"""P2-21 a3: C1 (the Zhang deck with the ONE-pair table) against P2-20's legs/Z1/run, byte for byte.

Compared: every file under every NNNNNcell plotfile of the reference (Header, Cell_H, Level_0/Cell_D_*) and
thermo.dat. (P2-20 pruned Z1's node plotfiles, PRUNED.txt, so there is no node reference; metadata, diff.patch and
*.visit carry the argv / git state and are not compared.) Also: the P2-21 columns must be ABSENT from C1's header.
Usage: a3_compare.py > A3.txt (exit 1 unless every compared file is identical and the columns are absent)
"""
import filecmp
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(HERE, "..", "p2_20_beam_deposition", "legs", "Z1", "run")
NEW = os.path.join(HERE, "legs", "C1", "run")


def main():
    pfs = sorted(p for p in os.listdir(REF) if p.endswith("cell"))
    n, bad = 0, []
    for pf in pfs:
        for dp, _ds, fs in os.walk(os.path.join(REF, pf)):
            for fn in fs:
                a = os.path.join(dp, fn)
                b = os.path.join(NEW, os.path.relpath(a, REF))
                n += 1
                if not (os.path.exists(b) and filecmp.cmp(a, b, shallow=False)):
                    bad.append(os.path.relpath(a, REF))
    th = filecmp.cmp(os.path.join(REF, "thermo.dat"), os.path.join(NEW, "thermo.dat"), shallow=False)
    hdr = open(os.path.join(NEW, "thermo.dat")).readline().split("\t")
    cols = [c.strip() for c in hdr if c.strip().startswith(("beam_tau", "beam_alpha_max", "beam_col_err", "bc_"))]
    print("a3: C1 vs p2_20_beam_deposition/legs/Z1/run")
    print("  %d cell plotfiles, %d files: %d identical%s" % (len(pfs), n, n - len(bad),
                                                            ("; DIFFER: " + ", ".join(bad[:6])) if bad else ""))
    print("  thermo.dat byte-identical: %s; P2-21 columns in C1's header: %s" % (th, cols or "none"))
    ok = not bad and n > 0 and th and not cols
    print("a3 %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
