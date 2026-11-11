#!/usr/bin/env python3
"""P2-29 b2 -- THE DELETION DID NOT MOVE A PORTED DECK: ZQ and hu_conduction (granite2, sandstone2) re-run on the
stage-2 binary (S/legs2/) against their stage-1 runs (S/legs/), byte for byte: thermo.dat and every file of every
plotfile (Header, Cell_H, Cell_D_*). The run's `metadata` (it embeds the build's git diff), the .visit lists and
diff.patch are not compared. Usage: b2.py   (exit 1 on any difference)"""
import filecmp
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PAIRS = [("ZQ", os.path.join(HERE, "legs", "ZQ", "run"), os.path.join(HERE, "legs2", "ZQ", "run")),
         ("hu_conduction/granite2", os.path.join(HERE, "legs", "hu_new", "hu_conduction", "granite2"),
          os.path.join(HERE, "legs2", "hu_conduction", "granite2")),
         ("hu_conduction/sandstone2", os.path.join(HERE, "legs", "hu_new", "hu_conduction", "sandstone2"),
          os.path.join(HERE, "legs2", "hu_conduction", "sandstone2"))]
SKIP = ("metadata", "diff.patch", "celloutput.visit", "nodeoutput.visit")


def files(root):
    out = []
    for dp, _d, fs in os.walk(root):
        for f in fs:
            if f not in SKIP:
                out.append(os.path.relpath(os.path.join(dp, f), root))
    return sorted(out)


bad = 0
for name, a, b in PAIRS:
    fa, fb = files(a), files(b)
    only = sorted(set(fa) ^ set(fb))
    diff = [f for f in sorted(set(fa) & set(fb)) if not filecmp.cmp(os.path.join(a, f), os.path.join(b, f),
                                                                    shallow=False)]
    ok = not only and not diff and fa
    bad += not ok
    print("  %s  %-26s %d files compared; %s; %s" % ("ok  " if ok else "FAIL", name, len(fa),
                                                    "none on one side only" if not only else "one side only: %s" % only[:4],
                                                    "all byte-identical" if not diff else "DIFFER: %s" % diff[:4]))
print("RESULT b2: %s" % ("PASS (the deletion moved no ported deck)" if not bad else "FAIL (%d)" % bad))
sys.exit(1 if bad else 0)
