#!/usr/bin/env python3
"""P2-30: take a pre-image S/<name>/ of the class's source and of the six source-reading tests (refuses to overwrite).
Layout: <name>/mmwspalling.cc, <name>/Integrator/MMWSpalling.H, <name>/Integrator/MMWSpalling/*.H, <name>/tests/<test>/test
Usage: preimage.py <name>"""
import glob
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
SIX = ("lateral_flux_lock", "retired_keys", "buckle_lock", "profile_depth", "face_heat", "face_removal")


def main():
    d = os.path.join(HERE, sys.argv[1])
    if os.path.exists(d):
        sys.exit("preimage: %s exists; not overwriting" % d)
    os.makedirs(os.path.join(d, "Integrator", "MMWSpalling"))
    shutil.copy2(os.path.join(ROOT, "src", "mmwspalling.cc"), d)
    shutil.copy2(os.path.join(ROOT, "src", "Integrator", "MMWSpalling.H"), os.path.join(d, "Integrator"))
    hs = sorted(glob.glob(os.path.join(ROOT, "src", "Integrator", "MMWSpalling", "*.H")))
    for h in hs:
        shutil.copy2(h, os.path.join(d, "Integrator", "MMWSpalling"))
    for t in SIX:
        os.makedirs(os.path.join(d, "tests", t))
        shutil.copy2(os.path.join(ROOT, "tests", "MMWSpalling", "unit", t, "test"), os.path.join(d, "tests", t))
    u = os.path.join(ROOT, "tests", "MMWSpalling", "unit", "_src.py")
    if os.path.exists(u):
        shutil.copy2(u, os.path.join(d, "tests"))
    print("pre-image %s: mmwspalling.cc, MMWSpalling.H, %d fragment(s), %d tests" % (sys.argv[1], len(hs), len(SIX)))


if __name__ == "__main__":
    main()
