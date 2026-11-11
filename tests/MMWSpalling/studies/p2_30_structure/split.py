#!/usr/bin/env python3
"""P2-30 stage 1 (deliverable 1): cut src/Integrator/MMWSpalling.H into consecutive slices, each #included in place.

A slice is the run of lines from the first comment line (or AMREX_GPU_* line) above its first member to its last
member's closing brace. Nothing is reordered, renamed or edited: the parent with each include line replaced by its slice
is the original file, byte for byte (locks.py bytes proves it). The two access specifiers stay in the parent.

Reads S/src_pre_1 (the pre-image), writes the live src/. Refuses unless the live file is still the pre-image.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
PRE = os.path.join(HERE, "src_pre_1", "Integrator", "MMWSpalling.H")
LIVE = os.path.join(ROOT, "src", "Integrator", "MMWSpalling.H")
DIR = os.path.join(ROOT, "src", "Integrator", "MMWSpalling")
# (slice, first member's first line, last line) -- last None: the line before the next slice's start, blanks stripped.
# Anchors are the packet's, re-resolved at z0.
CUT = [("Parse.H", 81, 1197),
       ("Model.H", 1338, None),
       ("Heat.H", 1537, None),
       ("Flame.H", 2006, None),
       ("FaceList.H", 2919, None),
       ("Ledger.H", 3418, None),
       ("BeamProps.H", 3683, 4473),
       ("ParseRemoval.H", 4484, 4735),
       ("Fields.H", 4746, None),
       ("SurfacePatch.H", 5703, None),
       ("State.H", 6172, 6820)]


def walk_up(L, first):
    """the slice's first line: above the member, every comment / AMREX_GPU_* / template< line and any blank separating
    comment blocks from it, stopping at code"""
    s = first
    i = first - 1
    while i >= 1:
        t = L[i - 1].strip()
        if t.startswith("//") or t.startswith("AMREX_GPU") or t.startswith("template<"):
            s = i
            i -= 1
        elif t == "":
            j = i - 1
            while j >= 1 and L[j - 1].strip() == "":
                j -= 1
            if j >= 1 and L[j - 1].strip().startswith("//"):
                i = j
            else:
                break
        else:
            break
    return s


def main():
    text = open(PRE).read()
    if open(LIVE).read() != text:
        sys.exit("split: the live MMWSpalling.H is not the pre-image; not splitting")
    L = text.split("\n")
    rng = []
    for nm, first, last in CUT:
        rng.append([nm, walk_up(L, first), last])
    for i, r in enumerate(rng):
        if r[2] is None:
            k = rng[i + 1][1] - 1
            while L[k - 1].strip() == "":
                k -= 1
            r[2] = k
    # refuse overlap and out-of-order slices
    for a, b in zip(rng, rng[1:]):
        if not a[1] <= a[2] < b[1]:
            sys.exit("split: bad ranges %s %s" % (a, b))
    out, cur = [], 1
    for nm, s, e in rng:
        out.extend(L[cur - 1:s - 1])
        out.append('    #include "Integrator/MMWSpalling/%s"' % nm)
        cur = e + 1
    out.extend(L[cur - 1:])
    for nm, s, e in rng:
        open(os.path.join(DIR, nm), "w").write("\n".join(L[s - 1:e]) + "\n")
        print("  %-16s lines %5d - %5d  (%4d)" % (nm, s, e, e - s + 1))
    open(LIVE, "w").write("\n".join(out))
    print("  parent: %d lines" % len(out))


if __name__ == "__main__":
    main()
