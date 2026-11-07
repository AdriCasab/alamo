#!/usr/bin/env python3
"""P2-23 -- RESULTS.txt: the offline z0, the unit and gate logs, the #19 test on the four legs, and the e0 verdict.
Usage: reduce.py > RESULTS.txt (reads only files the chain wrote; no plotfile)."""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def txt(f):
    p = os.path.join(HERE, f)
    return open(p).read().rstrip() if os.path.exists(p) else "NOT ON DISK"


def gate(f, pat):
    return re.search(pat, txt(f), re.M) is not None


def main():
    out = []
    for f in ("Z0.txt", "Z1_face_beam.log", "Z3_retired_keys.log", "GATE_tree.txt", "A1.txt", "A2_beam_void_closure.txt",
              "A2_beam.log", "A2_beam_tau.log", "TEST_all.log"):
        t = txt(f)
        if f in ("Z3_retired_keys.log", "A2_beam.log", "A2_beam_tau.log"):
            t = "\n".join(t.splitlines()[-3:])
        out += ["== %s" % f, t, ""]
    g = {
        "z1, z2": gate("Z1_face_beam.log", r"^PASS$"),
        "z3": gate("Z3_retired_keys.log", r"^PASS"),
        "a0": gate("GATE_tree.txt", r"vs p2_18: 0 of 2639 changed, 0 only in new, 0 only in ref"),
        "a1": gate("A1.txt", r"^a1 PASS"),
        "a2": gate("A2_beam_void_closure.txt", r"^a2 beam_void_closure PASS") and gate("A2_beam.log", r"PASS")
              and gate("A2_beam_tau.log", r"^PASS$"),
        "b1": gate("TEST_all.log", r"PASS  b1"),
    }
    T = txt("TEST_all.log")
    p1 = "P1 HOLDS" in T
    p2 = "P2 HOLDS" in T
    out.append("e0 VERDICT")
    out.append("   the face list runs without a flame and nothing frozen moved: %s (%s)" % (
        "YES" if all(g.values()) else "NO", ", ".join("%s %s" % (k, "PASS" if v else "FAIL") for k, v in g.items())))
    out.append("   P1 (B0: no thin flake before melt, every fired line withdrawn, bk_lb >= 2): %s; P2 (P: the same rule "
               "removes under a shallow beam): %s" % ("holds" if p1 else "does NOT hold", "holds" if p2 else
                                                       "does NOT hold"))
    out.append("   the packet's reading: %s" % (
        "P1 + P2 => the port's first statement stands (no thin-flake spallation under this beam on this rule)"
        if p1 and p2 else
        "P1 without P2 => INERT by the packet's mapping (read the control's buckle log for WHY before accepting it)"
        if p1 else "P2 without P1 => a finding (a thin flake under the deep beam)" if p2 else
        "neither => see the rows"))
    print("\n".join(out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
