#!/usr/bin/env python3
"""P2-24 -- RESULTS.txt: the offline z0, the unit and gate logs, the #19 test on the four legs, and the e0 verdict.
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
    for f in ("Z0.txt", "Z1_face_melt.log", "Z2_face_vapor.log", "Z3_retired_keys.log", "Z3_latent_heat.log",
              "GATE_tree.txt", "A1.txt", "A2_beam_void_closure.txt", "A2_beam.log", "A2_beam_tau.log",
              "A2_face_beam.log", "CHECK_SHARED.txt", "TEST_all.log"):
        t = txt(f)
        if f in ("Z3_retired_keys.log", "Z3_latent_heat.log", "A2_beam.log", "A2_beam_tau.log"):
            t = "\n".join(t.splitlines()[-3:])
        out += ["== %s" % f, t, ""]
    g = {
        "z1": gate("Z1_face_melt.log", r"^PASS$"),
        "z2": gate("Z2_face_vapor.log", r"^PASS$"),
        "z3": gate("Z3_retired_keys.log", r"^PASS") and gate("Z3_latent_heat.log", r"PASS"),
        "a0": gate("GATE_tree.txt", r"vs p2_18: 0 of 2639 changed, 0 only in new, 0 only in ref"),
        "a1": gate("A1.txt", r"^a1 PASS"),
        "a2": gate("A2_beam_void_closure.txt", r"^a2 beam_void_closure PASS") and gate("A2_beam.log", r"PASS")
              and gate("A2_beam_tau.log", r"^PASS$") and gate("A2_face_beam.log", r"^PASS$"),
        "b1": gate("TEST_all.log", r"PASS  b1 \(i\)") and gate("TEST_all.log", r"PASS  b1 \(ii\)")
              and gate("TEST_all.log", r"PASS  b1 \(iii\)"),
        "b2": gate("TEST_all.log", r"PASS  b2 "),
    }
    out.append("e0 VERDICT")
    out.append("   the melt is a state, the beam sees it, the crack rule stops in and over it, vapour spends once at the "
               "surface only, and nothing frozen moved: %s (%s)" % (
                   "YES" if all(g.values()) else "NO", ", ".join("%s %s" % (k, "PASS" if v else "FAIL")
                                                              for k, v in g.items())))
    out.append("   m1-m7 are REPORTS (the first post-melt statements about #19): TEST_all.log above")
    print("\n".join(out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
