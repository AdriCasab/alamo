#!/usr/bin/env python3
"""P2-2d ord 5 — scoring, with (b1)'s POSITIVE CONTROL.

⚠⚠ WHY THE CONTROL EXISTS. P2-2b scored `(b1)` twice as "0 of 2639" from a
`join -j 2 ... -o 0,1.1,2.1` that is INVALID on BSD join: it errored, printed
nothing, and the counter read zero differences from an EMPTY STREAM. A zero
from a passing comparison and a zero from a failed command are
indistinguishable in the output. So every identity score here is
 (i) taken by a sorted `diff` of the two manifests, and
 (ii) accompanied by a control that must find a KNOWN difference.
A control that reads 0 means the method is broken, whatever the headline says.

⚠⚠ P2-2f (a2): THE LOCAL `manifest_diff` IS DELETED. It counted `difflib`
DIFF LINES, and a file whose hash changed contributes BOTH its `-` and its `+`
line, so every number this script ever printed was exactly 2x: it reported
1726/1726/52 where the truth is 863/863/26. The control above was sound; the
counter behind it was not, and the fix shipped in `p2_2c_meier/manifest_cmp.py`
left this script -- and the `SCORE_p2_2d.txt` it backs -- still on the broken
one. Two implementations of one check, one of them wrong, is not an acceptable
end state, so there is now ONE: `compare()` is imported from `manifest_cmp`
and counts CHANGED PATHS.

Usage: score_p2_2d.py > SCORE_p2_2d.txt
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
W = os.path.join(ROOT, "tests", "MMWSpalling", "studies", "tree_sweep_0922",
                 "witness")


sys.path.insert(0, os.path.join(ROOT, "tests", "MMWSpalling", "studies",
                                "p2_2c_meier"))
from manifest_cmp import compare as _compare  # noqa: E402


def manifest_diff(a, b):
    """(changed_paths, n_ref, path_coverage_symmetric_difference).

    ⚠⚠ CHANGED PATHS, via the single implementation in `manifest_cmp.compare`.
    The deleted local version counted difflib `+`/`-` lines and doubled every
    figure. Signature and return shape are kept so the call sites and the
    printed labels below did not have to move -- only the label wording, which
    said "differing lines" and now says what it counts."""
    r = _compare(a, b)
    if r is None:
        return None, None, None
    changed, only_a, only_b, n_a, _n_b = r
    return changed, n_a, only_a + only_b


def main():
    print("P2-2d SCORE")
    print("=" * 72)
    fails = []

    def need(ok, msg):
        print("  %s  %s" % ("ok  " if ok else "FAIL", msg))
        if not ok:
            fails.append(msg)

    print("\n(b1) TREE IDENTITY AT THE SHIPPED DEFAULTS")
    d, n, cov = manifest_diff("p2_2b_final", "p2_2d")
    if d is None:
        need(False, "(b1) a manifest is missing")
    else:
        print("      p2_2b_final vs p2_2d: %d CHANGED PATHS of %d" % (d, n))
        print("      path coverage symmetric difference: %d" % cov)
        need(d == 0, "(b1) 0 changed paths of %d at the shipped defaults" % n)
        need(cov == 0, "(b1) both manifests cover exactly the same paths")
        # ⚠⚠ THE POSITIVE CONTROL -- the method must FIND a known difference.
        print("\n      POSITIVE CONTROL (the method must find real differences):")
        any_found = False
        for old in ("base_d2r0c", "d2t", "d2aa", "p2_2a_r7"):
            dc, nc, _ = manifest_diff(old, "p2_2d")
            if dc is None:
                continue
            print("        %-12s vs p2_2d: %d changed paths of %d" % (old, dc, nc))
            if dc > 0:
                any_found = True
        need(any_found,
             "(b1) the comparison METHOD detects known differences -- without "
             "this the row passes on an empty stream, which it did twice")
        # p2_2a_r7 SHOULD be 0 too (P2-2b was byte-identical), so it is not a
        # control; the three D2-era manifests are.

    print("\n  ⚠ the binary this was taken on is recorded beside the manifest;")
    print("    the guard line in TREE_p2_2d.log asserts it never moved.")
    g = os.path.join(HERE, "TREE_p2_2d.log")
    if os.path.exists(g):
        for ln in open(g):
            if "GUARD" in ln:
                print("    %s" % ln.rstrip())

    print("\n" + "=" * 72)
    if fails:
        print("FAILURES:")
        for f in fails:
            print("  " + f)
        print("\nFAIL")
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
