#!/usr/bin/env python3
"""P2-2c a1 — the identity comparison, with the doubling fixed.

⚠⚠ THE DEFECT THIS REPLACES. `score_p2_2d.py`'s `manifest_diff` ran
`difflib.unified_diff` over the two sorted manifests and counted every line
starting `+` or `-`. A file whose hash CHANGED contributes BOTH its old line
(`-`) and its new line (`+`), so every count was exactly 2x. It reported
1726/1726/52 where the truth is 863/863/26 -- and at that scale it would have
called P2-2a's 26 strip-leg files 52. The witness scripts only write manifests
and were never at fault.

⇒ count CHANGED PATHS: join the two manifests on the path and count paths whose
hash differs, plus paths present in only one. That is the quantity every
identity claim in this programme has meant.

⚠ The manifest is MD5, not SHA-256. A SHA-256 scorer reports every file as
differing (it did, once).

Usage:  manifest_cmp.py <ref-label> <new-label> [--control <known-different>]
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
W = os.path.join(ROOT, "tests", "MMWSpalling", "studies", "tree_sweep_0922",
                 "witness")


def load(label):
    p = os.path.join(W, label + ".md5")
    if not os.path.exists(p):
        return None, p
    d = {}
    for ln in open(p).read().splitlines():
        if not ln.strip():
            continue
        parts = ln.split(None, 1)
        if len(parts) == 2:
            d[parts[1]] = parts[0]
    return d, p


def compare(ref, new):
    """(changed, only_ref, only_new, n_ref, n_new) -- CHANGED PATHS, not diff
    lines."""
    a, _ = load(ref)
    b, _ = load(new)
    if a is None or b is None:
        return None
    common = set(a) & set(b)
    changed = sum(1 for k in common if a[k] != b[k])
    return (changed, len(set(a) - set(b)), len(set(b) - set(a)),
            len(a), len(b))


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    ref, new = sys.argv[1], sys.argv[2]
    ctrl = None
    if "--control" in sys.argv:
        ctrl = sys.argv[sys.argv.index("--control") + 1]

    fails = []

    def need(ok, msg):
        print("  %s  %s" % ("ok  " if ok else "FAIL", msg))
        if not ok:
            fails.append(msg)

    print("P2-2c IDENTITY COMPARISON (changed paths, not diff lines)")
    # ⚠ the gate must PRINT THE PIN IT COMPARED AGAINST. A bare "0 differing"
    # does not say 0 against what, and two of this programme's identity claims
    # were unreadable for exactly that reason.
    _, pa = load(ref)
    _, pb = load(new)
    print("  reference pin : %s" % os.path.relpath(pa, ROOT))
    print("  new manifest  : %s" % os.path.relpath(pb, ROOT))
    r = compare(ref, new)
    if r is None:
        need(False, "both manifests exist")
        return 1
    ch, only_a, only_b, na, nb = r
    print("  entries: %d ref / %d new" % (na, nb))
    print("  CHANGED PATHS      : %d" % ch)
    print("  only in reference  : %d" % only_a)
    print("  only in new        : %d" % only_b)
    need(ch == 0 and only_a == 0 and only_b == 0,
         "identity: 0 changed paths of %d, and no path in only one" % na)

    if ctrl:
        print("\n  POSITIVE CONTROL against a known-different pin: %s" % ctrl)
        rc = compare(ctrl, new)
        if rc is None:
            need(False, "the control manifest exists")
        else:
            cch, ca, cb, cn, _ = rc
            print("    changed paths %d of %d (only-ref %d, only-new %d)"
                  % (cch, cn, ca, cb))
            # ⚠⚠ THE CONTROL MUST FAIL THE IDENTITY TEST. If it passes, the
            # comparison is blind and the headline 0 means nothing -- which is
            # exactly how a broken `join` once reported 0 from an empty stream.
            need(cch > 0,
                 "the control DOES differ (%d changed paths), so the method "
                 "can see differences at all" % cch)

    print()
    if fails:
        print("FAIL")
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
