#!/usr/bin/env python3
"""P2-27 -- one tree-gate bar: NEW.md5 vs REF.md5 (witness/<label>.md5), path by path.

  usage: tree_bar.py NEW REF [--only-new FILE] [--only-ref FILE]
PASS iff no common path changed, the paths only in NEW are exactly FILE's list (default: none), and the paths only in
REF are exactly FILE's list (default: none). FILE: one path per line (as hash.sh writes them, relative to
tests/MMWSpalling), '#' comments. Exit 0 on PASS.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
W = os.path.join(HERE, "..", "tree_sweep_0922", "witness")


def load(label):
    m = {}
    for ln in open(os.path.join(W, label + ".md5")):
        h, p = ln.rstrip("\n").split(None, 1)
        m[p.strip()] = h
    return m


def listed(path):
    if not path:
        return set()
    return {l.split("#")[0].strip() for l in open(path) if l.split("#")[0].strip()}


def opt(name):
    return sys.argv[sys.argv.index(name) + 1] if name in sys.argv else None


def main():
    new, ref = sys.argv[1], sys.argv[2]
    N, R = load(new), load(ref)
    want_new, want_ref = listed(opt("--only-new")), listed(opt("--only-ref"))
    ch = sorted(p for p in N if p in R and N[p] != R[p])
    on = {p for p in N if p not in R}
    orf = {p for p in R if p not in N}
    print("tree bar %s vs %s: %d / %d paths; %d common changed; only in new %d (declared %d); only in ref %d "
          "(declared %d)" % (new, ref, len(N), len(R), len(ch), len(on), len(want_new), len(orf), len(want_ref)))
    for p in ch[:20]:
        print("  CHANGED   %s" % p)
    for tag, got, want in (("only-new", on, want_new), ("only-ref", orf, want_ref)):
        for p in sorted(got - want)[:20]:
            print("  UNDECLARED %s %s" % (tag, p))
        for p in sorted(want - got)[:20]:
            print("  MISSING    %s %s (declared, not seen)" % (tag, p))
    ok = not ch and on == want_new and orf == want_ref
    print("BAR: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
