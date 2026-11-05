#!/usr/bin/env python3
"""P2-9 rev 2: compare a tree-sweep manifest with reference manifests.
usage: tree_compare.py NEW REF [REF ...] [--declared FILE] [--list]
Every manifest is tree_sweep_0922/witness/<label>.md5 (`md5  path` lines). For each
reference: files whose hash differs, files only in one of the two. With
--declared, every changed path must start with a prefix listed in FILE (one per line,
'#' comments). Prints one line per reference; exit 0."""
import os
import sys

W = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tree_sweep_0922", "witness")


def load(label):
    m = {}
    for ln in open(os.path.join(W, label + ".md5")):
        h, p = ln.rstrip("\n").split(None, 1)
        m[p.strip()] = h
    return m


def main():
    a = [x for x in sys.argv[1:] if not x.startswith("--")]
    declared = None
    if "--declared" in sys.argv:
        f = sys.argv[sys.argv.index("--declared") + 1]
        a.remove(f)
        declared = [l.split("#")[0].strip() for l in open(f) if l.split("#")[0].strip()]
    new, refs = a[0], a[1:]
    N = load(new)
    for r in refs:
        R = load(r)
        ch = sorted(p for p in N if p in R and N[p] != R[p])
        only_n = sorted(p for p in N if p not in R)
        only_r = sorted(p for p in R if p not in N)
        line = "  %s vs %s: %d of %d changed, %d only in new, %d only in ref" % (
            new, r, len(ch), len(N), len(only_n), len(only_r))
        if declared is not None:
            out = [p for p in ch if not any(p.startswith(d) for d in declared)]
            line += "; outside the declared set: %d" % len(out)
            for p in out[:10]:
                line += "\n      UNDECLARED %s" % p
        print(line)
        if "--list" in sys.argv:
            tops = {}
            for p in ch:
                k = "/".join(p.split("/")[:4])
                tops[k] = tops.get(k, 0) + 1
            for k, v in sorted(tops.items()):
                print("      %4d  %s" % (v, k))


if __name__ == "__main__":
    main()
