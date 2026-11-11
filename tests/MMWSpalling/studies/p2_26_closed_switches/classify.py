#!/usr/bin/env python3
"""P2-26 b1 -- classify every path the stage-2 tree sweep moved against p2_18.md5.

Bar (the packet's b1, restated for dev2d): every differing path lies inside MOVERS_PRE.txt (the 233 paths of
beam_void_closure + dev2d, written by z0.py before any edit); nothing only on one side; beam_void_closure's thermo.dat
is among them and its header lost exactly `sn_recession_superseded` (dev2d writes no thermo.dat -- z0 section 3);
every differing path is classified; the paths of those two legs that did NOT move are listed. The old bytes are the
archive pre2_out/ (the p2_26a sweep's outputs of the two legs, equal to p2_18's hashes -- checked here).
Classes:
  HDR   thermo.dat: the column sn_recession_superseded left the header (C2's conditional registration); the shared
        columns are compared row by row and the first differing row is reported (the firing moved with the
        /cos(theta) sampling)
  CSV   an event log / cluster log: first differing row, row counts old -> new
  PLT   a plotfile fab (Level_0/Cell_D_*): the field moved because the surface-normal criterion's /cos(theta)
        sampling (site A, the column scan) and the h_col/cos recession are gone
Usage: classify.py NEW_LABEL [REF_LABEL]   (default p2_26 p2_18) > GATE_b1.txt; exit 1 on a failed bar
"""
import hashlib
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
T = os.path.join(ROOT, "tests", "MMWSpalling")
W = os.path.join(T, "studies", "tree_sweep_0922", "witness")
PRE = os.path.join(HERE, "pre2_out")


def load(label):
    m = {}
    for ln in open(os.path.join(W, label + ".md5")):
        h, p = ln.rstrip("\n").split(None, 1)
        m[p.strip()] = h
    return m


def md5(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest()


def rows(p):
    return open(p, errors="replace").read().split("\n")


def main():
    new = sys.argv[1] if len(sys.argv) > 1 else "p2_26"
    ref = sys.argv[2] if len(sys.argv) > 2 else "p2_18"
    N, R = load(new), load(ref)
    movers = [ln[3:].strip() for ln in open(os.path.join(HERE, "MOVERS_PRE.txt")).read().split("\n")[1:] if ln.strip()]
    t0 = {ln[3:].strip() for ln in open(os.path.join(HERE, "MOVERS_PRE.txt")).read().split("\n")[1:]
          if ln.startswith("T0 ")}
    ch = sorted(p for p in N if p in R and N[p] != R[p])
    only_n = sorted(p for p in N if p not in R)
    only_r = sorted(p for p in R if p not in N)
    bad = []
    print("P2-26 b1 -- %s vs %s: %d of %d differ, %d only in %s, %d only in %s" % (
        new, ref, len(ch), len(N), len(only_n), new, len(only_r), ref))
    outside = [p for p in ch if p not in movers]
    if outside or only_n or only_r:
        bad.append("paths outside the pre-registered 233, or on one side only")
    for p in outside[:20]:
        print("  ** OUTSIDE THE PRE-REGISTERED MOVERS: %s" % p)
    # the archive is the old bytes
    arch_ok = all(md5(os.path.join(PRE, p)) == R[p] for p in movers if os.path.exists(os.path.join(PRE, p)))
    n_arch = sum(1 for p in movers if os.path.exists(os.path.join(PRE, p)))
    print("  archive pre2_out/: %d of %d mover paths present, every one hashing to %s: %s"
          % (n_arch, len(movers), ref, arch_ok))
    if not (arch_ok and n_arch == len(movers)):
        bad.append("archive incomplete or not the reference bytes")
    print("\nclassified differing paths (%d):" % len(ch))
    cls = {"HDR": 0, "CSV": 0, "PLT": 0, "OTHER": 0}
    for p in ch:
        old, cur = os.path.join(PRE, p), os.path.join(T, p)
        if p.endswith("thermo.dat"):
            ro, rn = rows(old), rows(cur)
            ho, hn = ro[0].split(), rn[0].split()
            gone = [c for c in ho if c not in hn]
            added = [c for c in hn if c not in ho]
            idx = [(ho.index(c), hn.index(c)) for c in hn if c in ho]
            first = None
            for k in range(1, min(len(ro), len(rn))):
                a, b = ro[k].split(), rn[k].split()
                if not a or not b:
                    continue
                if any(a[i] != b[j] for i, j in idx):
                    first = k
                    break
            ok = gone == ["sn_recession_superseded"] and not added
            cls["HDR"] += 1
            print("  HDR  %-70s header: -%s +%s; rows %d -> %d; shared columns first differ at row %s%s"
                  % (p, gone, added, len(ro) - 1, len(rn) - 1, first, "" if ok else "  ** UNEXPECTED HEADER"))
            if not ok:
                bad.append("thermo header change is not exactly the C2 column")
        elif p.endswith(".csv"):
            ro, rn = rows(old), rows(cur)
            first = next((k for k in range(min(len(ro), len(rn))) if ro[k] != rn[k]), min(len(ro), len(rn)))
            cls["CSV"] += 1
            print("  CSV  %-70s rows %d -> %d; first differing line %d%s"
                  % (p, len(ro) - 1, len(rn) - 1, first, " (header)" if first == 0 else ""))
            if first == 0:
                bad.append("an event-log header moved")
        elif "/Level_0/Cell_D_" in p:
            cls["PLT"] += 1
            step = re.search(r"/(\d+)(cell|node)/Level_0", p)
            print("  PLT  %-70s step %s %s%s" % (p, step.group(1) if step else "?", step.group(2) if step else "",
                                                 "  ** a t = 0 plotfile moved" if p in t0 else ""))
            if p in t0:
                bad.append("a t = 0 plotfile moved")
        else:
            cls["OTHER"] += 1
            print("  ???  %s" % p)
            bad.append("an unclassified path")
    thermo = "unit/beam_void_closure/output/thermo.dat"
    if thermo not in ch:
        bad.append("beam_void_closure's thermo.dat did not move (the column must leave)")
    still = [p for p in movers if p not in ch]
    print("\nclass counts: %s" % cls)
    print("\nmover paths that did NOT move (%d of %d; the t = 0 plotfiles are expected here):" % (len(still), len(movers)))
    for p in still:
        print("  %s %s" % ("T0" if p in t0 else "  ", p))
    by_leg = {}
    for p in ch:
        k = "beam_void_closure" if "beam_void_closure" in p else "dev2d"
        by_leg[k] = by_leg.get(k, 0) + 1
    print("\nmoved per leg: %s; C3 movers (any melt-aware test path): %d" % (
        by_leg, sum(1 for p in ch if not ("beam_void_closure" in p or "/dev2d/" in p))))
    print("\nb1 %s%s" % ("PASS" if not bad else "FAIL", "" if not bad else ": " + "; ".join(sorted(set(bad)))))
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
