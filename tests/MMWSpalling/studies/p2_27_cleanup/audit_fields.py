#!/usr/bin/env python3
"""P2-27 L9 / d2 -- THE FIELD AUDIT of stage 4 (the grimsel-inert scalars refused; a law phase's t = 0 property
fields filled from the law): two face-gate lists, every differing plotfile read component by component.

  usage: audit_fields.py <old> <new> [--declared MOVERS_4_face_t0.txt]
PASS iff the differing paths are EXACTLY the declared list (the t = 0 plotfile's Cell_D + Cell_H of every grimsel leg;
nothing only on one side), every differing plotfile keeps its component list, and the components that differ are a
subset of {kappa_phase, beta_phase, E_phase, mu_phase}. Prints, per leg, the differing fields with max |d| and the
old -> new value range of each (the scalar vs the law at T_ref).
"""
import importlib.util
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("ac", os.path.join(HERE, "audit_components.py"))
ac = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ac)
ALLOWED = {"kappa_phase", "beta_phase", "E_phase", "mu_phase"}


def fields(pf):
    _, names = ac.read_header(pf)
    out = {n: [] for n in names}
    cache = {}
    for f, off in ac.read_cell_h(pf):
        if f not in cache:
            cache[f] = ac.read_fabs(pf, f)
        npts, nc, data, box = cache[f][off]
        for c, n in enumerate(names):
            out[n].append(np.array(data[c]))
    return {n: np.concatenate(v) for n, v in out.items()}


def main():
    old, new = sys.argv[1], sys.argv[2]
    decl = sys.argv[sys.argv.index("--declared") + 1] if "--declared" in sys.argv else \
        os.path.join(HERE, "MOVERS_4_face_t0.txt")
    want = {l.strip() for l in open(decl) if l.strip() and not l.startswith("#")}
    A = ac.load_list(os.path.join(ac.W, "face_%s.md5" % old))
    B = ac.load_list(os.path.join(ac.W, "face_%s.md5" % new))
    # --ignore-new LEGS: legs added to the gate after <old> ran (P2-27 stage 4 rev 2 ran face_gate.py with the stage-5
    # FEET legs already listed); their paths are not part of this bar
    ign = sys.argv[sys.argv.index("--ignore-new") + 1].split(",") if "--ignore-new" in sys.argv else []
    B = {p: h for p, h in B.items() if not any(p.startswith(l + "/") for l in ign) or p in A}
    only = sorted(set(A) ^ set(B))
    diff = {p for p in set(A) & set(B) if A[p] != B[p]}
    # thermo.dat: the wall-clock timers differ on every run; every other column must not (face_gate.py cmp's rule)
    timers = {p for p in diff if p.endswith("thermo.dat")
              and ac.thermo_same(os.path.join(ac.W, "face_" + old, p), os.path.join(ac.W, "face_" + new, p))}
    diff -= timers
    bad = ["only on one side: %s" % p for p in only[:10]]
    bad += ["UNDECLARED differing path: %s" % p for p in sorted(diff - want)[:20]]
    bad += ["declared but identical: %s" % p for p in sorted(want - diff)[:20]]
    pfs = sorted({p.split("/Level_0/")[0] for p in diff if "/Level_0/" in p})
    print("field audit face_%s -> face_%s: %d / %d paths, %d differ (declared %d; + %d thermo.dat by the timers only),"
          " %d plotfiles" % (old, new, len(A), len(B), len(diff), len(want), len(timers), len(pfs)))
    seen = set()
    for r in pfs:
        a = os.path.join(ac.W, "face_" + old, r)
        b = os.path.join(ac.W, "face_" + new, r)
        _, na = ac.read_header(a)
        _, nb = ac.read_header(b)
        if na != nb:
            bad.append("%s: the component list changed" % r)
            continue
        fa, fb = fields(a), fields(b)
        moved = []
        for n in na:
            if not np.array_equal(fa[n].view(np.uint64), fb[n].view(np.uint64)):
                d = np.max(np.abs(fa[n] - fb[n]))
                moved.append("%s(max|d| %.4g; %.6g..%.6g -> %.6g..%.6g)"
                             % (n, d, fa[n].min(), fa[n].max(), fb[n].min(), fb[n].max()))
                seen.add(n)
                if n not in ALLOWED:
                    bad.append("%s: %s moved (not one of the four law-filled fields)" % (r, n))
        print("  %-36s %s" % (r, " ".join(moved) if moved else "(no field moved: Cell_H only?)"))
    print("  fields that moved, over all legs: %s (allowed %s)" % (sorted(seen), sorted(ALLOWED)))
    for b in bad[:30]:
        print("  FAIL  %s" % b)
    print("FIELD AUDIT: %s" % ("PASS" if not bad else "FAIL (%d)" % len(bad)))
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
