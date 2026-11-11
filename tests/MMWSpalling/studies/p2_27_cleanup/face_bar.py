#!/usr/bin/env python3
"""P2-27 e2 -- one face-gate bar: face_<new>.md5 vs face_<old>.md5 when NEW legs were added.

  usage: face_bar.py <old> <new> --new-legs FEET_S,FEET_0
PASS iff every path of <old> is in <new> and identical (a thermo.dat: every column but the wall-clock timers), and
every path only in <new> lies under one of the named new legs (each of which must contribute files).
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
W = os.path.join(HERE, "..", "tree_sweep_0922", "witness")
TIMERS = {"fr_t_removal", "fr_t_step", "fl_t_path"}


def load(label):
    out = {}
    for ln in open(os.path.join(W, "face_%s.md5" % label)):
        h, p = ln.rstrip("\n").split("  ", 1)
        out[p] = h
    return out


def thermo_same(a, b):
    ra = [ln.rstrip("\n").split("\t") for ln in open(a) if ln.strip()]
    rb = [ln.rstrip("\n").split("\t") for ln in open(b) if ln.strip()]
    if len(ra) != len(rb) or [h.strip() for h in ra[0]] != [h.strip() for h in rb[0]]:
        return False
    keep = [i for i, h in enumerate(ra[0]) if h.strip() not in TIMERS]
    return all([x[i] for i in keep] == [y[i] for i in keep] for x, y in zip(ra, rb))


def main():
    old, new = sys.argv[1], sys.argv[2]
    legs = sys.argv[sys.argv.index("--new-legs") + 1].split(",") if "--new-legs" in sys.argv else []
    A, B = load(old), load(new)
    gone = sorted(set(A) - set(B))
    added = sorted(set(B) - set(A))
    changed, timers = [], 0
    for p in sorted(set(A) & set(B)):
        if A[p] == B[p]:
            continue
        if p.endswith("thermo.dat") and thermo_same(os.path.join(W, "face_" + old, p), os.path.join(W, "face_" + new, p)):
            timers += 1
            continue
        changed.append(p)
    stray = [p for p in added if not any(p.startswith(l + "/") for l in legs)]
    per = {l: sum(1 for p in added if p.startswith(l + "/")) for l in legs}
    print("face bar face_%s -> face_%s: %d / %d paths; %d of the old changed (+%d thermo timers only); %d only in old;"
          " %d new (%s); %d new outside the named legs" % (old, new, len(A), len(B), len(changed), timers, len(gone),
                                                          len(added), per, len(stray)))
    for p in (changed + gone + stray)[:20]:
        print("  FAIL  %s" % p)
    ok = not changed and not gone and not stray and all(n > 0 for n in per.values())
    print("BAR: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
