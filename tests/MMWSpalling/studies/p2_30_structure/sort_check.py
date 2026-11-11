#!/usr/bin/env python3
"""P2-30 a4 and a5 bookkeeping.
  sort_check.py sort     a4: every candidate has a verdict and a reason; lines removed from src/ == lines added to
                         CODE_HISTORY.md (multisets, non-blank); every flagged line justified in COMMENT_FLAGS.txt
  sort_check.py trim <a> <b>  a5: the six logs of label a vs b differ only by the clause lines S/LOCK_TRIM.txt deletes, each
                         with a reason and a covering row; reports the lock lines removed
Exit 1 on a failure."""
import collections
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "tests", "MMWSpalling", "unit"))
import _src  # noqa: E402
sys.path.insert(0, HERE)
import locks  # noqa: E402
import apply_sort  # noqa: E402

HIST = os.path.join(ROOT, "docs", "project", "CODE_HISTORY.md")


def sort():
    bad = 0
    cand = apply_sort.candidates()
    tab = apply_sort.sort_table()
    miss = sorted(set(cand) - set(tab))
    empty = [k for k, v in tab.items() if not v["reason"]]
    print("  %d candidates; %d verdicts; missing %s; without a reason %s" % (len(cand), len(tab), miss, empty))
    bad += 1 if miss or empty else 0
    pre = os.path.join(HERE, "src_pre_2", "Integrator")
    names = ["MMWSpalling.H"] + [os.path.join("MMWSpalling", n) for n in sorted(os.listdir(os.path.join(pre, "MMWSpalling")))]
    old, new = collections.Counter(), collections.Counter()
    for n in names:
        old.update(s for s in open(os.path.join(pre, n)).read().split("\n") if s.strip())
        new.update(s for s in _src.read("Integrator/" + n).split("\n") if s.strip())
    removed = old - new
    added_src = new - old
    fences, inside = collections.Counter(), False
    for s in open(HIST).read().split("\n"):
        if s.startswith("```"):
            inside = not inside
        elif inside and s.strip():
            fences[s] += 1
    pointer = [s for s in added_src if "CODE_HISTORY.md" in s]
    other_added = {k: v for k, v in added_src.items() if k not in pointer}
    print("  lines removed from src/ %d; lines in CODE_HISTORY.md %d; identical as multisets: %s; pointer lines added %d; other lines added to src/ %d" % (
        sum(removed.values()), sum(fences.values()), removed == fences, len(pointer), len(other_added)))
    bad += 0 if removed == fences and not other_added and len(pointer) <= 1 else 1
    moved = [(k, n) for k, n in removed.items()]
    nc = [k for k, _ in moved if not k.lstrip().startswith("//")]
    print("  non-comment lines among those removed: %d" % len(nc))
    bad += 1 if nc else 0
    # the flag list, as apply_sort computed it at the sort's own state
    just = {}
    for ln in open(os.path.join(HERE, "COMMENT_FLAGS.txt")):
        if ln.startswith("#") or not ln.strip():
            continue
        f = ln.split(" | ")
        just[f[0].strip()] = f[-1].strip()
    moved_rows = [ln.rstrip("\n").split("\t") for ln in open(os.path.join(HERE, "COMMENT_MOVED.txt"))]
    flagged = [r for r in moved_rows if re.search(r"\bmust\b|\bnever\b|\bNOT\b|⚠", r[2])]
    print("  flagged moved lines (beside must / never / NOT / the warning sign): %d; justified lines in COMMENT_FLAGS.txt: %d" % (len(flagged), len(just)))
    for r in flagged:
        print("    %s %s  %s" % (r[0], r[1], r[2].strip()[:110]))
    bad += 0 if len(flagged) == len(just) else 1
    print("VERDICT", "PASS" if not bad else "FAIL")
    return 1 if bad else 0


def trim(a, b):
    table = []
    for ln in open(os.path.join(HERE, "LOCK_TRIM.txt")):
        if ln.startswith("#") or not ln.strip():
            continue
        f = [x.strip() for x in ln.split(" | ")]
        table.append(f)
    print("  %d clause line(s) declared deleted; each row has a reason and a covering row: %s" % (
        len(table), all(len(f) == 4 and f[2] and f[3] for f in table)))
    bad = 0 if all(len(f) == 4 and f[2] and f[3] for f in table) else 1
    declared = collections.Counter((f[0], f[1]) for f in table)
    seen = collections.Counter()
    for name, _ in locks.SIX:
        la = open(os.path.join(HERE, "logs", a, name + ".log")).read().splitlines()
        lb = open(os.path.join(HERE, "logs", b, name + ".log")).read().splitlines()
        ca, cb = collections.Counter(la), collections.Counter(lb)
        gone = ca - cb
        new = cb - ca
        for ln, n in gone.items():
            key = (name, re.sub(r"^\s*(ok|PASS|FAIL)\s+", "", ln).strip())
            seen[key] += n
        bad += 1 if new else 0
        print("  %-20s %4d -> %4d lines; removed %d, added %d" % (name, len(la), len(lb), sum(gone.values()), sum(new.values())))
        for ln in list(new)[:5]:
            print("    ADDED: %s" % ln[:160])
    undeclared = seen - declared
    unused = declared - seen
    for k in undeclared:
        print("  UNDECLARED removal: %s %s" % k)
    for k in unused:
        print("  declared but still present: %s %s" % k)
    bad += len(undeclared) + len(unused)
    print("VERDICT", "PASS" if not bad else "FAIL")
    return 1 if bad else 0


def main():
    if sys.argv[1] == "sort":
        return sort()
    if sys.argv[1] == "trim":
        return trim(sys.argv[2], sys.argv[3])
    return 2


if __name__ == "__main__":
    sys.exit(main())
