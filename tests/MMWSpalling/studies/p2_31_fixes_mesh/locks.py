#!/usr/bin/env python3
"""P2-31: the six tests that read the source, their source clauses only (P2-30's locks.py pattern), plus the binary
copy and the pre-images every code stage takes before it edits.
  locks.py capture <label>        each test's output to S/logs/<label>/<test>.log; exit 1 if any test exits non-zero
  locks.py cmp <a> <b> [allow]    line multisets per log: lines only in <b> must match the regex `allow` (the stage's
                                  declared new clause lines), lines only in <a> are reported and fail the row
  locks.py preimage <stage> <rel>...   copy each repo-relative file to S/src_pre_<stage>/<rel> (and the binary to
                                  S/bins/pre_<stage>), refusing to overwrite an existing pre-image (SR 9)
"""
import collections
import hashlib
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
U = os.path.join(ROOT, "tests", "MMWSpalling", "unit")
PY = "/Users/tzetze20/Desktop/code/.venv/bin/python"
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
SIX = [("retired_keys", ["--source"]), ("lateral_flux_lock", ["--check"]), ("buckle_lock", []), ("profile_depth", []),
       ("face_heat", ["--source"]), ("face_removal", ["--source"])]


def capture(label):
    d = os.path.join(HERE, "logs", label)
    os.makedirs(d, exist_ok=True)
    bad = 0
    for name, args in SIX:
        r = subprocess.run([PY, os.path.join(U, name, "test")] + args, cwd=ROOT, capture_output=True, text=True)
        open(os.path.join(d, name + ".log"), "w").write(r.stdout + r.stderr)
        tail = (r.stdout.strip().splitlines() or ["(no output)"])[-1]
        print("  %-20s rc %d, %4d lines; %s" % (name, r.returncode, len(r.stdout.splitlines()), tail[:90]))
        bad += 1 if r.returncode else 0
    return 1 if bad else 0


def cmp_logs(a, b, allow=None):
    bad = 0
    for name, _ in SIX:
        la = open(os.path.join(HERE, "logs", a, name + ".log")).read().splitlines()
        lb = open(os.path.join(HERE, "logs", b, name + ".log")).read().splitlines()
        ca, cb = collections.Counter(la), collections.Counter(lb)
        gone, new = ca - cb, cb - ca
        odd = [ln for ln in new if not (allow and re.search(allow, ln))]
        print("  %-20s %4d -> %4d lines; only in %s %d, only in %s %d (declared %d, undeclared %d)" % (
            name, len(la), len(lb), a, sum(gone.values()), b, sum(new.values()), sum(new.values()) - len(odd), len(odd)))
        for ln in list(gone)[:6]:
            print("      - %s" % ln[:170])
        for ln in list(new)[:6]:
            print("      + %s" % ln[:170])
        bad += 1 if (gone or odd) else 0
    print("  %s" % ("every lock prints what it printed, plus only the declared lines" if not bad
                    else "%d lock log(s) changed outside the declaration" % bad))
    return 1 if bad else 0


def preimage(stage, rels):
    d = os.path.join(HERE, "src_pre_%s" % stage)
    for rel in rels:
        dst = os.path.join(d, rel)
        if os.path.exists(dst):
            print("  exists, kept: %s" % os.path.relpath(dst, ROOT))
            continue
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copy2(os.path.join(ROOT, rel), dst)
    bd = os.path.join(HERE, "bins")
    os.makedirs(bd, exist_ok=True)
    gi = os.path.join(bd, ".gitignore")
    if not os.path.exists(gi):
        open(gi, "w").write("# the pre-stage binaries the two-sided clauses run (3.5 MB each): local only\n*\n!.gitignore\n")
    b = os.path.join(bd, "pre_%s" % stage)
    if not os.path.exists(b):
        shutil.copy2(BIN, b)
    h = hashlib.sha256(open(b, "rb").read()).hexdigest()
    open(os.path.join(HERE, "BINARY_pre_%s.txt" % stage), "w").write("%s  bins/pre_%s\n" % (h, stage))
    print("  pre-image %s: %d file(s); binary bins/pre_%s sha256 %s" % (stage, len(rels), stage, h[:16]))
    return 0


def main():
    if sys.argv[1] == "capture":
        return capture(sys.argv[2])
    if sys.argv[1] == "cmp":
        return cmp_logs(sys.argv[2], sys.argv[3], sys.argv[4] if len(sys.argv) > 4 else None)
    if sys.argv[1] == "preimage":
        return preimage(sys.argv[2], sys.argv[3:])
    return 2


if __name__ == "__main__":
    sys.exit(main())
