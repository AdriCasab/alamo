#!/usr/bin/env python3
"""P2-30: the six tests that read the source, their source clauses only.
  locks.py capture <label>   runs each, its output to S/logs/<label>/<test>.log, exit 1 if any test exits non-zero
  locks.py cmp <a> <b>       every line of every log identical (the a2 / a5 bar); exit 1 on a difference
  locks.py bytes <pre>       a2's second bar: _src's reading of the class equals the pre-image byte for byte
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
U = os.path.join(ROOT, "tests", "MMWSpalling", "unit")
PY = "/Users/tzetze20/Desktop/code/.venv/bin/python"
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
        nl = len(r.stdout.splitlines())
        print("  %-20s rc %d, %4d lines; %s" % (name, r.returncode, nl, tail[:90]))
        bad += 1 if r.returncode else 0
    return 1 if bad else 0


def cmp_logs(a, b):
    bad, total = 0, 0
    for name, _ in SIX:
        la = open(os.path.join(HERE, "logs", a, name + ".log")).read().splitlines()
        lb = open(os.path.join(HERE, "logs", b, name + ".log")).read().splitlines()
        total += len(la)
        diff = [(i, x, y) for i, (x, y) in enumerate(zip(la, lb)) if x != y]
        ok = not diff and len(la) == len(lb)
        print("  %-20s %4d / %4d lines, %s" % (name, len(la), len(lb), "identical" if ok else "DIFFERENT"))
        for i, x, y in diff[:6]:
            print("    line %d\n      a: %s\n      b: %s" % (i + 1, x[:200], y[:200]))
        bad += 0 if ok else 1
    print("  %d lines over the six logs; %s" % (total, "every line identical" if not bad else "%d log(s) differ" % bad))
    return 1 if bad else 0


def check_bytes(pre):
    sys.path.insert(0, U)
    import _src
    d = os.path.join(HERE, pre, "Integrator")
    src = _src.tree()
    old = open(os.path.join(d, "MMWSpalling.H")).read()
    new = _src.header(src)
    print("  MMWSpalling.H in the pre-image %d bytes; _src.header() %d bytes; %s" % (
        len(old), len(new), "byte-identical" if old == new else "DIFFERENT"))
    fr = open(os.path.join(d, "MMWSpalling", "FaceRemoval.H")).read()
    sh = open(os.path.join(d, "MMWSpalling", "Shared.H")).read()
    full = old.replace('    #include "Integrator/MMWSpalling/Shared.H"\n', sh, 1).replace(
        '    #include "Integrator/MMWSpalling/FaceRemoval.H"\n', fr, 1)
    ctext = _src.text(src)
    print("  the class with every fragment in place: pre-image %d bytes; _src.text() %d bytes; %s; parent file now %d bytes" % (
        len(full), len(ctext), "byte-identical" if full == ctext else "DIFFERENT", len(src[_src.PARENT])))
    return 0 if old == new and full == ctext else 1


def main():
    if sys.argv[1] == "capture":
        return capture(sys.argv[2])
    if sys.argv[1] == "cmp":
        return cmp_logs(sys.argv[2], sys.argv[3])
    if sys.argv[1] == "bytes":
        return check_bytes(sys.argv[2])
    return 2


if __name__ == "__main__":
    sys.exit(main())
