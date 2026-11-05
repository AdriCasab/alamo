#!/usr/bin/env python3
"""P2-17 -- the legs. Serial, -np 4, under caffeinate; per leg legs/<leg>/{argv.sh, ALLOW.txt, DECK.txt, BINARY.txt,
run.log, rc, .done} (SR 8); an existing leg dir is archived to .old<k> first (SR 9).

  CAL1        studies/p2_16_span_rule/legs/CAL/argv.sh verbatim (plot_file re-pointed), on the edited deck; every
              plotfile kept (P2-16 CAL's set)
  <F>_S<s>    the Hoeser face-list ladder: studies/p2_16_span_rule/legs/HOS_R/argv.sh token for token, with
              stop_time=30.0, plot_file re-pointed and surface_patch.h_expr = h_expr(S) (P2-11's run_p2_11.py,
              imported: the G1 formula with C = S x 7272.3958, %.6f) replacing the base's -- the only heat change.
              Form F: O adds spall.buckle_span=0.0 (the pre-P2-17 form), N adds nothing (production: the rock
              file's 0.012), A adds spall.buckle_span=0.0132 (the alternative span). Plotfiles at 15 s and 30 s
              kept (00000 removed).
S lives in argv only; it is never written into an input.

Usage: run.py LEG [LEG ...] | --argv LEG ... | --list
"""
import glob
import hashlib
import importlib.machinery
import importlib.util
import os
import shlex
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
STUD = os.path.join(ROOT, "tests", "MMWSpalling", "studies")
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
LEGS_DIR = os.path.join(HERE, "legs")
P16 = os.path.join(STUD, "p2_16_span_rule", "legs")
SPAN_KEY = "spall.buckle_span"
FORMS = {"O": ["%s=0.0" % SPAN_KEY], "N": [], "A": ["%s=0.0132" % SPAN_KEY]}
S_GRID = {"O": (0.20, 0.1636, 0.14, 0.12, 0.10, 0.087, 0.07),
          "N": (0.20, 0.1636, 0.14, 0.12, 0.10, 0.087, 0.07),
          "A": (0.14, 0.12, 0.10, 0.087, 0.07)}


def _load(name, path):
    sys.argv, keep = [path], sys.argv
    ld = importlib.machinery.SourceFileLoader(name, path)
    m = importlib.util.module_from_spec(importlib.util.spec_from_loader(name, ld))
    ld.exec_module(m)
    sys.argv = keep
    return m


R11 = _load("run11_p217", os.path.join(STUD, "p2_11_heat", "run_p2_11.py"))


def name_of(form, S):
    return "%s_S%03d" % (form, round(1000 * S))


def ladder():
    return {name_of(f, S): (f, S) for f in ("O", "N", "A") for S in S_GRID[f]}


def all_legs():
    return ["CAL1"] + list(ladder())


def leg_argv(name):
    d = os.path.join(LEGS_DIR, name)
    if name == "CAL1":
        toks = shlex.split(open(os.path.join(P16, "CAL", "argv.sh")).read())
        out = [t for t in toks if not t.startswith("plot_file=")] + ["plot_file=" + os.path.join(d, "run")]
        return out, out[out.index(BIN) + 1] if BIN in out else out[7], []
    form, S = ladder()[name]
    toks = shlex.split(open(os.path.join(P16, "HOS_R", "argv.sh")).read())
    rep = {"stop_time": "30.0", "plot_file": os.path.join(d, "run"), "surface_patch.h_expr": R11.h_expr(S)}
    out, seen = [], set()
    for t in toks:
        k = t.split("=", 1)[0] if "=" in t else None
        if k in rep:
            out.append("%s=%s" % (k, rep[k]))
            seen.add(k)
        else:
            out.append(t)
    assert seen == set(rep), "the base argv lacks %s" % (set(rep) - seen)
    out += FORMS[form]
    deck = out[out.index(BIN) + 1]
    return out, deck, ([SPAN_KEY] if FORMS[form] else [])


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def run(name):
    d = os.path.join(LEGS_DIR, name)
    if os.path.isdir(d):                      # SR 9
        k = 1
        while os.path.isdir(d + ".old%d" % k):
            k += 1
        shutil.move(d, d + ".old%d" % k)
    os.makedirs(d)
    a, deck, allow = leg_argv(name)
    open(os.path.join(d, "argv.sh"), "w").write(" ".join(shlex.quote(x) for x in a) + "\n")
    open(os.path.join(d, "ALLOW.txt"), "w").write(" ".join(allow) + "\n")
    open(os.path.join(d, "DECK.txt"), "w").write(deck + "\n")
    open(os.path.join(d, "BINARY.txt"), "w").write(sha(BIN) + "  bin/mmwspalling-3d-g++\n")
    with open(os.path.join(d, "run.log"), "w") as log:
        rc = subprocess.call(a, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    open(os.path.join(d, "rc"), "w").write("%d\n" % rc)
    open(os.path.join(d, ".done"), "w").write("")
    if name != "CAL1":                        # keep 15 s and 30 s; drop the t = 0 plotfile
        for p in glob.glob(os.path.join(d, "run", "00000" + "*")):
            shutil.rmtree(p)
    print("  %-7s rc=%d" % (name, rc), flush=True)
    return rc


def main(args):
    if args[:1] == ["--list"]:
        print(" ".join(all_legs()))
        return 0
    if args[:1] == ["--argv"]:
        for n in args[1:]:
            a, deck, allow = leg_argv(n)
            print(n, " ".join(shlex.quote(x) for x in a[7:]), "| allow:", allow)
        return 0
    subprocess.Popen(["caffeinate", "-dis", "-w", str(os.getpid())])   # SR 11
    bad = 0
    for n in args:
        bad += run(n) != 0
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
