#!/usr/bin/env python3
"""P2-19 -- the legs. Serial, -np 4, under caffeinate; per leg legs/<leg>/{argv.sh, BINARY.txt, run.log, rc, .done}
(SR 8); an existing leg dir is archived to .old<k> first (SR 9). Every plotfile is kept (CAL's 98 files are a3's).

  T00_2 / T00_1 / T00_05   unit/face_removal's input + mesh_keys("2"/"1"/"05") + surface_patch.face_list=1 +
                           stop_time=8.0 (the face gate's T00_1_k1 form, at three meshes)
  CAL                      studies/p2_16_span_rule/legs/CAL/argv.sh verbatim, plot_file re-pointed (60 s)
  HOS                      studies/p2_16_span_rule/legs/HOS_R/argv.sh token for token, plot_file re-pointed (15 s).
                           ⚠ Since P2-17 input_lab's rock file carries spall.buckle_span = 0.012, so these tokens now
                           run with B ON; HOS_R ran B off. c2's like-for-like reference is P2-17's N_S164 (the same
                           tokens, B on, 30 s), cut by step at 15 s (reduce.py).
  C0H150                   studies/p2_14a_cref_face/legs/C0_h150/argv.sh verbatim, plot_file re-pointed (200 s)

Usage: run.py [LEG ...] (default: all, in that order) | --argv LEG ...
"""
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
FR_TEST = os.path.join(ROOT, "tests", "MMWSpalling", "unit", "face_removal", "test")
VERBATIM = {"CAL": os.path.join(STUD, "p2_16_span_rule", "legs", "CAL", "argv.sh"),
            "HOS": os.path.join(STUD, "p2_16_span_rule", "legs", "HOS_R", "argv.sh"),
            "C0H150": os.path.join(STUD, "p2_14a_cref_face", "legs", "C0_h150", "argv.sh")}
ORDER = ["T00_2", "T00_1", "T00_05", "CAL", "HOS", "C0H150"]


def _fr():
    keep = sys.argv
    sys.argv = [FR_TEST]
    ld = importlib.machinery.SourceFileLoader("fr_p219", FR_TEST)
    m = importlib.util.module_from_spec(importlib.util.spec_from_loader("fr_p219", ld))
    ld.exec_module(m)
    sys.argv = keep
    return m


def leg_argv(name):
    d = os.path.join(LEGS_DIR, name)
    if name in VERBATIM:
        toks = shlex.split(open(VERBATIM[name]).read())
        return [t for t in toks if not t.startswith("plot_file=")] + ["plot_file=" + os.path.join(d, "run")]
    fr = _fr()
    m = {"T00_2": "2", "T00_1": "1", "T00_05": "05"}[name]
    kv = dict(fr.mesh_keys(m), **{"surface_patch.face_list": "1", "stop_time": "8.0",
                                  "plot_file": os.path.join(d, "run")})
    return ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4", BIN, fr.INPUT] + \
        ["%s=%s" % (k, v) for k, v in kv.items()]


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
    a = leg_argv(name)
    open(os.path.join(d, "argv.sh"), "w").write(" ".join(shlex.quote(x) for x in a) + "\n")
    open(os.path.join(d, "BINARY.txt"), "w").write(sha(BIN) + "  bin/mmwspalling-3d-g++\n")
    with open(os.path.join(d, "run.log"), "w") as log:
        rc = subprocess.call(a, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    open(os.path.join(d, "rc"), "w").write("%d\n" % rc)
    open(os.path.join(d, ".done"), "w").write("")
    print("  %-7s rc=%d" % (name, rc), flush=True)
    return rc


def main(args):
    if args[:1] == ["--argv"]:
        for n in args[1:]:
            print(n, " ".join(shlex.quote(x) for x in leg_argv(n)[7:]))
        return 0
    subprocess.Popen(["caffeinate", "-dis", "-w", str(os.getpid())])   # SR 11
    bad = 0
    for n in (args or ORDER):
        bad += run(n) != 0
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
