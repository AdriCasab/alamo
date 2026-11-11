#!/usr/bin/env python3
"""P2-30: the two event/log writers the face gate does not reach with a data row, run on the current binary.
  V  unit/face_vapor's deck (its U3 leg)            -- the vapour row writer (regime 2); the packet's named cover (L6)
  B  unit/face_buckle's F_S + spall.buckle_log_csv=1 -- the buckle-log ROW writer (F:1011): both gate legs that write
     the file end before any line fires, so the gate hashes its header alone (found at z0; an added cover, declared)
Usage: covers.py <label>   writes S/covers/<label>/{V,B}/run*  and prints each log's row count and md5
       covers.py cmp <a> <b>   the logs byte for byte; exit 1 on a difference"""
import hashlib
import importlib.machinery
import importlib.util
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
T = os.path.join(ROOT, "tests", "MMWSpalling")
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
MPI = ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4"]
LOGS = {"V": ["run_removal_events.csv"], "B": ["run_removal_events.csv", "run_buckle_log.csv"]}


def load(path, name):
    loader = importlib.machinery.SourceFileLoader(name, path)
    spec = importlib.util.spec_from_loader(name, loader)
    m = importlib.util.module_from_spec(spec)
    loader.exec_module(m)
    return m


def run(d, deck, extra):
    if os.path.isdir(d):
        shutil.rmtree(d)
    os.makedirs(d)
    argv = MPI + [BIN, deck] + extra + ["plot_file=" + os.path.join(d, "run")]
    open(os.path.join(d, "argv.sh"), "w").write(" ".join(argv) + "\n")
    with open(os.path.join(d, "run.log"), "w") as log:
        rc = subprocess.call(argv, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    open(os.path.join(d, "rc"), "w").write("%d\n" % rc)
    return rc


def summary(label):
    for k, files in LOGS.items():
        for f in files:
            p = os.path.join(HERE, "covers", label, k, f)
            if not os.path.exists(p):
                print("  %s %-26s ABSENT" % (k, f))
                continue
            rows = open(p).read().splitlines()
            reg = {}
            if "removal_events" in f:
                for r in rows[1:]:
                    c = r.split(",")[3]
                    reg[c] = reg.get(c, 0) + 1
            print("  %s %-26s %6d data rows %s md5 %s" % (k, f, len(rows) - 1, reg or "", hashlib.md5(open(p, "rb").read()).hexdigest()))


def main():
    if sys.argv[1] == "cmp":
        bad = 0
        for k, files in LOGS.items():
            for f in files:
                a = open(os.path.join(HERE, "covers", sys.argv[2], k, f), "rb").read()
                b = open(os.path.join(HERE, "covers", sys.argv[3], k, f), "rb").read()
                same = a == b
                print("  %s %-26s %d bytes, %s" % (k, f, len(a), "byte-identical" if same else "DIFFERENT"))
                bad += 0 if same else 1
        return 1 if bad else 0
    label = sys.argv[1]
    fb = load(os.path.join(T, "unit", "face_buckle", "test"), "fb_cover")
    rcV = run(os.path.join(HERE, "covers", label, "V"), os.path.join(T, "unit", "face_vapor", "input"), [])
    rcB = run(os.path.join(HERE, "covers", label, "B"), fb.INPUT, list(fb.LEGS["F_S"]) + ["stop_time=8.0", "spall.buckle_log_csv=1"])
    print("covers %s on binary %s: V rc %d, B rc %d" % (label, hashlib.sha256(open(BIN, "rb").read()).hexdigest()[:16], rcV, rcB))
    summary(label)
    return 0 if rcV == 0 and rcB == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
