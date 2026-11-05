#!/usr/bin/env python3
"""P2-16 -- the legs. Serial, -np 4, under caffeinate; per leg legs/<leg>/{argv.sh, BINARY.txt, ALLOW.txt, DECK.txt,
run.log, rc, .done} (SR 8); an existing leg dir is archived to .old<k> first (SR 9). Plotfiles: only the last is kept
except CAL and HOS_R (packet guardrail). P2-15's leg helpers (kt, BT, TKEYS, hos_tokens, sha, prune) are imported.

  retirement checks (a3): CAL, HOS_R (P2-13 argv verbatim), F0_k1 (unit/face_removal B1 + face_list=1), KB150_0
                          (P2-15 argv verbatim)
  sensitivity (c1):       KB150_m30, KB160_m30 (P2-15 KB argv + stress_relief=R_30 + stop_time=30.0), CAL_m30
                          (CAL + R_30); R_30 is read from R_BETHEL.txt at launch (bethel_R.py)

Usage: run.py LEG [LEG ...] | --argv LEG ... | UNITS [name ...]
"""
import importlib.machinery
import importlib.util
import os
import shlex
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
T = os.path.join(ROOT, "tests", "MMWSpalling")


def load(name, path):
    ld = importlib.machinery.SourceFileLoader(name, path)
    m = importlib.util.module_from_spec(importlib.util.spec_from_loader(name, ld))
    ld.exec_module(m)
    return m


R15 = load("run15", os.path.join(T, "studies", "p2_15_buckle", "run.py"))
BIN, CAL_IN, MPI = R15.BIN, R15.CAL_IN, R15.MPI
LEGS_DIR = os.path.join(HERE, "legs")
FR_IN = os.path.join(T, "unit", "face_removal", "input")
B1 = ["amr.n_cell=24 24 64", "timestep=0.008", "amr.max_grid_size=16", "amr.plot_int=125"]
KEEP_ALL_PF = {"CAL", "HOS_R"}
RKEY = "spallation.sp.stress_relief"


def r30():
    return open(os.path.join(HERE, "R_BETHEL.txt")).read().split()[0]


def legs():
    """name -> (deck kind, tokens, allowed shared keys, needs R_30)"""
    return {
        "CAL": ("CAL", [], [], False),
        "HOS_R": ("HOS", [], [], False),
        "F0_k1": ("FR", B1 + ["surface_patch.face_list=1"], [], False),
        "KB150_0": ("CAL", R15.kt(150) + R15.BT + ["stop_time=200.0"], R15.TKEYS, False),
        "KB150_m30": ("CAL", R15.kt(150) + R15.BT + ["stop_time=30.0"], R15.TKEYS, True),
        "KB160_m30": ("CAL", R15.kt(160) + R15.BT + ["stop_time=30.0"], R15.TKEYS, True),
        "CAL_m30": ("CAL", [], [], True),
    }


def leg_argv(name):
    kind, toks, allow, rr = legs()[name]
    d = os.path.join(LEGS_DIR, name)
    if kind == "CAL":
        deck, base = CAL_IN, []
    elif kind == "HOS":
        deck, base = R15.hos_tokens()
    else:
        deck, base = FR_IN, []
    toks, allow = base + list(toks), list(allow)
    if rr:
        toks.append("%s=%s" % (RKEY, r30()))
        allow.append(RKEY)
    return MPI + [BIN, deck] + toks + ["plot_file=" + os.path.join(d, "run")], deck, allow


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
    open(os.path.join(d, "BINARY.txt"), "w").write(R15.sha(BIN) + "  bin/mmwspalling-3d-g++\n")
    with open(os.path.join(d, "run.log"), "w") as log:
        rc = subprocess.call(a, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    open(os.path.join(d, "rc"), "w").write("%d\n" % rc)
    open(os.path.join(d, ".done"), "w").write("")
    if name not in KEEP_ALL_PF:
        R15.prune(d)
    print("  %-10s rc=%d" % (name, rc), flush=True)
    return rc


def UNITS(only=None):
    R10 = load("r10_p216", os.path.join(T, "studies", "p2_10_retire", "run_p2_10.py"))
    R10.HERE = HERE
    R10.units("p2_16", only)
    return 0


def main(args):
    if args[:1] == ["--argv"]:
        for n in args[1:]:
            a, deck, allow = leg_argv(n)
            print(n, " ".join(shlex.quote(x) for x in a[7:]), "| allow:", allow)
        return 0
    subprocess.Popen(["caffeinate", "-dis", "-w", str(os.getpid())])   # SR 11
    if args[:1] == ["UNITS"]:
        return UNITS(args[1:] or None)
    bad = 0
    for n in args:
        bad += run(n) != 0
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
