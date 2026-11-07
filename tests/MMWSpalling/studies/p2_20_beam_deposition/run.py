#!/usr/bin/env python3
"""P2-20 -- the legs. Serial, -np 4, under caffeinate; per leg legs/<leg>/{argv.sh, BINARY.txt, run.log, rc, .done}
(SR 8); an existing leg dir is archived to .old<k> first (SR 9). Every plotfile is kept.

  Z1        the edited Zhang deck (#28c schedule, implicit, dt 2 s, 150 s, omega0 0.02 as Zhang)
  Z1w       Z1 + beam.omega0=0.01414 (Oglesby's own spot: 1/e power diameter 20 mm => 1/e^2 radius 14.1 mm)
  Z2        Z1 + beam.eps_r=1.0 beam.tan_delta=<z0's string for 14 1/m at 1.070687e-2> (Zhang's recorded inputs:
            alpha 14, R 0; the nearest double -- this pair has no bit-exact string, Z0.txt)
  Z1x       Z1 + solver=explicit timestep=0.03 stop_time=50.0 amr.plot_int=333 (explicit dt witness, b3)
  Z3a/b/c   Z1 + the #19 geometry (two stacked discs: 0.10 x 0.10 x 0.05 m, 128 x 128 x 64, beam.z0 0.05), 2 kW
            for 400 s, beam.omega0 = 0.015 / 0.020 / 0.025 (the spot is not recorded: bracket), stop 400, plot every
            10 s
  CAL_S     the Kant calibration deck, stop_time=0.016 (z3: check_shared reads the rock file's two new keys)
  HOS_S     the face gate's HOS_S tokens (face_gate.from_record(HOS_R, stop 2.0)) (z3, over input_lab)

Usage: run.py LEG [LEG ...] | --argv LEG ...
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
T = os.path.join(ROOT, "tests", "MMWSpalling")
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
LEGS_DIR = os.path.join(HERE, "legs")
ZDECK = os.path.join(T, "validation", "zhang", "zhang_oglesby", "input")
CAL = os.path.join(T, "validation", "kant", "grimsel_calibration", "input")
MPI = ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4"]


def load(name, path):
    keep = sys.argv
    sys.argv = [path]
    ld = importlib.machinery.SourceFileLoader(name, path)
    m = importlib.util.module_from_spec(importlib.util.spec_from_loader(name, ld))
    ld.exec_module(m)
    sys.argv = keep
    return m


Z0 = load("z0_run220", os.path.join(HERE, "z0.py"))
Z3 = ["geometry.prob_hi=0.10 0.10 0.05", "amr.n_cell=128 128 64", "beam.z0=0.05", "beam.schedule.t_end=400.0",
      "beam.schedule.P=2000.0", "stop_time=400.0", "amr.plot_int=5"]
EXTRA = {
    "Z1": [],
    "Z1w": ["beam.omega0=0.01414"],
    "Z2": ["beam.eps_r=1.0", "beam.tan_delta=" + Z0.fixture_string(14.0, "1.070687e-2")[0]],
    "Z1x": ["solver=explicit", "timestep=0.03", "stop_time=50.0", "amr.plot_int=333"],
    "Z3a": Z3 + ["beam.omega0=0.015"],
    "Z3b": Z3 + ["beam.omega0=0.020"],
    "Z3c": Z3 + ["beam.omega0=0.025"],
}


def leg_argv(name):
    d = os.path.join(LEGS_DIR, name)
    pf = "plot_file=" + os.path.join(d, "run")
    if name in EXTRA:
        return MPI + [BIN, ZDECK] + EXTRA[name] + [pf]
    if name == "CAL_S":
        return MPI + [BIN, CAL, "stop_time=0.016", pf]
    if name == "HOS_S":
        fg = load("fg_p220", os.path.join(T, "studies", "tree_sweep_0922", "witness", "face_gate.py"))
        return fg.from_record(fg.HOS_R, d, {"stop_time": "2.0", "amr.plot_int": "125"})
    raise KeyError(name)


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
    print("  %-6s rc=%d" % (name, rc), flush=True)
    return rc


def main(args):
    if args[:1] == ["--argv"]:
        for n in args[1:]:
            print(n, " ".join(shlex.quote(x) for x in leg_argv(n)[7:]))
        return 0
    subprocess.Popen(["caffeinate", "-dis", "-w", str(os.getpid())])   # SR 11
    bad = 0
    for n in args:
        bad += run(n) != 0
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
