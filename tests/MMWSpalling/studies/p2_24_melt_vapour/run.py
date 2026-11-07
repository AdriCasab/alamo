#!/usr/bin/env python3
"""P2-24 -- the legs. Serial, -np 4, under caffeinate; per leg legs/<leg>/{argv.sh, BINARY.txt, run.log, rc, .done}
(SR 8); an existing leg dir is archived to .old<k> first (SR 9).

D = validation/oglesby/input as edited by P2-24: the production rock WITH the melt set (Zhang's, LIT), the beam's melt
keys (x4, R_melt 0.3), the vapour rule on, Fig. 34's whole trace (158 x 10 s from beam-on, then 0), 1.5625 mm, dt 0.2 s,
1,800 s (1,580 s of heating, 220 s of cooling).
  M0    D as written
  M100  D + beam.melt_factor=100.0 (the DOE model's "two orders of magnitude": the recorded alternative)
  M0F   D + "beam.tan_delta_T=295.15 1.0 1268.15 1.0" (two equal pairs: a flat alpha 14 on the sum path, x4 on top)
  M1    D + amr.n_cell=128 128 64 amr.max_grid_size=32 timestep=0.05 amr.plot_int=400 amr.thermo.plot_int=10
        weibull.V0=4.76837158203125e-10 stop_time=170.0 (0.78 mm through the first ~120 s after onset)
  SMOKE_M0  M0 on 32 x 32 x 16, dt 0.5 s, the whole 1,800 s (the witnesses and the test's code paths only)
Usage: run.py LEG [LEG ...] | --argv LEG ...
"""
import hashlib
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
D = os.path.join(T, "validation", "oglesby", "input")
MPI = ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4"]
EXTRA = {
    "M0": [],
    "M100": ["beam.melt_factor=100.0"],
    "M0F": ["beam.tan_delta_T=295.15 1.0 1268.15 1.0"],
    "M1": ["amr.n_cell=128 128 64", "amr.max_grid_size=32", "timestep=0.05", "amr.plot_int=400",
           "amr.thermo.plot_int=10", "weibull.V0=4.76837158203125e-10", "stop_time=170.0"],
    # ADDENDUM (advisor, the user's request 2026-11-07, after the review): the melt heat-transport SENSITIVITY, two
    # declared alternatives on M0, one key each, scored on the same observables (no number adjusted); run on the
    # post-P2-24 binary that adds grimsel.lambda.cap_theta (default -1 = off), with M0R = M0 re-run on that binary as
    # the identity witness (thermo.dat byte for byte against legs/M0, the P2-24 binary 98a261d3 in bin_p2_24/).
    #   LAM  lambda(T) held at its 1000 degC (1273.15 K) value 2000/1563 = 1.280 W/mK above it (Zhang's cap, the
    #        rock file's recorded alternative), instead of the fit's 0.60 W/mK at 3000 K
    #   EPS  the surface emissivity 0.5 at every temperature (the rock file's 0.8); both keys so no step-down
    "M0R": [],
    "LAM": ["microstructure.phase0.grimsel.lambda.cap_theta=1000.0"],
    "EPS": ["losses.epsilon=0.5", "losses.epsilon_high=0.5"],
    # the smoke before the campaign: M0 on 32 x 32 x 16 (3.125 mm; the limit 0.97 s), the whole 1,800 s
    "SMOKE_M0": ["amr.n_cell=32 32 16", "amr.max_grid_size=16", "weibull.V0=3.0517578125e-8", "timestep=0.5",
                 "amr.plot_int=100000"],
}


def leg_argv(name):
    return MPI + [BIN, D] + EXTRA[name] + ["plot_file=" + os.path.join(LEGS_DIR, name, "run")]


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
    print("  %-10s rc=%d" % (name, rc), flush=True)
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
