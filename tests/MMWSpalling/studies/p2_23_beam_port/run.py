#!/usr/bin/env python3
"""P2-23 -- the legs. Serial, -np 4, under caffeinate; per leg legs/<leg>/{argv.sh, BINARY.txt, run.log, rc, .done}
(SR 8); an existing leg dir is archived to .old<k> first (SR 9).

D = validation/oglesby/input (the production rock and form, the patch OFF, the beam on #19 at 12.88 mm, the measured
ramp, B on, 1.5625 mm, dt 0.2 s, 80 s).
  B0    D as written
  P     D + the POSITIVE CONTROL: beam.eps_r=1.0 beam.tan_delta=0.6970986507425017 beam.lambda=2.19e-3 (alpha 2000 1/m
        exactly in the code's association, R 0) "beam.tan_delta_T=295.15 1.0" (flat) beam.schedule.t_end=20.0
        beam.schedule.P=150.0 timestep=0.016 stop_time=20.0 amr.plot_int=250
  B0F   D + "beam.tan_delta_T=295.15 1.0 1268.15 1.0" (two equal pairs: flat alpha 14 on the sum path)
  B1    D + amr.n_cell=128 128 64 amr.max_grid_size=32 timestep=0.05 amr.plot_int=400 amr.thermo.plot_int=10
        + weibull.V0=4.76837158203125e-10 (DECLARED: the cell volume at 0.78 mm, check_shared (iv)'s rule -- the
        packet's B1 argv omits it, which would rescale the flaw per cell by (V/V0)^(1/m))
  SMOKE_<leg>  the leg on 32 x 32 x 16 (V0 3.0517578125e-8, dt 0.5 s; P dt 0.1 s) for 10 s, the witnesses only
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
    "B0": [],
    "P": ["beam.eps_r=1.0", "beam.tan_delta=0.6970986507425017", "beam.lambda=2.19e-3", "beam.tan_delta_T=295.15 1.0",
          "beam.schedule.t_end=20.0", "beam.schedule.P=150.0", "timestep=0.016", "stop_time=20.0",
          "amr.plot_int=250"],
    # ADDENDUM (advisor, the user's request 2026-11-06, after the implementer's takeaway 1): the robust positive
    # control at 400 W (2.7x P's flux: 1.54 MW/m2 peak absorbed, ~460 K/s in the top cell), dt 8 ms, 20 s
    "P400": ["beam.eps_r=1.0", "beam.tan_delta=0.6970986507425017", "beam.lambda=2.19e-3", "beam.tan_delta_T=295.15 1.0",
             "beam.schedule.t_end=20.0", "beam.schedule.P=400.0", "timestep=0.008", "stop_time=20.0",
             "amr.plot_int=500"],
    "B0F": ["beam.tan_delta_T=295.15 1.0 1268.15 1.0"],
    "B1": ["amr.n_cell=128 128 64", "amr.max_grid_size=32", "timestep=0.05", "amr.plot_int=400",
           "amr.thermo.plot_int=10", "weibull.V0=4.76837158203125e-10"],
}
SMOKE = ["amr.n_cell=32 32 16", "amr.max_grid_size=16", "weibull.V0=3.0517578125e-8", "stop_time=10.0",
         "amr.plot_int=100000"]


def extra(name):
    if name.startswith("SMOKE_"):
        base = name[6:]
        e = list(EXTRA[base]) + SMOKE
        return e + (["timestep=0.1"] if base == "P" else ["timestep=0.5"])
    return EXTRA[name]


def leg_argv(name):
    return MPI + [BIN, D] + extra(name) + ["plot_file=" + os.path.join(LEGS_DIR, name, "run")]


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
