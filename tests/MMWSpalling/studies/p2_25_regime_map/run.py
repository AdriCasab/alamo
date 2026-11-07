#!/usr/bin/env python3
"""P2-25 -- the legs. Serial, -np 4, under caffeinate; per leg legs/<leg>/{argv.sh, BINARY.txt, run.log, rc, .done}
(SR 8); an existing leg dir is archived to .old<k> first (SR 9).

D = studies/p2_25_regime_map/input (the #19 deck, flat loss on the sum path, one-window schedule); D19 =
validation/oglesby/input (P400R only). The argv below is LITERAL: the map's stops come from z0's 1-D column
(STOPS.json), pasted here once; `run.py --check` rebuilds the table from STOPS.json and the packet's rules and asserts it
equals this one.
  D<d>_Q<q> (42)  beam.tan_delta=<td> (none at 71 mm) beam.schedule.P=<P> beam.schedule.t_end=<stop> stop_time=<stop>
                  timestep=<dt>
  P400R           D19 + P2-23's P400 argv with the one-pair table made two pairs (the melt keys abort on the scalar path)
  K05_/K2_D<d>_Q2 (14)  D<d>_Q2 + spallation.sp.kic.scale=0.5 / 2.0 (declared brackets, SR 1)
  W6_/W19_D<1|4>_Q2 (4) D<d>_Q2 with beam.omega0 6.44 / 19.32 mm and P x (w/12.88 mm)^2 (q held)
  F_D<1|2|4>_Q2, F_D2_Q1, F_D2_Q4 (5)  the 0.78 mm twin: amr.n_cell=128 128 64 amr.max_grid_size=32 weibull.V0 =
                  (0.78125 mm)^3, dt from the rule at dz/2, stop 20 s (60 s at q 1)
Usage: run.py LEG [LEG ...] | --argv LEG ... | --check | --list
"""
import hashlib
import json
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
D = os.path.join(HERE, "input")
D19 = os.path.join(T, "validation", "oglesby", "input")
MPI = ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4"]
P400R = ["beam.eps_r=1.0", "beam.tan_delta=0.6970986507425017", "beam.lambda=2.19e-3",
         "beam.tan_delta_T=295.15 1.0 1268.15 1.0", "beam.schedule.t_end=20.0", "beam.schedule.P=400.0",
         "timestep=0.008", "stop_time=20.0", "amr.plot_int=500"]
EXTRA = {
    'D0p5_Q0p3': ['beam.tan_delta=1.52111', 'beam.schedule.P=91.6', 'beam.schedule.t_end=77.0', 'stop_time=77.0', 'timestep=0.06'],
    'D1_Q0p3': ['beam.tan_delta=0.760555', 'beam.schedule.P=91.6', 'beam.schedule.t_end=92.0', 'stop_time=92.0', 'timestep=0.08'],
    'D2_Q0p3': ['beam.tan_delta=0.380278', 'beam.schedule.P=91.6', 'beam.schedule.t_end=121.0', 'stop_time=121.0', 'timestep=0.1'],
    'D4_Q0p3': ['beam.tan_delta=0.190139', 'beam.schedule.P=91.6', 'beam.schedule.t_end=180.0', 'stop_time=180.0', 'timestep=0.1'],
    'D8_Q0p3': ['beam.tan_delta=0.0950694', 'beam.schedule.P=91.6', 'beam.schedule.t_end=299.0', 'stop_time=299.0', 'timestep=0.2'],
    'D16_Q0p3': ['beam.tan_delta=0.0475347', 'beam.schedule.P=91.6', 'beam.schedule.t_end=554.0', 'stop_time=554.0', 'timestep=0.2'],
    'D71_Q0p3': ['beam.schedule.P=91.6', 'beam.schedule.t_end=600.0', 'stop_time=600.0', 'timestep=0.2'],
    'D0p5_Q0p5': ['beam.tan_delta=1.52111', 'beam.schedule.P=152.7', 'beam.schedule.t_end=60.0', 'stop_time=60.0', 'timestep=0.03'],
    'D1_Q0p5': ['beam.tan_delta=0.760555', 'beam.schedule.P=152.7', 'beam.schedule.t_end=60.0', 'stop_time=60.0', 'timestep=0.04'],
    'D2_Q0p5': ['beam.tan_delta=0.380278', 'beam.schedule.P=152.7', 'beam.schedule.t_end=60.0', 'stop_time=60.0', 'timestep=0.07'],
    'D4_Q0p5': ['beam.tan_delta=0.190139', 'beam.schedule.P=152.7', 'beam.schedule.t_end=80.0', 'stop_time=80.0', 'timestep=0.1'],
    'D8_Q0p5': ['beam.tan_delta=0.0950694', 'beam.schedule.P=152.7', 'beam.schedule.t_end=139.0', 'stop_time=139.0', 'timestep=0.2'],
    'D16_Q0p5': ['beam.tan_delta=0.0475347', 'beam.schedule.P=152.7', 'beam.schedule.t_end=263.0', 'stop_time=263.0', 'timestep=0.2'],
    'D71_Q0p5': ['beam.schedule.P=152.7', 'beam.schedule.t_end=600.0', 'stop_time=600.0', 'timestep=0.2'],
    'D0p5_Q1': ['beam.tan_delta=1.52111', 'beam.schedule.P=305.3', 'beam.schedule.t_end=60.0', 'stop_time=60.0', 'timestep=0.01'],
    'D1_Q1': ['beam.tan_delta=0.760555', 'beam.schedule.P=305.3', 'beam.schedule.t_end=60.0', 'stop_time=60.0', 'timestep=0.02'],
    'D2_Q1': ['beam.tan_delta=0.380278', 'beam.schedule.P=305.3', 'beam.schedule.t_end=60.0', 'stop_time=60.0', 'timestep=0.03'],
    'D4_Q1': ['beam.tan_delta=0.190139', 'beam.schedule.P=305.3', 'beam.schedule.t_end=60.0', 'stop_time=60.0', 'timestep=0.05'],
    'D8_Q1': ['beam.tan_delta=0.0950694', 'beam.schedule.P=305.3', 'beam.schedule.t_end=60.0', 'stop_time=60.0', 'timestep=0.1'],
    'D16_Q1': ['beam.tan_delta=0.0475347', 'beam.schedule.P=305.3', 'beam.schedule.t_end=107.0', 'stop_time=107.0', 'timestep=0.2'],
    'D71_Q1': ['beam.schedule.P=305.3', 'beam.schedule.t_end=525.0', 'stop_time=525.0', 'timestep=0.2'],
    'D0p5_Q2': ['beam.tan_delta=1.52111', 'beam.schedule.P=610.7', 'beam.schedule.t_end=60.0', 'stop_time=60.0', 'timestep=0.009'],
    'D1_Q2': ['beam.tan_delta=0.760555', 'beam.schedule.P=610.7', 'beam.schedule.t_end=60.0', 'stop_time=60.0', 'timestep=0.01'],
    'D2_Q2': ['beam.tan_delta=0.380278', 'beam.schedule.P=610.7', 'beam.schedule.t_end=60.0', 'stop_time=60.0', 'timestep=0.01'],
    'D4_Q2': ['beam.tan_delta=0.190139', 'beam.schedule.P=610.7', 'beam.schedule.t_end=60.0', 'stop_time=60.0', 'timestep=0.02'],
    'D8_Q2': ['beam.tan_delta=0.0950694', 'beam.schedule.P=610.7', 'beam.schedule.t_end=60.0', 'stop_time=60.0', 'timestep=0.05'],
    'D16_Q2': ['beam.tan_delta=0.0475347', 'beam.schedule.P=610.7', 'beam.schedule.t_end=60.0', 'stop_time=60.0', 'timestep=0.1'],
    'D71_Q2': ['beam.schedule.P=610.7', 'beam.schedule.t_end=223.0', 'stop_time=223.0', 'timestep=0.2'],
    'D0p5_Q4': ['beam.tan_delta=1.52111', 'beam.schedule.P=1221.3', 'beam.schedule.t_end=20.0', 'stop_time=20.0', 'timestep=0.004'],
    'D1_Q4': ['beam.tan_delta=0.760555', 'beam.schedule.P=1221.3', 'beam.schedule.t_end=20.0', 'stop_time=20.0', 'timestep=0.006'],
    'D2_Q4': ['beam.tan_delta=0.380278', 'beam.schedule.P=1221.3', 'beam.schedule.t_end=20.0', 'stop_time=20.0', 'timestep=0.008'],
    'D4_Q4': ['beam.tan_delta=0.190139', 'beam.schedule.P=1221.3', 'beam.schedule.t_end=20.0', 'stop_time=20.0', 'timestep=0.01'],
    'D8_Q4': ['beam.tan_delta=0.0950694', 'beam.schedule.P=1221.3', 'beam.schedule.t_end=20.0', 'stop_time=20.0', 'timestep=0.02'],
    'D16_Q4': ['beam.tan_delta=0.0475347', 'beam.schedule.P=1221.3', 'beam.schedule.t_end=22.0', 'stop_time=22.0', 'timestep=0.05'],
    'D71_Q4': ['beam.schedule.P=1221.3', 'beam.schedule.t_end=101.0', 'stop_time=101.0', 'timestep=0.2'],
    'D0p5_Q8': ['beam.tan_delta=1.52111', 'beam.schedule.P=2442.7', 'beam.schedule.t_end=20.0', 'stop_time=20.0', 'timestep=0.002'],
    'D1_Q8': ['beam.tan_delta=0.760555', 'beam.schedule.P=2442.7', 'beam.schedule.t_end=20.0', 'stop_time=20.0', 'timestep=0.003'],
    'D2_Q8': ['beam.tan_delta=0.380278', 'beam.schedule.P=2442.7', 'beam.schedule.t_end=20.0', 'stop_time=20.0', 'timestep=0.004'],
    'D4_Q8': ['beam.tan_delta=0.190139', 'beam.schedule.P=2442.7', 'beam.schedule.t_end=20.0', 'stop_time=20.0', 'timestep=0.007'],
    'D8_Q8': ['beam.tan_delta=0.0950694', 'beam.schedule.P=2442.7', 'beam.schedule.t_end=20.0', 'stop_time=20.0', 'timestep=0.01'],
    'D16_Q8': ['beam.tan_delta=0.0475347', 'beam.schedule.P=2442.7', 'beam.schedule.t_end=20.0', 'stop_time=20.0', 'timestep=0.02'],
    'D71_Q8': ['beam.schedule.P=2442.7', 'beam.schedule.t_end=47.0', 'stop_time=47.0', 'timestep=0.1'],
    'K05_D0p5_Q2': ['beam.tan_delta=1.52111', 'beam.schedule.P=610.7', 'beam.schedule.t_end=60.0', 'stop_time=60.0', 'timestep=0.009', 'spallation.sp.kic.scale=0.5'],
    'K05_D1_Q2': ['beam.tan_delta=0.760555', 'beam.schedule.P=610.7', 'beam.schedule.t_end=60.0', 'stop_time=60.0', 'timestep=0.01', 'spallation.sp.kic.scale=0.5'],
    'K05_D2_Q2': ['beam.tan_delta=0.380278', 'beam.schedule.P=610.7', 'beam.schedule.t_end=60.0', 'stop_time=60.0', 'timestep=0.01', 'spallation.sp.kic.scale=0.5'],
    'K05_D4_Q2': ['beam.tan_delta=0.190139', 'beam.schedule.P=610.7', 'beam.schedule.t_end=60.0', 'stop_time=60.0', 'timestep=0.02', 'spallation.sp.kic.scale=0.5'],
    'K05_D8_Q2': ['beam.tan_delta=0.0950694', 'beam.schedule.P=610.7', 'beam.schedule.t_end=60.0', 'stop_time=60.0', 'timestep=0.05', 'spallation.sp.kic.scale=0.5'],
    'K05_D16_Q2': ['beam.tan_delta=0.0475347', 'beam.schedule.P=610.7', 'beam.schedule.t_end=60.0', 'stop_time=60.0', 'timestep=0.1', 'spallation.sp.kic.scale=0.5'],
    'K05_D71_Q2': ['beam.schedule.P=610.7', 'beam.schedule.t_end=223.0', 'stop_time=223.0', 'timestep=0.2', 'spallation.sp.kic.scale=0.5'],
    'K2_D0p5_Q2': ['beam.tan_delta=1.52111', 'beam.schedule.P=610.7', 'beam.schedule.t_end=60.0', 'stop_time=60.0', 'timestep=0.009', 'spallation.sp.kic.scale=2.0'],
    'K2_D1_Q2': ['beam.tan_delta=0.760555', 'beam.schedule.P=610.7', 'beam.schedule.t_end=60.0', 'stop_time=60.0', 'timestep=0.01', 'spallation.sp.kic.scale=2.0'],
    'K2_D2_Q2': ['beam.tan_delta=0.380278', 'beam.schedule.P=610.7', 'beam.schedule.t_end=60.0', 'stop_time=60.0', 'timestep=0.01', 'spallation.sp.kic.scale=2.0'],
    'K2_D4_Q2': ['beam.tan_delta=0.190139', 'beam.schedule.P=610.7', 'beam.schedule.t_end=60.0', 'stop_time=60.0', 'timestep=0.02', 'spallation.sp.kic.scale=2.0'],
    'K2_D8_Q2': ['beam.tan_delta=0.0950694', 'beam.schedule.P=610.7', 'beam.schedule.t_end=60.0', 'stop_time=60.0', 'timestep=0.05', 'spallation.sp.kic.scale=2.0'],
    'K2_D16_Q2': ['beam.tan_delta=0.0475347', 'beam.schedule.P=610.7', 'beam.schedule.t_end=60.0', 'stop_time=60.0', 'timestep=0.1', 'spallation.sp.kic.scale=2.0'],
    'K2_D71_Q2': ['beam.schedule.P=610.7', 'beam.schedule.t_end=223.0', 'stop_time=223.0', 'timestep=0.2', 'spallation.sp.kic.scale=2.0'],
    'W6_D1_Q2': ['beam.tan_delta=0.760555', 'beam.schedule.P=152.7', 'beam.schedule.t_end=60.0', 'stop_time=60.0', 'timestep=0.01', 'beam.omega0=0.00644'],
    'W6_D4_Q2': ['beam.tan_delta=0.190139', 'beam.schedule.P=152.7', 'beam.schedule.t_end=60.0', 'stop_time=60.0', 'timestep=0.02', 'beam.omega0=0.00644'],
    'W19_D1_Q2': ['beam.tan_delta=0.760555', 'beam.schedule.P=1374.0', 'beam.schedule.t_end=60.0', 'stop_time=60.0', 'timestep=0.01', 'beam.omega0=0.01932'],
    'W19_D4_Q2': ['beam.tan_delta=0.190139', 'beam.schedule.P=1374.0', 'beam.schedule.t_end=60.0', 'stop_time=60.0', 'timestep=0.02', 'beam.omega0=0.01932'],
    'F_D1_Q2': ['beam.tan_delta=0.760555', 'beam.schedule.P=610.7', 'beam.schedule.t_end=20.0', 'stop_time=20.0', 'timestep=0.008', 'amr.n_cell=128 128 64', 'amr.max_grid_size=32', 'weibull.V0=4.76837158203125e-10'],
    'F_D2_Q2': ['beam.tan_delta=0.380278', 'beam.schedule.P=610.7', 'beam.schedule.t_end=20.0', 'stop_time=20.0', 'timestep=0.01', 'amr.n_cell=128 128 64', 'amr.max_grid_size=32', 'weibull.V0=4.76837158203125e-10'],
    'F_D4_Q2': ['beam.tan_delta=0.190139', 'beam.schedule.P=610.7', 'beam.schedule.t_end=20.0', 'stop_time=20.0', 'timestep=0.02', 'amr.n_cell=128 128 64', 'amr.max_grid_size=32', 'weibull.V0=4.76837158203125e-10'],
    'F_D2_Q1': ['beam.tan_delta=0.380278', 'beam.schedule.P=305.3', 'beam.schedule.t_end=60.0', 'stop_time=60.0', 'timestep=0.02', 'amr.n_cell=128 128 64', 'amr.max_grid_size=32', 'weibull.V0=4.76837158203125e-10'],
    'F_D2_Q4': ['beam.tan_delta=0.380278', 'beam.schedule.P=1221.3', 'beam.schedule.t_end=20.0', 'stop_time=20.0', 'timestep=0.007', 'amr.n_cell=128 128 64', 'amr.max_grid_size=32', 'weibull.V0=4.76837158203125e-10'],
}


def leg_argv(name):
    if name == "P400R":
        return MPI + [BIN, D19] + P400R + ["plot_file=" + os.path.join(LEGS_DIR, name, "run")]
    return MPI + [BIN, D] + EXTRA[name] + ["plot_file=" + os.path.join(LEGS_DIR, name, "run")]


def rebuild():
    """the table from STOPS.json and the packet's rules (the generator, kept so the literal table can be audited)"""
    J = json.load(open(os.path.join(HERE, "STOPS.json")))
    Dk, Qk = ["0p5", "1", "2", "4", "8", "16", "71"], ["0p3", "0p5", "1", "2", "4", "8"]

    def base(d, q):
        a = [("beam.tan_delta", J["td"][d])] if J["td"][d] else []
        st = J["stops"]["D%s_Q%s" % (d, q)]
        return a + [("beam.schedule.P", "%.1f" % J["P"][q]), ("beam.schedule.t_end", "%.1f" % st),
                    ("stop_time", "%.1f" % st), ("timestep", "%g" % J["dt"]["D%s_Q%s" % (d, q)])]

    def merge(a, upd):
        x = dict(a)
        x.update(upd)
        return list(x.items())
    L = {}
    for q in Qk:
        for d in Dk:
            L["D%s_Q%s" % (d, q)] = base(d, q)
    for k, s in (("K05", "0.5"), ("K2", "2.0")):
        for d in Dk:
            L["%s_D%s_Q2" % (k, d)] = base(d, "2") + [("spallation.sp.kic.scale", s)]
    for w, om, P in (("W6", "0.00644", "152.7"), ("W19", "0.01932", "1374.0")):
        for d in ("1", "4"):
            L["%s_D%s_Q2" % (w, d)] = merge(base(d, "2"), {"beam.omega0": om, "beam.schedule.P": P})
    mesh = {"amr.n_cell": "128 128 64", "amr.max_grid_size": "32", "weibull.V0": "4.76837158203125e-10"}
    for d, dt in (("1", "0.008"), ("2", "0.01"), ("4", "0.02")):
        L["F_D%s_Q2" % d] = merge(base(d, "2"), dict(mesh, **{"stop_time": "20.0", "beam.schedule.t_end": "20.0",
                                                            "timestep": dt}))
    for q, dt, st in (("1", "0.02", "60.0"), ("4", "0.007", "20.0")):
        L["F_D2_Q%s" % q] = merge(base("2", q), dict(mesh, **{"stop_time": st, "beam.schedule.t_end": st,
                                                             "timestep": dt}))
    return {k: ["%s=%s" % kv for kv in v] for k, v in L.items()}


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
    print("  %-14s rc=%d" % (name, rc), flush=True)
    return rc


def main(args):
    if args[:1] == ["--check"]:
        ok = rebuild() == EXTRA
        print("run.py --check: the literal table %s the table rebuilt from STOPS.json (%d legs + P400R)" % (
            "EQUALS" if ok else "DIFFERS FROM", len(EXTRA)))
        return 0 if ok else 1
    if args[:1] == ["--list"]:
        print(" ".join(EXTRA))
        return 0
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
