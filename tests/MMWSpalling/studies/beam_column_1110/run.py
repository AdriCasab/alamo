#!/usr/bin/env python3
"""Beam column resolution study (2026-11-10, run directly at the user's request; not a packet).

QUESTION: how small must the cell be, relative to the beam's absorption depth D, before the things the regime map
turns on stop moving -- the first spall firing, the melt onset, the top cell's temperature, the melt depth?

FIXTURE: P2-25's regime deck (studies/p2_25_regime_map/input: the #19 deck, flat loss, one-window schedule) cut down
to a COLUMN of 4 x 4 cubic cells, 50 mm deep, insulated sides, under a beam 100x wider (omega0 1.288 m) at 1e4 x the
power, so the absorbed centre flux q = (1 - R) 2P/(pi w0^2) is the map's and the heating is laterally uniform (1-D).
Random flaws are OFF (weibull.enabled = 0: the fixed 8 um flaw), so every column is the same and the only thing that
differs between the legs of a series is the cell. The time step is 0.0008 s on EVERY leg (stable at the finest cell:
dz^2 / (6 kappa_max) = 0.93 ms at 0.098 mm, kappa_max 1.71e-6 m2/s), so the series is a pure spatial comparison.

LEGS per (D, q, N):  A = the model as it is (removal on);  H = heat only (spall.enabled = 0, vapor.enabled = 0, events
CSV off): the melt onset and the temperatures undisturbed by removal. Before A's first removal A and H are the same
run, so A's first event row is the spall onset of both. (The first launch called the heat-only leg B and left vapour
on; Parse.H refuses vapour without spall on the face list, rc 6 on all of them. H is that leg with vapour off; the
B directories were deleted.)

Usage: run.py --list | --argv LEG | LEG [LEG ...] | --all [-j K]
"""
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
DECK = os.path.join(ROOT, "tests", "MMWSpalling", "studies", "p2_25_regime_map", "input")
LEGS = os.path.join(HERE, "legs")
DT = "0.0008"
LZ = 0.05
NZ = [16, 32, 64, 128, 256, 512]                      # dz 3.125, 1.5625, 0.78125, 0.390625, 0.1953125, 0.09765625 mm
TD = {"0p5": "1.52111", "1": "0.760555", "2": "0.380278", "4": "0.190139", "8": "0.0950694"}   # P2-25's ladder
PW = {"2": 610.7, "4": 1221.3}                        # P2-25's powers at w0 12.88 mm (absorbed centre q 2, 4 MW/m2)
KINDS = ("A", "H")
# not run (skip_watch.sh stopped them on the first launch): 41 and 82 cells per absorption depth, hours each
SKIP = ("D4_Q2_N512_A", "D4_Q2_N512_H", "D8_Q2_N512_A", "D8_Q2_N512_H")
POINTS = [("0p5", "2", 25.0), ("1", "2", 25.0), ("2", "2", 25.0), ("4", "2", 40.0), ("8", "2", 40.0), ("2", "4", 12.0)]


def leg_keys(d, q, stop, n, kind):
    dz = LZ / n
    k = {"geometry.prob_hi": "%.10g %.10g %.10g" % (4 * dz, 4 * dz, LZ), "amr.n_cell": "4 4 %d" % n,
         "amr.max_grid_size": "32", "beam.x0": "%.10g" % (2 * dz), "beam.y0": "%.10g" % (2 * dz),
         "surface_patch.x0": "%.10g" % (2 * dz), "surface_patch.y0": "%.10g" % (2 * dz),
         "beam.omega0": "1.288", "beam.schedule.P": "%.6g" % (PW[q] * 1.0e4), "beam.tan_delta": TD[d],
         "beam.schedule.t_end": "%g" % stop, "stop_time": "%g" % stop, "timestep": DT,
         "weibull.enabled": "0", "amr.plot_int": "-1", "amr.thermo.plot_int": "5"}
    if kind == "H":
        k.update({"spall.enabled": "0", "vapor.enabled": "0", "spall.removal_events_csv": "0"})
    return k


def legs():
    L = {}
    for d, q, stop in POINTS:
        for n in NZ:
            for kind in KINDS:
                L["D%s_Q%s_N%d_%s" % (d, q, n, kind)] = leg_keys(d, q, stop, n, kind)
    return L


def argv(name, np_="1", extra=None, out=None):
    k = dict(legs()[name])
    k.update(extra or {})
    d = out or os.path.join(LEGS, name)
    return ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", np_, BIN, DECK] + \
        ["%s=%s" % kv for kv in k.items()] + ["plot_file=" + os.path.join(d, "run")]


def run(name, np_="1", extra=None, out=None):
    d = out or os.path.join(LEGS, name)
    if os.path.exists(os.path.join(d, "rc")) and open(os.path.join(d, "rc")).read().strip() == "0":
        return name, 0, "kept"
    os.makedirs(d, exist_ok=True)
    a = argv(name, np_, extra, out)
    open(os.path.join(d, "argv.sh"), "w").write(" ".join("'%s'" % t if " " in t else t for t in a) + "\n")
    with open(os.path.join(d, "run.log"), "w") as f:
        rc = subprocess.call(a, stdout=f, stderr=subprocess.STDOUT, cwd=ROOT)
    open(os.path.join(d, "rc"), "w").write("%d\n" % rc)
    return name, rc, "ran"


if __name__ == "__main__":
    a = sys.argv[1:]
    if a[:1] == ["--list"]:
        print("\n".join(legs()))
    elif a[:1] == ["--argv"]:
        print(" ".join(argv(a[1])))
    else:
        j = int(a[a.index("-j") + 1]) if "-j" in a else 6
        names = list(legs()) if "--all" in a else [x for x in legs() if x.endswith("_H")] if "--heat" in a \
            else [x for x in a if x in legs()]
        names = [x for x in names if x not in SKIP]
        names.sort(key=lambda s: int(s.split("_N")[1].split("_")[0]))          # the cheap cells first
        with ThreadPoolExecutor(j) as ex:
            for name, rc, how in ex.map(run, names):
                print("%-22s rc %d (%s)" % (name, rc, how), flush=True)
