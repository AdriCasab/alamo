#!/usr/bin/env python3
"""P2-15 -- the CAL-deck legs (fixture F's legs are unit/face_buckle's). Serial, -np 4, under caffeinate; per leg
legs/<leg>/{argv.sh, BINARY.txt, ALLOW.txt, DECK.txt, run.log, rc, .done} (SR 8). Plotfiles: only the last pair is
kept (packet guardrail), except CAL and HOS_R (a4 compares every plotfile). A *_S leg reads SPAN.txt at launch.

Usage: run.py LEG [LEG ...] | --argv LEG ... | UNITS
"""
import glob
import hashlib
import os
import re
import shlex
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
CALD = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "kant", "grimsel_calibration")
CAL_IN = os.path.join(CALD, "input")
LEGS_DIR = os.path.join(HERE, "legs")
MPI = ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4"]
# Kant 2016 Table 4 (h_fl W/m2K, theta_f degC); 140/160 linear in SOD (studies/law_form/RATIO.txt)
KT = {50: (766, 1450), 80: (1009, 1410), 110: (873, 1360), 130: (747, 1250), 140: (690, 1175),
      150: (633, 1100), 160: (576, 1025)}
GT = ["microstructure.phase0.grimsel.lambda.a=1923.5", "microstructure.phase0.grimsel.rho=2594"]
BT = ["microstructure.phase0.grimsel.lambda.a=1540.6", "microstructure.phase0.grimsel.rho=2566"]
TKEYS = ["microstructure.phase0.grimsel.lambda.a", "microstructure.phase0.grimsel.rho"]
KEEP_ALL_PF = {"CAL", "HOS_R"}


def kt(sod):
    h, th = KT[sod]
    return ["surface_patch.h_expr=%.1f" % h, "surface_patch.T_flame_expr=%.2f" % (th + 273.15)]


def span():
    return open(os.path.join(HERE, "SPAN.txt")).read().split()[0]


def legs():
    L = {"CAL": ("CAL", [], [], False), "HOS_R": ("HOS", [], [], False)}
    for sod, stop in ((130, None), (140, 200), (150, 200), (160, 200)):
        L["KG%d_0" % sod] = ("CAL", kt(sod) + GT + ([] if stop is None else ["stop_time=200.0"]), TKEYS, False)
    for sod in (150, 160):
        L["KB%d_0" % sod] = ("CAL", kt(sod) + BT + ["stop_time=200.0"], TKEYS, False)
    for sod, stop in ((50, None), (80, None), (110, None), (130, None), (140, 200), (150, 200), (160, 200)):
        L["KG%d_S" % sod] = ("CAL", kt(sod) + GT + ([] if stop is None else ["stop_time=200.0"]), TKEYS, True)
    for sod in (150, 160):
        L["KB%d_S" % sod] = ("CAL", kt(sod) + BT + ["stop_time=200.0"], TKEYS, True)
    L["S_h150"] = ("CAL", ["surface_patch.h_expr=150.0", "stop_time=200.0"], [], True)
    L["S_h500"] = ("CAL", [], [], True)
    L["S_h1000"] = ("CAL", ["surface_patch.h_expr=1000.0", "timestep=0.004", "stop_time=30.0"], [], True)
    L["S_h2000"] = ("CAL", ["surface_patch.h_expr=2000.0", "timestep=0.004", "stop_time=20.0"], [], True)
    L["HOS_S"] = ("HOS", [], [], True)
    return L


def hos_tokens():
    toks = shlex.split(open(os.path.join(CALD, "legs", "HOS_R", "argv.sh")).read())
    i = next(k for k, t in enumerate(toks) if t.endswith("mmwspalling-3d-g++"))
    return toks[i + 1], [t for t in toks[i + 2:] if not t.startswith("plot_file=")]


def leg_argv(name):
    kind, toks, allow, sp = legs()[name]
    d = os.path.join(LEGS_DIR, name)
    deck, base = (CAL_IN, []) if kind == "CAL" else hos_tokens()
    toks = base + list(toks)
    allow = list(allow)
    if sp:
        toks.append("spall.buckle_span=" + span())
        allow.append("spall.buckle_span")
    return MPI + [BIN, deck] + toks + ["plot_file=" + os.path.join(d, "run")], deck, allow


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def prune(d):
    pf = sorted(glob.glob(os.path.join(d, "run", "[0-9]" * 5 + "cell")))
    for p in pf[:-1]:
        shutil.rmtree(p)
        shutil.rmtree(p[:-4] + "node", ignore_errors=True)


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
    if name not in KEEP_ALL_PF:
        prune(d)
    print("  %-9s rc=%d" % (name, rc), flush=True)
    return rc


def UNITS(only=None):
    import importlib.machinery
    import importlib.util
    stud = os.path.join(ROOT, "tests", "MMWSpalling", "studies")
    ld = importlib.machinery.SourceFileLoader("r10_p215", os.path.join(stud, "p2_10_retire", "run_p2_10.py"))
    R10 = importlib.util.module_from_spec(importlib.util.spec_from_loader("r10_p215", ld))
    ld.exec_module(R10)
    R10.HERE = HERE
    R10.units("p2_15", only)
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
    for n in args:
        run(n)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
