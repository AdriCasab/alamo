#!/usr/bin/env python3
"""P2-14a -- the legs not run by unit/face_cref (that unit runs F0_k0/k1, F2_k0/k1, F164_k1, T35_0/2, N2_0/2 on the
new binary; the reducer reads its output). Serial, -np 4 unless named, under caffeinate; per leg legs/<leg>/{argv.sh,
BINARY.txt, ALLOW.txt, run.log, rc, .done} (SR 8). A *_pre leg runs bin_pre/ (sha 0a96a095...): rows a6, a7, a8 only.

Usage: run.py LEG [LEG ...] | --argv LEG ... | ALL | UNITS
"""
import hashlib
import os
import shlex
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
BIN_PRE = os.path.join(HERE, "bin_pre", "mmwspalling-3d-g++")
F_IN = os.path.join(ROOT, "tests", "MMWSpalling", "unit", "face_removal", "input")
CALD = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "kant", "grimsel_calibration")
CAL_IN = os.path.join(CALD, "input")
LAB = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "hoeser", "input_lab")
LEGS_DIR = os.path.join(HERE, "legs")
MPI = ["mpirun", "--oversubscribe", "--bind-to", "none"]
B1 = ["amr.n_cell=24 24 64", "timestep=0.008", "amr.max_grid_size=16", "amr.plot_int=125"]
B2 = ["amr.n_cell=12 12 32", "timestep=0.008", "amr.max_grid_size=8", "amr.plot_int=125"]
B05 = ["amr.n_cell=48 48 128", "timestep=0.004", "amr.max_grid_size=32", "amr.plot_int=250"]
D2 = ["spallation.sp.constraint_depth=0.002"]
K0, K1 = ["surface_patch.face_list=0"], ["surface_patch.face_list=1"]
D_KEY = "spallation.sp.constraint_depth"

# name: (kind, tokens, pre?, np)   kind F = fixture F; P13:<leg> = P2-13's recorded argv; CAL = the deck
LEGS = {
    "F0_k0_pre": ("F", B1 + K0, True, 4), "F0_k1_pre": ("F", B1 + K1, True, 4),
    "F2_k0_pre": ("F", B1 + K0 + D2, True, 4), "F2_k1_pre": ("F", B1 + K1 + D2, True, 4),
    "F2_k1_np1": ("F", B1 + K1 + D2, False, 1),
    "F2_k1_g8": ("F", B1 + K1 + D2 + ["amr.max_grid_size=8"], False, 4),
    "M2_k1": ("F", B2 + K1 + D2, False, 4), "M05_k1": ("F", B05 + K1 + D2, False, 4),
    "CAL": ("P13:CAL", [], False, 4), "CAL_det": ("P13:CAL_det", [], False, 4),
    "CAL_h2000": ("P13:CAL_h2000", [], False, 4), "CAL_2mm": ("P13:CAL_2mm", [], False, 4),
    "CAL_p27": ("P13:CAL_p27", [], False, 4), "CAL_p48": ("P13:CAL_p48", [], False, 4),
    "HOS_R": ("P13:HOS_R", [], False, 4),
    "CAL_p27_mu": ("P13:CAL_p27", ["microstructure.phase0.mu=12.0e9"], False, 4),
    "CAL_p27_mu_pre": ("P13:CAL_p27", ["microstructure.phase0.mu=12.0e9"], True, 4),
    "C0_h150": ("CAL", ["surface_patch.h_expr=150.0", "stop_time=200.0"], False, 4),
    "C2_h150": ("CAL", ["surface_patch.h_expr=150.0", "stop_time=200.0"] + D2, False, 4),
    "C2_h500": ("CAL", D2, False, 4),
    "C0_h1000": ("CAL", ["surface_patch.h_expr=1000.0", "timestep=0.004", "stop_time=30.0"], False, 4),
    "C2_h1000": ("CAL", ["surface_patch.h_expr=1000.0", "timestep=0.004", "stop_time=30.0"] + D2, False, 4),
    "C2_h2000": ("CAL", ["surface_patch.h_expr=2000.0", "timestep=0.004", "stop_time=20.0"] + D2, False, 4),
    "HOS_d2": ("P13:HOS_R", D2, False, 4),
}
ORDER = ["F2_k0_pre", "F2_k1_pre", "F0_k0_pre", "F0_k1_pre", "CAL_p27_mu_pre",
         "CAL", "CAL_det", "CAL_h2000", "CAL_2mm", "CAL_p27", "CAL_p48", "HOS_R", "CAL_p27_mu",
         "F2_k1_np1", "F2_k1_g8", "M2_k1", "M05_k1",
         "C0_h150", "C2_h150", "C2_h500", "C0_h1000", "C2_h1000", "C2_h2000", "HOS_d2"]


def p13_tokens(leg):
    """P2-13's recorded argv after its input file, without plot_file"""
    toks = shlex.split(open(os.path.join(CALD, "legs", leg, "argv.sh")).read())
    i = next(k for k, t in enumerate(toks) if t.endswith("mmwspalling-3d-g++"))
    deck = toks[i + 1]
    return deck, [t for t in toks[i + 2:] if not t.startswith("plot_file=")]


def p13_allow(leg):
    p = os.path.join(CALD, "legs", leg, "ALLOW.txt")
    return open(p).read().split() if os.path.exists(p) else []


def leg_argv(name):
    kind, toks, pre, npr = LEGS[name]
    d = os.path.join(LEGS_DIR, name)
    allow = []
    if kind == "F":
        deck, base = F_IN, []
    elif kind == "CAL":
        deck, base = CAL_IN, []
    else:
        deck, base = p13_tokens(kind[4:])
        allow = p13_allow(kind[4:])
    toks = base + list(toks)
    if any(t.startswith(D_KEY + "=") for t in toks):
        allow.append(D_KEY)
    if any(t.startswith("microstructure.phase0.mu=") for t in toks):
        allow.append("microstructure.phase0.mu")
    a = MPI + ["-np", str(npr), BIN_PRE if pre else BIN, deck] + toks + ["plot_file=" + os.path.join(d, "run")]
    return a, deck, allow


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
    open(os.path.join(d, "BINARY.txt"), "w").write(sha(a[6]) + "  " + os.path.relpath(a[6], ROOT) + "\n")
    with open(os.path.join(d, "run.log"), "w") as log:
        rc = subprocess.call(a, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    open(os.path.join(d, "rc"), "w").write("%d\n" % rc)
    open(os.path.join(d, ".done"), "w").write("")
    print("  %-16s rc=%d" % (name, rc), flush=True)
    return rc


def UNITS(only=None):
    import importlib.machinery
    import importlib.util
    stud = os.path.join(ROOT, "tests", "MMWSpalling", "studies")
    ld = importlib.machinery.SourceFileLoader("r10_p214a", os.path.join(stud, "p2_10_retire", "run_p2_10.py"))
    R10 = importlib.util.module_from_spec(importlib.util.spec_from_loader("r10_p214a", ld))
    ld.exec_module(R10)
    R10.HERE = HERE
    R10.units("p2_14a", only)
    return 0


def main(args):
    if args[:1] == ["--argv"]:
        for n in args[1:]:
            a, deck, allow = leg_argv(n)
            print(n, " ".join(shlex.quote(x) for x in a[4:]), "| allow:", allow)
        return 0
    subprocess.Popen(["caffeinate", "-dis", "-w", str(os.getpid())])   # SR 11
    if args[:1] == ["UNITS"]:
        return UNITS(args[1:] or None)
    names = ORDER if args == ["ALL"] else args
    for n in names:
        run(n)                       # an rc is a row's input (a8 expects a non-zero one), not a stop
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
