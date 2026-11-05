#!/usr/bin/env python3
"""P2-13 -- the Kant 2017 Grimsel calibration deck: the legs, serial, under caffeinate.

Usage: run.py LEG [LEG ...]        (or --argv LEG: print the command and exit)
Legs (the deck's `input` plus exactly these argv tokens; ACTIVE_STEP.md P2-13 "Legs"):
  CAL        (none)                                              60 s
  CAL_det    weibull.enabled=0               (declared --allow)  60 s
  CAL_h2000  surface_patch.h_expr=2000.0 timestep=0.004          20 s
  CAL_p27    confining.p=27.0e6                                  60 s
  CAL_p48    confining.p=48.0e6                                  60 s
  CAL_2mm    amr.n_cell=12 12 32 weibull.V0=8.0e-9               60 s
  CAL_05mm   amr.n_cell=48 48 128 weibull.V0=1.25e-10            60 s
  RECAL_6 / RECAL_10 / RECAL_x   weibull.a0_ig = spallation.sp.a0 = 6 / 10 um / RECAL_X.txt (declared --allow)
  HOS_R      validation/hoeser/input_lab (re-pointed) + G1's argv minus every shared key; stop 15 s
  SMOKE      the deck, one step (a2's check_shared run)
  SMOKE_<leg>  any leg above, stopped at 10 s (the witness: firings > 0)
Each leg: legs/<leg>/{argv.sh, BINARY.txt, run.log, rc, .done, run/, run_removal_events.csv}  (SR 8).
"""
import hashlib
import os
import shlex
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", "..", ".."))
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
DECK = os.path.join(HERE, "input")
LAB = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "hoeser", "input_lab")
G1 = os.path.join(ROOT, "tests", "MMWSpalling", "studies", "p2_11_heat", "legs", "G1", "argv.sh")
SHARED = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "shared")
LEGS_DIR = os.path.join(HERE, "legs")
MPI = ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4"]

FLAW = lambda a: ["weibull.a0_ig=%s" % a, "spallation.sp.a0=%s" % a]
LEGS = {  # name: (tokens, the --allow keys check_shared must be given)
    "CAL":       ([], []),
    "CAL_det":   (["weibull.enabled=0"], ["weibull.enabled"]),
    "CAL_h2000": (["surface_patch.h_expr=2000.0", "timestep=0.004", "stop_time=20.0"], []),
    "CAL_p27":   (["confining.p=27.0e6"], []),
    "CAL_p48":   (["confining.p=48.0e6"], []),
    "CAL_2mm":   (["amr.n_cell=12 12 32", "weibull.V0=8.0e-9"], []),
    "CAL_05mm":  (["amr.n_cell=48 48 128", "weibull.V0=1.25e-10"], []),
    "RECAL_6":   (FLAW("6.0e-6"), ["weibull.a0_ig", "spallation.sp.a0"]),
    "RECAL_10":  (FLAW("10.0e-6"), ["weibull.a0_ig", "spallation.sp.a0"]),
}


def shared_keys():
    ks = set()
    for f in ("rock_grimsel.inputs", "model_form_flame.inputs"):
        for ln in open(os.path.join(SHARED, f)):
            ln = ln.split("#", 1)[0]
            if "=" in ln:
                ks.add(ln.split("=", 1)[0].strip())
    return ks


def hos_r_tokens(d):
    """G1's argv minus every token now in a shared file; stop 15 s; plot_file re-pointed"""
    toks = shlex.split(open(G1).read())
    i = toks.index(LAB)
    sk, out, dropped = shared_keys(), [], []
    for t in toks[i + 1:]:
        k = t.split("=", 1)[0]
        if k in sk:
            dropped.append(t)
        elif k == "stop_time":
            out.append("stop_time=15.0")
        elif k == "plot_file":
            out.append("plot_file=" + os.path.join(d, "run"))
        else:
            out.append(t)
    return out, dropped


def leg_argv(name):
    """(argv, allow-keys, notes)"""
    d = os.path.join(LEGS_DIR, name)
    if name == "HOS_R":
        toks, dropped = hos_r_tokens(d)
        return MPI + [BIN, LAB] + toks, [], dropped
    smoke = name.startswith("SMOKE_")
    base = name[6:] if smoke else name
    if base == "SMOKE":
        toks, allow = ["stop_time=0.008"], []
    elif base == "RECAL_x":
        a = open(os.path.join(HERE, "RECAL_X.txt")).read().split()[0]
        toks, allow = FLAW(a), ["weibull.a0_ig", "spallation.sp.a0"]
    else:
        toks, allow = LEGS[base]
        toks = list(toks)
    if smoke:
        toks = [t for t in toks if not t.startswith("stop_time=")] + ["stop_time=10.0"]
    return MPI + [BIN, DECK] + toks + ["plot_file=" + os.path.join(d, "run")], allow, []


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def run(name):
    d = os.path.join(LEGS_DIR, name)
    if os.path.isdir(d):                     # SR 9: archive, never overwrite
        k = 1
        while os.path.isdir(d + ".old%d" % k):
            k += 1
        shutil.move(d, d + ".old%d" % k)
    os.makedirs(d)
    a, allow, dropped = leg_argv(name)
    open(os.path.join(d, "argv.sh"), "w").write(" ".join(shlex.quote(x) for x in a) + "\n")
    open(os.path.join(d, "ALLOW.txt"), "w").write(" ".join(allow) + "\n")
    if dropped:
        open(os.path.join(d, "DROPPED.txt"), "w").write("\n".join(dropped) + "\n")
    open(os.path.join(d, "BINARY.txt"), "w").write(sha(BIN) + "  bin/mmwspalling-3d-g++\n")
    with open(os.path.join(d, "run.log"), "w") as log:
        rc = subprocess.call(a, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    open(os.path.join(d, "rc"), "w").write("%d\n" % rc)
    if rc == 0:
        open(os.path.join(d, ".done"), "w").write("")
    print("  %-14s rc=%d" % (name, rc), flush=True)
    return rc


def UNITS():
    """a4: every kept unit test on this binary (P2-10's `units` harness, output HERE), compared with
    studies/p2_12_retire_path/restore/UNIT_p2_12r exactly as restore/gates_restore.py UNITS compares"""
    import importlib.machinery
    import importlib.util
    stud = os.path.join(ROOT, "tests", "MMWSpalling", "studies")
    ld = importlib.machinery.SourceFileLoader("r10_p213", os.path.join(stud, "p2_10_retire", "run_p2_10.py"))
    R10 = importlib.util.module_from_spec(importlib.util.spec_from_loader("r10_p213", ld))
    ld.exec_module(R10)
    R10.HERE = HERE
    label = "p2_13"
    R10.units(label)
    A = os.path.join(stud, "p2_12_retire_path", "restore", "UNIT_p2_12r")
    B = os.path.join(HERE, "UNIT_" + label)
    rc = lambda d: {l.split()[0]: l.split()[1] for l in open(os.path.join(d, "rc.txt")).read().splitlines()[1:]}
    ra, rb = rc(A), rc(B)
    L = ["a4 -- every kept unit test with input_lab re-pointed vs P2-12r's UNIT_p2_12r: rc and verdict lines",
         "=" * 100]
    bad = []
    for n in sorted(set(ra) | set(rb)):
        va = R10.verdict_lines(os.path.join(A, n + ".log")) if n in ra else None
        vb = R10.verdict_lines(os.path.join(B, n + ".log")) if n in rb else None
        same = ra.get(n) == rb.get(n) and va == vb
        L.append("  %-28s rc %5s -> %5s   verdict lines %s" % (n, ra.get(n, "-"), rb.get(n, "-"),
                 "identical (%d)" % len(vb) if same else "DIFFER"))
        if not same:
            bad.append(n)
    L.append("a4 (every test): %s" % ("PASS" if not bad else "FAIL: " + ", ".join(bad)))
    txt = "\n".join(L)
    print(txt)
    open(os.path.join(HERE, "UNITS_%s_vs_p2_12r.txt" % label), "w").write(txt + "\n")
    return 0 if not bad else 1


def main(args):
    if args[:1] == ["--argv"]:
        for n in args[1:]:
            a, allow, dropped = leg_argv(n)
            print(n, " ".join(shlex.quote(x) for x in a[len(MPI) + 1:]), "| allow:", allow,
                  "| dropped %d shared tokens" % len(dropped) if dropped else "")
        return 0
    subprocess.Popen(["caffeinate", "-dis", "-w", str(os.getpid())])   # SR 11
    if args == ["UNITS"]:
        return UNITS()
    bad = 0
    for n in args:
        bad += run(n) != 0
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
