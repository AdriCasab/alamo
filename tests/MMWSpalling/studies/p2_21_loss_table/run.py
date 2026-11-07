#!/usr/bin/env python3
"""P2-21 -- the legs. Serial, -np 4, under caffeinate; per leg legs/<leg>/{argv.sh, BINARY.txt, run.log, rc, .done}
(SR 8); an existing leg dir is archived to .old<k> first (SR 9).

Z = validation/zhang/zhang_oglesby/input as edited (its FILE line includes the 17-pair Hartlieb table). Every beam
leg adds amr.node.any=0 amr.plot_int=50 (every 100 s; the per-step central-column columns carry the records).
  C1        Z + beam.tan_delta_T="295.15 1.0" amr.node.any=1 amr.plot_int=5 (ONE pair => flat: a3, byte-identical to
            p2_20_beam_deposition/legs/Z1/run)
  T1a/b/c   Z + test #19 (0.10 x 0.10 x 0.05 m, 128 x 128 x 64, beam.z0 0.05, 400 s) with beam.omega0 = 0.013 /
            0.0165 / 0.020, at Fig. 34's MEASURED incident power (ramp_schedule: 1.46 kW at beam-on to ~2.1 kW,
            40 intervals of 10 s) -- the USER'S DECISION 2026-11-06, replacing the packet's constant 2 kW on every
            T-leg (T0, M2, M3.5, T1d, T1x inherit it)
  T0        T1c + beam.tan_delta_T="295.15 1.0 1268.15 1.0" (two EQUAL pairs: the sum path on a flat coefficient)
  M2, M3.5  T1c with the table's ratio column scaled to its maximum m: r_m = 1 + (r - 1)(m - 1)/(r_max - 1)
  T1d       T1b + the 100 mm block (0.10 x 0.10 x 0.10 m, 128 x 128 x 128, beam.z0 0.10)
  T1x       T1b + solver=explicit timestep=0.03 amr.plot_int=333 stop_time=max(60, 1.25 x T1b's onset) (b3;
            T1b's onset is read from its thermo.dat, so T1b must have run)
  T2        Z + beam.omega0=0.0141 (#28c at Oglesby's spot) with amr.plot_int=25 (DECLARED deviation from the
            packet's 50: c9 reads 50 / 100 / 150 s, which plot_int 50 at dt 2 s would not write)
  SMOKE_<leg>  the leg on 32 x 32 x 16 (T1-type legs; 64 x 64 x 32 for T1d) for 10 steps (T1x: 50 steps), the
            witness diagnostics only

Usage: run.py LEG [LEG ...] | --argv LEG ...
"""
import hashlib
import math
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
TABLE = os.path.join(T, "validation", "shared", "loss_table_hartlieb.inputs")
MPI = ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4"]
BEAM = ["amr.node.any=0", "amr.plot_int=50"]
FIG34 = os.path.join(HERE, "doe", "FIG34.txt")


def ramp_schedule(t_stop=400.0, dt_w=10.0):
    """Fig. 34's MEASURED incident power on test #19 (doe/FIG34.txt: the digitised lower panel, 10 s windows) as a
    piecewise-constant beam.schedule from beam-on: interval i = [10 i, 10 (i+1)) s carries the window centred
    10 i + 5 s after beam-on; beam-on is taken 5 s before the first window with both traces (215 s => 210 s).
    THE USER'S DECISION (2026-11-06, after z0 found the packet's '~2 kW rising to ~2.3 kW' is the level AFTER
    onset; the trace reads 1.46 kW at beam-on, 1.96 kW just before the collapse): this ramp REPLACES the
    packet's constant 2 kW on every T-leg. Returns (t_end list [s], P list [W])."""
    lit = []
    for ln in open(FIG34):
        v = ln.split()
        if len(v) == 4 and v[0].isdigit():
            try:
                P = float(v[1])
            except ValueError:
                continue
            if P == P and P > 0.5:
                lit.append((float(v[0]), 1e3 * P))
    t_on = lit[0][0] - 0.5 * dt_w
    n = int(round(t_stop / dt_w))
    by_c = {round(t - t_on - 0.5 * dt_w): P for t, P in lit}
    t_end, Ps = [], []
    for i in range(n):
        key = round(i * dt_w)
        if key not in by_c:
            raise RuntimeError("FIG34.txt has no window %g s after beam-on" % (i * dt_w + 0.5 * dt_w))
        t_end.append((i + 1) * dt_w)
        Ps.append(by_c[key])
    return t_end, Ps


def ramp_tokens(t_stop=400.0):
    te, P = ramp_schedule(t_stop)
    return ["beam.schedule.t_end=" + " ".join("%g" % x for x in te), "beam.schedule.P=" + " ".join("%g" % x for x in P)]


T19 = ["geometry.prob_hi=0.10 0.10 0.05", "amr.n_cell=128 128 64", "beam.z0=0.05"] + ramp_tokens() + \
      ["stop_time=400.0"]
BLOCK100 = ["geometry.prob_hi=0.10 0.10 0.10", "amr.n_cell=128 128 128", "beam.z0=0.10"]


def table_pairs():
    """the 17 (T [K], r) pairs as the deck reads them (the one beam.tan_delta_T line of the shared file)"""
    for ln in open(TABLE):
        ln = ln.split("#", 1)[0]
        if ln.strip().startswith("beam.tan_delta_T"):
            v = [float(x) for x in ln.split("=", 1)[1].split()]
            return list(zip(v[0::2], v[1::2]))
    raise RuntimeError("no beam.tan_delta_T in " + TABLE)


def scaled(m):
    """the declared scan scaling on the ratio column only: r_m = 1 + (r - 1)(m - 1)/(r_max - 1)"""
    P = table_pairs()
    rmax = max(r for _, r in P)
    return [(t, 1.0 + (r - 1.0) * (m - 1.0) / (rmax - 1.0)) for t, r in P]


def table_token(pairs):
    return "beam.tan_delta_T=" + " ".join("%r %r" % (t, r) for t, r in pairs)


def onset_time(leg):
    """the first step row with a finite bc_z_melt, as time after beam-on (the step's start of the first lit row)"""
    p = os.path.join(LEGS_DIR, leg, "run", "thermo.dat")
    L = open(p).read().splitlines()
    h = [c.strip() for c in L[0].split("\t")]
    rows = [[float(x) for x in r.split("\t")] for r in L[1:] if r.strip()]
    rows = [r for r in rows if r[h.index("time")] > 0.0]
    dt = rows[1][0] - rows[0][0]
    t_on = next(r[h.index("time")] for r in rows if r[h.index("beam_P_inc")] > 0) - dt
    for r in rows:
        if math.isfinite(r[h.index("bc_z_melt")]):
            return r[h.index("time")] - t_on
    return None


def extra(name):
    if name == "C1":
        return ["beam.tan_delta_T=295.15 1.0", "amr.node.any=1", "amr.plot_int=5"]
    w = {"T1a": "0.013", "T1b": "0.0165", "T1c": "0.020"}
    if name in w:
        return T19 + ["beam.omega0=" + w[name]] + BEAM
    if name == "T0":
        return extra("T1c") + ["beam.tan_delta_T=295.15 1.0 1268.15 1.0"]
    if name in ("M2", "M3.5"):
        return extra("T1c") + [table_token(scaled(float(name[1:])))]
    if name == "T1d":
        return extra("T1b") + BLOCK100
    if name == "T1x":
        t_on = onset_time("T1b")
        if t_on is None:
            raise SystemExit("T1x needs T1b's onset, and T1b never melted")
        stop = max(60.0, math.ceil(1.25 * t_on))
        return extra("T1b") + ["solver=explicit", "timestep=0.03", "amr.plot_int=333", "stop_time=%g" % stop]
    if name == "T2":
        return ["beam.omega0=0.0141", "amr.node.any=0", "amr.plot_int=25"]
    if name.startswith("SMOKE_"):
        base = name[6:]
        if base == "T1x":
            e = extra("T1b") + ["solver=explicit", "timestep=0.03", "amr.plot_int=333", "stop_time=1.5"]
        else:
            e = extra(base) + ["stop_time=20.0"]
        if base == "T1d":
            e += ["amr.n_cell=32 32 32"]
        elif base not in ("C1", "T2"):
            e += ["amr.n_cell=32 32 16"]
        return e + ["amr.plot_int=100000"]
    raise KeyError(name)


def leg_argv(name):
    return MPI + [BIN, ZDECK] + extra(name) + ["plot_file=" + os.path.join(LEGS_DIR, name, "run")]


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
