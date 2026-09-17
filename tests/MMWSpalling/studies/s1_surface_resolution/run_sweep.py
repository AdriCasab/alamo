#!/usr/bin/env python3
"""S1 surface-resolution sweep driver (study, not a regression).

Runs the case matrix of ACTIVE_STEP S1 on 1-D columns (1 x 1 x N, depth
250 mm) by CLI key=value overrides of the base `input`, one MPI rank per run,
several runs in parallel.

  run_sweep.py                 full matrix -> output/<case>_dz<mm>
  run_sweep.py --smoke         one dz (1 mm) per case, stop <= 10 s
                               -> output_smoke/, plus the 1-vs-4-rank check
  run_sweep.py --rank-check    only the 1-vs-4-rank identity check
  run_sweep.py --cases F1 R-face --dz 0.5   subset
  run_sweep.py --jobs 8        parallel runs (default 6)
  run_sweep.py --force         rerun even if a run is marked done

Stop time per case = min(60 s, 200 mm / v_hand), v_hand = the packet's hand
estimate, so the depth scan never reaches the column bottom -- except where
STOP_OVERRIDE replaces it with the measured time-to-200-mm (see below).
"""
import argparse
import concurrent.futures as cf
import glob
import math
import os
import re
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = HERE.split(os.sep + os.path.join("tests", "MMWSpalling") + os.sep)[0]
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
INPUT = os.path.join(HERE, "input")

DEPTH = 0.250          # column depth [m]
DRILL_MAX = 0.200      # never drill further than this at the hand ROP
BEAM_RADIUS = 0.1      # uniform disk radius >> column
H_ROBIN = 1.0e4

DZ_ALL = [2.0, 1.0, 0.5, 0.25, 0.125, 0.0625]   # mm

# name: (source dict, dz list [mm], hand ROP [mm/s])
CASES = {
    "F1":          (dict(kind="beam", q=2.06e6, alpha=2000.0),   DZ_ALL[:5], 1.83),
    "F1s":         (dict(kind="beam", q=2.06e6, alpha=32000.0),  DZ_ALL,     1.83),
    "R-cell":      (dict(kind="robin", form="cell", T_gas=1000.0), DZ_ALL[:5], 3.6),
    "R-face":      (dict(kind="robin", form="face", T_gas=1000.0), DZ_ALL,     1.67),
    "M-2000":      (dict(kind="beam", q=15.0e6, alpha=2000.0),   DZ_ALL[:5], 13.0),
    "M-100":       (dict(kind="beam", q=15.0e6, alpha=100.0),    DZ_ALL[:5], 13.0),
    "R-face-1400": (dict(kind="robin", form="face", T_gas=1400.0), [2.0, 0.5, 0.125], 5.2),
}


# Measured stop times [s]. The packet's hand ROPs assumed <dT_flake> ~ 519 K
# (E1 Meier 2D); in these columns the criterion fires at a top-cell T of
# ~656 K (<dT_rem> ~ 363 K), so R-cell, M-*, and the finest R-face runs drilled
# the full 250 mm before the hand stop time. These overrides are the earliest
# time (over the case's dz set, from the first sweep's h_col CSV) at which
# 200 mm had been removed, floored to 0.1 s. A shorter stop_time leaves the
# earlier trajectory bit-identical.
STOP_OVERRIDE = {"R-cell": 31.7, "R-face": 51.5, "M-2000": 10.6, "M-100": 10.9,
                 "R-face-1400": 27.6}


def case_stop(case):
    return STOP_OVERRIDE.get(case, stop_time(CASES[case][2]))


def dz_tag(dz_mm):
    s = f"{dz_mm:g}".replace(".", "p")
    return f"dz{s}"


def stop_time(v_mm_s):
    return math.floor(min(60.0, DRILL_MAX * 1e3 / v_mm_s) * 10.0) / 10.0


def overrides(case, dz_mm, stop, plot_file, max_grid_size=None):
    src, _, _ = CASES[case]
    dz = dz_mm * 1e-3
    N = int(round(DEPTH / dz))
    mgs = max_grid_size if max_grid_size else max(32, N // 4)
    c = 0.5 * dz
    ov = [f"plot_file={plot_file}",
          f"stop_time={stop}",
          f"amr.n_cell=1 1 {N}",
          f"geometry.prob_hi={dz!r} {dz!r} {DEPTH!r}",
          f"amr.max_grid_size={mgs}"]
    if src["kind"] == "beam":
        P0 = src["q"] * math.pi * BEAM_RADIUS ** 2
        ov += [f"beam.P0={P0!r}", f"beam.radius={BEAM_RADIUS!r}",
               f"beam.x0={c!r}", f"beam.y0={c!r}", f"beam.z0={DEPTH!r}",
               f"beam.alpha_expr={src['alpha']!r}",
               f"microstructure.phase0.alpha_attenuation={src['alpha']!r}"]
    else:
        ov += ["beam.P0=0.0",
               "surface_patch.enabled=1", "surface_patch.mode=convective_flame",
               "surface_patch.face=zhi", f"surface_patch.robin_form={src['form']}",
               f"surface_patch.x0={c!r}", f"surface_patch.y0={c!r}",
               f"surface_patch.radius={BEAM_RADIUS!r}",
               f"surface_patch.h_conv={H_ROBIN!r}",
               f"surface_patch.T_flame={src['T_gas']!r}"]
    return ov


def run_one(case, dz_mm, stop, plot_file, nranks=1, force=False, max_grid_size=None):
    done = plot_file + ".done"
    if os.path.exists(done) and not force:
        return case, dz_mm, "skip (done)", 0.0
    for p in (plot_file, plot_file + "_h_col_events.csv", plot_file + "_clusters.csv"):
        if os.path.isdir(p):
            shutil.rmtree(p)
        elif os.path.exists(p):
            os.remove(p)
    os.makedirs(os.path.dirname(plot_file), exist_ok=True)
    cmd = [BIN, INPUT] + overrides(case, dz_mm, stop, plot_file, max_grid_size)
    if nranks > 1:
        cmd = ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", str(nranks)] + cmd
    t0 = time.time()
    with open(plot_file + ".log", "w") as log:
        log.write(" ".join(cmd) + "\n")
        log.flush()
        rc = subprocess.run(cmd, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT).returncode
    wall = time.time() - t0
    if rc == 0:
        with open(done, "w") as f:
            f.write(f"{wall:.1f}\n")
        return case, dz_mm, "ok", wall
    return case, dz_mm, f"FAILED rc={rc}", wall


def read_plotfile_fields(pf):
    import numpy as np
    import yt
    yt.funcs.mylog.setLevel(50)
    ds = yt.load(pf)
    cg = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    out = {"__time__": np.array([float(ds.current_time)])}
    for f in ds.field_list:
        if f[0] == "boxlib":
            out[f[1]] = np.asarray(cg[f])
    return out


def rank_check(outroot, force):
    """Same case/dz/box layout on 1 and 4 ranks: every cell plotfile field,
    thermo.dat and the h_col CSV must be identical."""
    import numpy as np
    case, dz_mm, stop = "F1s", 0.5, 5.0
    N = int(round(DEPTH / (dz_mm * 1e-3)))
    runs = {}
    for nr in (1, 4):
        pf = os.path.join(outroot, f"rankcheck_np{nr}")
        r = run_one(case, dz_mm, stop, pf, nranks=nr, force=force, max_grid_size=N // 4)
        print(f"  rank check np={nr}: {r[2]} ({r[3]:.1f} s)")
        if not r[2].startswith(("ok", "skip")):
            return False
        runs[nr] = pf
    p1 = sorted(glob.glob(os.path.join(runs[1], "*cell")))
    p4 = sorted(glob.glob(os.path.join(runs[4], "*cell")))
    same = [os.path.basename(p) for p in p1] == [os.path.basename(p) for p in p4] and len(p1) > 0
    nfields = 0
    for a, b in zip(p1, p4):
        fa, fb = read_plotfile_fields(a), read_plotfile_fields(b)
        if fa.keys() != fb.keys():
            same = False
            break
        for k in fa:
            nfields += 1
            if not np.array_equal(fa[k], fb[k]):
                print(f"  rank check: {os.path.basename(a)} field {k} differs "
                      f"(max |d| = {np.max(np.abs(fa[k] - fb[k])):.3e})")
                same = False
    with open(runs[1] + "_h_col_events.csv", "rb") as f1, \
         open(runs[4] + "_h_col_events.csv", "rb") as f4:
        if f1.read() != f4.read():
            print("  rank check: h_col CSV differs")
            same = False
    # thermo.dat: every column must print identically except ledger_err, a
    # round-off residual (~1e-16) whose value depends on the MPI summation
    # order; it is only required to stay at round-off on both.
    t1 = [l.split() for l in open(os.path.join(runs[1], "thermo.dat"))]
    t4 = [l.split() for l in open(os.path.join(runs[4], "thermo.dat"))]
    if len(t1) != len(t4) or t1[0] != t4[0]:
        print("  rank check: thermo.dat layout differs")
        same = False
    else:
        ie = t1[0].index("ledger_err")
        for r1, r4 in zip(t1[1:], t4[1:]):
            if [v for i, v in enumerate(r1) if i != ie] != [v for i, v in enumerate(r4) if i != ie] \
               or abs(float(r1[ie])) > 1e-12 or abs(float(r4[ie])) > 1e-12:
                print(f"  rank check: thermo.dat row t = {r1[0]} differs")
                same = False
                break
    print(f"  rank check ({case}, dz = {dz_mm} mm, {N // 4}-cell boxes, t = {stop} s): "
          f"{len(p1)} plotfiles, {nfields} field comparisons + thermo.dat + h_col CSV -> "
          f"{'IDENTICAL' if same else 'DIFFERENT'}")
    return same


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--rank-check", action="store_true")
    ap.add_argument("--cases", nargs="*", default=None)
    ap.add_argument("--dz", nargs="*", type=float, default=None)
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if not os.path.exists(BIN):
        sys.exit(f"missing binary {BIN}; build first")

    outroot = os.path.join(HERE, "output_smoke" if (args.smoke or args.rank_check) else "output")
    os.makedirs(outroot, exist_ok=True)
    ok = True

    if not args.rank_check:
        jobs = []
        for case in (args.cases or CASES):
            _, dzs, v = CASES[case]
            if args.smoke:
                dzs = [1.0]
            if args.dz:
                dzs = [d for d in dzs if d in args.dz]
            stop = case_stop(case)
            if args.smoke:
                stop = min(stop, 10.0)
            for dz_mm in dzs:
                jobs.append((case, dz_mm, stop, os.path.join(outroot, f"{case}_{dz_tag(dz_mm)}")))
        print(f"{len(jobs)} runs, {args.jobs} in parallel, 1 rank each -> {outroot}")
        # Longest (finest) runs first so they are not left for last.
        jobs.sort(key=lambda j: (j[1], -j[2]))
        with cf.ThreadPoolExecutor(max_workers=args.jobs) as ex:
            futs = [ex.submit(run_one, c, d, s, p, 1, args.force) for c, d, s, p in jobs]
            for fu in cf.as_completed(futs):
                case, dz_mm, status, wall = fu.result()
                print(f"  {case:12s} dz = {dz_mm:<7g} mm  {status}  ({wall:.1f} s)", flush=True)
                if not status.startswith(("ok", "skip")):
                    ok = False

    if args.smoke or args.rank_check:
        ok = rank_check(outroot, args.force) and ok

    print("SWEEP OK" if ok else "SWEEP HAD FAILURES")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
