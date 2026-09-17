#!/usr/bin/env python3
"""S1b: firing temperature under the production (Meier) Weibull block, over dz.

Study, not a regression. Reuses the S1 column base input
(../s1_surface_resolution/input) and S1's `overrides()` for the dz / F1 beam
keys; nothing is copied. Common: F1 beam (q = 2.06 MW/m^2, alpha = 2000/m),
timestep = 1 ms, stop_time = 60 s, 1 rank, spall.removal_events_csv = 1,
follow_mask + ledger on (base).

Overshoot per step (09-15b rule, dT = q dt/(rho Cp dz)): 0.47 K at 2 mm,
0.95 K at 1 mm, 1.9 K at 0.5 mm, 3.8 K at 0.25 mm -- all <= 5 K.

  run.py                 all cases -> output/<case>_dz<mm>
  run.py --meier-keys    add the W-def 2 mm variants with one Meier key each
  run.py --cases W-def --dz 2 --jobs 4 --force
"""
import argparse
import concurrent.futures as cf
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
S1 = os.path.join(os.path.dirname(HERE), "s1_surface_resolution")
sys.path.insert(0, S1)
from run_sweep import BIN, INPUT, ROOT, dz_tag, overrides  # noqa: E402

STOP = 60.0
MEIER_WEIBULL = ["weibull.enabled=1", "weibull.a0_gb=20.0e-6", "weibull.a0_ig=4.0e-6",
                 "weibull.m=20.0", "weibull.seed=12345"]


def v0_cell(dz_mm):
    dz = dz_mm * 1e-3
    return [f"weibull.V0={dz * dz * dz!r}"]


# name: (dz list [mm], function dz -> extra CLI keys)
CASES = {
    "W-def":  ([2.0, 1.0, 0.5, 0.25], lambda dz: MEIER_WEIBULL),
    "W-cell": ([2.0, 1.0, 0.5, 0.25], lambda dz: MEIER_WEIBULL + v0_cell(dz)),
    "W-seed": ([2.0, 0.5],            lambda dz: MEIER_WEIBULL + ["weibull.seed=777"]),
    "S-4p4":  ([2.0],                 lambda dz: ["weibull.enabled=0", "spallation.sp.a0=4.436e-6"]),
}
# Meier 2D keys that differ from the S1 column, one at a time on W-def 2 mm
# (item 8 (ii): only needed if T_fire at 2 mm misses 812 K by > 15 K).
MEIER_KEYS = {
    "W-def+sn":   ([2.0], lambda dz: MEIER_WEIBULL + ["spall.surface_normal=1"]),
    "W-def+melt0": ([2.0], lambda dz: MEIER_WEIBULL + ["spall.melt_aware=0"]),
    "W-def+coh":  ([2.0], lambda dz: MEIER_WEIBULL + ["spall.flake_coherence_length=0.004"]),
}


def run_one(case, dz_mm, extra, plot_file, force):
    done = plot_file + ".done"
    if os.path.exists(done) and not force:
        return case, dz_mm, "skip (done)", 0.0
    for p in (plot_file, plot_file + "_h_col_events.csv", plot_file + "_removal_events.csv"):
        if os.path.isdir(p):
            shutil.rmtree(p)
        elif os.path.exists(p):
            os.remove(p)
    os.makedirs(os.path.dirname(plot_file), exist_ok=True)
    cmd = ([BIN, INPUT] + overrides("F1", dz_mm, STOP, plot_file)
           + ["timestep=0.001", "spall.removal_events_csv=1"] + extra)
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


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cases", nargs="*", default=None)
    ap.add_argument("--dz", nargs="*", type=float, default=None)
    ap.add_argument("--meier-keys", action="store_true")
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    if not os.path.exists(BIN):
        sys.exit(f"missing binary {BIN}; build first")
    table = dict(CASES)
    if args.meier_keys:
        table.update(MEIER_KEYS)
    names = args.cases or list(table)
    out = os.path.join(HERE, "output")
    jobs = []
    for c in names:
        dzs, fn = table[c]
        for dz in dzs:
            if args.dz and dz not in args.dz:
                continue
            jobs.append((c, dz, fn(dz), os.path.join(out, f"{c}_{dz_tag(dz)}")))
    jobs.sort(key=lambda j: j[1])
    print(f"{len(jobs)} runs, {args.jobs} in parallel -> {out}")
    ok = True
    with cf.ThreadPoolExecutor(max_workers=args.jobs) as ex:
        for fu in cf.as_completed([ex.submit(run_one, c, d, e, p, args.force) for c, d, e, p in jobs]):
            c, d, status, wall = fu.result()
            print(f"  {c:12s} dz = {d:<5g} mm  {status}  ({wall:.1f} s)", flush=True)
            ok = ok and status.startswith(("ok", "skip"))
    print("RUNS OK" if ok else "RUNS HAD FAILURES")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
