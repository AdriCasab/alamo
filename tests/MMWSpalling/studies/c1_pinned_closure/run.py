#!/usr/bin/env python3
"""C1 study: pinned-surface Robin closure (surface_patch.robin_form = pinned).

Study, not a regression (no `test`). Nothing is copied: Part 1 reuses the S1
column base input + S1 `overrides()` (R-face source keys, then
robin_form=pinned); Part 2 overrides the Meier 2D dev harness input by CLI.

  run.py --part 1        1-D A/B against S1 (1 rank per run)  -> output/p1_*
  run.py --part 2        2-D machinery check (4 ranks per run) -> output/p2_*
  run.py --part 1 --cases P-1000 --dz 2 --force --jobs 4

Time step: the overshoot rule dT_step = q dt/(rho Cp dz) <= 5 K with q the
closed-form flux at the expected firing T; dt = min(1 ms, limit), halved
until it satisfies the rule. Plotfiles every 1 s (Part 1) / 2 s (Part 2).
"""
import argparse
import concurrent.futures as cf
import math
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
STUDIES = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(STUDIES, "s1_surface_resolution"))
from run_sweep import BIN, INPUT as S1_INPUT, ROOT, dz_tag, overrides  # noqa: E402

MEIER_INPUT = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "meier",
                           "sp_meier_pilot", "input_2d_dev")
RHOCP = 2750.0 * 790.0
SIG, EPS, HC, TA = 5.670374e-8, 0.8, 10.0, 293.15
OUT = os.path.join(HERE, "output")


def q_closed(h, T_gas, T_fire):
    """Pinned-surface closed form q = h(T_gas - T_f) - eps sig (T_f^4 - T_a^4) - h_c (T_f - T_a)."""
    return h * (T_gas - T_fire) - EPS * SIG * (T_fire ** 4 - TA ** 4) - HC * (T_fire - TA)


def dt_rule(q, dz):
    dt = 1.0e-3
    while q * dt / (RHOCP * dz) > 5.0:
        dt *= 0.5
    return dt


# ---------------------------------------------------------------- Part 1
# S1 R-cases: degenerate a = 20 um (base), p = 1 MPa, losses on; T_fire ~ 660 K.
T_FIRE_1D = 660.0
PART1 = {
    # name: (S1 source case, dz list, extra keys)
    "P-1000":        ("R-face",      [2.0, 1.0, 0.5, 0.25], ["surface_patch.pinned_idle_cycles=2.0"]),
    "P-1400":        ("R-face-1400", [2.0, 0.5],            ["surface_patch.pinned_idle_cycles=2.0"]),
    "P-1000-idle1":  ("R-face",      [2.0],                 ["surface_patch.pinned_idle_cycles=1.0"]),
    "P-1000-idle4":  ("R-face",      [2.0],                 ["surface_patch.pinned_idle_cycles=4.0"]),
}
T_GAS_1D = {"R-face": 1000.0, "R-face-1400": 1400.0}


def part1_job(name, dz_mm):
    src, _, extra = PART1[name]
    T_gas = T_GAS_1D[src]
    q = q_closed(1.0e4, T_gas, T_FIRE_1D)
    v = q / (RHOCP * (T_FIRE_1D - TA))            # m/s at the closed form
    stop = math.floor(min(60.0, 0.150 / v) * 10.0) / 10.0   # <= 150 mm of recession
    dz = dz_mm * 1e-3
    dt = dt_rule(q, dz)
    pf = os.path.join(OUT, f"p1_{name}_{dz_tag(dz_mm)}")
    cmd = [BIN, S1_INPUT] + overrides(src, dz_mm, stop, pf) + [
        "surface_patch.robin_form=pinned", f"timestep={dt!r}",
        f"amr.plot_int={int(round(1.0 / dt))}", f"amr.thermo.plot_int={max(1, int(round(0.01 / dt)))}",
        "spall.removal_events_csv=1"] + extra
    meta = dict(part=1, name=name, dz=dz_mm, dt=dt, stop=stop, q_cf=q,
                overshoot=q * dt / (RHOCP * dz), nranks=1)
    return pf, cmd, meta


# ---------------------------------------------------------------- Part 2
# Meier 2D harness, Weibull V0 = V_cell (user decision 2026-09-15): T_fire ~ 822 K.
T_FIRE_2D = 822.0
PART2 = {
    "M-5k-1200":  (5.0e3, 1200.0, [2.0]),
    "M-10k-1000": (1.0e4, 1000.0, [2.0, 1.0]),
    "M-20k-1400": (2.0e4, 1400.0, [2.0]),
}


def part2_job(name, dz_mm):
    h, T_gas, _ = PART2[name]
    q = q_closed(h, T_gas, T_FIRE_2D)
    v = q / (RHOCP * (T_FIRE_2D - TA))
    stop = math.floor(min(40.0, 0.080 / v) * 10.0) / 10.0
    dz = dz_mm * 1e-3
    dt = dt_rule(q, dz)
    pf = os.path.join(OUT, f"p2_{name}_{dz_tag(dz_mm)}")
    ov = [f"plot_file={pf}", f"stop_time={stop}", f"timestep={dt!r}",
          f"amr.plot_int={int(round(2.0 / dt))}", "amr.thermo.int=1",
          f"amr.thermo.plot_int={max(1, int(round(0.01 / dt)))}",
          "beam.P0=0.0",
          "surface_patch.enabled=1", "surface_patch.mode=convective_flame",
          "surface_patch.robin_form=pinned", "surface_patch.pinned_idle_cycles=2.0",
          "surface_patch.face=zhi", "surface_patch.x0=0.07", "surface_patch.y0=0.004",
          "surface_patch.radius=0.025", f"surface_patch.h_conv={h!r}",
          f"surface_patch.T_flame={T_gas!r}",
          "surface.follow_mask=1", "energy_ledger.enabled=1", "spall.removal_events_csv=1",
          f"weibull.V0={dz * dz * dz!r}"]
    if dz_mm != 2.0:
        ov += [f"amr.n_cell={int(round(0.14 / dz))} {int(round(0.008 / dz))} {int(round(0.120 / dz))}"]
    cmd = ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4", BIN, MEIER_INPUT] + ov
    meta = dict(part=2, name=name, dz=dz_mm, dt=dt, stop=stop, q_cf=q,
                overshoot=q * dt / (RHOCP * dz), nranks=4, h=h, T_gas=T_gas)
    return pf, cmd, meta


def run_one(pf, cmd, meta, force):
    done = pf + ".done"
    if os.path.exists(done) and not force:
        return meta, "skip (done)", 0.0
    for p in (pf, pf + "_h_col_events.csv", pf + "_removal_events.csv"):
        if os.path.isdir(p):
            shutil.rmtree(p)
        elif os.path.exists(p):
            os.remove(p)
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    with open(pf + ".log", "w") as log:
        log.write(" ".join(cmd) + "\n")
        log.flush()
        rc = subprocess.run(cmd, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT).returncode
    wall = time.time() - t0
    if rc == 0:
        with open(done, "w") as f:
            f.write(repr(meta) + "\n")
        return meta, "ok", wall
    return meta, f"FAILED rc={rc}", wall


def jobs_for(part, names=None, dzs=None):
    table, fn = (PART1, part1_job) if part == 1 else (PART2, part2_job)
    out = []
    for name in (names or table):
        if name not in table:
            continue
        dz_list = table[name][1] if part == 1 else table[name][2]
        for dz in dz_list:
            if dzs and dz not in dzs:
                continue
            out.append(fn(name, dz))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--part", type=int, choices=(1, 2), required=True)
    ap.add_argument("--cases", nargs="*", default=None)
    ap.add_argument("--dz", nargs="*", type=float, default=None)
    ap.add_argument("--jobs", type=int, default=None)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    if not os.path.exists(BIN):
        sys.exit(f"missing binary {BIN}; build first")
    jobs = jobs_for(args.part, args.cases, args.dz)
    nj = args.jobs or (6 if args.part == 1 else 2)
    print(f"{len(jobs)} runs, {nj} in parallel -> {OUT}")
    for pf, _, m in jobs:
        print(f"  {os.path.basename(pf)}: dt = {m['dt'] * 1e3:g} ms, stop = {m['stop']} s, "
              f"closed-form q = {m['q_cf'] / 1e6:.3f} MW/m^2, overshoot {m['overshoot']:.2f} K/step")
    ok = True
    with cf.ThreadPoolExecutor(max_workers=nj) as ex:
        for fu in cf.as_completed([ex.submit(run_one, pf, c, m, args.force) for pf, c, m in jobs]):
            m, status, wall = fu.result()
            print(f"  {m['name']:14s} dz = {m['dz']:<5g} mm  {status}  ({wall:.1f} s)", flush=True)
            ok = ok and status.startswith(("ok", "skip"))
    print("RUNS OK" if ok else "RUNS HAD FAILURES")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
