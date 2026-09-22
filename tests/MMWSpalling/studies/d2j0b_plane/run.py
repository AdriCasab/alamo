#!/usr/bin/env python3
"""D2j-0b harness: the D2j-0 hot disk under a descending exclusion plane.
Key-only; the D2j-0 binary. See CRITERION.md (read-only, hashed before launch).

  run.py --cases C0_2mm C1_2mm C2_2mm C3_2mm C1d_2mm [--jobs 3] [--force]
  run.py --cases C0_1mm C1_1mm C2_1mm [--jobs 3]
  run.py --list
  run.py --kill NAME              binary path + plot_file match (never pkill -f)

Imports studies/d2j0_hotdisk/run.py by path for run_one / kill / keep_awake /
binary_sha (not edited). Everything patch-related is replaced here:
  * quarter domain 0.08 x 0.08 x 0.15 m, surface at z = 0.15 m;
  * prescribed plane: nozzle_z0 0.170 m, nozzle_feed 4.3056e-4 m/s (1.55 m/h),
    nozzle_collision_radius 0.002 m;
  * h_expr / T_flame_expr in s only (the parser's if(cond, a, b)); r enters
    only through surface_patch.radius = 0.030 (h = 0 outside the patch);
  * dt 16 ms at 2 mm, 8 ms at 1 mm (C1d: 8 ms at 2 mm), blocking_factor 1,
    mgs 20 / 40, stop 250 s, plotfiles every 50 s.
"""
import argparse
import concurrent.futures as cf
import importlib.util
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
D2J0 = os.path.join(os.path.dirname(HERE), "d2j0_hotdisk")
_spec = importlib.util.spec_from_file_location("d2j0_run", os.path.join(D2J0, "run.py"))
J = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(J)

OUT = os.path.join(HERE, "output")
J.OUT = OUT                      # run_one's makedirs and kill()'s token
INPUT = J.INPUT                  # the D2j-0 input_hotdisk, unchanged

LXY, LZ = 0.08, 0.15
R_DISK = 0.030
H_IN, T_IN = 700.0, 1600.0
Z0, FEED, R_COLL = 0.170, 4.3056e-4, 0.002
STOP, PLOT_S = 250.0, 50.0
MESH = {"2mm": 2e-3, "1mm": 1e-3}
DT = {"2mm": 16e-3, "1mm": 8e-3}
MGS = {"2mm": 20, "1mm": 40}

# (h above the plane, T above the plane); below the plane always (H_IN, T_IN)
WALL = {
    "C0": (H_IN, T_IN),          # no exclusion
    "C1": (1.0e-3, T_IN),        # exclusion emulated (parser forbids h <= 0)
    "C2": (100.0, 1000.0),       # exposed rock held warm (~780 K), no spalling
    "C3": (100.0, T_IN),         # exhaust-like wall, exploratory, 2 mm only
    "C1d": (1.0e-3, T_IN),       # C1 at dt 8 ms, 2 mm only
}
CASES = ["C0_2mm", "C1_2mm", "C2_2mm", "C3_2mm", "C1d_2mm", "C0_1mm", "C1_1mm", "C2_1mm",
         # exploratory (CRITERION_addendum.md): the Step 20 surface-normal kernel on
         "C1n_2mm", "C0n_2mm", "C1n_1mm"]


def exprs(case):
    h_up, T_up = WALL[case]
    if h_up == H_IN:
        h = f"{H_IN!r}"
    else:
        h = f"if(s>0.0,{H_IN!r},{h_up!r})"
    T = f"{T_IN!r}" if T_up == T_IN else f"if(s>0.0,{T_IN!r},{T_up!r})"
    return h, T


def job(name):
    case, mesh = name.split("_")
    dz = MESH[mesh]
    dt = 8e-3 if case == "C1d" else DT[mesh]
    sn = case.endswith("n")
    if sn:
        case = case[:-1]
    nxy, nz = int(round(LXY / dz)), int(round(LZ / dz))
    h, T = exprs(case)
    q = J.q_net()
    k = [f"stop_time={STOP!r}", f"timestep={dt!r}", f"amr.plot_int={int(round(PLOT_S / dt))}",
         f"amr.n_cell={nxy} {nxy} {nz}", f"geometry.prob_hi={LXY!r} {LXY!r} {LZ!r}",
         f"amr.max_grid_size={MGS[mesh]}", f"weibull.V0={dz ** 3!r}", "amr.blocking_factor=1",
         f"surface_patch.radius={R_DISK!r}",
         f'surface_patch.h_expr="{h}"', f'surface_patch.T_flame_expr="{T}"',
         "surface_patch.nozzle_descent=prescribed",
         f"surface_patch.nozzle_z0={Z0!r}", f"surface_patch.nozzle_feed={FEED!r}",
         f"surface_patch.nozzle_collision_radius={R_COLL!r}"]
    if sn:
        k.append("spall.surface_normal=1")
    pf = os.path.join(OUT, name)
    cmd = (["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4", J.BIN, INPUT,
            f"plot_file={pf}"] + k)
    meta = dict(name=name, case=case, mesh=mesh, dz=dz, nxy=nxy, nz=nz, lxy=LXY, lz=LZ, dt=dt,
                stop=STOP, plot_s=PLOT_S, max_grid_size=MGS[mesh], V0=dz ** 3, patch_r=R_DISK,
                h_in=H_IN, T_in=T_IN, h_up=WALL[case][0], T_up=WALL[case][1], z0=Z0, feed=FEED,
                r_coll=R_COLL, h_expr=h, T_expr=T, q_net=q, surface_normal=int(sn),
                overshoot=q * dt / (J.RHOCP * dz), binary_sha256=J.binary_sha())
    return pf, cmd, meta


def run_all(names, jobs, force):
    J.keep_awake()
    js = [job(n) for n in names]
    print(f"{len(js)} runs, {min(jobs, len(js))} in parallel -> {OUT}", flush=True)
    for _, _, m in js:
        print(f"  {m['name']}: h_up {m['h_up']:g} T_up {m['T_up']:g}; dz {m['dz'] * 1e3:g} mm, "
              f"{m['nxy']}²×{m['nz']}, dt {m['dt'] * 1e3:g} ms ({m['overshoot']:.2f} K/step), "
              f"stop {m['stop']:g} s, mgs {m['max_grid_size']}", flush=True)
    ok = True
    with cf.ThreadPoolExecutor(max_workers=min(jobs, len(js))) as ex:
        for fu in cf.as_completed([ex.submit(J.run_one, pf, c, m, force) for pf, c, m in js]):
            m, status, wall = fu.result()
            print(f"  {m['name']:8s} {status} at t = {m.get('t_end', 0):g} s "
                  f"({wall / 60:.1f} min)", flush=True)
            ok = ok and not status.startswith("FAILED")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", nargs="*")
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--kill")
    a = ap.parse_args()
    if a.kill:
        return J.kill(a.kill)
    if a.list:
        print(f"binary {J.binary_sha()[:8]}…  q_net(821 K) = {J.q_net() / 1e6:.3f} MW/m², "
              f"feed {FEED * 3600:.3f} m/h, s0 {Z0 - LZ:g} m")
        for n in CASES:
            _, cmd, m = job(n)
            print(f"{n:8s} dz {m['dz'] * 1e3:g} mm {m['nxy']}²×{m['nz']} dt {m['dt'] * 1e3:g} ms "
                  f"{m['overshoot']:.2f} K/step  h={m['h_expr']}  T={m['T_expr']}")
        return
    bad = [n for n in (a.cases or []) if n not in CASES]
    if not a.cases or bad:
        sys.exit(f"give --cases from {CASES}" + (f" (unknown: {bad})" if bad else ""))
    sys.exit(0 if run_all(a.cases, a.jobs, a.force) else 1)


if __name__ == "__main__":
    main()
