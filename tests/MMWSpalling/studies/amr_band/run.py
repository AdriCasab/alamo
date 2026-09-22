#!/usr/bin/env python3
"""AMR surface-band gate on the D2j-0b C1 hot disk (static disk under a
descending exclusion plane, quarter domain 0.08 x 0.08 x 0.15 m, no feet, no
jet). Same keys as d2j0b_plane C1; the AMR runs add a 2 mm base level with a
1 mm finest level on the surface band.

  run.py --cases U1_60 F1_60            identity gate (whole domain refined)
  run.py --cases U1_250 B1_250 [--jobs 2] band gate (20 mm band, regrid 2 s)
  run.py --list

Cases
  U1_<T>  uniform 1 mm, max_level 0, stop T s           (reference)
  U2_<T>  uniform 2 mm, max_level 0
  F1_<T>  2 mm base + level 1 covering the whole domain (must equal U1)
  B1_<T>  2 mm base + level 1 on the surface band only  (the real gate)

Binary: bin/mmwspalling-3d-g++-amr (built with POSTFIX=3d-g++-amr so the D2k
campaign's bin/mmwspalling-3d-g++ is never touched).
"""
import argparse
import ast
import concurrent.futures as cf
import hashlib
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
STUD = os.path.dirname(HERE)
_spec = importlib.util.spec_from_file_location(
    "d2j0b_run", os.path.join(STUD, "d2j0b_plane", "run.py"))
P = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(P)
J = P.J                                   # d2j0_hotdisk helpers (run_one, keep_awake, q_net)

BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++-amr")
INPUT = J.INPUT                           # input_hotdisk, unchanged
OUT = os.path.join(HERE, "output")
J.OUT = OUT                               # run_one's makedirs

LXY, LZ = P.LXY, P.LZ
R_DISK = P.R_DISK
H_IN, T_IN = P.H_IN, P.T_IN
Z0, FEED, R_COLL = P.Z0, P.FEED, P.R_COLL
PLOT_S = 50.0
DZ_FINE, DT_FINE = 1e-3, 8e-3
BAND_DEPTH, BAND_HEIGHT = 0.020, 0.004    # m, finest-level band below / above the surface
REGRID_S = 2.0                            # s between regrids (surface moves < 2 mm)
MGS = {"U1": 40, "U2": 20, "F1": 40, "B1": 40}


def binary_sha():
    h = hashlib.sha256()
    with open(BIN, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


def job(name):
    kind, stop = name.split("_")
    stop = float(stop)
    h = f"if(s>0.0,{H_IN!r},{1.0e-3!r})"      # C1: exclusion above the plane
    T = f"{T_IN!r}"
    amr = kind in ("F1", "B1")
    dz0 = 2e-3 if kind in ("U2", "F1", "B1") else DZ_FINE
    dz_fine = 2e-3 if kind == "U2" else DZ_FINE
    dt = 16e-3 if kind == "U2" else DT_FINE
    nxy0, nz0 = int(round(LXY / dz0)), int(round(LZ / dz0))
    nxy, nz = int(round(LXY / dz_fine)), int(round(LZ / dz_fine))
    k = [f"stop_time={stop!r}", f"timestep={dt!r}", f"amr.plot_int={int(round(PLOT_S / dt))}",
         f"amr.n_cell={nxy0} {nxy0} {nz0}", f"geometry.prob_hi={LXY!r} {LXY!r} {LZ!r}",
         f"amr.max_grid_size={MGS[kind]}", f"weibull.V0={dz_fine ** 3!r}", "amr.blocking_factor=1",
         f"surface_patch.radius={R_DISK!r}",
         f'surface_patch.h_expr="{h}"', f'surface_patch.T_flame_expr="{T}"',
         "surface_patch.nozzle_descent=prescribed",
         f"surface_patch.nozzle_z0={Z0!r}", f"surface_patch.nozzle_feed={FEED!r}",
         f"surface_patch.nozzle_collision_radius={R_COLL!r}",
         "energy_ledger.enabled=0"]          # single-level only; not a gate quantity here
    if amr:
        depth = 1.0 if kind == "F1" else BAND_DEPTH          # F1: refine everything
        height = 1.0 if kind == "F1" else BAND_HEIGHT
        k += ["amr.max_level=1", "amr.ref_ratio=2", "amr.nsubsteps=1",
              f"amr.regrid_int={int(round(REGRID_S / dt))}", "amr.n_error_buf=1",
              "amr.grid_eff=0.9",
              f"surface.refine_depth={depth!r}", f"surface.refine_height={height!r}",
              "hc.heat.refinement_threshold=1.0e30"]          # band tagging only
    pf = os.path.join(OUT, name)
    cmd = (["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4", BIN, INPUT,
            f"plot_file={pf}"] + k)
    q = J.q_net()
    # nxy/dz are the FINEST-level values: removal_events.csv columns are in
    # the surface-owning level's index space, so d2j0b_plane/analyze.py loads
    # an AMR run exactly like a uniform fine run.
    meta = dict(name=name, case="C1", kind=kind, mesh="1mm" if dz_fine == 1e-3 else "2mm",
                dz=dz_fine, nxy=nxy, nz=nz, dz_base=dz0, nxy_base=nxy0, nz_base=nz0,
                max_level=1 if amr else 0, band_depth=(depth if amr else 0.0),
                band_height=(height if amr else 0.0), regrid_s=(REGRID_S if amr else 0.0),
                lxy=LXY, lz=LZ, dt=dt, stop=stop, plot_s=PLOT_S, max_grid_size=MGS[kind],
                V0=dz_fine ** 3, patch_r=R_DISK, h_in=H_IN, T_in=T_IN, h_up=1.0e-3, T_up=T_IN,
                z0=Z0, feed=FEED, r_coll=R_COLL, h_expr=h, T_expr=T, q_net=q,
                overshoot=q * dt / (J.RHOCP * dz_fine), binary_sha256=binary_sha())
    return pf, cmd, meta


def run_all(names, jobs, force):
    J.keep_awake()
    js = [job(n) for n in names]
    print(f"{len(js)} runs, {min(jobs, len(js))} in parallel -> {OUT}", flush=True)
    for _, _, m in js:
        print(f"  {m['name']}: base {m['dz_base'] * 1e3:g} mm {m['nxy_base']}²×{m['nz_base']}, "
              f"max_level {m['max_level']}, band {m['band_depth'] * 1e3:g}/{m['band_height'] * 1e3:g} mm, "
              f"dt {m['dt'] * 1e3:g} ms, stop {m['stop']:g} s", flush=True)
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
    a = ap.parse_args()
    if a.list:
        print(f"binary {binary_sha()[:8]}…")
        for n in (a.cases or ["U1_60", "F1_60", "U1_250", "B1_250"]):
            _, cmd, m = job(n)
            print(f"  {n}: {' '.join(cmd[7:])}")
        return
    if not a.cases:
        ap.error("--cases required")
    sys.exit(0 if run_all(a.cases, a.jobs, a.force) else 1)


if __name__ == "__main__":
    main()
