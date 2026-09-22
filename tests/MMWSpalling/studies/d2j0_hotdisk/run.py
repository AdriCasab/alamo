#!/usr/bin/env python3
"""D2j-0 harness: the hot-disk edge test. Standalone — it shares no keys with
the d2c→d2i Meier chain and imports none of it.

  run.py --cases S4 S2 S1 [--jobs 3] [--force]
  run.py --cases T4 T2 T1 [--jobs 3]
  run.py --cases S05                      conditional 0.5 mm leg (~2-3 h)
  run.py --list                           the matrix, with the edge check
  run.py --kill NAME                      binary path + plot_file match

A fixed Robin disk on flat rock, quarter domain 0.08 x 0.08 x 0.10 m with the
axis at the xlo/ylo corner, no jet, no feet, no descent, idle clock off.

  sharp   S4/S2/S1/S05: h_conv = 700, T_flame = 1600 constant inside
          radius 0.030 — h drops 700 -> 0 across one cell.
  tapered T4/T2/T1:     h ramped by a raised cosine from 700 at r = 27.5 mm to
          0 at r = 32.5 mm (patch radius 0.0325 so the ramp is not truncated;
          the parser aborts on h <= 0 inside the patch, so the patch edge sits
          exactly where the ramp reaches zero). T_flame_expr = 1600 constant.
          nozzle_z0 = 0.15 m with no feed keeps s = z_n - z_face >= 0.05 m, so
          the s <= 0 rule and the collision guard never fire; h does not
          depend on s.

dt is fixed at 4 ms at every mesh so dt is not a confound in the sweep (the
explicit limit is dz^2/(6*alpha) = 0.24 s at 1 mm, and the overshoot rule
q*dt/(rho*cp*dz) gives 0.96 K/step at 1 mm, 1.9 K at 0.5 mm). weibull.V0 =
V_cell = dz^3, the scored convention, so it moves with the mesh.

Output: studies/d2j0_hotdisk/output/.
"""
import argparse
import ast
import concurrent.futures as cf
import hashlib
import os
import shutil
import signal
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
INPUT = os.path.join(HERE, "input_hotdisk")
OUT = os.path.join(HERE, "output")

LXY, LZ = 0.08, 0.10
H_CONV, T_FLAME, T_AMB = 700.0, 1600.0, 293.15
R_SHARP = 0.030                      # sharp patch radius
R_IN, R_OUT = 0.0275, 0.0325         # tapered ramp: full h inside, 0 at R_OUT
DT = 0.004                           # fixed at every mesh
STOP = 150.0
PLOT_S = 25.0
RHOCP, KAPPA = 2750.0 * 790.0, 1.5
MESH = {"4": 4e-3, "2": 2e-3, "1": 1e-3, "05": 0.5e-3}
MGS = {4e-3: 16, 2e-3: 20, 1e-3: 40, 0.5e-3: 40}
CASES = ["S4", "S2", "S1", "T4", "T2", "T1", "S05"]

# Raised cosine from 700 at r <= 27.5 mm to 0 at r >= 32.5 mm. r is clamped
# before the cosine, and nothing is divided by a max() (the AMReX parser
# evaluates a/max(x,num) as 0).
H_EXPR = (f"{H_CONV!r}*0.5*(1.0+cos(3.141592653589793*"
          f"(min({R_OUT!r},max({R_IN!r},r))-{R_IN!r})/{R_OUT - R_IN!r}))")


def h_of(r):
    """The same ramp in numpy, for the pre-flight edge check."""
    rc = np.clip(r, R_IN, R_OUT)
    return H_CONV * 0.5 * (1.0 + np.cos(np.pi * (rc - R_IN) / (R_OUT - R_IN)))


def binary_sha():
    h = hashlib.sha256()
    with open(BIN, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def q_net(Ts=821.0):
    """Net flux into a pinned surface at Ts: Robin minus radiation and ambient."""
    return (H_CONV * (T_FLAME - Ts) - 0.8 * 5.670374e-8 * (Ts ** 4 - T_AMB ** 4)
            - 10.0 * (Ts - T_AMB))


def job(name):
    tapered = name.startswith("T")
    dz = MESH[name[1:]]
    nxy, nz = int(round(LXY / dz)), int(round(LZ / dz))
    q = q_net()
    k = [f"stop_time={STOP!r}", f"timestep={DT!r}", f"amr.plot_int={int(round(PLOT_S / DT))}",
         f"amr.n_cell={nxy} {nxy} {nz}", f"geometry.prob_hi={LXY!r} {LXY!r} {LZ!r}",
         f"amr.max_grid_size={MGS[dz]}", f"weibull.V0={dz ** 3!r}",
         # 0.10 m at 4 mm is 25 cells, which is not divisible by the input's
         # blocking_factor 2. Set 1 at EVERY mesh rather than vary it across
         # the sweep or move the packet's geometry. Grid layout does not
         # change results (the removed_mf ghost-cell fix); checked by running
         # S2 at both and comparing (RESULTS §1).
         "amr.blocking_factor=1"]
    if tapered:
        k += [f"surface_patch.radius={R_OUT!r}", f'surface_patch.h_expr="{H_EXPR}"',
              f'surface_patch.T_flame_expr="{T_FLAME!r}"',
              "surface_patch.nozzle_z0=0.15", "surface_patch.nozzle_feed=0.0"]
    else:
        k += [f"surface_patch.radius={R_SHARP!r}", f"surface_patch.h_conv={H_CONV!r}",
              f"surface_patch.T_flame={T_FLAME!r}"]
    pf = os.path.join(OUT, name)
    cmd = (["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4", BIN, INPUT,
            f"plot_file={pf}"] + k)
    meta = dict(name=name, set="tapered" if tapered else "sharp", dz=dz, nxy=nxy, nz=nz,
                lxy=LXY, lz=LZ, dt=DT, stop=STOP, plot_s=PLOT_S, max_grid_size=MGS[dz],
                V0=dz ** 3, patch_r=R_OUT if tapered else R_SHARP, h_conv=H_CONV,
                T_flame=T_FLAME, q_net=q, overshoot=q * DT / (RHOCP * dz),
                dt_explicit=dz * dz / (6.0 * KAPPA / RHOCP), binary_sha256=binary_sha())
    return pf, cmd, meta


def edge_check(name):
    """The tapered patch aborts if any in-patch cell centre has h <= 0, which
    happens only if a centre lands exactly on R_OUT. Report the margin."""
    _, _, m = job(name)
    dz = m["dz"]
    c = (np.arange(m["nxy"]) + 0.5) * (m["lxy"] / m["nxy"])
    X, Y = np.meshgrid(c, c, indexing="ij")
    r = np.hypot(X, Y)
    inp = r <= m["patch_r"]
    if m["set"] == "sharp":
        return f"{int(inp.sum())} in-patch columns, h = {H_CONV:g} flat"
    hmin = float(h_of(r[inp]).min())
    return (f"{int(inp.sum())} in-patch columns, min h = {hmin:.3e} W/m²K "
            f"(must be > 0), closest centre to R_OUT: {abs(r[inp] - m['patch_r']).min() * 1e3:.4f} mm")


def run_one(pf, cmd, meta, force):
    done = pf + ".done"
    if os.path.exists(done) and not force:
        return meta, "skip (done)", 0.0
    for p in (pf, pf + "_removal_events.csv", pf + "_h_col_events.csv", done):
        if os.path.isdir(p):
            shutil.rmtree(p)
        elif os.path.exists(p):
            os.remove(p)
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    with open(pf + ".log", "w") as log:
        log.write(" ".join(cmd) + "\n")
        log.flush()
        proc = subprocess.Popen(cmd, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT,
                                start_new_session=True)
        rc = proc.wait()
    wall = time.time() - t0
    t_end = 0.0
    th = os.path.join(pf, "thermo.dat")
    if os.path.exists(th):
        a = np.loadtxt(th, skiprows=1, ndmin=2)
        if a.size:
            t_end = float(a[-1, 0])
    meta = dict(meta, wall_s=wall, status="ok" if rc == 0 else f"FAILED rc={rc}", t_end=t_end)
    if rc == 0:
        with open(done, "w") as f:
            f.write(repr(meta) + "\n")
    return meta, meta["status"], wall


def keep_awake():
    """Hold a no-sleep assertion for as long as this harness lives, so a run
    is never lost to the laptop sleeping. `caffeinate -w <pid>` exits on its
    own when we do, so there is nothing to clean up and no effect if the
    caller already wrapped us in one."""
    try:
        return subprocess.Popen(["caffeinate", "-dis", "-w", str(os.getpid())],
                                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except (OSError, ValueError) as e:
        print(f"  (no caffeinate: {e}; the host may sleep mid-run)", flush=True)
        return None


def run_all(names, jobs, force):
    keep_awake()
    js = [job(n) for n in names]
    print(f"{len(js)} runs, {min(jobs, len(js))} in parallel -> {OUT}", flush=True)
    for _, _, m in js:
        print(f"  {m['name']}: {m['set']}, dz {m['dz'] * 1e3:g} mm, {m['nxy']}²×{m['nz']}, "
              f"dt {m['dt'] * 1e3:g} ms ({m['overshoot']:.2f} K/step), stop {m['stop']:g} s, "
              f"V0 {m['V0']:.3g} m³, mgs {m['max_grid_size']}", flush=True)
    ok = True
    with cf.ThreadPoolExecutor(max_workers=min(jobs, len(js))) as ex:
        for fu in cf.as_completed([ex.submit(run_one, pf, c, m, force) for pf, c, m in js]):
            m, status, wall = fu.result()
            print(f"  {m['name']:4s} {status} at t = {m.get('t_end', 0):g} s "
                  f"({wall / 60:.1f} min)", flush=True)
            ok = ok and not status.startswith("FAILED")
    return ok


def meta_of(pf):
    return ast.literal_eval(open(pf + ".done").read())


def kill(name):
    """Match the binary path AND plot_file=<prefix>; a bare-name match would
    also kill shell waiters (the D2d lesson)."""
    tok = f"plot_file={os.path.join(OUT, name)}"
    ps = subprocess.run(["ps", "-axo", "pid=,command="], capture_output=True, text=True).stdout
    pids = []
    for line in ps.splitlines():
        pid, _, cmd = line.strip().partition(" ")
        argv = cmd.split()
        if BIN in argv and tok in argv and int(pid) != os.getpid():
            pids.append(int(pid))
    for p in pids:
        try:
            os.kill(p, signal.SIGTERM)
        except ProcessLookupError:
            pass
    print(f"SIGTERM to {len(pids)} processes of {name}: {pids}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", nargs="*")
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--kill")
    a = ap.parse_args()
    if a.kill:
        return kill(a.kill)
    if a.list:
        print(f"binary {binary_sha()[:8]}…  q_net(821 K) = {q_net() / 1e6:.3f} MW/m², "
              f"v = {q_net() / (RHOCP * (821.0 - T_AMB)) * 3600:.3f} m/h, "
              f"thermal length kappa/(rho cp v) = "
              f"{KAPPA / (RHOCP * q_net() / (RHOCP * (821.0 - T_AMB))) * 1e3:.2f} mm")
        for n in CASES:
            _, _, m = job(n)
            print(f"{n:4s} {m['set']:7s} dz {m['dz'] * 1e3:4g} mm {m['nxy']:3d}²×{m['nz']:<4d} "
                  f"dt {m['dt'] * 1e3:g} ms {m['overshoot']:.2f} K/step (explicit limit "
                  f"{m['dt_explicit']:.3f} s) V0 {m['V0']:.3g} mgs {m['max_grid_size']}")
            print("     ", edge_check(n))
        return
    bad = [n for n in (a.cases or []) if n not in CASES]
    if not a.cases or bad:
        sys.exit(f"give --cases from {CASES}" + (f" (unknown: {bad})" if bad else ""))
    sys.exit(0 if run_all(a.cases, a.jobs, a.force) else 1)


if __name__ == "__main__":
    main()
