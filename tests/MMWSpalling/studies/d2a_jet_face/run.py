#!/usr/bin/env python3
"""D2a study: impinging-jet face source (surface_patch.h_expr / T_flame_expr).

Study, not a regression (no `test`). CLI overrides only; no input is copied or
edited. Profiles and the Martin h_ref come from jet.py.

  run.py --part 1          2-D slab machinery (input_2d_dev, 4 ranks)  -> output/p1_*
  run.py --probe           3-D cost probe (10 s of J-M at 1.5 m/h)       -> output/probe
  run.py --part 2          3-D quarter domain (input_drilling, 4 ranks) -> output/p2_*
  run.py --part 1 --cases JM --force --jobs 2

Common: pinned closure, follow_mask, energy ledger, removal_events.csv,
weibull.V0 = V_cell (user decision), beam off. Time step: the overshoot rule
q dt/(rho Cp dz) <= 5 K with q the maximum stagnation closed-form flux (r = 0,
s over [2, 12] D, core T_gas, T_fire = 821 K), dt = DT_CAP halved until it
holds. No parameter is adjusted toward Meier's ROP or hole diameter.
"""
import argparse
import concurrent.futures as cf
import math
import os
import shutil
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import jet  # noqa: E402

ROOT = HERE.split(os.sep + os.path.join("tests", "MMWSpalling") + os.sep)[0]
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
MEIER = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "meier", "sp_meier_pilot")
OUT = os.path.join(HERE, "output")
RHOCP = 2750.0 * 790.0
SIG, EPS, HC, TA = 5.670374e-8, 0.8, 10.0, 293.15
D = jet.D
H_MARTIN = round(jet.martin_reference()["h_loc"])     # local h at 2.5 D, SOD 7 D [W/m^2K]
ANCHORS = {"JM": (float(H_MARTIN), 1500.0), "J5": (5.0e3, 1200.0), "J10": (1.0e4, 1000.0)}
DT_CAP = 4.0e-3
FEED_MH = {"f0": 0.0, "f15": 1.5, "f30": 3.0}


def q_closed(h, T_gas, T_f=jet.T_FIRE):
    return h * (T_gas - T_f) - EPS * SIG * (T_f ** 4 - TA ** 4) - HC * (T_f - TA)


def q_stag_max(h_ref, T_ref, T_ent, stag):
    s = np.linspace(2.0 * D, 12.0 * D, 101)
    return float(np.max(q_closed(jet.h_py(h_ref, 0.0, s, stag), jet.Tg_py(T_ref, s, T_ent))))


def dt_rule(q, dz, cap=DT_CAP):
    dt = cap
    while q * dt / (RHOCP * dz) > 5.0:
        dt *= 0.5
    return dt


def flame_keys(anchor, T_ent, stag, z_top, feed, extra_r_coll=None):
    h_ref, T_ref = ANCHORS[anchor]
    k = ["surface_patch.enabled=1", "surface_patch.mode=convective_flame",
         "surface_patch.robin_form=pinned", "surface_patch.face=zhi",
         f'surface_patch.h_expr="{jet.h_expr(h_ref, stag)}"',
         f'surface_patch.T_flame_expr="{jet.Tg_expr(T_ref, T_ent)}"',
         f"surface_patch.nozzle_z0={z_top + 7.0 * D!r}", f"surface_patch.nozzle_feed={feed!r}",
         "surface.follow_mask=1", "energy_ledger.enabled=1", "spall.removal_events_csv=1",
         "beam.P0=0.0", "amr.thermo.int=1"]
    if extra_r_coll is not None:
        k.append(f"surface_patch.nozzle_collision_radius={extra_r_coll!r}")
    return k


# ---------------------------------------------------------------- Part 1: 2-D slab
# name: (anchor, dz_mm, T_ent, stagnation, idle, dt_cap)
PART1 = {
    "JM":          ("JM",  2.0, jet.T_AMB, "flat", 2.0, DT_CAP),
    "J10":         ("J10", 2.0, jet.T_AMB, "flat", 2.0, DT_CAP),
    "JM_1mm":      ("JM",  1.0, jet.T_AMB, "flat", 2.0, DT_CAP),
    "J10_idle1":   ("J10", 2.0, jet.T_AMB, "flat", 1.0, DT_CAP),
    "JM_Tent600":  ("JM",  2.0, 600.0,     "flat", 2.0, DT_CAP),
    "JM_stag":     ("JM",  2.0, jet.T_AMB, "rise", 2.0, DT_CAP),
    "JM_dt1ms":    ("JM",  2.0, jet.T_AMB, "flat", 2.0, 1.0e-3),   # dt check for the 4 ms cap
}
STOP1 = 60.0


def part1_job(name):
    anchor, dz_mm, T_ent, stag, idle, cap = PART1[name]
    h_ref, T_ref = ANCHORS[anchor]
    dz = dz_mm * 1e-3
    q = q_stag_max(h_ref, T_ref, T_ent, stag)
    dt = dt_rule(q, dz, cap)
    pf = os.path.join(OUT, f"p1_{name}")
    ov = [f"plot_file={pf}", f"stop_time={STOP1}", f"timestep={dt!r}",
          f"amr.plot_int={int(round(10.0 / dt))}",
          f"amr.thermo.plot_int={max(1, int(round(0.02 / dt)))}",
          "surface_patch.x0=0.07", "surface_patch.y0=0.004", "surface_patch.radius=0.07",
          f"surface_patch.pinned_idle_cycles={idle!r}",
          f"weibull.V0={dz ** 3!r}"] + flame_keys(anchor, T_ent, stag, 0.120, 0.0)
    if dz_mm != 2.0:
        ov.append(f"amr.n_cell={int(round(0.14 / dz))} {int(round(0.008 / dz))} {int(round(0.120 / dz))}")
    cmd = ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4", BIN,
           os.path.join(MEIER, "input_2d_dev")] + ov
    meta = dict(part=1, name=name, anchor=anchor, h_ref=h_ref, T_ref=T_ref, T_ent=T_ent, stag=stag,
                idle=idle, dz=dz_mm, dt=dt, stop=STOP1, q_stag=q, overshoot=q * dt / (RHOCP * dz),
                feed=0.0, z_top=0.120, x0=0.07, y0=0.004, dim=2, nranks=4)
    return pf, cmd, meta


# ---------------------------------------------------------------- Part 2: 3-D quarter
# name: (anchor, feed key)
PART2 = {
    "JM_f15":  ("JM", "f15"),
    "J5_f15":  ("J5", "f15"),
    "J10_f15": ("J10", "f15"),
    "JM_f0":   ("JM", "f0"),
    "JM_f30":  ("JM", "f30"),
}
Z_TOP3, BOTTOM_MARGIN, STOP2_MAX = 0.20, 0.040, 300.0
R_COLL = 2.5 * D        # nozzle_collision_radius (user decision 2026-09-15)


def stop2(anchor, feed):
    """Stop at 300 s or when the face would reach 40 mm above the bottom. For
    feed > 0 the face leads the nozzle path by at most (s_lead - 7 D); s_lead
    is bounded by the stand-off where the stagnation closed form drops to
    rho Cp (T_fire - T0) feed (T_gas(s) 1/s tail, h at the 12 D clamp)."""
    if feed <= 0.0:
        return STOP2_MAX
    h_ref, T_ref = ANCHORS[anchor]
    s = np.linspace(7.0 * D, 40.0 * D, 4001)
    q = q_closed(jet.h_py(h_ref, 0.0, s), jet.Tg_py(T_ref, s))
    q_need = RHOCP * (jet.T_FIRE - TA) * feed
    s_lead = float(s[np.argmax(q < q_need)]) if np.any(q < q_need) else float(s[-1])
    lead = max(0.0, s_lead - 7.0 * D)
    t = (Z_TOP3 - BOTTOM_MARGIN - lead) / feed
    return float(min(STOP2_MAX, math.floor(t)))


def part2_job(name, stop_override=None):
    anchor, fk = PART2[name]
    h_ref, T_ref = ANCHORS[anchor]
    feed = FEED_MH[fk] / 3600.0
    dz = 2.0e-3
    q = q_stag_max(h_ref, T_ref, jet.T_AMB, "flat")
    dt = dt_rule(q, dz)
    stop = stop_override if stop_override is not None else stop2(anchor, feed)
    pf = os.path.join(OUT, "probe" if stop_override is not None else f"p2_{name}")
    ov = [f"plot_file={pf}", f"stop_time={stop!r}", f"timestep={dt!r}",
          f"amr.plot_int={int(round(100.0 / dt))}",
          f"amr.thermo.plot_int={max(1, int(round(0.1 / dt)))}",
          "amr.n_cell=60 60 100", "geometry.prob_hi=0.12 0.12 0.20",
          "bit.enabled=0",
          "surface_patch.x0=0.0", "surface_patch.y0=0.0", "surface_patch.radius=0.2",
          "surface_patch.pinned_idle_cycles=2.0",
          f"weibull.V0={dz ** 3!r}"] + flame_keys(anchor, jet.T_AMB, "flat", Z_TOP3, feed, R_COLL)
    cmd = ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4", BIN,
           os.path.join(MEIER, "input_drilling")] + ov
    meta = dict(part=2, name=name, anchor=anchor, h_ref=h_ref, T_ref=T_ref, T_ent=jet.T_AMB, stag="flat",
                idle=2.0, dz=2.0, dt=dt, stop=stop, q_stag=q, overshoot=q * dt / (RHOCP * dz),
                feed=feed, z_top=Z_TOP3, x0=0.0, y0=0.0, dim=3, nranks=4, r_coll=R_COLL)
    return pf, cmd, meta


# ---------------------------------------------------------------- parser check
# One step on the 3-D quarter geometry, feed 0, nozzle at z_top + s0: every
# column starts at T = 293.15 K with no pin, so each takes the face form.
# patch_P_robin (row t = dt) must equal sum_cols h(T_gas - T_f) dx dy with h,
# T_gas from the numpy forms in jet.py and T_f from the face balance, to the
# 6-digit thermo print. Guards against parser mis-evaluation (see jet.py).
# name: (anchor, T_ent, stagnation, s0 / D)
PCHECK = {
    "JM_s1.5": ("JM", jet.T_AMB, "flat", 1.5), "JM_s3": ("JM", jet.T_AMB, "flat", 3.0),
    "JM_s5": ("JM", jet.T_AMB, "flat", 5.0), "JM_s7": ("JM", jet.T_AMB, "flat", 7.0),
    "JM_s10": ("JM", jet.T_AMB, "flat", 10.0), "JM_s13": ("JM", jet.T_AMB, "flat", 13.0),
    "JM_rise_s7": ("JM", jet.T_AMB, "rise", 7.0), "JM_Tent600_s9": ("JM", 600.0, "flat", 9.0),
    "J5_s7": ("J5", jet.T_AMB, "flat", 7.0), "J10_s7": ("J10", jet.T_AMB, "flat", 7.0),
}
KAPPA = 1.5


def face_P_robin(anchor, T_ent, stag, s0, n=60, dx=2.0e-3, dz=2.0e-3):
    h_ref, T_ref = ANCHORS[anchor]
    xc = (np.arange(n) + 0.5) * dx
    X, Y = np.meshgrid(xc, xc, indexing="ij")
    r = np.hypot(X, Y)
    h = jet.h_py(h_ref, r, np.full_like(r, s0), stag)
    Tg = jet.Tg_py(T_ref, np.full_like(r, s0), T_ent)
    G, Tc = 2.0 * KAPPA / dz, TA
    lo, hi = np.full_like(r, TA), np.maximum(Tg, TA)
    for _ in range(200):                   # bisection on the monotone face residual
        Tf = 0.5 * (lo + hi)
        F = G * (Tf - Tc) - h * (Tg - Tf) + EPS * SIG * (Tf ** 4 - TA ** 4) + HC * (Tf - TA)
        hi = np.where(F > 0.0, Tf, hi)
        lo = np.where(F > 0.0, lo, Tf)
    Tf = 0.5 * (lo + hi)
    return float(np.sum(h * (Tg - Tf)) * dx * dx)


def pcheck_job(name):
    anchor, T_ent, stag, sD = PCHECK[name]
    dt = 1.0e-3
    pf = os.path.join(OUT, "pcheck", name)
    ov = [f"plot_file={pf}", f"stop_time={2.0 * dt!r}", f"timestep={dt!r}", "amr.plot_int=100000",
          "amr.thermo.plot_int=1", "amr.n_cell=60 60 100", "geometry.prob_hi=0.12 0.12 0.20",
          "bit.enabled=0", "surface_patch.x0=0.0", "surface_patch.y0=0.0", "surface_patch.radius=0.2",
          "surface_patch.pinned_idle_cycles=2.0", "weibull.V0=8e-09"] + \
        flame_keys(anchor, T_ent, stag, Z_TOP3 - 7.0 * D + sD * D, 0.0, R_COLL)
    cmd = ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4", BIN,
           os.path.join(MEIER, "input_drilling")] + ov
    return pf, cmd, dict(name=name, anchor=anchor, T_ent=T_ent, stag=stag, s0=sD * D)


def parser_check():
    os.makedirs(os.path.join(OUT, "pcheck"), exist_ok=True)
    ok = True
    for name in PCHECK:
        pf, cmd, m = pcheck_job(name)
        m, status, _ = run_one(pf, cmd, m, True)
        if status != "ok":
            print(f"  {name}: {status}")
            ok = False
            continue
        rows = [l.split() for l in open(os.path.join(pf, "thermo.dat")) if l.strip()]
        P_alamo = float(rows[2][rows[0].index("patch_P_robin")])        # row t = dt
        s_alamo = float(rows[2][rows[0].index("patch_min_standoff")])
        P_py = face_P_robin(m["anchor"], m["T_ent"], m["stag"], m["s0"])
        rel = abs(P_alamo / P_py - 1.0)
        good = rel <= 2.0e-5 and abs(s_alamo / m["s0"] - 1.0) <= 1.0e-5
        ok = ok and good
        print(f"  [{'PASS' if good else 'FAIL'}] {name}: patch_P_robin ALAMO {P_alamo:.6g} W, "
              f"numpy {P_py:.6g} W, rel {rel:.1e}; s = {s_alamo:g} m")
    print("PARSER CHECK PASS" if ok else "PARSER CHECK FAIL")
    return ok


def run_one(pf, cmd, meta, force):
    done = pf + ".done"
    if os.path.exists(done) and not force:
        return meta, "skip (done)", 0.0
    for p in (pf, pf + "_h_col_events.csv", pf + "_removal_events.csv", done):
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
    meta = dict(meta, wall=wall)
    if rc == 0:
        with open(done, "w") as f:
            f.write(repr(meta) + "\n")
        return meta, "ok", wall
    return meta, f"FAILED rc={rc}", wall


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--part", type=int, choices=(1, 2))
    ap.add_argument("--probe", action="store_true")
    ap.add_argument("--parser-check", action="store_true")
    ap.add_argument("--cases", nargs="*", default=None)
    ap.add_argument("--jobs", type=int, default=None)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    if not os.path.exists(BIN):
        sys.exit(f"missing binary {BIN}; build first")
    print(f"h_Martin (local, 2.5 D, SOD 7 D) = {H_MARTIN} W/m^2K")
    if args.parser_check:
        sys.exit(0 if parser_check() else 1)
    if args.probe:
        jobs = [part2_job("JM_f15", stop_override=10.0)]
        nj = 1
    elif args.part == 1:
        jobs = [part1_job(n) for n in (args.cases or PART1)]
        nj = args.jobs or 3
    elif args.part == 2:
        jobs = [part2_job(n) for n in (args.cases or PART2)]
        nj = args.jobs or 3
    else:
        sys.exit("give --part 1|2 or --probe")
    print(f"{len(jobs)} runs, {nj} in parallel -> {OUT}")
    for pf, _, m in jobs:
        print(f"  {os.path.basename(pf)}: dt = {m['dt'] * 1e3:g} ms, stop = {m['stop']} s, "
              f"stagnation q_max = {m['q_stag'] / 1e6:.3f} MW/m^2, overshoot {m['overshoot']:.2f} K/step")
    ok = True
    with cf.ThreadPoolExecutor(max_workers=nj) as ex:
        for fu in cf.as_completed([ex.submit(run_one, pf, c, m, args.force) for pf, c, m in jobs]):
            m, status, wall = fu.result()
            print(f"  {m['name']:12s} {status}  ({wall:.1f} s)", flush=True)
            if args.probe and status == "ok":
                steps = m["stop"] / m["dt"]
                cells = 60 * 60 * 100
                for n in PART2:
                    _, _, mm = part2_job(n)
                    est = wall / steps * (mm["stop"] / mm["dt"])
                    print(f"    estimate {n}: stop {mm['stop']} s, dt {mm['dt'] * 1e3:g} ms -> "
                          f"{est / 60:.1f} min ({cells * steps / wall:.3g} cell-steps/s in probe)")
            ok = ok and status.startswith(("ok", "skip"))
    print("RUNS OK" if ok else "RUNS HAD FAILURES")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
