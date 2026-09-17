#!/usr/bin/env python3
"""D2a2 study: wall-jet enthalpy closure (surface_patch.jet_closure = enthalpy).

Study, not a regression (no `test`). CLI overrides only; no input is copied or
edited. h(r, s) and the Martin h_ref come from walljet.py (D2a's jet.py at
D = 7.5 mm, SOD = 50 mm).

  run.py --part 1          2-D slab machinery (input_2d_dev, 4 ranks)  -> output/p1_*
  run.py --probe           3-D cost probe (10 s of J10_19, dt 2 ms)     -> output/probe
  run.py --part 2          3-D quarter domain (input_drilling, 4 ranks) -> output/p2_*
  run.py --part 1 --cases JM_19_f0 --force --jobs 2

Anchors: h_ref in {Martin (film at (T_nozzle + T_fire)/2), 5e3, 1e4} x
T_nozzle in {1900 K adiabatic, 1436 K chamber TC}: JM_19 JM_14 J5_19 J5_14
J10_19 J10_14.

mdot: the jet march must see the face area its mdot feeds.
  Part 2 quarter domain (axis at the xlo/ylo corner): mdot/4.
  Part 1 slab (8 mm thick through the axis, both sides): the strip area per dr
  is 2 w dr against the annulus 2 pi r dr, so mdot * w/(pi r_core) (the slab's
  share at the core edge). Beyond r_core the slab keeps too much gas (the
  annulus grows, the strip does not): Part 1 is machinery and trend only.

Stagnation: jet_stagnation = decay (baseline; T_stag = T_ent + (T_nozzle -
T_ent) min(1, 5 D/s_c)) with nozzle mass; treatments: jet_entrained_mass = 1,
jet_stagnation = nozzle.
Stop: 300 s (Part 2) / 60 s (Part 1), or when the centre face could reach
40 mm above the domain bottom:
  decay:  the centre leads the nozzle by at most s_eq + 10 mm, s_eq the
          stand-off where the centre's closed-form ROP equals the feed (feed 0:
          where T_stag falls to T_fire, so the centre stalls);
  nozzle: the centre outruns any feed; stop at its fastest closed-form rate.
Time step: q dt/(rho Cp dz) <= 5 K with q the stagnation maximum over the
stand-offs the centre can see (s >= SOD - 2 mm if the centre outruns the feed
at SOD, else s >= 5 D); dt = 4 ms halved until it holds.
No parameter is adjusted toward Meier's ROP, hole diameter or volume rate.
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
sys.path.insert(0, HERE)
import walljet as wj  # noqa: E402

ROOT = HERE.split(os.sep + os.path.join("tests", "MMWSpalling") + os.sep)[0]
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
MEIER = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "meier", "sp_meier_pilot")
OUT = os.path.join(HERE, "output")
D = wj.D
DT_CAP = 4.0e-3
FEED = 1.5 / 3600.0
R_COLL = 2.5 * D
SLAB_W = 0.008
ANCHORS = {f"J{h}_{t}": (h, t) for h in ("M", "5", "10") for t in ("19", "14")}


def anchor(name):
    hk, tk = ANCHORS[name]
    Tn = wj.T_NOZZLES[tk]
    return float(wj.h_anchor(hk, Tn)), Tn


def dt_rule(q, dz, cap=DT_CAP):
    dt = cap
    while q * dt / (wj.RHOCP * dz) > 5.0:
        dt *= 0.5
    return dt


S_MARGIN = 0.010


def s_lead(h, Tn, feed, stag):
    """Largest centre stand-off (decay) or None (nozzle: runaway)."""
    if stag == "nozzle":
        return None
    if feed <= 0.0:
        phi = (wj.T_FIRE - wj.TA) / (Tn - wj.TA)
        return wj.CORE_LEN * D / phi
    se = wj.s_equilibrium(h, Tn, feed * 3600.0)
    return max(se, wj.SOD)


def stop_rule(name, z_top, feed, z_floor, t_max, stag):
    h, Tn = anchor(name)
    sl = s_lead(h, Tn, feed, stag)
    if sl is None:
        rop = wj.centre_rop_max(h, Tn, "nozzle") / 3600.0
        return float(min(t_max, math.floor((z_top - z_floor) / rop)))
    room = z_top + wj.SOD - (sl + S_MARGIN) - z_floor
    if feed <= 0.0:
        return t_max if room >= 0.0 else 0.0
    return float(min(t_max, math.floor(room / feed)))


def q_max(name, feed, stag):
    h, Tn = anchor(name)
    if stag == "nozzle":
        return wj.q_stag_max(h, Tn, "nozzle")
    s_min = wj.SOD - 0.002 if wj.centre_rop(h, Tn, wj.SOD) > feed * 3600.0 else wj.CORE_LEN * D
    return wj.q_stag_max(h, Tn, "decay", s_min)


def flame_keys(h_ref, Tn, z_top, feed, mdot, cp=wj.CP, T_ent=wj.TA, jet_dr=None, prof=None,
               stag="decay", entrained=0):
    k = ["surface_patch.enabled=1", "surface_patch.mode=convective_flame",
         "surface_patch.robin_form=pinned", "surface_patch.face=zhi",
         "surface_patch.pinned_idle_cycles=2.0",
         f'surface_patch.h_expr="{wj.h_expr(h_ref)}"',
         f"surface_patch.nozzle_z0={z_top + wj.SOD!r}", f"surface_patch.nozzle_feed={feed!r}",
         f"surface_patch.nozzle_collision_radius={R_COLL!r}",
         "surface_patch.jet_closure=enthalpy", f"surface_patch.jet_mdot={mdot!r}",
         f"surface_patch.jet_cp={cp!r}", f"surface_patch.jet_T_nozzle={Tn!r}",
         f"surface_patch.jet_T_ent={T_ent!r}", f"surface_patch.jet_D={D!r}",
         f"surface_patch.jet_stagnation={stag}", f"surface_patch.jet_entrained_mass={entrained}",
         "surface.follow_mask=1", "energy_ledger.enabled=1", "spall.removal_events_csv=1",
         "beam.P0=0.0", "amr.thermo.int=1"]
    if jet_dr is not None:
        k.append(f"surface_patch.jet_dr={jet_dr!r}")
    if prof is not None:
        k.append(f"surface_patch.jet_profile_interval={prof!r}")
    return k


def mpirun(inp, ov):
    return ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4", BIN, inp] + ov


# ---------------------------------------------------------------- Part 1: 2-D slab
Z_TOP1, STOP1, BUDGET1 = 0.120, 60.0, 0.080
MDOT_SLAB = wj.MDOT * SLAB_W / (math.pi * wj.R_CORE)
# name -> (anchor, feed, cp factor, T_ent, jet_dr, stagnation, entrained)
PART1 = {f"{a}_{fk}": (a, fv, 1.0, wj.TA, None, "decay", 0)
         for a in ANCHORS for fk, fv in (("f0", 0.0), ("f15", FEED))}
PART1.update({
    "JM_19_f15_cp092":   ("JM_19", FEED, 0.92, wj.TA, None, "decay", 0),
    "JM_19_f15_cp108":   ("JM_19", FEED, 1.08, wj.TA, None, "decay", 0),
    "JM_19_f15_Tent600": ("JM_19", FEED, 1.0, 600.0, None, "decay", 0),
    "JM_19_f15_dr1":     ("JM_19", FEED, 1.0, wj.TA, 1.0e-3, "decay", 0),
    "JM_19_f15_noz":     ("JM_19", FEED, 1.0, wj.TA, None, "nozzle", 0),
    "JM_19_f15_ent":     ("JM_19", FEED, 1.0, wj.TA, None, "decay", 1),
})


def part1_job(name):
    a, feed, cpf, T_ent, jdr, stag, ent = PART1[name]
    h, Tn = anchor(a)
    dz = 2.0e-3
    q = q_max(a, feed, stag)
    dt = dt_rule(q, dz)
    stop = stop_rule(a, Z_TOP1, feed, Z_TOP1 - BUDGET1, STOP1, stag)
    pf = os.path.join(OUT, f"p1_{name}")
    ov = [f"plot_file={pf}", f"stop_time={stop!r}", f"timestep={dt!r}",
          f"amr.plot_int={int(round(10.0 / dt))}",
          f"amr.thermo.plot_int={max(1, int(round(0.02 / dt)))}",
          "surface_patch.x0=0.07", "surface_patch.y0=0.004", "surface_patch.radius=0.07",
          f"weibull.V0={dz ** 3!r}"] + \
        flame_keys(h, Tn, Z_TOP1, feed, MDOT_SLAB, wj.CP * cpf, T_ent, jdr, 5.0, stag, ent)
    meta = dict(part=1, name=name, anchor=a, h_ref=h, T_nozzle=Tn, feed=feed, cp=wj.CP * cpf,
                stag=stag, entrained=ent,
                T_ent=T_ent, jet_dr=jdr or dz, mdot=MDOT_SLAB, dz=2.0, dt=dt, stop=stop, q_stag=q,
                overshoot=q * dt / (wj.RHOCP * dz), z_top=Z_TOP1, x0=0.07, y0=0.004, dim=2,
                sector=None)
    return pf, mpirun(os.path.join(MEIER, "input_2d_dev"), ov), meta


# ---------------------------------------------------------------- Part 2: 3-D quarter
Z_TOP2, STOP2, BUDGET2 = 0.20, 300.0, 0.160
# name -> (anchor, stagnation, entrained)
PART2 = {a: (a, "decay", 0) for a in ANCHORS}
PART2.update({"J5_19_ent": ("J5_19", "decay", 1), "J10_19_ent": ("J10_19", "decay", 1),
              "J5_19_noz": ("J5_19", "nozzle", 0)})


def part2_job(name, stop_override=None):
    a, stag, ent = PART2[name]
    h, Tn = anchor(a)
    dz = 2.0e-3
    q = q_max(a, FEED, stag)
    dt = dt_rule(q, dz)
    stop = stop_override if stop_override is not None else \
        stop_rule(a, Z_TOP2, FEED, Z_TOP2 - BUDGET2, STOP2, stag)
    pf = os.path.join(OUT, "probe" if stop_override is not None else f"p2_{name}")
    ov = [f"plot_file={pf}", f"stop_time={stop!r}", f"timestep={dt!r}",
          f"amr.plot_int={int(round(20.0 / dt))}",
          f"amr.thermo.plot_int={max(1, int(round(0.1 / dt)))}",
          "amr.n_cell=60 60 100", "geometry.prob_hi=0.12 0.12 0.20", "bit.enabled=0",
          "surface_patch.x0=0.0", "surface_patch.y0=0.0", "surface_patch.radius=0.2",
          f"weibull.V0={dz ** 3!r}"] + \
        flame_keys(h, Tn, Z_TOP2, FEED, wj.MDOT / 4.0, prof=5.0, stag=stag, entrained=ent)
    meta = dict(part=2, name=name, anchor=a, h_ref=h, T_nozzle=Tn, feed=FEED, cp=wj.CP,
                stag=stag, entrained=ent,
                T_ent=wj.TA, jet_dr=dz, mdot=wj.MDOT / 4.0, dz=2.0, dt=dt, stop=stop, q_stag=q,
                overshoot=q * dt / (wj.RHOCP * dz), z_top=Z_TOP2, x0=0.0, y0=0.0, dim=3,
                sector=0.25)
    return pf, mpirun(os.path.join(MEIER, "input_drilling"), ov), meta


def run_one(pf, cmd, meta, force):
    done = pf + ".done"
    if os.path.exists(done) and not force:
        return meta, "skip (done)", 0.0
    for p in (pf, pf + "_h_col_events.csv", pf + "_removal_events.csv",
              pf + "_jet_profile.csv", done):
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
    ap.add_argument("--cases", nargs="*", default=None)
    ap.add_argument("--jobs", type=int, default=None)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    if not os.path.exists(BIN):
        sys.exit(f"missing binary {BIN}; build first")
    if args.probe:
        jobs, nj = [part2_job("J10_19", stop_override=10.0)], 1
    elif args.part == 1:
        jobs, nj = [part1_job(n) for n in (args.cases or PART1)], args.jobs or 3
    elif args.part == 2:
        jobs, nj = [part2_job(n) for n in (args.cases or PART2)], args.jobs or 3
    else:
        sys.exit("give --part 1|2 or --probe")
    print(f"{len(jobs)} runs, {nj} in parallel -> {OUT}")
    for pf, _, m in jobs:
        print(f"  {os.path.basename(pf)}: h_ref {m['h_ref']:g}, T_nozzle {m['T_nozzle']:g}, "
              f"mdot {m['mdot']:.4g}, dt {m['dt'] * 1e3:g} ms, stop {m['stop']} s, "
              f"q_stag {m['q_stag'] / 1e6:.2f} MW/m^2, overshoot {m['overshoot']:.2f} K/step")
    ok = True
    with cf.ThreadPoolExecutor(max_workers=nj) as ex:
        for fu in cf.as_completed([ex.submit(run_one, pf, c, m, args.force) for pf, c, m in jobs]):
            m, status, wall = fu.result()
            print(f"  {m['name']:18s} {status}  ({wall:.1f} s)", flush=True)
            if args.probe and status == "ok":
                steps = m["stop"] / m["dt"]
                for n in PART2:
                    _, _, mm = part2_job(n)
                    est = wall / steps * (mm["stop"] / mm["dt"])
                    print(f"    estimate {n}: stop {mm['stop']} s, dt {mm['dt'] * 1e3:g} ms -> "
                          f"{est / 60:.1f} min alone ({360000 * steps / wall:.3g} cell-steps/s)")
            ok = ok and status.startswith(("ok", "skip"))
    print("RUNS OK" if ok else "RUNS HAD FAILURES")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
