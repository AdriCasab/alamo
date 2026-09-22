#!/usr/bin/env python3
"""Timestep-cap convergence pair for the D2d mesh case R7b_2mm (2026-09-18).

R7b_2mm ran at the harness cap DT_CAP = 4 ms (overshoot 1.42 K/step). This
script runs the same case at 16 ms and 8 ms (5.66 and 2.83 K/step), bypassing
the halving rule, through the same run_one (watchdogs, .done metadata), so
score.py can compare the three. It waits for R7_1mm.done so the 1 mm campaign
run keeps its 4 cores.

  dt_cap.py [--dts 0.016 0.008] [--no-wait]
"""
import argparse, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run as R

ap = argparse.ArgumentParser()
ap.add_argument("--dts", nargs="+", type=float, default=[0.016, 0.008])
ap.add_argument("--no-wait", action="store_true")
ap.add_argument("--force", action="store_true")
args = ap.parse_args()

if not args.no_wait:
    while not os.path.exists(os.path.join(R.OUT, "R7_1mm.done")):
        time.sleep(30)
    time.sleep(30)  # let the campaign ranks exit

for dt in args.dts:
    R.dt_rule = (lambda q, dz, cap=None, _dt=dt: _dt)  # bypass DT_CAP and the 5 K halving
    name = f"R7b_2mm_dt{int(round(dt * 1e3))}"
    pf, cmd, meta = R.job(name, **R.MATRIX["R7b_2mm"])
    assert abs(meta["dt"] - dt) < 1e-12, meta["dt"]
    meta["dt_cap_override"] = dt
    meta["parent"] = "R7b_2mm"
    print(f"{name}: dt {dt * 1e3:g} ms, overshoot {meta['overshoot']:.2f} K/step, "
          f"plot_int {[k for k in cmd if k.startswith('amr.plot_int')][0]}", flush=True)
    t0 = time.time()
    m, status, wall = R.run_one(pf, cmd, meta, args.force)
    print(f"{name}: {status}, t_end {m.get('t_end', 0):.1f} s, wall {wall / 60:.1f} min", flush=True)
