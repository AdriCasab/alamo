#!/usr/bin/env python3
"""D2h (reviewer diagnostic, 2026-09-20): three key-only 2 mm / 1 mm pairs on the
D2f support rule (Q09: pads, q 0.9, no clip), 350 s each, to separate the
inward silent front seen in D2g-1 into "local / lateral" and "gas feedback".

  run.py --cases I0_1mm RF_1mm OL_1mm I0_2mm RF_2mm OL_2mm [--jobs 3] [--force]
  run.py --kill NAME       binary path + plot_file match (never pkill -f)
  run.py --list

Switches (each = Q09's key replacements plus ONE change; no key is given twice):
  I0  idle clock off:      pinned_idle_cycles 2.0 -> 1e9 (the state rule
                           T_s = min(T_f, T_pin) in SurfaceCellFlux remains)
  RF  recirculation fixed: jet_T_ent_mode exhaust -> rec_fixed, jet_T_rec_fixed
                           1600 K (L09 T_rec at 300 s: 1624 / 1642 K)
  OL  open loop:           RF + nozzle path prescribed at the L09_2mm mean feed
                           (nozzle_z0 0.45 m, nozzle_feed 4.188e-4 m/s = 1.508 m/h,
                           0-350 s fit); the feet keys are dropped (unused under
                           prescribed) and so is free-surface dilution, which the
                           code ties to nozzle_descent = feet. OL is therefore a
                           two-mesh comparison of its own configuration, not a
                           match to L09 in absolute terms.
Loads studies/d2f_support/run.py by path (which loads d2e's; neither edited).
Stop 350 s; stall (60 s) and bottom watchdogs on; no steady-stop; dt 16 ms at
2 mm, 8 ms at 1 mm; mgs 60. Plotfiles every 50 s at BOTH meshes (1 mm: 1.58 GB
each, 8 per run) so the late 1 mm fields exist this time.
Output: studies/d2h_openloop/output/.
"""
import argparse
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
D2F = os.path.join(os.path.dirname(HERE), "d2f_support")
_spec = importlib.util.spec_from_file_location("d2f_run", os.path.join(D2F, "run.py"))
F = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(F)
E, R = F.E, F.R

OUT = os.path.join(HERE, "output")
E.OUT = OUT
R.OUT = OUT

Q09 = dict(E.SWITCH["Q09"][0])
FEED = 4.188e-4          # m/s, L09_2mm nozzle_z fit 0-350 s (1.508 m/h)
Z0 = 0.45                # m, feet on the original surface (0.40) + 0.050 stand-off
T_REC = 1600.0           # K

I0 = dict(Q09); I0["surface_patch.pinned_idle_cycles=2.0"] = ["surface_patch.pinned_idle_cycles=1.0e9"]
RF = dict(Q09); RF["surface_patch.jet_T_ent_mode=exhaust"] = ["surface_patch.jet_T_ent_mode=rec_fixed",
                                                              f"surface_patch.jet_T_rec_fixed={T_REC!r}"]
OL = dict(RF)
OL["surface_patch.nozzle_descent=feet"] = ["surface_patch.nozzle_descent=prescribed",
                                           f"surface_patch.nozzle_z0={Z0!r}",
                                           f"surface_patch.nozzle_feed={FEED!r}",
                                           # under feet the collision guard defaults to foot_r_inner
                                           # (0.028 m); under prescribed it defaults to the patch
                                           # radius (0.2 m) and aborted both OL runs at 119.4 s when
                                           # the nozzle plane reached the untouched far surface.
                                           "surface_patch.nozzle_collision_radius=0.028"]
for k in ("surface_patch.foot_r_inner=0.028", "surface_patch.foot_r_outer=0.040",
          "surface_patch.foot_standoff=0.050", "surface_patch.foot_rule=pads",
          "surface_patch.foot_npads=3", "surface_patch.foot_pad_quantile=0.9",
          "surface_patch.foot_body_clearance=1",
          "surface_patch.jet_free_surface=1", "surface_patch.jet_fs_aspect=1.0",
          "surface_patch.jet_fs_exponent=1.0"):
    OL[k] = []

E.SWITCH["I0"] = (I0, "Q09 + pinned_idle_cycles = 1e9 (idle clock off)")
E.SWITCH["RF"] = (RF, f"Q09 + jet_T_ent_mode = rec_fixed, T_rec {T_REC:g} K")
E.SWITCH["OL"] = (OL, f"RF + nozzle prescribed z0 {Z0} m, feed {FEED:.4g} m/s; feet + free-surface keys dropped")
for p in ("I0", "RF", "OL"):
    E.STOP[p] = 350.0
CASES = [f"{p}_{m}" for p in ("I0", "RF", "OL") for m in ("2mm", "1mm")]

_job = E.job


def job(name, stop=None, grid=None, plot_s=None):
    # plotfiles every 50 s at both meshes (the d2e default writes 1 mm only at 0 and stop)
    return _job(name, stop=stop, grid=grid, plot_s=50.0 if plot_s is None else plot_s)


E.job = job


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", nargs="*")
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--kill")
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()
    if a.kill:
        return E.kill(a.kill)
    if a.list:
        for n in CASES:
            pf, cmd, m = job(n)
            print(f"{n:8s} {m['switch']}; dt {m['dt'] * 1e3:g} ms {m['K_per_step']:.2f} K/step "
                  f"stop {m['stop']:g} s {[c for c in cmd if c.startswith('amr.plot_int')][0]} "
                  f"mgs {m['max_grid_size']}")
            if n.startswith("OL"):
                print("   keys:", " ".join(c for c in cmd if c.startswith("surface_patch.nozzle") or "foot" in c or "jet_fs" in c or "free_surface" in c or "idle" in c or "T_ent_mode" in c or "rec_fixed" in c))
        return
    bad = [n for n in (a.cases or []) if n not in CASES]
    if not a.cases or bad:
        sys.exit(f"give --cases from {CASES}" + (f" (unknown: {bad})" if bad else ""))
    sys.exit(0 if E.run_all(a.cases, a.jobs, a.force) else 1)


if __name__ == "__main__":
    main()
