#!/usr/bin/env python3
"""D2o-0: prescribed floor/wall heating — does a cylindrical hole appear?

  run.py --list                       PRE-FLIGHT: zones, h and T_gas per zone,
                                      column counts (ACTIVE_STEP item 2)
  run.py --cases W_2mm W_1mm --jobs 2         primary, ~2.5 h (1 mm leg)
  run.py --cases Wz_2mm Wsens_2mm --jobs 2    variants, ~25 min
  run.py --kill NAME                  binary path + plot_file match

Binary: bin/mmwspalling-3d-g++ = 62dea451... Key-only; no source change, no build.

The field
---------
The Meier `LI0` geometry (feet, descent, 0.40 m box) with `jet_closure = none`
and the heating written by hand in two zones, switched at the floor/wall corner
INSIDE the band (not at the nozzle plane):

  FLOOR, s > 45 mm   the control's own J-M h(r, s) string, verbatim, and the
                     control's own radial gas profile
  WALL,  s <= 45 mm  h = 40 W/m^2K (the annulus duct estimate) at the control's
                     own exhaust temperature
  Wz variant         the skirt shadow, r > 40 mm inside the band, at ~zero

Measured on the control (d2i_gate/output/LI0_2mm) the switch lands at r ~ 40 mm
-- the skirt radius -- with the floor zone spanning r 1..41 mm and the wall zone
r 40..69 mm. `--list` prints the counts.

Three things the packet does not state, resolved from the named control run and
declared in PREDICTIONS.md rather than invented
------------------------------------------------------------------------------
1. **The wall gas temperature.** ACTIVE_STEP gives the wall coefficient but not
   its gas temperature. Taken as the control's own `jet_T_mouth`, mean over
   150-450 s = **1610.1 K**.
2. **The floor gas temperature.** With `jet_closure = none` there is no march, so
   T_gas must be prescribed for the floor too, and "unchanged" has to mean the
   control's own field. A quadratic fit to the control's `_jet_profile.csv`
   T_gas(r), averaged over 150-450 s, reproduces it to **3.4 K** over the floor
   zone: T = 1742.48 - 435.739 r - 42654 r^2  (r in m).
3. **The s <= 0 rule must be written by hand here.** Under `jet_closure =
   enthalpy` the code zeroes flux above the nozzle plane itself; under
   `jet_closure = none` it does NOT, and it strictly requires h > 0 and T_gas > 0
   at every in-patch column (MMWSpalling.H:1690). Without the guard the whole
   undisturbed surface out to the 0.2 m patch radius would be heated. Emulated
   as h = 1e-3 above the plane, the d2j0b_plane convention (the parser forbids
   h <= 0).
"""
import argparse
import concurrent.futures as cf
import importlib.util
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)


def load(name, path):
    sp = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


B = load("d2j0b_run", os.path.join(STUD, "d2j0b_plane", "run.py"))   # J: the runner
J = B.J
I = load("d2i_run", os.path.join(STUD, "d2i_gate", "run.py"))        # the LI0 keys
AN = load("d2b_analyze", os.path.join(STUD, "d2b_feet_rop", "analyze.py"))

OUT = os.path.join(HERE, "output")
B.OUT = OUT
J.OUT = OUT
I.OUT = OUT
I.E.OUT = OUT
I.R.OUT = OUT

CONTROL = os.path.join(STUD, "d2i_gate", "output", "LI0_2mm")        # read-only reference
CONTROL_1MM = os.path.join(STUD, "d2i_gate", "output", "LI0_1mm")

STOP = 450.0
S_SWITCH = 0.045           # floor/wall corner [m]
R_SKIRT = 0.040            # skirt outer radius [m]
H_OFF = 1.0e-3             # emulated zero (parser forbids h <= 0)
T_WALL = 1610.1            # control jet_T_mouth, mean 150-450 s [K]
FLOOR_T = (1742.48, -435.739, -42654.0)     # a + b r + c r^2, fit to the control

# The annulus correlation. D_h = 2*gap = 0.020 m, u ~ 23 m/s, exhaust at ~1600 K
# (nu ~ 2.4e-4 m^2/s, k ~ 0.10 W/mK, Pr ~ 0.72) -> Re = u D_h / nu ~ 1.9e3.
#   laminar fully developed, annulus one side heated   Nu ~ 4.0  -> h ~ 20
#   laminar thermally developing (Hausen, D/L = 0.04)  Nu ~ 6.1  -> h ~ 30
#   turbulent Dittus-Boelter 0.023 Re^0.8 Pr^0.4       Nu ~ 9.2  -> h ~ 46
# ACTIVE_STEP states 40 W/m^2K as the duct estimate; that is the PRIMARY and is
# used as given, not re-derived. The pre-registered sensitivity is the other
# branch, laminar fully developed at 20 -- a factor 2, the spread the packet
# names. f is never tuned to make a prediction pass.
H_WALL = 40.0
H_WALL_SENS = 20.0

CASES = ["W_2mm", "W_1mm", "Wz_2mm", "Wsens_2mm"]
DROP = ("jet_closure", "jet_stagnation", "jet_entrained_mass", "jet_T_ent_mode",
        "jet_T_ent", "jet_mdot", "jet_cp", "jet_D", "jet_T_nozzle", "jet_core_length",
        "jet_bin_update", "jet_negative_flux", "jet_profile_interval", "jet_dr",
        "jet_decay_diameter", "jet_De_ref", "jet_free_surface", "jet_fs_aspect",
        "jet_fs_exponent", "jet_T_rec_fixed")


def floor_h_string(cmd):
    """The control's h_expr, verbatim — the floor field must be unchanged."""
    hit = [c for c in cmd if c.startswith("surface_patch.h_expr=")]
    if len(hit) != 1:
        sys.exit(f"expected exactly one h_expr in the control keys, got {len(hit)}")
    return hit[0].split("=", 1)[1].strip('"')


def floor_T_string():
    a, b, c = FLOOR_T
    # clamped so the unused branch can never go <= 0 (the parser evaluates it)
    return f"max(1400.0,{a!r}+{b!r}*r+{c!r}*r*r)"


def exprs(case, cmd):
    """(h_expr, T_flame_expr). Numeric arguments first in every min/max (the
    AMReX 25.12 rewrite trap); no a/max(x,num) (that form evaluates to 0)."""
    fh = floor_h_string(cmd)
    hw = H_WALL_SENS if case == "Wsens" else H_WALL
    if case == "Wz":
        # skirt shadow: r > R_SKIRT inside the band is reached only through the
        # slots -> ~zero. Inside R_SKIRT the band still grazes.
        wall = f"if(r>{R_SKIRT!r},{H_OFF!r},{hw!r})"
    else:
        wall = f"{hw!r}"
    band = f"if(s>{S_SWITCH!r},{fh},{wall})"
    h = f"if(s>0.0,{band},{H_OFF!r})"          # the s <= 0 rule, by hand
    T = f"if(s>{S_SWITCH!r},{floor_T_string()},{T_WALL!r})"
    return h, T


def h_py(case, r, s):
    """numpy twin of the h expression, for the pre-flight."""
    wj = AN.wj
    hw = H_WALL_SENS if case == "Wsens" else H_WALL
    floor = wj.h_py(1449.0, r, s, "power", 1.0)
    wall = np.full_like(np.asarray(r, dtype=float), hw)
    if case == "Wz":
        wall = np.where(np.asarray(r) > R_SKIRT, H_OFF, wall)
    out = np.where(np.asarray(s) > S_SWITCH, floor, wall)
    return np.where(np.asarray(s) > 0.0, out, H_OFF)


def T_py(r, s):
    a, b, c = FLOOR_T
    rr = np.asarray(r, dtype=float)
    floor = np.maximum(1400.0, a + b * rr + c * rr * rr)
    return np.where(np.asarray(s) > S_SWITCH, floor, T_WALL)


def job(name):
    case, mesh = name.rsplit("_", 1)
    _, cmd, meta = I.job(f"LI0_{mesh}")
    h, T = exprs(case, cmd)
    pf = os.path.join(OUT, name)
    out = []
    for c in cmd:
        key = c.split("=", 1)[0]
        if key.startswith("surface_patch.jet_") and key.split(".", 1)[1] in DROP:
            continue                                   # jet_closure = none (default)
        if c.startswith("plot_file="):
            out.append(f"plot_file={pf}")
        elif c.startswith("stop_time="):
            out.append(f"stop_time={STOP!r}")
        elif c.startswith("surface_patch.h_expr="):
            out.append(f'surface_patch.h_expr="{h}"')
        else:
            out.append(c)
    out.append(f'surface_patch.T_flame_expr="{T}"')
    left = [c for c in out if c.split("=", 1)[0].startswith("surface_patch.jet_")]
    if left:
        sys.exit(f"{name}: jet keys survive into a jet_closure = none run: {left}")
    if not any(c.startswith("surface_patch.h_expr=") for c in out):
        sys.exit(f"{name}: h_expr missing")
    meta = dict(meta, name=name, case=case, stop=STOP, h_wall=H_WALL_SENS if case == "Wsens"
                else H_WALL, s_switch=S_SWITCH, r_skirt=R_SKIRT, T_wall=T_WALL,
                skirt="zero" if case == "Wz" else "grazing", h_expr=h, T_expr=T,
                jet_closure="none")
    return pf, out, meta


def preflight():
    """ACTIVE_STEP item 2: zone boundaries, min/max h and T_gas per zone, and the
    column count in each. A silent zero or a half-height jump would make the run
    meaningless, so this is printed before anything launches."""
    print(f"binary {J.binary_sha()[:8]}…   stop {STOP:g} s   jet_closure = none")
    print(f"zones: FLOOR s > {S_SWITCH*1e3:g} mm (control h_expr + fitted T_gas(r));  "
          f"WALL s <= {S_SWITCH*1e3:g} mm (h = {H_WALL:g}, T = {T_WALL:g} K);  "
          f"above the plane s <= 0 -> h = {H_OFF:g} (emulated zero)")
    print(f"skirt radius {R_SKIRT*1e3:g} mm; sensitivity h_wall = {H_WALL_SENS:g}\n")
    run = AN.Run(CONTROL, AN.meta_of(CONTROL))
    for name in CASES:
        case, mesh = name.rsplit("_", 1)
        _, cmd, m = job(name)
        print(f"{name:10s} skirt {m['skirt']:8s} h_wall {m['h_wall']:5g}  dz {m['dz']*1e3:g} mm "
              f"dt {m['dt']*1e3:g} ms stop {m['stop']:g} s")
        if mesh != "2mm":
            continue
        for T in (150.0, 250.0, 350.0):
            j = int(np.argmin(np.abs(run.th["time"] - T)))
            zn = run.th["nozzle_z"][j]
            s = zn - (run.lz - run.depth(T))
            for lab, sel in (("above plane", s <= 0.0),
                             ("WALL 0<s<=45", (s > 0.0) & (s <= S_SWITCH)),
                             ("FLOOR s>45", s > S_SWITCH)):
                n = int(sel.sum())
                if n == 0:
                    print(f"    t={T:.0f}s {lab:13s}: 0 columns  ** EMPTY ZONE **")
                    continue
                hv, Tv = h_py(case, run.r[sel], s[sel]), T_py(run.r[sel], s[sel])
                print(f"    t={T:.0f}s {lab:13s}: {n:5d} cols  r {run.r[sel].min()*1e3:4.0f}"
                      f"–{run.r[sel].max()*1e3:3.0f} mm  h {hv.min():9.3g}–{hv.max():<9.3g} "
                      f"T_gas {Tv.min():6.1f}–{Tv.max():6.1f} K")
        print()
    print("control geometry is d2i_gate/output/LI0_2mm (read-only, not re-run).")


def run_all(names, jobs, force):
    J.keep_awake()
    js = [job(n) for n in names]
    print(f"binary {J.binary_sha()[:8]}…  {len(js)} runs, {min(jobs, len(js))} in parallel "
          f"-> {OUT}", flush=True)
    for _, _, m in js:
        print(f"  {m['name']:10s} skirt {m['skirt']}, h_wall {m['h_wall']:g}; "
              f"dz {m['dz']*1e3:g} mm, dt {m['dt']*1e3:g} ms, stop {m['stop']:g} s", flush=True)
    ok = True
    with cf.ThreadPoolExecutor(max_workers=min(jobs, len(js))) as ex:
        for fu in cf.as_completed([ex.submit(J.run_one, pf, c, m, force) for pf, c, m in js]):
            m, status, wall = fu.result()
            print(f"  {m['name']:10s} {status} at t = {m.get('t_end', 0):g} s "
                  f"({wall/60:.1f} min)", flush=True)
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
        return preflight()
    bad = [n for n in (a.cases or []) if n not in CASES]
    if not a.cases or bad:
        sys.exit(f"give --cases from {CASES}" + (f" (unknown: {bad})" if bad else ""))
    sys.exit(0 if run_all(a.cases, a.jobs, a.force) else 1)


if __name__ == "__main__":
    main()
