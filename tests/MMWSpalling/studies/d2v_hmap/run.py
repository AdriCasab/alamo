#!/usr/bin/env python3
"""D2v: the body-in-the-hole h map -- separate the funnel from the wall.

  run.py --cases M0 M1 M2 --jobs 3
  run.py --cases TF TP --jobs 2      item 3: the compiled-map trace, both branches
  run.py --list
  run.py --kill NAME                 never pkill -f

M0 is built from the FROZEN d2i_gate/LI0_2mm key set, imported by path from its
harness (which chains d2h -> d2f -> d2e -> d2c; none is edited).  NOT from
input_feet_d2c: those differ in two keys (pinned_idle_cycles 2.0 vs 1.0e9,
foot_body_clearance 1 vs 0), so an input_feet_d2c baseline would miss P1's bar
for CONFIGURATION reasons and the miss would be misattributed to physics.

The h map (item 1) -- ONE expression that BRANCHES ON s
------------------------------------------------------
The stale map decays as ~1/r and never reaches a floor: with radius = 0.2 m and
the quarter domain's corner at r = 168 mm, EVERY column is in the patch and
every column with s > 0 is excavated at 160-360 W/m^2K.  That is the funnel.

    floor (s > 0) :  h_imp(r,s)*w1(r) + h_ann*(1 - w1(r))
    wall  (s <= 0):  h_ann                       <- FLAT, deliberately
    both      x   :  branch*w2(r) + h_out*(1 - w2(r))

w1: 1 to r_imp = 40 mm, raised cosine to 0 by 56 mm.   r_imp IS the skirt
    radius foot_r_outer: starting the blend inside the 28-40 mm pad annulus
    costs the outer pad rock 10-20 % of its h, slows its recession and slows
    the feed -- D2o-0's jam mechanism.
w2: 1 to r_mouth = 60 mm, raised cosine to 0 by r_out = 76 mm.
Both blends are 16 mm = 8 coarse cells at dz = 2 mm: D2o-0b failed on a step,
D2t passed on 8 cells.  Taper width is a NUMERICAL parameter, never tuned.

The wall branch is FLAT h_ann and must stay flat.  h_imp clamps s to
[0.015, 0.09], so above the plane it returns its STRONGEST near-field profile;
blended, a column at r = 42-47 mm would receive 245-405 kW/m^2 against a floor
q_pin of ~520 kW/m^2 -- ~1 m/h of wall recession that would drive the hole out
to the hand-over radius regardless of physics.

jet_T_above -- declared, with the assumption stated, never fitted
----------------------------------------------------------------
The frozen run's own uncooled exhaust: jet_T_exhaust = 1689.5 K at 572 s
(1676.5-1689.7 over t >= 500 s).  Meier's water jacket is NOT in the model.
Its "186 W/K" is the WATER stream's mcp (160 L/h x 4186), and Meier does not
report the water dT, so the gas drop cannot be computed: it is 19 K for
dT_water = 5 K and 194 K for 50 K.  Declaring a drop would mean inventing
dT_water.  So jet_T_above = 1690 K, the measured uncooled value, which is
  * an UPPER bound on the physical gas temperature (no jacket), and
  * a LOWER bound relative to M1/M2's own exhaust, which will run hotter
    because the confined map extracts less power from the block top.
Both directions are stated; nothing is adjusted to a result.
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


D2I = load("d2i_run", os.path.join(STUD, "d2i_gate", "run.py"))
E = D2I.E                      # d2e: kill, rows_upto
R = D2I.R                      # d2c: run_one, and the stall/bottom watchdogs the
                               # frozen LI0_2mm run itself was driven by -- M0
                               # must keep them or it is not that run
OUT = os.path.join(HERE, "output")
for _m in (E, R, D2I):
    _m.OUT = OUT

FROZEN_LI0 = os.path.join(STUD, "d2i_gate", "output", "LI0_2mm")

STOP = 500.0
H_ANN = 142.0                  # declared; band 81 DB / 142 area-exp / 250 entry
H_OUT = 1.0e-3                 # outside the mouth (sealed wellhead)
T_ABOVE = 1690.0               # measured uncooled exhaust; see the docstring
R_IMP, R_ANN_R = 0.040, 0.056  # w1: skirt radius -> annulus
R_MOUTH, R_OUT = 0.060, 0.076  # w2: mouth -> outside
F_BELOW = 0.2                  # frozen stand-in, never scanned
TRACE_STOP, TRACE_STEPS = 0.48, 30

FEET_ONLY = ("surface_patch.foot_r_inner", "surface_patch.foot_r_outer",
             "surface_patch.foot_standoff", "surface_patch.foot_rule",
             "surface_patch.foot_npads", "surface_patch.foot_pad_quantile",
             "surface_patch.foot_body_clearance",
             "surface_patch.nozzle_collision_radius", "surface_patch.nozzle_feed_max",
             # jet_free_surface does NOT abort under prescribed descent -- it
             # silently misbehaves, which is worse: its test is
             # r > foot_r_outer (:2223) and foot_r_outer defaults to 0.0
             # (:5962) and is never parsed outside the feet branch, so EVERY
             # column would be treated as free surface.  Strip it explicitly.
             "surface_patch.jet_free_surface", "surface_patch.jet_fs_aspect",
             "surface_patch.jet_fs_exponent")


def cosramp(a, b, var="r"):
    """Raised cosine 1 -> 0 across [a, b]. Numeric first in every min/max (the
    AMReX parser evaluates a/max(x,num) as 0 and the frozen h_expr itself uses
    min(0.09,max(0.015,s))); the denominator is an inlined literal, never a
    max()."""
    return (f"0.5*(1.0+cos(3.141592653589793*"
            f"(min({b!r},max({a!r},{var}))-{a!r})/{b - a!r}))")


def h_expr_of(h_imp):
    w1, w2 = cosramp(R_IMP, R_ANN_R), cosramp(R_MOUTH, R_OUT)
    floor = f"(({h_imp})*({w1})+{H_ANN!r}*(1.0-({w1})))"
    wall = f"{H_ANN!r}"
    return f"(if(s>0.0,{floor},{wall})*({w2})+{H_OUT!r}*(1.0-({w2})))"


def w1_of(r):
    rc = np.clip(r, R_IMP, R_ANN_R)
    return 0.5 * (1.0 + np.cos(np.pi * (rc - R_IMP) / (R_ANN_R - R_IMP)))


def w2_of(r):
    rc = np.clip(r, R_MOUTH, R_OUT)
    return 0.5 * (1.0 + np.cos(np.pi * (rc - R_MOUTH) / (R_OUT - R_MOUTH)))


def h_imp_of(r, s):
    """The frozen input_feet impingement term, in numpy. A FILTER ONLY: the
    D2o-0 bug lived inside the AMReX parser, so this cannot prove the compiled
    map. That is what hmap_trace.py is for."""
    rr = np.maximum(0.01875, r)
    sc = np.minimum(0.09, np.maximum(0.015, s))
    b = 0.1 * (sc / 0.0075 - 6.0)
    num = (0.0075 / (2.0 * rr)) * (rr * rr + 2.0 * b * 0.0075 * rr
                                   - 1.1 * b * 0.0075 * 0.0075)
    den = (rr + b * 0.0075) ** 2
    return 1449.0 * (num / den) / 0.19763872491145218 * (0.09 / np.maximum(0.09, s)) ** 1.0


def h_map(r, s, above=False):
    """The full intended map, in numpy. Same caveat as h_imp_of."""
    w1, w2 = w1_of(r), w2_of(r)
    branch = np.where(above, H_ANN, h_imp_of(r, s) * w1 + H_ANN * (1.0 - w1))
    return branch * w2 + H_OUT * (1.0 - w2)


def frozen_keys():
    """The frozen LI0_2mm argv, verbatim."""
    _, cmd, meta = D2I.job("LI0_2mm")
    return list(cmd), dict(meta)


def _set(cmd, key, val):
    pre = key + "="
    hit = [n for n, c in enumerate(cmd) if c.startswith(pre)]
    if len(hit) > 1:
        raise SystemExit(f"{key} appears {len(hit)} times in the base keys")
    if hit:
        cmd[hit[0]] = pre + val
    else:
        cmd.append(pre + val)


def _drop(cmd, key):
    return [c for c in cmd if not c.startswith(key + "=")]


LEGS = ("M0", "M1", "M2", "TF", "TP")


def job(name):
    cmd, meta = frozen_keys()
    h_old = [c for c in cmd if c.startswith("surface_patch.h_expr=")][0]
    h_imp = h_old.split("=", 1)[1].strip('"')
    pf = os.path.join(OUT, name)
    _set(cmd, "plot_file", pf)
    _set(cmd, "stop_time", repr(STOP))
    new_map = h_expr_of(h_imp)

    if name != "M0":
        _set(cmd, "surface_patch.h_expr", f'"{new_map}"')
        _set(cmd, "surface_patch.side_face_flux", "1")
        _set(cmd, "surface_patch.side_face_factor", repr(F_BELOW))
        # MANDATORY on M1, not cosmetic: unset, the sentinel falls back to
        # side_face_factor = 0.2, and with h = 0 against 293 K the Newton solve
        # returns q_in < 0 -- M1 would COOL its above-plane faces.
        _set(cmd, "surface_patch.side_face_factor_above",
             "1.0" if name in ("M2", "TP") else "0.0")
    if name in ("M2", "TP"):
        _set(cmd, "surface_patch.jet_above_h", "1")
        _set(cmd, "surface_patch.jet_T_above", repr(T_ABOVE))

    if name in ("TF", "TP"):
        _set(cmd, "stop_time", repr(TRACE_STOP))
        _set(cmd, "amr.plot_int", str(TRACE_STEPS))
        _set(cmd, "surface_patch.jet_profile_interval", "0.016")
    if name == "TP":
        # The wall branch cannot be reached by a feet run: the nozzle rides
        # foot_standoff = 50 mm above the rock, so at t = 0 every column has
        # s = +50 mm.  Put the plane BELOW the original surface instead.
        for k in FEET_ONLY:
            cmd = _drop(cmd, k)
        _set(cmd, "surface_patch.nozzle_descent", "prescribed")
        _set(cmd, "surface_patch.nozzle_z0", "0.39")
        _set(cmd, "surface_patch.nozzle_feed", "0.0")
        _set(cmd, "surface_patch.nozzle_collision_radius", "0.0")

    meta.update(name=name, stop=float([c for c in cmd if c.startswith("stop_time=")][0]
                                      .split("=")[1]),
                h_imp=h_imp, h_expr=new_map if name != "M0" else h_imp,
                h_ann=H_ANN, h_out=H_OUT, r_imp=R_IMP, r_ann_r=R_ANN_R,
                r_mouth=R_MOUTH, r_out=R_OUT,
                side_face_flux=0 if name == "M0" else 1,
                side_face_factor=None if name == "M0" else F_BELOW,
                side_face_factor_above=(None if name == "M0"
                                        else (1.0 if name in ("M2", "TP") else 0.0)),
                jet_above_h=1 if name in ("M2", "TP") else 0,
                jet_T_above=T_ABOVE if name in ("M2", "TP") else None,
                trace=name in ("TF", "TP"), frozen_ref=FROZEN_LI0 if name == "M0" else None)
    return pf, cmd, meta


def keep_awake():
    """Hold a no-sleep assertion for as long as this harness lives. These runs
    are not checkpointed, so a host sleep loses one outright. Asserted here so
    it cannot be forgotten at the call site; harmless if already wrapped."""
    import subprocess
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
        print(f"  {m['name']:3s} stop {m['stop']:g} s, dt {m['dt'] * 1e3:g} ms, "
              f"side_flux {m['side_face_flux']}, f_above {m['side_face_factor_above']}, "
              f"jet_above_h {m['jet_above_h']}, T_above {m['jet_T_above']}", flush=True)
    ok = True
    with cf.ThreadPoolExecutor(max_workers=min(jobs, len(js))) as ex:
        for fu in cf.as_completed([ex.submit(R.run_one, pf, c, m, force) for pf, c, m in js]):
            m, status, wall = fu.result()
            print(f"  {m['name']:3s} {status} at t = {m.get('t_end', 0):g} s "
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
        return E.kill(a.kill)
    if a.list:
        print(f"h_ann {H_ANN:g}  h_out {H_OUT:g}  T_above {T_ABOVE:g} K")
        print(f"w1 {R_IMP*1e3:g}->{R_ANN_R*1e3:g} mm ({(R_ANN_R-R_IMP)/2e-3:g} coarse cells), "
              f"w2 {R_MOUTH*1e3:g}->{R_OUT*1e3:g} mm ({(R_OUT-R_MOUTH)/2e-3:g} cells)")
        for n in LEGS:
            _, cmd, m = job(n)
            print(f"\n{n}: stop {m['stop']:g} s  "
                  f"{' '.join(c for c in cmd if c.startswith(('surface_patch.side_face','surface_patch.jet_above','surface_patch.jet_T_above','surface_patch.nozzle_descent')))}")
        print("\nh(r) intended [W/m2K]:")
        rs = np.array([5, 20, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 90, 120]) * 1e-3
        print("  r [mm]     " + " ".join(f"{x*1e3:7.0f}" for x in rs))
        print("  floor s=50 " + " ".join(f"{h_map(r, 0.05):7.2f}" for r in rs))
        print("  wall       " + " ".join(f"{h_map(r, -0.01, True):7.2f}" for r in rs))
        return
    bad = [n for n in (a.cases or []) if n not in LEGS]
    if not a.cases or bad:
        sys.exit(f"give --cases from {list(LEGS)}" + (f" (unknown: {bad})" if bad else ""))
    sys.exit(0 if run_all(a.cases, a.jobs, a.force) else 1)


if __name__ == "__main__":
    main()
