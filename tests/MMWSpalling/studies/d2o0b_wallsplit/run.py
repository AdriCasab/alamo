#!/usr/bin/env python3
"""D2o-0b: prescribed floor/wall heating with the field bounded in radius.

  run.py --list                        G1 field-integrity gate (MUST pass first)
  run.py --cases W_2mm W_1mm --jobs 2  primary
  run.py --cases Wz_2mm Wsens_2mm --jobs 2
  run.py --kill NAME

Binary: bin/mmwspalling-3d-g++ = 62dea451... Key-only; no source change, no build.

Copied from studies/d2o0_wallsplit/run.py (frozen, not edited) with its three
defects fixed.

1. THE CORNER IS RADIAL, at the skirt
-------------------------------------
  FLOOR  r <= 40 mm, s > 0 : the control's own h(r, s) string, verbatim
  WALL   r >  40 mm, s > 0 : h = 40 W/m^2K (the annulus duct estimate)
  s <= 0                   : emulated zero (1e-3)
  r > 170 mm               : emulated zero (the hard outer bound)
Both zones share one T_gas(r).

D2o-0 put the corner at a STAND-OFF of 45 mm and it jammed the burner: the pads
rest on the 0.9 nearest-rank quantile of annulus heights, so ~10 % of annulus
rock sits above the pad plane, the 45 mm switch cut through that carrying
population, and those columns froze at 444-482 K and held the tool up. The
carriers live at r in [foot_r_inner, foot_r_outer] = [28, 40] mm, so **a radial
corner at 40 mm puts every carrier in the floor zone by construction.** A
stand-off corner provably cannot.

2. THE GAS FIELD IS BOUNDED IN RADIUS AND DECAYS TO AMBIENT
-----------------------------------------------------------
D2o-0 fitted T_gas on r <= 68 mm and guarded the parser with max(1400, ...),
which became a 1400 K blanket over everything beyond r = 85 mm. Here the fit is
made on the control's FULL tabulated range (0-168 mm) and decays to ambient:

  T(r) = 293.15 + 346.5465*exp(-(r/0.0671158)^5.89) + 1173.4772/(1+(r/0.1009858)^2)

max error 24 K against the control's own table (criterion 60 K), monotone over
0-250 mm, and T -> 293.15 K as r grows. **The parser guard is AMBIENT, not a
working gas temperature**, so it can never set a physical value -- that is the
exact defect being fixed.

3. RUN LONG ENOUGH TO SCORE A SHAFT
-----------------------------------
stop 600 s (D2o-0 used 450 s and never reached 10 cm because it jammed).
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


B = load("d2j0b_run", os.path.join(STUD, "d2j0b_plane", "run.py"))
J = B.J
I = load("d2i_run", os.path.join(STUD, "d2i_gate", "run.py"))
AN = load("d2b_analyze", os.path.join(STUD, "d2b_feet_rop", "analyze.py"))

OUT = os.path.join(HERE, "output")
B.OUT = OUT
J.OUT = OUT
I.OUT = OUT
I.E.OUT = OUT
I.R.OUT = OUT

CTRL_DIR = os.path.join(STUD, "d2i_gate", "output")
CONTROL = os.path.join(CTRL_DIR, "LI0_2mm")

STOP = 600.0
R_CORNER = 0.040           # the skirt = foot_r_outer; carriers live inside it
R_OUTER = 0.170            # hard outer bound on the prescribed field
H_OFF = 1.0e-3             # emulated zero (the parser forbids h <= 0)
T_AMB = 293.15             # the ONLY admissible guard value
TG = (346.5465, 0.0671158, 5.89, 1173.4772, 0.1009858)   # A, b, c, B, d
H_WALL, H_WALL_SENS = 40.0, 20.0

CASES = ["W_2mm", "W_1mm", "Wz_2mm", "Wsens_2mm"]
DROP = ("jet_closure", "jet_stagnation", "jet_entrained_mass", "jet_T_ent_mode",
        "jet_T_ent", "jet_mdot", "jet_cp", "jet_D", "jet_T_nozzle", "jet_core_length",
        "jet_bin_update", "jet_negative_flux", "jet_profile_interval", "jet_dr",
        "jet_decay_diameter", "jet_De_ref", "jet_free_surface", "jet_fs_aspect",
        "jet_fs_exponent", "jet_T_rec_fixed")
# G1 reference: the control's own march (studies/d2i_gate/output/LI0_2mm_jet_profile.csv),
# reproduced with this study's estimator before use -- 2710.7 W shallow / 1390.6 W deep,
# radial shares 41.9/29.9/22.9/5.3 % and 64.1/35.9/0/0 %.
G1_EDGES = (0.0, 0.040, 0.070, 0.120, 0.200)
G1_TIMES = (50.0, 150.0, 300.0)      # t = 25 s has no control plotfile (they are every 50 s)


def floor_h_string(cmd):
    hit = [c for c in cmd if c.startswith("surface_patch.h_expr=")]
    if len(hit) != 1:
        sys.exit(f"expected one h_expr in the control keys, got {len(hit)}")
    return hit[0].split("=", 1)[1].strip('"')


def T_string():
    A, b, c, Bc, d = TG
    core = (f"{T_AMB!r}+{A!r}*exp(-pow(max(1.0e-9,r)/{b!r},{c!r}))"
            f"+{Bc!r}/(1.0+(r/{d!r})*(r/{d!r}))")
    return f"max({T_AMB!r},{core})"        # guard is AMBIENT, never a working gas value


def exprs(case, cmd):
    """(h_expr, T_flame_expr). Numeric first in every min/max; no a/max(x,num)."""
    fh = floor_h_string(cmd)
    hw = H_WALL_SENS if case == "Wsens" else H_WALL
    wall = f"{H_OFF!r}" if case == "Wz" else f"{hw!r}"
    zoned = f"if(r>{R_CORNER!r},{wall},{fh})"                   # RADIAL corner
    bounded = f"if(r>{R_OUTER!r},{H_OFF!r},{zoned})"             # hard outer bound
    h = f"if(s>0.0,{bounded},{H_OFF!r})"                        # s <= 0 off
    return h, T_string()


def h_py(case, r, s):
    wj = AN.wj
    r = np.asarray(r, dtype=float)
    s = np.asarray(s, dtype=float)
    hw = H_WALL_SENS if case == "Wsens" else H_WALL
    wall = H_OFF if case == "Wz" else hw
    out = np.where(r > R_CORNER, wall, wj.h_py(1449.0, r, s, "power", 1.0))
    out = np.where(r > R_OUTER, H_OFF, out)
    return np.where(s > 0.0, out, H_OFF)


def T_py(r):
    A, b, c, Bc, d = TG
    r = np.maximum(1.0e-9, np.asarray(r, dtype=float))
    return np.maximum(T_AMB, T_AMB + A * np.exp(-(r / b) ** c) + Bc / (1.0 + (r / d) ** 2))


def job(name):
    case, mesh = name.rsplit("_", 1)
    _, cmd, meta = I.job(f"LI0_{mesh}")
    h, T = exprs(case, cmd)
    pf = os.path.join(OUT, name)
    out = []
    for c in cmd:
        key = c.split("=", 1)[0]
        if key.startswith("surface_patch.jet_") and key.split(".", 1)[1] in DROP:
            continue
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
    meta = dict(meta, name=name, case=case, stop=STOP, r_corner=R_CORNER, r_outer=R_OUTER,
                h_wall=H_WALL_SENS if case == "Wsens" else (0.0 if case == "Wz" else H_WALL),
                h_expr=h, T_expr=T, jet_closure="none")
    return pf, out, meta


# ------------------------------------------------------------------ G1
def control_march(t):
    """The control's own face power and radial split at time t, from its jet
    profile log -- the reference this field must not exceed."""
    a = np.genfromtxt(os.path.join(CTRL_DIR, "LI0_2mm_jet_profile.csv"),
                      delimiter=",", names=True)
    ts = np.unique(a["time"])
    tt = ts[int(np.argmin(np.abs(ts - t)))]
    m = a["time"] == tt
    tot = float(a["P_share"][m].sum())
    split = []
    for i in range(4):
        s = m & (a["r_lo"] >= G1_EDGES[i]) & (a["r_lo"] < G1_EDGES[i + 1])
        split.append(float(a["P_share"][s].sum()))
    return float(tt), tot, split


def surface_state(pf, want_t):
    """(t, T_surface, z_face) per column from the control's plotfile nearest to
    but not after want_t. Column order matches AN.Run.r (index j*nx + i)."""
    import yt
    names = sorted(d for d in os.listdir(pf) if d.endswith("cell") and d[:-4].isdigit())
    best = None
    for n in names:
        ds = yt.load(os.path.join(pf, n))
        if float(ds.current_time) <= want_t + 1e-6:
            best = (float(ds.current_time), ds)
    if best is None:
        return None
    t, ds = best
    nx, ny, nz = (int(v) for v in ds.domain_dimensions)
    dz = float(ds.domain_right_edge[2]) / nz
    g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    T = np.asarray(g["Temp"])
    solid = np.asarray(g["removed"]) < 0.5
    any_s = solid.any(axis=2)
    ktop = np.where(any_s, nz - 1 - np.argmax(solid[:, :, ::-1], axis=2), -1)
    ii, jj = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
    Ts = np.where(any_s, T[ii, jj, np.clip(ktop, 0, nz - 1)], np.nan)
    zf = np.where(any_s, (ktop + 1) * dz, np.nan)
    return t, Ts.T.ravel(), zf.T.ravel()


def g1(verbose=True):
    """The field-integrity gate. Evaluates the prescribed field on the CONTROL's
    own geometry and compares total face power and radial split against the
    control's own march at the same instant."""
    run = AN.Run(CONTROL, AN.meta_of(CONTROL))
    th = run.th
    ok = True
    L = []
    for want in G1_TIMES:
        st = surface_state(CONTROL, want)
        if st is None:
            continue
        t, Ts, zf = st
        j = int(np.argmin(np.abs(th["time"] - t)))
        zn = th["nozzle_z"][j]
        s = zn - zf
        good = np.isfinite(s) & np.isfinite(Ts)
        ss = np.where(good, s, -1.0)
        Tg = T_py(run.r)
        area = run.dz * run.dz

        def power(hv):
            q = np.where(good, np.maximum(0.0, hv * (Tg - Ts)), 0.0)
            tot = 4.0 * float(q.sum()) * area
            sp = [4.0 * float(q[(run.r >= G1_EDGES[i]) & (run.r < G1_EDGES[i + 1])].sum()) * area
                  for i in range(4)]
            return tot, sp

        # The comparison must use ONE estimator on BOTH fields. Reading the top
        # solid cell's temperature from a plotfile is not the pinned surface
        # temperature the solver uses (robin_form = pinned sets
        # T_s = min(T_face, T_pin)), so absolute powers from this estimator are
        # biased high. Evaluating the CONTROL's own field the same way makes the
        # bias cancel in the ratio, and the calibration is printed so the size of
        # it is visible.
        h_ctl = np.where(ss > 0.0, AN.wj.h_py(1449.0, run.r, ss, "power", 1.0), H_OFF)
        tot, split = power(h_py("W", run.r, ss))
        ctot_est, csplit_est = power(h_ctl)
        ct, ctot_log, csplit_log = control_march(t)
        far = (split[2] + split[3]) / tot * 100 if tot else 0.0
        cfar = (csplit_est[2] + csplit_est[3]) / ctot_est * 100 if ctot_est else 0.0
        far_ok = far <= cfar + 1.0            # 1 percentage point; the emulated
        ok = ok and far_ok                    # zero (1e-3) contributes a trace
        L.append((t, tot, ctot_est, ctot_log, split, csplit_est, far, cfar, far_ok))
    if verbose:
        print("\nG1 — field integrity (prescribed field on the CONTROL's geometry)\n")
        print("One estimator on both fields (see the note in g1()); the control's logged")
        print("march total is shown only to calibrate the estimator's bias.\n")
        print("| t [s] | prescribed | control field, same estimator | ratio | "
              "control logged | estimator bias | r>70mm share (mine / control) | ok |")
        print("|---|---|---|---|---|---|---|---|")
        for t, tot, cest, clog, sp, cs, far, cfar, fok in L:
            print(f"| {t:.0f} | {tot:.0f} W | {cest:.0f} W | **{tot/cest:.2f}×** | "
                  f"{clog:.0f} W | {cest/clog:.1f}× | {far:.1f} % / {cfar:.1f} % | "
                  f"{'ok' if fok else '**EXCEEDS**'} |")
        print("\n| t [s] | zone | prescribed [W] | control field, same estimator [W] |")
        print("|---|---|---|---|")
        for t, tot, cest, clog, sp, cs, *_ in L:
            for i in range(4):
                print(f"| {t:.0f} | r {G1_EDGES[i]*1e3:.0f}–{G1_EDGES[i+1]*1e3:.0f} mm | "
                      f"{sp[i]:.0f} | {cs[i]:.0f} |")
        print("\nT_gas(r) fit against the control's own table:")
        rt = np.array([0, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112,
                       120, 128, 136, 144, 152, 160, 168]) / 1000.0
        Tt = np.array([1824, 1813, 1781, 1736, 1685, 1638, 1576, 1440, 1273, 1149, 1051,
                       973, 908, 853, 806, 766, 733, 705, 681, 660, 641, 624], float)
        err = np.abs(T_py(rt) - Tt)
        print(f"  worst error {err.max():.1f} K at r = {rt[int(np.argmax(err))]*1e3:.0f} mm "
              f"(criterion 60 K) -> {'PASS' if err.max() <= 60 else 'FAIL'}")
        print(f"  T(0) = {T_py(0.0):.0f} K, T(170 mm) = {T_py(0.170):.0f} K, "
              f"T(250 mm) = {T_py(0.250):.0f} K, ambient {T_AMB:g} K "
              f"(guard is ambient; it never sets a physical value)")
        ok = ok and err.max() <= 60
    return ok


def preflight():
    print(f"binary {J.binary_sha()[:8]}…   stop {STOP:g} s   jet_closure = none")
    print(f"corner: RADIAL at r = {R_CORNER*1e3:g} mm (= foot_r_outer; carriers at "
          f"r 28–40 mm are all in the FLOOR zone by construction)")
    print(f"outer bound: h -> {H_OFF:g} beyond r = {R_OUTER*1e3:g} mm;  s <= 0 -> {H_OFF:g}")
    run = AN.Run(CONTROL, AN.meta_of(CONTROL))
    for name in CASES:
        case, mesh = name.rsplit("_", 1)
        _, cmd, m = job(name)
        print(f"\n{name:10s} h_wall {m['h_wall']:g}  dz {m['dz']*1e3:g} mm dt {m['dt']*1e3:g} ms "
              f"stop {m['stop']:g} s")
        if mesh != "2mm":
            continue
        for T in (50.0, 150.0, 300.0):
            j = int(np.argmin(np.abs(run.th["time"] - T)))
            zn = run.th["nozzle_z"][j]
            s = zn - (run.lz - run.depth(T))
            for lab, sel in (("s<=0 (off)", s <= 0.0),
                             ("FLOOR r<=40", (s > 0.0) & (run.r <= R_CORNER)),
                             ("WALL 40<r<=170", (s > 0.0) & (run.r > R_CORNER) & (run.r <= R_OUTER)),
                             ("beyond 170", (s > 0.0) & (run.r > R_OUTER))):
                n = int(sel.sum())
                if n == 0:
                    print(f"    t={T:.0f}s {lab:15s}: 0 columns")
                    continue
                hv = h_py(case, run.r[sel], s[sel])
                Tv = T_py(run.r[sel])
                print(f"    t={T:.0f}s {lab:15s}: {n:5d} cols  r {run.r[sel].min()*1e3:4.0f}"
                      f"–{run.r[sel].max()*1e3:3.0f} mm  h {hv.min():9.3g}–{hv.max():<9.3g} "
                      f"T_gas {Tv.min():6.1f}–{Tv.max():6.1f} K")
    passed = g1()
    print(f"\nG1 radial-share and T-fit gate: {'PASS' if passed else 'FAIL'}")
    return passed


def run_all(names, jobs, force):
    J.keep_awake()
    js = [job(n) for n in names]
    print(f"binary {J.binary_sha()[:8]}…  {len(js)} runs, {min(jobs, len(js))} in parallel "
          f"-> {OUT}", flush=True)
    for _, _, m in js:
        print(f"  {m['name']:10s} h_wall {m['h_wall']:g}; dz {m['dz']*1e3:g} mm, "
              f"stop {m['stop']:g} s", flush=True)
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
        sys.exit(0 if preflight() else 1)
    bad = [n for n in (a.cases or []) if n not in CASES]
    if not a.cases or bad:
        sys.exit(f"give --cases from {CASES}" + (f" (unknown: {bad})" if bad else ""))
    sys.exit(0 if run_all(a.cases, a.jobs, a.force) else 1)


if __name__ == "__main__":
    main()
