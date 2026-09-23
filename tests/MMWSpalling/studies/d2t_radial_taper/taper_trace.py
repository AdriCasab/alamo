#!/usr/bin/env python3
"""D2t Goal item 2: prove the COMPILED taper from a run. Criteria and method are
frozen in TAPER_TRACE.md (sha256 in TAPER_TRACE.sha256), written before any
temperature was read.

  taper_trace.py            score the two trace legs
  taper_trace.py --run      run them first (10 s each), then score

flame_h is never written anywhere, so the temperature field is the only witness
to h(r). Ten 16 ms steps from a uniform 293.15 K, before any spall, with every
column on the face form: the column enthalpy gain traces h(r) directly.
"""
import argparse
import ast
import importlib.util
import os
import subprocess
import sys

import numpy as np
import yt

yt.set_log_level(50)

HERE = os.path.dirname(os.path.abspath(__file__))
_sp = importlib.util.spec_from_file_location("d2t_run", os.path.join(HERE, "run.py"))
R = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(R)

OUT = R.OUT
RHO, CP, K = 2750.0, 790.0, 1.5
RHOCP = RHO * CP
T0, T_AMB, T_GAS = 293.15, 293.15, 1600.0
EPS, SIG, HC = 0.8, 5.670374e-8, 10.0
BIN = 2.0e-3                                   # fixed physical radial bin


def face_T(h, Tc, G):
    """The kernel's face-form surface balance (SurfaceCellFlux, :2648-2667),
    solved by bisection instead of its damped Newton. Same equation."""
    lo, hi = min(Tc, T_GAS, T_AMB), max(Tc, T_GAS, T_AMB)
    for _ in range(80):
        Tf = 0.5 * (lo + hi)
        F = (G * (Tf - Tc) - h * (T_GAS - Tf)
             + EPS * SIG * (Tf ** 4 - T_AMB ** 4) + HC * (Tf - T_AMB))
        if F > 0.0:
            hi = Tf
        else:
            lo = Tf
    return 0.5 * (lo + hi)


def column_E(h, dz, dt, nsteps, nz=12):
    """Column enthalpy gain per unit area after nsteps, for a given h. Same dz,
    dt, face form and explicit update as the kernel; the taper never enters."""
    T = np.full(nz, T0)
    G = 2.0 * K / dz
    for _ in range(nsteps):
        Tf = face_T(h, T[-1], G)
        q_in = G * (Tf - T[-1])
        lap = np.zeros(nz)
        lap[1:-1] = K * (T[2:] - 2.0 * T[1:-1] + T[:-2]) / dz ** 2
        lap[0] = K * (T[1] - T[0]) / dz ** 2
        lap[-1] = K * (T[-2] - T[-1]) / dz ** 2 + q_in / dz
        T = T + dt / RHOCP * lap
    return RHOCP * dz * float(np.sum(T - T0))


def invert(E, dz, dt, nsteps):
    """h from a measured column enthalpy gain. Monotone in h, so bisect."""
    if E <= 0.0:
        return 0.0
    lo, hi = 0.0, 1.0e5
    if column_E(hi, dz, dt, nsteps) < E:
        return float("nan")
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if column_E(mid, dz, dt, nsteps) < E:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def read(name):
    m = ast.literal_eval(open(os.path.join(OUT, name + ".done")).read())
    pfs = sorted(d for d in os.listdir(os.path.join(OUT, name)) if d.endswith("cell"))
    ds = yt.load(os.path.join(OUT, name, pfs[-1]))
    dims = ds.domain_dimensions
    cg = ds.covering_grid(0, left_edge=ds.domain_left_edge, dims=dims)
    T = np.array(cg["Temp"])
    t = float(ds.current_time)
    dz = float((ds.domain_right_edge[2] - ds.domain_left_edge[2]) / dims[2])
    nxy = dims[0]
    c = (np.arange(nxy) + 0.5) * dz
    r = np.hypot(*np.meshgrid(c, c, indexing="ij"))
    E = RHOCP * dz * (T - T0).sum(axis=2)
    dT = (T - T0).max(axis=2)
    nev = sum(1 for _ in open(os.path.join(OUT, name + "_removal_events.csv"))) - 1
    return m, t, dz, r, E, dT, nev, pfs[-1]


def profile(name):
    m, t, dz, r, E, dT, nev, pf = read(name)
    nsteps = int(round(t / m["dt"]))
    edges = np.arange(0.0, r.max() + BIN, BIN)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (r >= lo) & (r < hi)
        if not s.any():
            continue
        Em = float(E[s].mean())
        rows.append((0.5 * (lo + hi), int(s.sum()), Em,
                     invert(Em, dz, m["dt"], nsteps), float(dT[s].max())))
    return m, t, nsteps, dz, r, E, dT, nev, pf, np.array(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true")
    a = ap.parse_args()
    if a.run:
        subprocess.run([sys.executable, os.path.join(HERE, "run.py"),
                        "--cases", "TR_below", "TR_above", "--jobs", "2"], check=True)

    print("D2t item 2 — the compiled-taper thermal trace")
    print(f"criteria: TAPER_TRACE.md sha256 "
          f"{open(os.path.join(HERE, 'TAPER_TRACE.sha256')).read().split()[0][:16]}...")
    print(f"taper: h = 700*taper below the plane, 142*taper above; flat to "
          f"{R.R_FLAT * 1e3:g} mm, floor {R.H_FLOOR:g} by {R.R_TAPER * 1e3:g} mm; "
          f"patch radius {R.R_PATCH * 1e3:g} mm\n")

    P = {}
    for n in ("TR_below", "TR_above"):
        m, t, ns, dz, r, E, dT, nev, pf, rows = profile(n)
        P[n] = dict(m=m, t=t, ns=ns, dz=dz, r=r, E=E, dT=dT, nev=nev, pf=pf, rows=rows)
        print(f"{n}: {pf} at t = {t:g} s ({ns} steps of {m['dt'] * 1e3:g} ms), "
              f"dz {dz * 1e3:g} mm, nozzle_z0 {[c for c in R.job(n)[1] if 'nozzle_z0' in c][0].split('=')[1]}, "
              f"{nev} removal events")

    fails, notes = [], []

    # T8 (already confirmed before hashing; re-asserted here)
    for n in P:
        if P[n]["nev"] != 0:
            fails.append(f"T8 {n}: {P[n]['nev']} removal events, must be 0")

    b = P["TR_below"]["rows"]
    peak_ref = R.H_IN
    rc, nc, Ec, hc_, dTc = b.T
    tap = R.taper_of(rc)

    print(f"\nTR_below — the s > 0 branch, h = 700*taper")
    print(" r[mm]   n      E[J/m2]   h_meas   h/700   taper(r)   diff    maxdT[K]")
    for i in range(len(rc)):
        mark = "  <- patch edge" if abs(rc[i] - R.R_PATCH) < BIN else ""
        d = hc_[i] / peak_ref - (tap[i] if rc[i] <= R.R_PATCH else 0.0)
        print(f"{rc[i]*1e3:6.1f} {int(nc[i]):4d} {Ec[i]:12.4g} {hc_[i]:8.2f} "
              f"{hc_[i]/peak_ref:8.5f} {tap[i]:9.5f} {d:+8.5f} {dTc[i]:9.4g}{mark}")

    inn = rc < R.R_FLAT
    h_flat = float(hc_[inn].mean())
    t1 = abs(h_flat / peak_ref - 1.0)
    print(f"\nT1 self-calibration: h over r < {R.R_FLAT*1e3:g} mm = {h_flat:.2f} "
          f"W/m2K vs 700 known -> {t1*100:+.2f} %  (bar 3 %)")
    if t1 > 0.03:
        fails.append(f"T1: inversion returns {h_flat:.1f} at the known-700 plateau "
                     f"({t1*100:.1f} % off); the ESTIMATOR is void, not the taper")

    t2 = float(np.abs(hc_[inn] / peak_ref - 1.0).max())
    print(f"T2 flat inside r_flat: max |h/700 - 1| = {t2:.4f}  (bar 0.02)")
    if t2 > 0.02:
        fails.append(f"T2: {t2:.4f} > 0.02")

    ins = rc <= R.R_PATCH
    d = np.diff(hc_[ins])
    t3 = float(d.max()) if len(d) else 0.0
    print(f"T3 monotone to radius: largest rise bin-to-bin = {t3:+.4g} W/m2K "
          f"(bar {0.005*peak_ref:g})")
    if t3 > 0.005 * peak_ref:
        fails.append(f"T3: rises by {t3:.3g} W/m2K")

    t4 = float(np.abs(hc_[ins] / peak_ref - tap[ins]).max())
    i4 = int(np.argmax(np.abs(hc_[ins] / peak_ref - tap[ins])))
    print(f"T4 is it the intended profile: max |h/700 - taper| = {t4:.5f} at "
          f"r = {rc[ins][i4]*1e3:.1f} mm  (bar 0.02)")
    if t4 > 0.02:
        fails.append(f"T4: {t4:.5f} > 0.02 at r = {rc[ins][i4]*1e3:.1f} mm")

    flo = (rc >= R.R_TAPER) & (rc <= R.R_PATCH)
    t5 = float((hc_[flo] / peak_ref).max()) if flo.any() else 0.0
    print(f"T5 at the floor by radius: max h/700 over r >= r_taper = {t5:.3g}  (bar 0.01)")
    if t5 > 0.01:
        fails.append(f"T5: {t5:.3g} > 0.01")

    rr, dd = P["TR_below"]["r"], P["TR_below"]["dT"]
    out = rr > R.R_PATCH
    t6 = float(dd[out].max()) if out.any() else 0.0
    print(f"T6 no rise beyond radius: max dT over {int(out.sum())} columns with "
          f"r > {R.R_PATCH*1e3:g} mm = {t6:.3g} K  (bar 1e-3)")
    if t6 > 1.0e-3:
        fails.append(f"T6: {t6:.3g} K > 1e-3 -- the cliff MOVED, it was not removed (N2)")

    a_ = P["TR_above"]["rows"]
    ra, ha = a_[:, 0], a_[:, 3]
    sel = ra < R.R_FLAT
    ratio = ha[sel] / hc_[rc < R.R_FLAT]
    want = R.H_ABOVE / R.H_IN
    t7 = float(np.abs(ratio / want - 1.0).max())
    print(f"\nTR_above — the s <= 0 branch, h = 142*taper")
    print(f"T7 branch association: h_above/h_below over r < r_flat = "
          f"{ratio.mean():.6f} (want {want:.6f}), max deviation {t7*100:.2f} %  (bar 2 %)")
    for i in np.nonzero(sel)[0]:
        print(f"   r {ra[i]*1e3:5.1f} mm: h_above {ha[i]:8.3f}  h_below "
              f"{hc_[rc < R.R_FLAT][i]:8.2f}  ratio {ratio[i]:.6f}")
    if t7 > 0.02:
        fails.append(f"T7: ratio off by {t7*100:.1f} % -- the taper is not applied "
                     f"to both branches of the if()")

    print()
    if fails:
        print("*** STOP — the compiled taper is not the intended profile ***")
        for f in fails:
            print("  FAIL " + f)
        return 1
    print("ALL PASS (T1-T8): the compiled h(r) is the intended taper, on both "
          "branches of the s switch, and nothing is heated beyond the patch radius.")
    print("Cleared to launch the campaign.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
