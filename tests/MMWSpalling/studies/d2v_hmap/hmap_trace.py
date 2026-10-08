#!/usr/bin/env python3
"""D2v item 3: prove the COMPILED h map from a run, on BOTH branches of the s
switch. Criteria frozen in HMAP_TRACE.md (sha256 in HMAP_TRACE.sha256), written
before any temperature was read.

  hmap_trace.py            score TF (floor branch) and TP (wall branch)
  hmap_trace.py --run      run them first (~10 s each), then score
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
_sp = importlib.util.spec_from_file_location("d2v_run", os.path.join(HERE, "run.py"))
R = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(R)

OUT = R.OUT
RHO, CP, K = 2750.0, 790.0, 1.5
RHOCP = RHO * CP
T0, T_AMB = 293.15, 293.15
EPS, SIG, HC = 0.8, 5.670374e-8, 10.0
BIN = 2.0e-3
R_CLAMP = 0.01875                       # h_imp's own r clamp -> the known plateau


def face_T(h, Tc, G, T_gas):
    lo, hi = min(Tc, T_gas, T_AMB), max(Tc, T_gas, T_AMB)
    for _ in range(80):
        Tf = 0.5 * (lo + hi)
        F = (G * (Tf - Tc) - h * (T_gas - Tf)
             + EPS * SIG * (Tf ** 4 - T_AMB ** 4) + HC * (Tf - T_AMB))
        if F > 0.0:
            hi = Tf
        else:
            lo = Tf
    return 0.5 * (lo + hi)


def column_E(h, dz, dt, nsteps, T_gas, nz=12):
    """Column enthalpy gain per unit area: same dz, dt, face form and explicit
    update as the kernel. The MAP never enters -- only the response to a given h."""
    T = np.full(nz, T0)
    G = 2.0 * K / dz
    for _ in range(nsteps):
        q_in = G * (face_T(h, T[-1], G, T_gas) - T[-1])
        lap = np.zeros(nz)
        lap[1:-1] = K * (T[2:] - 2.0 * T[1:-1] + T[:-2]) / dz ** 2
        lap[0] = K * (T[1] - T[0]) / dz ** 2
        lap[-1] = K * (T[-2] - T[-1]) / dz ** 2 + q_in / dz
        T = T + dt / RHOCP * lap
    return RHOCP * dz * float(np.sum(T - T0))


def invert(E, dz, dt, nsteps, T_gas):
    if E <= 0.0 or T_gas <= T0:
        return 0.0
    lo, hi = 0.0, 1.0e5
    if column_E(hi, dz, dt, nsteps, T_gas) < E:
        return float("nan")
    for _ in range(70):
        mid = 0.5 * (lo + hi)
        if column_E(mid, dz, dt, nsteps, T_gas) < E:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def tgas_of(name, t, edges):
    """Per-bin T_gas at the plotfile time. TP has no marched column (its profile
    csv is header-only), so T_gas there is the constant jet_T_above."""
    p = os.path.join(OUT, name + "_jet_profile.csv")
    rows = np.loadtxt(p, delimiter=",", skiprows=1, ndmin=2)
    centres = 0.5 * (edges[:-1] + edges[1:])
    if rows.size == 0:
        return np.full(len(centres), R.T_ABOVE), "constant jet_T_above (no column marches)"
    ts = np.unique(rows[:, 0])
    tt = ts[np.argmin(np.abs(ts - t))]
    s = rows[:, 0] == tt
    r_lo, Tg = rows[s, 2], rows[s, 4]
    o = np.argsort(r_lo)
    out = np.interp(centres, r_lo[o], Tg[o])
    out[centres < r_lo.min()] = Tg[o][0]
    return out, f"per-bin march T_gas at t = {tt:g} s ({Tg.min():.0f}-{Tg.max():.0f} K)"


def profile(name):
    m = ast.literal_eval(open(os.path.join(OUT, name + ".done")).read())
    pfs = sorted(d for d in os.listdir(os.path.join(OUT, name)) if d.endswith("cell"))
    ds = yt.load(os.path.join(OUT, name, pfs[-1]))
    dims = ds.domain_dimensions
    cg = ds.covering_grid(0, left_edge=ds.domain_left_edge, dims=dims)
    T = np.array(cg["Temp"])
    t = float(ds.current_time)
    dx = float((ds.domain_right_edge[0] - ds.domain_left_edge[0]) / dims[0])
    dz = float((ds.domain_right_edge[2] - ds.domain_left_edge[2]) / dims[2])
    nxy = dims[0]
    c = (np.arange(nxy) + 0.5) * dx
    r = np.hypot(*np.meshgrid(c, c, indexing="ij"))
    E = RHOCP * dz * (T - T0).sum(axis=2)
    dT = (T - T0).max(axis=2)
    nev = sum(1 for _ in open(os.path.join(OUT, name + "_removal_events.csv"))) - 1
    edges = np.arange(0.0, r.max() + BIN, BIN)
    Tg, tg_note = tgas_of(name, t, edges)
    nsteps = int(round(t / m["dt"]))
    rows = []
    for k, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        s = (r >= lo) & (r < hi)
        if not s.any():
            continue
        Em = float(E[s].mean())
        rows.append((0.5 * (lo + hi), int(s.sum()), Em, Tg[k],
                     invert(Em, dz, m["dt"], nsteps, Tg[k]), float(dT[s].max())))
    return dict(m=m, t=t, nsteps=nsteps, dz=dz, r=r, dT=dT, nev=nev, pf=pfs[-1],
                rows=np.array(rows), tg_note=tg_note)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true")
    if ap.parse_args().run:
        subprocess.run([sys.executable, os.path.join(HERE, "run.py"),
                        "--cases", "TF", "TP", "--jobs", "2"], check=True)

    print("D2v item 3 — the compiled h map, BOTH branches of the s switch")
    print(f"criteria: HMAP_TRACE.md sha256 "
          f"{open(os.path.join(HERE, 'HMAP_TRACE.sha256')).read().split()[0][:16]}...")
    print(f"map: w1 {R.R_IMP*1e3:g}->{R.R_ANN_R*1e3:g} mm, h_ann {R.H_ANN:g}, "
          f"w2 {R.R_MOUTH*1e3:g}->{R.R_OUT*1e3:g} mm, h_out {R.H_OUT:g}\n")

    P = {n: profile(n) for n in ("TF", "TP")}
    for n, p in P.items():
        print(f"{n}: {p['pf']} at t = {p['t']:g} s ({p['nsteps']} steps), "
              f"{p['nev']} removal events; T_gas = {p['tg_note']}")
    fails = []
    for n, p in P.items():
        if p["nev"] != 0:
            fails.append(f"T10 {n}: {p['nev']} removal events, must be 0")

    # ---------------- TF: the floor branch ----------------
    p = P["TF"]
    rc, nc, Ec, Tg, hm, dTm = p["rows"].T
    intended = R.h_map(rc, 0.05, above=False)
    peak = float(R.h_map(np.array(0.0), 0.05, above=False))
    print(f"\n--- TF, the FLOOR branch (s = +50 mm), peak intended {peak:.1f} W/m2K ---")
    print("  r[mm]   n    T_gas[K]   h_meas   h_intended    diff    maxdT[K]")
    for i in range(len(rc)):
        if rc[i] > 0.12 and i % 4:
            continue
        print(f"{rc[i]*1e3:7.1f} {int(nc[i]):4d} {Tg[i]:9.1f} {hm[i]:9.2f} "
              f"{intended[i]:11.2f} {hm[i]-intended[i]:+8.2f} {dTm[i]:9.4g}")
    inn = rc < R_CLAMP
    h_flat = float(hm[inn].mean())
    t1 = abs(h_flat / 1449.0 - 1.0)
    print(f"\nT1 self-calibration: h over r < {R_CLAMP*1e3:g} mm = {h_flat:.1f} vs "
          f"1449 known -> {t1*100:+.2f} %  (bar 5 %)")
    if t1 > 0.05:
        fails.append(f"T1: inversion gives {h_flat:.0f} at the known-1449 plateau; "
                     f"the ESTIMATOR is void, not the map")
    d = np.diff(hm)
    t2 = float(d.max())
    print(f"T2 monotone: largest bin-to-bin rise = {t2:+.4g} W/m2K  (bar {0.01*peak:.1f})")
    if t2 > 0.01 * peak:
        fails.append(f"T2: rises by {t2:.3g}")
    t3 = float(np.abs(hm - intended).max()); i3 = int(np.argmax(np.abs(hm - intended)))
    print(f"T3 intended map: max |h_meas - h_intended| = {t3:.2f} at r = "
          f"{rc[i3]*1e3:.1f} mm  (bar {0.05*peak:.1f})")
    if t3 > 0.05 * peak:
        fails.append(f"T3: {t3:.1f} > {0.05*peak:.1f} at r = {rc[i3]*1e3:.1f} mm")
    ann = (rc >= R.R_ANN_R) & (rc <= R.R_MOUTH)
    t4 = float(np.abs(hm[ann] / R.H_ANN - 1.0).max()) if ann.any() else np.nan
    print(f"T4 annulus plateau {R.R_ANN_R*1e3:g}-{R.R_MOUTH*1e3:g} mm: h = "
          f"{hm[ann].mean():.1f} vs {R.H_ANN:g}, max dev {t4*100:.1f} %  (bar 15 %)")
    if t4 > 0.15:
        fails.append(f"T4: annulus off by {t4*100:.0f} %")
    flo = rc >= R.R_OUT
    t5 = float((hm[flo] / peak).max())
    print(f"T5 floor by r_out: max h/peak over r >= {R.R_OUT*1e3:g} mm = {t5:.3g}  (bar 0.01)")
    if t5 > 0.01:
        fails.append(f"T5: {t5:.3g}")
    out = p["r"] > R.R_OUT
    t6 = float(p["dT"][out].max())
    print(f"T6 THE FUNNEL WITNESS: max dT over {int(out.sum())} columns with r > "
          f"{R.R_OUT*1e3:g} mm = {t6:.4g} K  (bar 0.05).  The STALE map has "
          f"h = 362 W/m2K there.")
    if t6 > 0.05:
        fails.append(f"T6: {t6:.3g} K -- the mouth is still heated (N3)")

    # ---------------- TP: the wall branch ----------------
    p = P["TP"]
    rc, nc, Ec, Tg, hm, dTm = p["rows"].T
    intended = R.h_map(rc, -0.01, above=True)
    print(f"\n--- TP, the WALL branch (s = -10 mm), intended flat {R.H_ANN:g} then w2 ---")
    print("  r[mm]   n    T_gas[K]   h_meas   h_intended  blended-would-be   maxdT[K]")
    for i in range(len(rc)):
        if rc[i] > 0.10 and i % 4:
            continue
        bl = R.h_map(np.array(rc[i]), 0.015, above=False)
        print(f"{rc[i]*1e3:7.1f} {int(nc[i]):4d} {Tg[i]:9.1f} {hm[i]:9.2f} "
              f"{intended[i]:11.2f} {bl:17.1f} {dTm[i]:9.4g}")
    inn = rc < R.R_ANN_R
    t7 = float(np.abs(hm[inn] / R.H_ANN - 1.0).max())
    print(f"\nT7 wall branch flat: h over r < {R.R_ANN_R*1e3:g} mm = "
          f"{hm[inn].mean():.2f} vs {R.H_ANN:g}, max dev {t7*100:.1f} %  (bar 15 %)")
    if t7 > 0.15:
        fails.append(f"T7: wall branch off by {t7*100:.0f} %")
    bl = (rc >= R.R_IMP) & (rc <= R.R_ANN_R)
    t8 = float(hm[bl].max())
    print(f"T8 NOT blended: max h over {R.R_IMP*1e3:g}-{R.R_ANN_R*1e3:g} mm = "
          f"{t8:.1f}  (bar {1.3*R.H_ANN:.1f}); a blend would read "
          f"{float(R.h_map(np.array(0.044), 0.015, above=False)):.0f}-"
          f"{float(R.h_map(np.array(0.040), 0.015, above=False)):.0f}")
    if t8 > 1.3 * R.H_ANN:
        fails.append(f"T8: {t8:.0f} -- the wall branch is BLENDED, not flat")
    flo = rc >= R.R_OUT
    t9a = float(hm[flo].max() / R.H_ANN)
    out = p["r"] > R.R_OUT
    t9b = float(p["dT"][out].max())
    print(f"T9 floor/no rise beyond r_out: max h/h_ann = {t9a:.3g} (bar 0.01), "
          f"max dT = {t9b:.4g} K (bar 0.05)")
    if t9a > 0.01 or t9b > 0.05:
        fails.append(f"T9: {t9a:.3g} / {t9b:.3g} K")

    print()
    if fails:
        print("*** STOP — the compiled map is not the intended one ***")
        for f in fails:
            print("  FAIL " + f)
        return 1
    print("ALL PASS (T1-T10): the compiled h map is the intended one on BOTH "
          "branches of the s switch; the wall branch is flat h_ann and NOT "
          "blended; nothing is heated beyond the mouth.")
    print("Cleared to launch the campaign.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
