#!/usr/bin/env python3
"""D2w scoring: the 2x2 attribution. Criteria frozen in CRITERION.md (sha256 in
CRITERION.sha256), hashed before either new leg launched.

  analyze.py [--no-thermal]     --no-thermal skips the yt plotfile probe

The metric is IMPORTED from d2v_hmap/analyze.py, unchanged. One metric across
all four cells or the factorial is meaningless.

ROP is the slope of nozzle_z over 150-400 s -- burner descent, NEVER pit depth.
The burner rides the feet, not the deepest point; quoting pit depth as ROP is
the error D2v made.
"""
import argparse
import ast
import importlib.util
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)


def load_mod(name, path):
    sp = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


V = load_mod("d2v_an", os.path.join(STUD, "d2v_hmap", "analyze.py"))   # THE metric
W = load_mod("d2w_run", os.path.join(HERE, "run.py"))

OUT = {"M0": V.OUT, "M1": V.OUT, "A": W.OUT, "B": W.OUT}
CELLS = ("M0", "A", "B", "M1")
LAB = {"M0": "M0  map no  side no", "A": "A   map YES side no",
       "B": "B   map no  side YES", "M1": "M1  map YES side YES"}
ROP_WIN = (150.0, 400.0)
R_PAD = (0.028, 0.040)
R_IN = 0.076
V_MEIER_L, V_MEIER_M = 3.42, 0.5


def cell(n):
    d = OUT[n]
    m = ast.literal_eval(open(os.path.join(d, n + ".done")).read())
    dz, nxy = m["dz"], m["nxy"]
    c = (np.arange(nxy) + 0.5) * dz
    r = np.hypot(*np.meshgrid(c, c, indexing="ij"))
    a = np.loadtxt(os.path.join(d, n + "_removal_events.csv"), delimiter=",",
                   skiprows=1, ndmin=2)
    return m, r, a, d


def thermo(n, col):
    p = os.path.join(OUT[n], n, "thermo.dat")
    h = open(p).readline().split()
    if col not in h:
        return None
    a = np.loadtxt(p, skiprows=1, ndmin=2)
    a = a[a[:, 0] > 0]
    return a[:, 0], a[:, h.index(col)]


def rop(n):
    """Burner descent [m/h] = -slope of nozzle_z over 150-400 s."""
    t, z = thermo(n, "nozzle_z")
    s = (t >= ROP_WIN[0]) & (t <= ROP_WIN[1])
    return -np.polyfit(t[s], z[s], 1)[0] * 3600.0


def pmean(n, col):
    v = thermo(n, col)
    if v is None:
        return None
    t, x = v
    s = (t >= ROP_WIN[0]) & (t <= ROP_WIN[1])
    return float(x[s].mean())


def thermal_step(n):
    """Max lateral dT between a solid top cell and a TALLER solid neighbour,
    inside r = 76 mm, at the last plotfile. D2o-0b's cascade witness."""
    import yt
    yt.set_log_level(50)
    d = os.path.join(OUT[n], n)
    pfs = sorted(x for x in os.listdir(d) if x.endswith("cell"))
    ds = yt.load(os.path.join(d, pfs[-1]))
    dims = ds.domain_dimensions
    cg = ds.covering_grid(0, left_edge=ds.domain_left_edge, dims=dims)
    T = np.array(cg["Temp"]); rem = np.array(cg["removed"]) > 0.5
    dx = float((ds.domain_right_edge[0] - ds.domain_left_edge[0]) / dims[0])
    nx = dims[0]
    ktop = np.full((nx, nx), -1)
    for i in range(nx):
        for j in range(nx):
            s = np.nonzero(~rem[i, j, :])[0]
            if len(s):
                ktop[i, j] = s.max()
    c = (np.arange(nx) + 0.5) * dx
    r = np.hypot(*np.meshgrid(c, c, indexing="ij"))
    st = []
    for i in range(nx):
        for j in range(nx):
            k = ktop[i, j]
            if k < 0 or r[i, j] > R_IN:
                continue
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                ii, jj = i + di, j + dj
                if 0 <= ii < nx and 0 <= jj < nx and ktop[ii, jj] > k and not rem[ii, jj, k]:
                    st.append(T[i, j, k] - T[ii, jj, k])
    st = np.array(st) if st else np.array([0.0])
    return float(st.max()), float(np.median(st)), len(st), float(ds.current_time)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--no-thermal", action="store_true")
    no_th = ap.parse_args().no_thermal
    print("D2w — complete the 2x2: is the near-cylinder the MAP or the SIDE FACES?")
    print(f"criteria CRITERION.md sha256 "
          f"{open(os.path.join(HERE,'CRITERION.sha256')).read().split()[0][:16]}...")
    print("metric IMPORTED unchanged from d2v_hmap/analyze.py\n")
    missing = [n for n in CELLS if not os.path.exists(os.path.join(OUT[n], n + ".done"))]
    if missing:
        print(f"not finished: {missing}"); return 1

    C = {n: cell(n) for n in CELLS}
    T_COMMON = min(C[n][0]["t_end"] for n in CELLS)
    print(f"--- end time and depth reached (N3) ---")
    for n in CELLS:
        m, r, a, _ = C[n]
        zf, _ = V.zface(m, r, a, m["t_end"])
        print(f"  {LAB[n]}: status {m.get('status','?'):8s} t_end {m['t_end']:7.1f} s, "
              f"pit depth {(V.LZ - zf.min())*1e3:6.1f} mm")
    print(f"  -> common scoring time {T_COMMON:.1f} s\n")

    P, Z = {}, {}
    for n in CELLS:
        m, r, a, _ = C[n]
        zf, nv = V.zface(m, r, a, T_COMMON)
        Z[n] = (zf, nv); P[n] = V.profile(r, zf)

    # ---------- the ladder ----------
    print(f"--- diameter vs depth at t = {T_COMMON:.1f} s [mm] ---")
    D = {n: {d: V.diameter(*P[n][:2], d) * 1e3 for d in V.DEPTHS} for n in CELLS}
    print("depth  Meier*" + "".join(f"{n:>9s}" for n in CELLS))
    for d in V.DEPTHS:
        print(f"{d:5d} {V.MEIER[d]:6d}" + "".join(f"{D[n][d]:9.1f}" for n in CELLS))
    print("* a VISIBLE WIDTH = a LOWER BOUND on diameter, +/-5 mm; shape only.")

    # ---------- (a) ----------
    print("\n--- (a) SCORED — WHO OWNS THE FUNNEL: mean |Ø - Meier| over 10/25/50 mm ---")
    A_ = {}
    for n in CELLS:
        e = [abs(D[n][d] - V.MEIER[d]) for d in V.A_DEPTHS]
        A_[n] = float(np.mean(e))
        print(f"  {LAB[n]}: excess {' / '.join(f'{x:5.1f}' for x in e)}  -> mean {A_[n]:6.2f} mm")
    dA = abs(A_["A"] - A_["M1"])
    print(f"  P3: |A - M1| = {dA:.2f} mm  -> {'THE MAP OWNS THE FUNNEL' if dA <= 5.0 else 'the map does NOT fully own it'}"
          f"  (within 5 mm ⇒ map owns it; refuted if A > 15 mm: A = {A_['A']:.2f})")

    # ---------- (b) ----------
    print("\n--- (b) SCORED — WHO OWNS THE RATE: burner descent (nozzle_z slope, "
          "150-400 s) [m/h] ---")
    Rp = {n: rop(n) for n in CELLS}
    for n in CELLS:
        band = "IN Meier 1.3-1.6" if 1.3 <= Rp[n] <= 1.6 else "OUT of band"
        print(f"  {LAB[n]}: {Rp[n]:5.2f} m/h   {band}")
    okb = abs(Rp["B"] - Rp["M1"]) < abs(Rp["B"] - Rp["M0"])
    print(f"  bar |B - M1| < |B - M0|: {abs(Rp['B']-Rp['M1']):.2f} < "
          f"{abs(Rp['B']-Rp['M0']):.2f}  -> {'PASS (the side path owns the rate)' if okb else 'FAIL'}")
    print(f"  P2 (B > 2.5, refuted < 2.0): B = {Rp['B']:.2f} -> "
          f"{'HELD' if Rp['B'] > 2.5 else ('REFUTED' if Rp['B'] < 2.0 else 'partial')}")
    lo, up = Rp["A"] < 1.3, Rp["A"] > 2.0
    print(f"  P1 (A two-sided): A = {Rp['A']:.2f} -> "
          + ("LOWER BRANCH: the map alone COSTS rate" if lo else
             ("UPPER BRANCH: map is NOT rate-neutral" if up else
              "NEITHER branch: A is between 1.3 and 2.0, i.e. near M0")))

    # ---------- (c) ----------
    print("\n--- (c) WHO OWNS THE FLAT FLOOR: median radial step in z(r), r < 76 mm ---")
    S = {}
    for n in CELLS:
        rc, zc, nc = P[n]
        s_in = rc <= R_IN
        d = np.abs(np.diff(zc[s_in])) * 1e3
        S[n] = float(np.median(d))
        print(f"  {LAB[n]}: median {S[n]:5.2f} mm, max {d.max():6.1f} mm")
    closer = "A" if abs(S["A"] - S["M1"]) < abs(S["B"] - S["M1"]) else "B"
    print(f"  P4 (flat floor belongs to B): closer to M1's {S['M1']:.2f} is {closer} "
          f"-> {'HELD' if closer == 'B' else 'REFUTED — the MAP flattens the floor'}")

    # ---------- (d) ----------
    print("\n--- (d) ADDITIVITY: interaction = M1 - (A + B - M0) ---")
    ia = A_["M1"] - (A_["A"] + A_["B"] - A_["M0"])
    ib = Rp["M1"] - (Rp["A"] + Rp["B"] - Rp["M0"])
    print(f"  on (a): {A_['M1']:.2f} - ({A_['A']:.2f} + {A_['B']:.2f} - {A_['M0']:.2f}) "
          f"= {ia:+.2f} mm")
    print(f"  on (b): {Rp['M1']:.2f} - ({Rp['A']:.2f} + {Rp['B']:.2f} - {Rp['M0']:.2f}) "
          f"= {ib:+.2f} m/h")
    print("  A large interaction is itself a result: the mechanisms are not separable.")

    # ---------- (e) ----------
    print("\n--- (e) POWER, mean over 150-400 s [W] — THREE COLUMNS, NOT TWO ---")
    print("  cell                 tops   sides  s_loss   TOTAL  identity")
    for n in CELLS:
        tp = pmean(n, "patch_P_robin"); sf = pmean(n, "side_P_face")
        sl = pmean(n, "side_P_loss"); tot = pmean(n, "jet_P_face")
        f = lambda x: "  absent" if x is None else f"{x:8.1f}"
        idn = (tot - (tp or 0) - (sf or 0) - (sl or 0)) if tot is not None else None
        print(f"  {LAB[n]}: {f(tp)}{f(sf)}{f(sl)}{f(tot)}  "
              f"{'' if idn is None else f'{idn:+.2f} W'}")
    print("  jet_P_face is the TOTAL drawn from the stream, not the tops: the march "
          "gathers side faces too.")

    # ---------- (f) ----------
    print("\n--- (f) VOLUME ON BOTH BASES (the basis changes the verdict) ---")
    print("  cell                 vol[L]  pit[mm]  x(pit)  descent[mm]  x(descent)")
    for n in CELLS:
        m, r, a, _ = C[n]
        nv = Z[n][1]
        vol = float(nv.sum()) * m["dz"] ** 3 * 4.0 * 1e3
        pit = float(V.LZ - Z[n][0].min())
        t, z = thermo(n, "nozzle_z")
        desc = float(np.interp(0.0, [0, 1], [0, 0]) + (z[np.argmin(np.abs(t-0.0))] - np.interp(T_COMMON, t, z)))
        e_pit = V_MEIER_L * pit / V_MEIER_M
        e_des = V_MEIER_L * desc / V_MEIER_M if desc > 0 else np.nan
        print(f"  {LAB[n]}: {vol:7.3f} {pit*1e3:8.1f} {vol/e_pit:7.2f} "
              f"{desc*1e3:12.1f} {vol/e_des:11.2f}")

    # ---------- P1 witnesses ----------
    print("\n--- P1 WITNESSES (required): the pad ring, the feet, the wall ---")
    for n in CELLS:
        m, r, a, _ = C[n]
        s = (a[:, 0] >= ROP_WIN[0]) & (a[:, 0] <= ROP_WIN[1])
        rec = np.zeros_like(r)
        np.add.at(rec, (a[s, 1].astype(int), a[s, 2].astype(int)), a[s, 9])
        pad = (r >= R_PAD[0]) & (r < R_PAD[1])
        pr = rec[pad].mean() / (ROP_WIN[1] - ROP_WIN[0]) * 3600.0
        fs = pmean(n, "foot_stall_time"); fc = pmean(n, "foot_cols")
        fcar = pmean(n, "foot_carry_cols")
        print(f"  {LAB[n]}: pad ring 28-40 mm recession {pr:5.2f} m/h;  "
              f"foot_stall_time {fs:5.2f} s, foot_cols {fc:5.0f}, carry {fcar:5.1f}")
    if not no_th:
        print("  lateral dT at the wall (D2o-0b's cascade witness):")
        for n in CELLS:
            try:
                mx, md, cnt, tt = thermal_step(n)
                print(f"    {LAB[n]}: max {mx:6.1f} K, median {md:6.1f} K over {cnt} "
                      f"step faces (t = {tt:.0f} s)")
            except Exception as e:
                print(f"    {LAB[n]}: plotfile probe failed: {e}")

    # ---------- N witnesses ----------
    print("\n--- inert-table witnesses ---")
    for n in ("A", "B"):
        sf = thermo(n, "side_P_face")
        sa = thermo(n, "side_P_face_above")
        print(f"  {n}: side_P_face {'ABSENT (N1 ok)' if sf is None else f'present max {sf[1].max():.1f} W'}; "
              f"side_P_face_above " + ("ABSENT" if sa is None else
              f"max {sa[1].max():+.2f} min {sa[1].min():+.2f} W "
              f"{'(N4 ok: not cooling)' if sa[1].min() >= 0 else '*** NEGATIVE: COOLING (N4) ***'}"))
    hA = [c for c in W.job("A")[1] if c.startswith("surface_patch.h_expr=")][0]
    hB = [c for c in W.job("B")[1] if c.startswith("surface_patch.h_expr=")][0]
    print(f"  N5: A's map == M1's verbatim: {hA == W.map_token()}")
    print(f"  N2: B's map == M0's: {hB == [c for c in W.base_argv() if c.startswith('surface_patch.h_expr=')][0]}")
    for n in CELLS:
        le = thermo(n, "ledger_err")
        print(f"  {n}: max|ledger_err| = {np.abs(le[1]).max():.2e}")

    # ---------- verdict ----------
    print("\n=== VERDICT — the attribution ===")
    print(f"  FUNNEL (a): M0 {A_['M0']:.1f} -> A {A_['A']:.1f} | B {A_['B']:.1f} | M1 {A_['M1']:.1f} mm")
    print(f"  RATE   (b): M0 {Rp['M0']:.2f} -> A {Rp['A']:.2f} | B {Rp['B']:.2f} | M1 {Rp['M1']:.2f} m/h")
    print(f"  FLOOR  (c): M0 {S['M0']:.2f} -> A {S['A']:.2f} | B {S['B']:.2f} | M1 {S['M1']:.2f} mm")
    print(f"\n  (a) P3 {'HELD' if A_['A'] <= 15.0 else 'REFUTED'};  "
          f"(b) bar {'PASS' if okb else 'FAIL'};  "
          f"P4 {'HELD' if closer == 'B' else 'REFUTED'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
