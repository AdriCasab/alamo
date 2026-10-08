#!/usr/bin/env python3
"""D2x scoring. Criteria frozen in CRITERION.md (sha256 in CRITERION.sha256),
hashed before the first leg launched.

SCORING IS AT EQUAL BURNER DEPTH, NEVER AT EQUAL TIME. Each leg is sampled at
the first time its burner descent reaches the target, and only depths at least
50 mm above the pads are scored -- the pads start at the rock surface, so the
cap is (descent - 50 mm). Equal-time scoring across legs with different rates
compares unfinished pit against finished wall; that is the error D2w made.

`A` is the CONTROL THAT ALREADY PASSES. The bar is not "beat M0": a leg
succeeds only by closing A's wall-band narrowing without losing A's shape, rate,
pit lead or volume.
"""
import argparse
import ast
import filecmp
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
D2V_OUT = os.path.join(STUD, "d2v_hmap", "output")
D2W_OUT = os.path.join(STUD, "d2w_attribution", "output")
OUT = os.path.join(HERE, "output")

WHERE = {"M0": D2V_OUT, "M1": D2V_OUT, "A": D2W_OUT, "B": D2W_OUT,
         "C0": OUT, "C": OUT, "D": OUT, "E": OUT}
CELLS = ("M0", "A", "B", "M1", "C", "D", "E")
LAB = {"M0": "M0 stale map, no side", "A": "A  map, no side      ",
       "B": "B  stale map, f=0.2  ", "M1": "M1 map, f=0.2 all    ",
       "C": "C  map, f=0.2 r>=38  ", "D": "D  map, exact all    ",
       "E": "E  map, exact r>=38  "}
LZ = 0.4
TARGETS = ((150, (10, 25, 50, 75, 100)), (180, (10, 25, 50, 75, 100, 125)))
PRIMARY = 150
REF_D = 87.0                     # Meier lower bound at 100-125 mm
BAND = (1.3, 1.6)
V_L, V_M = 3.42, 0.5
RINGS = list(range(30, 50, 2))
ROP_DESC = (0.100, 0.150)        # measure descent slope over the same stretch


def have(n):
    return os.path.exists(os.path.join(WHERE[n], n + ".done"))


def thermo(n, col):
    p = os.path.join(WHERE[n], n, "thermo.dat")
    h = open(p).readline().split()
    if col not in h:
        return None
    a = np.loadtxt(p, skiprows=1, ndmin=2)
    a = a[a[:, 0] > 0]
    return a[:, 0], a[:, h.index(col)]


def descent(n):
    t, z = thermo(n, "nozzle_z")
    return t, z[0] - z


def t_at(n, d_mm):
    t, d = descent(n)
    return (float(np.interp(d_mm * 1e-3, d, t)) if d.max() >= d_mm * 1e-3 else None)


def cell(n):
    m = ast.literal_eval(open(os.path.join(WHERE[n], n + ".done")).read())
    dz, nxy = m["dz"], m["nxy"]
    c = (np.arange(nxy) + 0.5) * dz
    r = np.hypot(*np.meshgrid(c, c, indexing="ij"))
    a = np.loadtxt(os.path.join(WHERE[n], n + "_removal_events.csv"), delimiter=",",
                   skiprows=1, ndmin=2)
    return m, r, a


def shape(n, d_mm, depths):
    """(a), the ladder and the geometry at equal burner descent."""
    t = t_at(n, d_mm)
    if t is None:
        return None
    m, r, a = cell(n)
    zf, nv = V.zface(m, r, a, t)
    rc, zc, _ = V.profile(r, zf)
    cap = d_mm - 50.0
    use = [d for d in depths if d <= cap]
    D = {d: V.diameter(rc, zc, d) * 1e3 for d in depths}
    e = [abs(D[d] - V.MEIER[d]) for d in use]
    pit = (LZ - zf.min()) * 1e3
    vol = float(nv.sum()) * m["dz"] ** 3 * 4.0 * 1e3
    return dict(t=t, D=D, used=use, a=float(np.mean(e)) if e else None,
                pit=pit, lead=pit - d_mm, vol=vol,
                x_pit=vol / (V_L * pit * 1e-3 / V_M),
                x_desc=vol / (V_L * d_mm * 1e-3 / V_M))


def rop(n):
    """Burner descent [m/h] over the SAME physical stretch on every leg."""
    t, d = descent(n)
    s = (d >= ROP_DESC[0]) & (d <= ROP_DESC[1])
    if s.sum() < 3:
        return float("nan")
    return np.polyfit(t[s], d[s], 1)[0] * 3600.0


def rings(n, win=(150.0, 400.0)):
    m, r, a = cell(n)
    s = (a[:, 0] >= win[0]) & (a[:, 0] <= win[1])
    rec = np.zeros_like(r)
    np.add.at(rec, (a[s, 1].astype(int), a[s, 2].astype(int)), a[s, 9])
    return [float(rec[(r >= lo * 1e-3) & (r < (lo + 2) * 1e-3)].mean() * 1e3) for lo in RINGS]


def pmean(n, col, win=(150.0, 400.0)):
    v = thermo(n, col)
    if v is None:
        return None
    t, x = v
    s = (t >= win[0]) & (t <= win[1])
    return float(x[s].mean()) if s.any() else float(x.mean())


def main():
    ap = argparse.ArgumentParser(); ap.parse_args()
    print("D2x — close A's wall band without losing its rate")
    print(f"criteria CRITERION.md sha256 "
          f"{open(os.path.join(HERE,'CRITERION.sha256')).read().split()[0][:16]}...")
    print("metric imported unchanged from d2v_hmap/analyze.py; "
          "SCORED AT EQUAL BURNER DEPTH\n")
    missing = [n for n in CELLS if not have(n)]
    if missing:
        print(f"not finished: {missing}"); return 1

    # ---------- C0: the identity witness the tree gate cannot replace ----------
    print("--- C0: both new keys at their defaults, vs M1, BYTE FOR BYTE to 200 s ---")
    ok_c0 = True
    for lab, mine, ref, sep in (
            # 199.9, not 200.0: C0's stop_time lands exactly on a step, so its
            # final step at t = 199.99999999997618 is never produced while M1
            # (which ran on to 401 s) has it. C0 is a strict PREFIX of M1; the
            # comparison is over the interval both legs completed.
            ("removal_events.csv", os.path.join(OUT, "C0_removal_events.csv"),
             os.path.join(D2V_OUT, "M1_removal_events.csv"), ","),
            ("thermo.dat", os.path.join(OUT, "C0", "thermo.dat"),
             os.path.join(D2V_OUT, "M1", "thermo.dat"), None)):
        A_ = V.E.rows_upto(mine, 199.9, sep=sep)
        B_ = V.E.rows_upto(ref, 199.9, sep=sep)
        same = A_ == B_
        ok_c0 &= same
        print(f"  {lab:20s} {len(A_)-1:6d} rows vs {len(B_)-1:6d}  "
              f"{'BYTE-IDENTICAL' if same else '*** DIFFER ***'}")
    print(f"  -> new keys inert at their defaults on a side-path run: "
          f"{'PASS' if ok_c0 else 'FAIL'}")

    # ---------- the ladder at equal descent ----------
    S = {}
    for d_mm, depths in TARGETS:
        print(f"\n--- diameter vs depth at BURNER DESCENT {d_mm} mm "
              f"(depths capped at {d_mm-50} mm) [mm] ---")
        S[d_mm] = {n: shape(n, d_mm, depths) for n in CELLS}
        hdr = "depth Meier*" + "".join(f"{n:>8s}" for n in CELLS)
        print(hdr)
        for dep in depths:
            row = f"{dep:5d} {V.MEIER[dep]:6d}"
            for n in CELLS:
                s = S[d_mm][n]
                row += "     n/a" if s is None else f"{s['D'][dep]:8.1f}"
            print(row + ("" if dep <= d_mm - 50 else "   (below the plane: not scored)"))
        print("  (a)        " + "".join(
            "     n/a" if S[d_mm][n] is None or S[d_mm][n]["a"] is None
            else f"{S[d_mm][n]['a']:8.2f}" for n in CELLS))
        print("  t [s]      " + "".join(
            "     n/a" if S[d_mm][n] is None else f"{S[d_mm][n]['t']:8.1f}" for n in CELLS))

    # ---------- (c) rate ----------
    print("\n--- (c) burner descent over the 100-150 mm stretch [m/h], and pit lead ---")
    Rp = {n: rop(n) for n in CELLS}
    for n in CELLS:
        s150, s180 = S[150][n], S[180][n]
        band = "IN 1.3-1.6" if BAND[0] <= Rp[n] <= BAND[1] else "out of band"
        print(f"  {LAB[n]} {Rp[n]:5.2f} {band:12s} pit lead "
              f"{'n/a' if s150 is None else f'{s150[chr(108)+chr(101)+chr(97)+chr(100)]:5.0f} mm'} @150"
              f"{'' if s180 is None else f', {s180[chr(108)+chr(101)+chr(97)+chr(100)]:5.0f} mm @180'}")

    # ---------- (b) volume ----------
    print(f"\n--- (b) volume at the PRIMARY target (descent {PRIMARY} mm), BOTH bases ---")
    for n in CELLS:
        s = S[PRIMARY][n]
        if s is None:
            print(f"  {LAB[n]} n/a"); continue
        cand = (0.75 <= s["x_pit"] <= 1.25) or (0.75 <= s["x_desc"] <= 1.25)
        print(f"  {LAB[n]} {s['vol']:6.3f} L   pit {s['pit']:5.1f} mm -> {s['x_pit']:5.2f}x   "
              f"descent basis -> {s['x_desc']:5.2f}x   {'candidate' if cand else '-'}")

    # ---------- (d) rings ----------
    print("\n--- (d) REQUIRED WITNESS: ring-resolved recession 30-50 mm, 150-400 s "
          "[mm] (never averaged) ---")
    RG = {n: rings(n) for n in CELLS}
    print("  ring    " + "".join(f"{n:>8s}" for n in CELLS))
    for i, lo in enumerate(RINGS):
        print(f"  {lo:2d}-{lo+2:2d}   " + "".join(f"{RG[n][i]:8.1f}" for n in CELLS))
    print("  vs M0 % " + "".join(f"{'':>8s}" for n in CELLS))
    for i, lo in enumerate(RINGS):
        print(f"  {lo:2d}-{lo+2:2d}   " + "".join(
            f"{100*(RG[n][i]/RG['M0'][i]-1):+8.1f}" if RG["M0"][i] > 0 else "     n/a"
            for n in CELLS))

    # ---------- (e) power ----------
    print("\n--- (e) REQUIRED WITNESS: power split, mean 150-400 s [W], and the weight ---")
    print("  cell                   tops    sides  s_loss    TOTAL  identity   side_faces  w_mean")
    for n in CELLS:
        tp, sf = pmean(n, "patch_P_robin"), pmean(n, "side_P_face")
        sl, tot = pmean(n, "side_P_loss"), pmean(n, "jet_P_face")
        nf, wm = pmean(n, "side_faces"), pmean(n, "side_w_mean")
        f = lambda x, w=8: "  absent" if x is None else f"{x:{w}.1f}"
        idn = None if tot is None else tot - (tp or 0) - (sf or 0) - (sl or 0)
        print(f"  {LAB[n]}{f(tp)}{f(sf)}{f(sl)}{f(tot)}"
              f"{'' if idn is None else f'{idn:+9.3f}':>10s}"
              f"{f(nf,11)}{'   absent' if wm is None else f'{wm:8.4f}'}")

    # ---------- witnesses ----------
    print("\n--- inertness witnesses ---")
    for n in ("C", "D", "E"):
        ab = thermo(n, "side_P_face_above")
        le = thermo(n, "ledger_err")
        print(f"  {n}: side_P_face_above "
              + ("absent" if ab is None else f"max {ab[1].max():+.3e} min {ab[1].min():+.3e} W")
              + f";  max|ledger_err| {np.abs(le[1]).max():.2e}")
    nfC, nfM1 = pmean("C", "side_faces"), pmean("M1", "side_faces")
    spD, spM1 = pmean("D", "side_P_face"), pmean("M1", "side_P_face")
    print(f"  r_min reaches the block: side_faces C {nfC:.0f} vs M1 {nfM1:.0f} "
          f"({100*(nfC/nfM1-1):+.0f} %, bar -30 %)")
    print(f"  area_exact reaches:      side_P_face D {spD:.1f} vs M1 {spM1:.1f} W "
          f"({spD/spM1:.2f}x, bar 2x)")

    # ---------- the four bars ----------
    print("\n=== (a) THE FOUR BARS — A is the control that already passes ===")
    aA = S[PRIMARY]["A"]["a"]
    print(f"  reference A: (a) {aA:.2f} mm @150, Ø 84.4/82.7 @100/125 (descent 180), "
          f"rate {Rp['A']:.3f}, lead {S[PRIMARY]['A']['lead']:.0f} mm, "
          f"vol {S[PRIMARY]['A']['x_desc']:.2f}x")
    for n in ("C", "D", "E"):
        s180, s150 = S[180][n], S[PRIMARY][n]
        if s180 is None:
            b1 = None
        else:
            b1 = all(abs(s180["D"][d] - REF_D) <= 2.0 for d in (100, 125))
        b2 = s150 is not None and s150["a"] <= aA
        b3 = BAND[0] <= Rp[n] <= BAND[1]
        b4 = (s150 is not None and s150["lead"] <= S[PRIMARY]["A"]["lead"]
              and abs(s150["x_desc"] - 1.0) < abs(S[PRIMARY]["A"]["x_desc"] - 1.0))
        d100 = "n/a" if s180 is None else f"{s180['D'][100]:.1f}"
        d125 = "n/a" if s180 is None else f"{s180['D'][125]:.1f}"
        print(f"  {n}: 1) Ø@100/125 = {d100}/{d125} vs 87+-2 -> "
              f"{'n/a (never reached 180 mm)' if b1 is None else ('PASS' if b1 else 'FAIL')}"
              f"\n     2) (a) {s150['a']:.2f} <= {aA:.2f} -> {'PASS' if b2 else 'FAIL'}"
              f"   3) rate {Rp[n]:.2f} in 1.3-1.6 -> {'PASS' if b3 else 'FAIL'}"
              f"   4) lead/vol -> {'PASS' if b4 else 'FAIL'}"
              f"\n     ALL FOUR: {'PASS' if (b1 and b2 and b3 and b4) else 'FAIL'}")

    print("\n=== pre-registered predictions ===")
    wD = pmean("D", "side_w_mean")
    print(f"  P1 D rate > 4.0 (refuted < 3.5): {Rp['D']:.2f} -> "
          f"{'HELD' if Rp['D'] > 4.0 else ('REFUTED' if Rp['D'] < 3.5 else 'neither')}")
    print(f"  P2 side_w_mean on D in 0.55-0.75 (refuted outside 0.45-0.85): "
          f"{'absent' if wD is None else f'{wD:.4f}'} -> "
          + ("n/a" if wD is None else
             ('HELD' if 0.55 <= wD <= 0.75 else
              ('REFUTED' if (wD < 0.45 or wD > 0.85) else 'neither')))
          + "   <- THE DERIVED NUMBER THAT REPLACES f = 0.2")
    i42 = RINGS.index(42); i44 = RINGS.index(44)
    d42 = 100 * (RG["C"][i42] / RG["M0"][i42] - 1); d44 = 100 * (RG["C"][i44] / RG["M0"][i44] - 1)
    print(f"  P3 C restores 42-46 mm within 20 % of M0 (refuted below -50 %): "
          f"{d42:+.1f} / {d44:+.1f} % -> "
          f"{'HELD' if max(abs(d42),abs(d44)) <= 20 else ('REFUTED' if min(d42,d44) < -50 else 'partial')}")
    s180C = S[180]["C"]
    if s180C is None:
        print("  P4 C wall band: n/a (C never reached 180 mm descent)")
    else:
        g100, g125 = s180C["D"][100] - 84.4, s180C["D"][125] - 82.7
        print(f"  P4 C closes the wall band: Ø {s180C['D'][100]:.1f}/{s180C['D'][125]:.1f} "
              f"(gain {g100:+.1f}/{g125:+.1f} mm vs A) -> "
              f"{'HELD' if abs(s180C['D'][100]-REF_D)<=2 and abs(s180C['D'][125]-REF_D)<=2 else ('REFUTED' if g100 < 1 or g125 < 1 else 'partial')}")
    print(f"  P5 C rate in 1.3-1.9 (refuted > 2.4): {Rp['C']:.2f} -> "
          f"{'HELD' if 1.3 <= Rp['C'] <= 1.9 else ('REFUTED' if Rp['C'] > 2.4 else 'neither')}")
    print(f"  P6 E between C and D on descent: C {Rp['C']:.2f} | E {Rp['E']:.2f} | "
          f"D {Rp['D']:.2f} -> {'HELD' if min(Rp['C'],Rp['D']) <= Rp['E'] <= max(Rp['C'],Rp['D']) else 'REFUTED'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
