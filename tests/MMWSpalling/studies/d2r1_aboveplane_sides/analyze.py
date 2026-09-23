#!/usr/bin/env python3
"""D2r-1 scoring against CRITERION.md (hashed 33b27e47..., read-only, before
either 1 mm leg finished).

Metric definitions are IMPORTED from studies/d2j0b_plane/analyze.py and used
UNCHANGED -- the same estimator as the frozen A_f02 baseline, so the comparison
is like for like. The frozen baseline is read by path and never re-run.
"""
import importlib.util
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)


def load(name, path):
    sp = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


A = load("d2j0b_analyze", os.path.join(STUD, "d2j0b_plane", "analyze.py"))
BLOCKS = A.BLOCKS[:3]
OUT = os.path.join(HERE, "output")
FROZ = os.path.join(STUD, "d2l_scan", "output")


def rates(outdir, name):
    m = A.R.J.meta_of(os.path.join(outdir, name))
    nxy, dz = m["nxy"], m["dz"]
    c = (np.arange(nxy) + 0.5) * dz
    X, Y = np.meshgrid(c, c, indexing="ij")
    r = np.hypot(X, Y)
    a = np.loadtxt(os.path.join(outdir, name + "_removal_events.csv"),
                   delimiter=",", skiprows=1, ndmin=2)
    return {f"{int(lo)}-{int(hi)}": float(A.col_rate(r, a, (lo, hi))[r < 0.030].mean())
            for lo, hi in BLOCKS}


def thermo(name):
    f = os.path.join(OUT, name, "thermo.dat")
    n = open(f).readline().split()
    rows = np.loadtxt(f, skiprows=1, ndmin=2)
    k = {x: i for i, x in enumerate(n)}
    return {x: rows[:, k[x]] for x in n}


print("=" * 78)
print("D2r-1 -- above-plane SIDE faces: does the mesh gap close?")
print("=" * 78)

CASES = ["B_2mm", "B_1mm", "T_2mm", "T_1mm", "S_2mm"]
if os.path.exists(os.path.join(OUT, "S_1mm.done")):
    CASES.append("S_1mm")
R = {c: rates(OUT, c) for c in CASES}
F = {"A_f02": rates(FROZ, "A_f02"), "A_f02_1mm": rates(FROZ, "A_f02_1mm")}

print("\n## Disk-mean rate [m/h], r < 30 mm")
print(f"{'block':10s} {'B_2mm':>8s} {'B_1mm':>8s} {'gap':>8s} | "
      f"{'T_2mm':>8s} {'T_1mm':>8s} {'gap':>8s} | {'S_2mm':>8s}")
gb, gt = {}, {}
for b in R["B_2mm"]:
    gb[b] = R["B_1mm"][b] / R["B_2mm"][b] - 1.0
    gt[b] = R["T_1mm"][b] / R["T_2mm"][b] - 1.0
    print(f"{b:10s} {R['B_2mm'][b]:8.4f} {R['B_1mm'][b]:8.4f} {100*gb[b]:+7.1f}% | "
          f"{R['T_2mm'][b]:8.4f} {R['T_1mm'][b]:8.4f} {100*gt[b]:+7.1f}% | {R['S_2mm'][b]:8.4f}")

print("\n## Goal item 1 -- baseline validity (B vs the FROZEN A_f02)")
for mine, froz in (("B_2mm", "A_f02"), ("B_1mm", "A_f02_1mm")):
    d = max(abs(R[mine][b] / F[froz][b] - 1.0) for b in R[mine])
    print(f"   {mine} vs {froz}: max relative difference {100*d:.3e} %"
          f"   {'IDENTICAL' if d == 0.0 else 'DIFFERS'}")

print("\n## (a) THE SCORED ONE -- mesh gap <= 5 % in all three blocks")
a_ok = True
for b in gt:
    ok = abs(gt[b]) <= 0.05
    a_ok &= ok
    print(f"   {b:10s} baseline {100*gb[b]:+6.1f}%  ->  test {100*gt[b]:+6.1f}%   "
          f"({'PASS' if ok else 'FAIL'})   change {100*(abs(gt[b])-abs(gb[b])):+.1f} pts")
print(f"   => (a) {'PASS' if a_ok else 'FAIL'}")

print("\n## (b) THE SEPARATOR -- S must reproduce B (tops alone are inert)")
print("   S carries above-plane TOPS at h = 142 but forces the above-plane SIDE")
print("   factor to 0. So S - B isolates the tops, and T - S isolates the sides.")
for mesh in ("2mm", "1mm"):
    if f"S_{mesh}" not in R:
        print(f"   S_{mesh}: NOT RUN")
        continue
    for b in R[f"S_{mesh}"]:
        sb = R[f"S_{mesh}"][b] / R[f"B_{mesh}"][b] - 1.0
        ts = R[f"T_{mesh}"][b] / R[f"S_{mesh}"][b] - 1.0
        print(f"   S_{mesh} {b:10s} S/B-1 = {100*sb:+6.2f} %   (tops)   |   "
              f"T/S-1 = {100*ts:+6.2f} %   (SIDES)")
    if f"S_{mesh}" in R:
        g = {b: R[f"S_{mesh}"][b] for b in R[f"S_{mesh}"]}
        if mesh == "1mm":
            print("   S mesh gap (S_1mm vs S_2mm):")
            for b in g:
                print(f"      {b:10s} {100*(R['S_1mm'][b]/R['S_2mm'][b]-1):+6.1f} %"
                      f"   (baseline {100*gb[b]:+.1f} %, test {100*gt[b]:+.1f} %)")

print("\n## (c) THE WITNESS -- side_P_face_above > 0 on T, = 0 on S")
for c in [x for x in ("T_2mm", "T_1mm", "S_2mm", "S_1mm") if x in R]:
    t = thermo(c)
    ab = t.get("side_P_face_above", np.zeros(1))
    nc = t.get("side_cols_above", np.zeros(1))
    print(f"   {c:8s} max {ab.max():8.2f} W, mean {ab.mean():8.2f} W, faces max {nc.max():5.0f}"
          f"   -> {'ACTIVE' if ab.max() > 0 else 'zero (as intended)'}")

print("\n## (d) THE DRAIN -- wall_step_q vs D2o-0b's 302 kW/m2 (58 % of q_pin)")
for c in [x for x in ("T_2mm", "T_1mm", "S_2mm", "S_1mm") if x in R]:
    t = thermo(c)
    q, nf = t["wall_step_q"], t["wall_step_n"]
    print(f"   {c:8s} final {q[-1]/1e3:7.1f} kW/m2 ({100*q[-1]/520e3:4.0f} % of q_pin) "
          f"over {nf[-1]:5.0f} faces;  max {q.max()/1e3:7.1f}")
print("   (B carries no wall_step_q: the column is gated on side_face_factor_above so that")
print("    B stays byte-identical to the frozen A_f02 -- see CRITERION.md section 3.)")

print("\n## (e) THE CEILING -- T vs B at the same mesh, bound +15 %")
e_ok = True
for mesh in ("2mm", "1mm"):
    for b in R[f"T_{mesh}"]:
        ex = R[f"T_{mesh}"][b] / R[f"B_{mesh}"][b] - 1.0
        ok = ex <= 0.15
        e_ok &= ok
        print(f"   T_{mesh} {b:10s} {100*ex:+6.2f} %   {'PASS' if ok else 'FAIL'}")
print(f"   => (e) {'PASS' if e_ok else 'FAIL'}")

print("\n## ledger closure")
for c in CASES:
    print(f"   {c:8s} |ledger_err| max {np.abs(thermo(c)['ledger_err']).max():.2e}")

print("\n" + "=" * 78)
print(f"(a) {'PASS' if a_ok else 'FAIL'}    (e) {'PASS' if e_ok else 'FAIL'}")
print("=" * 78)
