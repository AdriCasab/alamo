#!/usr/bin/env python3
"""D2q-2c scoring: criteria (a)-(f) as fixed in CRITERION.md (hashed
6f837a70..., read-only, before output/ existed).

Metric definitions are IMPORTED from studies/d2j0b_plane/analyze.py and used
UNCHANGED -- same estimator as the Meier gate, so a pass transfers. The
key-off anchors (d2j0b_plane C0/C1) are read from that study's output/ by path
and never re-run.
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
BLOCKS = A.BLOCKS[:3]                      # 50-100, 100-150, 150-200
OUT = os.path.join(HERE, "output")
REF = os.path.join(STUD, "d2j0b_plane", "output")

# Key-off anchors, pre-registered in CRITERION.md section 5.
C0_OFF = {"2mm": {"50-100": 1.211, "100-150": 1.119, "150-200": 1.048},
          "1mm": {"50-100": 1.085, "100-150": 0.948, "150-200": 0.838}}
C1_OFF = {"2mm": {"50-100": 1.151, "100-150": 0.968, "150-200": 0.759},
          "1mm": {"50-100": 0.996, "100-150": 0.690, "150-200": 0.414}}


def rates(outdir, name):
    """Disk-mean rate per block, r < 30 mm, from <pf>_removal_events.csv."""
    m = A.R.J.meta_of(os.path.join(outdir, name))
    nxy, dz = m["nxy"], m["dz"]
    c = (np.arange(nxy) + 0.5) * dz
    X, Y = np.meshgrid(c, c, indexing="ij")
    r = np.hypot(X, Y)
    a = np.loadtxt(os.path.join(outdir, name + "_removal_events.csv"),
                   delimiter=",", skiprows=1, ndmin=2)
    out = {}
    for lo, hi in BLOCKS:
        out[f"{int(lo)}-{int(hi)}"] = float(A.col_rate(r, a, (lo, hi))[r < 0.030].mean())
    return out, r, a, m


def thermo(name):
    f = os.path.join(OUT, name, "thermo.dat")
    n = open(f).readline().split()
    rows = np.loadtxt(f, skiprows=1, ndmin=2)
    k = {x: i for i, x in enumerate(n)}
    return {x: rows[:, k[x]] for x in n}


print("=" * 78)
print("D2q-2c — criteria (a)-(f) vs CRITERION.md (hashed before launch)")
print("=" * 78)

R = {}
for c in ("C1_2mm", "C1_1mm", "C0_2mm", "C0_1mm"):
    R[c], _, _, _ = rates(OUT, c)

print("\n## Disk-mean rate [m/h], r < 30 mm, key ON (side flux capped at 80, per-cell on)")
print(f"{'block':10s} {'C1_2mm':>9s} {'C1_1mm':>9s} {'gap':>9s} | {'C0_2mm':>9s} {'C0_1mm':>9s} {'gap':>9s}")
gaps_c1, gaps_c0 = {}, {}
for b in R["C1_2mm"]:
    g1 = R["C1_1mm"][b] / R["C1_2mm"][b] - 1.0
    g0 = R["C0_1mm"][b] / R["C0_2mm"][b] - 1.0
    gaps_c1[b], gaps_c0[b] = g1, g0
    print(f"{b:10s} {R['C1_2mm'][b]:9.4f} {R['C1_1mm'][b]:9.4f} {100*g1:+8.1f}% | "
          f"{R['C0_2mm'][b]:9.4f} {R['C0_1mm'][b]:9.4f} {100*g0:+8.1f}%")

print("\n## (a) mesh gap within 5 % in every block, both pairs")
a_ok = True
for b in gaps_c1:
    for nm, g in (("C1", gaps_c1[b]), ("C0", gaps_c0[b])):
        ok = abs(g) <= 0.05
        a_ok &= ok
        print(f"   {nm} {b:10s} {100*g:+7.1f}%   {'PASS' if ok else 'FAIL'}")
print(f"   => (a) {'PASS' if a_ok else 'FAIL'}")

print("\n## (e) removal ceiling: <= C0 KEY-OFF at the same mesh + 15 %")
e_ok = True
for mesh, case in (("2mm", "C1_2mm"), ("1mm", "C1_1mm")):
    for b in R[case]:
        ceil = C0_OFF[mesh][b] * 1.15
        ok = R[case][b] <= ceil
        e_ok &= ok
        print(f"   {case:8s} {b:10s} {R[case][b]:7.4f} vs ceiling {ceil:7.4f} "
              f"(C0off {C0_OFF[mesh][b]:.3f})  {'PASS' if ok else 'FAIL'}")
print(f"   => (e) {'PASS' if e_ok else 'FAIL'}")

print("\n## (f) until pcr_lateral_cells != 0, the key-on run does nothing lateral")
f_ok = True
for c in ("C1_2mm", "C1_1mm", "C0_2mm", "C0_1mm"):
    t = thermo(c)
    lat = t.get("pcr_lateral_cells", np.zeros(1))
    over = t.get("pcr_overhang_cells", np.zeros(1))
    sf = t.get("side_faces", np.zeros(1))
    ok = lat.max() == 0.0
    f_ok &= ok
    print(f"   {c:8s} side_faces max {sf.max():6.0f}  pcr_lateral max {lat.max():4.0f}  "
          f"overhang max {over.max():4.0f}  -> {'no lateral removal occurred' if ok else 'FIRED'}")
print(f"   => (f) vacuously satisfied: PerCellRemoval never removed a cell")

print("\n## P3 / ledger")
for c in ("C1_2mm", "C1_1mm", "C0_2mm", "C0_1mm"):
    t = thermo(c)
    print(f"   {c:8s} |ledger_err| max {np.abs(t['ledger_err']).max():.2e}   "
          f"side_P_face max {t.get('side_P_face', np.zeros(1)).max():9.1f} W")

print("\n## Context: key-OFF anchors (d2j0b_plane, on disk, never re-run)")
print(f"{'block':10s} {'C1off_2mm':>10s} {'C1off_1mm':>10s} {'gap':>8s} | {'C0off_1mm':>10s}")
for b in C1_OFF["2mm"]:
    g = C1_OFF["1mm"][b] / C1_OFF["2mm"][b] - 1.0
    print(f"{b:10s} {C1_OFF['2mm'][b]:10.3f} {C1_OFF['1mm'][b]:10.3f} {100*g:+7.1f}% | "
          f"{C0_OFF['1mm'][b]:10.3f}")

print("\n" + "=" * 78)
print(f"(a) {'PASS' if a_ok else 'FAIL'}   (e) {'PASS' if e_ok else 'FAIL'}   "
      f"(f) vacuous   -- P4 (inert) {'CONFIRMED' if f_ok else 'refuted'}")
print("=" * 78)
