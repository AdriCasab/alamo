#!/usr/bin/env python3
"""D2z-1 scoring against CRITERION.md (sha256 12db5609fec1652c..., hashed
read-only with zero outputs on disk).

Metric definitions are IMPORTED UNCHANGED from studies/d2j0b_plane/analyze.py --
the same estimator the frozen A_f02 baseline was scored with. Frozen studies are
read by path and never re-run.
"""
import filecmp
import glob
import hashlib
import importlib.util
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)
OUT = os.path.join(HERE, "output")
FROZ = os.path.join(STUD, "d2l_scan", "output")
Z0 = os.path.join(STUD, "d2z0_increment", "output")
R_GATE = 0.025          # the gate mask; NOT 0.030 (that disc has the rim stall)
WIN = (50.0, 200.0)     # SINGLE pre-registered window (D2z-0 takeaway 5)


def load(name, path):
    sp = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


A = load("d2j0b_analyze", os.path.join(STUD, "d2j0b_plane", "analyze.py"))


def geom(d, n):
    m = A.R.J.meta_of(os.path.join(d, n))
    nxy, dz = m["nxy"], m["dz"]
    c = (np.arange(nxy) + 0.5) * dz
    X, Y = np.meshgrid(c, c, indexing="ij")
    return np.hypot(X, Y), dz


def ev(d, n):
    return np.loadtxt(os.path.join(d, n + "_removal_events.csv"),
                      delimiter=",", skiprows=1, ndmin=2)


def cells_voided(d, n):
    """Cells voided inside the gate disc over the single window."""
    r, _ = geom(d, n)
    a = ev(d, n)
    i, j = a[:, 1].astype(int), a[:, 2].astype(int)
    m = (r[i, j] < R_GATE) & (a[:, 10] > 0) & (a[:, 0] >= WIN[0]) & (a[:, 0] < WIN[1])
    return float(a[m, 10].sum())


def agg_rate(d, n):
    r, _ = geom(d, n)
    return float(A.col_rate(r, ev(d, n), WIN)[r < R_GATE].mean())


def thermo(d, n):
    f = os.path.join(d, n, "thermo.dat")
    h = open(f).readline().split()
    a = np.loadtxt(f, skiprows=1, ndmin=2)
    return h, {c: a[:, i] for i, c in enumerate(h)}


h = hashlib.sha256(open(os.path.join(HERE, "CRITERION.md"), "rb").read()).hexdigest()
print("=" * 78)
print("D2z-1 -- whole-cell removal unconditional, spall.removal_increment deleted")
print(f"criteria CRITERION.md sha256 {h[:16]}...")
print("metric imported UNCHANGED from d2j0b_plane/analyze.py; gate r < 25 mm; "
      "SINGLE window 50-200 s")
print("=" * 78)

# --------------------------------------------------- (c0) byte-identity vs D2z-0
print("\n--- (c0) SCORED: byte-identity against D2z-0 (same physics, key removed) ---")
c0_ok = True
for old, new in (("F2", "W2"), ("F1", "W1")):
    same_log = filecmp.cmp(os.path.join(Z0, old + "_removal_events.csv"),
                           os.path.join(OUT, new + "_removal_events.csv"), shallow=False)
    ga = sorted(glob.glob(os.path.join(Z0, old, "*cell", "Level_0", "Cell_D_*")))
    gb = sorted(glob.glob(os.path.join(OUT, new, "*cell", "Level_0", "Cell_D_*")))
    bad = [x for x, y in zip(ga, gb) if not filecmp.cmp(x, y, shallow=False)]
    ha, da = thermo(Z0, old)
    hb, db = thermo(OUT, new)
    shared = [c for c in ha if c in hb]
    numdiff = [c for c in shared if not np.array_equal(da[c], db[c])]
    ok = same_log and ga and len(ga) == len(gb) and not bad and not numdiff
    c0_ok &= bool(ok)
    print(f"  {new} vs {old}: removal log {'IDENTICAL' if same_log else 'DIFFERS'}; "
          f"{len(gb)} Cell_D, {len(bad)} differing; "
          f"thermo shared {len(shared)} cols, {len(numdiff)} differing")
    print(f"      header delta: added {[c for c in hb if c not in ha]}, "
          f"removed {[c for c in ha if c not in hb]}")
print(f"  (c0) {'PASS' if c0_ok else 'FAIL'}")

# ------------------------------------------------------- (c) cells voided
print("\n--- (c) SCORED: cells voided inside r < 25 mm, 50-200 s  [bar +0.2..+1.0 %] ---")
print(f"{'mesh':6s} {'frozen f=0.2':>13s} {'whole-cell':>11s} {'change':>9s}")
c_ok, changes = True, []
for mesh, fn, mn in (("2mm", "A_f02", "W2"), ("1mm", "A_f02_1mm", "W1")):
    a, b = cells_voided(FROZ, fn), cells_voided(OUT, mn)
    ch = 100.0 * (b / a - 1.0)
    changes.append(ch)
    good = 0.2 <= ch <= 1.0
    c_ok &= good
    print(f"{mesh:6s} {a:13.0f} {b:11.0f} {ch:+8.2f}% {'ok' if good else 'OUT'}")
print(f"  (c) {'PASS' if c_ok else 'FAIL'}")

# ------------------------------------------------------- (d) reported only
print("\n--- (d) REPORTED, not gated: aggregate rate and mesh gap (50-200 s) ---")
rates = {}
for mesh, fn, mn in (("2mm", "A_f02", "W2"), ("1mm", "A_f02_1mm", "W1")):
    r0, r1 = agg_rate(FROZ, fn), agg_rate(OUT, mn)
    rates[mesh] = (r0, r1)
    print(f"  {mesh}: {r0:.5f} -> {r1:.5f} m/h  ({100*(r1/r0-1):+.3f} %)   "
          f"[D2z-0 reference {'+0.352' if mesh=='2mm' else '+0.438'} %]")
g0 = 100.0 * (rates["1mm"][0] / rates["2mm"][0] - 1.0)
g1 = 100.0 * (rates["1mm"][1] / rates["2mm"][1] - 1.0)
print(f"  aggregate mesh gap: {g0:+.2f}% -> {g1:+.2f}%   shift {g1-g0:+.3f} pts "
      f"[D2z-0 reference +0.085]")

# ------------------------------------------------------- (b) plumbing witnesses
print("\n--- plumbing witnesses (criterion b, fourth bullet) ---")
# D2z-1 Goal 4 (revised): there is NO firings_per_void column -- under
# whole-cell it is 1.000 by construction, it carried no information, and it
# broke jet_enthalpy's "jet_* are the last eight columns". The invariant is
# enforced by the n_voided == 1 assertion in Removal.H and read back here from
# the removal log, which is the stronger statement anyway: an event IS a
# removal.
p3_ok = True
for d, n in ((OUT, "W2"), (OUT, "W1")):
    hh, tt = thermo(d, n)
    a_ = ev(d, n)
    r1 = a_[a_[:, 3] == 1]
    nv = r1[:, 10]
    rows, voids = r1.shape[0], float(nv.sum())
    ok = rows > 0 and rows == voids and set(np.unique(nv).tolist()) == {1.0}
    p3_ok &= bool(ok)
    print(f"  {n}: {rows} regime-1 rows, {voids:.0f} cells voided, "
          f"ratio {rows/max(voids,1):.4f}; n_voided values "
          f"{sorted(set(np.unique(nv).tolist()))}; "
          f"firings_per_void column absent: {'firings_per_void' not in tt}")

# ------------------------------------------------------- Stage C
print("\n--- Stage C: re-baselined anchors (written under THIS study's output/) ---")
print("  W2 / W1 ARE the re-baselined d2l_scan A_f02 pair (same argv, whole-cell).")
for n in ("RB_C1_2mm", "RB_D"):
    p = os.path.join(OUT, n + ".done")
    print(f"  {n}: {'done -- ' + open(p).read()[:72] + '...' if os.path.exists(p) else 'NOT RUN'}")

print("\n=== pre-registered predictions ===")
print(f"  P2 cells voided +0.2..+1.0 % at both meshes: "
      f"[{min(changes):+.2f}, {max(changes):+.2f}] -> {'HELD' if c_ok else 'REFUTED'}")
print(f"  P3 firings_per_void == 1 and n_voided == 1 everywhere: "
      f"{'HELD' if p3_ok else 'REFUTED'}")
print("  P1 (all tests), P4 (message + witness), P5 (re-pin audit): see RESULTS.md")
