#!/usr/bin/env python3
"""D2z-0 scoring against CRITERION.md (sha256 ce81f371962f6588..., hashed
read-only before the first scored leg).

Metric definitions are IMPORTED UNCHANGED from studies/d2j0b_plane/analyze.py --
the same estimator the frozen A_f02 baseline was scored with, so the comparison
is like for like. Frozen studies are read by path and never re-run.
"""
import csv
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
R_GATE = 0.025          # the gate mask; NOT 0.030 (that disc has the rim stall)


def load(name, path):
    sp = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


A = load("d2j0b_analyze", os.path.join(STUD, "d2j0b_plane", "analyze.py"))
BLOCKS = A.BLOCKS[:3]
BKEY = [f"{int(lo)}-{int(hi)}" for lo, hi in BLOCKS]


def rates(outdir, name, rmax=R_GATE):
    m = A.R.J.meta_of(os.path.join(outdir, name))
    nxy, dz = m["nxy"], m["dz"]
    c = (np.arange(nxy) + 0.5) * dz
    X, Y = np.meshgrid(c, c, indexing="ij")
    r = np.hypot(X, Y)
    a = np.loadtxt(os.path.join(outdir, name + "_removal_events.csv"),
                   delimiter=",", skiprows=1, ndmin=2)
    return {f"{int(lo)}-{int(hi)}": float(A.col_rate(r, a, (lo, hi))[r < rmax].mean())
            for lo, hi in BLOCKS}


def thermo(outdir, name):
    f = os.path.join(outdir, name, "thermo.dat")
    n = open(f).readline().split()
    rows = np.loadtxt(f, skiprows=1, ndmin=2)
    return {x: rows[:, i] for i, x in enumerate(n)}


def firing_rows(outdir, name):
    p = os.path.join(outdir, name + "_removal_events.csv")
    out = []
    with open(p) as fh:
        for d in csv.DictReader(fh):
            if int(d["regime"]) == 1:
                out.append((float(d["h_applied"]), int(d["n_voided"])))
    return out


def same_bytes(a_dir, a, b_dir, b):
    for rel in ("thermo.dat",):
        pa, pb = os.path.join(a_dir, a, rel), os.path.join(b_dir, b, rel)
        if not filecmp.cmp(pa, pb, shallow=False):
            return False, rel
    pa = os.path.join(a_dir, a + "_removal_events.csv")
    pb = os.path.join(b_dir, b + "_removal_events.csv")
    if not filecmp.cmp(pa, pb, shallow=False):
        return False, "_removal_events.csv"
    ga = sorted(glob.glob(os.path.join(a_dir, a, "*cell", "Level_0", "Cell_D_*")))
    gb = sorted(glob.glob(os.path.join(b_dir, b, "*cell", "Level_0", "Cell_D_*")))
    if not ga or len(ga) != len(gb):
        return False, f"plotfile count {len(ga)} vs {len(gb)}"
    for x, y in zip(ga, gb):
        if not filecmp.cmp(x, y, shallow=False):
            return False, os.path.basename(os.path.dirname(os.path.dirname(x)))
    return True, f"{len(ga)} Cell_D + thermo.dat + removal log"


h = hashlib.sha256(open(os.path.join(HERE, "CRITERION.md"), "rb").read()).hexdigest()
print("=" * 78)
print("D2z-0 -- retire the half-cell removal increment")
print(f"criteria CRITERION.md sha256 {h[:16]}...")
print("metric imported UNCHANGED from d2j0b_plane/analyze.py; gate mask r < 25 mm")
print("=" * 78)

have_f1 = os.path.exists(os.path.join(OUT, "F1.done"))

# ---------------------------------------------------------------- (d) identity
print("\n--- (d) IDENTITY: H2 (key absent) vs frozen d2l_scan/A_f02 ---")
ok, what = same_bytes(OUT, "H2", FROZ, "A_f02")
print(f"  {'PASS' if ok else 'FAIL'}  byte-for-byte: {what}")
print("  tree identity 2639/2639 with the key unset: see PREFLIGHT.md "
      "(witness/d2z0.md5 vs witness/tree0922.md5, guarded)")

# ------------------------------------------------------------------- (a) rate
print(f"\n--- (a) SCORED: disk-mean rate, r < 25 mm [m/h], and the change ---")
froz = {"2mm": rates(FROZ, "A_f02"), "1mm": rates(FROZ, "A_f02_1mm")}
mine = {"2mm": rates(OUT, "F2")}
if have_f1:
    mine["1mm"] = rates(OUT, "F1")
print(f"{'mesh':6s} {'block':10s} {'frozen f=0.2':>13s} {'full_cell':>11s} {'change':>9s}  bar +0.2..+0.6 %")
a_ok = True
for mesh in ("2mm", "1mm"):
    if mesh not in mine:
        print(f"{mesh:6s} (not run)")
        continue
    for b in BKEY:
        d = 100.0 * (mine[mesh][b] / froz[mesh][b] - 1.0)
        good = 0.2 <= d <= 0.6
        a_ok &= good
        print(f"{mesh:6s} {b:10s} {froz[mesh][b]:13.4f} {mine[mesh][b]:11.4f} "
              f"{d:+8.3f}% {'ok' if good else 'OUT'}")
print(f"  (a) {'PASS' if a_ok and have_f1 else ('FAIL' if have_f1 else 'PARTIAL -- F1 not run')}")

print("\n  r < 30 mm alongside (NEVER scored -- includes the 25-30 mm rim stall):")
for mesh, nm in (("2mm", "F2"), ("1mm", "F1")):
    if mesh == "1mm" and not have_f1:
        continue
    f30 = rates(FROZ, "A_f02" if mesh == "2mm" else "A_f02_1mm", 0.030)
    m30 = rates(OUT, nm, 0.030)
    ch = [f"{100*(m30[b]/f30[b]-1.0):+.2f}%" for b in BKEY]
    print(f"    {mesh}: {ch}")

# --------------------------------------------------------------- (b) mesh gap
print("\n--- (b) SCORED: the mesh gap does not move (r < 25 mm) ---")
ref = {b: 100.0 * (froz["1mm"][b] / froz["2mm"][b] - 1.0) for b in BKEY}
if have_f1:
    got = {b: 100.0 * (mine["1mm"][b] / mine["2mm"][b] - 1.0) for b in BKEY}
    b_ok = True
    print(f"{'block':10s} {'f=0.2 ref':>10s} {'full_cell':>10s} {'shift':>8s}  bar < 0.2 pts")
    for b in BKEY:
        s = got[b] - ref[b]
        good = abs(s) < 0.2
        b_ok &= good
        print(f"{b:10s} {ref[b]:+9.2f}% {got[b]:+9.2f}% {s:+7.3f}  {'ok' if good else 'OUT'}")
    print(f"  (b) {'PASS' if b_ok else 'FAIL'}")
else:
    print(f"  reference gaps {[f'{ref[b]:+.2f}%' for b in BKEY]}; F1 not run -> (b) not scored")

# --------------------------------------------------------------- (c) plumbing
print("\n--- (c) SCORED: plumbing ---")
print(f"{'leg':6s} {'firings':>9s} {'voids':>7s} {'firings/void':>13s} "
      f"{'n_voided set':>14s} {'h/dz':>10s} {'thermo fpv':>11s} {'thermo ric':>11s}")
c_ok = True
legs = [("A_f02", FROZ, "frozen"), ("F2", OUT, "full")]
if have_f1:
    legs += [("A_f02_1mm", FROZ, "frozen"), ("F1", OUT, "full")]
for nm, d, kind in legs:
    rows = firing_rows(d, nm)
    dz = A.R.J.meta_of(os.path.join(d, nm))["dz"]
    v = sum(n for _, n in rows)
    fpv = len(rows) / v if v else float("nan")
    nvs = sorted({n for _, n in rows})
    hh = np.array([x for x, _ in rows]) / dz
    T = thermo(d, nm)
    if "firings_per_void" in T:
        x = T["firings_per_void"][T["firings_per_void"] > 0]
        y = T["removal_increment_cells"][T["removal_increment_cells"] > 0]
        tf = f"{x.max():.5f}" if x.size else "-"
        tr = f"{y.min():.5f}" if y.size else "-"
    else:
        tf = tr = "absent"
    print(f"{nm:6s} {len(rows):9d} {v:7d} {fpv:13.4f} {str(nvs):>14s} "
          f"{hh.mean():10.5f} {tf:>11s} {tr:>11s}")
    if kind == "full":
        good = (abs(fpv - 1.0) < 1e-12 and nvs == [1] and np.all(hh == 1.0)
                and tf != "absent" and abs(float(tf) - 1.0) < 1e-12)
        c_ok &= good
print(f"  (c) {'PASS' if c_ok and have_f1 else ('FAIL' if not c_ok else 'PARTIAL -- F1 not run')}")
print("  NOTE: ACTIVE_STEP states the baseline as 1.977; that is MEIER A's value.")
print("        On the arbiter this packet runs it is 1.957 (2 mm) / 1.962 (1 mm).")
print("        The bar is the destination (exactly 1.000), which is unaffected.")

# -------------------------------------------------------------- predictions
print("\n=== pre-registered predictions ===")
if have_f1:
    ch = [100.0 * (mine[m][b] / froz[m][b] - 1.0) for m in ("2mm", "1mm") for b in BKEY]
    print(f"  P1 rate change in +0.2..+0.6 % at both meshes: "
          f"[{min(ch):+.3f}, {max(ch):+.3f}] -> {'HELD' if all(0.2 <= c <= 0.6 for c in ch) else 'REFUTED'}")
    sh = [abs(got[b] - ref[b]) for b in BKEY]
    print(f"  P2 mesh-gap shift < 0.2 pts: max {max(sh):.3f} -> "
          f"{'HELD' if max(sh) < 0.2 else 'REFUTED'}")
print(f"  P3 firings_per_void -> exactly 1.000: {'HELD' if c_ok else 'REFUTED'}")
print(f"  P4 H2 byte-identical to frozen A_f02: {'HELD' if ok else 'REFUTED'}")
