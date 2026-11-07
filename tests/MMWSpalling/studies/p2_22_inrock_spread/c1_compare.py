#!/usr/bin/env python3
"""P2-22: C1 (the Zhang deck, one-pair table, eps_r 5.02) against P2-20's legs/Z1/run. EXPECTED TO MOVE (the in-rock
spread changes every eps_r > 1 run); this records the byte verdict AND the magnitude from thermo.dat (column-wise
max relative difference over the common rows). Usage: c1_compare.py > C1_vs_Z1.txt
"""
import filecmp, math, os
HERE = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(HERE, "..", "p2_20_beam_deposition", "legs", "Z1", "run")
NEW = os.path.join(HERE, "legs", "C1", "run")

def table(p):
    L = open(p).read().splitlines()
    h = [c.strip() for c in L[0].split("\t")]
    rows = [[float(x) for x in r.split("\t")] for r in L[1:] if r.strip()]
    return h, rows

pfs = sorted(p for p in os.listdir(REF) if p.endswith("cell"))
n, bad = 0, 0
for pf in pfs:
    for dp, _ds, fs in os.walk(os.path.join(REF, pf)):
        for fn in fs:
            a = os.path.join(dp, fn); b = os.path.join(NEW, os.path.relpath(a, REF)); n += 1
            if not (os.path.exists(b) and filecmp.cmp(a, b, shallow=False)): bad += 1
ha, ra = table(os.path.join(REF, "thermo.dat")); hb, rb = table(os.path.join(NEW, "thermo.dat"))
print("C1 vs p2_20 Z1: %d plotfile files, %d differ (byte identity EXPECTED to fail: eps_r 5.02 => in-rock spread)" % (n, bad))
m = min(len(ra), len(rb))
print("thermo.dat %d common rows; per-column max relative difference (skipping timers):" % m)
for c in ha:
    if c not in hb or c in ("time", "step") or c.endswith(("_t_path", "_t_removal", "_t_step")) or c.startswith("fl_t"):
        continue
    ia, ib = ha.index(c), hb.index(c); worst = 0.0; at = None
    for k in range(m):
        x, y = ra[k][ia], rb[k][ib]
        if not (math.isfinite(x) and math.isfinite(y)): continue
        d = abs(y - x) / max(abs(x), 1e-300)
        if d > worst: worst, at = d, (ra[k][ha.index("time")], x, y)
    if worst > 1e-12:
        print("   %-24s %9.3e  (t %g s: %g -> %g)" % (c, worst, *at))
