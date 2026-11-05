#!/usr/bin/env python3
"""Reduce the what-if legs with P2-7's accepted per-leg reduction (reduce_p2_7.one), plus rate-vs-time and exponents."""
import math, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
P27 = os.path.join(HERE, "..", "p2_7_hoeser_heat")
sys.path.insert(0, P27); sys.path.insert(0, HERE)
import reduce_p2_7 as R
from run_whatif import SCALES
SC = dict(HG=1.0, **SCALES)
R.legdir = lambda n: os.path.join(P27, "legs", "HG") if n == "HG" else os.path.join(HERE, "legs", n)
ORDER = ["HG", "s0581", "s150", "s087", "s100", "s058", "s034"]
HER_V, HER_EMAX = 12.9, (41.54, 47.39)
res = {}
print("What-if: P2-7 HG with the heat scaled by S (README.md). Quarter domain x4 on volume and energy.")
print("=" * 118)
print(f"{'leg':6s} {'S':>7s} | {'V cm3':>8s} {'/her':>6s} | {'E_abs kJ':>8s} {'flake%':>6s} | {'fireT':>6s} | {'axis mm':>7s} {'r17.5':>6s} {'ratio':>6s} | {'onset s':>7s} | rate by 5 s window [cm3/min]")
for n in ORDER:
    d = R.legdir(n)
    if not os.path.exists(os.path.join(d, ".done")) and n != "HG":
        print(f"{n:6s} not done"); continue
    r = R.one(n); ev = R.events(d)
    win = [0.0] * 6; t_on = None
    for e in ev:
        t = R._f(e, "time"); nv = R._i(e, "n_voided")
        if nv > 0 and t_on is None: t_on = t
        win[min(5, int(t // 5.0))] += nv * r["g"]["V_cell"] * 4e6 * 12.0
    r["win"], r["t_on"] = win, t_on; res[n] = r
    ratio = r["dr"] if r["axis"] > 0 else float("nan")
    print(f"{n:6s} {SC[n]:7.4f} | {r['V']:8.3f} {r['V']/HER_V:6.2f} | {r['E_src']/1e3:8.2f} {100*r['share']:6.1f} | {r['fire_med']:6.1f} |"
          f" {1e3*r['axis']:7.2f} {1e3*r['d175']:6.2f} {ratio:6.3f} | {('%.2f' % t_on) if t_on is not None else '  none':>7s} | "
          + " ".join(f"{w:6.1f}" for w in win))
print("-" * 118)
print("FLUX EXPONENT n = ln(V_hi/V_lo) / ln(S_hi/S_lo)   (her n_q 1.44-1.83, 47 MPa series)")
for a, b in [("HG", "s0581"), ("s150", "s087"), ("s100", "s058"), ("s058", "s034")]:
    if a in res and b in res and res[b]["V"] > 0:
        print(f"  {a:6s}/{b:6s}  V {res[a]['V']:8.3f} / {res[b]['V']:8.3f}   n = {math.log(res[a]['V']/res[b]['V'])/math.log(SC[a]/SC[b]):.2f}"
              f"   (on absorbed energy: {math.log(res[a]['V']/res[b]['V'])/math.log(res[a]['E_src']/res[b]['E_src']):.2f})")
    elif a in res and b in res:
        print(f"  {a:6s}/{b:6s}  low leg removed nothing")
print(f"her 0 MPa: V {HER_V} cm3; absorbable at most {HER_EMAX[0]}-{HER_EMAX[1]} kJ (P2-7 BOUND, 850 K); her 47 MPa bowl ratio 0.61")
