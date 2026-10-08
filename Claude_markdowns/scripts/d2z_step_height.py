#!/usr/bin/env python3
"""D2z step-height estimator — settles WHICH delta the packet's table measured.

READ-ONLY. Reads frozen removal_events.csv only; re-runs nothing, writes
nothing into any study. Arbiter runs are a QUARTER domain with the hot disk at
the ORIGIN corner (x0 = y0 = 0, col_i/col_j 0..14 at 2 mm) -- centring on the
domain middle silently yields zero columns.

Continuous surface height per column, in cells measured DOWNWARD, linear
between that column's own void events (flip-interval CV 0.01-0.11, so the
interpolation is sound). For an adjacent pair with model k_top drop n:

  TRUE  delta = |S_taller - S_shorter|        <- both columns' own phase
  PKT   delta = (n - 1) + phi_lower           <- ACTIVE_STEP D2z Step 1,
                                                 phi_col = clamp((t-t_last)/dt_mean,0,1)

RESULT (2026-09-24, blocks 50-100 / 150-200 s, r < 25 mm):

  leg           notes' delta_birth   TRUE (this)      PKT (this)
  f=0.2  2 mm     0.048 / 0.068      0.046 / 0.065    0.005 / 0.005
  f=0    2 mm     0.171 / 0.411      0.168 / 0.399    0.005 / 0.005
  f=0.2  1 mm     0.056 / 0.089      0.055 / 0.088    0.011 / 0.011

TRUE reproduces the frozen measurement to <= 0.003 on all six; PKT is ~10x low
and is ~0 BY CONSTRUCTION -- a step is born the instant the deeper column
flips, so phi_lower ~ 0 exactly then. PKT measures time since the deeper
column flipped, not the step height.

TIME-WEIGHTED means of TRUE (what an in-code thermo `phi_mean` must read,
since it is a population mean over live faces, not an at-birth sample):

  f=0.2  2 mm   0.132 / 0.156 / 0.169   (blocks 50-100 / 100-150 / 150-200)
  f=0.2  1 mm   0.193 / 0.260 / 0.319

=> 0.13-0.17 on the cured 2 mm floor is essentially the f = 0.2 that D2k/D2l
found empirically, which is WHY f = 0.2 works at 2 mm; at 1 mm the true value
is 0.19-0.32, so a fixed f would need ~0.25-0.3. That is the closure's actual
value-add, and the reason the gate has to be run at 1 mm.
"""
import csv, sys, math
from collections import defaultdict

def load(path, dx, n_cell, r_max, x0=None, y0=None):
    if x0 is None: x0 = 0.5 * dx * n_cell
    if y0 is None: y0 = 0.5 * dx * n_cell
    flips = defaultdict(list)          # (i,j) -> [(t, n_voided), ...]
    with open(path) as f:
        for row in csv.DictReader(f):
            nv = int(row["n_voided"])
            if nv <= 0: continue
            i, j = int(row["col_i"]), int(row["col_j"])
            x = (i + 0.5) * dx; y = (j + 0.5) * dx
            if math.hypot(x - x0, y - y0) > r_max: continue
            flips[(i, j)].append((float(row["time"]), nv))
    for v in flips.values(): v.sort()
    return flips

def build(flips):
    """cols -> (times[], level_after[], dt_mean_running[])  level in cells down."""
    out = {}
    for c, ev in flips.items():
        if len(ev) < 2: continue
        ts = [e[0] for e in ev]
        lv, s = [], 0
        for _, nv in ev:
            s += nv; lv.append(s)
        out[c] = (ts, lv)
    return out

def S_of(ts, lv, t):
    """continuous level (cells down) at time t, linear between flips."""
    if t <= ts[0]:  return 0.0
    if t >= ts[-1]: return float(lv[-1])
    lo, hi = 0, len(ts) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if ts[mid] <= t: lo = mid
        else: hi = mid
    f = (t - ts[lo]) / (ts[hi] - ts[lo])
    return lv[lo] + f * (lv[hi] - lv[lo])

def k_of(ts, lv, t):
    """discrete model level (cells down) at time t."""
    if t < ts[0]: return 0
    lo, hi = 0, len(ts)
    while hi - lo > 0:
        mid = (lo + hi) // 2
        if ts[mid] <= t: lo = mid + 1
        else: hi = mid
    return lv[lo - 1] if lo > 0 else 0

def phi_of(ts, lv, t):
    """ACTIVE_STEP phi_col = clamp((t - t_last)/dt_mean, 0, 1), dt_mean = that
    column's running mean interval (>=2 removals, else fallback 1.0)."""
    if t < ts[0]: return 1.0, True
    lo, hi = 0, len(ts)
    while hi - lo > 0:
        mid = (lo + hi) // 2
        if ts[mid] <= t: lo = mid + 1
        else: hi = mid
    m = lo - 1
    if m < 1: return 1.0, True            # fewer than two removals
    dt_mean = (ts[m] - ts[0]) / m
    if dt_mean <= 0: return 1.0, True
    return min(1.0, max(0.0, (t - ts[m]) / dt_mean)), False

def run(path, dx, n_cell, r_max, blocks, label, sample_dt=0.05, x0=None):
    cols = build(load(path, dx, n_cell, r_max, x0=x0, y0=x0))
    keys = set(cols)
    pairs = []
    for (i, j) in keys:
        for d in ((1, 0), (0, 1)):
            o = (i + d[0], j + d[1])
            if o in keys: pairs.append(((i, j), o))
    print(f"\n=== {label}  ({len(cols)} cols, {len(pairs)} adjacent pairs, r<{r_max*1e3:.0f} mm) ===")
    print(f"{'block':>12} {'step%':>7} {'>=2c%':>7} {'TRUE birth':>11} {'PKT birth':>10} "
          f"{'TRUE t-wtd':>11} {'PKT t-wtd':>10} {'fallback%':>10}")
    for (t0, t1) in blocks:
        tw_true, tw_pkt, b_true, b_pkt = [], [], [], []
        n_samp = n_step = n_multi = n_fb = 0
        prev_n = {}
        t = t0
        while t <= t1:
            for (a, b) in pairs:
                ta, la = cols[a]; tb, lb = cols[b]
                ka, kb = k_of(ta, la, t), k_of(tb, lb, t)
                n = abs(ka - kb)
                n_samp += 1
                if n >= 1:
                    n_step += 1
                    if n >= 2: n_multi += 1
                    Sa, Sb = S_of(ta, la, t), S_of(tb, lb, t)
                    true_d = abs(Sa - Sb)
                    # lower = the column that is DEEPER (larger level down)
                    low = a if ka > kb else b
                    tl, ll = cols[low]
                    ph, fb = phi_of(tl, ll, t)
                    if fb: n_fb += 1
                    pkt_d = (n - 1) + ph
                    tw_true.append(true_d); tw_pkt.append(pkt_d)
                    if prev_n.get((a, b), 0) == 0:      # birth of a step
                        b_true.append(true_d); b_pkt.append(pkt_d)
                prev_n[(a, b)] = n
            t += sample_dt
        med = lambda v: (sorted(v)[len(v)//2] if v else float('nan'))
        mean = lambda v: (sum(v)/len(v) if v else float('nan'))
        print(f"{t0:5.0f}-{t1:<6.0f} {100*n_step/max(1,n_samp):7.1f} "
              f"{100*n_multi/max(1,n_step):7.1f} "
              f"{med(b_true):11.3f} {med(b_pkt):10.3f} "
              f"{med(tw_true):11.3f} {med(tw_pkt):10.3f} "
              f"{100*n_fb/max(1,n_step):10.1f}")
        print(f"{'':12} {'':7} {'':7} {mean(b_true):11.3f} {mean(b_pkt):10.3f} "
              f"{mean(tw_true):11.3f} {mean(tw_pkt):10.3f}   <- MEANS")

if __name__ == "__main__":
    B = [(50, 100), (100, 150), (150, 200)]
    D = "/Users/tzetze20/amr_tools/alamo/tests/MMWSpalling/studies"
    run(f"{D}/d2l_scan/output/A_f02_removal_events.csv", 2e-3, 40, 0.025, B, "f=0.2  2 mm (cured)", x0=0.0)
    run(f"{D}/d2j0b_plane/output/C1_2mm_removal_events.csv", 2e-3, 40, 0.025, B, "f=0    2 mm (cascade)", x0=0.0)
    run(f"{D}/d2l_scan/output/A_f02_1mm_removal_events.csv", 1e-3, 80, 0.025, B, "f=0.2  1 mm (cured)", x0=0.0)
