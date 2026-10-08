#!/usr/bin/env python3
"""D2z item 0a -- the OFFLINE both-column delta estimator. READ-ONLY.

  S_i / dz = n_i + phi_i ,  phi_i = clamp((t - t_last_i)/dt_i, 0, 1)
  delta    = | n + phi_S - phi_T |        n = n_S - n_T  (shorter S is deeper)

This is the reference criterion (d) scores the LIVE delta_mean against, so it
uses the SAME trailing-interval rule the live code must use -- the column's own
previous interval, never the next one in hindsight.

⚠ (n-1) + phi_lower is WRONG: it reads one column's phase and a step is born the
instant the deeper column flips, so it returns ~0 by construction.
⚠ Arbiter runs are a QUARTER domain with the disk at the ORIGIN corner
(x0 = y0 = 0); centring on the domain middle silently yields zero columns.
⚠ Filter n_voided > 0.5. Trivially true under whole-cell, but it is the filter
that separated 0.034 from 0.716.
"""
import csv, math, sys
from collections import defaultdict


def voids(path, dx, r_max, x0=0.0):
    """(i,j) -> sorted [t] of VOID events inside the gate disc."""
    ev = defaultdict(list)
    with open(path) as f:
        for r in csv.DictReader(f):
            if float(r["n_voided"]) <= 0.5:          # the load-bearing filter
                continue
            i, j = int(r["col_i"]), int(r["col_j"])
            if math.hypot((i + 0.5) * dx - x0, (j + 0.5) * dx - x0) > r_max:
                continue
            ev[(i, j)].append(float(r["time"]))
    for v in ev.values():
        v.sort()
    return {c: v for c, v in ev.items() if len(v) >= 2}


def state(ts, t):
    """(n_removed, phi, fallback) with phi from the TRAILING interval."""
    lo, hi = 0, len(ts)
    while lo < hi:
        m = (lo + hi) // 2
        if ts[m] <= t: lo = m + 1
        else: hi = m
    n = lo                       # voids so far
    if n < 2:
        return n, 1.0, True      # no interval yet -> phi = 1 -> delta = cell count
    dt = ts[n - 1] - ts[n - 2]   # trailing, not hindsight
    if dt <= 0.0:
        return n, 1.0, True
    return n, min(1.0, max(0.0, (t - ts[n - 1]) / dt)), False


def run(path, dx, r_max, blocks, label, sample_dt=0.05):
    cols = voids(path, dx, r_max)
    keys = set(cols)
    pairs = [(c, (c[0] + d[0], c[1] + d[1])) for c in keys for d in ((1, 0), (0, 1))
             if (c[0] + d[0], c[1] + d[1]) in keys]
    print(f"\n=== {label}  ({len(cols)} cols, {len(pairs)} adjacent pairs, "
          f"r < {r_max*1e3:.0f} mm) ===")
    print(f"{'block':>12} {'step%':>7} {'>=2c%':>7} {'birth med':>10} {'birth mean':>11} "
          f"{'t-wtd med':>10} {'t-wtd MEAN':>11} {'fallback%':>10}")
    out = {}
    for (t0, t1) in blocks:
        tw, birth = [], []
        n_s = n_step = n_multi = n_fb = 0
        prev = {}
        t = t0
        while t <= t1:
            for (a, b) in pairs:
                na, pa, fa = state(cols[a], t)
                nb, pb, fb = state(cols[b], t)
                n_s += 1
                nd = na - nb
                if nd == 0:
                    prev[(a, b)] = 0
                    continue
                n_step += 1
                if abs(nd) >= 2: n_multi += 1
                if fa or fb: n_fb += 1
                d = abs(nd + (pa - pb)) if nd > 0 else abs(-nd + (pb - pa))
                tw.append(d)
                if prev.get((a, b), 0) == 0:
                    birth.append(d)
                prev[(a, b)] = nd
            t += sample_dt
        med = lambda v: (sorted(v)[len(v)//2] if v else float("nan"))
        mean = lambda v: (sum(v)/len(v) if v else float("nan"))
        key = f"{int(t0)}-{int(t1)}"
        out[key] = dict(step=n_step/max(1, n_s), multi=n_multi/max(1, n_step),
                        birth_med=med(birth), birth_mean=mean(birth),
                        tw_med=med(tw), tw_mean=mean(tw), fb=n_fb/max(1, n_step))
        o = out[key]
        print(f"{key:>12} {100*o['step']:7.1f} {100*o['multi']:7.1f} "
              f"{o['birth_med']:10.3f} {o['birth_mean']:11.3f} "
              f"{o['tw_med']:10.3f} {o['tw_mean']:11.3f} {100*o['fb']:10.2f}")
    return out


if __name__ == "__main__":
    D = "tests/MMWSpalling/studies/d2z1_delete/output"
    B = [(50, 100), (100, 150), (150, 200), (50, 200)]
    run(f"{D}/W2_removal_events.csv", 2e-3, 0.025, B, "W2  2 mm  whole-cell (re-baselined)")
    run(f"{D}/W1_removal_events.csv", 1e-3, 0.025, B, "W1  1 mm  whole-cell (re-baselined)")
