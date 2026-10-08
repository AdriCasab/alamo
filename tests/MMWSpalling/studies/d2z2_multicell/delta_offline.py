#!/usr/bin/env python3
"""D2z-2 offline both-column delta estimator. READ-ONLY.

  S_i / dz = n_i + phi_i ,  phi_i = clamp((t - t_last_i)/dt_i, 0, 1)
  delta    = n + phi_S - phi_T ,  n = nvoid_S - nvoid_T  (S is deeper/shorter)

⚠⚠ THIS IS NOT D2z's delta_offline.py. Criterion (b) compares this to the LIVE
`delta_mean`, and D2z's version specified its population DIFFERENTLY from the
live diagnostic, which is why that comparison failed for population reasons
rather than physics ones. The two differences, both fixed here:

  1. D2z's `voids()` DROPPED every column with fewer than two removals
     (`len(v) >= 2`). The live `ColPhi` KEEPS them, at phi = 1. So never-fired
     and once-fired columns are included here, at phi = 1, and a column that
     never appears in the event log at all is included too -- it is the TALLEST
     column there is (nvoid = 0), and on Meier those sit right next to active
     ones well inside any scoring radius.
  2. D2z required BOTH columns inside the radius. The live gate is on the
     column that owns the riser -- the TALLER one. The radius is applied to the
     taller column here, the same way.

Absolute k_top is never needed: every column starts at the same flat surface,
so k_top = k0 - nvoid and n = nvoid_S - nvoid_T.

⚠ Arbiter runs are a QUARTER domain with the disk at the ORIGIN corner
(x0 = y0 = 0). Meier is centred; pass the centre explicitly.
⚠ Filter n_voided > 0.5 on the EVENT rows: a row with no voided cell is a
firing that removed nothing and does not move the column.
"""
import csv
import math
import sys
from collections import defaultdict


def load(path, nxy, dx, x0=0.0, y0=0.0):
    """(i,j) -> sorted [t] of VOID events; EVERY column in the grid is present,
    including the ones that never fired (empty list)."""
    ev = {(i, j): [] for i in range(nxy) for j in range(nxy)}
    with open(path) as f:
        for r in csv.DictReader(f):
            if float(r["n_voided"]) <= 0.5:
                continue
            c = (int(r["col_i"]), int(r["col_j"]))
            if c in ev:
                ev[c].append(float(r["time"]))
    for v in ev.values():
        v.sort()
    rad = {c: math.hypot((c[0] + 0.5) * dx - x0, (c[1] + 0.5) * dx - y0) for c in ev}
    return ev, rad


def state(ts, t):
    """(n_voided_so_far, phi, fallback) with phi from the TRAILING interval.

    ⚠ n < 2 returns phi = 1, matching ColPhi -- NOT a dropped column."""
    lo, hi = 0, len(ts)
    while lo < hi:
        m = (lo + hi) // 2
        if ts[m] <= t:
            lo = m + 1
        else:
            hi = m
    n = lo
    if n < 2:
        return n, 1.0, True
    dt = ts[n - 1] - ts[n - 2]
    if dt <= 0.0:
        return n, 1.0, True
    return n, min(1.0, max(0.0, (t - ts[n - 1]) / dt)), False


def run(path, nxy, dx, r_max, blocks, label, x0=0.0, y0=0.0, sample_dt=0.25,
        quiet=False, require_tall_fired=False):
    """`require_tall_fired` drops pairs whose TALLER column has fewer than two
    removals.

    ⚠ MEASURED, not chosen: the live population contains no such pair. Live
    `delta_fallback_cols` is exactly 0 on every closure leg, and ColPhi sets
    that flag for any column with nvoid < 2, so the live diagnostic never sees
    one. The packet's instruction to include them on both sides is right about
    the PHI CONVENTION (phi = 1) and wrong about INCLUSION: a never-fired column
    at the patch edge sits beside a column that has receded 40 cells, so the
    pair carries delta = 40 and it is the crater WALL, not a phase lag. Including
    them moved the R2w r<30 offline mean from 0.540 to 3.123 while the live
    diagnostic read 0.504.

    Both populations are reported by analyze.py; the pre-registered one is
    scored. See RESULTS.md.
    """
    ev, rad = load(path, nxy, dx, x0, y0)
    pairs = [(c, (c[0] + d[0], c[1] + d[1]))
             for c in ev for d in ((1, 0), (0, 1)) if (c[0] + d[0], c[1] + d[1]) in ev]
    if not quiet:
        print(f"\n=== {label}  ({len(ev)} cols, {len(pairs)} adjacent pairs, "
              f"taller column r < {r_max*1e3:.0f} mm) ===")
        print(f"{'block':>12} {'step%':>7} {'>=2c%':>7} {'birth med':>10} "
              f"{'t-wtd med':>10} {'t-wtd MEAN':>11} {'fallback%':>10} {'maxdelta':>9}")
    out = {}
    for (t0, t1) in blocks:
        tw, birth = [], []
        n_s = n_step = n_multi = n_fb = 0
        prev = {}
        t = t0
        while t <= t1:
            for (a, b) in pairs:
                na, pa, fa = state(ev[a], t)
                nb, pb, fb = state(ev[b], t)
                n_s += 1                      # EVERY sample, so step% is real
                # more voids => deeper => SHORTER. n = nvoid_S - nvoid_T.
                if na == nb:
                    prev[(a, b)] = 0
                    continue
                if na > nb:
                    cS, cT, nn, pS, pT, fS, fT = a, b, na - nb, pa, pb, fa, fb
                else:
                    cS, cT, nn, pS, pT, fS, fT = b, a, nb - na, pb, pa, fb, fa
                # the radius is applied to the TALLER column, as the live gate is
                if rad[cT] > r_max:
                    prev[(a, b)] = 0
                    continue
                if require_tall_fired and len(ev[cT]) < 2:
                    prev[(a, b)] = 0
                    continue
                n_step += 1
                if nn >= 2:
                    n_multi += 1
                if fS or fT:
                    n_fb += 1
                d = nn + pS - pT
                if d < 0.0:
                    d = 0.0
                tw.append(d)
                if prev.get((a, b), 0) == 0:
                    birth.append(d)
                prev[(a, b)] = nn
            t += sample_dt
        med = lambda v: (sorted(v)[len(v) // 2] if v else float("nan"))
        mean = lambda v: (sum(v) / len(v) if v else float("nan"))
        key = f"{int(t0)}-{int(t1)}"
        out[key] = dict(step=n_step / max(1, n_s), multi=n_multi / max(1, n_step),
                        birth_med=med(birth), tw_med=med(tw), tw_mean=mean(tw),
                        fb=n_fb / max(1, n_step), n=n_step,
                        dmax=(max(tw) if tw else float("nan")))
        if not quiet:
            o = out[key]
            print(f"{key:>12} {100*o['step']:7.1f} {100*o['multi']:7.1f} "
                  f"{o['birth_med']:10.3f} {o['tw_med']:10.3f} {o['tw_mean']:11.3f} "
                  f"{100*o['fb']:10.2f} {o['dmax']:9.3f}")
    return out


if __name__ == "__main__":
    print(__doc__)
    print("Import `run` from analyze.py; this module has no standalone campaign.")
