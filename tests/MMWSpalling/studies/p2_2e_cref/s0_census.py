#!/usr/bin/env python3
"""P2-2e deliverable 1(b) -- the FROZEN-COLUMN CENSUS, P1's definition.

P1's definition, verbatim from the packet: a column at r < 40 mm is FROZEN at
time t if its last regime-1 firing precedes t - 100 s WHILE >= 50 % of its 8
neighbours fired in the last 100 s. Scored at t = 200 / 350 / 500 s.

⚠ Self-contained on purpose. `reduce_price.py`'s `frozen_floor` uses
silent_s = 30 / active_s = 50, not P1's 100/100, and carries a module-level
2 mm `V_CELL`. The definition is the packet's, implemented here against the
event log directly.
⚠ Columns are indexed by (col_i, col_j); the radius uses the leg's own dx,
parsed from its argv.sh (never a module constant).

Usage: s0_census.py <leg-dir> [<leg-dir> ...]
"""
import csv
import os
import shlex
import sys

R_FLOOR = 1          # regime label: a vertical (floor) firing
R_ACTIVE = 40.0e-3   # the active radius, P1's 40 mm
SILENT = 100.0       # a column is silent if its last firing precedes t - 100 s
NEIGH_FRAC = 0.5     # >= 50 % of the 8 neighbours must have fired recently
TIMES = (200.0, 350.0, 500.0)


def dx_of(leg):
    toks = shlex.split(open(os.path.join(leg, "argv.sh")).read())
    kv, i = {}, 0
    while i < len(toks):
        t = toks[i]
        if "=" in t:
            k, v = t.split("=", 1)
            i += 1
            while i < len(toks) and "=" not in toks[i]:
                v += " " + toks[i]
                i += 1
            kv[k] = v
        else:
            i += 1
    nx = int(kv["amr.n_cell"].split()[0])
    lo = float(kv["geometry.prob_lo"].split()[0])
    hi = float(kv["geometry.prob_hi"].split()[0])
    return (hi - lo) / nx


def last_fired(leg):
    """{(i,j): last regime-1 firing time}, and the max event time."""
    p = os.path.join(leg, "run_removal_events.csv")
    last, tmax = {}, 0.0
    with open(p) as f:
        for r in csv.DictReader(f):
            try:
                reg = int(float(r["regime"] or 0))
                t = float(r["time"])
            except (ValueError, KeyError):
                continue
            if t > tmax:
                tmax = t
            if reg != R_FLOOR:
                continue
            key = (int(float(r["col_i"])), int(float(r["col_j"])))
            if t > last.get(key, -1.0):
                last[key] = t
    return last, tmax


def census(last, dx, t):
    """(n_frozen, n_pop, [frozen keys]) at time t."""
    pop = [k for k in last
           if dx * ((k[0] + 0.5) ** 2 + (k[1] + 0.5) ** 2) ** 0.5 < R_ACTIVE]
    frozen = []
    for (i, j) in pop:
        if last[(i, j)] >= t - SILENT:
            continue                       # not silent
        nb = [(i + a, j + b) for a in (-1, 0, 1) for b in (-1, 0, 1)
              if not (a == 0 and b == 0)]
        have = [n for n in nb if n in last]
        if not have:
            continue
        recent = sum(1 for n in have if last[n] >= t - SILENT)
        if recent >= NEIGH_FRAC * len(have):
            frozen.append((i, j))
    return len(frozen), len(pop), frozen


def report(leg, label=None):
    P = print
    dx = dx_of(leg)
    last, tmax = last_fired(leg)
    nm = label or os.path.basename(leg.rstrip("/"))
    P(f"  {nm}:  dx {dx*1e3:.3f} mm (from argv.sh), {len(last)} columns with a"
      f" regime-1 firing, last event t = {tmax:.3f} s")
    out = {}
    for t in TIMES:
        # ⚠ A leg that stops at 500 s has its last event a fraction of a step
        # earlier (S0: 499.984). Scoring t = 500 against it is correct -- the
        # window is 100 s wide -- and the archive's own P1 census had the same
        # tie (its note: "t_max 499.98 vs t_stop 500.0", a 3-column difference).
        # Only skip a time the leg genuinely never reached.
        if t > tmax + 30.0:
            P(f"    t = {t:5.0f} s  -- the leg stops at {tmax:.1f} s,"
              f" not scored")
            continue
        nf, npop, _ = census(last, dx, t)
        out[t] = (nf, npop)
        P(f"    t = {t:5.0f} s  frozen {nf:4d} of {npop:4d} active columns"
          f" (r < {R_ACTIVE*1e3:.0f} mm) = {100.0*nf/max(npop,1):5.2f} %")
    return out


def main():
    print("P2-2e (a1 b) -- the FROZEN-COLUMN CENSUS, P1's definition")
    print("=" * 94)
    print(f"  frozen at t: last regime-1 firing < t - {SILENT:.0f} s AND"
          f" >= {NEIGH_FRAC*100:.0f} % of its 8 neighbours fired in the last"
          f" {SILENT:.0f} s")
    print(f"  population: columns at r < {R_ACTIVE*1e3:.0f} mm that have ever"
          f" fired regime 1")
    for leg in sys.argv[1:]:
        report(leg)
    print("=" * 94)


if __name__ == "__main__":
    main()
