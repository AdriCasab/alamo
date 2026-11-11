#!/usr/bin/env python3
"""Beam column resolution study -- the tables. Usage: reduce.py > RESULTS.txt   (run.py's docstring has the fixture)

Timing: a thermo row stamped t is the step that ENDED at t (written every 5th step, 0.004 s); an events row stamps the
step's START. Columns read (thermo): bc_T_peak (the column's hottest CELL value, K); bc_z_melt (finite once a cell
reaches T_m = 1510.5 K -> the melt onset); bc_z_liq (centre depth of the deepest fully liquid cell; + dz/2 = the liquid
layer in whole cells); bc_ell_p (the profile's own e-fold depth from the surface temperature, m); fr_smax_floor (the
largest onset score); bk_n_eval / bk_n_block (lines that met the onset criterion this step / of those, withdrawn by the
buckling condition). Events CSV (leg A): one row per fired line (regime 1), 16 identical columns.
"""
import csv
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run as R

NCOL = 16
DMM = {"0p5": 0.5, "1": 1.0, "2": 2.0, "4": 4.0, "8": 8.0}


def thermo(leg):
    p = os.path.join(R.LEGS, leg, "run", "thermo.dat")
    rcp = os.path.join(R.LEGS, leg, "rc")
    if not (os.path.exists(rcp) and open(rcp).read().strip() == "0" and os.path.exists(p)):
        return None, None
    L = open(p).read().splitlines()
    h = [x.strip() for x in L[0].split("\t")]
    x = np.array([[float(v) for v in r.split("\t")] for r in L[1:] if r.strip()])
    return {n: i for i, n in enumerate(h)}, x[x[:, 0] > 0.0]


def first(x, mask):
    return float(x[np.argmax(mask), 0]) if mask.any() else None


def at(x, c, col, t):
    if t is None or t > x[-1, 0] + 1e-9:
        return None
    return float(x[np.argmin(np.abs(x[:, 0] - t)), c[col]])


def heat(leg, dz):
    c, x = thermo(leg)
    if c is None:
        return None
    ev = (x[:, c["bk_n_eval"]] > 0) | (x[:, c["fr_smax_floor"]] >= 1.0)
    tc = first(x, ev)
    o = {"t_crit": tc, "wd": None, "t_melt": first(x, np.isfinite(x[:, c["bc_z_melt"]])),
         "t_liq": first(x, np.isfinite(x[:, c["bc_z_liq"]])), "T05": at(x, c, "bc_T_peak", 0.5),
         "T10": at(x, c, "bc_T_peak", 1.0), "Tstop": float(x[-1, c["bc_T_peak"]]), "stop": float(x[-1, 0])}
    if tc is not None:
        i = np.argmax(ev)
        o["wd"] = bool(x[i, c["bk_n_block"]] > 0)
        o["ellp"] = float(x[i, c["bc_ell_p"]])
        o["Tc"] = float(x[i, c["bc_T_peak"]])
        blk = x[:, c["bk_n_block"]] > 0
        ok = (x[:, c["bk_n_eval"]] > 0) & ~blk
        o["t_pass"] = first(x, ok)                     # the first step a fired line is NOT withdrawn
    for tag, dtm in (("L5", 5.0), ("L10", 10.0)):
        v = at(x, c, "bc_z_liq", o["t_melt"] + dtm) if o["t_melt"] is not None else None
        o[tag] = (v + 0.5 * dz) if (v is not None and math.isfinite(v)) else (0.0 if v is not None else None)
    return o


def model(leg, dz):
    c, x = thermo(leg)
    if c is None:
        return None
    ev = [r for r in csv.DictReader(open(os.path.join(R.LEGS, leg, "run_removal_events.csv"))) if r["regime"] == "1"]
    tm = first(x, np.isfinite(x[:, c["bc_z_melt"]]))
    o = {"n": len(ev) / NCOL, "t_melt": tm, "blk": float(np.sum(x[:, c["bk_n_block"]]) > 0),
         "vap": float(np.nansum(x[:, c["bc_n_vapor"]]))}
    if ev:
        e = ev[0]
        o.update(t=float(e["time"]), T_top=float(e["T_top"]), T_surf=float(e["T_surf"]), ell=float(e["ell_eff"]),
                 bk=float(e["bk"]), ellp=float(e["ell_p"]),
                 depth=sum(float(r["h_applied"]) for r in ev) / NCOL,
                 n_pre=sum(1 for r in ev if tm is None or float(r["time"]) < tm) / NCOL)
    nb = o.get("n_pre", 0.0)
    o["cls"] = "S" if nb >= 3 else "St" if nb >= 1 else "M" if tm is not None else "B" if o["blk"] else "C"
    return o


def f(v, fmt="%.3f", sc=1.0):
    return "-" if v is None or (isinstance(v, float) and not math.isfinite(v)) else fmt % (sc * v)


def richardson(vals):
    """the limit from the three finest values that exist, when they move one way with a shrinking change"""
    v = [a for a in vals if a is not None]
    if len(v) < 3:
        return None
    a, b, c = v[-3:]
    d1, d2 = b - a, c - b
    if d2 == 0.0:
        return c
    if d1 * d2 <= 0 or abs(d2) >= abs(d1):
        return None
    return c + d2 / (d1 / d2 - 1.0)


def main():
    p = print
    p("BEAM COLUMN RESOLUTION STUDY (binary %s)" % open(os.path.join(R.HERE, "BINARY.txt")).read().strip())
    p("fixture, legs and keys: run.py. dt 0.0008 s on every leg; 16 identical columns; random flaws off (the 8 um flaw).")
    p("A = the model as it is; H = heat only (no removal: the criterion is still evaluated and logged in thermo).")
    summary = []
    for d, q, stop in R.POINTS:
        D = DMM[d]
        p("")
        p("=" * 150)
        p("absorption depth D = %g mm, absorbed flux %s MW/m2, stop %g s" % (D, q, stop))
        p("  HEAT ONLY (H)")
        p("  %5s %7s %6s | %8s %8s %9s %8s | %8s %8s %8s | %8s %8s | %8s %8s %9s" % (
            "N", "dz mm", "dz/D", "t_crit s", "T_top K", "withdrawn", "t_pass s", "T(0.5) K", "T(1.0) K", "ell_p mm",
            "t_melt s", "t_liq s", "liq+5 mm", "liq+10mm", "T(stop) K"))
        H, A = {}, {}
        for n in R.NZ:
            dz = R.LZ / n
            H[n] = h = heat("D%s_Q%s_N%d_H" % (d, q, n), dz)
            if h is None:
                p("  %5d %7.4f %6.3f | not run" % (n, 1e3 * dz, 1e3 * dz / D))
                continue
            p("  %5d %7.4f %6.3f | %8s %8s %9s %8s | %8s %8s %8s | %8s %8s | %8s %8s %9s" % (
                n, 1e3 * dz, 1e3 * dz / D, f(h["t_crit"]), f(h.get("Tc"), "%.1f"),
                "-" if h["wd"] is None else ("YES" if h["wd"] else "no"), f(h.get("t_pass")),
                f(h["T05"], "%.1f"), f(h["T10"], "%.1f"), f(h.get("ellp"), "%.3f", 1e3), f(h["t_melt"]), f(h["t_liq"]),
                f(h["L5"], "%.3f", 1e3), f(h["L10"], "%.3f", 1e3), f(h["Tstop"], "%.0f")))
        p("  THE MODEL AS IT IS (A): the first fired line that was NOT withdrawn, firings per column, depth removed per column")
        p("  %5s %7s %6s | %8s %8s %8s %8s %7s %8s | %8s %9s %9s %8s | %5s" % (
            "N", "dz mm", "dz/D", "t_fire s", "T_top K", "T_surf K", "ell mm", "bk", "ell_p mm", "firings", "pre-melt",
            "depth mm", "t_melt s", "class"))
        for n in R.NZ:
            dz = R.LZ / n
            A[n] = a = model("D%s_Q%s_N%d_A" % (d, q, n), dz)
            if a is None:
                p("  %5d %7.4f %6.3f | not run" % (n, 1e3 * dz, 1e3 * dz / D))
                continue
            p("  %5d %7.4f %6.3f | %8s %8s %8s %8s %7s %8s | %8s %9s %9s %8s | %5s" % (
                n, 1e3 * dz, 1e3 * dz / D, f(a.get("t")), f(a.get("T_top"), "%.1f"), f(a.get("T_surf"), "%.1f"),
                f(a.get("ell"), "%.3f", 1e3), f(a.get("bk")), f(a.get("ellp"), "%.3f", 1e3), f(a["n"], "%.1f"),
                f(a.get("n_pre"), "%.1f"), f(a.get("depth"), "%.3f", 1e3), f(a["t_melt"]), a["cls"]))
        for key, lab in (("t_crit", "onset criterion met (H)"), ("t_melt", "melt onset (H)"), ("T10", "top cell T at 1 s (H)")):
            vals = [H[n][key] if H[n] else None for n in R.NZ]
            lim = richardson(vals)
            if lim is None:
                p("  %-26s no limit from the three finest (missing, not one-way, or not shrinking): %s" % (
                    lab, " ".join(f(v) for v in vals)))
                continue
            errs = [None if v is None else 100.0 * (v - lim) / lim for v in vals]
            p("  %-26s limit %.4g; error against it per cell: %s  [%%]" % (lab, lim, "  ".join(f(e, "%+.1f") for e in errs)))
            for n, e in zip(R.NZ, errs):
                if e is not None:
                    summary.append((key, D, q, 1e3 * R.LZ / n / D, e))
    p("")
    p("=" * 150)
    p("THE RULE: error against the extrapolated limit vs cells per absorption depth (all points)")
    for key, lab in (("t_crit", "onset criterion met"), ("t_melt", "melt onset"), ("T10", "top cell T at 1 s")):
        p("  %s" % lab)
        rows = sorted([s for s in summary if s[0] == key], key=lambda s: -s[3])
        for lo, hi in ((2.0, 99), (1.0, 2.0), (0.5, 1.0), (0.25, 0.5), (0.125, 0.25), (0.0, 0.125)):
            b = [s for s in rows if lo <= s[3] < hi]
            if b:
                p("    dz/D in [%.3f, %s): %2d legs, |error| max %.1f %%, median %.1f %%   (%s)" % (
                    lo, ("%.3f" % hi) if hi < 99 else "inf", len(b), max(abs(s[4]) for s in b),
                    float(np.median([abs(s[4]) for s in b])),
                    ", ".join("D%g/Q%s %.2f: %+.1f" % (s[1], s[2], s[3], s[4]) for s in b)))


if __name__ == "__main__":
    main()
