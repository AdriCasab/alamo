#!/usr/bin/env python3
"""S1 surface-resolution study: analysis -> RESULTS.md.

  analyze.py                  reads output/, writes RESULTS.md
  analyze.py --out output_smoke --results RESULTS_smoke.md

Everything above the `<!-- DISCUSSION -->` marker in RESULTS.md is
regenerated; the hand-written discussion below it is preserved.

Per run (plot_file P of run_sweep.py):
  P/thermo.dat            ledger_* (every 0.1 s)
  P_h_col_events.csv      EDGE-TRIGGERED: one row per column's FIRST firing at
                          each new top cell (time, h_col); repeat firings at
                          the same top cell are not logged
  P/NNNNNcell             plotfiles every 2 s (Temp, H, Lambda_L, removed, phi)

Surface history: level-set height z_c(k_top) + phi(k_top) at each plotfile
(k_top from `removed`); the scoring window and its halves are snapped to
plotfile times. (The CSV cannot give the surface: sum(h_col) misses the
unlogged repeat firings.)

Voided cells keep their H/T/Lambda frozen (no reset; Removal.H), so the final
plotfile holds every removed cell's state AT THE STEP IT WAS VOIDED:
  * <dT_rem> = mean H of cells voided in the window / (rho*Cp)
               (H = 0 at T_ref = T0 = 293.15 K; includes latent heat)
  * top-cell T at firing = frozen T of the cells voided in the window. When
               CSV rows == voided cells (`1-cell` column: no event voided two
               cells) every voided cell was the top solid cell at its firing
               step, so this is exact; otherwise it also samples sub-top cells.
  * flake = h_col of CSV rows in the window.
"""
import argparse
import glob
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from run_sweep import CASES, DEPTH, dz_tag, H_ROBIN  # noqa: E402

RHO, CP, K = 2750.0, 790.0, 1.5
KAPPA = K / (RHO * CP)
T0 = 293.15
BOTTOM_EXCL = 0.020
SURFACE_FLUX_CASES = ("F1s", "R-cell", "R-face", "R-face-1400")
MARK = "<!-- DISCUSSION -->"


def load_thermo(path):
    with open(path) as f:
        names = f.readline().split()
    rows = np.loadtxt(path, skiprows=1, ndmin=2)
    return {n: rows[:, i] for i, n in enumerate(names)}


def load_events(path):
    if not os.path.exists(path):
        return np.zeros(0), np.zeros(0)
    d = np.genfromtxt(path, delimiter=",", names=True, ndmin=1)
    if d.size == 0:
        return np.zeros(0), np.zeros(0)
    t, h = np.atleast_1d(d["time"]), np.atleast_1d(d["h_col"])
    keep = h > 0.0
    return t[keep], h[keep]


def load_pf(pf):
    import yt
    yt.funcs.mylog.setLevel(50)
    ds = yt.load(pf)
    cg = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    col = lambda f: np.asarray(cg[f])[0, 0, :]
    return dict(t=float(ds.current_time), T=col("Temp"), H=col("H"), LL=col("Lambda_L"),
                rem=col("removed") > 0.5, phi=col("phi"))


def pct(a, q):
    return float(np.percentile(a, q)) if len(a) else float("nan")


def analyze_run(case, dz_mm, P):
    dz = dz_mm * 1e-3
    N = int(round(DEPTH / dz))
    zc = (np.arange(N) + 0.5) * dz
    r = dict(case=case, dz=dz_mm, ok=False)
    if not os.path.isdir(P) or not os.path.exists(P + ".done"):
        r["status"] = "missing"
        return r
    th = load_thermo(os.path.join(P, "thermo.dat"))
    te, he = load_events(P + "_h_col_events.csv")
    pfs = sorted(glob.glob(os.path.join(P, "*cell")),
                 key=lambda p: int(re.search(r"(\d+)cell$", p).group(1)))
    pdat = [load_pf(p) for p in pfs]
    last = pdat[-1]
    t_end = last["t"]
    r.update(ok=True, t_end=t_end, n_events=len(te))
    r["ledger_err"] = float(np.max(np.abs(th["ledger_err"])))

    # Surface (level-set height) at each plotfile: z_c(k_top) + phi(k_top),
    # k_top from the `removed` mask.
    def z_surface(d):
        solid = np.nonzero(~d["rem"])[0]
        if len(solid) == 0:
            return 0.0
        k = solid.max()
        return zc[k] + d["phi"][k]
    tp = np.array([d["t"] for d in pdat])
    zs = np.array([z_surface(d) for d in pdat])
    nvoid = np.array([int(d["rem"].sum()) for d in pdat])
    r["depth_mm"] = (DEPTH - zs[-1]) * 1e3
    r["n_voided"] = int(nvoid[-1])
    # Every CSV row is a column's first firing at a new top cell; rows ==
    # voided cells means no event voided more than one cell, so every voided
    # cell was the top solid cell at its firing step.
    r["single_cell"] = (len(te) == nvoid[-1])
    if nvoid[-1] == 0:
        r["status"] = "no removal"
        return r

    # First removal from the CSV (first firing whose cell was voided): the
    # first row's time, since row n <-> the n-th voided cell.
    t_first = float(te[0]) if len(te) else float(tp[np.argmax(nvoid > 0)])
    t0_req = max(t_first + 10.0, 0.3 * t_end)
    if t0_req >= 0.8 * t_end:            # short runs (M-*): see RESULTS header
        t0_req = max(t_first + 2.0, 0.3 * t_end)
    i0 = int(np.searchsorted(tp, t0_req - 1e-9))
    i1 = len(tp) - 1
    if i0 >= i1:
        i0 = max(0, i1 - 1)
    im = i0 + int(np.argmin(np.abs(tp[i0:i1 + 1] - 0.5 * (tp[i0] + tp[i1]))))
    t0, t1, tm = tp[i0], tp[i1], tp[im]
    r.update(t_first=t_first, t0=t0, t1=t1)
    mh = 3600.0
    r["rop"] = (zs[i0] - zs[i1]) / (t1 - t0) * mh
    r["rop_a"] = (zs[i0] - zs[im]) / (tm - t0) * mh if tm > t0 else float("nan")
    r["rop_b"] = (zs[im] - zs[i1]) / (t1 - tm) * mh if t1 > tm else float("nan")

    # Cells voided inside the window (mask difference), frozen at firing.
    cells_w = last["rem"] & ~pdat[i0]["rem"]
    r["n_cells_w"] = int(cells_w.sum())
    Hw = last["H"][cells_w]
    Tw = last["T"][cells_w]
    r["dT_rem"] = float(np.mean(Hw) / (RHO * CP)) if len(Hw) else float("nan")
    r.update(Ttop_mean=float(np.mean(Tw)) if len(Tw) else float("nan"),
             Ttop_p10=pct(Tw, 10), Ttop_p90=pct(Tw, 90), n_top=len(Tw))
    win = (te > t0) & (te <= t1)
    hw = he[win]
    r["n_ev_w"] = int(win.sum())
    r.update(h_med=pct(hw, 50) * 1e3, h_p10=pct(hw, 10) * 1e3, h_p90=pct(hw, 90) * 1e3)

    # Absorbed flux (window mean) and losses.
    A = dz * dz
    tt = th["time"]
    sw = (tt > t0) & (tt <= t1)
    src = CASES[case][0]
    Pabs = th["ledger_P_robin"] if src["kind"] == "robin" else th["ledger_P_beam"]
    r["q_abs"] = float(np.mean(Pabs[sw]) / A) if sw.any() else float("nan")
    r["q_loss"] = float(np.mean(th["ledger_P_loss"][sw]) / A) if sw.any() else float("nan")

    # Implied true-surface T from the top-cell average (exponential layer).
    v = r["rop"] / mh
    if case in SURFACE_FLUX_CASES and v > 0 and np.isfinite(r["Ttop_mean"]):
        delta = KAPPA / v
        x = dz / delta
        r["delta_mm"] = delta * 1e3
        r["Ts_implied"] = T0 + (r["Ttop_mean"] - T0) * x / (1.0 - math.exp(-x))

    # Max T and melt, over all plotfiles, excluding the bottom 20 mm.
    Tmax_solid, melt_t, melt_depth = -np.inf, None, []
    keep = zc > BOTTOM_EXCL
    for d in pdat:
        solid = (~d["rem"]) & keep
        if solid.any():
            Tmax_solid = max(Tmax_solid, float(d["T"][solid].max()))
            ktop = np.nonzero(~d["rem"])[0].max()
            z_face = zc[ktop] + 0.5 * dz
            m = solid & (d["LL"] > 0.0)
            if m.any():
                if melt_t is None:
                    melt_t = d["t"]
                melt_depth.extend(list(z_face - zc[m]))
    voided = last["rem"] & keep
    r["Tmax_solid"] = Tmax_solid
    r["Tmax_void"] = float(last["T"][voided].max()) if voided.any() else float("nan")
    r["melt_solid"] = (melt_t, (min(melt_depth) * 1e3, max(melt_depth) * 1e3) if melt_depth else None)
    r["melt_removed_cells"] = int(np.count_nonzero(voided & (last["LL"] > 0.0)))
    r["missing_cols_max"] = float(np.max(th["ledger_surface_missing_cols"]))
    r["status"] = "ok"
    return r


def f(x, fmt=".3g"):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "-"
    return format(x, fmt)


def trend(vals):
    v = [x for x in vals if np.isfinite(x)]
    if len(v) < 3:
        return "too few points"
    d = np.diff(v)
    if np.all(d > 0):
        return "monotone increasing as dz decreases"
    if np.all(d < 0):
        return "monotone decreasing as dz decreases"
    return "non-monotone"


def convergence(rows, key, label):
    rs = [r for r in rows if r.get("status") == "ok" and np.isfinite(r.get(key, np.nan))]
    rs.sort(key=lambda r: -r["dz"])
    if len(rs) < 2:
        return f"- {label}: not enough converged runs"
    a, b = rs[-2], rs[-1]
    ch = (b[key] - a[key]) / a[key] if a[key] != 0 else float("nan")
    tag = "< 10% (converged)" if abs(ch) < 0.10 else ">= 10% (NOT converged)"
    return (f"- {label}: {a[key]:.4g} (dz = {a['dz']:g} mm) -> {b[key]:.4g} (dz = {b['dz']:g} mm), "
            f"change {ch * 100:+.1f}% {tag}; {trend([r[key] for r in rs])}")


def write_results(out, results_path, all_rows):
    L = []
    L.append("# S1 surface-resolution study: results\n")
    L.append(f"Generated by `analyze.py` from `{os.path.relpath(out, HERE)}/`. "
             "1-D columns (1 x 1 x N cubic cells, 250 mm deep, Neumann-0), sp_weibull "
             "top_cell spallation with melt_aware, losses on (eps = 0.8, h_c = 10), "
             "surface.follow_mask = 1, energy ledger on, dt = 1 ms, 1 rank per run.\n")
    L.append("**Scoring window:** from max(first removal + 10 s, 0.3 * t_end) to t_end, "
             "snapped to the 2 s plotfile times (first removal = first h_col CSV row). "
             "For runs too short for that (start later than 0.8 * t_end: the ~11 s M-* runs) "
             "the start is max(first removal + 2 s, 0.3 * t_end). ROP = level-set "
             "recession between window plotfiles; halves split at the plotfile nearest "
             "the midpoint. Top-cell T at firing = frozen T of all cells voided in the "
             "window (n = cells); exact top-cell sample when `1-cell` = yes (h_col CSV rows "
             "== voided cells, i.e. no event voided two cells). <dT_rem> = mean frozen "
             "H/(rho Cp) of the same cells. Flake = h_col of CSV rows in the window "
             "(first firing at each new top cell). T_s implied = T0 + (T_top - T0) x/(1 - e^-x), "
             "x = dz/delta, delta = kappa/ROP (surface-flux cases). Melt columns exclude the "
             "bottom 20 mm; `solid` = Lambda_L > 0 in a still-solid cell of any 2 s "
             "plotfile (first time, depth range below the mask top face); `removed` = "
             "voided cells whose frozen Lambda_L > 0.\n")
    for case, (src, dzs, v) in CASES.items():
        rows = [r for r in all_rows if r["case"] == case]
        if not rows:
            continue
        desc = (f"beam q = {src['q'] / 1e6:g} MW/m^2, alpha = {src['alpha']:g}/m"
                if src["kind"] == "beam" else
                f"Robin h = {H_ROBIN:g} W/m^2K, T_gas = {src['T_gas']:g} K, robin_form = {src['form']}")
        L.append(f"\n## {case}: {desc}\n")
        L.append("| dz mm | t_end s | window s | depth mm | ROP m/h (1st / 2nd half) | "
                 "<dT_rem> K | <q_abs> MW/m^2 | <q_loss> kW/m^2 | top-cell T at firing K "
                 "mean [p10, p90] (n) | T_s implied K | flake h mm median [p10, p90] | "
                 "flake cells median | Tmax K (solid snapshots / voided at firing) | melt solid: t s, depth mm | melt removed cells | "
                 "max abs ledger_err | 1-cell |")
        L.append("|" + "---|" * 18)
        for r in sorted(rows, key=lambda r: -r["dz"]):
            if r.get("status") != "ok":
                L.append(f"| {r['dz']:g} | {r.get('status')} |" + " |" * 16)
                continue
            ms = r["melt_solid"]
            melt = "none" if ms[0] is None else f"{ms[0]:g}, {ms[1][0]:.3g}-{ms[1][1]:.3g}"
            L.append(
                f"| {r['dz']:g} | {r['t_end']:.1f} | {r['t0']:.1f}-{r['t1']:.1f} | {r['depth_mm']:.1f} | "
                f"{r['rop']:.3g} ({r['rop_a']:.3g} / {r['rop_b']:.3g}) | {f(r['dT_rem'], '.0f')} | "
                f"{r['q_abs'] / 1e6:.3g} | {r['q_loss'] / 1e3:.3g} | "
                f"{f(r['Ttop_mean'], '.0f')} [{f(r['Ttop_p10'], '.0f')}, {f(r['Ttop_p90'], '.0f')}] ({r['n_top']}) | "
                f"{f(r.get('Ts_implied'), '.0f')} | "
                f"{f(r['h_med'])} [{f(r['h_p10'])}, {f(r['h_p90'])}] | {f(r['h_med'] / r['dz'], '.2g')} | "
                f"{r['Tmax_solid']:.0f} / {f(r['Tmax_void'], '.0f')} | {melt} | {r['melt_removed_cells']} | "
                f"{r['ledger_err']:.1e} | {'yes' if r['single_cell'] else 'no'} |")
        L.append("\nConvergence (two finest dz):\n")
        L.append(convergence(rows, "rop", "ROP [m/h]"))
        L.append(convergence(rows, "h_med", "median flake [mm]"))
        L.append(convergence(rows, "Ttop_mean", "top-cell T at firing [K]"))
        L.append(convergence(rows, "q_abs", "<q_abs> [W/m^2]"))
        if case in SURFACE_FLUX_CASES:
            rs = sorted([r for r in rows if r.get("status") == "ok"], key=lambda r: r["dz"])
            if rs:
                fin = rs[0]
                L.append(f"- implied T_s vs resolved top-cell T at the finest dz ({fin['dz']:g} mm): "
                         + ", ".join(f"dz {r['dz']:g}: {f(r.get('Ts_implied'), '.0f')}"
                                     for r in sorted(rs, key=lambda r: -r["dz"]))
                         + f" vs {f(fin['Ttop_mean'], '.0f')} K (delta = kappa/v = "
                         f"{f(fin.get('delta_mm'), '.3g')} mm at the finest dz)")
    text = "\n".join(L) + "\n"
    discussion = f"\n{MARK}\n\n## Discussion\n\n(to be written)\n"
    if os.path.exists(results_path):
        old = open(results_path).read()
        if MARK in old:
            discussion = "\n" + MARK + old.split(MARK, 1)[1]
    with open(results_path, "w") as fh:
        fh.write(text + discussion)
    print(f"wrote {results_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="output")
    ap.add_argument("--results", default="RESULTS.md")
    args = ap.parse_args()
    out = os.path.join(HERE, args.out)
    rows = []
    for case, (_, dzs, _) in CASES.items():
        for dz_mm in dzs:
            P = os.path.join(out, f"{case}_{dz_tag(dz_mm)}")
            if not os.path.isdir(P):
                continue
            r = analyze_run(case, dz_mm, P)
            rows.append(r)
            print(f"{case:12s} dz={dz_mm:<7g} {r.get('status')}"
                  + (f"  ROP {r['rop']:.3g} m/h  dT_rem {r['dT_rem']:.0f} K  Ttop {r['Ttop_mean']:.0f} K  "
                     f"h_med {r['h_med']:.3g} mm  q_abs {r['q_abs'] / 1e6:.3g} MW/m2  "
                     f"ledger {r['ledger_err']:.1e}  1-cell {r['single_cell']}  depth {r['depth_mm']:.0f} mm"
                     if r.get("status") == "ok" else ""), flush=True)
    write_results(out, os.path.join(HERE, args.results), rows)


if __name__ == "__main__":
    main()
