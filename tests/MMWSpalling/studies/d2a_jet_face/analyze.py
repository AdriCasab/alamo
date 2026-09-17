#!/usr/bin/env python3
"""D2a analysis -> RESULTS.md (tables regenerate above <!-- DISCUSSION -->; the
hand-written discussion below the marker is preserved) and depth-profile PNGs.

All recession numbers come from removal_events.csv (C1 takeaway 8): a column's
level-set depth at time t is the cumulative h_applied of its rows with time <= t
(the phi shift). Bin-mean depth over radial bins of 1 D (r of the column centre
from the jet axis) is sampled on a 0.5 s grid for the profiles. A bin's window
ROP is the mean over its columns of the C1 estimator: the least-squares slope of
the column's cumulative h_applied at its event rows inside the window (>= 3
rows spanning >= half the window; otherwise the depth secant over the window).
A/B comparisons (mesh, dt, sensitivities) add a band ROP: the centre-bin
depth increase over the time it takes to go from 20 % to 85 % of the pair's
smaller end centre depth. It removes the ignition-time offset (the time to first
firing is face-form and dt/mesh dependent), which shifts a feed-0 decelerating
face along its curve and biases fixed-time windows.

Windows: W_mid = [t_end/3, 2 t_end/3], W_end = [2 t_end/3, t_end] (t_end = the
run's stop time). Feed runs are scored on W_end; W_mid is shown for steadiness.
W_all = [t_end/6, t_end] is also used for the A/B comparisons (mesh, dt,
sensitivities): a 20 s window holds only ~6 firings of ~dz/2 per column in the
2-D runs.

Definitions:
  R_h (feed > 0)  outer edge of the contiguous run of bins from the axis whose
                  W_end ROP >= 0.9 v_feed.
  R_h (feed = 0)  r where the end bin-mean depth = 1/2 the centre-bin depth.
  centre s(t)     z_n(t) - mask top face of the column nearest the axis
                  (k_top - n_voided + 1)*dz after its last event row.
  T_fire          regime-1 T_top of columns with r <= max(R_h, 1 D) in W_end.
  face power      W_end mean of thermo patch_P_robin (x4 in the quarter domain).
  hole flux       mean over columns with r <= R_h of the pinned closed form
                  h(r,s)(T_gas(s) - T_pin) - eps sig (T_pin^4 - T_a^4) - h_c(T_pin - T_a)
                  at t_end (s from the column's mask top, T_pin its last regime-1
                  T_top; 0 for columns that never fired).
  branch shares   W_end thermo pinned_cols_* / (pinned + face), domain totals.
  idle proxy      per region, share of a column's consecutive regime-1 gaps
                  longer than idle * t_cell, t_cell = rho Cp dz (T_pin - T_a)/q_pin,
                  q_pin = h(r,s)(T_gas(s) - T_pin) (the BuildPinEffective rule).
  rim             columns with r < R_h (hole) that never fired.
  out of range    share of patch area with r > 7.5 D, r < 2.5 D, or s outside
                  [2, 12] D at t_end (Martin validity), all patch columns.
  max face slope  max |d(depth)/dr| between adjacent bins at t_end, as an angle.
"""
import ast
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import jet  # noqa: E402
from run import (ANCHORS, D, EPS, HC, OUT, PART1, PART2, RHOCP, SIG, TA, H_MARTIN,  # noqa: E402
                 q_closed)

MARK = "<!-- DISCUSSION -->"
P_BURNER = 38.0e3
T0 = TA
DT_GRID = 0.5


def meta_of(pf):
    try:
        return ast.literal_eval(open(pf + ".done").read())
    except Exception:
        return None


def thermo(pf):
    lines = [l.split() for l in open(os.path.join(pf, "thermo.dat")) if l.strip()]
    a = np.array(lines[1:], dtype=float)
    return {n: a[:, i] for i, n in enumerate(lines[0])}


def events(pf):
    ev = np.genfromtxt(pf + "_removal_events.csv", delimiter=",", names=True, ndmin=1)
    return ev


def geometry(m):
    dz = m["dz"] * 1e-3
    if m["dim"] == 2:
        nx, ny, nz = int(round(0.14 / dz)), int(round(0.008 / dz)), int(round(0.120 / dz))
    else:
        nx, ny, nz = 60, 60, 100
    xc = (np.arange(nx) + 0.5) * dz
    yc = (np.arange(ny) + 0.5) * dz
    X, Y = np.meshgrid(xc, yc, indexing="xy")          # [j, i]
    r = np.hypot(X - m["x0"], Y - m["y0"])
    return dz, nx, ny, nz, r


def fit(t, y, t0, t1):
    sel = (t >= t0) & (t <= t1)
    if sel.sum() < 3:
        return float("nan")
    return float(np.polyfit(t[sel], y[sel], 1)[0])


def analyze(pf, m):
    dz, nx, ny, nz, r = geometry(m)
    ncol = nx * ny
    th = thermo(pf)
    ev = events(pf)
    t_end = float(m["stop"])
    feed = m["feed"]
    z_n = lambda t: m["z_top"] + 7.0 * D - feed * t
    tg = np.arange(0.0, t_end + 1e-9, DT_GRID)
    W_mid, W_end = (t_end / 3.0, 2.0 * t_end / 3.0), (2.0 * t_end / 3.0, t_end)
    W_all = (t_end / 6.0, t_end)

    col = (ev["col_j"].astype(int) * nx + ev["col_i"].astype(int)) if ev.size else np.zeros(0, int)
    order = np.lexsort((ev["time"], col)) if ev.size else np.zeros(0, int)
    ev, col = ev[order], col[order]
    starts = np.searchsorted(col, np.arange(ncol + 1))

    depth_g = np.zeros((ncol, tg.size))       # level-set depth per column on the grid
    k_after = np.full(ncol, nz - 1)           # mask top at t_end
    T_pin = np.zeros(ncol)
    fired = np.zeros(ncol, bool)
    gaps_long = np.zeros(ncol)
    gaps_all = np.zeros(ncol)
    rflat = r.ravel()
    for c in range(ncol):
        e = ev[starts[c]:starts[c + 1]]
        if e.size == 0:
            continue
        cum = np.cumsum(e["h_applied"])
        idx = np.searchsorted(e["time"], tg, side="right") - 1
        depth_g[c] = np.where(idx >= 0, cum[np.maximum(idx, 0)], 0.0)
        k_after[c] = int(e["k_top"][-1] - e["n_voided"][-1])
        f = e[e["regime"] == 1]
        if f.size:
            fired[c] = True
            T_pin[c] = f["T_top"][-1]
            # idle proxy on W_end gaps
            sel = (f["time"][1:] >= W_end[0])
            if sel.sum():
                t_prev, Tp = f["time"][:-1][sel], f["T_top"][:-1][sel]
                kk = np.interp(t_prev, e["time"], e["k_top"] - e["n_voided"])
                s = z_n(t_prev) - (kk + 1) * dz - 0.0
                s = s - (m["z_top"] - nz * dz)          # z_face = PLO_z + (k+1) dz, PLO_z = z_top - nz dz
                h = jet.h_py(m["h_ref"], rflat[c], s, m["stag"])
                qp = h * (jet.Tg_py(m["T_ref"], s, m["T_ent"]) - Tp)
                tcell = np.where(qp > 0, RHOCP * dz * (Tp - TA) / np.where(qp > 0, qp, 1.0), 0.0)
                gap = f["time"][1:][sel] - t_prev
                gaps_long[c] = np.sum((qp <= 0) | (gap > m["idle"] * tcell))
                gaps_all[c] = sel.sum()

    # radial bins of 1 D
    nb = int(math.ceil(rflat.max() / D))
    b = np.minimum((rflat / D).astype(int), nb - 1)
    counts = np.bincount(b, minlength=nb)
    have = counts > 0
    bin_depth = np.zeros((nb, tg.size))
    for k in range(tg.size):
        bin_depth[:, k] = np.bincount(b, weights=depth_g[:, k], minlength=nb) / np.maximum(counts, 1)
    # C1 estimator: per-column least-squares slope of cumulative h_applied at
    # its event rows inside the window (>= 3 rows spanning >= half the window),
    # else the depth secant over the window; bin ROP = mean over the bin's columns.
    def col_rops(W):
        out = np.zeros(ncol)
        i0, i1 = int(round(W[0] / DT_GRID)), int(round(W[1] / DT_GRID))
        for c in range(ncol):
            e = ev[starts[c]:starts[c + 1]]
            sel = (e["time"] > W[0]) & (e["time"] <= W[1]) if e.size else np.zeros(0, bool)
            if sel.sum() >= 3 and np.ptp(e["time"][sel]) >= 0.5 * (W[1] - W[0]):
                out[c] = np.polyfit(e["time"][sel], np.cumsum(e["h_applied"])[sel], 1)[0]
            else:
                out[c] = (depth_g[c, i1] - depth_g[c, i0]) / (W[1] - W[0])
        return out * 3600.0
    def bin_mean(v):
        return np.where(have, np.bincount(b, weights=v, minlength=nb) / np.maximum(counts, 1), np.nan)
    rop_mid = bin_mean(col_rops(W_mid))
    rop_end = bin_mean(col_rops(W_end))
    rop_all = bin_mean(col_rops(W_all))
    rb = (np.arange(nb) + 0.5) * D
    d_end = bin_depth[:, -1]

    if feed > 0.0:
        ok = rop_end >= 0.9 * feed * 3600.0
        n_ok = int(np.argmin(ok)) if not ok.all() else nb
        R_h = n_ok * D
    else:
        c0 = d_end[0]
        R_h = float("nan")
        for i in range(1, nb):
            if have[i] and d_end[i] <= 0.5 * c0:
                R_h = float(np.interp(0.5 * c0, [d_end[i], d_end[i - 1]], [rb[i], rb[i - 1]]))
                break

    # centre column stand-off (mask) vs t
    c_axis = int(np.argmin(rflat))
    e = ev[starts[c_axis]:starts[c_axis + 1]]
    plo_z = m["z_top"] - nz * dz
    kg = np.full(tg.size, nz - 1.0)
    if e.size:
        idx = np.searchsorted(e["time"], tg, side="right") - 1
        kg = np.where(idx >= 0, (e["k_top"] - e["n_voided"])[np.maximum(idx, 0)], nz - 1)
    s_c = z_n(tg) - (plo_z + (kg + 1) * dz)
    s_c_ls = z_n(tg) - (m["z_top"] - depth_g[c_axis])
    ds_end = fit(tg, s_c_ls, *W_end)

    # temperatures
    R_T = max(R_h if np.isfinite(R_h) else 0.0, D)
    fsel = (ev["regime"] == 1) & (ev["time"] >= W_end[0]) & (rflat[col] <= R_T) if ev.size else np.zeros(0, bool)
    Tf = ev["T_top"][fsel] if ev.size else np.zeros(0)
    Tstats = (float(np.mean(Tf)), float(np.percentile(Tf, 10)), float(np.percentile(Tf, 90))) if Tf.size else (np.nan,) * 3

    # power, hole flux
    tw = (th["time"] >= W_end[0]) & (th["time"] <= W_end[1])
    mult = 4.0 if m["dim"] == 3 else 1.0
    P_face = float(np.mean(th["patch_P_robin"][tw])) * mult
    s_end = z_n(t_end) - (plo_z + (k_after + 1) * dz)
    h_end = jet.h_py(m["h_ref"], rflat, s_end, m["stag"])
    Tg_end = jet.Tg_py(m["T_ref"], s_end, m["T_ent"])
    q_pin = np.where(fired, h_end * (Tg_end - T_pin) - EPS * SIG * (T_pin ** 4 - TA ** 4) - HC * (T_pin - TA), 0.0)
    hole = rflat <= (R_h if np.isfinite(R_h) else 0.0)
    q_hole = float(np.mean(q_pin[hole])) if hole.any() else float("nan")

    # branch shares
    tot = th["pinned_cols_pinned"][tw] + th["pinned_cols_face"][tw]
    share = {k: float(np.sum(th[f"pinned_cols_{k}"][tw]) / max(np.sum(tot), 1.0))
             for k in ("pinned", "face_unset", "face_idle", "face_qneg", "face_minrule")}
    centre = rflat < 2.5 * D
    ring = (rflat >= max(R_h - 2.0 * D, 2.5 * D)) & (rflat <= R_h) if np.isfinite(R_h) else np.zeros_like(centre)
    idle_c = float(np.sum(gaps_long[centre]) / max(np.sum(gaps_all[centre]), 1.0))
    idle_r = float(np.sum(gaps_long[ring]) / max(np.sum(gaps_all[ring]), 1.0)) if ring.any() else float("nan")

    rim = int(np.sum(hole & ~fired))
    r_fired_max = float(rflat[fired].max()) if fired.any() else 0.0
    oor = (rflat > 7.5 * D) | (rflat < 2.5 * D) | (s_end < 2.0 * D) | (s_end > 12.0 * D)
    oor_r_hi = float(np.mean(rflat > 7.5 * D))
    oor_r_lo = float(np.mean(rflat < 2.5 * D))
    oor_s = float(np.mean((s_end < 2.0 * D) | (s_end > 12.0 * D)))
    in_hole_oor = float(np.mean(oor[hole])) if hole.any() else float("nan")
    slope = np.degrees(np.arctan(np.nanmax(np.abs(np.diff(d_end[have])) / D))) if have.sum() > 1 else float("nan")

    snaps = [t_end / 4.0, t_end / 2.0, 3.0 * t_end / 4.0, t_end]
    prof = {t: bin_depth[:, int(round(t / DT_GRID))] for t in snaps}
    return dict(m=m, nb=nb, rb=rb, have=have, rop_mid=rop_mid, rop_end=rop_end, rop_all=rop_all, prof_c=bin_depth[0], d_end=d_end, prof=prof,
                R_h=R_h, tg=tg, s_c=s_c, s_c_ls=s_c_ls, ds_end=ds_end, Tstats=Tstats, P_face=P_face,
                q_hole=q_hole, share=share, idle_c=idle_c, idle_r=idle_r, rim=rim, r_fired_max=r_fired_max,
                oor=(oor_r_hi, oor_r_lo, oor_s, in_hole_oor), slope=slope,
                ledger=float(np.max(np.abs(th["ledger_err"]))), pmin=th["patch_min_standoff"],
                tth=th["time"], n_fired=int(fired.sum()), ncol=ncol, depth_c=depth_g[c_axis])


def fmt(x, f="{:.2f}"):
    return "—" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f.format(x)


def run_table(res, part):
    L = []
    L.append("| run | h_ref [W/m²K] | T_ref [K] | T_ent | stag | idle | dz | dt [ms] | feed [m/h] | "
             "centre ROP W_mid / W_end (W_all) [m/h] | R_h [mm] | centre s end [mm] (ds/dt W_end mm/s) | "
             "T_fire mean (p10–p90) [K] | face P [kW] (% of 38) | hole flux [MW/m²] | ledger |")
    L.append("|" + "---|" * 16)
    for k, a in res.items():
        m = a["m"]
        pct = f"{a['P_face'] / 1e3:.2f} ({100 * a['P_face'] / P_BURNER:.1f} %)" if m["dim"] == 3 else f"{a['P_face']:.1f} W (slab)"
        L.append(f"| {k} | {m['h_ref']:.0f} | {m['T_ref']:.0f} | {m['T_ent']:.0f} | {m['stag']} | {m['idle']:g} | "
                 f"{m['dz']:g} | {m['dt'] * 1e3:g} | {m['feed'] * 3600:g} | "
                 f"{fmt(a['rop_mid'][0])} / {fmt(a['rop_end'][0])} ({fmt(a['rop_all'][0])}) | {fmt(a['R_h'] * 1e3, '{:.1f}')} | "
                 f"{a['s_c'][-1] * 1e3:.1f} ({fmt(a['ds_end'] * 1e3, '{:+.3f}')}) | "
                 f"{fmt(a['Tstats'][0], '{:.1f}')} ({fmt(a['Tstats'][1], '{:.0f}')}–{fmt(a['Tstats'][2], '{:.0f}')}) | "
                 f"{pct} | {fmt(a['q_hole'] / 1e6, '{:.3f}')} | {a['ledger']:.1e} |")
    return L


def check_table(res):
    L = ["| run | pinned | face: unset | idle | q ≤ 0 | min rule | idle-gap share centre / ring | "
         "rim never-fired (hole cols) | fired cols / all, max r fired [mm] | out of range: r>7.5D / r<2.5D / s∉[2,12]D / in hole | "
         "max face slope [°] | min patch s (thermo) [mm] |",
         "|" + "---|" * 12]
    for k, a in res.items():
        sh = a["share"]
        o = a["oor"]
        L.append(f"| {k} | {sh['pinned']:.3f} | {sh['face_unset']:.3f} | {sh['face_idle']:.3f} | {sh['face_qneg']:.3f} | "
                 f"{sh['face_minrule']:.3f} | {a['idle_c']:.3f} / {fmt(a['idle_r'], '{:.3f}')} | {a['rim']} | "
                 f"{a['n_fired']}/{a['ncol']}, {a['r_fired_max'] * 1e3:.1f} | "
                 f"{o[0]:.2f} / {o[1]:.2f} / {o[2]:.2f} / {fmt(o[3])} | {fmt(a['slope'], '{:.1f}')} | "
                 f"{np.min(a['pmin'][a['tth'] > 0]) * 1e3:.1f} |")
    return L


def profile_table(k, a, max_bins=None):
    nb = a["nb"] if max_bins is None else min(a["nb"], max_bins)
    idx = [i for i in range(nb) if a["have"][i]]
    step = max(1, len(idx) // 12)
    idx = idx[::step]
    L = [f"**{k}**: bin-mean level-set depth [mm] vs r (bin centres, 1 D bins), and W_end ROP [m/h]", "",
         "| t [s] | " + " | ".join(f"{a['rb'][i] / D:.1f} D" for i in idx) + " |",
         "|---|" + "---|" * len(idx)]
    for t, p in a["prof"].items():
        L.append(f"| {t:g} | " + " | ".join(f"{p[i] * 1e3:.1f}" for i in idx) + " |")
    L.append("| ROP W_end | " + " | ".join(fmt(a["rop_end"][i]) for i in idx) + " |")
    L.append("")
    return L


def closed_form_table():
    """Analytic pinned closed-form ROP at SOD = 7 D vs r (D2b feet question)."""
    rs = np.array([0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5]) * D
    L = ["| anchor | " + " | ".join(f"r = {x / D:.1f} D ({x * 1e3:.0f} mm)" for x in rs) + " | r where ROP_cf = 1.5 m/h [mm] |",
         "|---|" + "---|" * (len(rs) + 1)]
    for an, (h_ref, T_ref) in ANCHORS.items():
        s = 7.0 * D
        q = q_closed(jet.h_py(h_ref, rs, s), jet.Tg_py(T_ref, s))
        rop = q / (RHOCP * (jet.T_FIRE - T0)) * 3600.0
        rr = np.linspace(0.0, 0.3, 30001)
        qq = q_closed(jet.h_py(h_ref, rr, s), jet.Tg_py(T_ref, s)) / (RHOCP * (jet.T_FIRE - T0)) * 3600.0
        r15 = float(rr[np.argmax(qq < 1.5)]) if np.any(qq < 1.5) else float("nan")
        L.append(f"| {an} ({h_ref:.0f}, {T_ref:.0f} K) | " + " | ".join(f"{x:.2f}" for x in rop) +
                 f" | {fmt(r15 * 1e3, '{:.0f}')} |")
    return L


def plots(res, name):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return None
    n = len(res)
    fig, axs = plt.subplots(1, n, figsize=(3.2 * n, 3.2), squeeze=False)
    for ax, (k, a) in zip(axs[0], res.items()):
        for t, p in a["prof"].items():
            ax.plot(a["rb"][a["have"]] * 1e3, -p[a["have"]] * 1e3, label=f"{t:g} s")
        ax.set_title(k, fontsize=9)
        ax.set_xlabel("r [mm]")
        ax.grid(alpha=0.3)
    axs[0][0].set_ylabel("face z − z_top [mm]")
    axs[0][0].legend(fontsize=7)
    fig.tight_layout()
    out = os.path.join(HERE, f"{name}.png")
    fig.savefig(out, dpi=110)
    plt.close(fig)
    return os.path.basename(out)


def plot_pair(res, keys, name, title):
    """J-M vs J-10 side by side: face profiles at 4 times, centre depth and
    centre stand-off vs t, W_end ROP vs r, and the imposed h(r) / T_gas(s)."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return None
    keys = [k for k in keys if k in res]
    if len(keys) < 2:
        return None
    cols = {keys[0]: "tab:red", keys[1]: "tab:blue"}
    fig, axs = plt.subplots(2, 3, figsize=(13.5, 7.5))
    for j, k in enumerate(keys):
        a = res[k]
        ax = axs[0][j]
        for t, p in a["prof"].items():
            ax.plot(a["rb"][a["have"]] * 1e3, -p[a["have"]] * 1e3, label=f"{t:g} s")
        m = a["m"]
        ax.set_title(f"{k}: h_ref {m['h_ref']:.0f} W/m²K, T_ref {m['T_ref']:.0f} K", fontsize=10)
        ax.set_xlabel("r [mm]")
        ax.set_ylabel("face z − z_top [mm]")
        if np.isfinite(a["R_h"]):
            ax.axvline(a["R_h"] * 1e3, color="k", ls=":", lw=1, label=f"R_h {a['R_h'] * 1e3:.1f} mm")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    ax = axs[0][2]
    for k in keys:
        a = res[k]
        ax.plot(a["tg"], a["prof_c"] * 1e3, color=cols[k], label=f"{k} centre depth")
        ax.plot(a["tg"], a["s_c_ls"] * 1e3, color=cols[k], ls="--", label=f"{k} centre stand-off s")
    ax.set_xlabel("t [s]")
    ax.set_ylabel("[mm]")
    ax.set_title("centre bin depth and stand-off", fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    ax = axs[1][0]
    for k in keys:
        a = res[k]
        ax.plot(a["rb"][a["have"]] * 1e3, a["rop_end"][a["have"]], "o-", ms=3, color=cols[k], label=f"{k} W_end")
        ax.plot(a["rb"][a["have"]] * 1e3, a["rop_mid"][a["have"]], "o--", ms=3, color=cols[k], alpha=0.5,
                label=f"{k} W_mid")
        if a["m"]["feed"] > 0:
            ax.axhline(a["m"]["feed"] * 3600, color="k", ls=":", lw=1)
    ax.set_xlabel("r [mm]")
    ax.set_ylabel("ROP fit [m/h]")
    ax.set_title("window ROP vs r (removal_events fits)", fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    ax = axs[1][1]
    rr = np.linspace(0.0, 0.1, 400)
    for k in keys:
        m = res[k]["m"]
        for sD, ls in ((7.0, "-"), (12.0, "--")):
            ax.plot(rr * 1e3, jet.h_py(m["h_ref"], rr, sD * D, m["stag"]) / 1e3, ls, color=cols[k],
                    label=f"{k} s = {sD:g} D")
    ax.set_yscale("log")
    ax.set_xlabel("r [mm]")
    ax.set_ylabel("h [kW/m²K]")
    ax.set_title("imposed h(r, s)", fontsize=10)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8)
    ax = axs[1][2]
    ss = np.linspace(1.0 * D, 20.0 * D, 400)
    for k in keys:
        m = res[k]["m"]
        ax.plot(ss * 1e3, jet.Tg_py(m["T_ref"], ss, m["T_ent"]), color=cols[k], label=f"{k} T_gas(s)")
    ax.axhline(jet.T_FIRE, color="k", ls=":", lw=1, label="T_fire 821 K")
    ax.set_xlabel("s [mm]")
    ax.set_ylabel("T [K]")
    ax.set_title("imposed T_gas(s)", fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    out = os.path.join(HERE, f"{name}.png")
    fig.savefig(out, dpi=110)
    plt.close(fig)
    return os.path.basename(out)


def main():
    m0 = jet.martin_reference()
    res1, res2 = {}, {}
    for n in PART1:
        pf = os.path.join(OUT, f"p1_{n}")
        m = meta_of(pf)
        if m:
            res1[n] = analyze(pf, m)
    for n in PART2:
        pf = os.path.join(OUT, f"p2_{n}")
        m = meta_of(pf)
        if m:
            res2[n] = analyze(pf, m)

    L = ["# D2a: impinging-jet face source — results", "",
         "Generated by `analyze.py` (tables above the discussion marker regenerate). "
         "Definitions are in the `analyze.py` docstring. No parameter was adjusted toward Meier's ROP or hole diameter.", "",
         "## Reference arithmetic (jet.py)", "",
         f"- Martin (1977) at T_film = {m0['T_film']:.0f} K (air, Sutherland): μ = {m0['mu']:.3e} Pa·s, "
         f"k = {m0['k']:.4f} W/m·K, Pr = {m0['Pr']:.3f}; Re = 4ṁ/(πDμ) = {m0['Re']:.3e}.",
         f"- F(Re) = {m0['F']:.1f}, G(2.5 D, 7 D) = {m0['G']:.4f} → Nu_avg = {m0['Nu_avg']:.1f}, "
         f"h_avg(2.5 D) = {m0['h_avg']:.0f} W/m²K, **h_loc(2.5 D, 7 D) = {m0['h_loc']:.0f} W/m²K = J-M h_ref ({H_MARTIN})**.",
         f"- h_loc derivation check: area average reproduces G over 2.5–7.5 D (H/D = 2, 7, 12) to {jet.verify_hloc():.1e}.",
         "- Parser check (`run.py --parser-check`): compiled expressions vs numpy at s0 = 1.5–13 D, all anchors, "
         "both sensitivities — see README.", "",
         "Closed-form pinned ROP at SOD = 7 D (T_fire = 821 K, no conduction losses; analytic, information for D2b):", ""]
    L += closed_form_table() + [""]

    if res1:
        L += ["## Part 1 — 2-D slab machinery (input_2d_dev, feed 0, 60 s)", ""]
        L += run_table(res1, 1) + [""] + check_table(res1) + [""]
        for k, a in res1.items():
            L += profile_table(k, a, max_bins=11)
        png = plots(res1, "part1_profiles")
        if png:
            L += [f"![Part 1 profiles]({png})", ""]
        png = plot_pair(res1, ["JM", "J10"], "part1_JM_vs_J10",
                        "Part 1 (2-D slab, feed 0, 2 mm): J-M (Martin) vs J-10 (Meier Ch. 7)")
        if png:
            L += [f"![Part 1 J-M vs J-10]({png})", ""]
        # comparisons
        def pct(a, b):
            return 100.0 * (a / b - 1.0)

        def band(A, B):
            dmin = min(A["d_end"][0], B["d_end"][0])
            d1, d2 = 0.20 * dmin, 0.85 * dmin
            out = []
            for X in (A, B):
                dc = np.maximum.accumulate(X["prof_c"])
                t1, t2 = np.interp([d1, d2], dc + 1e-12 * np.arange(dc.size), X["tg"])
                out.append((d2 - d1) / (t2 - t1) * 3600.0)
            return out[0], out[1], d1, d2
        cmp = []
        if "JM" in res1 and "JM_1mm" in res1:
            A, B = res1["JM"], res1["JM_1mm"]
            ba, bb, d1, d2 = band(A, B)
            cmp.append(f"- **2 → 1 mm (J-M):** band ROP ({d1 * 1e3:.1f}→{d2 * 1e3:.1f} mm) {ba:.3f} → {bb:.3f} m/h "
                       f"({pct(bb, ba):+.1f} %); centre ROP fit W_all {A['rop_all'][0]:.3f} → {B['rop_all'][0]:.3f} m/h "
                       f"({pct(B['rop_all'][0], A['rop_all'][0]):+.1f} %; W_mid {pct(B['rop_mid'][0], A['rop_mid'][0]):+.1f} %, "
                       f"W_end {pct(B['rop_end'][0], A['rop_end'][0]):+.1f} %); T_fire {pct(B['Tstats'][0], A['Tstats'][0]):+.2f} % "
                       f"(5 % bar); end centre depth {A['d_end'][0] * 1e3:.1f} → {B['d_end'][0] * 1e3:.1f} mm "
                       f"({pct(B['d_end'][0], A['d_end'][0]):+.1f} %); R_h {A['R_h'] * 1e3:.1f} → {B['R_h'] * 1e3:.1f} mm.")
        if "JM" in res1 and "JM_dt1ms" in res1:
            A, B = res1["JM"], res1["JM_dt1ms"]
            ba, bb, d1, d2 = band(A, B)
            cmp.append(f"- **dt 4 → 1 ms (J-M):** band ROP {pct(bb, ba):+.2f} %; centre ROP fit W_all {pct(B['rop_all'][0], A['rop_all'][0]):+.2f} % "
                       f"(W_mid {pct(B['rop_mid'][0], A['rop_mid'][0]):+.2f} %, W_end {pct(B['rop_end'][0], A['rop_end'][0]):+.2f} %); "
                       f"end centre depth {pct(B['d_end'][0], A['d_end'][0]):+.2f} %; T_fire {pct(B['Tstats'][0], A['Tstats'][0]):+.2f} %; "
                       f"R_h {A['R_h'] * 1e3:.1f} → {B['R_h'] * 1e3:.1f} mm.")
        for s, lab in (("J10_idle1", "idle 2 → 1 (J-10)"), ("JM_Tent600", "T_ent 293 → 600 K (J-M)"),
                       ("JM_stag", "flat → rising stagnation (J-M)")):
            base = "J10" if s.startswith("J10") else "JM"
            if base in res1 and s in res1:
                A, B = res1[base], res1[s]
                ba, bb, d1, d2 = band(A, B)
                cmp.append(f"- **{lab}:** band ROP ({d1 * 1e3:.1f}→{d2 * 1e3:.1f} mm) {ba:.2f} → {bb:.2f} m/h "
                           f"({pct(bb, ba):+.1f} %); centre ROP fit W_all {pct(B['rop_all'][0], A['rop_all'][0]):+.1f} % "
                           f"(W_end {pct(B['rop_end'][0], A['rop_end'][0]):+.1f} %); end centre depth "
                           f"{A['d_end'][0] * 1e3:.1f} → {B['d_end'][0] * 1e3:.1f} mm; R_h {A['R_h'] * 1e3:.1f} → "
                           f"{B['R_h'] * 1e3:.1f} mm; idle share {A['share']['face_idle']:.3f} → {B['share']['face_idle']:.3f}; "
                           f"face P {A['P_face']:.1f} → {B['P_face']:.1f} W.")
        L += cmp + [""]

    if res2:
        L += ["## Part 2 — 3-D quarter domain (input_drilling, 60×60×100 at 2 mm, axis at the xlo/ylo corner)", ""]
        L += run_table(res2, 2) + [""] + check_table(res2) + [""]
        for k, a in res2.items():
            L += profile_table(k, a)
        png = plots(res2, "part2_profiles")
        if png:
            L += [f"![Part 2 profiles]({png})", ""]
        png = plot_pair(res2, ["JM_f15", "J10_f15"], "part2_JM_vs_J10",
                        "Part 2 (3-D quarter, feed 1.5 m/h, 2 mm): J-M (Martin) vs J-10 (Meier Ch. 7)")
        if png:
            L += [f"![Part 2 J-M vs J-10]({png})", ""]
        L += ["Centre stand-off s(t) (mask top) [mm]:", "",
              "| run | " + " | ".join(f"{t:g} s" for t in (10, 30, 60, 100, 150, 200, 250, 300)) + " |",
              "|---|" + "---|" * 8]
        for k, a in res2.items():
            vals = []
            for t in (10, 30, 60, 100, 150, 200, 250, 300):
                i = int(round(t / DT_GRID))
                vals.append(f"{a['s_c'][i] * 1e3:.1f}" if i < a["tg"].size else "—")
            L.append(f"| {k} | " + " | ".join(vals) + " |")
        L.append("")

    old = open(os.path.join(HERE, "RESULTS.md")).read() if os.path.exists(os.path.join(HERE, "RESULTS.md")) else ""
    disc = old[old.index(MARK):] if MARK in old else MARK + "\n\n## Discussion\n\n(to be written)\n"
    with open(os.path.join(HERE, "RESULTS.md"), "w") as f:
        f.write("\n".join(L) + "\n" + disc)
    print("wrote RESULTS.md;", len(res1), "part-1 runs,", len(res2), "part-2 runs")


if __name__ == "__main__":
    main()
