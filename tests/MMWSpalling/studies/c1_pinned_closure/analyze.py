#!/usr/bin/env python3
"""C1 analysis -> RESULTS.md (tables regenerate above <!-- DISCUSSION -->; the
hand-written discussion below the marker is preserved).

Part 1 (1-D columns): ROP from the level set at the 1 s plotfiles over the S1
window (max(first removal + 10 s, 0.3 t_end) -> t_end, snapped); <q_abs> =
window-mean ledger_P_robin / A; T_fire = removal_events.csv regime-1 T_top in
the window; <dT_rem> = frozen H/(rho Cp) of cells voided in the window;
closed form q_cf(T_fire) = h(T_gas - T_f) - eps sig (T_f^4 - T_a^4) - h_c(T_f - T_a),
ROP_cf = q_cf/(rho Cp (T_f - T_0)); pinned fraction = window thermo rows'
pinned_cols_pinned / (pinned + face).

Part 2 (Meier 2D slab): same window rule on the 2 s plotfiles; per-column
level-set recession; footprint = columns inside the 25 mm patch disk; centre =
the columns nearest x0; ring = footprint columns with |x - x0| > 20 mm.
Active-pin fraction per column group is reconstructed exactly from
removal_events.csv with BuildPinEffective's rule (pin set and
time - pin_t < idle * t_cell), sampled every 10 ms in the window; the pinned
BRANCH additionally needs T_face > T_pin, which only the global thermo
counts record.

ROP fit: least-squares slope of each column's cumulative h_applied
(removal_events.csv) over the window, sub-cell resolution. The plotfile ROP is
a whole-cell staircase sampled at 1-2 s and carries up to dz/(window) noise
(e.g. +-6 % for 2 mm over 20 s).
"""
import glob
import importlib.util
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
STUDIES = os.path.dirname(HERE)
S1 = os.path.join(STUDIES, "s1_surface_resolution")
sys.path.insert(0, S1)
sys.path.insert(0, HERE)
_spec = importlib.util.spec_from_file_location("s1_analyze", os.path.join(S1, "analyze.py"))
_s1 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_s1)
load_pf, load_thermo, pct = _s1.load_pf, _s1.load_thermo, _s1.pct
from run_sweep import DEPTH, dz_tag  # noqa: E402
from run import (PART1, PART2, T_GAS_1D, OUT, RHOCP, TA, q_closed,  # noqa: E402
                 part1_job, part2_job)

T0 = 293.15
MARK = "<!-- DISCUSSION -->"
P_BURNER = 38.0e3


def meta_of(pf):
    try:
        return eval(open(pf + ".done").read())
    except Exception:
        return None


def window(tp, t_first):
    t_end = tp[-1]
    t0 = max(t_first + 10.0, 0.3 * t_end)
    if t0 >= 0.8 * t_end:
        t0 = max(t_first + 2.0, 0.3 * t_end)
    i0 = int(np.searchsorted(tp, t0 - 1e-9))
    i1 = len(tp) - 1
    if i0 >= i1:
        i0 = max(0, i1 - 1)
    im = i0 + int(np.argmin(np.abs(tp[i0:i1 + 1] - 0.5 * (tp[i0] + tp[i1]))))
    return i0, im, i1


def fit_rop(ev, t0, t1):
    """Level-set recession rate [m/h] of one column from removal_events.csv:
    least-squares slope of cumulative h_applied vs time over rows in (t0, t1].
    Sub-cell resolution (the phi shift is continuous), unlike the whole-cell
    plotfile staircase."""
    if ev.size == 0:
        return 0.0
    cum = np.cumsum(ev["h_applied"])
    sel = (ev["time"] > t0) & (ev["time"] <= t1)
    if sel.sum() < 3:
        return 0.0
    return float(np.polyfit(ev["time"][sel], cum[sel], 1)[0]) * 3600.0


def dlnrop(T_gas, T_f):
    return -1.0 / (T_gas - T_f) - 1.0 / (T_f - T0)


# ------------------------------------------------------------------ Part 1
def analyze_p1(name, dz_mm):
    pf = os.path.join(OUT, f"p1_{name}_{dz_tag(dz_mm)}")
    m = meta_of(pf)
    r = dict(name=name, dz=dz_mm)
    if m is None:
        r["status"] = "missing"
        return r
    dz = dz_mm * 1e-3
    N = int(round(DEPTH / dz))
    zc = (np.arange(N) + 0.5) * dz
    pfs = sorted(glob.glob(os.path.join(pf, "*cell")), key=lambda p: int(re.search(r"(\d+)cell$", p).group(1)))
    pd = [load_pf(p) for p in pfs]
    tp = np.array([d["t"] for d in pd])

    def zs(d):
        k = np.nonzero(~d["rem"])[0].max()
        return zc[k] + d["phi"][k]
    z = np.array([zs(d) for d in pd])
    ev = np.genfromtxt(pf + "_removal_events.csv", delimiter=",", names=True, ndmin=1)
    fire = ev[ev["regime"] == 1]
    i0, im, i1 = window(tp, float(fire["time"][0]))
    t0, tm, t1 = tp[i0], tp[im], tp[i1]
    th = load_thermo(os.path.join(pf, "thermo.dat"))
    tt = th["time"]
    sw = (tt > t0) & (tt <= t1)
    fw = fire[(fire["time"] > t0) & (fire["time"] <= t1)]
    w = pd[-1]["rem"] & ~pd[i0]["rem"]
    T_gas = T_GAS_1D[PART1[name][0]]
    Tf = float(np.mean(fw["T_top"]))
    q_cf = q_closed(1.0e4, T_gas, Tf)
    npin, nface = th["pinned_cols_pinned"][sw].sum(), th["pinned_cols_face"][sw].sum()
    r.update(status="ok", dt=m["dt"], overshoot=m["overshoot"], t0=t0, t1=t1,
             depth=(DEPTH - z[-1]) * 1e3,
             rop=(z[i0] - z[i1]) / (t1 - t0) * 3600, rop_a=(z[i0] - z[im]) / (tm - t0) * 3600,
             rop_b=(z[im] - z[i1]) / (t1 - tm) * 3600,
             q_abs=float(np.mean(th["ledger_P_robin"][sw])) / (dz * dz),
             Tf=Tf, Tf10=pct(fw["T_top"], 10), Tf90=pct(fw["T_top"], 90),
             dT_rem=float(np.mean(pd[-1]["H"][w]) / RHOCP),
             q_cf=q_cf, rop_cf=q_cf / (RHOCP * (Tf - T0)) * 3600,
             pin_frac=npin / max(npin + nface, 1.0),
             ledger=float(np.max(np.abs(th["ledger_err"]))), T_gas=T_gas,
             rop_fit=fit_rop(ev, t0, t1), rop_fit_a=fit_rop(ev, t0, 0.5 * (t0 + t1)),
             rop_fit_b=fit_rop(ev, 0.5 * (t0 + t1), t1))
    return r


# ------------------------------------------------------------------ Part 2
def load_pf3(p):
    import yt
    yt.funcs.mylog.setLevel(50)
    ds = yt.load(p)
    cg = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    return dict(t=float(ds.current_time), rem=np.asarray(cg["removed"]) > 0.5,
                phi=np.asarray(cg["phi"]), H=np.asarray(cg["H"]),
                dx=np.array((ds.domain_width / ds.domain_dimensions).to_value(), float))


def analyze_p2(name, dz_mm):
    pf = os.path.join(OUT, f"p2_{name}_{dz_tag(dz_mm)}")
    m = meta_of(pf)
    r = dict(name=name, dz=dz_mm)
    if m is None:
        r["status"] = "missing"
        return r
    h, T_gas = m["h"], m["T_gas"]
    pfs = sorted(glob.glob(os.path.join(pf, "*cell")), key=lambda p: int(re.search(r"(\d+)cell$", p).group(1)))
    pd = [load_pf3(p) for p in pfs]
    tp = np.array([d["t"] for d in pd])
    dx = pd[0]["dx"]
    nx, ny, nz = pd[0]["rem"].shape
    x = (np.arange(nx) + 0.5) * dx[0]
    y = (np.arange(ny) + 0.5) * dx[1]
    zc = (np.arange(nz) + 0.5) * dx[2]
    ztop = nz * dx[2]
    X, Y = np.meshgrid(x, y, indexing="ij")
    foot = (X - 0.07) ** 2 + (Y - 0.004) ** 2 <= 0.025 ** 2
    xd = np.abs(x - 0.07)
    centre_i = np.nonzero(xd <= xd.min() + 1e-12)[0]
    ring = foot & (np.abs(X - 0.07) > 0.020)
    centre = np.zeros_like(foot)
    centre[centre_i, :] = True

    def surf(d):
        solid = ~d["rem"]
        k = nz - 1 - np.argmax(solid[:, :, ::-1], axis=2)
        ii, jj = np.indices(k.shape)
        return zc[k] + d["phi"][ii, jj, k]
    zs = np.array([surf(d) for d in pd])                 # (nt, nx, ny)
    ev = np.genfromtxt(pf + "_removal_events.csv", delimiter=",", names=True, ndmin=1)
    fire = ev[ev["regime"] == 1]
    i0, im, i1 = window(tp, float(fire["time"][0]))
    t0, tm, t1 = tp[i0], tp[im], tp[i1]
    rate = lambda a, b: (zs[a] - zs[b]) / (tp[b] - tp[a]) * 3600.0
    V = dx.prod()
    Ly = ny * dx[1]
    vol = np.array([d["rem"].sum() * V for d in pd])
    fw = fire[(fire["time"] > t0) & (fire["time"] <= t1)]
    th = load_thermo(os.path.join(pf, "thermo.dat"))
    tt = th["time"]
    sw = (tt > t0) & (tt <= t1)
    npin, nface = th["pinned_cols_pinned"][sw].sum(), th["pinned_cols_face"][sw].sum()
    q_mean = float(np.mean(th["pinned_q_mean"][sw]))

    # Active-pin reconstruction (BuildPinEffective rule) per column group.
    dz = dx[2]
    idle = 2.0
    samples = np.arange(t0, t1, 0.01)
    act = {"centre": [0, 0], "ring": [0, 0]}
    cols = {}
    for row in fire:
        cols.setdefault((int(row["col_i"]), int(row["col_j"])), []).append((row["time"], row["T_top"]))
    for (ci, cj), lst in cols.items():
        grp = "centre" if centre[ci, cj] else ("ring" if ring[ci, cj] else None)
        if grp is None:
            continue
        ft = np.array([a for a, _ in lst])
        fT = np.array([b for _, b in lst])
        k = np.searchsorted(ft, samples, side="right") - 1
        ok = k >= 0
        Tp = np.where(ok, fT[np.maximum(k, 0)], 0.0)
        tpn = np.where(ok, ft[np.maximum(k, 0)], -1e30)
        qp = h * (T_gas - Tp)
        tcell = np.where(qp > 0, RHOCP * dz * (Tp - TA) / np.where(qp > 0, qp, 1.0), 0.0)
        active = ok & (qp > 0) & (samples - tpn < idle * tcell)
        act[grp][0] += int(active.sum())
        act[grp][1] += samples.size
    for grp in ("centre", "ring"):
        n_ring_cols = int((ring if grp == "ring" else centre).sum())
        act[grp][1] = max(act[grp][1], 1)
    # columns that never fired count as inactive for all samples
    n_centre, n_ring = int(centre.sum()), int(ring.sum())
    frac_c = act["centre"][0] / max(samples.size * n_centre, 1)
    frac_r = act["ring"][0] / max(samples.size * n_ring, 1)

    tmid = 0.5 * (t0 + t1)
    fits = {}
    for (ci, cj) in zip(*np.nonzero(foot)):
        sub = ev[(ev["col_i"] == ci) & (ev["col_j"] == cj)]
        fits[(ci, cj)] = (fit_rop(sub, t0, t1), fit_rop(sub, t0, tmid), fit_rop(sub, tmid, t1))
    fc = [fits[k] for k in fits if centre[k]]
    ff = list(fits.values())
    never = int(sum(1 for k in fits if not np.any((ev["col_i"] == k[0]) & (ev["col_j"] == k[1]))))
    Tf = float(np.mean(fw["T_top"]))
    q_cf = q_closed(h, T_gas, Tf)
    prof_x = x[(x >= 0.07 - 0.035) & (x <= 0.07 + 0.035)]
    depth_end = (ztop - zs[-1].mean(axis=1)) * 1e3
    r.update(status="ok", dt=m["dt"], overshoot=m["overshoot"], t0=t0, t1=t1, stop=m["stop"],
             vrate=(vol[i1] - vol[i0]) / (t1 - t0) / Ly * 1e6,           # mm^2/s per unit slab depth
             rop_c=float(rate(i0, i1)[centre].mean()), rop_c_a=float(rate(i0, im)[centre].mean()),
             rop_c_b=float(rate(im, i1)[centre].mean()),
             rop_f=float(rate(i0, i1)[foot].mean()), rop_f_a=float(rate(i0, im)[foot].mean()),
             rop_f_b=float(rate(im, i1)[foot].mean()),
             depth_c=float((ztop - zs[-1][centre]).mean() * 1e3),
             Tf=Tf, Tf10=pct(fw["T_top"], 10), Tf90=pct(fw["T_top"], 90), n_fire=int(fw.size),
             pin_frac=npin / max(npin + nface, 1.0), act_c=frac_c, act_r=frac_r,
             q_mean=q_mean, P=q_mean * math.pi * 0.025 ** 2,
             q_cf=q_cf, rop_cf=q_cf / (RHOCP * (Tf - T0)) * 3600,
             ledger=float(np.max(np.abs(th["ledger_err"]))), h=h, T_gas=T_gas,
             fit_c=[float(np.mean([f[i] for f in fc])) for i in range(3)],
             fit_f=[float(np.mean([f[i] for f in ff])) for i in range(3)],
             n_foot=len(ff), never=never,
             profile=[(float(xx), float(dd)) for xx, dd in zip(x, depth_end) if abs(xx - 0.07) <= 0.035])
    return r


def main():
    L = ["# C1: pinned-surface Robin closure (`robin_form = pinned`)\n",
         "Generated by `analyze.py` (see its docstring for every definition). "
         "The closed form is q_cf = h(T_gas - T_f) - eps sig (T_f^4 - T_a^4) - h_c (T_f - T_a) "
         "at the run's own mean T_fire, ROP_cf = q_cf / (rho Cp (T_f - 293.15 K)).\n",
         "## Part 1: 1-D columns vs S1\n",
         "S1 column (1 x 1 x N, 250 mm), degenerate a = 20 um, confining 1 MPa, losses on "
         "(eps 0.8, h_c 10), follow_mask + ledger, `pinned_idle_cycles` = 2 unless noted, h = 1e4 W/m^2K. "
         "Window max(first removal + 10 s, 0.3 t_end) -> t_end on 1 s plotfiles.\n",
         "| case | dz mm | dt ms (K/step) | window s | depth mm | ROP m/h (1st / 2nd) | <q_abs> MW/m^2 | "
         "ROP fit m/h (1st / 2nd) | T_fire K mean [p10, p90] | <dT_rem> K | q_cf MW/m^2 | ROP_cf m/h | ROP fit/ROP_cf | pinned frac | "
         "d ln ROP/dT_f %/K | max abs ledger_err |",
         "|" + "---|" * 16]
    p1 = []
    for name, (_, dzs, _) in PART1.items():
        for dz in dzs:
            r = analyze_p1(name, dz)
            p1.append(r)
            if r["status"] != "ok":
                L.append(f"| {name} | {dz:g} | {r['status']} |" + " |" * 13)
                continue
            print(f"P1 {name} {dz}: ROP {r['rop']:.2f} cf {r['rop_cf']:.2f} q {r['q_abs']/1e6:.3f} Tf {r['Tf']:.1f}", flush=True)
            L.append(f"| {name} | {dz:g} | {r['dt']*1e3:g} ({r['overshoot']:.2f}) | {r['t0']:.0f}-{r['t1']:.0f} | "
                     f"{r['depth']:.1f} | {r['rop']:.2f} ({r['rop_a']:.2f} / {r['rop_b']:.2f}) | "
                     f"{r['q_abs']/1e6:.3f} | {r['rop_fit']:.2f} ({r['rop_fit_a']:.2f} / {r['rop_fit_b']:.2f}) | "
                     f"{r['Tf']:.1f} [{r['Tf10']:.1f}, {r['Tf90']:.1f}] | {r['dT_rem']:.1f} | "
                     f"{r['q_cf']/1e6:.3f} | {r['rop_cf']:.2f} | {r['rop_fit']/r['rop_cf']:.3f} | {r['pin_frac']:.3f} | "
                     f"{dlnrop(r['T_gas'], r['Tf'])*100:.3f} | {r['ledger']:.1e} |")
    seq = sorted([r for r in p1 if r["name"] == "P-1000" and r["status"] == "ok"], key=lambda r: -r["dz"])
    if seq:
        L.append("\nP-1000 change per halving: " + "; ".join(
            f"{a['dz']:g} -> {b['dz']:g} mm: ROP fit {(b['rop_fit']-a['rop_fit'])/a['rop_fit']*100:+.2f}% "
            f"(plotfile {(b['rop']-a['rop'])/a['rop']*100:+.2f}%), "
            f"q {(b['q_abs']-a['q_abs'])/a['q_abs']*100:+.2f}%" for a, b in zip(seq, seq[1:])))
    base = next((r for r in p1 if r["name"] == "P-1000" and r["dz"] == 2.0 and r["status"] == "ok"), None)
    if base:
        parts = []
        for nm in ("P-1000-idle1", "P-1000-idle4"):
            rr = next((r for r in p1 if r["name"] == nm and r["status"] == "ok"), None)
            if rr:
                parts.append(f"{nm}: ROP fit {(rr['rop_fit']-base['rop_fit'])/base['rop_fit']*100:+.2f}% vs idle 2 "
                             f"(q {(rr['q_abs']-base['q_abs'])/base['q_abs']*100:+.2f}%)")
        L.append("\nIdle-cycles sensitivity at 2 mm: " + "; ".join(parts))

    L += ["\n## Part 2: 2-D machinery check on the Meier dev harness\n",
          "`input_2d_dev` by CLI: beam off; pinned Robin patch x0 = 0.07, y0 = 0.004, radius 25 mm; "
          "follow_mask + ledger + removal_events; `weibull.V0` = V_cell; `pinned_idle_cycles` = 2; 4 ranks. "
          "Hole diameter = footprint by construction; this is not a Meier ROP comparison.\n",
          "| case | dz mm | dt ms (K/step) | stop s | window s | vol. rate / slab depth mm^2/s | centre ROP m/h (1st / 2nd) | "
          "footprint ROP m/h (1st / 2nd) | centre ROP fit m/h (1st / 2nd) | footprint ROP fit m/h (1st / 2nd) | "
          "never-fired footprint cols | centre depth mm | T_fire K mean [p10, p90] (n) | q_cf MW/m^2 | ROP_cf m/h | "
          "pinned frac (thermo) | active pin centre / ring | <q> MW/m^2 | P kW (% of 38 kW) | max abs ledger_err |",
          "|" + "---|" * 20]
    p2 = []
    for name, (_, _, dzs) in PART2.items():
        for dz in dzs:
            r = analyze_p2(name, dz)
            p2.append(r)
            if r["status"] != "ok":
                L.append(f"| {name} | {dz:g} | {r['status']} |" + " |" * 17)
                continue
            print(f"P2 {name} {dz}: centre ROP {r['rop_c']:.2f} foot {r['rop_f']:.2f} cf {r['rop_cf']:.2f} Tf {r['Tf']:.1f}", flush=True)
            L.append(f"| {name} | {dz:g} | {r['dt']*1e3:g} ({r['overshoot']:.2f}) | {r['stop']:g} | {r['t0']:.0f}-{r['t1']:.0f} | "
                     f"{r['vrate']:.2f} | {r['rop_c']:.2f} ({r['rop_c_a']:.2f} / {r['rop_c_b']:.2f}) | "
                     f"{r['rop_f']:.2f} ({r['rop_f_a']:.2f} / {r['rop_f_b']:.2f}) | "
                     f"{r['fit_c'][0]:.2f} ({r['fit_c'][1]:.2f} / {r['fit_c'][2]:.2f}) | "
                     f"{r['fit_f'][0]:.2f} ({r['fit_f'][1]:.2f} / {r['fit_f'][2]:.2f}) | {r['never']}/{r['n_foot']} | "
                     f"{r['depth_c']:.1f} | "
                     f"{r['Tf']:.1f} [{r['Tf10']:.1f}, {r['Tf90']:.1f}] ({r['n_fire']}) | {r['q_cf']/1e6:.3f} | "
                     f"{r['rop_cf']:.2f} | {r['pin_frac']:.3f} | {r['act_c']:.3f} / {r['act_r']:.3f} | "
                     f"{r['q_mean']/1e6:.3f} | {r['P']/1e3:.2f} ({r['P']/P_BURNER*100:.1f}%) | {r['ledger']:.1e} |")
    L.append("\nEnd-of-run depth profile (mm below the initial top, mean over y) vs x - x0 (mm):\n")
    for r in p2:
        if r["status"] != "ok":
            continue
        step = max(1, int(round(4.0 / (r["dz"]))))
        pts = r["profile"][::step]
        L.append(f"- {r['name']} dz {r['dz']:g}: " + ", ".join(f"{(xx-0.07)*1e3:+.0f}: {d:.1f}" for xx, d in pts))
    a = next((r for r in p2 if r["name"] == "M-10k-1000" and r["dz"] == 2.0 and r["status"] == "ok"), None)
    b = next((r for r in p2 if r["name"] == "M-10k-1000" and r["dz"] == 1.0 and r["status"] == "ok"), None)
    if a and b:
        L.append(f"\nM-10k-1000 mesh check 2 -> 1 mm: T_fire {a['Tf']:.1f} -> {b['Tf']:.1f} K "
                 f"({(b['Tf']-T0-(a['Tf']-T0))/(a['Tf']-T0)*100:+.2f}% of dT_fire); centre ROP "
                 f"{a['rop_c']:.2f} -> {b['rop_c']:.2f} m/h ({(b['rop_c']-a['rop_c'])/a['rop_c']*100:+.2f}%; fit "
                 f"{a['fit_c'][0]:.2f} -> {b['fit_c'][0]:.2f} m/h, {(b['fit_c'][0]-a['fit_c'][0])/a['fit_c'][0]*100:+.2f}%); "
                 f"footprint ROP {a['rop_f']:.2f} -> {b['rop_f']:.2f} m/h ({(b['rop_f']-a['rop_f'])/a['rop_f']*100:+.2f}%; fit "
                 f"{a['fit_f'][0]:.2f} -> {b['fit_f'][0]:.2f} m/h, {(b['fit_f'][0]-a['fit_f'][0])/a['fit_f'][0]*100:+.2f}%)")

    text = "\n".join(L) + "\n"
    res = os.path.join(HERE, "RESULTS.md")
    disc = f"\n{MARK}\n\n## Discussion\n\n(to be written)\n"
    if os.path.exists(res) and MARK in open(res).read():
        disc = "\n" + MARK + open(res).read().split(MARK, 1)[1]
    open(res, "w").write(text + disc)
    print("wrote", res)


if __name__ == "__main__":
    main()
