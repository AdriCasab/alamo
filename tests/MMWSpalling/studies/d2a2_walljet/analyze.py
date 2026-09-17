#!/usr/bin/env python3
"""D2a2 analysis -> RESULTS.md and PNGs.

RESULTS.md layout: the pre-registration (above <!-- RESULTS -->) and the
hand-written discussion (below <!-- DISCUSSION -->) are preserved; the tables
between the two markers are regenerated.

Recession comes from removal_events.csv only (C1 takeaway 8): a column's depth
below the original top at time t is the cumulative h_applied of its rows with
time <= t. Definitions:
  window W        [t_end/2, t_end] (t_end = the run's stop time).
  column ROP      least-squares slope of cumulative h_applied at the column's
                  event rows in W (>= 3 rows spanning >= W/2), else the depth
                  secant over W (the C1/D2a estimator).
  ring ROP        mean column ROP over a radial ring of 1 D (r of the column
                  centre from the jet axis); "ROP @40" is the ring holding 40 mm
                  (37.5-45 mm).
  r_wall(d)       max r of a column whose depth at t_end exceeds d (the packet's
                  wall profile, d = depth below the original top); hole
                  diameter = 2 r_wall. In the quarter domain r >= 120 mm is
                  sampled only off-axis: flagged "dom" (domain-limited).
  r_fire          max r of a column that ever fired (regime 1).
  centre          the column nearest the axis: depth, ROP in W, stand-off
                  s = z_n(t_end) - mask top face at t_end.
  powers          thermo at the end and W means; quarter-domain values x4
                  (full jet); jet_P_face against jet_P_cap and against Meier's
                  2.83 kW rock-side removal power. floor = patch_P_robin -
                  jet_P_face (power not drawn from the nozzle stream: T_ent
                  floor supply or negative q on hot columns).
  T_fire          regime-1 T_top of all columns in W: mean, p10, p90.
  branch shares   W sums of thermo pinned_cols_* / (pinned + face).
  rim             columns with r < r_wall(2 mm) that never fired.
  max slope       max |d depth/dr| between adjacent 1 D rings at t_end, as an angle.
  out of range    share of in-hole columns (depth > 2 mm) with r outside
                  [2.5, 7.5] D or s outside [2, 12] D at t_end (Martin validity).
  closed form @40 pinned q = h(0.04, s_ring)(T_gas(40 mm) - 821) - losses and its
                  ROP, T_gas from the last profile row, s_ring = mean ring s at t_end.
  steady          centre ROP in [t_end/2, 3t_end/4] vs [3t_end/4, t_end] and
                  r_wall(10 mm) at 3t_end/4 vs t_end; "centre s steady from" =
                  the first time after which the centre stand-off stays within
                  +-3 mm of its W mean.
  hole Ø (ref)    2 x the outer edge of the contiguous run of 1 D rings from
                  the axis whose W ROP >= 0.9 x 1.5 m/h (the reference table's
                  removal criterion applied to the simulation).
  reference       walljet.ref_table (flat face at SOD, with losses) for the
                  run's anchor and treatment; the packet's table is in section 0.
  nozzle budget   jet_mdot cp (T_nozzle - T_ent) = jet_P_face + jet_P_exhaust +
                  jet_P_decay (x4 in the quarter domain).
  early window E  Part 2: [40, 110] s (or [t_end/2, t_end] for runs shorter than
                  60 s), before the r = 40 mm ring rises above the descending
                  nozzle plane; t_flank = first 5 s grid time at which that
                  ring's mean face is at or above z_nozzle (after which the
                  no-flux rule switches it off). Ring ROP @40 in E is the
                  number compared with the reference table.
T_gas(r) comes from <plot_file>_jet_profile.csv (bin inlet values).
"""
import ast
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import walljet as wj  # noqa: E402
from run import ANCHORS, OUT, PART1, PART2, part1_job, part2_job  # noqa: E402

D = wj.D
MARK_R, MARK_D = "<!-- RESULTS -->", "<!-- DISCUSSION -->"
DEPTHS = (0.002, 0.010, 0.020, 0.040)


KILLED = {"p2_J10_19_ent": "stopped by user decision (cost)"}


def meta_of(pf):
    """.done metadata; for a run that ended early (nozzle-collision abort or
    stopped by hand) the job metadata with stop = its last thermo time and the
    reason in meta['ended']."""
    try:
        return ast.literal_eval(open(pf + ".done").read())
    except Exception:
        pass
    name = os.path.basename(pf)
    if not os.path.exists(os.path.join(pf, "thermo.dat")):
        return None
    try:
        part, n = int(name[1]), name[3:]
        m = (part1_job if part == 1 else part2_job)(n)[2]
    except Exception:
        return None
    rows = [l.split() for l in open(os.path.join(pf, "thermo.dat")) if l.strip()]
    t_last = float(rows[-1][0])
    log = open(pf + ".log", errors="replace").read()
    if "nozzle collision" in log:
        i = log.index("nozzle collision")
        m["ended"] = "nozzle collision at column " + log[i:].split("(i, j) = ")[1].split(")")[0] + ")"
    elif name in KILLED:
        m["ended"] = KILLED[name]
    else:
        return None
    m["planned_stop"] = m["stop"]
    m["stop"] = t_last
    return m


def thermo(pf):
    lines = [l.split() for l in open(os.path.join(pf, "thermo.dat")) if l.strip()]
    a = np.array(lines[1:], dtype=float)
    return {n: a[:, i] for i, n in enumerate(lines[0])}


def csv(path):
    return np.genfromtxt(path, delimiter=",", names=True, ndmin=1)


def geometry(m):
    dz = m["dz"] * 1e-3
    if m["dim"] == 2:
        nx, ny, nz = 70, 4, 60
    else:
        nx, ny, nz = 60, 60, 100
    xc = (np.arange(nx) + 0.5) * dz
    yc = (np.arange(ny) + 0.5) * dz
    X, Y = np.meshgrid(xc, yc, indexing="xy")          # [j, i]
    return dz, nx, ny, nz, np.hypot(X - m["x0"], Y - m["y0"]).ravel()


def analyze(pf, m, t_c=None):
    dz, nx, ny, nz, r = geometry(m)
    ncol = nx * ny
    th, ev = thermo(pf), csv(pf + "_removal_events.csv")
    prof = csv(pf + "_jet_profile.csv")
    t_end = float(m["stop"])
    W = (0.5 * t_end, t_end)
    scale = 1.0 / m["sector"] if m["sector"] else 1.0
    z_n = lambda t: m["z_top"] + wj.SOD - m["feed"] * t
    plo_z = m["z_top"] - nz * dz

    col = ev["col_j"].astype(int) * nx + ev["col_i"].astype(int)
    order = np.lexsort((ev["time"], col))
    ev, col = ev[order], col[order]
    st = np.searchsorted(col, np.arange(ncol + 1))

    def depth_at(c, t):
        e = ev[st[c]:st[c + 1]]
        if e.size == 0:
            return 0.0
        i = np.searchsorted(e["time"], t, side="right")
        return float(np.sum(e["h_applied"][:i]))

    def col_rop(c, w):
        e = ev[st[c]:st[c + 1]]
        sel = (e["time"] > w[0]) & (e["time"] <= w[1])
        if sel.sum() >= 3 and np.ptp(e["time"][sel]) >= 0.5 * (w[1] - w[0]):
            return float(np.polyfit(e["time"][sel], np.cumsum(e["h_applied"])[sel], 1)[0]) * 3600.0
        return (depth_at(c, w[1]) - depth_at(c, w[0])) / (w[1] - w[0]) * 3600.0

    depth = np.array([depth_at(c, t_end) for c in range(ncol)])
    depth34 = np.array([depth_at(c, 0.75 * t_end) for c in range(ncol)])
    ktop = np.full(ncol, nz - 1)
    fired = np.zeros(ncol, bool)
    for c in range(ncol):
        e = ev[st[c]:st[c + 1]]
        if e.size:
            ktop[c] = int(e["k_top"][-1] - e["n_voided"][-1])
            fired[c] = bool(np.any(e["regime"] == 1))
    s_end = z_n(t_end) - (plo_z + (ktop + 1) * dz)
    rop = np.array([col_rop(c, W) for c in range(ncol)])

    nb = int(math.ceil(r.max() / D))
    b = np.minimum((r / D).astype(int), nb - 1)
    cnt = np.bincount(b, minlength=nb)
    have = cnt > 0
    ring = lambda v: np.where(have, np.bincount(b, weights=v, minlength=nb) / np.maximum(cnt, 1), np.nan)
    ring_rop, ring_depth = ring(rop), ring(depth)
    rb = (np.arange(nb) + 0.5) * D
    i40 = int(0.040 / D)

    def r_wall(dep, d):
        sel = dep > d
        return float(r[sel].max()) if sel.any() else 0.0

    rw = {d: r_wall(depth, d) for d in DEPTHS}
    # common-time snapshot (runs stop at different times)
    common = None
    if t_c is not None:
        dc = np.array([depth_at(c, t_c) for c in range(ncol)])
        fc = np.zeros(ncol, bool)
        for c in range(ncol):
            e = ev[st[c]:st[c + 1]]
            fc[c] = bool(np.any((e["regime"] == 1) & (e["time"] <= t_c))) if e.size else False
        ic = int(np.searchsorted(th["time"], t_c))
        ic = min(ic, th["time"].size - 1)
        common = dict(t=t_c, rw2=r_wall(dc, 0.002), rw10=r_wall(dc, 0.010),
                      r_fire=float(r[fc].max()) if fc.any() else 0.0,
                      Pf=th["jet_P_face"][ic] * scale, Tex=th["jet_T_exhaust"][ic],
                      c_depth=dc[int(np.argmin(r))])
    rw10_34 = r_wall(depth34, 0.010)
    ca = int(np.argmin(r))
    # centre stand-off series (mask top after each event)
    tg = np.arange(0.0, t_end + 1e-9, 1.0)
    e = ev[st[ca]:st[ca + 1]]
    kk = np.full(tg.size, nz - 1.0)
    if e.size:
        i = np.searchsorted(e["time"], tg, side="right") - 1
        kk = np.where(i >= 0, (e["k_top"] - e["n_voided"])[np.maximum(i, 0)], nz - 1.0)
    s_c_t = z_n(tg) - (plo_z + (kk + 1) * dz)
    sW = s_c_t[tg >= W[0]].mean()
    off = np.abs(s_c_t - sW) > 0.003
    steady_from = float(tg[np.nonzero(off)[0][-1] + 1]) if off.any() and np.nonzero(off)[0][-1] + 1 < tg.size \
        else (0.0 if not off.any() else float("nan"))
    c_rop1 = col_rop(ca, (0.5 * t_end, 0.75 * t_end))
    c_rop2 = col_rop(ca, (0.75 * t_end, t_end))

    # early window E and the flank-crossing time of the 40 mm ring
    E = (40.0, min(110.0, t_end)) if t_end > 60.0 else (0.5 * t_end, t_end)
    ringc = np.nonzero(b == i40)[0]
    t_flank = None
    for tt in np.arange(0.0, t_end + 1e-9, 5.0):
        face = m["z_top"] - np.mean([depth_at(c, tt) for c in ringc])
        if z_n(tt) - face <= 0.0:
            t_flank = float(tt)
            break
    rop40_E = float(np.mean([(depth_at(c, E[1]) - depth_at(c, E[0])) / (E[1] - E[0]) * 3600.0
                             for c in ringc]))
    cropE = (depth_at(ca, E[1]) - depth_at(ca, E[0])) / (E[1] - E[0]) * 3600.0
    dE = np.array([depth_at(c, E[1]) for c in range(ncol)])
    rw10_E = r_wall(dE, 0.010)
    inE = (th["time"] > E[0]) & (th["time"] <= E[1])
    ptimes = np.unique(prof["time"])
    tpE = ptimes[np.argmin(np.abs(ptimes - E[1]))]
    rr, TT = prof["r_lo"][prof["time"] == tpE] + 0.5 * m["jet_dr"], prof["T_gas"][prof["time"] == tpE]
    T40_E = float(np.interp(0.040, rr, TT, right=np.nan))
    early = dict(E=E, t_flank=t_flank, rop40=rop40_E, c_rop=cropE, rw10=rw10_E,
                 Tstag=th["jet_T_stag"][inE].mean(), s_c=th["jet_s_c"][inE].mean(),
                 Pf=th["jet_P_face"][inE].mean() * scale, Pdec=th["jet_P_decay"][inE].mean() * scale,
                 T40=T40_E, tpE=tpE)
    inW = (th["time"] > W[0]) & (th["time"] <= t_end)
    last = -1
    Pf, cap = th["jet_P_face"] * scale, th["jet_P_cap"] * scale
    fire = ev[(ev["regime"] == 1) & (ev["time"] > W[0])]
    Tf = fire["T_top"] if fire.size else np.array([np.nan])
    tot = th["pinned_cols_pinned"][inW].sum() + th["pinned_cols_face"][inW].sum()
    share = lambda k: th[k][inW].sum() / tot if tot > 0 else float("nan")
    rim = int(np.sum((r < rw[0.002]) & ~fired))
    dd = np.diff(ring_depth[have])
    slope = math.degrees(math.atan(np.nanmax(np.abs(dd)) / D)) if dd.size else 0.0
    hole = depth > 0.002
    oor = hole & ((r < 2.5 * D) | (r > 7.5 * D) | (s_end < 2.0 * D) | (s_end > 12.0 * D))

    # T_gas(r) profiles
    times = np.unique(prof["time"])
    pick = [times[min(len(times) - 1, int(round(q * (len(times) - 1))))] for q in (0.0, 1 / 3, 2 / 3, 1.0)]
    pick = list(dict.fromkeys(pick))
    profiles = {t: (prof["r_lo"][prof["time"] == t] + 0.5 * m["jet_dr"], prof["T_gas"][prof["time"] == t])
                for t in pick}
    tl = times[-1]
    rl, Tl = profiles[tl]
    T40 = float(np.interp(0.040, rl, Tl))
    Trc = float(np.interp(wj.R_CORE, rl, Tl))
    s40 = float(np.mean(s_end[b == i40])) if np.any(b == i40) else wj.SOD
    rop40_cf = float(wj.rop_closed(wj.h_py(m["h_ref"], 0.040, s40), T40))
    okr = (ring_rop >= 0.9 * 1.5) & have
    n_ok = int(np.argmin(okr)) if not okr.all() else nb
    hole_ref = 2.0 * n_ok * D
    ref = wj.ref_table(m["h_ref"], m["T_nozzle"], m.get("stag", "nozzle"), bool(m.get("entrained", 0)))
    P_noz = th["jet_P_cap"] * 0.0 + (m["mdot"] * m["cp"] * (m["T_nozzle"] - m["T_ent"]))
    q40 = float(wj.h_py(m["h_ref"], 0.040, s40) * (T40 - wj.T_FIRE) - wj.loss(wj.T_FIRE))

    return dict(
        early=early, s_c_t=(tg, s_c_t), steady_from=steady_from, sW=sW, hole_ref=hole_ref, ref=ref,
        Tstag=th["jet_T_stag"][last], Tstag_W=th["jet_T_stag"][inW].mean(),
        s_c=th["jet_s_c"][last], Pdec=th["jet_P_decay"][last] * scale,
        Pdec_W=th["jet_P_decay"][inW].mean() * scale, Pnoz=float(P_noz[last]) * scale,
        m=m, t_end=t_end, rw=rw, common=common, rw10_34=rw10_34, r_fire=float(r[fired].max()) if fired.any() else 0.0,
        dom=(m["dim"] == 3), c_depth=depth[ca], c_rop=rop[ca], c_rop1=c_rop1, c_rop2=c_rop2,
        c_s=s_end[ca], Pf_end=Pf[last], Pf_W=Pf[inW].mean(), cap=cap[last],
        Pex=th["jet_P_exhaust"][last] * scale, Tex=th["jet_T_exhaust"][last],
        Tex_W=th["jet_T_exhaust"][inW].mean(), reach=th["jet_r_reach"][last],
        floor=(th["patch_P_robin"][inW] - th["jet_P_face"][inW]).mean() * scale,
        Tf=(np.nanmean(Tf), np.nanpercentile(Tf, 10), np.nanpercentile(Tf, 90)),
        sh=dict(pin=share("pinned_cols_pinned"), unset=share("pinned_cols_face_unset"),
                idle=share("pinned_cols_face_idle"), qneg=share("pinned_cols_face_qneg"),
                minrule=share("pinned_cols_face_minrule")),
        ledger=float(np.max(np.abs(th["ledger_err"]))), rim=rim, slope=slope,
        oor=float(oor.sum() / max(1, hole.sum())), n_hole=int(hole.sum()),
        rb=rb, ring_rop=ring_rop, ring_depth=ring_depth, have=have, rop40=ring_rop[i40],
        rop40_cf=rop40_cf, q40=q40, T40=T40, Trc=Trc, profiles=profiles, s40=s40,
        vol_rate=float(np.sum(rop) / 3600.0 * dz * dz * scale) if m["dim"] == 3 else float("nan"))


def fr(x, d=1):
    return "—" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{d}f}"


def rwall_txt(a, d):
    v = a["rw"][d] * 1e3
    tag = " dom" if a["dom"] and v >= 118.0 else ""
    return f"{v:.0f}{tag}"


def table_main(res, part):
    L = []
    L.append("| run | h_ref | T_noz | treatment | t_end | centre depth | centre ROP (W) | centre s (W mean) | "
             "r_wall 2/10/20/40 mm | hole Ø (10 mm) | hole Ø (ref crit.) | r_fire | ring ROP @40 | "
             "closed form @40 (q) | reference @40 / Ø | steady? |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for n, a in res.items():
        m = a["m"]
        stdy = ("centre ROP {:.1f}→{:.1f} m/h; r_wall(10) {:.0f}→{:.0f} mm; centre s steady from {} s"
                .format(a["c_rop1"], a["c_rop2"], a["rw10_34"] * 1e3, a["rw"][0.010] * 1e3,
                        fr(a["steady_from"], 0)))
        tr = m.get("stag", "?") + (" +mass" if m.get("entrained") else "")
        te = f"{a['t_end']:.0f} s" + (f" of {m['planned_stop']:.0f} s, {m['ended']}" if m.get("ended") else "")
        L.append(f"| {n} | {m['h_ref']:.0f} | {m['T_nozzle']:.0f} | {tr} | {te} | "
                 f"{a['c_depth'] * 1e3:.1f} mm | {a['c_rop']:.2f} m/h | {a['c_s'] * 1e3:.0f} "
                 f"({a['sW'] * 1e3:.0f}) mm | "
                 f"{' / '.join(rwall_txt(a, d) for d in DEPTHS)} | {2 * a['rw'][0.010] * 1e3:.0f} mm | "
                 f"{a['hole_ref'] * 1e3:.0f} mm | "
                 f"{a['r_fire'] * 1e3:.0f} mm | {fr(a['rop40'], 2)} m/h | {fr(a['rop40_cf'], 2)} m/h "
                 f"({a['q40'] / 1e6:.2f} MW/m², s {a['s40'] * 1e3:.0f} mm) | "
                 f"{a['ref']['rop40']:.2f} m/h / {a['ref']['hole'] * 1e3:.0f} mm | {stdy} |")
    return L


def table_power(res):
    L = ["| run | jet_P_face end (W mean) | jet_P_cap | face / cap | face / 2.83 kW | "
         "jet_T_stag end (W mean) | jet_s_c | jet_P_decay (W mean) | nozzle budget | "
         "jet_P_exhaust | jet_T_exhaust (W mean) | jet_r_reach | floor supply | T_gas(r_core) / T_gas(40 mm) |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for n, a in res.items():
        L.append(f"| {n} | {a['Pf_end'] / 1e3:.2f} ({a['Pf_W'] / 1e3:.2f}) kW | {a['cap'] / 1e3:.2f} kW | "
                 f"{a['Pf_end'] / a['cap']:.2f} | {a['Pf_end'] / wj.P_MEIER_ROCK:.1f}× | "
                 f"{a['Tstag']:.0f} ({a['Tstag_W']:.0f}) K | {a['s_c'] * 1e3:.0f} mm | "
                 f"{a['Pdec'] / 1e3:.2f} ({a['Pdec_W'] / 1e3:.2f}) kW | {a['Pnoz'] / 1e3:.2f} kW | "
                 f"{a['Pex'] / 1e3:.2f} kW | {a['Tex']:.0f} ({a['Tex_W']:.0f}) K | "
                 f"{a['reach'] * 1e3:.0f} mm | {a['floor']:.1f} W | {a['Trc']:.0f} / {a['T40']:.0f} K |")
    return L


def table_quality(res):
    L = ["| run | T_fire mean [p10, p90] | pinned | face: unset / idle / qneg / minrule | rim | "
         "ledger max | max slope | out of Martin range | vol. rate |",
         "|---|---|---|---|---|---|---|---|---|"]
    for n, a in res.items():
        s = a["sh"]
        L.append(f"| {n} | {a['Tf'][0]:.1f} [{a['Tf'][1]:.1f}, {a['Tf'][2]:.1f}] K | {s['pin']:.2f} | "
                 f"{s['unset']:.2f} / {s['idle']:.2f} / {s['qneg']:.2f} / {s['minrule']:.3f} | {a['rim']} | "
                 f"{a['ledger']:.1e} | {a['slope']:.0f}° | {a['oor']:.2f} of {a['n_hole']} | "
                 f"{fr(a['vol_rate'] * 1e6, 2)} cm³/s |")
    return L


def table_common(res):
    a0 = next(iter(res.values()))
    if a0["common"] is None:
        return []
    L = [f"At the common time t = {a0['common']['t']:g} s (the shortest run's stop):", "",
         "| run | centre depth | r_fire | r_wall(2 mm) | r_wall(10 mm) | jet_P_face | jet_T_exhaust |",
         "|---|---|---|---|---|---|---|"]
    for n, a in res.items():
        c = a["common"]
        L.append(f"| {n} | {c['c_depth'] * 1e3:.1f} mm | {c['r_fire'] * 1e3:.0f} mm | {c['rw2'] * 1e3:.0f} mm | "
                 f"{c['rw10'] * 1e3:.0f} mm | {c['Pf'] / 1e3:.2f} kW | {c['Tex']:.0f} K |")
    return L


def table_sens(res):
    base = res.get("JM_19_f15")
    if base is None:
        return []
    L = ["| run | centre ROP (W) | ring ROP @40 | r_wall(10 mm) | jet_P_face | jet_T_exhaust | T_gas(40 mm) |",
         "|---|---|---|---|---|---|---|"]
    for n in ("JM_19_f15", "JM_19_f15_cp092", "JM_19_f15_cp108", "JM_19_f15_Tent600", "JM_19_f15_dr1",
              "JM_19_f15_noz", "JM_19_f15_ent"):
        a = res.get(n)
        if a is None:
            continue
        rel = lambda k: "" if n == "JM_19_f15" or not base[k] else f" ({a[k] / base[k] - 1:+.1%})"
        L.append(f"| {n} | {a['c_rop']:.2f}{rel('c_rop')} | {fr(a['rop40'], 2)}{rel('rop40')} | "
                 f"{a['rw'][0.010] * 1e3:.0f} mm | {a['Pf_end']:.0f} W{rel('Pf_end')} | {a['Tex']:.0f} K | "
                 f"{a['T40']:.0f} K |")
    return L


def table_early(res):
    L = ["| run | treatment | window E | t_flank (40 mm ring above nozzle) | ring ROP @40 (E) | reference @40 | "
         "centre ROP (E) | jet_s_c (E mean) | jet_T_stag (E mean) | T_gas(40 mm) | jet_P_face (E mean) | "
         "jet_P_decay (E mean) | hole Ø (10 mm) at E end | reference Ø |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for n, a in res.items():
        e, m = a["early"], a["m"]
        tr = m.get("stag", "?") + (" +mass" if m.get("entrained") else "")
        tf = "not reached" if e["t_flank"] is None else f"{e['t_flank']:.0f} s"
        L.append(f"| {n} | {tr} | {e['E'][0]:.0f}–{e['E'][1]:.0f} s | {tf} | {e['rop40']:.2f} m/h | "
                 f"{a['ref']['rop40']:.2f} m/h | {e['c_rop']:.2f} m/h | {e['s_c'] * 1e3:.0f} mm | "
                 f"{e['Tstag']:.0f} K | {fr(e['T40'], 0)} K ({e['tpE']:.0f} s) | {e['Pf'] / 1e3:.2f} kW | "
                 f"{e['Pdec'] / 1e3:.2f} kW | {2 * e['rw10'] * 1e3:.0f} mm | {a['ref']['hole'] * 1e3:.0f} mm |")
    return L


def table_treat(res):
    L = ["| run | ring ROP @40 (W) | hole Ø (ref crit.) | hole Ø (10 mm) | r_fire | centre ROP (W) | "
         "centre s (W mean) | jet_T_stag (W mean) | jet_P_face (W mean) |",
         "|---|---|---|---|---|---|---|---|---|"]
    for base, alts in (("J5_19", ("J5_19_ent", "J5_19_noz")), ("J10_19", ("J10_19_ent",))):
        for n in (base,) + alts:
            a = res.get(n)
            if a is None:
                continue
            L.append(f"| {n} | {fr(a['rop40'], 2)} m/h | {a['hole_ref'] * 1e3:.0f} mm | "
                     f"{2 * a['rw'][0.010] * 1e3:.0f} mm | {a['r_fire'] * 1e3:.0f} mm | {a['c_rop']:.2f} m/h | "
                     f"{a['sW'] * 1e3:.0f} mm | {a['Tstag_W']:.0f} K | {a['Pf_W'] / 1e3:.2f} kW |")
    return L


def plots(res, tag):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print("no matplotlib:", e)
        return []
    names = [n for n in res if n in ANCHORS or n.endswith(("_f0", "_f15"))]
    names = [n for n in names if n.split("_f")[0] in ANCHORS] if tag == "part1" else \
        [n for n in res if n in ANCHORS or n.endswith(("_ent", "_noz"))]
    fig, axs = plt.subplots(1, 4, figsize=(21, 4.8))
    for n in names:
        a = res[n]
        t = sorted(a["profiles"])[-1]
        rr, TT = a["profiles"][t]
        axs[0].plot(rr * 1e3, TT, label=f"{n} ({t:.0f} s)")
        axs[1].plot(a["rb"][a["have"]] * 1e3, -a["ring_depth"][a["have"]] * 1e3, label=n)
        axs[2].plot(a["rb"][a["have"]] * 1e3, a["ring_rop"][a["have"]], label=n)
        axs[3].plot(a["s_c_t"][0], a["s_c_t"][1] * 1e3, label=n)
    axs[0].axhline(wj.T_FIRE, color="k", lw=0.8, ls=":")
    axs[0].set(xlabel="r [mm]", ylabel="T_gas (bin inlet) [K]", title="T_gas(r) at t_end")
    axs[1].set(xlabel="r [mm]", ylabel="ring-mean depth [mm]", title="face at t_end")
    axs[2].axhline(1.5, color="k", lw=0.8, ls=":")
    axs[2].set(xlabel="r [mm]", ylabel="ring ROP in W [m/h]", title="ring ROP", yscale="log")
    axs[3].axhline(wj.CORE_LEN * D * 1e3, color="k", lw=0.8, ls=":")
    axs[3].set(xlabel="t [s]", ylabel="centre stand-off [mm]", title="centre s(t)")
    for ax in axs[:3]:
        ax.axvline(40.0, color="grey", lw=0.6, ls="--")
        ax.grid(alpha=0.3)
    axs[2].legend(fontsize=7)
    fig.tight_layout()
    f1 = os.path.join(HERE, f"{tag}_profiles.png")
    fig.savefig(f1, dpi=110)
    plt.close(fig)
    base = [n for n in names if n in ANCHORS or n.split("_f")[0] in ANCHORS][:6]
    if tag == "part1":
        base = [n for n in names if n.endswith("_f15")][:6]
    fig, axs = plt.subplots(2, 3, figsize=(15, 8), sharex=True, sharey=True)
    for ax, n in zip(axs.ravel(), base):
        a = res[n]
        for t in sorted(a["profiles"]):
            rr, TT = a["profiles"][t]
            ax.plot(rr * 1e3, TT, label=f"{t:.0f} s")
        ax.axhline(wj.T_FIRE, color="k", lw=0.8, ls=":")
        ax.set_title(n)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
    for ax in axs[-1]:
        ax.set_xlabel("r [mm]")
    for ax in axs[:, 0]:
        ax.set_ylabel("T_gas [K]")
    fig.tight_layout()
    f2 = os.path.join(HERE, f"{tag}_Tgas_times.png")
    fig.savefig(f2, dpi=110)
    plt.close(fig)
    return [f1, f2]


def main():
    out = []
    for part, names, tag in ((1, PART1, "part1"), (2, PART2, "part2")):
        res = {}
        metas = {n: meta_of(os.path.join(OUT, f"p{part}_{n}")) for n in names}
        metas = {n: mm for n, mm in metas.items() if mm is not None}
        if not metas:
            continue
        t_c = min(float(mm["stop"]) for mm in metas.values())
        for n, mm in metas.items():
            res[n] = analyze(os.path.join(OUT, f"p{part}_{n}"), mm, t_c)
        if not res:
            continue
        figs = plots(res, tag)
        out.append(f"## {'1. Part 1: 2-D slab (machinery and trend only)' if part == 1 else '2. Part 2: 3-D quarter domain, feed 1.5 m/h'}")
        out.append("")
        if part == 1:
            out.append(f"Slab jet share mdot·w/(π r_core) = {res[next(iter(res))]['m']['mdot']:.4g} kg/s; "
                       "powers are slab values. No number from this part enters the conclusions.")
        else:
            out.append("Quarter domain with mdot/4; powers are ×4 (full jet).")
        out.append("")
        out += table_main(res, part) + [""] + table_common(res) + [""] + table_power(res) + [""] \
            + table_quality(res) + [""]
        if part == 1:
            s = table_sens(res)
            if s:
                out += ["Sensitivities (JM_19, feed 1.5 m/h):", ""] + s + [""]
        else:
            out += ["**Early window (the closure's numbers; read these, not W, for the decay runs — "
                    "see §3):**", ""] + table_early(res) + [""]
            out += ["Stagnation treatments (1900 K), window W:", ""] + table_treat(res) + [""]
        out += [f"![{os.path.basename(f)}]({os.path.basename(f)})" for f in figs] + [""]
        # profiles table
        out.append("T_gas(r) [K] at the profile times (bin inlet, interpolated):")
        out.append("")
        out.append("| run | t | r = 10 | 18.75 | 30 | 40 | 60 | 80 | 100 mm |")
        out.append("|---|---|---|---|---|---|---|---|---|")
        for n, a in res.items():
            if not (n in ANCHORS or n.endswith("_f15")):
                continue
            for t in sorted(a["profiles"]):
                rr, TT = a["profiles"][t]
                vals = [f"{np.interp(x, rr, TT, right=np.nan):.0f}" for x in (0.010, wj.R_CORE, 0.030, 0.040,
                                                                               0.060, 0.080, 0.100)]
                out.append(f"| {n} | {t:.0f} s | " + " | ".join(vals) + " |")
        out.append("")
        meta_lines = [f"{n}: dt {a['m']['dt'] * 1e3:g} ms, overshoot {a['m']['overshoot']:.2f} K/step, "
                      + (f"wall {a['m']['wall'] / 60:.1f} min" if 'wall' in a['m'] else "ended early")
                      for n, a in res.items()]
        out.append("Run settings: " + "; ".join(meta_lines) + ".")
        out.append("")
    path = os.path.join(HERE, "RESULTS.md")
    txt = open(path).read()
    head = txt.split(MARK_R)[0]
    tail = txt.split(MARK_D, 1)[1] if MARK_D in txt else "\n"
    open(path, "w").write(head + MARK_R + "\n\n" + "\n".join(out) + "\n" + MARK_D + tail)
    print("\n".join(out))


if __name__ == "__main__":
    main()
