#!/usr/bin/env python3
"""D2b analysis -> RESULTS.md (between <!-- RESULTS --> and <!-- DISCUSSION -->;
the pre-registration above and the discussion below are preserved) + PNGs.

  analyze.py            all stages found in output/
  analyze.py --stage A  print the Stage A decision only

Definitions (quarter domain; powers and volumes x4 = full hole):
  depth(col, t)   original top (lz) minus the mask top face after the column's
                  last event at or before t (k_top - n_voided + 1)*dz.
  ROP             -d(nozzle_z)/dt, least squares over the steady window;
                  cross-check: mean over annulus columns of the least-squares
                  slope of cumulative h_applied in the window.
  steady window   the earliest t0 such that [t0, t_end] (>= W_min: 150 s for B/C,
                  60 s for D) splits into 50 s sub-windows whose descent-rate fits
                  are all within +-10 % of their mean, and |d jet_s_c/dt| < 0.02
                  mm/s over it; t0 on a 10 s grid. If none: "not steady", and the
                  last max(W_min, t_end/2) s are used, flagged.
  r_wall(z)       max r of a column with depth > z at t_end, z = depth below the
                  original top; drilled depth = the feet depth lz - foot_z(t_end).
                  depth-mean Ø = mean over z in [0, drilled depth] of 2 r_wall.
  block profile   final depth per column: frozen (face above the nozzle) as
                  drilled; active projected with its radius bin's window rate
                  v(r) until the nozzle passes it. Hole Ø is scored over
                  0-0.5 m of this profile (depth-mean and minimum).
  in hole         r <= R_HOLE = 46.5 mm (Ø 93); the scored volume rate and
                  the rock-side face power are taken there.
  volume rate     (a) d/dt of sum_z pi r_wall(z)^2 dz from wall profiles every
                  10 s in the window (the packet's); (b) removed volume rate
                  4 * sum h_applied dx dy / window (every removal, pit included).
  R               rho Cp (T_fire - 293.15) * (b) / (4 * mean jet_P_face); per 1 D
                  ring with the ring's removal and the ring's P_bin (jet profile).
  T_fire          regime-1 T_top in the window: mean, p10, p90, p99.
  needles         columns whose end depth exceeds every in-domain 4-neighbour's
                  by > 2 dz.
  rim             never-fired columns with r < r_wall(2 mm).
  out of range    in-hole columns (depth > 2 mm) with r outside [2.5, 7.5] D or s
                  outside [2, 12] D at t_end.
"""
import argparse
import ast
import glob
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import feetmodel as fm  # noqa: E402

wj = fm.wj
D = wj.D
OUT = os.path.join(HERE, "output")
MARK_R, MARK_D = "<!-- RESULTS -->", "<!-- DISCUSSION -->"
TFIRE_BAND_HI, OVERHEAT = 823.0, 50.0
MDOT = wj.MDOT
P_MEIER = wj.P_MEIER_ROCK
BAND_ROP = (1.04, 1.92)
R_HOLE = 0.0465   # Ø 93 mm: Meier's 3.42 L over the 0.5 m block
BAND_D = (85.0, 93.0)
BAND_V = (1.98, 2.96)
BAND_DT = (500.0, 560.0)


def meta_of(pf):
    try:
        return ast.literal_eval(open(pf + ".done").read())
    except Exception:
        return None


def thermo(pf):
    lines = [l.split() for l in open(os.path.join(pf, "thermo.dat")) if l.strip()]
    a = np.array(lines[1:], dtype=float)
    return {n: a[:, i] for i, n in enumerate(lines[0])}


def csv(path):
    return np.genfromtxt(path, delimiter=",", names=True, ndmin=1)


class Run:
    def __init__(self, pf, m):
        self.pf, self.m = pf, m
        self.name = m["name"]
        self.th = thermo(pf)
        self.ev = csv(pf + "_removal_events.csv")
        self.prof = csv(pf + "_jet_profile.csv") if os.path.exists(pf + "_jet_profile.csv") else None
        self.dz = m["dz"]
        self.nx = m["nxy"]
        self.lz = m["lz"]
        self.t_end = float(self.th["time"][-1])
        xc = (np.arange(self.nx) + 0.5) * self.dz
        X, Y = np.meshgrid(xc, xc, indexing="xy")
        self.r = np.hypot(X, Y).ravel()
        self.ncol = self.nx * self.nx
        col = self.ev["col_j"].astype(int) * self.nx + self.ev["col_i"].astype(int)
        order = np.lexsort((self.ev["time"], col))
        self.ev, self.col = self.ev[order], col[order]
        self.st = np.searchsorted(self.col, np.arange(self.ncol + 1))
        self.ann = np.nonzero((self.r >= 0.028) & (self.r <= 0.040))[0]
        self.fired = np.array([np.any(self.ev["regime"][self.st[c]:self.st[c + 1]] == 1)
                               for c in range(self.ncol)])

    def depth(self, t):
        """mask depth of every column at time t"""
        d = np.zeros(self.ncol)
        for c in range(self.ncol):
            e = self.ev[self.st[c]:self.st[c + 1]]
            if e.size == 0:
                continue
            i = np.searchsorted(e["time"], t, side="right") - 1
            if i >= 0:
                d[c] = self.lz - (e["k_top"][i] - e["n_voided"][i] + 1) * self.dz
        return d

    def happ(self, c, w):
        e = self.ev[self.st[c]:self.st[c + 1]]
        sel = (e["time"] > w[0]) & (e["time"] <= w[1])
        return e[sel]

    def col_fit(self, c, w):
        e = self.ev[self.st[c]:self.st[c + 1]]
        cum = np.cumsum(e["h_applied"])
        sel = (e["time"] > w[0]) & (e["time"] <= w[1])
        if sel.sum() >= 3 and np.ptp(e["time"][sel]) >= 0.5 * (w[1] - w[0]):
            return float(np.polyfit(e["time"][sel], cum[sel], 1)[0])
        i0 = np.searchsorted(e["time"], w[0], side="right")
        i1 = np.searchsorted(e["time"], w[1], side="right")
        return float(np.sum(e["h_applied"][i0:i1])) / (w[1] - w[0])

    def rate(self, w):
        th = self.th
        sel = (th["time"] >= w[0]) & (th["time"] <= w[1])
        return -float(np.polyfit(th["time"][sel], th["nozzle_z"][sel], 1)[0])

    def steady(self, w_min):
        th = self.th
        best = None
        for t0 in np.arange(0.0, self.t_end - w_min + 1e-9, 10.0):
            n = int((self.t_end - t0) // 50.0)
            if n < 2:
                continue
            subs = [self.rate((self.t_end - 50.0 * (k + 1), self.t_end - 50.0 * k)) for k in range(n)]
            mu = np.mean(subs)
            if mu <= 0.0 or np.max(np.abs(np.array(subs) / mu - 1.0)) > 0.10:
                continue
            sel = th["time"] >= t0
            dsc = abs(float(np.polyfit(th["time"][sel], th["jet_s_c"][sel], 1)[0]))
            if dsc >= 2e-5:
                continue
            best = (t0, subs, dsc)
            break
        return best

    def r_wall(self, d, z):
        sel = d > z
        return float(self.r[sel].max()) if sel.any() else 0.0

    def ring_of(self, c):
        return int(self.r[c] / D)


def analyze(run, kind):
    m, th = run.m, run.th
    w_min = 60.0 if kind == "D" else 150.0
    if kind == "A":
        W = (min(50.0, run.t_end / 3.0), min(150.0, run.t_end))
        steady = None
    else:
        s = run.steady(w_min)
        if s:
            W = (s[0], run.t_end)
            steady = s
        else:
            W = (max(0.0, run.t_end - max(w_min, run.t_end / 2.0)), run.t_end)
            steady = None
    inW = (th["time"] > W[0]) & (th["time"] <= W[1])
    dur = W[1] - W[0]
    res = dict(run=run, W=W, steady=steady, status=m.get("status", "?"), t_end=run.t_end)
    # ROP
    rop = run.rate(W) * 3600.0 if dur > 0 else float("nan")
    ring_rop = np.mean([run.col_fit(c, W) for c in run.ann]) * 3600.0
    res.update(rop=rop, ring_rop=ring_rop)
    selW = (th["time"] >= W[0]) & (th["time"] <= W[1])
    res["sc_drift"] = float(np.polyfit(th["time"][selW], th["jet_s_c"][selW], 1)[0]) if dur > 0 else float("nan")
    res["Pf_drift"] = 4.0 * float(np.polyfit(th["time"][selW], th["jet_P_face"][selW], 1)[0]) if dur > 0 else float("nan")
    if steady:
        res["subs"] = [x * 3600.0 for x in steady[1]]
    else:
        n = max(1, int(dur // 50.0))
        res["subs"] = [run.rate((W[1] - 50.0 * (k + 1), W[1] - 50.0 * k)) * 3600.0
                       for k in range(n) if W[1] - 50.0 * (k + 1) >= W[0] - 1e-9]
    # temperatures
    fire = run.ev[(run.ev["regime"] == 1) & (run.ev["time"] > W[0]) & (run.ev["time"] <= W[1])]
    Tf = fire["T_top"] if fire.size else np.array([np.nan])
    res["Tf"] = (np.nanmean(Tf), np.nanpercentile(Tf, 10), np.nanpercentile(Tf, 90), np.nanpercentile(Tf, 99),
                 np.nanmax(Tf))
    fcol = fire["col_j"].astype(int) * run.nx + fire["col_i"].astype(int) if fire.size else np.zeros(0, int)
    ring_fire = np.isin(fcol, run.ann)
    res["Tf_ring_max"] = float(np.max(fire["T_top"][ring_fire])) if ring_fire.any() else float("nan")
    dT = res["Tf"][0] - 293.15
    # removal and power
    hsum = np.sum(run.ev["h_applied"][(run.ev["time"] > W[0]) & (run.ev["time"] <= W[1])])
    vdot = 4.0 * hsum * run.dz * run.dz / dur
    Pf = 4.0 * float(np.mean(th["jet_P_face"][inW]))
    res.update(vdot=vdot, P_face=Pf, R=wj.RHOCP * dT * vdot / Pf if Pf > 0 else float("nan"), dT=dT)
    # per ring R
    nb = int(math.ceil(run.r.max() / D))
    rings = np.minimum((run.r / D).astype(int), nb - 1)
    vr = np.zeros(nb)
    ev_w = run.ev[(run.ev["time"] > W[0]) & (run.ev["time"] <= W[1])]
    cw = ev_w["col_j"].astype(int) * run.nx + ev_w["col_i"].astype(int)
    np.add.at(vr, rings[cw], ev_w["h_applied"])
    vr *= 4.0 * run.dz * run.dz / dur
    pr = np.zeros(nb)
    if run.prof is not None:
        ps = run.prof[(run.prof["time"] > W[0]) & (run.prof["time"] <= W[1])]
        nt = max(1, np.unique(ps["time"]).size)
        rb = np.minimum(((ps["r_lo"] + 0.5 * run.dz) / D).astype(int), nb - 1)
        np.add.at(pr, rb, ps["P_bin"])
        pr *= 4.0 / nt
    with np.errstate(divide="ignore", invalid="ignore"):
        Rr = np.where(pr > 0, wj.RHOCP * dT * vr / pr, np.nan)
    ring_bins = sorted(set(rings[run.ann]))
    res.update(R_ring=Rr, P_ring=pr, V_ring=vr, ring_bins=ring_bins,
               R_foot=float(np.nansum(vr[ring_bins]) * wj.RHOCP * dT / np.nansum(pr[ring_bins]))
               if np.nansum(pr[ring_bins]) > 0 else float("nan"))
    # geometry at t_end
    d_end = run.depth(run.t_end)
    z_feet = run.lz - float(th["foot_z"][-1])
    zg = np.arange(0.0, max(z_feet, run.dz), run.dz / 2.0)
    rw = np.array([run.r_wall(d_end, z) for z in zg])
    res.update(d_end=d_end, z_feet=z_feet, zg=zg, rw=rw,
               dmean=2.0 * float(np.mean(rw)) if rw.size else float("nan"),
               dmin=2.0 * float(np.min(rw)) if rw.size else float("nan"),
               c_depth=float(d_end[int(np.argmin(run.r))]),
               dom_limited=bool(rw.size and rw.max() >= 0.99 * run.m["lxy"]))
    # D2b scoring inside the hole (set before the Stage B results): Meier's
    # volume rate and rock-side power are ROP x nominal hole area, so score
    # removal and face power inside r <= R_HOLE (Ø 93 = 3.42 L / 0.5 m)
    ev_h = ev_w[run.r[cw] <= R_HOLE]
    res["vdot_hole"] = 4.0 * float(np.sum(ev_h["h_applied"])) * run.dz * run.dz / dur
    if run.prof is not None:
        ph = ps[(ps["r_lo"] + 0.5 * run.dz) <= R_HOLE]
        res["P_hole"] = 4.0 * float(np.sum(ph["P_bin"])) / nt
    else:
        res["P_hole"] = float("nan")
    # Ø over the 0.5 m block: every column's final depth. A column whose face
    # is above the nozzle (depth <= nozzle depth n) is frozen and keeps its
    # depth. An active column recedes at v(r) until the nozzle, descending at
    # v, passes it: t* = (d - n) / (v - v(r)), final depth d + v(r) t*
    # (never frozen if v(r) >= v). v(r) = window recession rate of the
    # column's radius bin (width dz); v = burner ROP. r_wall(z) = max r of a
    # column whose final depth exceeds z.
    v = rop / 3600.0
    rb = (run.r / run.dz).astype(int)
    vcol = np.array([run.col_fit(c, W) for c in range(run.ncol)])
    nbin = rb.max() + 1
    cnt = np.bincount(rb, minlength=nbin)
    vbin = (np.bincount(rb, weights=vcol, minlength=nbin) / np.maximum(cnt, 1))[rb]
    n_depth = run.lz - float(th["nozzle_z"][-1])
    active = d_end > n_depth
    with np.errstate(divide="ignore", invalid="ignore"):
        tstar = np.where(vbin < v, (d_end - n_depth) / (v - vbin), np.inf)
        final = np.where(active, np.where(np.isfinite(tstar), d_end + vbin * tstar, np.inf), d_end)
    zb = np.arange(0.0, fm.BLOCK, run.dz / 2.0)
    rw_b = np.array([run.r_wall(final, z) for z in zb])
    res.update(d_final=final, zb=zb, rw_b=rw_b,
               dmean_block=2.0 * float(np.mean(rw_b)), dmin_block=2.0 * float(np.min(rw_b)),
               vdot_block=v * float(np.mean(math.pi * rw_b ** 2)),
               z93=float(zb[np.argmax(rw_b <= R_HOLE)]) if np.any(rw_b <= R_HOLE) else float("nan"),
               z85=float(zb[np.argmax(rw_b <= 0.0425)]) if np.any(rw_b <= 0.0425) else float("nan"))
    # volume rate from wall profiles
    ts = np.arange(W[0], W[1] + 1e-9, max(10.0, dur / 10.0))
    vols = []
    for t in ts:
        dt_ = run.depth(t)
        zmax = dt_.max()
        zz = np.arange(0.0, zmax, run.dz)
        vols.append(float(np.sum([math.pi * run.r_wall(dt_, z) ** 2 * run.dz for z in zz])))
    res["vdot_wall"] = float(np.polyfit(ts, vols, 1)[0]) if len(ts) >= 2 else float("nan")
    # checks
    d2 = run.r_wall(d_end, 0.002)
    rim = [c for c in range(run.ncol) if run.r[c] < d2 and not run.fired[c]]
    res["rim"] = len(rim)
    res["rim_ring"] = int(np.sum(np.isin(rim, run.ann)))
    nx = run.nx
    needles = []
    D2 = d_end.reshape(nx, nx)
    for j in range(nx):
        for i in range(nx):
            nb4 = [D2[jj, ii] for jj, ii in ((j, i - 1), (j, i + 1), (j - 1, i), (j + 1, i))
                   if 0 <= ii < nx and 0 <= jj < nx]
            if nb4 and D2[j, i] - max(nb4) > 2.0 * run.dz + 1e-12:
                needles.append(j * nx + i)
    res["needles"] = len(needles)
    res["needles_ring"] = int(np.sum(np.isin(needles, run.ann)))
    rmean = np.array([np.mean(d_end[rings == b]) if np.any(rings == b) else np.nan for b in range(nb)])
    dd = np.diff(rmean[np.isfinite(rmean)])
    res["slope"] = math.degrees(math.atan(np.nanmax(np.abs(dd)) / D)) if dd.size else 0.0
    res["ring_mean_depth"] = rmean
    tot = th["pinned_cols_pinned"][inW].sum() + th["pinned_cols_face"][inW].sum()
    res["sh"] = {k: th[k][inW].sum() / tot for k in ("pinned_cols_pinned", "pinned_cols_face_unset",
                                                     "pinned_cols_face_idle", "pinned_cols_face_qneg",
                                                     "pinned_cols_face_minrule")} if tot > 0 else {}
    res["ledger"] = float(np.max(np.abs(th["ledger_err"])))
    s_end = th["nozzle_z"][-1] - (run.lz - d_end)
    hole = d_end > 0.002
    oor = hole & ((run.r < 2.5 * D) | (run.r > 7.5 * D) | (s_end < 2 * D) | (s_end > 12 * D))
    res["oor"] = float(oor.sum() / max(1, hole.sum()))
    for k in ("jet_s_c", "jet_T_stag", "jet_T_rec", "jet_T_exhaust", "jet_P_cap"):
        if k in th:
            res[k] = float(np.mean(th[k][inW]))
    res["gap"] = float(np.mean(th["nozzle_z"][inW] - th["foot_z"][inW]))
    res["P_cap"] = 4.0 * res.get("jet_P_cap", float("nan"))
    res["P_net"] = MDOT * wj.CP * (m["T_nozzle"] - res.get("jet_T_exhaust", float("nan")))
    res["carry"] = (float(np.min(th["foot_carry_cols"][inW])), float(np.mean(th["foot_carry_cols"][inW])))
    # longest hold: simulated time without a nozzle_z decrease (from nozzle_z, so
    # it works for runs made before foot_stall_time replaced foot_stalled_steps)
    live = th["time"] > 0
    tz, zz = th["time"][live], th["nozzle_z"][live]
    t_last, hold = tz[0], 0.0
    for k in range(1, tz.size):
        if zz[k] < zz[k - 1]:
            t_last = tz[k]
        hold = max(hold, tz[k] - t_last)
    res["hold"] = float(hold)
    res["hold_col"] = float(np.max(th["foot_stall_time"])) if "foot_stall_time" in th else float("nan")
    # thermo.dat is written every few steps; count solver steps from dt
    res["steps"] = int(round(th["time"][-1] / m["dt"])) if "dt" in m else int(np.sum(th["time"] > 0))
    return res


def stage_a_decision(res):
    lines = ["| treatment | R (hole) | R (foot rings) | T_fire p99 / max | ring T_top max | needles (ring) | "
             "never-fired in hole (ring) | passes 1-3 | ring clean |",
             "|---|---|---|---|---|---|---|---|---|"]
    passing, ring_clean = [], []
    for w in ("A1", "A2", "A3"):
        r = res.get(f"A_{w}")
        if r is None:
            continue
        p1 = r["R"] <= 1.05
        p2 = r["Tf"][3] <= TFIRE_BAND_HI + OVERHEAT
        p3 = r["needles"] == 0 and r["rim"] == 0
        rc = (r["R_foot"] <= 1.05 and not (r["Tf_ring_max"] > TFIRE_BAND_HI + OVERHEAT)
              and r["rim_ring"] == 0 and r["needles_ring"] == 0)
        if p1 and p2 and p3:
            passing.append(w)
        if rc:
            ring_clean.append(w)
        lines.append(f"| {w} | {r['R']:.3f} | {r['R_foot']:.3f} | {r['Tf'][3]:.0f} / {r['Tf'][4]:.0f} K | "
                     f"{r['Tf_ring_max']:.0f} K | {r['needles']} ({r['needles_ring']}) | {r['rim']} ({r['rim_ring']}) | "
                     f"{'yes' if p1 and p2 and p3 else 'no'} ({'R' if not p1 else ''}{'T' if not p2 else ''}"
                     f"{'A' if not p3 else ''}) | {'yes' if rc else 'no'} |")
    order = ["A2", "A3", "A1"]
    if passing:
        choice = next(w for w in order if w in passing)
        why = f"passes rules 1-3 (passers: {', '.join(passing)}); order A2 > A3 > A1"
    elif ring_clean:
        choice = next(w for w in order if w in ring_clean)
        why = (f"no treatment passes 1-3; tie-break (rule 5): ring-clean treatments {', '.join(ring_clean)}; "
               f"order A2 > A3 > A1; its violations are carried as a caveat")
    else:
        choice, why = None, "every treatment violates inside the ring: STOP (rule 6)"
    return lines, choice, why


def hand_rwall(hk, Tn, mode, z):
    roots, _, _ = fm.equilibria(float(wj.h_anchor(hk, Tn)), Tn, mode, True)
    if not roots:
        return None
    prof = fm.depth_profile(float(wj.h_anchor(hk, Tn)), roots[0])
    return np.interp(z, prof["z"], prof["r_wall"]), roots[0], prof


def fmt(x, d=2):
    return "—" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{d}f}"


def band(x, b):
    return "PASS" if b[0] <= x <= b[1] else "FAIL"


def run_table(res):
    L = ["| run | status / t_end | window | ROP (burner) | ring columns | sub-fits (50 s) | T_fire mean [p10, p90] | "
         "ΔT_fire | jet_s_c | T_stag | T_rec | gap | face power ×4 | cap ×4 | mdot cp (T_noz − T_exh) | "
         "face power inside Ø 93 ×4 (vs 2.83 kW) | R | removed vol. rate | wall vol. rate | depth-mean Ø (drilled) | min Ø | feet depth | "
         "centre depth |",
         "|" + "---|" * 23]
    for n, r in res.items():
        W = r["W"]
        L.append(
            f"| {n} | {r['status']} / {r['t_end']:.0f} s | {W[0]:.0f}–{W[1]:.0f} s{' (steady)' if r['steady'] else ' (not steady)'} | "
            f"{fmt(r['rop'])} m/h | {fmt(r['ring_rop'])} m/h | {', '.join(f'{x:.2f}' for x in r['subs'])} | "
            f"{r['Tf'][0]:.1f} [{r['Tf'][1]:.1f}, {r['Tf'][2]:.1f}] K | {r['dT']:.0f} K | "
            f"{fmt(r.get('jet_s_c', float('nan')) * 1e3, 0)} mm | {fmt(r.get('jet_T_stag'), 0)} K | "
            f"{fmt(r.get('jet_T_rec'), 0)} K | {r['gap'] * 1e3:.1f} mm | {r['P_face'] / 1e3:.2f} kW | "
            f"{r['P_cap'] / 1e3:.2f} kW | {r['P_net'] / 1e3:.2f} kW | {r['P_hole'] / 1e3:.2f} kW ({r['P_hole'] / P_MEIER:.2f}×) | {r['R']:.3f} | "
            f"{r['vdot'] * 1e6:.2f} cm³/s | {r['vdot_wall'] * 1e6:.2f} cm³/s | {r['dmean'] * 1e3:.1f} mm"
            f"{' (dom)' if r['dom_limited'] else ''} | {r['dmin'] * 1e3:.1f} mm | {r['z_feet'] * 1e3:.0f} mm | "
            f"{r['c_depth'] * 1e3:.0f} mm |")
    return L


def check_table(res):
    L = ["| run | pinned | face: unset / idle / qneg / minrule | foot_carry_cols min / mean | longest hold [s] | "
         "rim (ring) | needles (ring) | R foot rings | T_fire p99 | ledger max | max slope | out of Martin range |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for n, r in res.items():
        s = r["sh"]
        L.append(f"| {n} | {s.get('pinned_cols_pinned', float('nan')):.2f} | "
                 f"{s.get('pinned_cols_face_unset', 0):.2f} / {s.get('pinned_cols_face_idle', 0):.2f} / "
                 f"{s.get('pinned_cols_face_qneg', 0):.2f} / {s.get('pinned_cols_face_minrule', 0):.3f} | "
                 f"{r['carry'][0]:.0f} / {r['carry'][1]:.1f} | {r['hold']:.1f} | "
                 f"{r['rim']} ({r['rim_ring']}) | {r['needles']} ({r['needles_ring']}) | {r['R_foot']:.3f} | "
                 f"{r['Tf'][3]:.0f} K | {r['ledger']:.1e} | {r['slope']:.0f}° | {r['oor']:.2f} |")
    return L


def scored_table(r, hand, mesh=None):
    rop, dm, dmin, v = r["rop"], r["dmean"] * 1e3, r["dmin"] * 1e3, r["vdot_wall"] * 1e6
    hz = hand if hand is not None else dict(v=float("nan"), dm_drill=float("nan"), dm_block=float("nan"),
                                            dmin=float("nan"), vdot=float("nan"))
    L = ["| criterion | target | simulation | verdict | hand prediction |", "|---|---|---|---|---|"]
    L.append(f"| ROP | [1.04, 1.92] m/h | {fmt(rop)} m/h (ring columns {fmt(r['ring_rop'])}) | "
             f"{band(rop, BAND_ROP)} | {fmt(hz['v'])} m/h{'' if hand is not None else ' (stall predicted)'} |")
    db, dbmin, vh = r["dmean_block"] * 1e3, r["dmin_block"] * 1e3, r["vdot_hole"] * 1e6
    L.append(f"| hole Ø, depth-mean over 0–0.5 m | 85–93 mm | {db:.1f} mm (frozen columns as drilled, active "
             f"columns projected with the run's own v(r); Ø ≤ 93 from {fmt(r['z93'] * 1e3, 0)} mm, "
             f"≤ 85 from {fmt(r['z85'] * 1e3, 0)} mm); drilled depth only: {dm:.1f} mm | "
             f"{band(db, BAND_D)} | {fmt(hz['dm_block'], 1)} mm; {fmt(hz['dm_drill'], 1)} mm over the drilled depth |")
    L.append(f"| hole Ø, minimum over 0–0.5 m | ≥ 80 mm | {dbmin:.1f} mm (drilled depth: {dmin:.1f} mm) | "
             f"{'PASS' if dbmin >= 80.0 else 'FAIL'} | {fmt(hz['dmin'], 1)} mm |")
    L.append(f"| volume rate inside Ø 93 | [1.98, 2.96] cm³/s | {vh:.2f} cm³/s; whole top {r['vdot'] * 1e6:.2f}, "
             f"wall profile {v:.2f}, ROP × block-mean area {r['vdot_block'] * 1e6:.2f} | "
             f"{band(vh, BAND_V)} | {fmt(hz['vdot'])} cm³/s (π r_wall² at the ring rate) |")
    L.append(f"| ΔT_fire | 500–560 K | {r['dT']:.0f} K | {band(r['dT'], BAND_DT)} | 528 K |")
    flat = max(abs(x / np.mean(r['subs']) - 1.0) for x in r['subs']) if r['subs'] else float("nan")
    L.append(f"| ROP flatness | sub-fits within ±10 % | max {flat * 100:.1f} %; jet_s_c drift {r['sc_drift'] * 1e3:.2f} mm/s (limit 0.02); "
             f"face power ×4 drift {r['Pf_drift']:+.0f} W/s ({'steady' if r['steady'] else 'no steady window'}) | "
             f"{'PASS' if r['steady'] else 'FAIL'} | — |")
    L.append(f"| ledger | round-off | {r['ledger']:.1e} | {'PASS' if r['ledger'] < 1e-10 else 'FAIL'} | — |")
    if mesh:
        L.append(f"| mesh 2 → 1 mm | ring ROP and T_fire within 5 % | {mesh['text']} | {mesh['verdict']} | — |")
    else:
        L.append("| mesh 2 → 1 mm | within 5 % | not evaluated | — | — |")
    return L


def hand_numbers(hk, Tn, mode, entrained, z_drill):
    h = float(wj.h_anchor(hk, Tn))
    roots, _, _ = fm.equilibria(h, Tn, mode, bool(entrained))
    if not roots:
        return None
    eq = roots[0]
    prof = fm.depth_profile(h, eq)
    zsel = prof["z"] <= z_drill
    dm_drill = float(np.mean(2.0 * prof["r_wall"][zsel])) * 1e3 if zsel.any() else float("nan")
    # hand volume rate: the hole deepens at v with the deep wall radius; at depth
    # z the wall is r_wall(z), fixed once frozen, so the removed volume per unit
    # burner descent tends to pi r_q^2 plus the frozen collar (included once).
    vdot = math.pi * float(prof["r_wall"][-1]) ** 2 * eq["v"] / 3600.0 * 1e6
    return dict(v=eq["v"], dm_drill=dm_drill, dm_block=fm.depth_mean(prof, 0.0, fm.BLOCK) * 1e3,
                dmin=2.0 * float(prof["r_wall"].min()) * 1e3, vdot=vdot, prof=prof, eq=eq)


def plots(res, tag):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return []
    fig, axs = plt.subplots(1, 4, figsize=(21, 4.8))
    for n, r in res.items():
        th = r["run"].th
        live = th["time"] > 0  # the t = 0 row is written before the feet are placed
        axs[0].plot(th["time"][live], (r["run"].lz + 0.05 - th["nozzle_z"][live]) * 1e3, label=n)
        axs[1].plot(th["time"][live], th["jet_s_c"][live] * 1e3, label=n)
        axs[2].plot(r["rw"] * 1e3, -r["zg"] * 1e3, label=n)
        if r["run"].prof is not None:
            p = r["run"].prof
            tl = np.unique(p["time"])[-1]
            s = p["time"] == tl
            axs[3].plot((p["r_lo"][s] + 0.5 * r["run"].dz) * 1e3, p["T_gas"][s], label=f"{n} ({tl:.0f} s)")
    axs[0].set(xlabel="t [s]", ylabel="burner descent [mm]", title="nozzle descent")
    axs[1].set(xlabel="t [s]", ylabel="jet_s_c [mm]", title="centre stand-off")
    axs[2].axvline(40.0, color="k", lw=0.6, ls=":")
    axs[2].set(xlabel="r_wall [mm]", ylabel="depth [mm]", title="wall profile at t_end")
    axs[3].set(xlabel="r [mm]", ylabel="T_gas [K]", title="T_gas(r), last profile")
    for ax in axs:
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
    fig.tight_layout()
    f = os.path.join(HERE, f"{tag}.png")
    fig.savefig(f, dpi=110)
    plt.close(fig)
    return [f]


def load(prefix):
    out = {}
    for d in sorted(glob.glob(os.path.join(OUT, prefix + "*.done"))):
        pf = d[:-5]
        m = meta_of(pf)
        if m and os.path.exists(os.path.join(pf, "thermo.dat")):
            out[m["name"]] = Run(pf, m)
    return out


def mesh_result(res_d):
    """Stage D verdict (1 mm vs 2 mm, trimmed) into MESH["result"]; returns (d1, d2)."""
    d1 = next((r for n, r in res_d.items() if "1mm" in n), None)
    d2 = next((r for n, r in res_d.items() if "2mm" in n), None)
    if d1 and d2:
        a = d1["ring_rop"] / d2["ring_rop"] - 1.0
        b = d1["Tf"][0] / d2["Tf"][0] - 1.0
        c = d1["rop"] / d2["rop"] - 1.0
        verdict = "PASS" if abs(a) <= 0.05 and abs(b) <= 0.05 and d1["steady"] and d2["steady"] else "FAIL"
        MESH["result"] = dict(text=f"ring ROP {a * 100:+.1f} %, T_fire {b * 100:+.2f} %, burner ROP "
                                   f"{c * 100:+.1f} % (1 mm vs 2 mm, trimmed; windows "
                                   f"{d1['W'][0]:.0f}–{d1['W'][1]:.0f} / {d2['W'][0]:.0f}–{d2['W'][1]:.0f} s; "
                                   f"steady {bool(d1['steady'])} / {bool(d2['steady'])})",
                              verdict=verdict)
    return d1, d2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=("A",))
    args = ap.parse_args()
    runs_a = load("A_")
    res_a = {n: analyze(r, "A") for n, r in runs_a.items()}
    a_lines, choice, why = stage_a_decision(res_a)
    if args.stage == "A":
        print("\n".join(a_lines))
        print(f"CHOICE: {choice} -- {why}")
        print("\n".join(run_table(res_a)))
        print("\n".join(check_table(res_a)))
        for n, r in res_a.items():
            print(n, "R per ring:", ", ".join(f"{i}:{x:.2f}" for i, x in enumerate(r["R_ring"]) if np.isfinite(x)))
        return
    out = []
    if res_a:
        out += ["## 1. Stage A: wall treatment", ""] + a_lines + ["",
                f"**Choice: {choice}**: {why}.", ""]
        out += run_table(res_a) + [""] + check_table(res_a) + [""]
        out += plots(res_a, "stageA") and ["![stageA](stageA.png)", ""]
    # the mesh verdict (Stage D) is needed by the Stage B scored table
    runs_d = load("D_")
    if runs_d:
        mesh_result({n: analyze(r, "D") for n, r in runs_d.items()})
    for stage, title in (("B_", "## 2. Stage B: scored runs"), ("D_", "## 3. Stage D: mesh check"),
                         ("C_", "## 4. Stage C: sensitivities")):
        runs = load(stage)
        if not runs:
            continue
        res = {n: analyze(r, stage[0]) for n, r in runs.items()}
        out += [title, ""] + run_table(res) + [""] + check_table(res) + [""]
        if stage in ("B_", "C_"):
            for n, r in res.items():
                hk = r["run"].m["anchor"]
                hn = hand_numbers(hk, r["run"].m["T_nozzle"], r["run"].m["mode"], r["run"].m["entrained"],
                                  r["z_feet"])
                mesh = MESH.get("result")
                tag = (" (scored)" if hk == "M" else " (bracket implication)") if stage == "B_" \
                    else " (sensitivity, verdict table)"
                out += [f"### {n}{tag}", ""]
                if stage == "C_":
                    mesh = None
                out += scored_table(r, hn, mesh if hk == "M" else None) + [""]
                zz = [0.01, 0.02, 0.05, 0.1, 0.15, 0.2]
                out.append("Wall profile r_wall(z) [mm], simulation | hand: " + "; ".join(
                    f"z {z * 1e3:.0f}: {np.interp(z, r['zg'], r['rw']) * 1e3 if z <= r['zg'][-1] else float('nan'):.0f}"
                    f" | {np.interp(z, hn['prof']['z'], hn['prof']['r_wall']) * 1e3 if hn else float('nan'):.0f}"
                    for z in zz))
                out.append("")
        if stage == "D_":
            d1, d2 = mesh_result(res)
            if d1 and d2:
                out += [f"Mesh: {MESH['result']['text']} → **{MESH['result']['verdict']}**; pinned share "
                        f"{d2['sh'].get('pinned_cols_pinned', 0):.2f} (2 mm) vs "
                        f"{d1['sh'].get('pinned_cols_pinned', 0):.2f} (1 mm); unset-face share "
                        f"{d2['sh'].get('pinned_cols_face_unset', 0):.2f} vs "
                        f"{d1['sh'].get('pinned_cols_face_unset', 0):.2f}.", ""]
        f = plots(res, "stage" + stage[0])
        if f:
            out += [f"![stage{stage[0]}](stage{stage[0]}.png)", ""]
    path = os.path.join(HERE, "RESULTS.md")
    txt = open(path).read()
    head = txt.split(MARK_R)[0]
    tail = txt.split(MARK_D, 1)[1] if MARK_D in txt else "\n"
    open(path, "w").write(head + MARK_R + "\n\n" + "\n".join(out) + "\n" + MARK_D + tail)
    print("\n".join(out))


MESH = {}

if __name__ == "__main__":
    main()
