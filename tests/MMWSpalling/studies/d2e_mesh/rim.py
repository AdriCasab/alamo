#!/usr/bin/env python3
"""D2e rim reconstruction from <pf>_removal_events.csv and <pf>/thermo.dat.

  rim.py RUN_PREFIX [RUN_PREFIX ...] [--blocks 0 50 100 ...] [--check-dt 10]

Per run:
  * rebuilds the per-column mask-top staircase k_top(t): every column starts
    at the domain top (nz - 1); each event sets k = k_top(event) - n_voided
    (and the event's k_top is checked against the running value);
  * applies FeetDescent's pad rule to it (annulus [28, 40] mm, 3 equal sectors
    of the observed azimuth range, nearest rank ceil(q n) - 1 with q = 0.9 on
    z_face = z_lo + (k_top + 1) dz, mean over the pads) and reports the
    reconstructed foot_z against the thermo foot_z at every check time
    (acceptance <= 0.5 mm; the worst is reported), at the thermo rows' own
    times; under foot_rule = mean (P-c) the rule is the annulus mean;
  * per block (default 50 s): the burner rate (least-squares fit of
    nozzle_z over the block, t > 0) and the annulus-mean face recession
    (every column in 28-40 mm), both mm/min; face recession per 2 mm radial band over
    28-60 mm, summed h_applied (the phi shift) per column, averaged over the
    band's columns, split into spall (regimes 1-3) and clip (regime 4), in
    mm/min; the share of each band's columns standing at their pad height
    at the block end, and the band of each pad's rank-selected column; the
    removal rate in the 40-60 mm band (sum h_applied dx dy, x4 for the full
    circle) in cm^3/s.
Geometry (dz, nxy, nz, lz) comes from <pf>.done. Event/thermo alignment: an
event at time t_e is applied to the thermo row with time t if t_e <= t - dt/2
(the log stamps a step's removal with its start time; thermo stamps its end).
"""
import argparse
import ast
import math
import os

import numpy as np

R_IN, R_OUT, NPADS, QUANT = 0.028, 0.040, 3, 0.9
BANDS = np.arange(0.028, 0.0601, 0.002)          # 2 mm bands, 28-60 mm
FLANK = (0.040, 0.060)


def meta_of(pf):
    return ast.literal_eval(open(pf + ".done").read())


def load_events(pf):
    a = np.loadtxt(pf + "_removal_events.csv", delimiter=",", skiprows=1, ndmin=2)
    # time,col_i,col_j,regime,k_top,T_top,a_f,Sp_top,h_scan,h_applied,n_voided
    o = np.argsort(a[:, 0], kind="stable")
    return a[o]


def load_thermo(pf, cols):
    with open(os.path.join(pf, "thermo.dat")) as f:
        head = f.readline().split()
    a = np.loadtxt(os.path.join(pf, "thermo.dat"), skiprows=1, ndmin=2)
    return {c: a[:, head.index(c)] for c in cols}


class Geometry:
    def __init__(self, m):
        self.dz, self.n, self.nz, self.lz = m["dz"], m["nxy"], m["nz"], m["lz"]
        self.dx = m["lxy"] / m["nxy"]
        c = (np.arange(self.n) + 0.5) * self.dx
        X, Y = np.meshgrid(c, c, indexing="ij")          # [i, j]
        self.r = np.hypot(X, Y)
        self.az = np.arctan2(Y, X)
        ann = (self.r * self.r >= R_IN ** 2) & (self.r * self.r <= R_OUT ** 2)
        self.ann = ann
        a = self.az[ann]
        amin, amax = a.min(), a.max()
        pad = np.zeros(self.r.shape, int)
        if amax > amin:
            pad = np.minimum(NPADS - 1, np.floor((self.az - amin) / (amax - amin) * NPADS).astype(int))
        self.pad = np.where(ann, pad, -1)
        self.band = np.digitize(self.r, BANDS) - 1        # band b = [BANDS[b], BANDS[b+1])
        self.band[(self.r < BANDS[0]) | (self.r >= BANDS[-1])] = -1
        self.nb = len(BANDS) - 1
        self.rule = "mean" if "foot_rule = mean" in m.get("switch", "") else "pads"

    def zface(self, k):
        return (k + 1) * self.dz                          # z_lo = 0

    def pads(self, k):
        """(foot_z, per-pad height, per-pad rank-selected (i, j) columns).
        foot_rule = mean (D2e P-c): one pad, the mean over the annulus."""
        z = self.zface(k)
        if self.rule == "mean":
            h = float(np.mean(z[self.ann]))
            return h, [h] * NPADS, [(-1, -1)] * NPADS
        hs, sel = [], []
        for p in range(NPADS):
            m = self.pad == p
            zi = z[m]
            idx = np.argsort(zi, kind="stable")
            rk = max(0, min(int(math.ceil(QUANT * zi.size)) - 1, zi.size - 1))
            hs.append(zi[idx[rk]])
            ii, jj = np.nonzero(m)
            sel.append((ii[idx[rk]], jj[idx[rk]]))
        return float(np.mean(hs)), hs, sel


def analyse(pf, edges, check_dt=10.0):
    m = meta_of(pf)
    g = Geometry(m)
    dt = m["dt"]
    ev = load_events(pf)
    th = load_thermo(pf, ["time", "foot_z", "nozzle_z"])
    t_end = float(th["time"][-1])
    # a run stops at the last step before stop_time (e.g. 299.904 s): keep
    # an edge within one thermo interval of t_end, clipped to t_end
    edges = [min(e, t_end) for e in edges if e <= t_end + 0.5]
    k = np.full((g.n, g.n), g.nz - 1, int)
    rec_sp = np.zeros((g.n, g.n))
    rec_cl = np.zeros((g.n, g.n))
    # check at the thermo rows nearest to every check_dt (their own times:
    # rows fall every 0.096 s at 16 ms, not on the 10 s marks)
    checks = np.unique([th["time"][int(np.argmin(np.abs(th["time"] - t)))]
                        for t in np.arange(check_dt, t_end + 1e-9, check_dt)])
    snap_t = sorted(set(list(checks) + [e for e in edges if e > 0]))
    snaps = {}
    bad_k = 0
    ie = 0
    ne = ev.shape[0]
    for t in snap_t:
        cut = t - 0.5 * dt
        while ie < ne and ev[ie, 0] <= cut:
            i, j, reg = int(ev[ie, 1]), int(ev[ie, 2]), int(ev[ie, 3])
            if int(ev[ie, 4]) != k[i, j]:
                bad_k += 1
            k[i, j] = int(ev[ie, 4]) - int(ev[ie, 10])
            if reg == 4:
                rec_cl[i, j] += ev[ie, 9]
            else:
                rec_sp[i, j] += ev[ie, 9]
            ie += 1
        snaps[round(t, 6)] = (k.copy(), rec_sp.copy(), rec_cl.copy())
    # foot_z check
    worst, worst_t, rows = 0.0, None, 0
    for t in checks:
        j = int(np.argmin(np.abs(th["time"] - t)))
        fz, _, _ = g.pads(snaps[round(t, 6)][0])
        d = abs(fz - th["foot_z"][j])
        rows += 1
        if d > worst:
            worst, worst_t = d, t
    # blocks
    blocks = []
    prev = (np.full((g.n, g.n), g.nz - 1, int), np.zeros((g.n, g.n)), np.zeros((g.n, g.n)))
    for a, b in zip(edges[:-1], edges[1:]):
        kb, sb, cb = snaps[round(b, 6)]
        if a > 0:
            _, sa, ca = snaps[round(a, 6)]
        else:
            sa, ca = prev[1], prev[2]
        mins = (b - a) / 60.0
        sp_band, cl_band, at_h = [], [], []
        _, hs, sel = g.pads(kb)
        z = g.zface(kb)
        padh = np.full(z.shape, np.nan)
        for p in range(NPADS):
            padh[g.pad == p] = hs[p]
        for bi in range(g.nb):
            mb = g.band == bi
            sp_band.append(float(np.mean(sb[mb] - sa[mb])) * 1e3 / mins)
            cl_band.append(float(np.mean(cb[mb] - ca[mb])) * 1e3 / mins)
            ma = mb & g.ann
            at_h.append(float(np.mean(np.isclose(z[ma], padh[ma]))) if ma.any() else float("nan"))
        sel_band = [int(g.band[i, j]) if i >= 0 else -1 for i, j in sel]
        ann_mean = float(np.mean((sb + cb - sa - ca)[g.ann])) * 1e3 / mins
        st = (th["time"] >= max(a, 1e-12)) & (th["time"] <= b)
        burner = -np.polyfit(th["time"][st], th["nozzle_z"][st], 1)[0] * 6e4 if st.sum() > 2 else np.nan
        mf = (g.r >= FLANK[0]) & (g.r < FLANK[1])
        vol = float(np.sum((sb + cb - sa - ca)[mf])) * g.dx * g.dx * 4.0 / (b - a) * 1e6
        blocks.append(dict(a=a, b=b, sp=sp_band, cl=cl_band, at_h=at_h, sel_band=sel_band,
                           flank_cm3s=vol, pad_h=hs, ann_mean=ann_mean, burner=float(burner)))
    return dict(name=m["name"], dz=g.dz, worst=worst, worst_t=worst_t, checks=rows, bad_k=bad_k,
                blocks=blocks, geom=g, n_events=ne)


def band_label(bi):
    return f"{BANDS[bi] * 1e3:.0f}–{BANDS[bi + 1] * 1e3:.0f}"


def report(res):
    g = res["geom"]
    L = [f"### {res['name']} (dz {res['dz'] * 1e3:g} mm, {res['n_events']} events)",
         f"foot_z reconstruction: worst |Δ| = {res['worst'] * 1e3:.3f} mm at t = {res['worst_t']} s "
         f"over {res['checks']} check times ({'PASS' if res['worst'] <= 5e-4 + 1e-12 else 'FAIL'}, "
         f"≤ 0.5 mm); event k_top vs running staircase mismatches: {res['bad_k']}",
         "",
         "Rim (38–40 mm band) and pad, per block:",
         "",
         "| block | rim spall | rim clip | rim total [mm/min] | annulus mean | burner | "
         "pad-setting columns in band | share of the 38–40 band at H | 40–60 mm removal [cm³/s] |",
         "|---|---|---|---|---|---|---|---|---|"]
    rim = int(np.nonzero(np.isclose(BANDS[:-1], 0.038))[0][0])
    for bl in res["blocks"]:
        sb = ", ".join(band_label(b) if b >= 0 else "—" for b in bl["sel_band"])
        L.append(f"| {bl['a']:.0f}–{bl['b']:.0f} s | {bl['sp'][rim]:.1f} | {bl['cl'][rim]:.1f} | "
                 f"{bl['sp'][rim] + bl['cl'][rim]:.1f} | {bl['ann_mean']:.1f} | {bl['burner']:.1f} | "
                 f"{sb} | {bl['at_h'][rim]:.2f} | "
                 f"{bl['flank_cm3s']:.2f} |")
    L += ["", "Recession by 2 mm band, spall + clip [mm/min] (pad shape = the 28–40 mm columns):", ""]
    hdr = "| block | " + " | ".join(band_label(b) for b in range(g.nb)) + " |"
    L += [hdr, "|---|" + "---|" * g.nb]
    for bl in res["blocks"]:
        L.append(f"| {bl['a']:.0f}–{bl['b']:.0f} s | "
                 + " | ".join(f"{bl['sp'][b]:.0f}+{bl['cl'][b]:.0f}" for b in range(g.nb)) + " |")
    L += ["", "Share of each annulus band's columns standing at their pad height (block end):", ""]
    ab = [b for b in range(g.nb) if BANDS[b] < R_OUT - 1e-9]
    L += ["| block | " + " | ".join(band_label(b) for b in ab) + " |", "|---|" + "---|" * len(ab)]
    for bl in res["blocks"]:
        L.append(f"| {bl['a']:.0f}–{bl['b']:.0f} s | " + " | ".join(f"{bl['at_h'][b]:.2f}" for b in ab) + " |")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--blocks", nargs="*", type=float, default=[0, 50, 100, 150, 200, 250, 300])
    ap.add_argument("--check-dt", type=float, default=10.0)
    a = ap.parse_args()
    for pf in a.runs:
        print(report(analyse(os.path.abspath(pf), a.blocks, a.check_dt)))
        print()


if __name__ == "__main__":
    main()
