#!/usr/bin/env python3
"""D2j-0b analysis: the metrics fixed in CRITERION.md (v_c, per-band t_stop /
t_cross / depth-below-plane at stop, disk-mean rate by 50 s block and the mesh
gap, r_front), and the P1-P4 scoring. Writes analysis.md and d2j0b_bands.png.

  analyze.py            every finished run
  analyze.py --quiet
"""
import argparse
import importlib.util
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_sp = importlib.util.spec_from_file_location("d2j0b_run", os.path.join(HERE, "run.py"))
R = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(R)

OUT = R.OUT
EDGES = np.arange(0.018, 0.0301, 0.002)      # 2 mm bands 18-30 mm
BLOCKS = [(50.0, 100.0), (100.0, 150.0), (150.0, 200.0), (200.0, 250.0)]
FRONT_T = (100.0, 150.0, 200.0, 250.0)
D2J0_DISKMEAN = {2e-3: 1.1331, 1e-3: 1.0142}  # D2j-0 sharp, 50-150 s, r < 30 mm [m/h]
V_C_REF = 1.58


def have(n):
    return os.path.exists(os.path.join(OUT, n + ".done"))


def load(n):
    m = R.J.meta_of(os.path.join(OUT, n))
    nxy, dz = m["nxy"], m["dz"]
    c = (np.arange(nxy) + 0.5) * dz
    X, Y = np.meshgrid(c, c, indexing="ij")
    r = np.hypot(X, Y)
    a = np.loadtxt(os.path.join(OUT, n + "_removal_events.csv"), delimiter=",",
                   skiprows=1, ndmin=2)
    t_end = m["t_end"]
    return m, r, a, t_end


def col_rate(r, a, win):
    s = (a[:, 0] >= win[0]) & (a[:, 0] <= win[1])
    v = np.zeros_like(r)
    np.add.at(v, (a[s, 1].astype(int), a[s, 2].astype(int)), a[s, 9])
    return v / (win[1] - win[0]) * 3600.0


def per_column_history(m, a, i, j):
    """(t, cumulative recession) for one column, cell-quantised face height."""
    s = (a[:, 1] == i) & (a[:, 2] == j)
    t = a[s, 0]
    cum = np.cumsum(a[s, 9])
    return t, cum


def band_timing(m, r, a, t_end, lo, hi):
    """Median over the band's columns of t_stop, t_cross, and depth below the
    plane at t_stop. t_cross = first time the (cell-quantised) face is above the
    plane; the face is z_top - (voided cells)*dz, plane z0 - feed t."""
    dz, lz, z0, feed = m["dz"], m["lz"], m["z0"], m["feed"]
    cols = np.argwhere((r >= lo) & (r < hi))
    t_stop, t_cross, dbelow = [], [], []
    for i, j in cols:
        s = (a[:, 1] == i) & (a[:, 2] == j)
        t = a[s, 0]
        nv = np.cumsum(a[s, 10])                 # cells voided so far -> face drops nv*dz
        if len(t) == 0:
            ts, tc, db = 0.0, 0.0, np.nan       # never fired: face at the surface from t = 0
            # crossing when z0 - feed t < lz  ->  t > (z0 - lz)/feed
            tc = (z0 - lz) / feed
            db = -(z0 - feed * ts - lz)
        else:
            ts = t[-1] if t[-1] < t_end - 5.0 else t_end
            # face height after each event; plane at that time
            zf = lz - nv * dz
            zn = z0 - feed * t
            above = np.nonzero(zf > zn)[0]
            if len(above):
                # crossing happened between events: the face was constant since the
                # previous event, so solve for the time the plane fell below it
                k = above[0]
                zf_prev = lz - (nv[k - 1] if k > 0 else 0.0) * dz
                tc = (z0 - zf_prev) / feed
            else:
                # still below at the last event; does the plane pass the final face before t_end?
                tc_final = (z0 - zf[-1]) / feed
                tc = tc_final if tc_final <= t_end else np.inf
            db = (z0 - feed * ts) - (lz - nv[-1] * dz)   # s at t_stop (>0 below the plane)
        t_stop.append(ts)
        t_cross.append(tc)
        dbelow.append(db)
    return (float(np.median(t_stop)), float(np.median(t_cross)), float(np.nanmedian(dbelow)),
            len(cols))


def r_front(r, a, t, v_c):
    v = col_rate(r, a, (t - 25.0, t))
    edges = np.arange(0.0, 0.0301, 0.002)
    cen = 0.5 * (edges[:-1] + edges[1:])
    b = np.digitize(r, edges) - 1
    prof = np.array([v[b == k].mean() if np.any(b == k) else np.nan for k in range(len(cen))])
    k = len(cen) - 1
    if not (np.isfinite(prof[k]) and prof[k] < 0.1 * v_c):
        return np.nan
    while k - 1 >= 0 and np.isfinite(prof[k - 1]) and prof[k - 1] < 0.1 * v_c:
        k -= 1
    return float(edges[k])


def metrics(n):
    m, r, a, t_end = load(n)
    v = col_rate(r, a, (50.0, min(250.0, t_end)))
    v_c = float(v[r < 0.006].mean())
    d = dict(name=n, case=m["case"], dz=m["dz"], t_end=t_end, v_c=v_c,
             disk50_150=float(col_rate(r, a, (50.0, 150.0))[r < 0.030].mean()),
             blocks={}, bands={}, front={})
    for b in BLOCKS:
        if b[1] <= t_end + 1e-6:
            d["blocks"][b] = float(col_rate(r, a, b)[r < 0.030].mean())
    for lo, hi in zip(EDGES[:-1], EDGES[1:]):
        d["bands"][(round(lo * 1e3), round(hi * 1e3))] = band_timing(m, r, a, t_end, lo, hi)
    for t in FRONT_T:
        if t <= t_end + 1e-6:
            d["front"][t] = r_front(r, a, t, v_c)
    return d


def fmt_t(x):
    return "∞" if not np.isfinite(x) else f"{x:.0f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quiet", action="store_true")
    ap.parse_args()
    names = [n for n in R.CASES if have(n)]
    D = {n: metrics(n) for n in names}
    L = ["# D2j-0b analysis (metrics per CRITERION.md)", ""]
    L += ["## Runs", "", "| run | dz | t_end | v_c [m/h] | disk-mean 50–150 s [m/h] | D2j-0 ref | "
          + " | ".join(f"R {int(b[0])}–{int(b[1])}" for b in BLOCKS) + " |",
          "|---|---|---|---|---|---|" + "---|" * len(BLOCKS)]
    for n in names:
        d = D[n]
        ref = D2J0_DISKMEAN[d["dz"]]
        L.append(f"| {n} | {d['dz']*1e3:g} mm | {d['t_end']:.0f} s | {d['v_c']:.3f} | "
                 f"{d['disk50_150']:.3f} | {ref:.3f} ({(d['disk50_150']/ref-1)*100:+.1f} %) | "
                 + " | ".join(f"{d['blocks'][b]:.3f}" if b in d["blocks"] else "—" for b in BLOCKS) + " |")
    # mesh gaps per case
    L += ["", "## Mesh gap g = R_1mm/R_2mm − 1 per 50 s block", "",
          "| case | " + " | ".join(f"{int(b[0])}–{int(b[1])}" for b in BLOCKS) + " |",
          "|---|" + "---|" * len(BLOCKS)]
    gaps = {}
    for case in ("C0", "C1", "C2"):
        n2, n1 = f"{case}_2mm", f"{case}_1mm"
        if n2 in D and n1 in D:
            g = {b: D[n1]["blocks"][b] / D[n2]["blocks"][b] - 1.0
                 for b in BLOCKS if b in D[n1]["blocks"] and b in D[n2]["blocks"]}
            gaps[case] = g
            L.append(f"| {case} | " + " | ".join(f"{g[b]*100:+.1f} %" if b in g else "—" for b in BLOCKS) + " |")
    # band timing
    L += ["", "## Band timing (median over the band's columns): t_stop / t_cross(own) / Δ_own = t_cross − t_stop / "
          "depth below plane at stop [mm] / Δ_outer = t_stop − t_cross(next band out)", ""]
    for n in names:
        d = D[n]
        L += [f"### {n}", "", "| band [mm] | n | t_stop | t_cross | Δ_own | below plane at stop | Δ_outer |",
              "|---|---|---|---|---|---|---|"]
        keys = sorted(d["bands"].keys(), reverse=True)
        outer_cross = None
        for k in keys:
            ts, tc, db, nc = d["bands"][k]
            d_own = tc - ts if np.isfinite(tc) else np.nan
            d_out = ts - outer_cross if outer_cross is not None and np.isfinite(outer_cross) else np.nan
            stopped = ts < d["t_end"] - 1e-6
            L.append(f"| {k[0]}–{k[1]} | {nc} | {ts:.0f}{'' if stopped else ' (firing)'} | {fmt_t(tc)} | "
                     f"{'—' if not np.isfinite(d_own) else f'{d_own:+.0f} s'} | "
                     f"{'—' if not np.isfinite(db) else f'{db*1e3:.0f}'} | "
                     f"{'—' if not np.isfinite(d_out) else f'{d_out:+.0f} s'} |")
            outer_cross = tc
        L.append("")
    L += ["## r_front(t) [mm] (inner edge of the stalled run reaching 30 mm, trailing 25 s, < 0.1 v_c)", "",
          "| run | " + " | ".join(f"{int(t)} s" for t in FRONT_T) + " |", "|---|" + "---|" * len(FRONT_T)]
    for n in names:
        d = D[n]
        L.append(f"| {n} | " + " | ".join(("none" if not np.isfinite(d["front"][t]) else f"{d['front'][t]*1e3:.0f}")
                                          if t in d["front"] else "—" for t in FRONT_T) + " |")

    # ---- scoring ----
    L += ["", "## Predictions scored", ""]
    p1 = all(abs(D[n]["v_c"] / V_C_REF - 1.0) <= 0.05 for n in names)
    p1b = all(abs(D[n]["disk50_150"] / D2J0_DISKMEAN[D[n]["dz"]] - 1.0) <= 0.03
              for n in names if D[n]["case"] == "C0")
    L.append(f"- **P1** v_c within 5 % of {V_C_REF} in every run: {'yes' if p1 else 'NO'}; C0 disk-mean 50–150 s "
             f"within 3 % of D2j-0: {'yes' if p1b else 'NO'} — **{'HOLDS' if p1 and p1b else 'FAILS'}**.")

    def stops_early(n):
        d = D[n]
        out = []
        for k, (ts, tc, db, nc) in d["bands"].items():
            if k[1] <= 28 and ts < d["t_end"] - 1e-6 and np.isfinite(tc) and tc - ts > 10.0:
                out.append(k)
        return out

    for case, want in (("C1", True), ("C0", False), ("C2", False)):
        for mesh in ("2mm", "1mm"):
            n = f"{case}_{mesh}"
            if n in D:
                e = stops_early(n)
                inner_stopped = [k for k, (ts, tc, db, nc) in D[n]["bands"].items()
                                 if k[1] <= 28 and ts < D[n]["t_end"] - 1e-6]
                L.append(f"- {n}: bands inside 28 mm that stopped: {inner_stopped or 'none'}; "
                         f"with Δ_own > 10 s: {e or 'none'}.")
    c1 = [f"C1_{m}" for m in ("2mm", "1mm") if f"C1_{m}" in D]
    c0 = [f"C0_{m}" for m in ("2mm", "1mm") if f"C0_{m}" in D]
    p2_c1 = bool(c1) and all(len(stops_early(n)) >= 2 for n in c1)
    p2_c0 = bool(c0) and all(not [k for k, (ts, tc, db, nc) in D[n]["bands"].items()
                                  if k[1] <= 28 and ts < D[n]["t_end"] - 1e-6] for n in c0)
    L.append(f"- **P2** C1 cascades (≥ 2 bands inside 28 mm stop with Δ_own > 10 s, both meshes): "
             f"{'yes' if p2_c1 else 'NO'}; C0 no inner band stops: {'yes' if p2_c0 else 'NO'} — "
             f"**{'HOLDS' if p2_c1 and p2_c0 else 'FAILS'}** (runs present: {c1 + c0}).")
    if "C1" in gaps and "C0" in gaps:
        g1 = [gaps["C1"][b] for b in BLOCKS if b in gaps["C1"]]
        g0 = [gaps["C0"][b] for b in BLOCKS if b in gaps["C0"]]
        mono = all(abs(g1[i + 1]) >= abs(g1[i]) - 0.005 for i in range(len(g1) - 1))
        p3 = mono and g1[-1] < -0.20 and all(abs(x - g0[0]) <= 0.04 for x in g0)
        L.append(f"- **P3** C1 gap grows monotonically to < −20 %: {mono and g1[-1] < -0.20} "
                 f"(last {g1[-1]*100:+.1f} %); C0 gap standing within ±4 pp of its first block: "
                 f"{all(abs(x - g0[0]) <= 0.04 for x in g0)} — **{'HOLDS' if p3 else 'FAILS'}**.")
    if "C2" in gaps:
        c2 = [f"C2_{m}" for m in ("2mm", "1mm") if f"C2_{m}" in D]
        g2 = [gaps["C2"][b] for b in BLOCKS if b in gaps["C2"]]
        no_casc = all(len(stops_early(n)) == 0 for n in c2)
        standing = all(abs(x - g2[0]) <= 0.04 for x in g2)
        L.append(f"- **P4** C2 no band stops with Δ_own > 10 s: {no_casc}; gap standing: {standing} "
                 f"(blocks {' '.join(f'{x*100:+.1f}%' for x in g2)}) — **{'HOLDS' if no_casc and standing else 'FAILS'}**.")
    if "C1_2mm" in D and "C1d_2mm" in D:
        dts = max(abs(D["C1_2mm"]["bands"][k][0] - D["C1d_2mm"]["bands"][k][0]) for k in D["C1_2mm"]["bands"])
        drr = max(abs(D["C1d_2mm"]["blocks"][b] / D["C1_2mm"]["blocks"][b] - 1.0)
                  for b in BLOCKS if b in D["C1_2mm"]["blocks"] and b in D["C1d_2mm"]["blocks"])
        L.append(f"- **dt** C1d vs C1 (2 mm): max |Δt_stop| {dts:.1f} s (limit 5), max block-rate change "
                 f"{drr*100:.2f} % (limit 2) — **{'HOLDS' if dts <= 5.0 and drr <= 0.02 else 'FAILS'}**.")

    txt = "\n".join(L)
    open(os.path.join(HERE, "analysis.md"), "w").write(txt + "\n")
    print(txt)
    fig(D)


def fig(D):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    names = list(D)
    f, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    for n in names:
        d = D[n]
        ks = sorted(d["bands"].keys())
        cen = [0.5 * (k[0] + k[1]) for k in ks]
        ax[0].plot(cen, [d["bands"][k][0] for k in ks], "o-", ms=4, label=f"{n} t_stop")
        ax[0].plot(cen, [min(d["bands"][k][1], 300.0) for k in ks], "x--", ms=4, color=ax[0].lines[-1].get_color())
    ax[0].set(xlabel="band centre r [mm]", ylabel="time [s]", title="t_stop (o) and own t_cross (x, capped 300)")
    ax[0].grid(alpha=0.3)
    ax[0].legend(fontsize=7)
    for n in names:
        d = D[n]
        bs = sorted(d["blocks"])
        ax[1].plot([0.5 * (b[0] + b[1]) for b in bs], [d["blocks"][b] for b in bs], "o-", ms=4, label=n)
    ax[1].set(xlabel="block centre [s]", ylabel="disk-mean rate r < 30 mm [m/h]", title="disk-mean recession by block")
    ax[1].grid(alpha=0.3)
    ax[1].legend(fontsize=7)
    f.suptitle("D2j-0b: hot disk under a descending exclusion plane", y=0.99)
    f.tight_layout()
    f.savefig(os.path.join(HERE, "d2j0b_bands.png"), dpi=110)
    plt.close(f)
    print("wrote d2j0b_bands.png")


if __name__ == "__main__":
    main()
