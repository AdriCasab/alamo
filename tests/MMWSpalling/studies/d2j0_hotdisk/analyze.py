#!/usr/bin/env python3
"""D2j-0 analysis: the metrics frozen in PREDICTIONS.md, P1–P4, the
sharp-vs-tapered discriminator, and d2j0_vprofile.png. Writes analysis.md.

Standalone, like run.py — no import of the d2c→d2i chain.
"""
import argparse
import importlib.util
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_sp = importlib.util.spec_from_file_location("d2j0_run", os.path.join(HERE, "run.py"))
R = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(R)

OUT = R.OUT
WIN = (50.0, 150.0)                       # the metric window
EDGES = np.arange(0.0, 0.0401, 0.002)     # 2 mm annuli, 0-40 mm
CEN = 0.5 * (EDGES[:-1] + EDGES[1:])
R_STALL_MAX = 0.030                       # f_stall counts columns with r < 30 mm
SHARP = ["S4", "S2", "S1", "S05"]
TAPER = ["T4", "T2", "T1"]


def have(n):
    return os.path.exists(os.path.join(OUT, n + ".done"))


def geom(n):
    m = R.meta_of(os.path.join(OUT, n))
    c = (np.arange(m["nxy"]) + 0.5) * (m["lxy"] / m["nxy"])
    X, Y = np.meshgrid(c, c, indexing="ij")
    return m, np.hypot(X, Y)


def recession(n, win):
    """Per-column recession rate in m/h over `win`, and the column radii."""
    m, r = geom(n)
    a = np.loadtxt(os.path.join(OUT, n + "_removal_events.csv"), delimiter=",",
                   skiprows=1, ndmin=2)
    s = (a[:, 0] >= win[0]) & (a[:, 0] <= win[1])
    v = np.zeros_like(r)
    np.add.at(v, (a[s, 1].astype(int), a[s, 2].astype(int)), a[s, 9])
    return m, r, v / (win[1] - win[0]) * 3600.0


def binned(r, v):
    b = np.digitize(r, EDGES) - 1
    return np.array([v[b == i].mean() if np.any(b == i) else np.nan
                     for i in range(len(CEN))])


def cross(prof, level):
    """Outermost radius where the binned profile falls through `level`,
    by linear interpolation between bin centres."""
    ok = np.isfinite(prof)
    x, y = CEN[ok], prof[ok]
    hit = np.nonzero((y[:-1] >= level) & (y[1:] < level))[0]
    if not len(hit):
        return np.nan
    i = hit[-1]
    if y[i] == y[i + 1]:
        return float(x[i])
    return float(x[i] + (x[i + 1] - x[i]) * (y[i] - level) / (y[i] - y[i + 1]))


def metrics(n):
    m, r, v = recession(n, WIN)
    prof = binned(r, v)
    v_c = float(v[r < 0.006].mean())
    d = dict(name=n, dz=m["dz"], set=m["set"], v_c=v_c, prof=prof,
             r_half=cross(prof, 0.5 * v_c), w_edge=np.nan,
             f_stall=float(np.mean(v[r < R_STALL_MAX] < 0.1 * v_c)),
             wall=m["wall_s"] / 60.0, t_end=m["t_end"], ncol=int((r < R_STALL_MAX).sum()))
    r90, r10 = cross(prof, 0.9 * v_c), cross(prof, 0.1 * v_c)
    d["r90"], d["r10"] = r90, r10
    d["w_edge"] = r10 - r90
    for t in (50.0, 100.0, 150.0):
        d[f"rsi{int(t)}"] = r_stall_inner(n, t, v_c)
    return d


def r_stall_inner(n, t, v_c):
    """Inner edge of the contiguous run of bins reaching r = 30 mm with mean
    v < 0.1 v_c over the trailing window [t-50, t]. NaN if there is none."""
    _, r, v = recession(n, (t - 50.0, t))
    prof = binned(r, v)
    i30 = int(np.searchsorted(CEN, R_STALL_MAX)) - 1
    if not (np.isfinite(prof[i30]) and prof[i30] < 0.1 * v_c):
        return np.nan
    i = i30
    while i - 1 >= 0 and np.isfinite(prof[i - 1]) and prof[i - 1] < 0.1 * v_c:
        i -= 1
    return float(EDGES[i])


def table(rows):
    L = ["| run | dz | v_c [m/h] | r_half [mm] | r(0.9) [mm] | r(0.1) [mm] | **w_edge** [mm] | "
         "w_edge/dz | f_stall | r_stall_inner 50/100/150 s [mm] | wall |",
         "|---|---|---|---|---|---|---|---|---|---|---|"]
    for d in rows:
        rs = " / ".join("none" if not np.isfinite(d[f"rsi{t}"]) else f"{d[f'rsi{t}']*1e3:.0f}"
                        for t in (50, 100, 150))
        L.append(f"| {d['name']} | {d['dz']*1e3:g} mm | {d['v_c']:.3f} | {d['r_half']*1e3:.2f} | "
                 f"{d['r90']*1e3:.2f} | {d['r10']*1e3:.2f} | **{d['w_edge']*1e3:.2f}** | "
                 f"{d['w_edge']/d['dz']:.2f} | {d['f_stall']*100:.1f} % | {rs} | {d['wall']:.1f} min |")
    return L


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quiet", action="store_true")
    ap.parse_args()
    L = [f"# D2j-0 analysis (window {WIN[0]:.0f}–{WIN[1]:.0f} s, 2 mm annuli)", ""]
    sh = [metrics(n) for n in SHARP if have(n)]
    ta = [metrics(n) for n in TAPER if have(n)]

    L += ["## Sharp set (h = 700 inside r = 30 mm, drops to 0 across one cell)", ""] + table(sh)
    if ta:
        L += ["", "## Tapered set (raised cosine, 700 at 27.5 mm to 0 at 32.5 mm)", ""] + table(ta)

    L += ["", "## v(r)/v_c by 2 mm annulus", "",
          "| r [mm] | " + " | ".join(d["name"] for d in sh + ta) + " |",
          "|---|" + "---|" * len(sh + ta)]
    for i, rc in enumerate(CEN):
        vals = " | ".join("—" if not np.isfinite(d["prof"][i]) else f"{d['prof'][i]/d['v_c']:.3f}"
                          for d in sh + ta)
        L.append(f"| {rc*1e3:.0f} | {vals} |")

    # --- P1-P4 ---
    L += ["", "## Predictions scored", ""]
    prim = [d for d in sh if d["dz"] >= 1e-3]
    vc = [d["v_c"] for d in prim]
    spread = (max(vc) - min(vc)) / np.mean(vc) if vc else float("nan")
    p1 = spread <= 0.05
    L.append(f"- **P1** v_c across 4/2/1 mm = " + ", ".join(f"{x:.3f}" for x in vc)
             + f" m/h, spread **{spread*100:.1f} %** (limit 5 %) — **{'HOLDS' if p1 else 'FAILS'}**.")
    if not p1:
        L.append("  - P1 failing means the configuration is wrong, not the edge. P2–P4 are not "
                 "scored; see PREDICTIONS.md.")
    fs = [d["f_stall"] for d in prim]
    any_stall = all(d["f_stall"] > 0 for d in prim)
    grows = all(fs[i] <= fs[i + 1] for i in range(len(fs) - 1))  # 4 mm -> 1 mm order
    p2 = any_stall and grows
    L.append(f"- **P2** stalled annulus at every mesh: {'yes' if any_stall else 'NO'}; f_stall "
             f"4→1 mm = " + ", ".join(f"{x*100:.1f} %" for x in fs)
             + f" ({'grows' if grows else 'does NOT grow'} as dx falls) — "
             f"**{'HOLDS' if p2 else 'FAILS'}**.")
    w2 = next((d["w_edge"] for d in prim if abs(d["dz"] - 2e-3) < 1e-9), np.nan)
    w1 = next((d["w_edge"] for d in prim if abs(d["dz"] - 1e-3) < 1e-9), np.nan)
    ch = abs(w1 / w2 - 1.0) if np.isfinite(w1) and np.isfinite(w2) and w2 else np.nan
    p3 = np.isfinite(ch) and ch > 0.20
    L.append(f"- **P3** w_edge 2 mm → 1 mm: {w2*1e3:.2f} → {w1*1e3:.2f} mm, change "
             f"**{ch*100:.1f} %** (> 20 % = not converged) — **{'HOLDS' if p3 else 'FAILS'}**.")
    for d in prim:
        pass
    d1 = next((d for d in prim if abs(d["dz"] - 1e-3) < 1e-9), None)
    d4 = next((d for d in prim if abs(d["dz"] - 4e-3) < 1e-9), None)

    def inward(d):
        a, b = d["rsi100"], d["rsi150"]
        if not np.isfinite(b):
            return None
        if not np.isfinite(a):
            return True                      # appeared where there was none
        return b < a - 1e-9

    i1, i4 = (inward(d1) if d1 else None), (inward(d4) if d4 else None)
    p4 = bool(i1) and not bool(i4)
    L.append(f"- **P4** r_stall_inner moves inward 100 → 150 s: 1 mm {i1}, 4 mm {i4} — "
             f"**{'HOLDS' if p4 else 'FAILS'}**.")

    if ta:
        tw = {d["dz"]: d["w_edge"] for d in ta}
        sw = {d["dz"]: d["w_edge"] for d in sh}
        tch = abs(tw.get(1e-3, np.nan) / tw.get(2e-3, np.nan) - 1.0)
        L += ["", "## Sharp-vs-tapered discriminator", "",
              f"- sharp w_edge 2 → 1 mm: {sw.get(2e-3, np.nan)*1e3:.2f} → "
              f"{sw.get(1e-3, np.nan)*1e3:.2f} mm ({ch*100:.1f} %)",
              f"- tapered w_edge 2 → 1 mm: {tw.get(2e-3, np.nan)*1e3:.2f} → "
              f"{tw.get(1e-3, np.nan)*1e3:.2f} mm ({tch*100:.1f} %)",
              f"- tapered f_stall 4→1 mm: " + ", ".join(f"{d['f_stall']*100:.1f} %" for d in ta)]

    txt = "\n".join(L)
    open(os.path.join(HERE, "analysis.md"), "w").write(txt + "\n")
    print(txt)
    fig(sh, ta)
    print(f"\nS05 condition: w_edge moves {ch*100:.1f} % from 2 mm to 1 mm -> "
          f"{'RUN S05' if p3 else 'do not run S05'}")


def fig(sh, ta):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    f, ax = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)
    for a, rows, ttl, redge in ((ax[0], sh, "sharp: h = 700 inside r = 30 mm", [30.0]),
                                (ax[1], ta, "tapered: cosine 27.5 → 32.5 mm", [27.5, 32.5])):
        for d in rows:
            a.plot(CEN * 1e3, d["prof"] / d["v_c"], "o-", ms=4,
                   label=f"{d['name']} (dz {d['dz']*1e3:g} mm)")
        for x in redge:
            a.axvline(x, color="k", ls=":", lw=0.8)
        a.axhline(0.1, color="r", lw=0.6, ls="--")
        a.axhline(0.9, color="r", lw=0.6, ls="--")
        a.set(xlabel="r [mm]", title=ttl, xlim=(0, 40))
        a.grid(alpha=0.3)
        a.legend(fontsize=8)
    ax[0].set_ylabel("v(r) / v_c")
    f.suptitle("D2j-0: recession profile at a lateral hot/cold edge, 50–150 s", y=0.99)
    f.tight_layout()
    f.savefig(os.path.join(HERE, "d2j0_vprofile.png"), dpi=110)
    plt.close(f)
    print("wrote d2j0_vprofile.png")


if __name__ == "__main__":
    main()
