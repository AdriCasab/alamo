#!/usr/bin/env python3
"""D2e pair analysis (Goal 6): tables, the PREDICTIONS rows and figures.

  analyze_pairs.py [--out pairs.md] [--no-figs]

For every pair with both meshes done (output/<P>_2mm, output/<P>_1mm):
  * the matched window 150 s - min(t_end) with score.mesh (text) and its
    burner-ROP gap;
  * 50 s block ROP (run.py blocks) and the 1 mm / 2 mm gap per block;
  * rim.py per mesh (rim spall/clip, annulus mean, burner, 40-60 mm removal);
  * jet_r_reach, T_rec, s_c and jet_P_face at 100 / 200 / 300 s; max |ledger_err|.
Then each PREDICTIONS.md row is evaluated mechanically (held / refuted /
not run) and the decision rule is applied. Operational definitions fixed here,
before the pair outputs existed (the hashed PREDICTIONS.md wording is
unchanged):
  * P-a "early gap <= 3 %": |gap| <= 3 % in both the 0-50 and 50-100 s blocks;
  * P-a "burner slower than B0 at 2 mm": Pa_2mm block ROP below B0_2mm's in
    every block to 200 s;
  * P-a "rim = spall only": no regime-4 row in either Pa run;
  * P-b "no growth": gap(250-300) <= gap(50-100) + 3 points;
  * P-c "150-300 s gap <= 5 %": |matched-window gap| <= 5 %;
  * P-c "burner follows the annulus mean": at both meshes, over 150-300 s,
    |burner - annulus mean| < |burner - rim-band total| (rim.py);
  * B0 gap: matched-window gap >= +30 %.
Figures: d2e_blocks.png, d2e_rim.png, d2e_loop.png.
"""
import argparse
import importlib.util
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rim  # noqa: E402

_spec = importlib.util.spec_from_file_location("d2e_run", os.path.join(HERE, "run.py"))
E = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(E)
sys.path.insert(0, E.D2C)
import score as sc  # noqa: E402

OUT = E.OUT
EDGES = [0.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0]
PAIRS = ["B0", "Pa", "Pb", "Pc"]
RIM_I = int(np.nonzero(np.isclose(rim.BANDS[:-1], 0.038))[0][0])


def have(n):
    return os.path.exists(os.path.join(OUT, n + ".done"))


def thermo(pf):
    with open(os.path.join(pf, "thermo.dat")) as f:
        h = f.readline().split()
    a = np.loadtxt(os.path.join(pf, "thermo.dat"), skiprows=1, ndmin=2)
    return {c: a[:, i] for i, c in enumerate(h)}


def window_rop(th, a, b):
    s = (th["time"] >= a) & (th["time"] <= b) & (th["time"] > 0)
    return -float(np.polyfit(th["time"][s], th["nozzle_z"][s], 1)[0]) * 3600.0


def n_clip(pf):
    ev = rim.load_events(pf)
    return int(np.sum(ev[:, 3] == 4))


def pair_data(p):
    d = {"pair": p}
    for mesh in ("2mm", "1mm"):
        pf = os.path.join(OUT, f"{p}_{mesh}")
        th = thermo(pf)
        t_end = float(th["time"][-1])
        edges = [min(e, t_end) for e in EDGES if e <= t_end + 0.5]
        d[mesh] = dict(pf=pf, th=th, t_end=t_end, blocks=E.blocks(pf, edges)[0], edges=edges,
                       rim=rim.analyse(pf, edges), clip=n_clip(pf),
                       rim_late=rim.analyse(pf, [150.0, min(300.0, t_end)]) if t_end >= 250 else None,
                       led=float(np.max(np.abs(th["ledger_err"]))) if "ledger_err" in th else np.nan,
                       meta=rim.meta_of(pf))
    t1 = min(d["2mm"]["t_end"], d["1mm"]["t_end"])
    d["W"] = (150.0, t1)
    d["rop2"], d["rop1"] = (window_rop(d[m]["th"], 150.0, t1) for m in ("2mm", "1mm"))
    d["gap"] = d["rop1"] / d["rop2"] - 1.0
    n = min(len(d["2mm"]["blocks"]), len(d["1mm"]["blocks"]))
    d["bgap"] = [d["1mm"]["blocks"][i] / d["2mm"]["blocks"][i] - 1.0 for i in range(n)]
    d["mesh_txt"] = sc.mesh(d["2mm"]["pf"], d["1mm"]["pf"])
    return d


def at(th, c, t):
    return float(np.interp(t, th["time"], th[c])) if t <= th["time"][-1] + 1e-6 else np.nan


def pair_table(d):
    L = [f"### {d['pair']}: {E.SWITCH[d['pair']][1]}", ""]
    L.append(f"Matched window {d['W'][0]:.0f}–{d['W'][1]:.0f} s: 2 mm {d['rop2']:.3f} m/h, 1 mm "
             f"{d['rop1']:.3f} m/h, **gap {d['gap'] * 100:+.1f} %**.")
    L += ["", "| block | 2 mm ROP | 1 mm ROP | gap |", "|---|---|---|---|"]
    for i, g in enumerate(d["bgap"]):
        a, b = EDGES[i], EDGES[i + 1]
        L.append(f"| {a:.0f}–{b:.0f} s | {d['2mm']['blocks'][i]:.3f} | {d['1mm']['blocks'][i]:.3f} | "
                 f"{g * 100:+.1f} % |")
    L += ["", "| quantity | t | 2 mm | 1 mm |", "|---|---|---|---|"]
    for c, lab, f in (("jet_r_reach", "reach [mm]", 1e3), ("jet_T_rec", "T_rec [K]", 1.0),
                      ("jet_s_c", "s_c [mm]", 1e3), ("jet_P_face", "P_face [W]", 1.0)):
        for t in (100.0, 200.0, 300.0):
            v2, v1 = at(d["2mm"]["th"], c, t), at(d["1mm"]["th"], c, t)
            if np.isnan(v2) and np.isnan(v1):
                continue
            L.append(f"| {lab} | {t:.0f} s | {v2 * f:.1f} | {v1 * f:.1f} |")
    L += ["", "| mesh | rim spall + clip by block [mm/min] | annulus mean | burner | 40–60 mm [cm³/s] | "
          "regime-4 rows | max abs ledger_err |", "|---|---|---|---|---|---|---|"]
    for m in ("2mm", "1mm"):
        bl = d[m]["rim"]["blocks"]
        L.append(f"| {m} | " + ", ".join(f"{x['sp'][RIM_I]:.1f}+{x['cl'][RIM_I]:.1f}" for x in bl)
                 + " | " + ", ".join(f"{x['ann_mean']:.1f}" for x in bl)
                 + " | " + ", ".join(f"{x['burner']:.1f}" for x in bl)
                 + " | " + ", ".join(f"{x['flank_cm3s']:.2f}" for x in bl)
                 + f" | {d[m]['clip']} | {d[m]['led']:.1e} |")
    L += ["", "`score.mesh`:", "", d["mesh_txt"], ""]
    return L


def verdicts(D, gate_ok, ident_ok):
    rows, res = [], {}
    if "B0" in D:
        d = D["B0"]
        ok = d["gap"] >= 0.30
        res["B0_gap"] = ok
        rows.append(("B0", "2 mm byte-identical to R7b_2mm_dt16 to 300 s",
                     "see §4 identity table", "held" if ident_ok else "REFUTED"))
        rows.append(("B0", "1 mm step gate (blocks 2 %, s_c/T_rec 1 %)", "see §4 gate table",
                     "held" if gate_ok else "REFUTED"))
        rows.append(("B0", "gap over 150–300 s ≥ +30 %", f"{d['gap'] * 100:+.1f} %",
                     "held" if ok else "REFUTED"))
    if "Pa" in D:
        d = D["Pa"]
        e = d["bgap"][:2]
        ok1 = all(abs(g) <= 0.03 for g in e)
        rows.append(("P-a", "early gap (0–50, 50–100 s) ≤ 3 %",
                     ", ".join(f"{g * 100:+.1f} %" for g in e), "held" if ok1 else "REFUTED"))
        res["Pa_early"] = ok1
        if "B0" in D:
            b0 = D["B0"]["2mm"]["blocks"]
            pa = d["2mm"]["blocks"]
            ok2 = all(pa[i] < b0[i] for i in range(min(len(pa), len(b0), 4)))
            rows.append(("P-a", "burner slower than B0 at 2 mm (every block to 200 s)",
                         ", ".join(f"{pa[i]:.3f} vs {b0[i]:.3f}" for i in range(min(len(pa), 4))),
                         "held" if ok2 else "REFUTED"))
        ok3 = d["2mm"]["clip"] == 0 and d["1mm"]["clip"] == 0
        rows.append(("P-a", "rim = spall only (no regime-4 rows)",
                     f"{d['2mm']['clip']} / {d['1mm']['clip']} rows", "held" if ok3 else "REFUTED"))
    if "Pb" in D:
        d = D["Pb"]
        if len(d["bgap"]) >= 6:
            ok = d["bgap"][5] <= d["bgap"][1] + 0.03
            rows.append(("P-b", "gap(250–300) ≤ gap(50–100) + 3 points",
                         f"{d['bgap'][5] * 100:+.1f} % vs {d['bgap'][1] * 100:+.1f} % + 3",
                         "held" if ok else "REFUTED"))
            res["Pb_nogrowth"] = ok
    if "Pc" in D:
        d = D["Pc"]
        ok = abs(d["gap"]) <= 0.05
        rows.append(("P-c", "150–300 s gap ≤ 5 %", f"{d['gap'] * 100:+.1f} %", "held" if ok else "REFUTED"))
        fol, txt = True, []
        for m in ("2mm", "1mm"):
            x = d[m]["rim_late"]["blocks"][0]
            rimt = x["sp"][RIM_I] + x["cl"][RIM_I]
            f = abs(x["burner"] - x["ann_mean"]) < abs(x["burner"] - rimt)
            fol = fol and f
            txt.append(f"{m}: burner {x['burner']:.1f}, annulus {x['ann_mean']:.1f}, rim {rimt:.1f}")
        rows.append(("P-c", "burner follows the annulus mean (150–300 s, both meshes)", "; ".join(txt),
                     "held" if fol else "REFUTED"))
    for p, row in (("P-a", "Pa"), ("P-b", "Pb"), ("P-c", "Pc")):
        if row not in D:
            rows.append((p, "all", "not run", "not tested"))
    return rows, res


def figures(D):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for p, d in D.items():
        mids = [(EDGES[i] + EDGES[i + 1]) / 2 for i in range(len(d["bgap"]))]
        ax.plot(mids, [g * 100 for g in d["bgap"]], "o-", label=f"{p}: {E.SWITCH[p][1]}")
    ax.axhline(0, color="k", lw=0.6)
    ax.axhspan(-5, 5, color="g", alpha=0.1, label="±5 %")
    ax.set(xlabel="block centre [s]", ylabel="1 mm / 2 mm burner ROP − 1 [%]",
           title="mesh gap per 50 s block")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "d2e_blocks.png"), dpi=110)
    plt.close(fig)

    fig, axs = plt.subplots(1, len(D), figsize=(4.2 * len(D), 4), squeeze=False)
    xb = [(rim.BANDS[b] + rim.BANDS[b + 1]) / 2 * 1e3 for b in range(len(rim.BANDS) - 1)]
    for ax, (p, d) in zip(axs[0], D.items()):
        for m, st in (("2mm", "-"), ("1mm", "--")):
            bl = d[m]["rim"]["blocks"]
            sp = np.mean([x["sp"] for x in bl[1:]], axis=0)
            cl = np.mean([x["cl"] for x in bl[1:]], axis=0)
            ax.plot(xb, sp, "b" + st, label=f"{m} spall")
            ax.plot(xb, cl, "r" + st, label=f"{m} clip")
        ax.axvspan(28, 40, color="k", alpha=0.06)
        ax.set(xlabel="r [mm]", ylabel="recession [mm/min] (mean of blocks after 50 s)", title=p)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "d2e_rim.png"), dpi=110)
    plt.close(fig)

    fig, axs = plt.subplots(1, 3, figsize=(15, 4.4))
    for p, d in D.items():
        for m, st in (("2mm", "-"), ("1mm", "--")):
            th = d[m]["th"]
            s = th["time"] > 0
            axs[0].plot(th["time"][s], th["jet_T_rec"][s], st, label=f"{p} {m}")
            axs[1].plot(th["time"][s], th["jet_s_c"][s] * 1e3, st, label=f"{p} {m}")
            axs[2].plot(th["time"][s], th["jet_r_reach"][s] * 1e3, st, label=f"{p} {m}")
    for ax, yl in zip(axs, ("T_rec [K]", "s_c [mm]", "jet_r_reach [mm]")):
        ax.set(xlabel="t [s]", ylabel=yl)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "d2e_loop.png"), dpi=110)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "pairs.md"))
    ap.add_argument("--no-figs", action="store_true")
    a = ap.parse_args()
    D = {p: pair_data(p) for p in PAIRS if have(f"{p}_2mm") and have(f"{p}_1mm")}
    import io
    import contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        ident_ok = E.identity() if have("B0_2mm") else False
    ident_txt = buf.getvalue()
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        gate_ok = E.gate() if have("B0_1mm") else False
    gate_txt = buf.getvalue()
    L = ["## B0 identity (2 mm) and 1 mm step gate", "", ident_txt, gate_txt]
    for p, d in D.items():
        L += pair_table(d)
    rows, res = verdicts(D, gate_ok, ident_ok)
    L += ["## PREDICTIONS rows", "", "| pair | expectation | measurement | verdict |", "|---|---|---|---|"]
    L += [f"| {a} | {b} | {c} | {v} |" for a, b, c, v in rows]
    txt = "\n".join(L)
    open(a.out, "w").write(txt + "\n")
    print(txt)
    if D and not a.no_figs:
        figures(D)
        print("\nwrote d2e_blocks.png, d2e_rim.png, d2e_loop.png")


if __name__ == "__main__":
    main()
