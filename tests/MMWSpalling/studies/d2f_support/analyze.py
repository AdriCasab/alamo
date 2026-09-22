#!/usr/bin/env python3
"""D2f analysis (Goal 5): the Q10 pair, the q-sensitivity, the PREDICTIONS
rows, the decision outcome and the figures.

  analyze.py [--out analysis.md] [--no-figs]

Reuses studies/d2e_mesh/analyze_pairs.py (pair tables, loaded by path and
pointed at this study's output), d2e rim.py (unchanged; imported) and
d2c score.py (az_rwall / load_run, for the hole at the run end). Operational
definitions (fixed with PREDICTIONS.md, before the runs):
  * block / window ROP: least-squares slope of nozzle_z (t > 0), m/h;
  * E1: |1 mm / 2 mm - 1| <= 5 % over 150 s - min(t_end) AND in every 50 s
    block after 50 s (i.e. 50-100 ... 250-300 s);
  * E2: the three 2 mm ROPs over 150-300 s within +-5 % of their mean;
  * E3: ROP(q 0.8) > ROP(q 0.9) > ROP(q 1.0) at 2 mm;
  * E4: 2 mm ROP within +-27 % of the hand value; E5: s_c at 300 s within
    +-20 mm of the hand value;
  * E6: run.py --identity.
Hole at the run end (direction only): score.py's azimuth-mean Ø(z) (9
sectors) from the removal-event depth map at t_end, Ø at 25 and 50 mm and the
minimum to the feet depth.
"""
import argparse
import contextlib
import importlib.util
import io
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)


def load(name, path):
    sp = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


F = load("d2f_run", os.path.join(HERE, "run.py"))
E = F.E
OUT = F.OUT
A = load("d2e_analyze", os.path.join(STUD, "d2e_mesh", "analyze_pairs.py"))
A.E = E                               # the d2e run module carrying the Q switches
A.OUT = OUT
rim = A.rim
sc = A.sc
EDGES = A.EDGES
QS = {"Q10": 1.0, "Q09": 0.9, "Q08": 0.8}
_rim_analyse = rim.analyse


def rim_analyse_q(pf, edges, check_dt=10.0):
    """rim.analyse with the run's own foot_pad_quantile (rim.QUANT is 0.9;
    D2e output unchanged: the wrapper restores it)."""
    m = re.search(r"foot_pad_quantile = ([\d.]+)", rim.meta_of(pf).get("switch", ""))
    old = rim.QUANT
    rim.QUANT = float(m.group(1)) if m else old
    try:
        return _rim_analyse(pf, edges, check_dt)
    finally:
        rim.QUANT = old


rim.analyse = rim_analyse_q
D2E_OUT = os.path.join(STUD, "d2e_mesh", "output")


def have(n):
    return os.path.exists(os.path.join(OUT, n + ".done"))


def hand():
    """Hand numbers parsed from the hashed PREDICTIONS.md table."""
    rows = {}
    for line in open(os.path.join(HERE, "PREDICTIONS.md")):
        m = re.match(r"\| (1\.0|0\.9|0\.8) \| ([\d.]+) \| ([\d.]+) \| ([\d.]+) \| ([\d.]+) / ([\d.]+) \|", line)
        if m:
            rows[float(m.group(1))] = dict(R_Q=float(m.group(2)), steady=float(m.group(3)),
                                           rop=float(m.group(4)), sc=float(m.group(5)), Trec=float(m.group(6)))
    return rows


def hole_at_end(pf):
    run = sc.load_run(pf)
    th = run.th
    t = float(th["time"][-1])
    feet = run.lz - float(np.interp(t, th["time"], th["foot_z"]))
    d_w = run.depth(t)
    zg = np.arange(0.0, max(feet, run.dz), run.dz / 2.0)
    D = np.array([2.0 * sc.az_rwall(run, d_w, z) for z in zg]) * 1e3
    at = lambda z_mm: float(np.interp(z_mm * 1e-3, zg, D))
    return dict(t=t, d25=at(25.0), d50=at(50.0), dmin=float(D.min()), zmin=float(zg[int(np.argmin(D))]) * 1e3,
                feet=feet * 1e3)


def q_row(name):
    pf = os.path.join(OUT, name)
    th = A.thermo(pf)
    t_end = float(th["time"][-1])
    rop = A.window_rop(th, 150.0, t_end)
    s = (th["time"] >= 150.0) & (th["time"] <= t_end)
    r = rim.analyse(pf, [150.0, t_end])
    bl = r["blocks"][0]
    return dict(name=name, t_end=t_end, rop=rop, sc=A.at(th, "jet_s_c", t_end), Trec=A.at(th, "jet_T_rec", t_end),
                carry=float(np.mean(th["foot_carry_cols"][s])) if "foot_carry_cols" in th else np.nan,
                stall=float(np.max(th["foot_stall_time"])) if "foot_stall_time" in th else np.nan,
                sp=bl["sp"], cl=bl["cl"], ann=bl["ann_mean"], burner=bl["burner"],
                sel=bl["sel_band"], worst=r["worst"], hole=hole_at_end(pf), meta=rim.meta_of(pf))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "analysis.md"))
    ap.add_argument("--no-figs", action="store_true")
    a = ap.parse_args()
    H = hand()
    L = []
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        ident = F.identity() if have("Q09_2mm") else False
    L += ["## Harness identity (E6)", "", buf.getvalue()]
    # Q10 pair (and Q09 if both meshes exist)
    pairs = {p: A.pair_data(p) for p in QS if have(f"{p}_2mm") and have(f"{p}_1mm")}
    for p, d in pairs.items():
        L += A.pair_table(d)
    # q-sensitivity at 2 mm
    Qr = {p: q_row(f"{p}_2mm") for p in QS if have(f"{p}_2mm")}
    nb = len(rim.BANDS) - 1
    ann_b = [b for b in range(nb) if rim.BANDS[b] < rim.R_OUT - 1e-9]
    L += ["## q-sensitivity (2 mm, 150 s – end)", "",
          "| q | ROP [m/h] | hand | s_c / T_rec at end | hand s_c / T_rec | mean foot_carry_cols | max stall [s] | "
          "burner / annulus mean [mm/min] | pad-setting bands | foot_z recon. |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for p, x in Qr.items():
        q = QS[p]
        h = H.get(q, {})
        sb = ", ".join(rim.band_label(b) if b >= 0 else "—" for b in x["sel"])
        L.append(f"| {q:.1f} | {x['rop']:.3f} | {h.get('rop', float('nan')):.3f} | {x['sc'] * 1e3:.1f} / "
                 f"{x['Trec']:.0f} | {h.get('sc', float('nan')):.1f} / {h.get('Trec', float('nan')):.0f} | "
                 f"{x['carry']:.1f} | {x['stall']:.1f} | {x['burner']:.1f} / {x['ann']:.1f} | {sb} | "
                 f"{x['worst'] * 1e3:.3f} mm |")
    L += ["", "Rim-zone spall by 2 mm band over 150 s – end [mm/min] (no clip in any D2f run):", "",
          "| q | " + " | ".join(rim.band_label(b) for b in ann_b) + " |", "|---|" + "---|" * len(ann_b)]
    for p, x in Qr.items():
        L.append(f"| {QS[p]:.1f} | " + " | ".join(f"{x['sp'][b] + x['cl'][b]:.1f}" for b in ann_b) + " |")
    L += ["", "Hole at the run end (direction only; azimuth-mean Ø, score.py):", "",
          "| q | t | Ø 25 mm | Ø 50 mm | min Ø to feet depth (at depth) | feet depth |", "|---|---|---|---|---|---|"]
    for p, x in Qr.items():
        o = x["hole"]
        L.append(f"| {QS[p]:.1f} | {o['t']:.0f} s | {o['d25']:.1f} | {o['d50']:.1f} | {o['dmin']:.1f} "
                 f"({o['zmin']:.0f} mm) | {o['feet']:.0f} mm |")
    # PREDICTIONS rows
    rows = []
    rows.append(("E6", "Q09_2mm byte-identical to d2e Pa_2mm to 200 s", "see identity table",
                 "held" if ident else "REFUTED"))
    e1 = None
    if "Q10" in pairs:
        d = pairs["Q10"]
        late = d["bgap"][1:]
        e1 = abs(d["gap"]) <= 0.05 and all(abs(g) <= 0.05 for g in late)
        rows.append(("E1 (gate i)", "Q10 1 mm / 2 mm within 5 % over 150–300 s and every block after 50 s",
                     f"window {d['gap'] * 100:+.1f} %; blocks " + ", ".join(f"{g * 100:+.1f}" for g in d["bgap"]) + " %",
                     "held" if e1 else "REFUTED"))
    else:
        rows.append(("E1 (gate i)", "Q10 converges", "Q10_1mm not run", "not tested"))
    e2 = None
    if len(Qr) == 3:
        rops = np.array([Qr[p]["rop"] for p in QS])
        mean = rops.mean()
        spread = float(np.max(np.abs(rops / mean - 1.0)))
        e2 = spread <= 0.05
        rows.append(("E2 (gate ii)", "2 mm ROP for q 0.8 / 0.9 / 1.0 within ±5 % of the mean",
                     f"{Qr['Q08']['rop']:.3f} / {Qr['Q09']['rop']:.3f} / {Qr['Q10']['rop']:.3f} m/h, mean "
                     f"{mean:.3f}, spread ±{spread * 100:.1f} %; {(Qr['Q10']['rop'] - Qr['Q08']['rop']) / 2:+.3f} m/h "
                     f"per +0.1 q", "held" if e2 else "REFUTED"))
        e3 = Qr["Q08"]["rop"] > Qr["Q09"]["rop"] > Qr["Q10"]["rop"]
        rows.append(("E3", "ROP falls as q rises", f"{Qr['Q08']['rop']:.3f} > {Qr['Q09']['rop']:.3f} > "
                     f"{Qr['Q10']['rop']:.3f}?", "held" if e3 else "REFUTED"))
        for p in QS:
            q = QS[p]
            h = H[q]
            ok4 = abs(Qr[p]["rop"] / h["rop"] - 1.0) <= 0.27
            ok5 = abs(Qr[p]["sc"] * 1e3 - h["sc"]) <= 20.0
            rows.append(("E4", f"q {q:.1f} ROP within ±27 % of hand {h['rop']:.3f}",
                         f"{Qr[p]['rop']:.3f} m/h ({(Qr[p]['rop'] / h['rop'] - 1) * 100:+.1f} %)",
                         "held" if ok4 else "REFUTED"))
            rows.append(("E5", f"q {q:.1f} s_c at end within ±20 mm of hand {h['sc']:.1f}",
                         f"{Qr[p]['sc'] * 1e3:.1f} mm", "held" if ok5 else "REFUTED"))
    L += ["", "## PREDICTIONS rows", "", "| id | expectation | measurement | verdict |", "|---|---|---|---|"]
    L += [f"| {a_} | {b} | {c} | {v} |" for a_, b, c, v in rows]
    if e2 is None:
        oc = "not decidable (2 mm legs missing)"
    elif not e2:
        oc = "(c): q controls the rate"
    elif e1 is None:
        oc = "(ii) holds; (i) not tested"
    else:
        oc = "(a): score q = 1.0" if e1 else "(b): score q = 0.9 (jolt allowance)"
    L += ["", f"**Decision outcome: {oc}**"]
    txt = "\n".join(L)
    open(a.out, "w").write(txt + "\n")
    print(txt)
    if not a.no_figs:
        figures(pairs, Qr, H)


def figures(pairs, Qr, H):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if "Q10" in pairs:
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ref = {}
        for p in ("B0", "Pa"):
            if all(os.path.exists(os.path.join(D2E_OUT, f"{p}_{m}.done")) for m in ("2mm", "1mm")):
                b2 = E.blocks(os.path.join(D2E_OUT, f"{p}_2mm"), EDGES)[0]
                b1 = E.blocks(os.path.join(D2E_OUT, f"{p}_1mm"), EDGES)[0]
                ref[p] = [x1 / x2 - 1.0 for x1, x2 in zip(b1, b2) if np.isfinite(x1) and np.isfinite(x2)]
        series = dict(ref)
        series["Q10"] = pairs["Q10"]["bgap"]
        if "Q09" in pairs:
            series["Q09"] = pairs["Q09"]["bgap"]
        lab = {"B0": "D2e B0 (q 0.9 + clip)", "Pa": "D2e P-a (q 0.9, no clip)", "Q10": "Q10 (q 1.0, no clip)",
               "Q09": "Q09 (q 0.9, no clip)"}
        for k, g in series.items():
            mids = [(EDGES[i] + EDGES[i + 1]) / 2 for i in range(len(g))]
            ax.plot(mids, [x * 100 for x in g], "o-", label=lab[k])
        ax.axhline(0, color="k", lw=0.6)
        ax.axhspan(-5, 5, color="g", alpha=0.1, label="±5 %")
        ax.set(xlabel="block centre [s]", ylabel="1 mm / 2 mm burner ROP − 1 [%]", title="mesh gap per 50 s block")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(os.path.join(HERE, "d2f_blocks.png"), dpi=110)
        plt.close(fig)
        fig, ax = plt.subplots(figsize=(6.5, 4.5))
        xb = [(rim.BANDS[b] + rim.BANDS[b + 1]) / 2 * 1e3 for b in range(len(rim.BANDS) - 1)]
        for m, st in (("2mm", "-"), ("1mm", "--")):
            bl = pairs["Q10"][m]["rim"]["blocks"][1:]
            ax.plot(xb, np.mean([x["sp"] for x in bl], axis=0), "b" + st, label=f"{m} spall (mean of blocks after 50 s)")
        ax.axvspan(28, 40, color="k", alpha=0.06)
        ax.set(xlabel="r [mm]", ylabel="recession [mm/min]", title="Q10 rim by band")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(os.path.join(HERE, "d2f_rim.png"), dpi=110)
        plt.close(fig)
    if Qr:
        fig, ax = plt.subplots(figsize=(6.5, 4.5))
        qs = sorted(H)
        ax.plot(qs, [H[q]["rop"] for q in qs], "s--", color="gray", label="hand, transient 150–300 s")
        ax.plot(qs, [H[q]["steady"] for q in qs], "^:", color="gray", label="hand, steady")
        pts = sorted((QS[p], x["rop"]) for p, x in Qr.items())
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "o-", label="simulation 2 mm, 150–300 s")
        if "Q10" in pairs:
            ax.plot([1.0], [pairs["Q10"]["rop1"]], "D", label="simulation 1 mm, Q10")
        ax.set(xlabel="foot_pad_quantile q", ylabel="burner ROP [m/h]", title="ROP vs q (diagnostic, not scored)")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(HERE, "d2f_q.png"), dpi=110)
        plt.close(fig)


if __name__ == "__main__":
    main()
