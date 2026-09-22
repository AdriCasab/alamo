#!/usr/bin/env python3
"""D2g Stage 1 analysis: the CRITERION.md verdict, per-block diagnostics and
d2g1_blocks.png. Writes analysis.md.

Reuses d2f analyze.py (loaded by path; its rim.QUANT wrapper, d2e
analyze_pairs thermo / window_rop and rim.py), unchanged.
"""
import contextlib
import importlib.util
import io
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)


def load(name, path):
    sp = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


G = load("d2g_run", os.path.join(HERE, "run.py"))
FA = load("d2f_analyze", os.path.join(STUD, "d2f_support", "analyze.py"))
A, rim = FA.A, FA.rim            # rim.analyse is FA's q-aware wrapper
OUT = G.OUT
D2F_OUT = os.path.join(STUD, "d2f_support", "output")
TOL = 0.05


def blocks_100(t_m):
    ed, a = [], 150.0
    while a + 100.0 <= t_m + 1e-9:
        ed.append((a, a + 100.0))
        a += 100.0
    if t_m - a >= 50.0:
        ed.append((a, t_m))
    return ed


def blocks_50(t_m):
    ed, a = [], 0.0
    while a + 50.0 <= t_m + 1e-9:
        ed.append((a, a + 50.0))
        a += 50.0
    if t_m - a >= 25.0:
        ed.append((a, t_m))
    return ed


def gaps(th2, th1, ed):
    out = []
    for a, b in ed:
        r2, r1 = A.window_rop(th2, a, b), A.window_rop(th1, a, b)
        out.append((a, b, r2, r1, r1 / r2 - 1.0))
    return out


def diag(pf, th, ed):
    """Per-block rim/flank (rim.py) and thermo means."""
    r = rim.analyse(pf, [ed[0][0]] + [b for _, b in ed])
    ri = int(np.nonzero(np.isclose(rim.BANDS[:-1], 0.038))[0][0])
    rows = []
    for (a, b), bl in zip(ed, r["blocks"]):
        s = (th["time"] >= a) & (th["time"] <= b)
        pin = th["pinned_cols_pinned"][s] / np.maximum(th["pinned_cols_pinned"][s] + th["pinned_cols_face"][s], 1e-30)
        rows.append(dict(a=a, b=b, rim_sp=bl["sp"][ri], ann=bl["ann_mean"], flank=bl["flank_cm3s"],
                         Trec=float(np.mean(th["jet_T_rec"][s])), sc=float(np.mean(th["jet_s_c"][s])) * 1e3,
                         reach=float(np.mean(th["jet_r_reach"][s])) * 1e3, pin=float(np.mean(pin)),
                         led=float(np.max(np.abs(th["ledger_err"][s])))))
    return rows, r["worst"]


def main():
    L = []
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        ident = G.identity()
    L += ["## Identity (rows t <= 299.9 s vs d2f Q09)", "", buf.getvalue()]
    pf2, pf1 = (os.path.join(OUT, n) for n in ("L09_2mm", "L09_1mm"))
    th2, th1 = A.thermo(pf2), A.thermo(pf1)
    m2, m1 = rim.meta_of(pf2), rim.meta_of(pf1)
    t_m = min(float(th2["time"][-1]), float(th1["time"][-1]))
    L += ["## Runs", "", "| run | status | t_end | wall |", "|---|---|---|---|"]
    for n, m in (("L09_2mm", m2), ("L09_1mm", m1)):
        L.append(f"| {n} | {m['status']} | {m['t_end']:.1f} s | {m['wall_s'] / 60:.1f} min |")
    r2, r1 = A.window_rop(th2, 150.0, t_m), A.window_rop(th1, 150.0, t_m)
    gw = r1 / r2 - 1.0
    b100 = gaps(th2, th1, blocks_100(t_m))
    b50 = gaps(th2, th1, blocks_50(t_m))
    ok_a = abs(gw) <= TOL
    ok_b = all(abs(g) <= TOL for *_, g in b100)
    mids = np.array([(a + b) / 2 for a, b, *_ in b100])
    gv = np.array([g for *_, g in b100]) * 100
    slope = float(np.polyfit(mids, gv, 1)[0]) * 100.0 if len(b100) >= 2 else float("nan")
    L += ["", f"## Pair (t_m = {t_m:.1f} s)", "",
          f"Matched window 150–{t_m:.0f} s: 2 mm {r2:.3f} m/h, 1 mm {r1:.3f} m/h, gap **{gw * 100:+.1f} %** "
          f"({'within' if ok_a else 'OUTSIDE'} ±5 %).", "",
          "| 100 s block | 2 mm ROP | 1 mm ROP | gap |", "|---|---|---|---|"]
    L += [f"| {a:.0f}–{b:.0f} s | {x2:.3f} | {x1:.3f} | {g * 100:+.1f} % |" for a, b, x2, x1, g in b100]
    L += ["", f"Growth diagnostic: {slope:+.2f} points per 100 s (100 s block gaps vs mid-time).", "",
          "| 50 s block | 2 mm ROP | 1 mm ROP | gap |", "|---|---|---|---|"]
    L += [f"| {a:.0f}–{b:.0f} s | {x2:.3f} | {x1:.3f} | {g * 100:+.1f} % |" for a, b, x2, x1, g in b50]
    ed = [(a, b) for a, b, *_ in b100]
    d2, w2 = diag(pf2, th2, ed)
    d1, w1 = diag(pf1, th1, ed)
    L += ["", "## Per 100 s block diagnostics (2 mm / 1 mm)", "",
          "| block | rim spall 38–40 [mm/min] | annulus mean [mm/min] | 40–60 mm removal [cm³/s] | mean T_rec [K] | "
          "mean s_c [mm] | mean reach [mm] | pinned share | max abs ledger_err |",
          "|---|---|---|---|---|---|---|---|---|"]
    for x, y in zip(d2, d1):
        L.append(f"| {x['a']:.0f}–{x['b']:.0f} s | {x['rim_sp']:.1f} / {y['rim_sp']:.1f} | {x['ann']:.1f} / "
                 f"{y['ann']:.1f} | {x['flank']:.2f} / {y['flank']:.2f} | {x['Trec']:.0f} / {y['Trec']:.0f} | "
                 f"{x['sc']:.1f} / {y['sc']:.1f} | {x['reach']:.1f} / {y['reach']:.1f} | {x['pin']:.2f} / "
                 f"{y['pin']:.2f} | {x['led']:.1e} / {y['led']:.1e} |")
    L += ["", f"rim.py foot_z reconstruction worst |Δ|: 2 mm {w2 * 1e3:.3f} mm, 1 mm {w1 * 1e3:.3f} mm.", ""]
    verdict = "PASS" if (ok_a and ok_b) else "FAIL"
    L += [f"**Verdict: {verdict}** — (a) window {gw * 100:+.1f} % {'ok' if ok_a else 'fails'}; (b) 100 s blocks "
          + ", ".join(f"{g * 100:+.1f}" for *_, g in b100) + f" % {'all ok' if ok_b else 'fail'}. "
          f"Identity: {'PASS' if ident else 'FAIL'}."]
    txt = "\n".join(L)
    open(os.path.join(HERE, "analysis.md"), "w").write(txt + "\n")
    print(txt)
    fig(b50, b100)


def fig(b50, b100):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    f, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.plot([(a + b) / 2 for a, b, *_ in b50], [g * 100 for *_, g in b50], "o-", ms=4, label="L09, 50 s blocks")
    ax.plot([(a + b) / 2 for a, b, *_ in b100], [g * 100 for *_, g in b100], "s-", ms=7, lw=2,
            label="L09, 100 s blocks (criterion)")
    q2, q1 = A.thermo(os.path.join(D2F_OUT, "Q09_2mm")), A.thermo(os.path.join(D2F_OUT, "Q09_1mm"))
    q = gaps(q2, q1, blocks_50(299.9))
    ax.plot([(a + b) / 2 for a, b, *_ in q], [g * 100 for *_, g in q], "x--", color="gray",
            label="D2f Q09, 50 s blocks (to 300 s)")
    ax.axhline(0, color="k", lw=0.6)
    ax.axhspan(-5, 5, color="g", alpha=0.1, label="±5 %")
    ax.set(xlabel="block centre [s]", ylabel="1 mm / 2 mm burner ROP − 1 [%]",
           title="D2g-1: q 0.9 mesh gap to the 0.40 m bottom stop")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    f.tight_layout()
    f.savefig(os.path.join(HERE, "d2g1_blocks.png"), dpi=110)
    plt.close(f)
    print("wrote d2g1_blocks.png")


if __name__ == "__main__":
    main()
