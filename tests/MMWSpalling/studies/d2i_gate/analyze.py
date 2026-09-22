#!/usr/bin/env python3
"""D2i analysis: the CRITERION.md verdict, per-block diagnostics, the silent
fraction table, the idle sensitivity trio and d2i_blocks.png. Writes
analysis.md.

Reuses d2g_mesh/analyze.py (its blocks_100 / blocks_50 / gaps / diag, which in
turn use d2f's q-aware rim wrapper, d2e analyze_pairs and rim.py). Nothing
frozen is edited.
"""
import contextlib
import importlib.util
import io
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)


def load(name, path):
    sp = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


I = load("d2i_run", os.path.join(HERE, "run.py"))
GA = load("d2g_analyze", os.path.join(STUD, "d2g_mesh", "analyze.py"))
A, rim = GA.A, GA.rim
OUT = I.OUT
D2G_OUT = os.path.join(STUD, "d2g_mesh", "output")
D2H_OUT = os.path.join(STUD, "d2h_openloop", "output")
TOL = 0.05
SNAPS = (250.0, 350.0, 450.0, 550.0)
SBANDS = np.arange(0.024, 0.0521, 0.004)         # 4 mm bands, 24-52 mm (D2h §2)
RHOCP = I.R.wj.RHOCP


def silent(pf, times):
    """D2h §2: per 4 mm band, the share of columns whose time since the last
    removal event exceeds 2 t_cell, t_cell = rho cp dz 528 / 3e5; plus the
    firings per pad-annulus column over the preceding 100 s."""
    m = rim.meta_of(pf)
    g = rim.Geometry(m)
    t_cell = RHOCP * g.dz * 528.0 / 3.0e5
    ev = rim.load_events(pf)
    band = np.digitize(g.r, SBANDS) - 1
    band[(g.r < SBANDS[0]) | (g.r >= SBANDS[-1])] = -1
    i, j = ev[:, 1].astype(int), ev[:, 2].astype(int)
    rows = []
    for t in times:
        last = np.zeros((g.n, g.n))
        s = ev[:, 0] <= t
        np.maximum.at(last, (i[s], j[s]), ev[s, 0])
        sil = (t - last) > 2.0 * t_cell
        fr = [100.0 * float(np.mean(sil[band == b])) if np.any(band == b) else float("nan")
              for b in range(len(SBANDS) - 1)]
        w = s & (ev[:, 0] > t - 100.0) & g.ann[i, j]
        rows.append((t, fr, int(w.sum()) / max(int(g.ann.sum()), 1)))
    return rows, t_cell


def th_block(th, a, b, col):
    s = (th["time"] >= a) & (th["time"] <= b)
    return float(np.mean(th[col][s]))


def main():
    L = []
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        ident = I.identity()
    L += ["## 1. Identity (rows t <= 349.9 s vs d2h I0)", "", buf.getvalue()]

    pf2, pf1 = (os.path.join(OUT, n) for n in ("LI0_2mm", "LI0_1mm"))
    th2, th1 = A.thermo(pf2), A.thermo(pf1)
    m2, m1 = rim.meta_of(pf2), rim.meta_of(pf1)
    t_m = min(float(th2["time"][-1]), float(th1["time"][-1]))
    L += ["## 2. Runs", "", "| run | status | t_end | wall | dt | plotfiles |", "|---|---|---|---|---|---|"]
    for n, m in (("LI0_2mm", m2), ("LI0_1mm", m1)):
        npf = len([d for d in os.listdir(os.path.join(OUT, n)) if d.endswith("cell")])
        L.append(f"| {n} | {m['status']} | {m['t_end']:.1f} s | {m['wall_s'] / 60:.1f} min | "
                 f"{m['dt'] * 1e3:g} ms | {npf} |")

    r2, r1 = A.window_rop(th2, 150.0, t_m), A.window_rop(th1, 150.0, t_m)
    gw = r1 / r2 - 1.0
    b100 = GA.gaps(th2, th1, GA.blocks_100(t_m))
    b50 = GA.gaps(th2, th1, GA.blocks_50(t_m))
    ok_a = abs(gw) <= TOL
    ok_b = all(abs(g) <= TOL for *_, g in b100)
    mids = np.array([(a + b) / 2 for a, b, *_ in b100])
    gv = np.array([g for *_, g in b100]) * 100
    slope = float(np.polyfit(mids, gv, 1)[0]) * 100.0 if len(b100) >= 2 else float("nan")
    L += ["", f"## 3. The pair (t_m = {t_m:.1f} s)", "",
          f"Matched window 150–{t_m:.0f} s: 2 mm {r2:.3f} m/h, 1 mm {r1:.3f} m/h, gap "
          f"**{gw * 100:+.1f} %** ({'within' if ok_a else 'OUTSIDE'} 5 %).", "",
          "| 100 s block | 2 mm ROP | 1 mm ROP | gap |", "|---|---|---|---|"]
    L += [f"| {a:.0f}–{b:.0f} s | {x2:.3f} | {x1:.3f} | {g * 100:+.1f} % |" for a, b, x2, x1, g in b100]
    L += ["", f"Growth diagnostic: {slope:+.2f} points per 100 s (100 s block gaps vs mid-time).", "",
          "| 50 s block | 2 mm ROP | 1 mm ROP | gap |", "|---|---|---|---|"]
    L += [f"| {a:.0f}–{b:.0f} s | {x2:.3f} | {x1:.3f} | {g * 100:+.1f} % |" for a, b, x2, x1, g in b50]

    ed = [(a, b) for a, b, *_ in b100]
    d2, w2 = GA.diag(pf2, th2, ed)
    d1, w1 = GA.diag(pf1, th1, ed)
    L += ["", "## 4. Per 100 s block (2 mm / 1 mm)", "",
          "| block | rim spall 38–40 [mm/min] | annulus mean [mm/min] | 40–60 mm removal [cm³/s] | "
          "pinned share | idle cols | T_rec [K] | s_c [mm] | reach [mm] | P_face [W] | max abs ledger_err |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]
    for x, y, (a, b) in zip(d2, d1, ed):
        L.append(f"| {a:.0f}–{b:.0f} s | {x['rim_sp']:.1f} / {y['rim_sp']:.1f} | {x['ann']:.1f} / "
                 f"{y['ann']:.1f} | {x['flank']:.2f} / {y['flank']:.2f} | {x['pin']:.3f} / {y['pin']:.3f} | "
                 f"{th_block(th2, a, b, 'pinned_cols_face_idle'):.0f} / "
                 f"{th_block(th1, a, b, 'pinned_cols_face_idle'):.0f} | "
                 f"{x['Trec']:.0f} / {y['Trec']:.0f} | {x['sc']:.1f} / {y['sc']:.1f} | "
                 f"{x['reach']:.1f} / {y['reach']:.1f} | "
                 f"{th_block(th2, a, b, 'jet_P_face'):.0f} / {th_block(th1, a, b, 'jet_P_face'):.0f} | "
                 f"{x['led']:.1e} / {y['led']:.1e} |")
    L += ["", f"rim.py foot_z reconstruction worst |Δ|: 2 mm {w2 * 1e3:.3f} mm, 1 mm {w1 * 1e3:.3f} mm.", ""]

    times = [t for t in SNAPS if t <= t_m]
    L += ["## 5. Silent fraction [%] by 4 mm band, and pad firings per column per 100 s", "",
          "| run | t | " + " | ".join(f"{a * 1e3:.0f}–{b * 1e3:.0f}" for a, b in zip(SBANDS[:-1], SBANDS[1:]))
          + " | fir/col |", "|---|---|" + "---|" * (len(SBANDS) - 1) + "---|"]
    for n, pf in (("LI0 2 mm", pf2), ("LI0 1 mm", pf1)):
        rows, tc = silent(pf, times)
        for t, fr, fc in rows:
            L.append(f"| {n} | {t:.0f} s | " + " | ".join(f"{v:.0f}" for v in fr) + f" | {fc:.1f} |")
        L.append(f"| | t_cell = {tc:.2f} s | " + " | " * (len(SBANDS) - 1) + " |")
    L += ["", "Silence in the 44–52 mm bands is the `s <= 0` rule above the nozzle plane and does "
          "not count against the gate (CRITERION.md).", ""]

    L += ["## 6. Idle sensitivity at 1 mm: burner ROP 250–350 s, idle columns, silence at 350 s", "",
          "| pinned_idle_cycles | run | ROP 250–350 s | vs idle off | idle cols at 350 s | "
          + " | ".join(f"{a * 1e3:.0f}–{b * 1e3:.0f}" for a, b in zip(SBANDS[:-1], SBANDS[1:])) + " |",
          "|---|---|---|---|---|" + "---|" * (len(SBANDS) - 1)]
    trio = [("2", os.path.join(D2G_OUT, "L09_1mm")), ("20", os.path.join(OUT, "I20_1mm")),
            ("1e9 (off)", os.path.join(D2H_OUT, "I0_1mm"))]
    vals = [A.window_rop(A.thermo(pf), 250.0, 350.0) for _, pf in trio]
    for (lab, pf), v in zip(trio, vals):
        th = A.thermo(pf)
        idl = th_block(th, 349.0, 350.0, "pinned_cols_face_idle")
        fr = silent(pf, [350.0])[0][0][1]
        L.append(f"| {lab} | {os.path.basename(pf)} | {v:.3f} m/h | {(v / vals[-1] - 1) * 100:+.1f} % | "
                 f"{idl:.0f} | " + " | ".join(f"{x:.0f}" for x in fr) + " |")

    verdict = "PASS" if (ok_a and ok_b) else "FAIL"
    L += ["", f"## 7. Verdict: {verdict}", "",
          f"(a) matched window {gw * 100:+.1f} % {'ok' if ok_a else 'fails'}; (b) 100 s blocks "
          + ", ".join(f"{g * 100:+.1f}" for *_, g in b100) + f" % {'all ok' if ok_b else 'fail'}. "
          f"Identity: {'PASS' if ident else 'FAIL'}."]
    txt = "\n".join(L)
    open(os.path.join(HERE, "analysis.md"), "w").write(txt + "\n")
    print(txt)
    fig(b50, b100, t_m)


def fig(b50, b100, t_m):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    f, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.plot([(a + b) / 2 for a, b, *_ in b50], [g * 100 for *_, g in b50], "o-", ms=4,
            color="C0", label="LI0 (idle off), 50 s blocks")
    ax.plot([(a + b) / 2 for a, b, *_ in b100], [g * 100 for *_, g in b100], "s-", ms=7, lw=2,
            color="C0", label="LI0, 100 s blocks (criterion)")
    for lab, d, t, c in (("D2g L09 (idle on), 50 s", D2G_OUT, "L09", "gray"),
                         ("D2h I0 (idle off, 350 s), 50 s", D2H_OUT, "I0", "C2")):
        try:
            a2, a1 = A.thermo(os.path.join(d, f"{t}_2mm")), A.thermo(os.path.join(d, f"{t}_1mm"))
            tt = min(float(a2["time"][-1]), float(a1["time"][-1]))
            q = GA.gaps(a2, a1, GA.blocks_50(tt))
            ax.plot([(x + y) / 2 for x, y, *_ in q], [g * 100 for *_, g in q], "x--", color=c, label=lab)
        except OSError:
            pass
    ax.axhline(0, color="k", lw=0.6)
    ax.axhspan(-5, 5, color="g", alpha=0.1, label="±5 %")
    ax.set(xlabel="block centre [s]", ylabel="1 mm / 2 mm burner ROP − 1 [%]",
           title="D2i: support rule with the idle clock off, to the 0.40 m bottom stop")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    f.tight_layout()
    f.savefig(os.path.join(HERE, "d2i_blocks.png"), dpi=110)
    plt.close(f)
    print("wrote d2i_blocks.png")


if __name__ == "__main__":
    main()
