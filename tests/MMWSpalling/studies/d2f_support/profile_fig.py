#!/usr/bin/env python3
"""d2f_profile.png: azimuth-mean hole Ø(z) at the run end (300 s) for the
three 2 mm legs, D2e B0 (clip on) and Q10_1mm, against Meier Fig. 8.8's
visible width (a lower bound for the finished 500 mm hole). Direction only."""
import importlib.util
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sp = importlib.util.spec_from_file_location("d2f_an", os.path.join(HERE, "analyze.py"))
A = importlib.util.module_from_spec(sp)
sp.loader.exec_module(A)
sc = A.sc


def prof(pf):
    run = sc.load_run(pf)
    th = run.th
    t = float(th["time"][-1])
    feet = run.lz - float(np.interp(t, th["time"], th["foot_z"]))
    noz = run.lz - float(np.interp(t, th["time"], th["nozzle_z"]))
    zg = np.arange(0.0, feet, run.dz / 2.0)
    d_w = run.depth(t)
    return zg * 1e3, np.array([2.0 * sc.az_rwall(run, d_w, z) for z in zg]) * 1e3, noz * 1e3


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.4, 6.4))
    f88 = sc.load_fig()
    ax.plot([w for _, w, _ in f88], [z for z, _, _ in f88], "k-o", ms=3,
            label="Meier Fig. 8.8 (visible width, finished hole)")
    runs = [("Q10 2 mm (q 1.0)", "output/Q10_2mm", "-"), ("Q09 2 mm (q 0.9)", "output/Q09_2mm", "-"),
            ("Q08 2 mm (q 0.8)", "output/Q08_2mm", "-"), ("Q10 1 mm", "output/Q10_1mm", "--"),
            ("D2e B0 2 mm (q 0.9 + clip)", os.path.join(A.D2E_OUT, "B0_2mm"), ":")]
    noz = None
    for lab, pf, st in runs:
        pf = os.path.join(HERE, pf) if not os.path.isabs(pf) else pf
        z, D, noz_ = prof(pf)
        noz = noz_ if noz is None else noz
        ax.plot(D, z, st, label=lab)
    ax.axhline(noz, color="gray", lw=0.6, ls="--")
    ax.text(150, noz - 3, "nozzle plane at 300 s (walls above: final)", fontsize=7, color="gray")
    ax.axvline(80, color="r", lw=0.6, ls=":", label="burner Ø 80")
    ax.set(xlabel="hole Ø [mm]", ylabel="depth below the original surface [mm]", ylim=(140, 0),
           title="hole at 300 s (direction only, not scored)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "d2f_profile.png"), dpi=110)
    print("wrote d2f_profile.png")


if __name__ == "__main__":
    main()
