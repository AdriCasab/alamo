#!/usr/bin/env python3
# Compare the Rossi spall-event depth distribution from the elastic MLMG run
# (output/rossi_h_col_events.csv) vs the MLMG-free local-thermoelastic run
# (output/rossi_no_mlmg_h_col_events.csv). The h_col events CSV is the actual
# Rossi validation observable (Step 16c). If the two agree, the local
# thermoelastic stress reproduces the MLMG result without the solve.
#
#   python3 tests/MMWSpalling/validation/rossi/sp_rossi_damage_profile/compare_mlmg_vs_local.py

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(TEST_DIR, "output")
RUNS = [
    ("MLMG (el.type=static)",        os.path.join(OUT, "rossi_h_col_events.csv"),         "C0"),
    ("local TE (el.type=disable)",   os.path.join(OUT, "rossi_no_mlmg_h_col_events.csv"), "C3"),
]
# Rossi depth bins: 94 um wide (test uses ~dz/256-refined h_col)
BIN_W = 94.0   # um
BINS = np.arange(0.0, 1000.0 + BIN_W, BIN_W)


def load(path):
    if not os.path.isfile(path) or os.path.getsize(path) < 80:
        return None
    d = np.genfromtxt(path, delimiter=",", names=True)
    if d.size == 0:
        return None
    return d


def main():
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    print(f"{'run':28} {'events':>7} {'first_t':>8} {'peak_bin_um':>12} "
          f"{'h_col_max_um':>13} {'GB_frac':>8} {'Sp_max':>7}")
    for name, path, color in RUNS:
        d = load(path)
        if d is None:
            print(f"{name:28} -- no events / missing csv --")
            continue
        h_um = np.atleast_1d(d["h_col"]) * 1e6
        t = np.atleast_1d(d["time"])
        is_gb = np.atleast_1d(d["is_gb_at_top"])
        sp = np.atleast_1d(d["Sp_field_at_top"])
        counts, edges = np.histogram(h_um, bins=BINS)
        centres = 0.5 * (edges[:-1] + edges[1:])
        peak_idx = int(np.argmax(counts))
        peak_bin = centres[peak_idx]
        gb_frac = float(np.mean(is_gb))
        print(f"{name:28} {len(h_um):7d} {t.min():8.1f} {peak_bin:12.1f} "
              f"{h_um.max():13.1f} {gb_frac:8.3f} {sp.max():7.3f}")

        axes[0].step(centres, counts, where="mid", color=color, label=name, lw=2)
        axes[1].plot(t, h_um, ".", color=color, label=name, alpha=0.4, ms=4)

    ax = axes[0]
    ax.axvspan(100, 200, color="green", alpha=0.10, label="Rossi peak 100-200 um")
    ax.axvline(520, color="0.5", ls=":", lw=1, label="Rossi falloff ~520 um")
    ax.set_xlabel("spall-event depth h_col [um]")
    ax.set_ylabel("event count per 94 um bin")
    ax.set_title("Depth distribution (Rossi observable)")
    ax.set_xlim(0, 1000); ax.legend(fontsize=8)

    ax = axes[1]
    ax.set_xlabel("time [s]"); ax.set_ylabel("event depth h_col [um]")
    ax.set_title("Event depth vs time")
    ax.legend(fontsize=8)

    fig.tight_layout()
    out_png = os.path.join(OUT, "rossi_mlmg_vs_local.png")
    fig.savefig(out_png, dpi=120)
    print(f"\nWrote {out_png}")


if __name__ == "__main__":
    main()
