#!/usr/bin/env python3
# Side-by-side comparison: hu_spall_onset_mmwbeam (damage_law) vs
# hu_spall_onset_mmwbeam_lefm (sp_weibull LEFM). Identical MMW beam,
# Voronoi granite, BCs, losses. Different spallation criterion.
#
# Writes output/compare_damage_vs_lefm.png with 6 panels:
#   1. T_max(t)              — heating sanity (should be similar)
#   2. removed cells(t)      — drilling progression
#   3. drill_depth(t) [mm]   — surface advance
#   4. damage_law D_max(t)  +  LEFM Sp_field max(t) — criterion drive
#   5. h_spall_field max(t) [mm]                   — LEFM only
#   6. spall_event count(t) per snapshot           — per-step firing

import glob
import os
import re
import sys

import numpy as np
import yt

yt.funcs.mylog.setLevel(40)

ROOT = os.path.dirname(os.path.abspath(__file__)).split(os.sep + os.path.join("tests", "MMWSpalling") + os.sep)[0]
DAMAGE_DIR = os.path.join(ROOT, "tests", "MMWSpalling", "extensions", "hu_spall_onset_mmwbeam",
                          "output", "granite2")
LEFM_DIR   = os.path.join(ROOT, "tests", "MMWSpalling",
                          "hu_spall_onset_mmwbeam_lefm", "output", "granite2")
OUT_PNG    = os.path.join(ROOT, "tests", "MMWSpalling",
                          "hu_spall_onset_mmwbeam_lefm", "output",
                          "compare_damage_vs_lefm.png")

Z_TOP = 0.1   # m
PHI_HI = 0.1  # m, upper z bound for drill-depth computation


def plotfiles_cell(outdir):
    pat = re.compile(r"(\d+)cell$")
    out = []
    for p in glob.glob(os.path.join(outdir, "*cell")):
        m = pat.search(os.path.basename(p))
        if m:
            out.append((int(m.group(1)), p))
    out.sort()
    return [p for _, p in out]


def to_numpy(a):
    return np.asarray(a.to_value() if hasattr(a, "to_value") else a)


def summarize(ds, has_lefm=False):
    ad = ds.all_data()
    fields = {f[1] for f in ds.field_list}
    T = to_numpy(ad["Temp"])
    rem = to_numpy(ad["removed"]) > 0.5
    n_rem = int(np.count_nonzero(rem))
    if n_rem > 0:
        z = to_numpy(ad["z"])
        drill = float(PHI_HI - np.min(z[rem]))
    else:
        drill = 0.0
    out = {
        "t":          float(ds.current_time),
        "T_max":      float(np.max(T)),
        "n_removed":  n_rem,
        "drill":      drill,
        "D_max":      float(np.max(to_numpy(ad["D"]))) if "D" in fields else 0.0,
        "spall_ev":   int(np.count_nonzero(to_numpy(ad["spall_event"]) > 0.5))
                      if "spall_event" in fields else 0,
        "vapor_ev":   int(np.count_nonzero(to_numpy(ad["vapor_event"]) > 0.5))
                      if "vapor_event" in fields else 0,
    }
    if has_lefm:
        out["Sp_max"] = float(np.max(to_numpy(ad["Sp_field"]))) if "Sp_field" in fields else 0.0
        out["hs_max"] = float(np.max(to_numpy(ad["h_spall_field"]))) if "h_spall_field" in fields else 0.0
    return out


print("Loading damage_law plotfiles...")
dl_files = plotfiles_cell(DAMAGE_DIR)
print(f"  {len(dl_files)} plotfiles in {DAMAGE_DIR}")

print("Loading LEFM plotfiles...")
lef_files = plotfiles_cell(LEFM_DIR)
print(f"  {len(lef_files)} plotfiles in {LEFM_DIR}")

dl_series  = [summarize(yt.load(f), has_lefm=False) for f in dl_files]
lef_series = [summarize(yt.load(f), has_lefm=True)  for f in lef_files]

# --- summary table ---
print()
print(f"{'metric':<28} {'damage_law':>14} {'LEFM':>14}")
print("-" * 60)
print(f"{'final removed cells':<28} {dl_series[-1]['n_removed']:>14d} {lef_series[-1]['n_removed']:>14d}")
print(f"{'final drill depth [mm]':<28} {dl_series[-1]['drill']*1e3:>14.3f} {lef_series[-1]['drill']*1e3:>14.3f}")
print(f"{'peak T_max [K]':<28} "
      f"{max(s['T_max'] for s in dl_series):>14.1f} "
      f"{max(s['T_max'] for s in lef_series):>14.1f}")
print(f"{'final T_max [K]':<28} "
      f"{dl_series[-1]['T_max']:>14.1f} "
      f"{lef_series[-1]['T_max']:>14.1f}")
print(f"{'peak D_max':<28} "
      f"{max(s['D_max'] for s in dl_series):>14.4f} "
      f"{max(s['D_max'] for s in lef_series):>14.4f}")
print(f"{'peak Sp_field max':<28} {'N/A':>14} "
      f"{max(s.get('Sp_max', 0.0) for s in lef_series):>14.4f}")
print(f"{'peak h_spall_field [mm]':<28} {'N/A':>14} "
      f"{max(s.get('hs_max', 0.0) for s in lef_series)*1e3:>14.3f}")
spall_total_dl  = sum(s['spall_ev'] for s in dl_series)
spall_total_lef = sum(s['spall_ev'] for s in lef_series)
print(f"{'cum spall_event snapshots':<28} {spall_total_dl:>14d} {spall_total_lef:>14d}")

# Onset time: first plotfile with >= 1 removed cell
def onset(series):
    for s in series:
        if s['n_removed'] > 0:
            return s['t']
    return None

print(f"{'onset (first removal) [s]':<28} "
      f"{onset(dl_series) if onset(dl_series) else 'N/A':>14} "
      f"{onset(lef_series) if onset(lef_series) else 'N/A':>14}")

# Mean rate of penetration (final / final-t)
final_t = dl_series[-1]['t']
print(f"{'mean RoP [mm/s]':<28} "
      f"{dl_series[-1]['drill']*1e3/final_t:>14.4f} "
      f"{lef_series[-1]['drill']*1e3/final_t:>14.4f}")

# --- plot ---
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 3, figsize=(15, 8.5))

    t_dl  = [s['t'] for s in dl_series]
    t_lef = [s['t'] for s in lef_series]

    # 1: T_max
    ax[0, 0].plot(t_dl,  [s['T_max'] for s in dl_series],  'C3o-', label='damage_law', lw=1.6, ms=3)
    ax[0, 0].plot(t_lef, [s['T_max'] for s in lef_series], 'C0s--', label='LEFM',       lw=1.6, ms=3)
    ax[0, 0].axhline(3233.0, color='gray', ls=':', lw=1, alpha=0.7)
    ax[0, 0].annotate('T_vap_lo (3233 K)', (1.0, 3270), color='gray', fontsize=8)
    ax[0, 0].set_xlabel('t [s]'); ax[0, 0].set_ylabel('max Temp [K]')
    ax[0, 0].set_title('Peak temperature (MMW beam, identical)')
    ax[0, 0].legend(loc='lower right'); ax[0, 0].grid(alpha=0.3)

    # 2: removed cells (cumulative)
    ax[0, 1].plot(t_dl,  [s['n_removed'] for s in dl_series],  'C3o-', label='damage_law', lw=1.6, ms=3)
    ax[0, 1].plot(t_lef, [s['n_removed'] for s in lef_series], 'C0s--', label='LEFM',       lw=1.6, ms=3)
    ax[0, 1].set_xlabel('t [s]'); ax[0, 1].set_ylabel('cumulative removed cells')
    ax[0, 1].set_title('Drilling progression')
    ax[0, 1].legend(loc='lower right'); ax[0, 1].grid(alpha=0.3)

    # 3: drill depth
    ax[0, 2].plot(t_dl,  [s['drill']*1e3 for s in dl_series],  'C3o-', label='damage_law', lw=1.6, ms=3)
    ax[0, 2].plot(t_lef, [s['drill']*1e3 for s in lef_series], 'C0s--', label='LEFM',       lw=1.6, ms=3)
    ax[0, 2].set_xlabel('t [s]'); ax[0, 2].set_ylabel('drill depth [mm]')
    ax[0, 2].set_title('Surface advance')
    ax[0, 2].legend(loc='lower right'); ax[0, 2].grid(alpha=0.3)

    # 4: criterion drive (D vs Sp)
    ax[1, 0].plot(t_dl,  [s['D_max'] for s in dl_series],  'C3o-', label='damage_law D_max', lw=1.6, ms=3)
    ax[1, 0].set_xlabel('t [s]'); ax[1, 0].set_ylabel('damage D_max', color='C3')
    ax[1, 0].tick_params(axis='y', labelcolor='C3')
    ax[1, 0].axhline(0.51, color='C3', ls=':', lw=1, alpha=0.6)
    ax[1, 0].grid(alpha=0.3)
    ax4_r = ax[1, 0].twinx()
    ax4_r.plot(t_lef, [s.get('Sp_max', 0.0) for s in lef_series], 'C0s--', label='LEFM Sp_max', lw=1.6, ms=3)
    ax4_r.set_ylabel('LEFM Sp_field max', color='C0')
    ax4_r.tick_params(axis='y', labelcolor='C0')
    ax4_r.axhline(1.0, color='C0', ls=':', lw=1, alpha=0.6)
    ax[1, 0].set_title('Criterion drive (D vs Sp)')

    # 5: h_spall_field for LEFM
    ax[1, 1].plot(t_lef, [s.get('hs_max', 0.0)*1e3 for s in lef_series], 'C0s--', lw=1.6, ms=3)
    ax[1, 1].set_xlabel('t [s]'); ax[1, 1].set_ylabel('h_spall_field max [mm]')
    ax[1, 1].set_title('LEFM running-max crack depth (per column)')
    ax[1, 1].grid(alpha=0.3)

    # 6: spall_event per-step firing count
    ax[1, 2].plot(t_dl,  [s['spall_ev'] for s in dl_series],  'C3o-', label='damage_law', lw=1.6, ms=3)
    ax[1, 2].plot(t_lef, [s['spall_ev'] for s in lef_series], 'C0s--', label='LEFM',       lw=1.6, ms=3)
    ax[1, 2].set_xlabel('t [s]'); ax[1, 2].set_ylabel('spall_event count per snapshot')
    ax[1, 2].set_title('Per-snapshot spall firings')
    ax[1, 2].legend(loc='upper right'); ax[1, 2].grid(alpha=0.3)

    fig.suptitle('Side-by-side: damage_law vs LEFM (sp_weibull) under identical MMW beam '
                 '(P0=10 kW, ω0=2 cm)', fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(OUT_PNG, dpi=120)
    print(f"\nWrote: {OUT_PNG}")
except ImportError:
    print("(matplotlib not installed - skipping PNG plot)")
