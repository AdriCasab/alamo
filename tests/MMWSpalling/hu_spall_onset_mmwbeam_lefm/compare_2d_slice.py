#!/usr/bin/env python3
# Side-by-side 2D slice: damage_law vs LEFM (sp_weibull) at t=60s.
# xz mid-plane slice (y = 0.05 m). Both runs share the MMW beam axis at
# (x, y) = (0.05, 0.05), beam entering the +z face at z = 0.1 m.
#
# Layout (2 rows x 3 cols):
#   Row 1 (damage_law): Temperature [K], damage D, removed mask
#   Row 2 (LEFM):       Temperature [K], Sp_field,   removed mask
#
# Writes output/compare_2d_slice.png

import glob
import os
import re

import numpy as np
import yt

yt.funcs.mylog.setLevel(40)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
DL_DIR  = os.path.join(ROOT, "tests", "MMWSpalling", "hu_spall_onset_mmwbeam",
                       "output", "granite2")
LEF_DIR = os.path.join(ROOT, "tests", "MMWSpalling",
                       "hu_spall_onset_mmwbeam_lefm", "output", "granite2")
OUT_PNG = os.path.join(ROOT, "tests", "MMWSpalling",
                       "hu_spall_onset_mmwbeam_lefm", "output",
                       "compare_2d_slice.png")

NX = NY = NZ = 32
PLO = 0.0
PHI = 0.1
DX  = (PHI - PLO) / NX

# midplane y index (cells centred at y = (j+0.5)*dx; pick j closest to y=0.05)
J_MID = int(round((0.05 / DX) - 0.5))   # → 15

# y-bounds for the thin-slab cut (one cell thick around the midplane)
Y_LO = PLO + (J_MID    ) * DX
Y_HI = PLO + (J_MID + 1) * DX


def last_cell_plotfile(outdir):
    pat = re.compile(r"(\d+)cell$")
    cands = []
    for p in glob.glob(os.path.join(outdir, "*cell")):
        m = pat.search(os.path.basename(p))
        if m:
            cands.append((int(m.group(1)), p))
    return sorted(cands)[-1][1]


def to_numpy(a):
    return np.asarray(a.to_value() if hasattr(a, "to_value") else a)


def xz_slice(ds, field):
    """Extract an xz slice at y = J_MID midplane (single cell thick)."""
    ad = ds.all_data()
    x = to_numpy(ad["x"])
    y = to_numpy(ad["y"])
    z = to_numpy(ad["z"])
    f = to_numpy(ad[field])
    mask = (y > Y_LO) & (y < Y_HI)
    xs = x[mask]; zs = z[mask]; fs = f[mask]
    # Bin to NX x NZ grid
    ix = np.clip(((xs - PLO) / DX).astype(int), 0, NX - 1)
    iz = np.clip(((zs - PLO) / DX).astype(int), 0, NZ - 1)
    grid = np.full((NZ, NX), np.nan)
    grid[iz, ix] = fs
    return grid  # row 0 = z=0 (bottom), row NZ-1 = z=0.1 (top)


print(f"midplane y index J_MID = {J_MID} (y ∈ [{Y_LO:.4f}, {Y_HI:.4f}] m)")

# Load the final-time plotfile in each case
dl_path  = last_cell_plotfile(DL_DIR)
lef_path = last_cell_plotfile(LEF_DIR)
ds_dl  = yt.load(dl_path)
ds_lef = yt.load(lef_path)
print(f"damage_law: {os.path.basename(dl_path)}  t = {float(ds_dl.current_time):.2f} s")
print(f"LEFM:       {os.path.basename(lef_path)} t = {float(ds_lef.current_time):.2f} s")

# Extract slices
Temp_dl   = xz_slice(ds_dl,  "Temp")
D_dl      = xz_slice(ds_dl,  "D")
rem_dl    = xz_slice(ds_dl,  "removed")

Temp_lef  = xz_slice(ds_lef, "Temp")
Sp_lef    = xz_slice(ds_lef, "Sp_field")
hs_lef    = xz_slice(ds_lef, "h_spall_field")
rem_lef   = xz_slice(ds_lef, "removed")

# Convert removed → binary mask (1 = removed, 0 = solid)
rem_dl_mask  = (rem_dl  > 0.5).astype(float)
rem_lef_mask = (rem_lef > 0.5).astype(float)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize

fig, ax = plt.subplots(2, 3, figsize=(13.5, 8.5))

EXT = (PLO * 1e3, PHI * 1e3, PLO * 1e3, PHI * 1e3)  # mm
T_lo, T_hi = 280.0, max(np.nanmax(Temp_dl), np.nanmax(Temp_lef))

def add_beam_marker(a):
    a.axvline(50.0, color='w', ls=':', lw=0.8, alpha=0.6)
    a.text(50.5, 95, 'beam axis', color='w', fontsize=7, alpha=0.8)

# Row 1: damage_law
im00 = ax[0, 0].imshow(Temp_dl, origin='lower', extent=EXT, cmap='inferno',
                      vmin=T_lo, vmax=T_hi, aspect='equal')
ax[0, 0].set_title('damage_law: Temp [K]')
ax[0, 0].set_ylabel('z [mm] (surface at top)')
plt.colorbar(im00, ax=ax[0, 0], fraction=0.046, pad=0.04)
add_beam_marker(ax[0, 0])

im01 = ax[0, 1].imshow(D_dl, origin='lower', extent=EXT, cmap='magma',
                      vmin=0, vmax=1, aspect='equal')
ax[0, 1].set_title('damage_law: D')
plt.colorbar(im01, ax=ax[0, 1], fraction=0.046, pad=0.04)
add_beam_marker(ax[0, 1])

im02 = ax[0, 2].imshow(rem_dl_mask, origin='lower', extent=EXT, cmap='Greys_r',
                      vmin=0, vmax=1, aspect='equal')
ax[0, 2].set_title(f'damage_law: removed (black=solid, white=removed)\n'
                   f'cum cells removed: {int(np.sum(to_numpy(ds_dl.all_data()["removed"])>0.5))}')
plt.colorbar(im02, ax=ax[0, 2], fraction=0.046, pad=0.04)
add_beam_marker(ax[0, 2])

# Row 2: LEFM
im10 = ax[1, 0].imshow(Temp_lef, origin='lower', extent=EXT, cmap='inferno',
                      vmin=T_lo, vmax=T_hi, aspect='equal')
ax[1, 0].set_title('LEFM: Temp [K]')
ax[1, 0].set_xlabel('x [mm]'); ax[1, 0].set_ylabel('z [mm] (surface at top)')
plt.colorbar(im10, ax=ax[1, 0], fraction=0.046, pad=0.04)
add_beam_marker(ax[1, 0])

# Sp_field with overlay of h_spall_field contour
im11 = ax[1, 1].imshow(Sp_lef, origin='lower', extent=EXT, cmap='viridis',
                      vmin=0, vmax=max(1.0, np.nanmax(Sp_lef)), aspect='equal')
ax[1, 1].set_title('LEFM: Sp_field (= K_I / K_Ic)')
ax[1, 1].set_xlabel('x [mm]')
cb11 = plt.colorbar(im11, ax=ax[1, 1], fraction=0.046, pad=0.04)
cb11.ax.axhline(1.0, color='r', lw=1)  # firing threshold
add_beam_marker(ax[1, 1])

im12 = ax[1, 2].imshow(rem_lef_mask, origin='lower', extent=EXT, cmap='Greys_r',
                      vmin=0, vmax=1, aspect='equal')
ax[1, 2].set_title(f'LEFM: removed (black=solid, white=removed)\n'
                   f'cum cells removed: {int(np.sum(to_numpy(ds_lef.all_data()["removed"])>0.5))}')
ax[1, 2].set_xlabel('x [mm]')
plt.colorbar(im12, ax=ax[1, 2], fraction=0.046, pad=0.04)
add_beam_marker(ax[1, 2])

fig.suptitle(f'damage_law vs LEFM (sp_weibull) at t = {float(ds_dl.current_time):.1f} s — '
             f'xz mid-plane slice (y = {(Y_LO+Y_HI)/2*1e3:.1f} mm), MMW beam P0 = 10 kW',
             fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(OUT_PNG, dpi=120)
print(f"\nWrote: {OUT_PNG}")
