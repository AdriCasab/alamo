#!/usr/bin/env python3
"""Vertical section through the hole axis, per plotfile, for S0 / W1 / B1L, with Meier's
Fig. 8.8 outline and the burner drawn in; frames -> GIF. Reads only `Temp` and `removed`."""
import os, re, sys, math, struct
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from PIL import Image

ROOT = "/Users/tzetze20/amr_tools/alamo"
ST = ROOT + "/tests/MMWSpalling/studies"
LEGS = [("S0  (floor x0.743, wall OFF)", ST + "/p2_4_body/legs/S0"),
        ("W1  (wall ON, march)", ST + "/p2_3_wall/legs/W1"),
        ("B1L (wall ON, march + cold body)", ST + "/p2_4_body/legs/B1L")]
NX, NY, NZ, DX = 60, 60, 300, 0.002
Z_TOP = 0.5
DT = 0.016
MEIER = [(10, 96), (25, 92), (50, 88), (75, 87), (85, 87), (100, 87), (125, 87), (150, 86),
         (175, 85), (200, 85), (250, 84), (300, 84), (350, 83), (400, 82), (450, 80), (480, 78)]
HDR = re.compile(rb"FAB \(\((\d+), \(\[?[^)]*\)\),\(\d+, \(\d+ \d+ \d+ \d+ \d+ \d+ \d+ \d+\)\)\)\(\((-?\d+),(-?\d+),(-?\d+)\) \((-?\d+),(-?\d+),(-?\d+)\) \((\d+),(\d+),(\d+)\)\) (\d+)\n")

def var_index(pf, name):
    L = open(os.path.join(pf, "Header")).read().split("\n")
    n = int(L[1]); names = L[2:2 + n]
    return names.index(name), n

def load(pf, comps):
    """-> dict name -> np.array (NX,NY,NZ)"""
    idx = {c: var_index(pf, c)[0] for c in comps}
    ncomp = var_index(pf, comps[0])[1]
    out = {c: np.zeros((NX, NY, NZ)) for c in comps}
    lvl = os.path.join(pf, "Level_0")
    for fn in sorted(x for x in os.listdir(lvl) if x.startswith("Cell_D_")):
        blob = open(os.path.join(lvl, fn), "rb").read()
        pos = 0
        while True:
            m = HDR.search(blob, pos)
            if not m: break
            lo = tuple(int(m.group(i)) for i in (2, 3, 4)); hi = tuple(int(m.group(i)) for i in (5, 6, 7))
            nc = int(m.group(11))
            nx, ny, nz = (hi[0]-lo[0]+1), (hi[1]-lo[1]+1), (hi[2]-lo[2]+1)
            npts = nx*ny*nz
            start = m.end()
            for c in comps:
                ci = idx[c]
                a = np.frombuffer(blob[start + ci*npts*8: start + (ci+1)*npts*8], dtype="<f8")
                a = a.reshape((nz, ny, nx)).transpose(2, 1, 0)   # i fastest
                out[c][lo[0]:hi[0]+1, lo[1]:hi[1]+1, lo[2]:hi[2]+1] = a
            pos = start + nc*npts*8
    return out

def plotfiles(leg):
    r = os.path.join(leg, "run")
    pfs = sorted(x for x in os.listdir(r) if x.endswith("cell"))
    return [(int(re.search(r"(\d+)cell", x).group(1)) * DT, os.path.join(r, x)) for x in pfs]

def nozzle_z(leg, t):
    head = None; best = None
    for line in open(os.path.join(leg, "run", "thermo.dat")):
        f = line.split()
        if head is None: head = f; continue
        try: tt = float(f[0])
        except ValueError: continue
        if best is None or abs(tt - t) < abs(best[0] - t):
            best = (tt, float(f[head.index("nozzle_z")]))
    return best[1] if best else None

def profile(removed):
    """equivalent-area diameter per level [mm] vs depth [mm] (quarter x4)."""
    depth, dia = [], []
    for k in range(NZ):
        c = removed[:, :, k].sum()
        if c > 0:
            A = 4.0 * c * DX * DX
            depth.append((NZ - 1 - k) * DX * 1e3 + 1.0); dia.append(2e3 * math.sqrt(A / math.pi))
    return depth, dia

def main():
    times = sorted({t for _, leg in LEGS for t, _ in plotfiles(leg)})
    times = [t for t in times if t > 0]
    data = {}   # (leg, t) -> arrays (cached per leg: last available <= t)
    pfmap = {leg: dict(plotfiles(leg)) for _, leg in LEGS}
    out = ROOT + "/Claude_markdowns/movies/frames"; os.makedirs(out, exist_ok=True)
    frames = []
    cmap = plt.get_cmap("inferno").copy()
    for fi, t in enumerate(times):
        fig = plt.figure(figsize=(17, 7.6), dpi=100, facecolor="white")
        gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 1.15], wspace=0.25, left=0.04, right=0.935, top=0.88, bottom=0.09)
        axp = fig.add_subplot(gs[0, 3])
        for li, (label, leg) in enumerate(LEGS):
            ax = fig.add_subplot(gs[0, li])
            avail = [tt for tt in pfmap[leg] if tt <= t and tt > 0]
            if not avail:
                ax.set_axis_off(); continue
            tt = max(avail); held = tt < t - 1e-6
            key = (leg, tt)
            if key not in data:
                data[key] = load(pfmap[leg][tt], ["Temp", "removed"])
            T = data[key]["Temp"][:, 0, :]; R = data[key]["removed"][:, 0, :] > 0.5
            # mirror about x = 0 (quarter -> half section); rows: depth increasing downward
            Tm = np.concatenate([T[::-1, :], T], axis=0)          # x from -120..120 mm
            Rm = np.concatenate([R[::-1, :], R], axis=0)
            img = Tm.T[::-1, :]                                    # depth down the page: k high = top
            msk = Rm.T[::-1, :]
            rgba = cmap((np.clip(img, 293, 900) - 293) / (900 - 293))
            rgba[msk] = (0.84, 0.90, 0.96, 1.0)                   # void = pale blue
            ax.imshow(rgba, extent=[-120, 120, 600, 0], aspect="equal", interpolation="nearest")
            # Meier outline
            dM = [d for d, _ in MEIER]; rM = [w / 2 for _, w in MEIER]
            ax.plot(rM, dM, "w--", lw=1.2); ax.plot([-r for r in rM], dM, "w--", lw=1.2)
            # burner: nozzle plane, feet, Ø67 mantle above the flange, Ø80 skirt below
            nz = nozzle_z(leg, tt)
            if nz is not None:
                dn = (Z_TOP - nz) * 1e3                            # nozzle depth [mm]
                ax.add_patch(Rectangle((-33.5, -40), 67, max(0.0, dn + 40), fc=(0.55, 0.6, 0.65, 0.55), ec="k", lw=0.6))
                ax.add_patch(Rectangle((-40, dn), 80, 50, fc=(0.4, 0.45, 0.5, 0.45), ec="k", lw=0.6))
                ax.plot([-3.75, 3.75], [dn, dn], "k-", lw=2)
            ax.set_xlim(-120, 120); ax.set_ylim(500, -40)
            ax.set_xlabel("x [mm]  (quarter domain mirrored)"); 
            if li == 0: ax.set_ylabel("depth below the original surface [mm]")
            ax.set_title(label + ("   [held at %d s]" % tt if held else ""), fontsize=10)
            d, w = profile(data[key]["removed"])
            axp.plot(w, d, lw=1.6, label=label.split("(")[0].strip() + (" (held)" if held else ""))
        axp.plot([w for _, w in MEIER], [d for d, _ in MEIER], "k--", lw=1.5, label="Meier Fig. 8.8")
        axp.axvline(80, color="grey", lw=0.8, ls=":"); axp.text(81, 470, "Ø80 skirt", fontsize=8, color="grey")
        axp.set_xlim(0, 140); axp.set_ylim(500, 0); axp.grid(alpha=0.3)
        axp.set_xlabel("hole diameter, equivalent-area [mm]"); axp.set_ylabel("depth [mm]")
        axp.set_title("Ø(z) vs Meier", fontsize=10); axp.legend(fontsize=8, loc="lower right")
        fig.suptitle("Meier flame-jet drilling, Grimsel granite — t = %d s     (colour: rock temperature 293–900 K; pale blue: removed rock; white dashes: Meier Fig. 8.8)" % t, fontsize=12)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(293, 900)); sm.set_array([])
        cax = fig.add_axes([0.955, 0.12, 0.011, 0.74]); fig.colorbar(sm, cax=cax, label="rock temperature [K]")
        fn = os.path.join(out, "frame_%03d.png" % fi); fig.savefig(fn); plt.close(fig)
        frames.append(fn); print("frame", fi, "t =", t, flush=True)
    imgs = [Image.open(f).convert("P", palette=Image.ADAPTIVE, colors=256) for f in frames]
    gif = ROOT + "/Claude_markdowns/movies/meier_sections.gif"
    imgs[0].save(gif, save_all=True, append_images=imgs[1:], duration=900, loop=0)
    print("wrote", gif, "frames", len(imgs))

if __name__ == "__main__":
    main()
