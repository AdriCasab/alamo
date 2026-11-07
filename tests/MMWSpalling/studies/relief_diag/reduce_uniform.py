#!/usr/bin/env python3
"""Reduce the uniform-heating resolution ladder: centre-column sigma_op / sigma_1D = R * Sp_M / Sp_L at each
plot step for dx = 2, 1, 0.5 mm, with the skin depth ell = dx / ln((T_top - T0)/(T_1 - T0)) read from the
centre column's top two cells (the same two-cell estimate the production reconstruction uses)."""
import os, sys, math
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import reduce as RD
R = RD.R
T0 = 293.15
L_BLOCK = 0.024

def column(plot, name, dx):
    """values down the centre column, top first"""
    hdr = os.path.join(plot, "Header"); lines = open(hdr).read().split("\n")
    n = int(lines[1]); names = lines[2:2 + n]; t = float(lines[2 + n + 1]); c = names.index(name)
    ic = int(round(L_BLOCK / dx / 2)) - 1
    col = {}
    lvl = os.path.join(plot, "Level_0")
    for fn in sorted(os.listdir(lvl)):
        if not fn.startswith("Cell_D_"): continue
        for f in RD.fabcmp.read_fabs_boxed(os.path.join(lvl, fn)):
            lo, hi = f["lo"], f["hi"]
            if not (lo[0] <= ic <= hi[0] and lo[1] <= ic <= hi[1]): continue
            nx = hi[0] - lo[0] + 1; ny = hi[1] - lo[1] + 1
            comp = f["comp"][c]
            for k in range(lo[2], hi[2] + 1):
                col[k] = comp[(k - lo[2]) * nx * ny + (ic - lo[1]) * nx + (ic - lo[0])]
    ks = sorted(col, reverse=True)
    return [col[k] for k in ks], t

def main():
    out = ["uniform-heating resolution ladder: centre column, sigma_op / sigma_1D = R * Sp_M / Sp_L (exact answer 1)",
           "%-5s %-10s %7s | %9s %9s | %8s %8s | %8s" % ("dx", "step", "t [s]", "Sp_L", "Sp_M", "T_top", "ell/dx", "ratio")]
    for dxk, dx in (("2", 0.002), ("1", 0.001), ("05", 0.0005)):
        dl = os.path.join(HERE, "U", "Lc" + dxk, "run"); dm = os.path.join(HERE, "U", "M" + dxk, "run")
        if not (os.path.isdir(dl) and os.path.isdir(dm)): continue
        steps = sorted(set(p for p in os.listdir(dl) if p.endswith("cell")) & set(p for p in os.listdir(dm) if p.endswith("cell")))
        for s in steps:
            spl, t = column(os.path.join(dl, s), "Sp_field", dx)
            spm, _ = column(os.path.join(dm, s), "Sp_field", dx)
            T, _ = column(os.path.join(dl, s), "Temp", dx)
            ell = float("nan")
            if T[0] > T0 and T[1] > T0 and T[0] > T[1]:
                ell = dx / math.log((T[0] - T0) / (T[1] - T0))
            rat = R * spm[0] / spl[0] if abs(spl[0]) > 1e-12 else float("nan")
            out.append("%-5s %-10s %7.3f | %9.4f %9.4f | %8.2f %8.3f | %8.3f" % (dxk, s, t, spl[0], spm[0], T[0], ell / dx, rat))
    txt = "\n".join(out); print(txt)
    open(os.path.join(HERE, "RESULTS_uniform.txt"), "w").write(txt + "\n")

if __name__ == "__main__":
    main()
