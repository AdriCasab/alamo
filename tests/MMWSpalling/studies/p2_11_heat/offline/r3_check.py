#!/usr/bin/env python3
"""P2-11 extra (report only) -- why R3's offline LOCAL-Goldstein shares differ from the G1 run.
On P47n's OWN plotfiles (binary B), the conductance share within r < 10 mm and on floor faces for:
  (A) P47n as run: fl_h x fl_w from the plotfile;
  (B) the G1 formula applied the CODE's way: max over each cell's open faces of h(r_face, s_face), x fl_w;
  (C) R3's construction: g_loc at the CELL-CENTRE r, the TOP-face s for every cell, lid cells zeroed, x fl_w.
(B) vs (A) isolates the formula on a fixed shape; (C) vs (B) is R3's construction error; G1's measured share
(its own shape) vs (B) is the shape's feedback. Writes R3_CHECK.txt beside the study."""
import importlib.util
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("r11", os.path.join(HERE, "reduce_p2_11.py"))
r = importlib.util.module_from_spec(spec)
spec.loader.exec_module(r)
S_G1 = r.h_of(r.P11, "G1")[1]
CG1 = float("%.6f" % (S_G1 * r.H0G))


def h_g1(rr, s):
    if rr > 0.02 and s < 0.0633:
        return 0.001
    p = math.pow(rr / 0.0079, 1.285)
    return CG1 * 533 * (533 + 0.3575 * 44 * p) / math.pow(533 + 44 * p, 2) * math.pow(0.0632 / max(0.0632, s), 1.3185305455)


def main():
    mod, n = r.P9, "P47n"
    g, P = r.pfs_by_t(mod, n)
    zn = 0.2132
    L = [__doc__.strip().split("\n")[0], "=" * 110,
         "   t |  share r<10 mm: (A) P47n  (B) G1 code-way  (C) R3-way  G1 run | floor: (A)    (B)    (C)   G1 run"]
    for t in (15, 30, 60, 90, 120):
        f = mod.read_pf(mod.legdir(n), P[t], g, comps=("fr_f", "fl_h", "fl_w", "removed"))
        listed = f["fl_w"] > 0
        fo = r.faces(f["removed"])
        HB = np.zeros_like(f["fl_h"])
        HC = np.zeros_like(f["fl_h"])
        for (i, j, k) in np.argwhere(listed):
            best = 0.0
            for dr in range(6):
                if fo[dr][i, j, k]:
                    x, y, zf = r.face_pos(i, j, k, dr, g)
                    best = max(best, h_g1(math.sqrt(x * x + y * y), zn - zf))
            HB[i, j, k] = best
            xc, yc = (i + 0.5) * g["dx"], (j + 0.5) * g["dx"]
            rc = math.hypot(xc, yc)
            sf = zn - (k + 1) * g["dz"]
            HC[i, j, k] = 0.0 if (rc > 0.02 and sf < 0.0633) else h_g1(rc, sf)
        nz, _, zero = r.smoothed_nz(f["fr_f"])
        x = (np.arange(g["nx"]) + 0.5) * g["dx"]
        R = np.hypot(x[:, None], x[None, :])[:, :, None] * np.ones((1, 1, g["nz"]))
        fl = listed & (nz > math.cos(math.radians(45))) & ~zero
        out = []
        for H in (f["fl_h"], HB, HC):
            G = np.where(listed, H * f["fl_w"], 0.0)
            out.append((G[listed & (R < 0.010)].sum() / G.sum(), G[fl].sum() / G.sum()))
        m = r.conductance(r.P11, "G1", t)
        L.append(f"  {t:3d} |  {out[0][0]:.3f}  {out[1][0]:.3f}  {out[2][0]:.3f}  {m['axis']:.3f} | "
                 f"{out[0][1]:.3f}  {out[1][1]:.3f}  {out[2][1]:.3f}  {m['floor']:.3f}")
    txt = "\n".join(L)
    print(txt)
    open(os.path.join(HERE, "R3_CHECK.txt"), "w").write(txt + "\n")


if __name__ == "__main__":
    main()
