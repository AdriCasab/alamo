#!/usr/bin/env python3
"""P2-8 (b2 diagnosis) -- WHERE the zero-gradient (w = 1 fallback) faces are, from a
leg's final void mask, with the code's own rule rebuilt independently (solid
indicator; above the domain top = gas, elsewhere outside = nearest cell; central
difference; zero gradient => fallback). Classifies each fallback cell by geometry."""
import math, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import reduce_p2_8 as R
leg = sys.argv[1] if len(sys.argv) > 1 else "L0"
d = R.legdir(leg); g = R.geom_of(d); pf = R.plotfiles(d)[-1]
rem = R.read_pf(d, pf, g)["removed"]
nx, ny, nz = rem.shape
S = (rem <= 0.5).astype(float)
P = np.pad(S, 1, mode="edge"); P[:, :, -1] = 0.0          # above the top = gas
gx = P[2:, 1:-1, 1:-1] - P[:-2, 1:-1, 1:-1]
gy = P[1:-1, 2:, 1:-1] - P[1:-1, :-2, 1:-1]
gz = P[1:-1, 1:-1, 2:] - P[1:-1, 1:-1, :-2]
solid = S > 0.5
V = ~solid
op = np.zeros((6,) + S.shape, bool)
op[0, 1:] = V[:-1]; op[1, :-1] = V[1:]; op[2, :, 1:] = V[:, :-1]; op[3, :, :-1] = V[:, 1:]
op[4, :, :, 1:] = V[:, :, :-1]; op[5, :, :, :-1] = V[:, :, 1:]; op[5, :, :, -1] = True
op &= solid[None]
nface = op.sum(0)
zero = (gx == 0) & (gy == 0) & (gz == 0) & (nface > 0)
xc = (np.arange(nx) + .5) * g["dx"]; X, Y = np.meshgrid(xc, xc, indexing="ij"); Rr = np.hypot(X, Y)
roof = zero & op[4] & op[5]                      # void below AND open above: a one-cell plate
print(f"{leg} {pf}: surface cells {int((nface>0).sum())}, faces {int(nface.sum())}")
print(f"  zero-gradient (fallback) cells {int(zero.sum())}, their faces {int(nface[zero].sum())}"
      f" = {nface[zero].sum()/nface.sum()*100:.2f} % of faces")
print(f"  of which ONE-CELL PLATES (void below, open above): {int(roof.sum())} cells,"
      f" {int(nface[roof].sum())} faces")
ks = np.nonzero(zero)
rr = Rr[ks[0], ks[1]] * 1e3
print(f"  fallback cells: k {sorted(set(ks[2].tolist()))} (top k = {nz-1}),"
      f" r {rr.min():.1f}-{rr.max():.1f} mm")
