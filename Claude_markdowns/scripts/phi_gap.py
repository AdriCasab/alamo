#!/usr/bin/env python3
"""Test the phi/mask-gap hypothesis for the beam-energy leak.

For every lit column, take the topmost solid (removed<0.5) cell, read its phi
(the optical path BeamSource uses), form s_top = max(0, phi - dz/2) and the
implied column deposition fraction f = exp(-alpha * s_top).  If <f> over lit
columns matches the instantaneous energy closure, the leak is the mask/phi gap.
"""
import sys, glob, os, re
import numpy as np
import yt
yt.funcs.mylog.setLevel(50)

out   = sys.argv[1]
rbeam = float(sys.argv[2]); x0, y0 = float(sys.argv[3]), float(sys.argv[4])
alpha = 2000.0
pfs = sorted(glob.glob(os.path.join(out, "[0-9]*cell")),
             key=lambda p: int(re.search(r"(\d+)cell", p).group(1)))
print(f"{'t[s]':>6} {'ncol':>5} {'<phi_top>[mm]':>13} {'<s_top>[mm]':>11} {'<f_dep>':>8} {'min f':>6} {'max f':>6}  phi_top histogram (mm)")
for pf in pfs[1:]:
    ds = yt.load(pf); ad = ds.all_data(); t = float(ds.current_time)
    dx = np.array(ds.domain_width / ds.domain_dimensions, dtype=float); dz = dx[2]
    x = ad["boxlib","x"].v; y = ad["boxlib","y"].v; z = ad["boxlib","z"].v
    phi = ad["boxlib","phi"].v; rem = ad["boxlib","removed"].v > 0.5
    lit = ((x-x0)**2 + (y-y0)**2) <= rbeam**2
    cols = {}
    for n in np.where(lit & ~rem)[0]:
        key = (round(x[n],9), round(y[n],9))
        if key not in cols or z[n] > cols[key][0]:
            cols[key] = (z[n], phi[n])
    phit = np.array([v[1] for v in cols.values()])
    s_top = np.maximum(0.0, phit - 0.5*dz)
    f = np.exp(-alpha * s_top)
    hist, edges = np.histogram(phit*1e3, bins=[-1,0,0.5,1.0,1.5,2.0,2.5,3.0,10])
    print(f"{t:6.1f} {len(phit):5d} {phit.mean()*1e3:13.3f} {s_top.mean()*1e3:11.3f} {f.mean():8.3f} {f.min():6.3f} {f.max():6.3f}  "
          + " ".join(f"[{edges[i]:.1f},{edges[i+1]:.1f}):{hist[i]}" for i in range(len(hist))))
