#!/usr/bin/env python3
"""Close the energy budget on an MMWSpalling removal run from its plotfiles.

E_in(t)      = P_flux * A_lit * t              (uniform disk, only lit columns)
E_stored(t)  = sum over solid cells of (H - H0) * Vcell   (from H field)
E_stored_T   = sum over solid cells of rho*Cp*(T - T0) * Vcell  (cross-check)
E_removed(t) = E_in - E_stored - E_loss   (loss bounded, reported)
<dT_rem>     = E_removed / (V_removed * rho * Cp)   -- must be <= Tmax - T0
"""
import sys, glob, os, re
import numpy as np
import yt
yt.funcs.mylog.setLevel(50)

out   = sys.argv[1]
P0    = float(sys.argv[2])      # W
rbeam = float(sys.argv[3])      # m
x0, y0 = float(sys.argv[4]), float(sys.argv[5])
rho, Cp, T0 = 2750.0, 790.0, 293.15
eps, sigSB, hconv = 0.8, 5.670374e-8, 10.0

flux = P0 / (np.pi * rbeam**2)
pfs = sorted(glob.glob(os.path.join(out, "[0-9]*cell")),
             key=lambda p: int(re.search(r"(\d+)cell", p).group(1)))
print(f"flux = {flux/1e6:.3f} MW/m^2   rho*Cp = {rho*Cp:.3e} J/m^3K")
print(f"{'t[s]':>6} {'A_lit[cm2]':>10} {'E_in[kJ]':>9} {'E_stH[kJ]':>9} {'E_stT[kJ]':>9} "
      f"{'V_rem[cm3]':>10} {'E_rem[kJ]':>9} {'<dT_rem>':>8} {'Tmax-T0':>8} {'closure':>12} {'zc_top[mm]':>10}")
H0 = None
for pf in pfs:
    ds = yt.load(pf)
    ad = ds.all_data()
    t = float(ds.current_time)
    dx = np.array(ds.domain_width / ds.domain_dimensions, dtype=float)
    Vc = dx.prod()
    x = ad["boxlib", "x"].v; y = ad["boxlib", "y"].v; z = ad["boxlib", "z"].v
    T = ad["boxlib", "Temp"].v; H = ad["boxlib", "H"].v
    rem = ad["boxlib", "removed"].v > 0.5
    if H0 is None:
        H0 = np.median(H[~rem])
    lit = ((x - x0)**2 + (y - y0)**2) <= rbeam**2
    # lit columns = unique (x,y) with lit
    cols = {}
    for xi, yi, li in zip(x, y, lit):
        cols[(round(xi, 9), round(yi, 9))] = li
    A_lit = sum(cols.values()) * dx[0] * dx[1]
    E_in = flux * A_lit * t
    E_stH = np.sum((H[~rem] - H0)) * Vc
    E_stT = np.sum(rho * Cp * (T[~rem] - T0)) * Vc
    V_rem = rem.sum() * Vc
    # loss bound: top solid cell per column radiating/convecting; use Tmax as bound over lit area
    Tmax = T[~rem].max()
    E_loss_bound = (eps * sigSB * Tmax**4 + hconv * (Tmax - T0)) * A_lit * t
    E_rem = E_in - E_stH                       # inferred from balance
    E_remD = np.sum((H[rem] - H0)) * Vc        # DIRECT: void cells keep frozen H
    dT_rem = E_remD / (V_rem * rho * Cp) if V_rem > 0 else float("nan")
    closure = (E_stH + E_remD) / E_in if E_in > 0 else float("nan")
    # central column: highest removed z
    cen = np.argmin((x - x0)**2 + (y - y0)**2)
    xc, yc = x[cen], y[cen]
    colmask = (np.abs(x - xc) < 1e-9) & (np.abs(y - yc) < 1e-9)
    zsolid = z[colmask & ~rem]
    ztop = zsolid.max() if zsolid.size else float("nan")
    print(f"{t:6.1f} {A_lit*1e4:10.2f} {E_in/1e3:9.2f} {E_stH/1e3:9.2f} {E_stT/1e3:9.2f} "
          f"{V_rem*1e6:10.2f} {E_remD/1e3:9.2f} {dT_rem:8.0f} {Tmax-T0:8.0f} {closure:12.3f} {ztop*1e3:10.2f}"
          f"   (E_rem by balance {E_rem/1e3:.2f} kJ; loss bound {E_loss_bound/1e3:.2f} kJ)")
