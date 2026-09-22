"""Top-cell heat budget vs stall state. For every disk column at a plotfile time:
T_top, T_below (one cell down), T_wall = coldest solid lateral neighbour cell at
k_top (a neighbour whose column top is higher), q_down = k (T_top - T_below)/dz,
q_lat = k sum_walls (T_top - T_wall)/dx, and the firing state from the removal
log (fired within the last 10 s before the plotfile time). Prints medians per
2 mm band and a 2x2 check of the rule 'stalled <=> q_pin - q_down - q_lat <= 0
at T_top = T_pin'."""
import sys, numpy as np, yt
yt.set_log_level(50)
run, pfname = sys.argv[1], sys.argv[2]
K, H, TG, T_PIN = 1.5, 700.0, 1600.0, 821.0
ds = yt.load(f"output/{run}/{pfname}")
dims = ds.domain_dimensions
cg = ds.covering_grid(0, left_edge=ds.domain_left_edge, dims=dims)
T = np.array(cg["Temp"]); rem = np.array(cg["removed"]) > 0.5
dz = float((ds.domain_right_edge[2] - ds.domain_left_edge[2]) / dims[2]); dx = dz
t = float(ds.current_time); zn = 0.170 - 4.3056e-4 * t
nx = dims[0]
ktop = np.full((nx, nx), -1)
for i in range(nx):
    for j in range(nx):
        s = np.nonzero(~rem[i, j, :])[0]
        if len(s): ktop[i, j] = s.max()
c = (np.arange(nx) + 0.5) * dx
r = np.hypot(*np.meshgrid(c, c, indexing="ij"))
a = np.loadtxt(f"output/{run}_removal_events.csv", delimiter=",", skiprows=1, ndmin=2)
last = np.zeros((nx, nx))
np.maximum.at(last, (a[:, 1].astype(int), a[:, 2].astype(int)), a[:, 0])
q_pin = H * (TG - T_PIN) - (0.8 * 5.670374e-8 * (T_PIN ** 4 - 293.15 ** 4) + 10 * (T_PIN - 293.15))
rows = []
for i in range(nx):
    for j in range(nx):
        if r[i, j] >= 0.030: continue
        k = ktop[i, j]; zf = (k + 1) * dz; s = zn - zf
        Tt = T[i, j, k]; Tb = T[i, j, k - 1]
        walls = []
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            ii, jj = i + di, j + dj
            if ii < 0 or jj < 0 or ii >= nx or jj >= nx: continue
            if ktop[ii, jj] > k and not rem[ii, jj, k]: walls.append(T[ii, jj, k])
        Tw = min(walls) if walls else np.nan
        q_down = K * (T_PIN - Tb) / dz
        q_lat = K * sum(T_PIN - w for w in walls) / dx
        firing = last[i, j] > t - 10.0
        rows.append((r[i, j], s, Tt, Tb, Tw, len(walls), q_down, q_lat, firing))
R = np.array(rows, dtype=float)
print(f"{run} {pfname} t={t:.0f}s  q_pin(821 K) = {q_pin/1e3:.0f} kW/m², dx {dx*1e3:g} mm")
print("band[mm]  n  firing  s_med[mm]  T_top  T_below  T_wall  nwalls  q_down  q_lat  margin(q_pin-q_down-q_lat) [kW/m²]")
for lo in range(0, 30, 2):
    sel = (R[:, 0] >= lo * 1e-3) & (R[:, 0] < (lo + 2) * 1e-3)
    if not sel.any(): continue
    b = R[sel]
    med = np.nanmedian
    print(f"{lo:2d}-{lo+2:<3d} {sel.sum():3d}  {b[:,8].mean()*100:4.0f}%  {med(b[:,1])*1e3:7.1f}  {med(b[:,2]):5.0f}  {med(b[:,3]):6.0f}  {med(b[:,4]):6.0f}  {med(b[:,5]):4.1f}  {med(b[:,6])/1e3:6.0f} {med(b[:,7])/1e3:6.0f}  {(q_pin-med(b[:,6])-med(b[:,7]))/1e3:6.0f}")
below = R[R[:, 1] > 0]           # columns still under the disk flux
marg = q_pin - below[:, 6] - below[:, 7]
f = below[:, 8] > 0.5
print(f"columns below the plane: {len(below)}; rule margin<=0 -> stalled: "
      f"stalled&margin<=0 {int(((~f)&(marg<=0)).sum())}, stalled&margin>0 {int(((~f)&(marg>0)).sum())}, "
      f"firing&margin<=0 {int((f&(marg<=0)).sum())}, firing&margin>0 {int((f&(marg>0)).sum())}")
