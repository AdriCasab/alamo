"""Temperature at the step between a firing band and its stalled outer neighbour.
Row j = 0 (y = dz/2); for columns i in a radial range print k_top, T_top, the
top cell's outward lateral neighbour temperature at the same k (the wall cell),
the outer column's surface T, and the step height."""
import sys, numpy as np, yt
yt.set_log_level(50)
run, pfname = sys.argv[1], sys.argv[2]
ds = yt.load(f"output/{run}/{pfname}")
dims = ds.domain_dimensions
cg = ds.covering_grid(0, left_edge=ds.domain_left_edge, dims=dims)
T = np.array(cg["Temp"]); rem = np.array(cg["removed"])
dz = float((ds.domain_right_edge[2] - ds.domain_left_edge[2]) / dims[2])
lz = float(ds.domain_right_edge[2]); t = float(ds.current_time)
zn = 0.170 - 4.3056e-4 * t
def ktop(i, j):
    col = rem[i, j, :]; solid = np.nonzero(col <= 0.5)[0]
    return solid.max() if len(solid) else -1
j = 0
print(f"{run} {pfname} t={t:.1f}s plane z={zn*1e3:.1f} mm (dz {dz*1e3:g} mm)")
print(" i   r[mm]  k_top  z_face[mm]  s[mm]   T_top   T_wall(i+1,k_top)  T_surf(i+1)  step[mm]  T_1below")
i_lo, i_hi = int(sys.argv[3]), int(sys.argv[4])
for i in range(i_lo, i_hi + 1):
    k = ktop(i, j); k1 = ktop(i + 1, j)
    r = (i + 0.5) * dz; zf = (k + 1) * dz; s = zn - zf
    Tw = T[i + 1, j, k] if rem[i + 1, j, k] <= 0.5 else float('nan')
    print(f"{i:2d}  {r*1e3:5.1f}  {k:4d}   {zf*1e3:7.1f}  {s*1e3:6.1f}  {T[i,j,k]:6.0f}   {Tw:8.0f}            {T[i+1,j,k1]:6.0f}    {(k1-k)*dz*1e3:5.1f}    {T[i,j,k-1]:6.0f}")
