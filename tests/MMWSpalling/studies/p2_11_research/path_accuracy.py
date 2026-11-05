#!/usr/bin/env python3
"""P2-11 planning check: accuracy of graph shortest paths over SURFACE CELLS (binary staircase), neighbours within
radius rho cells (rho 1 = 6-connected, sqrt3 = 26-connected), Euclidean edge weights, vs the analytic surface distance. Shapes: flat floor (s = r), planes tilted
20/35/45 deg in x-z (s = in-plane distance), a flat-floored cylindrical bore R 12, depth 24 (s = r on the floor,
R + height up the wall, ignoring the corner). Source: the surface cell nearest the analytic stagnation point; its
initial distance = Euclidean distance from the point. Reports median / p95 / max |rel err| over cells with s >= 4 cells."""
import heapq, math, itertools
import numpy as np

def surface_cells(solid):
    pad = np.pad(solid, 1, constant_values=True)   # domain edges and symmetry planes are not gas
    out = np.zeros_like(solid)
    for a, b, c in [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]:
        nb = pad[1+a:pad.shape[0]-1+a, 1+b:pad.shape[1]-1+b, 1+c:pad.shape[2]-1+c]
        out |= solid & ~nb

    return out

def dijkstra(surf, src, d0, rho):
    offs = [o for o in itertools.product(range(-int(rho), int(rho)+1), repeat=3) if o != (0,0,0) and math.sqrt(sum(q*q for q in o)) <= rho + 1e-9]
    dist = np.full(surf.shape, np.inf); dist[src] = d0; pq = [(d0, src)]
    while pq:
        d, (i, j, k) = heapq.heappop(pq)
        if d > dist[i, j, k]: continue
        for a, b, c in offs:
            p, q, r = i+a, j+b, k+c
            if 0 <= p < surf.shape[0] and 0 <= q < surf.shape[1] and 0 <= r < surf.shape[2] and surf[p, q, r]:
                nd = d + math.sqrt(a*a+b*b+c*c)
                if nd < dist[p, q, r]: dist[p, q, r] = nd; heapq.heappush(pq, (nd, (p, q, r)))
    return dist

def report(name, solid, stag, exact_fn, rhos=(math.sqrt(3), 2, 3, 4)):
    surf = surface_cells(solid)
    # restrict "surface" to cells exposed upward/sideways to gas (not through the domain bottom/sides)
    idx = np.argwhere(surf)
    cen = idx + 0.5
    src = tuple(idx[np.argmin(((cen - stag) ** 2).sum(1))])
    d0 = math.sqrt(((np.array(src) + 0.5 - stag) ** 2).sum())
    ex = exact_fn(cen)
    for rho in rhos:
        dist = dijkstra(surf, src, d0, rho)
        dv = dist[tuple(idx.T)]
        out = []
        for lo in (4, 10):
            m = (ex >= lo) & np.isfinite(dv)
            e = np.abs(dv[m] / ex[m] - 1)
            out.append("s>=%2d: median %.3f p95 %.3f max %.3f (n %d)" % (lo, np.median(e), np.percentile(e, 95), e.max(), m.sum()))
        print("  %-26s rho %.2f: %s" % (name, rho, " | ".join(out)))

N = 40
X, Y, Z = np.meshgrid(np.arange(N) + .5, np.arange(N) + .5, np.arange(N) + .5, indexing="ij")
# flat floor, quarter domain, axis at the corner (0,0): exposed top at z = 20; sides of the domain are not gas
solid = Z < 20
stag = np.array([0.0, 0.0, 20.0])
surf_only_top = solid & ~(np.roll(solid, -1, 2))
report("flat floor (quarter)", solid, stag, lambda c: np.hypot(c[:, 0], c[:, 1]))
for th in (20, 35, 45):
    t = math.radians(th)
    n = np.array([math.sin(t), 0, math.cos(t)])
    solid = (X * n[0] + Z * n[2]) < (20 * n[2] + 10 * n[0])          # plane through (10, *, 20)
    stag = np.array([10.0, 0.0, 20.0])
    def ex(c, n=n):
        d = c - stag; d = d - np.outer(d @ n, n)                    # in-plane displacement
        return np.sqrt((d ** 2).sum(1))
    report("plane %d deg (quarter)" % th, solid, stag, ex)
# bore: R 12, depth 24 below z = 36, quarter domain
R, top, dep = 12.0, 36.0, 24.0
rr = np.hypot(X, Y)
solid = (Z < top) & ~((rr < R) & (Z > top - dep))
stag = np.array([0.0, 0.0, top - dep])
def exb(c):
    r = np.hypot(c[:, 0], c[:, 1])
    zc = c[:, 2]
    return np.where(zc > top - 1.0, R + dep + np.maximum(0.0, r - R),       # the original top, beyond the mouth
           np.where(r < R - 0.5, r, R + np.maximum(0.0, zc - (top - dep))))
report("bore R12 depth24 (quarter)", solid, stag, exb)
