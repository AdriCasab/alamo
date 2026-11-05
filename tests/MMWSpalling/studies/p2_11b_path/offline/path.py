#!/usr/bin/env python3
"""P2-11b (a0) -- THE SURFACE-PATH ROUTINE TO BE PORTED, and its accuracy on the research fixtures.

surface_path(nodes, sources, rho=4, dx=1.0):
  nodes    bool (nx, ny, nz): the graph's cells -- listed cells with at least one UNBLOCKED open face
  sources  {(i, j, k): d0}: start cells and their initial distance (X1: the cell under the +z face of each
           column in the march's innermost non-empty r-bin; d0 = that face centre's radius from (x0, y0))
  edges    every pair of nodes whose index offset o has |o| <= rho (rho = 4 cells), weight dx |o| (Euclidean)
  returns  dist (nx, ny, nz), inf where unreached
Plain binary-heap Dijkstra, exactly what the C++ port does (a fixed offset table, lazy deletion).

a0: on `studies/p2_11_research/path_accuracy.py`'s fixtures (flat floor, planes 20/35/45 deg, bore R12 depth 24;
quarter domains, single source = the surface cell nearest the analytic stagnation point, d0 its Euclidean distance),
|s_path / s_exact - 1| over cells with s_exact >= 10 cells: median <= 5 %, p95 <= 8 % on EVERY fixture.
Writes A0.txt beside the study.
"""
import heapq
import itertools
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def offsets(rho=4):
    R = int(math.floor(rho))
    return [(o, math.sqrt(o[0] * o[0] + o[1] * o[1] + o[2] * o[2]))
            for o in itertools.product(range(-R, R + 1), repeat=3)
            if o != (0, 0, 0) and o[0] * o[0] + o[1] * o[1] + o[2] * o[2] <= rho * rho + 1e-9]


def surface_path(nodes, sources, rho=4, dx=1.0):
    OFF = offsets(rho)
    nx, ny, nz = nodes.shape
    dist = np.full(nodes.shape, np.inf)
    pq = []
    for c, d0 in sources.items():
        if nodes[c] and d0 < dist[c]:
            dist[c] = d0
            heapq.heappush(pq, (d0, c))
    while pq:
        d, (i, j, k) = heapq.heappop(pq)
        if d > dist[i, j, k]:
            continue
        for (a, b, c), L in OFF:
            p, q, r = i + a, j + b, k + c
            if 0 <= p < nx and 0 <= q < ny and 0 <= r < nz and nodes[p, q, r]:
                nd = d + dx * L
                if nd < dist[p, q, r]:
                    dist[p, q, r] = nd
                    heapq.heappush(pq, (nd, (p, q, r)))
    return dist


# ------------------------------------------------------------------ a0 fixtures (path_accuracy.py's) ----
def surface_cells(solid):
    pad = np.pad(solid, 1, constant_values=True)   # domain edges and symmetry planes are not gas
    out = np.zeros_like(solid)
    for a, b, c in [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]:
        nb = pad[1 + a:pad.shape[0] - 1 + a, 1 + b:pad.shape[1] - 1 + b, 1 + c:pad.shape[2] - 1 + c]
        out |= solid & ~nb
    return out


def score(name, solid, stag, exact_fn, L):
    surf = surface_cells(solid)
    idx = np.argwhere(surf)
    cen = idx + 0.5
    src = tuple(int(v) for v in idx[np.argmin(((cen - stag) ** 2).sum(1))])
    d0 = math.sqrt(((np.array(src) + 0.5 - stag) ** 2).sum())
    dist = surface_path(surf, {src: d0}, rho=4)
    ex = exact_fn(cen)
    dv = dist[tuple(idx.T)]
    m = (ex >= 10) & np.isfinite(dv)
    e = np.abs(dv[m] / ex[m] - 1)
    med, p95 = float(np.median(e)), float(np.percentile(e, 95))
    ok = med <= 0.05 and p95 <= 0.08
    L.append(f"  {name:28s} s >= 10 cells: median {100 * med:5.2f} %  p95 {100 * p95:5.2f} %  max {100 * e.max():5.2f} %"
             f"  (n {m.sum()}, unreached {int((~np.isfinite(dv)).sum())})  {'ok' if ok else 'OUT'}")
    return ok


def main():
    L = ["P2-11b (a0) -- the surface path (Dijkstra, neighbours within 4 cells, Euclidean edges) vs the analytic "
         "surface distance", "=" * 110]
    N = 40
    X, Y, Z = np.meshgrid(np.arange(N) + .5, np.arange(N) + .5, np.arange(N) + .5, indexing="ij")
    ok = score("flat floor (quarter)", Z < 20, np.array([0.0, 0.0, 20.0]), lambda c: np.hypot(c[:, 0], c[:, 1]), L)
    for th in (20, 35, 45):
        t = math.radians(th)
        n = np.array([math.sin(t), 0, math.cos(t)])
        solid = (X * n[0] + Z * n[2]) < (20 * n[2] + 10 * n[0])
        stag = np.array([10.0, 0.0, 20.0])

        def ex(c, n=n, stag=stag):
            d = c - stag
            d = d - np.outer(d @ n, n)
            return np.sqrt((d ** 2).sum(1))
        ok &= score("plane %d deg (quarter)" % th, solid, stag, ex, L)
    R, top, dep = 12.0, 36.0, 24.0
    rr = np.hypot(X, Y)
    solid = (Z < top) & ~((rr < R) & (Z > top - dep))
    stag = np.array([0.0, 0.0, top - dep])

    def exb(c):
        r = np.hypot(c[:, 0], c[:, 1])
        zc = c[:, 2]
        return np.where(zc > top - 1.0, R + dep + np.maximum(0.0, r - R),
                        np.where(r < R - 0.5, r, R + np.maximum(0.0, zc - (top - dep))))
    ok &= score("bore R12 depth24 (quarter)", solid, stag, exb, L)
    L.append(f"a0 {'PASS' if ok else 'FAIL'} (bar: median <= 5 %, p95 <= 8 % on every fixture; research rho 4: "
             f"flat 1.0/2.5, planes 3.3/6.6, 3.3/5.8, 3.1/5.7, bore 4.5/6.6 %)")
    txt = "\n".join(L)
    print(txt)
    open(os.path.join(os.path.dirname(HERE), "A0.txt"), "w").write(txt + "\n")


if __name__ == "__main__":
    main()
