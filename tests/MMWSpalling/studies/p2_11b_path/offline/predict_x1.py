#!/usr/bin/env python3
"""P2-11b (a1, a2) -- THE BUILD DECISION, offline, before any source edit.

On G1's plotfiles (studies/p2_11_heat/legs/G1, binary C2) at 60 / 90 / 120 s, every open face of every listed cell
is evaluated at its own centre the code's way (reduce_p2_11.faces / face_pos), under
  G1  h = the deck's h_expr at (r_face, s_face)                       (what ran; reproduces fl_h, P2-11 a1)
  X1  h = C g_loc(s_path) stand(s_c), no lid; BLOCKED faces 0; unreached 0   (the packet's item 4)
      sources  the cell under the +z face of each column in the march's innermost non-empty r-bin (dr = dx,
               column-centre r, as JetEnthalpyMarch bins); d0 = that face centre's radius; H = s_c = mean of
               z_n - z_top over the bin
      barrier  (r > 0.02) and (z > 0.15) at the point 0.5 dx outside the face centre
      s_path   path.surface_path over listed cells with >= 1 unblocked open face, rho 4, Euclidean edges
Conductance per cell, PRIMARY = max over its faces of h x fl_w (P2-11 c2's definition and what a realised plotfile
gives: fl_h x fl_w, so c3 compares like with like); SECONDARY = (fl_w / n_faces) x sum of face h (the kernel's sum).
Both rules are scaled to G1's total conductance. Shares: within r < 10 mm (cell centre); the top-10-mm WALL band
(listed, smoothed n_z <= cos 45, cell centre <= 10 mm below the original top).
a1: X1/G1 within-10-mm ratio >= 1.10 at >= 2 of the 3 times => BUILD, else STOP.  Writes OFFLINE_X1.txt.
"""
import importlib.util
import math
import os
import sys
import time as _time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
S11 = os.path.dirname(HERE)
STUD = os.path.dirname(S11)
sys.path.insert(0, HERE)
import path as P  # noqa: E402

spec = importlib.util.spec_from_file_location("r11", os.path.join(STUD, "p2_11_heat", "reduce_p2_11.py"))
r11 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(r11)
D = 0.0079
ZTOP = 0.15


def g_loc(s):
    u = 44.0 * math.pow(s / D, 1.285)
    return 533.0 * (533.0 + 0.3575 * u) / (533.0 + u) ** 2


def stand(s):
    return math.pow(0.0632 / max(0.0632, s), 1.3185305455)


def barrier(x, y, z):
    return math.hypot(x, y) > 0.02 and z > 0.15


OUT = {0: (-1, 0, 0), 1: (1, 0, 0), 2: (0, -1, 0), 3: (0, 1, 0), 4: (0, 0, -1), 5: (0, 0, 1)}


def evaluate(t):
    mod, n = r11.P11, "G1"
    g, PF = r11.pfs_by_t(mod, n)
    d = mod.legdir(n)
    f = mod.read_pf(d, PF[t], g, comps=("fr_f", "fl_h", "fl_w", "removed"))
    hG1, S, _ = r11.h_of(mod, n)
    zn = float(r11.argv_kv(mod, n)["surface_patch.nozzle_z0"])
    dx, dz = g["dx"], g["dz"]
    listed = f["fl_w"] > 0
    fo = r11.faces(f["removed"])
    nfaces = sum(fo[dr].astype(int) for dr in range(6))
    # per-face geometry and the barrier
    blocked = {dr: np.zeros(listed.shape, bool) for dr in range(6)}
    for (i, j, k) in np.argwhere(listed):
        for dr in range(6):
            if fo[dr][i, j, k]:
                x, y, zf = r11.face_pos(i, j, k, dr, g)
                o = OUT[dr]
                if barrier(x + 0.5 * dx * o[0], y + 0.5 * dx * o[1], zf + 0.5 * dz * o[2]):
                    blocked[dr][i, j, k] = True
    unbl = sum((fo[dr] & ~blocked[dr]).astype(int) for dr in range(6))
    nodes = listed & (unbl > 0)
    # the march's innermost non-empty bin: column tops (top solid cell), column-centre r, dr = dx
    solid = ~(f["removed"] > 0.5)
    bins = {}
    for i in range(g["nx"]):
        for j in range(g["ny"]):
            ks = np.nonzero(solid[i, j, :])[0]
            if ks.size == 0:
                continue
            kt = int(ks.max())
            xc, yc = (i + 0.5) * dx, (j + 0.5) * dx
            rc = math.hypot(xc, yc)
            if rc * rc > 0.08 ** 2:
                continue
            bins.setdefault(int(math.floor(rc / dx)), []).append((i, j, kt, rc))
    b0 = bins[min(bins)]
    s_c = float(np.mean([zn - (kt + 1) * dz for (_, _, kt, _) in b0]))
    sources = {(i, j, kt): rc for (i, j, kt, rc) in b0}
    t0 = _time.time()
    dist = P.surface_path(nodes, sources, rho=4, dx=dx)
    tp = _time.time() - t0
    # per-face h under both rules -> per-cell max and sum
    Hmax = {"G1": np.zeros(listed.shape), "X1": np.zeros(listed.shape)}
    Hsum = {"G1": np.zeros(listed.shape), "X1": np.zeros(listed.shape)}
    unreached = 0
    ratio_sr = []
    for (i, j, k) in np.argwhere(listed):
        sp = dist[i, j, k]
        if nodes[i, j, k] and not np.isfinite(sp):
            unreached += 1
        for dr in range(6):
            if not fo[dr][i, j, k]:
                continue
            x, y, zf = r11.face_pos(i, j, k, dr, g)
            hg = hG1(math.sqrt(x * x + y * y), zn - zf)
            hx = 0.0 if (blocked[dr][i, j, k] or not np.isfinite(sp)) else g_loc(sp) * stand(s_c)
            for nm, h in (("G1", hg), ("X1", hx)):
                Hmax[nm][i, j, k] = max(Hmax[nm][i, j, k], h)
                Hsum[nm][i, j, k] += h
        if np.isfinite(sp) and (ZTOP - (k + 0.5) * dz) > 0.005:
            ratio_sr.append(sp / math.hypot((i + 0.5) * dx, (j + 0.5) * dx))
    nz, _, zero = r11.smoothed_nz(f["fr_f"])
    x = (np.arange(g["nx"]) + 0.5) * dx
    R = np.hypot(x[:, None], x[None, :])[:, :, None] * np.ones((1, 1, g["nz"]))
    depth_c = (ZTOP - (np.arange(g["nz"]) + 0.5) * dz)[None, None, :] * np.ones(listed.shape)
    floor = listed & (nz > math.cos(math.radians(45))) & ~zero
    band = listed & ~floor & (depth_c <= 0.010)
    ax10 = listed & (R < 0.010)
    res = {}
    for kind, H in (("max", Hmax), ("sum", Hsum)):
        for nm in ("G1", "X1"):
            G = H[nm] * f["fl_w"] if kind == "max" else H[nm] * f["fl_w"] / np.maximum(nfaces, 1)
            tot = G[listed].sum()
            res[(kind, nm)] = (G[ax10].sum() / tot, G[band].sum() / tot, G[floor].sum() / tot)
    # the G1 rule reproduces the plotfile (sanity: P2-11 a1)
    rel = np.abs(Hmax["G1"][listed] / f["fl_h"][listed] - 1)
    return dict(t=t, res=res, s_c=s_c, n_src=len(sources), n_blk=int(sum(blocked[dr].sum() for dr in range(6))),
                n_unr=unreached, n_nodes=int(nodes.sum()), n_listed=int(listed.sum()), tp=tp,
                sr=float(np.median(ratio_sr)) if ratio_sr else float("nan"),
                smax=float(np.nanmax(np.where(np.isfinite(dist), dist, np.nan))), g1ok=float((rel <= 1e-12).mean()))


def main():
    L = ["P2-11b (a1, a2) -- offline X1 vs G1 on G1's own plotfiles, per open face the code's way, equal total "
         "conductance", "=" * 116]
    rows = []
    for t in (60, 90, 120):
        e = evaluate(t)
        rows.append(e)
        L.append(f"  t {t:3d} s: listed {e['n_listed']}, graph nodes {e['n_nodes']}, sources {e['n_src']}, "
                 f"H = s_c {1e3 * e['s_c']:.1f} mm, blocked faces {e['n_blk']}, unreached {e['n_unr']}, "
                 f"s_path max {1e3 * e['smax']:.1f} mm, median s_path/r (cells > 5 mm deep) {e['sr']:.3f}; "
                 f"G1 rule vs plotfile fl_h within 1e-12 on {e['g1ok']:.4f}; Dijkstra {e['tp']:.1f} s (Python)")
    L.append("")
    L.append("  PRIMARY (per cell: max face h x fl_w -- P2-11 c2's definition, comparable with a realised fl_h x fl_w)")
    L.append("     t |  within 10 mm: G1    X1   X1/G1 | top-10-mm wall band: G1    X1   X1/G1 | floor: G1    X1")
    n_ok = 0
    for e in rows:
        a, b = e["res"][("max", "G1")], e["res"][("max", "X1")]
        r10 = b[0] / a[0]
        n_ok += r10 >= 1.10
        L.append(f"   {e['t']:3d} |        {a[0]:.3f}  {b[0]:.3f}  {r10:5.2f} |                    {a[1]:.3f}  {b[1]:.3f}"
                 f"  {b[1] / a[1]:5.2f} |   {a[2]:.3f}  {b[2]:.3f}")
    L.append("  SECONDARY (per cell: (fl_w / n_faces) x sum of face h -- the kernel's own sum)")
    for e in rows:
        a, b = e["res"][("sum", "G1")], e["res"][("sum", "X1")]
        L.append(f"   {e['t']:3d} |        {a[0]:.3f}  {b[0]:.3f}  {b[0] / a[0]:5.2f} |                    {a[1]:.3f}  {b[1]:.3f}"
                 f"  {b[1] / a[1]:5.2f} |   {a[2]:.3f}  {b[2]:.3f}")
    dec = "BUILD" if n_ok >= 2 else "STOP"
    L.append(f"\n  a1: X1/G1 within-10-mm ratio >= 1.10 at {n_ok} of 3 times (primary) => {dec}   (P0 predicted 1.15-1.45)")
    L.append(f"a1 {dec}")
    txt = "\n".join(L)
    print(txt)
    open(os.path.join(S11, "OFFLINE_X1.txt"), "w").write(txt + "\n")


if __name__ == "__main__":
    main()
