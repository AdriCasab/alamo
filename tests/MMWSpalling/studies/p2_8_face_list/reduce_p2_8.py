#!/usr/bin/env python3
"""P2-8 -- the reductions.

  reduce_p2_8.py GATE    b1 (K0 vs whatif_uptake/legs/s100, md5) + b2 (L0 witnesses)
                         -> GATE_b.txt
  reduce_p2_8.py MESH    b3 (L0c vs L0), c1 (L0 vs K0), c2 (L0 sectors)  -> MESH.txt
  reduce_p2_8.py C3      c3 (L0f vs L0)                                   -> C3.txt
  reduce_p2_8.py SSTAR   d1: S* from C15/C20/C25 (the packet's rule)      -> SSTAR.txt
  reduce_p2_8.py         every leg, per plotfile; d2-d4                   (stdout -> RESULTS.txt)

CONVENTIONS (each a standing rule or a prior estimator, named):
  ⚠ geometry from each leg's OWN argv.sh (SR 5 -- the meshes differ).
  ⚠ x4 on volume AND energy (SR 3); event-log volume and void-mask volume both.
  ⚠ cumulative ledger columns read at their value, per-step columns never summed
    (SR 7); `ledger_err_solid` is RELATIVE. b2's identities are the in-code DOUBLE
    residuals fl_id_march / fl_id_ledger (thermo's %g cannot resolve 1e-9).
  ⚠ the bore PER 30-degree SECTOR (SR 10): per level, each sector's max void-cell
    CENTRE radius (bore_profile.sector_radii's rule); a level's diameter is twice
    the sector MEAN (P2-4's convention), its spread max - min. A level index for a
    depth d below the original top is k = int((z_top - d - z_lo)/dz + 1e-9)
    (bore_profile's rule).
  ⚠ "max Ø" = the largest level diameter over all depths (the bulge); "mouth Ø" =
    the level just below the original surface (depth 0.5 dz).
  ⚠ "axis depth" = the median over columns with centre r < 2.5 mm of the column's
    CONTIGUOUS void run from the top (P2-7's estimator; an overhang's roof keeps a
    column's top solid, so this is the floor depth only where the roof is open).
    "axis depth (any)" = the deepest void cell in those columns.
  ⚠ r17.5/axis = mean of the 16.25/18.75 mm ring medians / the axis (P2-7).
"""
import csv
import math
import os
import re
import shlex
import statistics
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)
ROOT = os.path.abspath(os.path.join(STUD, "..", "..", ".."))
sys.path.insert(0, os.path.join(STUD, "p2_2_props"))
import fabcmp                                        # noqa: E402

BASE_LEG = os.path.join(STUD, "whatif_uptake", "legs", "s100")
H0G = 7272.3958
N_SEC = 3
HER_T = [15.0, 30.0, 60.0, 90.0, 120.0]
HER_DEPTH = [18.0, 31.0, 48.0, 61.0, 77.0]       # mm, Fig 5.10 (log 10-26 §6)
HER_DMAX = [40.0, 47.0, 56.0, 64.0, 73.0]        # mm, +-2
HER_V = {15.0: 16.0, 30.0: 33.0, 60.0: 81.0, 90.0: 135.0}   # cm3 (CSV `time` rows)
HER_MOUTH = 46.0


def legdir(n):
    return os.path.join(HERE, "legs", n)


def done(n):
    return os.path.exists(os.path.join(legdir(n), ".done"))


def geom_of(d):
    kv = {}
    for t in shlex.split(open(os.path.join(d, "argv.sh")).read()):
        if "=" in t:
            k, v = t.split("=", 1)
            kv[k] = v
    nx, ny, nz = (int(v) for v in kv["amr.n_cell"].split())
    lo = [float(v) for v in kv["geometry.prob_lo"].split()]
    hi = [float(v) for v in kv["geometry.prob_hi"].split()]
    return dict(nx=nx, ny=ny, nz=nz, dx=(hi[0] - lo[0]) / nx, dz=(hi[2] - lo[2]) / nz,
                lo=lo, hi=hi, dt=float(kv["timestep"]), tint=int(kv["amr.thermo.plot_int"]),
                st=float(kv["stop_time"]), fl=kv.get("surface_patch.face_list", "0"))


def thermo(d):
    L = open(os.path.join(d, "run", "thermo.dat")).read().split("\n")
    hdr = L[0].split()
    rows = [l.split() for l in L[1:] if l.strip()]
    return hdr, rows


def tcol(hdr, rows, c, skip0=True):
    i = hdr.index(c)
    it = hdr.index("time")
    return [float(r[i]) for r in rows if not (skip0 and float(r[it]) <= 0)]


def events(d):
    with open(os.path.join(d, "run_removal_events.csv")) as f:
        return list(csv.DictReader(f))


def plotfiles(d):
    run = os.path.join(d, "run")
    return sorted(x for x in os.listdir(run) if x.endswith("cell"))


def read_pf(d, pf, g, comps=("removed",)):
    """{comp: ndarray (nx, ny, nz)} for the named components."""
    lvl = os.path.join(d, "run", pf, "Level_0")
    hdrp = os.path.join(d, "run", pf, "Header")
    files = sorted(f for f in os.listdir(lvl) if f.startswith("Cell_D_"))
    out = {c: np.zeros((g["nx"], g["ny"], g["nz"])) for c in comps}
    names = None
    for fn in files:
        blob = open(os.path.join(lvl, fn), "rb").read()
        pos = 0
        while True:
            m = fabcmp.HDR.search(blob, pos)
            if not m:
                break
            lo = [int(m.group(i)) for i in (2, 3, 4)]
            hi = [int(m.group(i)) for i in (5, 6, 7)]
            nc = int(m.group(11))
            if names is None:
                names = fabcmp.var_names(hdrp, nc)
            shp = [hi[q] - lo[q] + 1 for q in range(3)]
            npts = shp[0] * shp[1] * shp[2]
            start = m.end()
            for c in comps:
                ci = names.index(c)
                a = np.frombuffer(blob, dtype="<f8", count=npts, offset=start + ci * npts * 8)
                a = a.reshape(shp[2], shp[1], shp[0]).transpose(2, 1, 0)
                out[c][lo[0]:hi[0] + 1, lo[1]:hi[1] + 1, lo[2]:hi[2] + 1] = a
            pos = start + nc * npts * 8
    return out


def t_of(pf, g):
    return int(pf[:5]) * g["dt"]


def shape(rem, g):
    """the cavity statistics of one void mask."""
    nx, ny, nz, dx, dz = g["nx"], g["ny"], g["nz"], g["dx"], g["dz"]
    void = rem > 0.5
    xc = (np.arange(nx) + 0.5) * dx
    yc = (np.arange(ny) + 0.5) * dx
    X, Y = np.meshgrid(xc, yc, indexing="ij")
    R = np.hypot(X, Y)
    A = np.arctan2(Y, X)
    SEC = np.minimum(N_SEC - 1, (A / (0.5 * math.pi) * N_SEC).astype(int))
    # contiguous void run from the top, per column
    top_solid = ~void[:, :, ::-1]
    run = np.argmax(top_solid, axis=2)
    allvoid = ~top_solid.any(axis=2)
    run = np.where(allvoid, nz, run)
    depth = run * dz
    anyd = np.zeros((nx, ny))
    for i in range(nx):
        for j in range(ny):
            ks = np.nonzero(void[i, j, :])[0]
            anyd[i, j] = (nz - ks.min()) * dz if ks.size else 0.0
    ring = lambda a, b: (R >= a) & (R < b)
    axis = float(np.median(depth[ring(0, 0.0025)]))
    axis_any = float(np.max(anyd[ring(0, 0.0025)]))
    d1625 = float(np.median(depth[ring(0.015, 0.0175)]))
    d1875 = float(np.median(depth[ring(0.0175, 0.020)]))
    r175 = 0.5 * (d1625 + d1875) / axis if axis > 0 else float("nan")

    def level_d(k):
        rs = []
        for s in range(N_SEC):
            m = void[:, :, k] & (SEC == s)
            rs.append(float(R[m].max()) if m.any() else None)
        got = [r for r in rs if r is not None]
        return rs, (2e3 * sum(got) / len(got) if got else 0.0)
    levels = []
    for k in range(nz - 1, -1, -1):
        if not void[:, :, k].any():
            continue
        rs, D = level_d(k)
        levels.append(((nz - k - 0.5) * dz, D, rs))
    dmax = max(levels, key=lambda x: x[1]) if levels else (0.0, 0.0, [None] * N_SEC)
    z_top = g["hi"][2]

    def at_depth(dmm):
        k = int((z_top - dmm * 1e-3 - g["lo"][2]) / dz + 1e-9)
        if not (0 <= k < nz):
            return None, None
        rs, D = level_d(k)
        got = [2e3 * r for r in rs if r is not None]
        return [None if r is None else 2e3 * r for r in rs], \
               (max(got) - min(got) if len(got) == N_SEC else None)
    mouth = level_d(nz - 1)[1]
    over = int((void[:, :, :-1] & ~void[:, :, 1:]).sum())
    return dict(axis=axis, axis_any=axis_any, r175=r175, dmax=dmax[1], dmax_depth=dmax[0],
                dmax_sec=dmax[2], mouth=mouth, at=at_depth, overhang=over,
                vol_mask=float(void.sum()) * dx * dx * dz * 4.0 * 1e6, levels=levels)


def plate_census(rem, g):
    """F3's fallback, rebuilt INDEPENDENTLY from a void mask with the code's rule
    (solid indicator; above the domain top = gas, elsewhere outside = nearest
    cell; central difference; zero gradient => w = 1), split into PLATES (rev 3:
    a fallback cell whose +-z neighbours are BOTH non-solid) and the rest.
    unit/face_list (a4) proved this rule equal to the code's census cell for cell."""
    S = (rem <= 0.5).astype(float)
    P = np.pad(S, 1, mode="edge")
    P[:, :, -1] = 0.0
    gx = P[2:, 1:-1, 1:-1] - P[:-2, 1:-1, 1:-1]
    gy = P[1:-1, 2:, 1:-1] - P[1:-1, :-2, 1:-1]
    gz = P[1:-1, 1:-1, 2:] - P[1:-1, 1:-1, :-2]
    solid = S > 0.5
    V = ~solid
    op = np.zeros((6,) + S.shape, bool)
    op[0, 1:] = V[:-1]; op[1, :-1] = V[1:]
    op[2, :, 1:] = V[:, :-1]; op[3, :, :-1] = V[:, 1:]
    op[4, :, :, 1:] = V[:, :, :-1]
    op[5, :, :, :-1] = V[:, :, 1:]; op[5, :, :, -1] = True
    op &= solid[None]
    nface = op.sum(0)
    zero = (gx == 0) & (gy == 0) & (gz == 0) & (nface > 0)
    # +-z both non-solid: below void (or nothing -- the base is solid by the mirror),
    # above void or the open top
    below = np.zeros_like(solid); below[:, :, 1:] = V[:, :, :-1]
    above = np.ones_like(solid);  above[:, :, :-1] = V[:, :, 1:]
    plate = zero & below & above
    nx = rem.shape[0]
    xc = (np.arange(nx) + 0.5) * g["dx"]
    X, Y = np.meshgrid(xc, xc, indexing="ij")
    R = np.hypot(X, Y)
    ks = np.nonzero(plate)
    rr = R[ks[0], ks[1]] * 1e3
    return dict(faces=int(nface.sum()), fb_cells=int(zero.sum()), fb_faces=int(nface[zero].sum()),
                pl_cells=int(plate.sum()), pl_faces=int(nface[plate].sum()),
                np_faces=int(nface[zero & ~plate].sum()),
                r=(float(rr.min()), float(rr.max())) if rr.size else None,
                k=(int(ks[2].min()), int(ks[2].max())) if rr.size else None)


def leg_summary(n, which=-1):
    d = legdir(n)
    g = geom_of(d)
    pfs = plotfiles(d)
    pf = pfs[which]
    t = t_of(pf, g)
    rem = read_pf(d, pf, g)["removed"]
    sh = shape(rem, g)
    ev = events(d)
    vc = g["dx"] * g["dx"] * g["dz"]
    # ⚠ the event log stamps a step's START (SR 7 iii): events before t
    nv = sum(int(float(e["n_voided"])) for e in ev if float(e["time"]) < t - 0.5 * g["dt"])
    nv6 = sum(int(float(e["n_voided"])) for e in ev
              if float(e["time"]) < t - 0.5 * g["dt"] and int(float(e["regime"])) == 6)
    hdr, rows = thermo(d)
    it = hdr.index("time")
    sel = [r for r in rows if 0 < float(r[it]) <= t + 1e-9]
    last = min(rows, key=lambda r: abs(float(r[it]) - t))
    row = {c: float(last[hdr.index(c)]) for c in hdr}
    led = max(abs(float(r[hdr.index("ledger_err_solid")])) for r in sel) if sel else float("nan")
    idm = max(float(r[hdr.index("fl_id_march")]) for r in sel) if "fl_id_march" in hdr and sel else None
    idl = max(float(r[hdr.index("fl_id_ledger")]) for r in sel) if "fl_id_ledger" in hdr and sel else None
    return dict(leg=n, g=g, t=t, pf=pf, sh=sh, V=nv * vc * 4e6, V6=nv6 * vc * 4e6,
                row=row, led=led, idm=idm, idl=idl, E_src=row["ledger_E_src"] * 4.0)


def fmt_sec(rs):
    return " ".join("  -  " if r is None else f"{r:5.1f}" for r in rs)


# ---------------------------------------------------------------- GATE ----
def md5(p):
    import hashlib
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def gate():
    out = []
    P = out.append
    ok_all = True
    P("P2-8 (b1, b2) -- the identity gate and the face-list witnesses")
    P("=" * 96)
    want = {"run_removal_events.csv": "58cb6ab598461b0af7e6fb5cbb60bec8",
            os.path.join("run", "thermo.dat"): "9fc066b02910cee89b03d7a0c557b0c1"}
    b1 = done("K0")
    for rel, w in want.items():
        got = md5(os.path.join(legdir("K0"), rel)) if done("K0") else "(no K0)"
        ref = md5(os.path.join(BASE_LEG, rel))
        same = got == w == ref
        b1 &= same
        P(f"  b1  {rel:24s} K0 {got}  packet {w}  s100 {ref}  {'IDENTICAL' if same else 'DIFFER'}")
    P(f"  b1  {'PASS' if b1 else 'FAIL'} -- the lid edit and the group-A deletion leave Hoeser"
      f" byte-identical with the key off" if b1 else "  b1  FAIL")
    ok_all &= b1
    if not (done("L0") and done("K0")):
        P("  b2  L0 or K0 not run")
        print("\n".join(out))
        return 1
    hk, rk = thermo(legdir("K0"))
    hl, rl = thermo(legdir("L0"))
    tk = [r for r in rk if float(r[hk.index("time")]) > 0][0]
    tl = [r for r in rl if float(r[hl.index("time")]) > 0][0]
    pk, pl = tk[hk.index("patch_P_robin")], tl[hl.index("patch_P_robin")]
    c1 = pk == pl
    P(f"  b2  step-1 patch_P_robin (printed token, t = {tl[hl.index('time')]}): K0 {pk}  L0 {pl}"
      f"  {'EQUAL' if c1 else 'DIFFER'}")
    body = [r for r in rl if float(r[hl.index("time")]) > 0]
    col = lambda c: [float(r[hl.index(c)]) for r in body]
    idm, idl = max(col("fl_id_march")), max(col("fl_id_ledger"))
    led = max(abs(x) for x in col("ledger_err_solid"))
    ltop = max(col("fl_lat_top"))
    mism = max(col("fl_h_top_mismatch"))
    nf = [a + b + c + d for a, b, c, d in zip(col("fl_n_top"), col("fl_n_up"),
                                              col("fl_n_side"), col("fl_n_down"))]
    fbr = [f / n if n > 0 else 0.0 for f, n in zip(col("fl_w_fallback"), nf)]
    i_fb = int(np.argmax(fbr))
    c2 = idm <= 1e-9 and idl <= 1e-9
    c3 = led <= 5e-3
    c4 = ltop == 0.0
    # rev 3: the bar EXCLUDES one-cell plates. The plate split needs the mask, and
    # L0 writes one mid-run mask (30 s), so the split is scored THERE; the
    # thermo-only total is printed every row for the record.
    g0 = geom_of(legdir("L0"))
    pfl = plotfiles(legdir("L0"))[-1]
    pc = plate_census(read_pf(legdir("L0"), pfl, g0)["removed"], g0)
    npf = pc["np_faces"] / pc["faces"] if pc["faces"] else 0.0
    c5 = npf <= 0.01
    c6 = mism == 0.0
    P(f"  b2  power identity, every row: max |jet_P_face - sum fl_P|/sum = {idm:.3e};"
      f" max |ledger_P_robin - sum fl_P|/sum = {idl:.3e}  (bar <= 1e-9)  {'ok' if c2 else 'MISS'}")
    P(f"  b2  ledger_err_solid max |.| {led:.3e} (bar <= 5e-3)  {'ok' if c3 else 'MISS'}")
    P(f"  b2  fl_lat_top max {ltop:.0f} over {len(body)} rows (bar 0 every row)  {'ok' if c4 else 'MISS'}")
    P(f"  b2  fl_w_fallback / faces (thermo, ALL fallback incl. plates): max {max(fbr)*100:.3f} % at"
      f" t = {float(body[i_fb][hl.index('time')]):g} s ({col('fl_w_fallback')[i_fb]:.0f} of {nf[i_fb]:.0f})")
    P(f"  b2  at the {t_of(pfl, g0):g} s mask: {pc['faces']} faces; fallback {pc['fb_cells']} cells /"
      f" {pc['fb_faces']} faces; PLATES {pc['pl_cells']} cells / {pc['pl_faces']} faces"
      f" ({pc['pl_faces']/pc['faces']*100:.2f} %), r {pc['r'][0]:.1f}-{pc['r'][1]:.1f} mm, k {pc['k'][0]}-{pc['k'][1]}"
      if pc["r"] else f"  b2  at the mask: no plates")
    P(f"  b2  NON-PLATE fallback {pc['np_faces']} faces = {npf*100:.3f} % (bar <= 1 %)  {'ok' if c5 else 'MISS'}")
    P(f"  b2  fl_h_top_mismatch max {mism:.0f} (the top face reproduces flame_h)  {'ok' if c6 else 'MISS'}")
    b2 = c1 and c2 and c3 and c4 and c5 and c6
    P(f"  b2  {'PASS' if b2 else 'FAIL'}")
    ok_all &= b2
    print("\n".join(out))
    open(os.path.join(HERE, "GATE_b.txt"), "w").write("\n".join(out) + "\n")
    return 0 if ok_all else 1


# ---------------------------------------------------------------- MESH ----
def mesh():
    out = []
    P = out.append
    P("P2-8 (b3, c1, c2) -- mesh twin, the key's effect at 30 s, isotropy")
    P("=" * 96)
    S = {n: leg_summary(n) for n in ("K0", "L0", "L0c") if done(n)}
    for n, s in S.items():
        sh = s["sh"]
        P(f"  {n:4s} t {s['t']:.1f} s  V(log) {s['V']:7.3f}  V(mask) {sh['vol_mask']:7.3f} cm3 x4"
          f"  lateral {s['V6']/max(s['V'],1e-30)*100:5.1f} %  axis {sh['axis']*1e3:5.1f} mm"
          f" (any {sh['axis_any']*1e3:5.1f})  r17.5/axis {sh['r175']:.3f}  max Ø {sh['dmax']:5.1f} mm"
          f" at {sh['dmax_depth']*1e3:4.1f} mm  mouth Ø {sh['mouth']:5.1f}  overhangs {sh['overhang']}")
    if "L0" in S and "L0c" in S:
        a, b = S["L0"], S["L0c"]
        dv = abs(b["V"] - a["V"]) / a["V"]
        dd = abs(b["sh"]["dmax"] - a["sh"]["dmax"])
        ok = dv <= 0.15 and dd <= 4.0
        P(f"  b3  L0c / L0 volume {b['V']:.3f} / {a['V']:.3f} = {dv*100:.1f} % apart (bar 15 %);"
          f" max Ø {b['sh']['dmax']:.1f} vs {a['sh']['dmax']:.1f} mm, {dd:.1f} mm apart (bar 4)"
          f"  {'PASS' if ok else 'FAIL'}")
    if "L0" in S and "K0" in S:
        a, k = S["L0"], S["K0"]
        P(f"  c1  max Ø L0 {a['sh']['dmax']:.1f} vs K0 {k['sh']['dmax']:.1f} mm (P1 42-50, K0 <= 41)")
        P(f"  c1  volume L0/K0 {a['V']/k['V']:.3f} (P1 0.95-1.50); lateral share L0"
          f" {a['V6']/a['V']*100:.1f} % (P1 8-35); K0 {k['V6']/k['V']*100:.1f} %")
        P(f"  c1  axis depth ratio L0/K0 {a['sh']['axis']/k['sh']['axis']:.3f} (P1 0.90-1.10);"
          f" r17.5/axis L0 {a['sh']['r175']:.3f} vs K0 {k['sh']['r175']:.3f} (P1 L0 0.50-0.64, BELOW K0)")
    if "L0" in S:
        for dmm in (5, 10):
            rs, sp = S["L0"]["sh"]["at"](dmm)
            P(f"  c2  L0 Ø by 30° sector at {dmm} mm: {fmt_sec(rs)} mm  spread "
              f"{'n/a' if sp is None else f'{sp:.2f}'} (bar <= 3)"
              f"  {'' if sp is None else ('PASS' if sp <= 3.0 else 'MISS')}")
    print("\n".join(out))
    open(os.path.join(HERE, "MESH.txt"), "w").write("\n".join(out) + "\n")


def c3():
    out = []
    P = out.append
    P("P2-8 (c3) -- L0f (0.5 mm) vs L0 (1 mm) at 30 s")
    P("=" * 96)
    a, f = leg_summary("L0"), leg_summary("L0f")
    for s in (a, f):
        sh = s["sh"]
        P(f"  {s['leg']:4s} V {s['V']:7.3f} cm3  max Ø {sh['dmax']:5.1f} mm  axis {sh['axis']*1e3:5.2f} mm"
          f"  mouth {sh['mouth']:5.1f}  lateral {s['V6']/s['V']*100:5.1f} %  overhangs {sh['overhang']}")
    dv = abs(f["V"] - a["V"]) / a["V"]
    dd = abs(f["sh"]["dmax"] - a["sh"]["dmax"])
    da = abs(f["sh"]["axis"] - a["sh"]["axis"]) * 1e3
    P(f"  c3  volume {dv*100:.2f} % (bar 5) | max Ø {dd:.2f} mm (bar 2) | axis depth {da:.2f} mm (bar 1.5)"
      f"  {'PASS' if (dv <= 0.05 and dd <= 2.0 and da <= 1.5) else 'MISS'}")
    print("\n".join(out))
    open(os.path.join(HERE, "C3.txt"), "w").write("\n".join(out) + "\n")


# ---------------------------------------------------------------- SSTAR ----
def sstar():
    out = []
    P = out.append
    P("P2-8 (d1) -- S* by the packet's rule: log V(30 s, x4) vs log S, bracketing pair, -> 33.0 cm3")
    P("=" * 96)
    pts = []
    for n, S in (("C15", 0.15), ("C20", 0.20), ("C25", 0.25), ("C10", 0.10), ("C30", 0.30)):
        src = n
        if n == "C10" and not done("C10") and done("L0"):
            # S 0.10 x H0G = 727.239580 = BASE_C: the C10 leg IS L0 (same key, same h_expr,
            # same binary); the extension point is read from L0, checked token for token
            import run_p2_8 as R
            a = R.argv_of("L0", "/X", probe=True)[0]
            b = shlex.split(open(os.path.join(legdir("L0"), "argv.sh")).read())
            b = ["plot_file=/X/run" if t.startswith("plot_file=") else t for t in b]
            assert R.c_of(0.10) == R.BASE_C and a == b, "L0 argv is not the C10 argv"
            src = "L0"
            P("  C10  = L0 (S 0.10 x H0G = BASE_C; argv token-for-token, binary f4c27b27...) -- not rerun")
        if done(src):
            s = leg_summary(src)
            pts.append((S, s["V"], n))
            P(f"  {n}  S {S:.2f}  V(30 s) {s['V']:.3f} cm3 x4")
    pts.sort()
    target = 33.0
    star = None
    for (s1, v1, n1), (s2, v2, n2) in zip(pts[:-1], pts[1:]):
        if (v1 - target) * (v2 - target) <= 0 and v1 > 0 and v2 > 0:
            x = math.log(s1) + (math.log(target) - math.log(v1)) * \
                (math.log(s2) - math.log(s1)) / (math.log(v2) - math.log(v1))
            star = round(math.exp(x), 4)
            P(f"  bracket {n1}-{n2}: S* = exp(interp) = {math.exp(x):.6f} -> {star:.4f}")
            break
    if star is None:
        P("  NOT BRACKETED -- the rule adds ONE point (S 0.10 or 0.30) and repeats; still outside => stop")
        print("\n".join(out))
        open(os.path.join(HERE, "SSTAR.txt"), "w").write("\n".join(out) + "\n")
        return 1
    P(f"  C* = S* x {H0G} = {star*H0G:.6f}   (P3 predicted S* 0.16-0.24)")
    P(f"S_STAR={star:.4f}")
    print("\n".join(out))
    open(os.path.join(HERE, "SSTAR.txt"), "w").write("\n".join(out) + "\n")
    return 0


# ---------------------------------------------------------------- RESULTS ----
def results():
    P = print
    P("P2-8 -- the face list on Hoeser: per leg and per plotfile")
    P("=" * 110)
    for n in ("K0", "L0", "L0c", "L0f", "C15", "C20", "C25", "P47", "P47off"):
        if not done(n):
            continue
        d = legdir(n)
        g = geom_of(d)
        P(f"\n  LEG {n}   n_cell {g['nx']}x{g['ny']}x{g['nz']}  dx {g['dx']*1e3:.2f} mm  face_list {g['fl']}"
          f"  stop {g['st']:g} s")
        P(f"    {'t[s]':>6s} {'V log':>8s} {'V mask':>8s} {'lat%':>5s} {'axis':>6s} {'any':>6s}"
          f" {'r17.5/ax':>8s} {'maxØ':>6s} {'@depth':>6s} {'mouth':>6s} {'over':>5s}"
          f" | {'n_top':>6s} {'n_up':>5s} {'n_side':>6s} {'n_down':>6s} {'fb%':>6s}"
          f" | {'led_solid':>9s} {'id_m':>8s} {'id_l':>8s}")
        for q in range(1, len(plotfiles(d))):
            s = leg_summary(n, q)
            sh, r = s["sh"], s["row"]
            nf = sum(r.get(f"fl_n_{c}", 0.0) for c in ("top", "up", "side", "down"))
            fb = (r.get("fl_w_fallback", 0.0) / nf * 100) if nf else float("nan")
            P(f"    {s['t']:6.1f} {s['V']:8.3f} {sh['vol_mask']:8.3f} {s['V6']/max(s['V'],1e-30)*100:5.1f}"
              f" {sh['axis']*1e3:6.1f} {sh['axis_any']*1e3:6.1f} {sh['r175']:8.3f} {sh['dmax']:6.1f}"
              f" {sh['dmax_depth']*1e3:6.1f} {sh['mouth']:6.1f} {sh['overhang']:5d}"
              f" | {r.get('fl_n_top',0):6.0f} {r.get('fl_n_up',0):5.0f} {r.get('fl_n_side',0):6.0f}"
              f" {r.get('fl_n_down',0):6.0f} {fb:6.2f}"
              f" | {s['led']:9.2e} {s['idm'] if s['idm'] is not None else float('nan'):8.1e}"
              f" {s['idl'] if s['idl'] is not None else float('nan'):8.1e}")
            for dmm in (5, 10, 30, 50):
                rs, sp = sh["at"](dmm)
                if rs is not None and any(x is not None for x in rs):
                    P(f"           Ø by sector at {dmm:2d} mm: {fmt_sec(rs)}  spread "
                      f"{'n/a' if sp is None else f'{sp:.2f}'}")
    # d2 / d3 / d4
    if done("P47"):
        d = legdir("P47")
        g = geom_of(d)
        P("\n" + "=" * 110)
        P("  d2  P47 vs Hoeser 47 MPa (Fig. 5.10; volumes from the CSV `time` rows)")
        P(f"    {'t':>5s} | {'V model':>8s} {'V her':>6s} | {'axis':>6s} {'her':>5s} {'%':>6s} |"
          f" {'maxØ':>6s} {'her':>5s} {'Δ':>6s} | {'mouth':>6s}")
        rows = {}
        for q in range(1, len(plotfiles(d))):
            s = leg_summary("P47", q)
            rows[round(s["t"], 3)] = s
        agree = True
        for t, dh, Dh in zip(HER_T, HER_DEPTH, HER_DMAX):
            s = rows.get(t)
            if s is None:
                continue
            sh = s["sh"]
            ax = sh["axis_any"] * 1e3
            pa = (ax - dh) / dh * 100
            dD = sh["dmax"] - Dh
            agree &= abs(pa) <= 15.0 and (t < 60 or abs(dD) <= 6.0)
            P(f"    {t:5.0f} | {s['V']:8.2f} {HER_V.get(t, float('nan')):6.1f} | {ax:6.1f} {dh:5.0f}"
              f" {pa:+6.1f} | {sh['dmax']:6.1f} {Dh:5.0f} {dD:+6.1f} | {sh['mouth']:6.1f}")
        if 30.0 in rows:
            v30 = rows[30.0]["V"]
            for t, her_r in ((60.0, 81.0 / 33.0), (90.0, 135.0 / 33.0)):
                if t in rows:
                    P(f"    V{int(t)}/V30 model {rows[t]['V']/v30:.2f}  hers {her_r:.2f}")
        if 15.0 in rows and 120.0 in rows:
            a, b = rows[15.0]["sh"], rows[120.0]["sh"]
            dD = b["dmax"] - a["dmax"]
            gr = (b["axis_any"] - a["axis_any"]) * 1e3 / dD if dD else float("inf")
            P(f"    depth growth / diameter growth 15->120 s: model {gr:.2f}  hers 1.79")
        P(f"    SHAPE AGREES (max Ø within 6 mm at 60/90/120 s AND axis within 15 % at every time):"
          f" {'YES' if agree else 'NO'}")
        P("    ⚠ axis depth here is the deepest void on the axis (`any`): her axis depth is the cavity's.")
    if done("P47off"):
        d = legdir("P47off")
        P("\n  d3  P47off max Ø at every plotfile (bar <= 41 mm)")
        mx = []
        for q in range(1, len(plotfiles(d))):
            s = leg_summary("P47off", q)
            mx.append((s["t"], s["sh"]["dmax"]))
        P("    " + "  ".join(f"{t:.0f}s {D:.1f}" for t, D in mx)
          + f"   -> {'PASS' if all(D <= 41.0 for _, D in mx) else 'MISS'}")
        if done("P47"):
            s47 = leg_summary("P47")
            P(f"    P47 - P47off max Ø at 120 s: {s47['sh']['dmax'] - mx[-1][1]:.1f} mm (P5 >= 10)")
    if done("P47"):
        s = leg_summary("P47")
        sh, r = s["sh"], s["row"]
        P("\n  d4  P47 at 120 s")
        P(f"    ledger_err_solid max {s['led']:.3e} (bar 5e-3); identity march {s['idm']:.2e}"
          f" ledger {s['idl']:.2e} (bar 1e-9); overhangs {sh['overhang']}")
        for dmm in (30, 50):
            rs, sp = sh["at"](dmm)
            P(f"    Ø by sector at {dmm} mm: {fmt_sec(rs) if rs else '-'}  spread "
              f"{'n/a' if sp is None else f'{sp:.2f}'} (bar 3)")
        d47 = legdir("P47"); g47 = geom_of(d47)
        for pf in plotfiles(d47)[1:]:
            pc = plate_census(read_pf(d47, pf, g47)["removed"], g47)
            P(f"    {t_of(pf, g47):5.0f} s plates {pc['pl_cells']} cells / {pc['pl_faces']} faces"
              f" ({pc['pl_faces']/pc['faces']*100:.2f} % of {pc['faces']})"
              + (f" r {pc['r'][0]:.1f}-{pc['r'][1]:.1f} mm k {pc['k'][0]}-{pc['k'][1]}" if pc["r"] else "")
              + f"; non-plate fallback {pc['np_faces']} ({pc['np_faces']/pc['faces']*100:.3f} %)")
        P(f"    direction census: top {r['fl_n_top']:.0f} up {r['fl_n_up']:.0f} side {r['fl_n_side']:.0f}"
          f" down {r['fl_n_down']:.0f}; fallback {r['fl_w_fallback']:.0f}")
        hdr, rows = thermo(legdir("P47"))
        for c in ("lat_shift_refuse", "lat_orphan", "lat_orphan_n"):
            if c in hdr:
                P(f"    {c}: sampled max {max(tcol(hdr, rows, c)):.0f}, sampled sum"
                  f" {sum(tcol(hdr, rows, c)):.0f}  (per-step counter, printed every {geom_of(legdir('P47'))['tint']} steps)")


def main():
    m = sys.argv[1] if len(sys.argv) > 1 else ""
    if m == "GATE":
        sys.exit(gate())
    if m == "MESH":
        return mesh()
    if m == "C3":
        return c3()
    if m == "SSTAR":
        sys.exit(sstar())
    results()


if __name__ == "__main__":
    main()
