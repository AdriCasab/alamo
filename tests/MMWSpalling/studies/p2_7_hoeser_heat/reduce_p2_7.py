#!/usr/bin/env python3
"""P2-7 deliverable 3 -- the reductions (criteria b1, c1-c8, witnesses).

  reduce_p2_7.py GATE     b1: A0 vs p2_5_hoeser/legs/A, byte for byte
                          (run_removal_events.csv and run/thermo.dat); also
                          writes GATE_b1.txt
  reduce_p2_7.py          every leg, then the criteria table  (-> RESULTS.txt)

⚠ GEOMETRY IS PARSED FROM EACH LEG'S OWN argv.sh (the LZ trap, SR 5): A0 is a
  100 mm block and HA/HG/HGd are 150 mm, so no module constant is safe.
⚠ Quarter domain (SR 3): volume AND energy are x4; a ratio of two quadrant
  quantities (flake share, price) is formed from x4 on BOTH sides.
⚠ `ledger_E_src` and `ledger_E_flake` are CUMULATIVE: read at their FINAL
  value, never summed (SR 7 i; amr.thermo.plot_int = 6). `ledger_err_solid` is
  RELATIVE (divided by ledger_E_src, SR 7 ii).
⚠ The depth profile is read from the VOID MASK at 30 s (the state). The event
  log is a cross-check through SR 4: its k_top is the top BEFORE the removal,
  so the post-removal top is k_top - n_voided.
⚠ The bore is reported PER 30-degree SECTOR (SR 10), as the max void-cell
  centre radius in the sector (bore_profile.sector_radii's rule, re-implemented
  here so its module-level constants are not imported).
"""
import csv
import hashlib
import math
import os
import shlex
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)
ROOT = os.path.abspath(os.path.join(STUD, "..", "..", ".."))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(STUD, "p2_2_props"))
import fabcmp                                   # noqa: E402
import energy_bound as EB                       # noqa: E402

A_P25 = os.path.join(STUD, "p2_5_hoeser", "legs", "A")
LEGS = ("A0", "HA", "HG", "HGd")
R_FLOOR = 1
RING = 0.0025
N_SEC = 3                    # 30-degree sectors across the quadrant
HER_V0 = 12.9                # cm3, 0 MPa 30 s (Fig 5.12)
HER_DR = None                # filled from energy_bound.bowl()


def legdir(n):
    return A_P25 if n == "A" else os.path.join(HERE, "legs", n)


def geom_of(d):
    toks = shlex.split(open(os.path.join(d, "argv.sh")).read())
    kv = {}
    for t in toks:
        if "=" in t:
            k, v = t.split("=", 1)
            kv[k] = v
    nx, ny, nz = (int(v) for v in kv["amr.n_cell"].split())
    lo = [float(v) for v in kv["geometry.prob_lo"].split()]
    hi = [float(v) for v in kv["geometry.prob_hi"].split()]
    dx, dz = (hi[0] - lo[0]) / nx, (hi[2] - lo[2]) / nz
    return dict(nx=nx, ny=ny, nz=nz, dx=dx, dz=dz, lo=lo, hi=hi,
                V_cell=dx * dx * dz, dt=float(kv["timestep"]),
                tint=int(kv["amr.thermo.plot_int"]),
                pint=int(kv["amr.plot_int"]), st=float(kv["stop_time"]),
                cref=kv.get("spallation.sp.constraint_depth"))


def thermo(d):
    L = open(os.path.join(d, "run", "thermo.dat")).read().split("\n")
    hdr = L[0].split()
    rows = [l.split() for l in L[1:] if l.strip()]
    out = {n: [] for n in hdr}
    for r in rows:
        if float(r[hdr.index("time")]) <= 0:
            continue
        for n in hdr:
            out[n].append(float(r[hdr.index(n)]))
    return out


def events(d):
    with open(os.path.join(d, "run_removal_events.csv")) as f:
        return list(csv.DictReader(f))


def _f(r, k):
    try:
        return float(r.get(k, "") or 0.0)
    except ValueError:
        return 0.0


def _i(r, k):
    try:
        return int(float(r.get(k, "") or 0))
    except ValueError:
        return 0


def plotfiles(d):
    run = os.path.join(d, "run")
    return sorted(x for x in os.listdir(run) if x.endswith("cell"))


def read_pf(d, pf):
    """(removed, Temp) dicts keyed (i, j, k) from one cell plotfile."""
    lvl = os.path.join(d, "run", pf, "Level_0")
    hdr = os.path.join(d, "run", pf, "Header")
    files = sorted(f for f in os.listdir(lvl) if f.startswith("Cell_D_"))
    first = fabcmp.read_fabs_boxed(os.path.join(lvl, files[0]))
    names = fabcmp.var_names(hdr, first[0]["ncomp"])
    ir, it = names.index("removed"), names.index("Temp")
    iT = names.index("T_surf_field") if "T_surf_field" in names else None
    rem, tmp, tsf = {}, {}, {}
    for fn in files:
        for fab in fabcmp.read_fabs_boxed(os.path.join(lvl, fn)):
            lo, hi = fab["lo"], fab["hi"]
            nxb, nyb, nzb = (hi[q] - lo[q] + 1 for q in range(3))
            cr, ct = fab["comp"][ir], fab["comp"][it]
            cs = fab["comp"][iT] if iT is not None else None
            for kk in range(nzb):
                for jj in range(nyb):
                    base = nxb * (jj + nyb * kk)
                    for ii in range(nxb):
                        o = ii + base
                        key = (lo[0] + ii, lo[1] + jj, lo[2] + kk)
                        rem[key] = cr[o]
                        tmp[key] = ct[o]
                        if cs is not None:
                            tsf[key] = cs[o]
    return rem, tmp, tsf


def col_depths(rem, g):
    """Per column: contiguous void cells from the top x dz (the recession)."""
    dep = {}
    for j in range(g["ny"]):
        for i in range(g["nx"]):
            n = 0
            for k in range(g["nz"] - 1, -1, -1):
                if rem.get((i, j, k), 0.0) > 0.5:
                    n += 1
                else:
                    break
            dep[(i, j)] = n * g["dz"]
    return dep


def rc(i, j, g):
    return g["dx"] * math.hypot(i + 0.5, j + 0.5)


def face_mean(rem, tmp, g, rmax=0.010):
    """P2-5 c12's estimator: mean T of the top solid cell over r < rmax.
    ⚠ This is the top cell's CENTRE temperature, dz/2 below the face."""
    vals = []
    for j in range(g["ny"]):
        for i in range(g["nx"]):
            if rc(i, j, g) >= rmax:
                continue
            for k in range(g["nz"] - 1, -1, -1):
                if rem.get((i, j, k), 0.0) <= 0.5:
                    vals.append(tmp[(i, j, k)])
                    break
    return (sum(vals) / len(vals), len(vals)) if vals else (float("nan"), 0)


def face_recon(rem, tsf, g, rmax=0.010):
    """The model's OWN face T: `T_surf_field` (MMWSpalling.H:7921, the P0
    reconstruction the criterion fires on) on the top solid cell over r < rmax.
    Written per step on top faces only; 0 where the reconstruction fell back,
    so only positive entries are averaged and the count is printed."""
    vals, nfb = [], 0
    for j in range(g["ny"]):
        for i in range(g["nx"]):
            if rc(i, j, g) >= rmax:
                continue
            for k in range(g["nz"] - 1, -1, -1):
                if rem.get((i, j, k), 0.0) <= 0.5:
                    v = tsf.get((i, j, k), 0.0)
                    if v > 0.0:
                        vals.append(v)
                    else:
                        nfb += 1
                    break
    return (sum(vals) / len(vals) if vals else float("nan")), len(vals), nfb


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def gate():
    out = []
    P = out.append
    P("P2-7 (b1) -- A0 vs p2_5_hoeser/legs/A, byte for byte")
    P("=" * 92)
    a0 = legdir("A0")
    if not os.path.exists(os.path.join(a0, ".done")):
        P("  A0 has no .done sentinel -- not run or failed")
        print("\n".join(out))
        return 1
    ok = True
    for rel in ("run_removal_events.csv", os.path.join("run", "thermo.dat")):
        p1, p2 = os.path.join(A_P25, rel), os.path.join(a0, rel)
        m1, m2 = md5(p1), md5(p2)
        same = open(p1, "rb").read() == open(p2, "rb").read()
        ok &= same
        P(f"  {rel:<24s} A  {m1}  ({os.path.getsize(p1)} B)")
        P(f"  {'':<24s} A0 {m2}  ({os.path.getsize(p2)} B)   "
          f"{'IDENTICAL' if same else 'DIFFER'}")
    P(f"  b1 {'PASS' if ok else 'FAIL'} -- "
      + ("the binary has not moved Hoeser since P2-5; A is a valid baseline"
         if ok else "the binary MOVED Hoeser since P2-5: re-baseline before ord 3"))
    print("\n".join(out))
    open(os.path.join(HERE, "GATE_b1.txt"), "w").write("\n".join(out) + "\n")
    return 0 if ok else 1


def one(n):
    d = legdir(n)
    g = geom_of(d)
    th = thermo(d)
    ev = events(d)
    pfs = plotfiles(d)
    H = g["hi"][2] - g["lo"][2]
    r = dict(leg=n, g=g, H=H)
    # ---- volume and energy --------------------------------------------
    nv = sum(_i(e, "n_voided") for e in ev)
    r["V"] = nv * g["V_cell"] * 4.0 * 1e6
    r["E_src"] = th["ledger_E_src"][-1] * 4.0
    r["E_fl"] = th["ledger_E_flake"][-1] * 4.0
    r["t_end"] = th["time"][-1]
    r["share"] = r["E_fl"] / r["E_src"]
    r["price"] = r["E_fl"] / (r["V"] * 1e-6) / 1e9
    r["err"] = max(abs(x) for x in th["ledger_err_solid"])
    # ---- the mask at the end -------------------------------------------
    rem, tmp, tsf = read_pf(d, pfs[-1])
    r["pf_end"] = pfs[-1]
    r["t_pf_end"] = int(pfs[-1][:5]) * g["dt"]
    dep = col_depths(rem, g)
    r["V_mask"] = sum(dep.values()) * g["dx"] ** 2 * 4.0 * 1e6
    # event-log cross-check of the per-column top (SR 4)
    top0 = g["nz"] - 1
    post = {}
    for e in ev:
        if _i(e, "regime") != R_FLOOR:
            continue
        key = (_i(e, "col_i"), _i(e, "col_j"))
        p = _i(e, "k_top") - _i(e, "n_voided")
        post[key] = min(post.get(key, top0), p)
    agree = sum(1 for key, p in post.items()
                if abs((top0 - p) * g["dz"] - dep[key]) < 0.5 * g["dz"])
    r["ev_cols"], r["ev_agree"] = len(post), agree
    # rings
    rings = {}
    for (i, j), dv in dep.items():
        b = int(rc(i, j, g) / RING)
        rings.setdefault(b, []).append(dv)
    r["rings"] = {b: (max(v), statistics.median(v), len(v))
                  for b, v in sorted(rings.items()) if b * RING < 0.0301}
    ax = r["rings"][0][1]
    d1625, d1875 = r["rings"][6][1], r["rings"][7][1]
    r["axis"] = ax
    r["axis00"] = dep[(0, 0)]
    r["d175"] = 0.5 * (d1625 + d1875)
    cen = [dv for (i, j), dv in dep.items() if 0.01625 <= rc(i, j, g) < 0.01875]
    r["d175c"] = statistics.median(cen)
    r["dr"] = r["d175"] / ax if ax > 0 else float("nan")
    r["drc"] = r["d175c"] / ax if ax > 0 else float("nan")
    r["deep_col"] = max(dep.values())
    kmin = min((k for (i, j, k), v in rem.items() if v > 0.5), default=None)
    r["deep_any"] = (g["hi"][2] - (g["lo"][2] + kmin * g["dz"])) if kmin is not None else 0.0
    # bore per sector at 10 mm (bore_profile's level rule)
    z = g["hi"][2] - 0.010
    k10 = int((z - g["lo"][2]) / g["dz"] + 1e-9)
    sec = [None] * N_SEC
    for j in range(g["ny"]):
        for i in range(g["nx"]):
            if rem.get((i, j, k10), 0.0) <= 0.5:
                continue
            x, y = (i + 0.5) * g["dx"], (j + 0.5) * g["dx"]
            s = min(N_SEC - 1, int(math.atan2(y, x) / (0.5 * math.pi) * N_SEC))
            rr = math.hypot(x, y)
            sec[s] = rr if sec[s] is None else max(sec[s], rr)
    r["k10"] = k10
    r["sec"] = [2e3 * x if x is not None else None for x in sec]
    # ---- criterion and face --------------------------------------------
    fl = [_f(e, "T_surf") for e in ev
          if _i(e, "regime") == R_FLOOR and _f(e, "T_surf") > 0]
    r["fire_med"] = statistics.median(fl) if fl else float("nan")
    r["n_fire"] = len(fl)
    faces = []
    for pf in pfs:
        t = int(pf[:5]) * g["dt"]
        if t <= 0:
            continue
        rm_, tp_, ts_ = (rem, tmp, tsf) if pf == pfs[-1] else read_pf(d, pf)
        m, cnt = face_mean(rm_, tp_, g)
        mr, nr, nfb = face_recon(rm_, ts_, g)
        faces.append((t, m, cnt, mr, nr, nfb))
    r["faces"] = faces
    r["face_tmean"] = sum(f[1] for f in faces) / len(faces)
    r["recon_tmean"] = sum(f[3] for f in faces) / len(faces)
    # ---- constraint-depth witnesses ------------------------------------
    if "cref_d" in th:
        r["cref_d"] = th["cref_d"][-1]
        r["cref_vert"] = sum(th["cref_n_vert"])
        r["cref_vert_max"] = max(th["cref_n_vert"])
        r["cref_clamp"] = sum(th["cref_n_clamp"])
        r["cref_neg"] = sum(th["cref_n_neg"])
        r["cref_lat"] = sum(th["cref_n_lat"])
    return r


def show(r):
    P = print
    g = r["g"]
    P("=" * 96)
    P(f"  LEG {r['leg']}   block {r['H']*1e3:.0f} mm   n_cell {g['nx']}x{g['ny']}x{g['nz']}"
      f"   dx {g['dx']*1e3:.3f} mm   dt {g['dt']}   thermo every {g['tint']} steps"
      f"   constraint_depth {g['cref'] or 'off'}")
    P(f"    run to t = {r['t_end']:.3f} s (thermo), end plotfile {r['pf_end']} = {r['t_pf_end']:.3f} s")
    P(f"    VOLUME x4   event log {r['V']:9.3f} cm3   void mask {r['V_mask']:9.3f} cm3"
      f"   ({abs(r['V']-r['V_mask'])/r['V']*100:.2f} %)   / her 12.9 = {r['V']/HER_V0:.2f}x")
    P(f"    ENERGY x4   absorbed ledger_E_src {r['E_src']/1e3:9.3f} kJ"
      f"   flakes ledger_E_flake {r['E_fl']/1e3:9.3f} kJ   share {r['share']*100:5.1f} %"
      f"   price {r['price']:.4f} GJ/m3")
    P(f"    LEDGER      ledger_err_solid max |.| {r['err']:.3e}  (relative, SR 7 ii)")
    P(f"    CRITERION   regime-1 firing T median {r['fire_med']:.2f} K  (n = {r['n_fire']})")
    P(f"    FACE        top-solid-cell mean at r < 10 mm (P2-5 c12 estimator), per plotfile:")
    P("                " + "  ".join(f"{f[0]:4.0f}s {f[1]:7.2f}" for f in r["faces"]))
    P(f"                time-mean over {len(r['faces'])} plotfiles {r['face_tmean']:.2f} K"
      f"   (the top cell's CENTRE, dz/2 below the face)")
    P(f"    FACE recon  T_surf_field (the criterion's reconstructed face T) on the same cells:")
    P("                " + "  ".join(f"{f[0]:4.0f}s {f[3]:7.2f}" for f in r["faces"]))
    P(f"                time-mean {r['recon_tmean']:.2f} K   (cells with a reconstruction per"
      f" plotfile: {', '.join(str(f[4]) for f in r['faces'])}; fallback "
      f"{', '.join(str(f[5]) for f in r['faces'])})")
    P(f"    PROFILE at {r['t_pf_end']:.0f} s from the void mask, 2.5 mm rings (max / median"
      f" depth [mm], columns):")
    for b, (mx, md, nc) in r["rings"].items():
        P(f"      r {b*2.5:5.1f}-{(b+1)*2.5:5.1f} mm   max {mx*1e3:7.2f}   median {md*1e3:7.2f}   n {nc:4d}")
    P(f"      axis (ring 0 median) {r['axis']*1e3:.2f} mm   column (0,0) {r['axis00']*1e3:.2f} mm")
    P(f"      depth at r = 17.5 mm  {r['d175']*1e3:.2f} mm (mean of the 16.25/18.75 ring medians)"
      f"  => ratio {r['dr']:.4f}")
    P(f"                            {r['d175c']*1e3:.2f} mm (median of columns at 16.25-18.75 mm)"
      f"  => ratio {r['drc']:.4f}")
    P(f"      event-log cross-check (SR 4: post-top = k_top - n_voided, regime 1): "
      f"{r['ev_agree']}/{r['ev_cols']} columns agree with the mask to 0.5 dz")
    P(f"    DEEPEST     column {r['deep_col']*1e3:.2f} mm; any void cell to {r['deep_any']*1e3:.2f} mm"
      f" of the {r['H']*1e3:.0f} mm block")
    s = [x for x in r["sec"] if x is not None]
    P(f"    BORE Ø at 10 mm (level k = {r['k10']}), per 30-degree sector: "
      + "  ".join("solid" if x is None else f"{x:6.2f}" for x in r["sec"])
      + (f"   spread {max(s)-min(s):.2f} mm" if s else ""))
    if "cref_d" in r:
        P(f"    CONSTRAINT DEPTH  cref_d {r['cref_d']:g} m   cref_n_vert sum {r['cref_vert']:.6g}"
          f" (per-step max {r['cref_vert_max']:.6g})   lat {r['cref_lat']:.6g}"
          f"   clamp {r['cref_clamp']:.6g}   neg {r['cref_neg']:.6g}")
        P(f"      ⚠ the cref_n_* columns are PER-STEP counters printed every"
          f" {g['tint']} steps: the sums are a SAMPLE, used only as > 0 witnesses (SR 7 i)")


def verdicts(R):
    P = print
    Vb, Sb, z175 = EB.bowl()
    her_dr = z175 / EB.BOWL_AXIS
    E850 = [sum(EB.E_terms(850.0, 12.9e-6, A, fl)) for A in (2.0e-3, 2.5e-3)
            for fl in ("const", "grimsel")]
    lo850, hi850 = min(E850), max(E850)
    P("=" * 96)
    P("  CRITERIA (ord 2-3)")
    out = {}
    if "A0" in R:
        a0 = R["A0"]
        P(f"  (A0 reference: V {a0['V']:.3f} cm3, E_src {a0['E_src']/1e3:.3f} kJ, share"
          f" {a0['share']*100:.1f} %, price {a0['price']:.4f} -- P2-5 A read 100.592 cm3, 1.1879)")
    HA, HG, HGd = R.get("HA"), R.get("HG"), R.get("HGd")
    if HA:
        ok = 90.0 <= HA["V"] <= 130.0
        P(f"  c1  HA volume {HA['V']:.3f} cm3 = {HA['V']/HER_V0:.2f}x her 12.9   bar 90-130"
          f"   {'PASS' if ok else 'MISS'}")
    if HA and HG:
        q = HG["V"] / HA["V"]
        P(f"  c2  HG / HA volume {HG['V']:.3f} / {HA['V']:.3f} = {q:.4f}   bar 0.80-0.95"
          f"   {'PASS' if 0.80 <= q <= 0.95 else 'MISS'}")
        P(f"      (energy ratio HG/HA {HG['E_src']/HA['E_src']:.4f}; the Goldstein mean over"
          f" Ø40 is 0.8819 x the uniform h)")
    for L in (HA, HG):
        if not L:
            continue
        rmin, rmax = L["E_src"] / hi850, L["E_src"] / lo850
        ok = rmin >= 2.0 and L["share"] >= 0.9
        P(f"  c3  {L['leg']:3s} absorbed {L['E_src']/1e3:.2f} kJ vs her 0 MPa E_max(850 K)"
          f" {lo850/1e3:.2f}-{hi850/1e3:.2f} kJ => {rmin:.2f}-{rmax:.2f}x (bar >= 2);"
          f" flake share {L['share']*100:.1f} % (bar >= 90)   {'PASS' if ok else 'MISS'}")
        for T, lab in ((L["fire_med"], "firing T"), (L["face_tmean"], "time-mean face T"),
                       (L["recon_tmean"], "reconstructed face T")):
            Es = [sum(EB.E_terms(T, 12.9e-6, A, "const")) for A in (2.0e-3, 2.5e-3)]
            P(f"        at the model's own {lab} {T:.1f} K her E_max is {min(Es)/1e3:.2f}-"
              f"{max(Es)/1e3:.2f} kJ => {L['E_src']/max(Es):.2f}-{L['E_src']/min(Es):.2f}x")
    for L, bar in ((HA, 0.9), (HG, 0.8)):
        if not L:
            continue
        P(f"  c4  {L['leg']:3s} depth ratio at r = 17.5 mm {L['dr']:.4f} (centred ring"
          f" {L['drc']:.4f}) vs her 47 MPa {her_dr:.4f}   bar >= {bar}"
          f"   {'PASS' if L['dr'] >= bar else 'MISS'}")
    if HG and HGd:
        q = HGd["V"] / HG["V"]
        P(f"  c5  HGd / HG volume {HGd['V']:.3f} / {HG['V']:.3f} = {q:.4f}   bar 0.90-1.05"
          f"   {'PASS' if 0.90 <= q <= 1.05 else 'MISS'}")
    for L in (HA, HG, HGd):
        if not L:
            continue
        dp = max(L["deep_col"], L["deep_any"])
        P(f"  c6  {L['leg']:3s} deepest {dp*1e3:.2f} mm of {L['H']*1e3:.0f}   bar <= 140"
          f"   {'PASS' if dp <= 0.140 + 1e-12 else 'MISS'}")
    for L in (HA, HG, HGd):
        if not L:
            continue
        okf = 850.0 <= L["fire_med"] <= 870.0
        okc = 700.0 <= L["face_tmean"] <= 870.0
        P(f"  c7  {L['leg']:3s} firing T {L['fire_med']:.2f} K (bar 850-870 {'ok' if okf else 'MISS'});"
          f" time-mean face {L['face_tmean']:.2f} K (bar 700-870 {'ok' if okc else 'MISS'})"
          f"   {'PASS' if okf and okc else 'MISS'}")
        P(f"        reconstructed face (T_surf_field) time-mean {L['recon_tmean']:.2f} K"
          f"   {'inside' if 700.0 <= L['recon_tmean'] <= 870.0 else 'OUTSIDE'} 700-870")
    for L in (HA, HG, HGd):
        if not L:
            continue
        s = [x for x in L["sec"] if x is not None]
        sp = (max(s) - min(s)) if len(s) == N_SEC else float("inf")
        ok = L["err"] <= 5e-3 and sp <= 3.0
        P(f"  c8  {L['leg']:3s} ledger_err_solid max {L['err']:.3e} (bar <= 5e-3); Ø at 10 mm by"
          f" sector {', '.join('solid' if x is None else f'{x:.2f}' for x in L['sec'])}"
          f" spread {sp:.2f} mm (bar <= 3)   {'PASS' if ok else 'MISS'}")
    if HGd and "cref_d" in HGd:
        ok = abs(HGd["cref_d"] - 0.002) < 1e-12 and HGd["cref_vert"] > 0
        P(f"  L3  HGd cref_d {HGd['cref_d']:g}, cref_n_vert sampled sum {HGd['cref_vert']:.6g}"
          f"   {'REACHED' if ok else 'NOT REACHED'}")
    P("=" * 96)


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "GATE":
        sys.exit(gate())
    names = sys.argv[1:] or list(LEGS)
    R = {}
    print("P2-7 -- Hoeser at measured heat: the legs")
    print(f"  legs dir {os.path.relpath(os.path.join(HERE, 'legs'), ROOT)};"
          f" every leg's argv.sh + DROPPED.txt beside its run (SR 8)")
    for n in names:
        if not os.path.exists(os.path.join(legdir(n), ".done")):
            print(f"  {n}: no .done sentinel -- skipped")
            continue
        R[n] = one(n)
        show(R[n])
    verdicts(R)


if __name__ == "__main__":
    main()
