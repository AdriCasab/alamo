#!/usr/bin/env python3
"""P2-5 rev 2, deliverable 5 -- the reductions.

⚠⚠ V_cell IS PARSED FROM THE LEG'S OWN argv.sh AND PRINTED. `reduce_price.py`
hard-codes `V_CELL = 8e-9` at module level (its line 51), which is 2 mm only, so
it is NOT imported here -- the packet bars it for a 1 mm leg and that is the LZ
trap (standing rule 5) in another costume. Its event-census LOGIC is mirrored:
every volume is counted from `n_voided`, one cell being dx*dy*dz.

⚠ Volumes are one QUADRANT and are reported x4 (standing rule 3). A RATIO of two
quadrant quantities is left unscaled (both sides or neither).
⚠ Per-step thermo columns are averaged or integrated, never summed as a total
(standing rule 7 i); `patch_P_robin` is per-step power, so the delivered energy
is sum(P)*dt_sample with dt_sample = timestep*amr.thermo.plot_int.
⚠ The bore/cavity is read from the VOID MASK, not the event log: the log's
`k_top` is the top BEFORE that removal (standing rule 4), and the mask is the
state. Reported PER AZIMUTH SECTOR (standing rule 10).

Usage: reduce_p2_5.py <leg> [<leg> ...]
"""
import csv
import math
import os
import shlex
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "tests", "MMWSpalling", "studies",
                                "p2_2_props"))
import fabcmp                                                   # noqa: E402
from p2_5_legs import LEGS, SHAPES, D_NOZ, Z_ROCK, s_lid, s_leg  # noqa: E402

R_FLOOR, R_LAT = 1, 6          # regime labels, from Removal.H
T_CLAMP = 873.15               # both Grimsel property tables floor here
MELT = 1400.0
N_SEC = 9                      # azimuth sectors across the quadrant
HER_V30 = 13.0                 # cm3 in 30 s at 0 MPa (26 cm3/min)


def geom_of(leg):
    """(nx, ny, nz, dx, dz, V_cell, dt, thermo_int, V0) from the LEG'S OWN argv
    -- the same file the leg was launched from, so it cannot drift."""
    toks = shlex.split(open(os.path.join(HERE, "legs", leg, "argv.sh")).read())
    kv, i = {}, 0
    while i < len(toks):
        t = toks[i]
        if "=" in t:
            k, v = t.split("=", 1)
            i += 1
            while i < len(toks) and "=" not in toks[i]:
                v += " " + toks[i]
                i += 1
            kv[k] = v
        else:
            i += 1
    nx, ny, nz = (int(v) for v in kv["amr.n_cell"].split())
    lo = [float(v) for v in kv["geometry.prob_lo"].split()]
    hi = [float(v) for v in kv["geometry.prob_hi"].split()]
    dx, dz = (hi[0] - lo[0]) / nx, (hi[2] - lo[2]) / nz
    return dict(nx=nx, ny=ny, nz=nz, dx=dx, dz=dz, lo=lo, hi=hi,
                V_cell=dx * dx * dz, dt=float(kv["timestep"]),
                tint=int(kv["amr.thermo.plot_int"]), V0=kv["weibull.V0"],
                core=float(kv["surface_patch.jet_core_length"]),
                z0=float(kv["surface_patch.nozzle_z0"]))


def thermo(leg):
    L = open(os.path.join(HERE, "legs", leg, "run", "thermo.dat")).read().split("\n")
    hdr = L[0].split()
    rows = [l.split() for l in L[1:] if l.strip()]
    out = {n: [] for n in hdr}
    for r in rows:
        if float(r[hdr.index("time")]) <= 0:
            continue              # the t = 0 row is all zeros pre-advance
        for n in hdr:
            out[n].append(float(r[hdr.index(n)]))
    return out


def events(leg):
    p = os.path.join(HERE, "legs", leg, "run_removal_events.csv")
    with open(p) as f:
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


def pctl(xs, q):
    if not xs:
        return float("nan")
    x = sorted(xs)
    return x[min(len(x) - 1, max(0, int(round(q * (len(x) - 1)))))]


def mask(leg, which=-1):
    """(removed, Temp, grid) from a cell plotfile. `which` indexes the sorted
    plotfile list (-1 = last)."""
    run = os.path.join(HERE, "legs", leg, "run")
    pfs = sorted(d for d in os.listdir(run) if d.endswith("cell"))
    if not pfs:
        return None
    lvl = os.path.join(run, pfs[which], "Level_0")
    hdr = os.path.join(run, pfs[which], "Header")
    files = sorted(f for f in os.listdir(lvl) if f.startswith("Cell_D_"))
    first = fabcmp.read_fabs_boxed(os.path.join(lvl, files[0]))
    names = fabcmp.var_names(hdr, first[0]["ncomp"])
    ir, it = names.index("removed"), names.index("Temp")
    g = geom_of(leg)
    rem = {}
    tmp = {}
    for fn in files:
        for fab in fabcmp.read_fabs_boxed(os.path.join(lvl, fn)):
            lo, hi = fab["lo"], fab["hi"]
            nxb = hi[0] - lo[0] + 1
            nyb = hi[1] - lo[1] + 1
            nzb = hi[2] - lo[2] + 1
            # ⚠ read_fabs_boxed returns {"comp": [list-of-float per
            # component]}, NOT a flat "data" array; AMReX FAB ordering is i
            # fastest, then j, then k.
            cr, ct = fab["comp"][ir], fab["comp"][it]
            for kk in range(nzb):
                for jj in range(nyb):
                    base = nxb * (jj + nyb * kk)
                    for ii in range(nxb):
                        o = ii + base
                        key = (lo[0] + ii, lo[1] + jj, lo[2] + kk)
                        rem[key] = cr[o]
                        tmp[key] = ct[o]
    return rem, tmp, g, pfs[which]


def cavity(leg):
    """Depth on the axis and width per azimuth sector, from the void mask."""
    m = mask(leg)
    if m is None:
        return None
    rem, tmp, g, tag = m
    nx, ny, nz, dx, dz = g["nx"], g["ny"], g["nz"], g["dx"], g["dz"]
    ztop = g["hi"][2]
    # per column: how many top cells are void => the recession depth
    depth = {}
    for j in range(ny):
        for i in range(nx):
            d = 0
            for k in range(nz - 1, -1, -1):
                if rem.get((i, j, k), 0.0) > 0.5:
                    d += 1
                else:
                    break
            depth[(i, j)] = d * dz
    axis = depth[(0, 0)]
    # width at the ORIGINAL surface and 5 mm below it, per azimuth sector
    def width_at(zoff):
        kk = int(round((zoff) / dz))      # cells below the original surface
        per = [0.0] * N_SEC
        for j in range(ny):
            for i in range(nx):
                if depth[(i, j)] <= kk * dz:
                    continue
                x, y = (i + 0.5) * dx, (j + 0.5) * dx
                r = math.hypot(x, y)
                sec = min(N_SEC - 1, int(math.atan2(y, x) / (math.pi / 2)
                                         * N_SEC))
                per[sec] = max(per[sec], r)
        return per
    return dict(axis=axis, w0=width_at(0.0), w5=width_at(0.005),
                depth=depth, g=g, tag=tag,
                vol_mask=sum(depth.values()) * dx * dx)


def one(leg):
    P = print
    g = geom_of(leg)
    dzl, st, K, z0, shape = LEGS[leg]
    sig = None if SHAPES[shape] is None else SHAPES[shape] * D_NOZ
    th = thermo(leg)
    ev = events(leg)
    dts = g["dt"] * g["tint"]

    P("=" * 94)
    P(f"  LEG {leg}   shape {shape}   K {g['core']:g}   nozzle_z0 {g['z0']}"
      f"   stop_time {st:g} s")
    P(f"  GEOMETRY PARSED FROM legs/{leg}/argv.sh (never a module constant):")
    P(f"    n_cell {g['nx']}x{g['ny']}x{g['nz']}   dx {g['dx']*1e3:.3f} mm"
      f"   dz {g['dz']*1e3:.3f} mm   V_cell {g['V_cell']:.4g} m3"
      f"   weibull.V0 {g['V0']}   dt {g['dt']} s   thermo every {dts:.3f} s")
    if abs(g["V_cell"] - float(g["V0"])) / g["V_cell"] > 0.01:
        P(f"    ⚠ weibull.V0 {g['V0']} != V_cell {g['V_cell']:.4g}")

    # ---- volumes: event log vs mask ------------------------------------
    nv = sum(_i(r, "n_voided") for r in ev)
    v_ev = nv * g["V_cell"] * 4.0 * 1e6          # cm3, x4 for the quarter
    cv = cavity(leg)
    v_mk = cv["vol_mask"] * 4.0 * 1e6 if cv else float("nan")
    rel = abs(v_ev - v_mk) / max(v_ev, 1e-30)
    P(f"\n  VOLUME at {st:g} s (x4, SR 3)")
    P(f"    event log  {v_ev:9.3f} cm3   ({nv} cells voided, all regimes)")
    P(f"    void mask  {v_mk:9.3f} cm3   (plotfile {cv['tag'] if cv else '-'})")
    P(f"    agreement  {rel*100:6.2f} %   (b1 bar <= 5 %)"
      f"   {'ok' if rel <= 0.05 else 'MISS'}")
    P(f"    her 0 MPa 30 s = {HER_V30} cm3  =>  model / Hoeser ="
      f" {v_ev/HER_V30:6.2f}x")

    # ---- energy ----------------------------------------------------------
    # ⚠⚠ THE DELIVERED ENERGY IS `ledger_E_src`, NOT int patch_P_robin dt.
    # `patch_P_robin` is the EXPRESSION PATCH's own Robin power -- the floor map
    # alone. `jet_P_face` additionally carries the side-face term the march
    # charges into it (the D2ab note at MMWSpalling.H:2530ff), and
    # `ledger_E_src` is the energy ledger's own cumulative source, i.e. the
    # authoritative input. On G05 they read 3.54 / 11.01 / 10.69 kJ per
    # quarter: scoring the flake share against the patch term alone gave
    # 124.6 %, which is how this was caught. `ledger_E_src` is CUMULATIVE, so
    # it is read at its final value and never summed (SR 7 i).
    E_src = th["ledger_E_src"][-1] * 4.0
    E_patch = sum(th["patch_P_robin"]) * dts * 4.0
    E_jet = sum(th["jet_P_face"]) * dts * 4.0
    E_fl = th.get("ledger_E_flake", [0.0])[-1] * 4.0
    V_m3 = v_ev * 1e-6                     # cm3 -> m3, already x4
    P(f"\n  ENERGY (x4; ledger_E_src and ledger_E_flake are CUMULATIVE)")
    P(f"    ledger_E_src  (authoritative input) {E_src/1e3:9.3f} kJ"
      f"   (mean {E_src/4.0/st:8.1f} W quarter)")
    P(f"    int jet_P_face dt  (floor + side)   {E_jet/1e3:9.3f} kJ")
    P(f"    int patch_P_robin dt  (FLOOR MAP only) {E_patch/1e3:6.3f} kJ"
      f"   = {E_patch/max(E_src,1e-30)*100:5.1f} % of the input")
    P(f"    ⚠ the rest is the TYPED 169.3 W/m2K side-face term on the cavity")
    P(f"      wall ({100-E_patch/max(E_src,1e-30)*100:.1f} %): the derived wall")
    P(f"      law is not ready at this bore, as the packet states, but the")
    P(f"      fallback is NOT inert -- it dominates the energy input.")
    P(f"    flakes        ledger_E_flake final  {E_fl/1e3:9.3f} kJ"
      f"   share {E_fl/max(E_src,1e-30)*100:5.1f} %  (c3 bar 60-85 %)")
    P(f"    price         {E_fl/max(V_m3,1e-30)/1e9:8.4f} GJ/m3"
      f"   (flake energy / removed volume; both sides x4 -- SR 3)")
    E_del = E_src
    cap = max(a / b for a, b in zip(th["jet_P_face"], th["jet_P_cap"]))
    P(f"    max jet_P_face/jet_P_cap {cap:.4f}   (b1 bar < 0.9)"
      f"   {'ok' if cap < 0.9 else 'MISS'}")

    # ---- the closure's standoff ------------------------------------------
    ts = th["jet_T_stag"]
    pm = th["patch_min_standoff"]
    first = next((i for i, t in enumerate(ts) if abs(t - 1866.0) > 0.5), None)
    P(f"\n  CLOSURE   jet_T_stag {ts[0]:.1f} -> {ts[-1]:.1f} K")
    if first is None:
        P(f"    phi stayed clamped at 1 for the whole leg"
          f" (patch_min_standoff {pm[0]:.5f} -> {pm[-1]:.5f} m,"
          f" threshold {g['core']*3.745e-3:.5f} m)")
    else:
        P(f"    jet_T_stag first leaves 1866 K at t = {th['time'][first]:.3f} s,"
          f" patch_min_standoff = {pm[first]:.5f} m"
          f" (threshold {g['core']*3.745e-3:.5f} m)")

    # ---- the criterion ----------------------------------------------------
    fl = [_f(r, "T_surf") for r in ev
          if _i(r, "regime") == R_FLOOR and _f(r, "T_surf") > 0]
    allt = [_f(r, "T_surf") for r in ev if _f(r, "T_surf") > 0]
    P(f"\n  CRITERION   regime-1 rows with T_surf > 0: n = {len(fl)}"
      f" of {len(ev)} events")
    if fl:
        P(f"    firing T  median {pctl(fl,0.5):8.2f} K   p10 {pctl(fl,0.10):8.2f}"
          f"   p90 {pctl(fl,0.90):8.2f}   (c6 bar 845 +- 10)")
    if allt:
        mx = max(allt)
        P(f"    max T_surf over ALL regimes {mx:8.2f} K"
          f"   MELT FLAG (> {MELT:.0f} K): {'YES' if mx > MELT else 'no'}")
    cl = [r for r in ev if _i(r, "regime") == R_FLOOR
          and _f(r, "T_surf") > T_CLAMP]
    nvf = sum(_i(r, "n_voided") for r in ev if _i(r, "regime") == R_FLOOR)
    nvc = sum(_i(r, "n_voided") for r in cl)
    P(f"    above the {T_CLAMP:.2f} K property-table clamp:"
      f" {len(cl)}/{len(fl)} by count, {nvc/max(nvf,1)*100:.1f} % by volume")

    # ---- the periphery share (c4) ----------------------------------------
    from derive_lab_map import r_half, h0_of, sigma_m
    rh = r_half(sigma_m(shape), h0_of(shape))
    if rh is None:
        P(f"\n  PERIPHERY (c4)   shape A has no periphery by construction"
          f" (plateau to {20.0:.0f} mm, q flat)")
    else:
        inb = outb = 0
        for r in ev:
            if _i(r, "regime") != R_FLOOR:
                continue
            i, j = _i(r, "col_i"), _i(r, "col_j")
            rr = g["dx"] * math.hypot(i + 0.5, j + 0.5)
            n = _i(r, "n_voided")
            if rr > rh:
                outb += n
            else:
                inb += n
        tot = inb + outb
        P(f"\n  PERIPHERY (c4)   r_0.5 = {rh*1e3:.2f} mm"
          f" (q(850 K) < 0.5 MW/m2 beyond it)")
        P(f"    regime-1 volume beyond r_0.5: {outb}/{tot} cells ="
          f" {outb/max(tot,1)*100:5.1f} %   (prediction <= 10 % on G05;"
          f" > 25 % => the cut-off is live)")

    # ---- the cavity, per sector (c5) -------------------------------------
    if cv:
        w0, w5 = cv["w0"], cv["w5"]
        nz0 = [w for w in w0 if w > 0]
        P(f"\n  CAVITY from the void mask, {N_SEC} azimuth sectors (SR 10)")
        Hdom = (cv["g"]["hi"][2] - cv["g"]["lo"][2])
        P(f"    axis depth              {cv['axis']*1e3:8.2f} mm"
          f"   (c5 bar 15-40)")
        if cv["axis"] > 0.8 * Hdom:
            P(f"    ⚠⚠ DEPTH GUARD: {cv['axis']*1e3:.1f} mm of a"
              f" {Hdom*1e3:.0f} mm block -- the cavity is within"
              f" {(Hdom-cv['axis'])*1e3:.0f} mm of the ADIABATIC base, so the")
            P(f"       late part of this leg is contaminated (ROADMAP:"
              f" 'a prescribed-feed leg can run into the bottom of the block')")
        if nz0:
            P(f"    width at the original surface  sector mean"
              f" {sum(nz0)/len(nz0)*2e3:8.2f} mm   spread"
              f" {(max(nz0)-min(nz0))*2e3:6.2f} mm   (c5 bar 20-35, spread <= 3)")
        nz5 = [w for w in w5 if w > 0]
        if nz5:
            P(f"    width 5 mm below it            sector mean"
              f" {sum(nz5)/len(nz5)*2e3:8.2f} mm   spread"
              f" {(max(nz5)-min(nz5))*2e3:6.2f} mm")
        # the undercut: void at the original surface OUTSIDE the Ø 40 opening
        und = [math.hypot((i+0.5)*g["dx"], (j+0.5)*g["dx"])
               for (i, j), d in cv["depth"].items()
               if d > 0 and math.hypot((i+0.5)*g["dx"], (j+0.5)*g["dx"]) > 0.020]
        P(f"    UNDERCUT under the lid (void beyond r = 20 mm):"
          f" {'none' if not und else f'to {max(und)*2e3:.2f} mm width'}"
          f"   (c13: her 44-64 mm at 47 MPa)")

    # ---- the Wilkinson row (c12) -----------------------------------------
    P(f"\n  WILKINSON ROW (c12)  time-mean exposed-surface T at r < 10 mm")
    run = os.path.join(HERE, "legs", leg, "run")
    pfs = sorted(d for d in os.listdir(run) if d.endswith("cell"))
    # ⚠ PICK BY TIME, NOT BY FRACTION. Fractions of the plotfile list are only
    # 10/20/30 s on a 30 s leg: on G05T120 they are 40/80/120 s, and the first
    # rev printed them under the 10/20/30 labels. amr.plot_int is in STEPS, so
    # the plotfile tag / plot_int * plot_int * dt recovers the time.
    pint = int(round(0.001 + (st / g["dt"]) / max(len(pfs) - 1, 1)))
    want = [t for t in (10.0, 20.0, 30.0) if t <= st + 1e-9]
    for tw in want:
        lbl = f"{tw:.0f} s"
        idx = min(len(pfs) - 1,
                  max(0, int(round(tw / g["dt"] / max(pint, 1)))))
        mm = mask(leg, idx)
        if mm is None:
            continue
        rem, tmp, gg, tag = mm
        vals = []
        for j in range(gg["ny"]):
            for i in range(gg["nx"]):
                if gg["dx"] * math.hypot(i + 0.5, j + 0.5) >= 0.010:
                    continue
                for k in range(gg["nz"] - 1, -1, -1):
                    if rem.get((i, j, k), 0.0) <= 0.5:
                        vals.append(tmp[(i, j, k)])
                        break
        if vals:
            P(f"    {lbl:>5s} ({tag}): mean {sum(vals)/len(vals):8.2f} K"
              f"   max {max(vals):8.2f} K   n = {len(vals)} surface cells")
    if fl:
        P(f"    firing T median {pctl(fl,0.5):.2f} K  (c12 bar: the time-mean"
          f" sits 50-150 K ABOVE it)")
    return dict(leg=leg, v=v_ev, E=E_del, Efl=E_fl, cv=cv, fl=fl,
                cap=cap, g=g, shape=shape)


def main():
    res = [one(l) for l in sys.argv[1:]]
    print("=" * 94)
    if len(res) > 1:
        print("  SUMMARY")
        print(f"    {'leg':>9s} {'V30 [cm3] x4':>13s} {'/Hoeser':>8s}"
              f" {'E_del [kJ]':>11s} {'flake %':>8s} {'firing T':>9s}")
        for r in res:
            print(f"    {r['leg']:>9s} {r['v']:13.3f} {r['v']/HER_V30:8.2f}x"
                  f" {r['E']/1e3:11.3f} {r['Efl']/max(r['E'],1e-30)*100:8.1f}"
                  f" {pctl(r['fl'],0.5) if r['fl'] else float('nan'):9.2f}")


if __name__ == "__main__":
    main()
