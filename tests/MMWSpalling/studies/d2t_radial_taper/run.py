#!/usr/bin/env python3
"""D2t: retire the mask edge with a smooth radial taper, and score SHAPE.

  run.py --cases P_2mm P_1mm PA_2mm PA_1mm --jobs 4
  run.py --cases TR_below TR_above --jobs 2    item 2: the compiled-taper trace
  run.py --list
  run.py --kill NAME                           never pkill -f

Imports studies/d2j0b_plane/run.py by path (FROZEN; never edited) and reuses its
arbiter keys.  The `TA` wall entry (h = 142 above the nozzle plane) is added IN
MEMORY, exactly as D2r-1 did, so the frozen file on disk is untouched.

The one change: `h_expr` gains a smooth radial taper
-----------------------------------------------------
D2j-0b/D2l/D2r-1 all ran a HARD mask: h at full strength at r = 29.9 mm and
zero at r = 30.1 mm.  The residual 6/8/14 % mesh gap is entirely that cliff
(core r < 25 mm is already converged to -0.2/+3.0/-1.8 % at f = 0.2).  This
packet removes the cliff instead of heating it:

    h(r, s) = taper(r) * if(s > 0, 700, 142)

    taper = F0 + (1-F0) * 0.5 * (1 + cos(pi * (clamp(r) - R_FLAT) / W))

a raised cosine (C1 at both ends, unlike a linear ramp whose kinks are two
weaker corners), full strength to R_FLAT = 18 mm, at its floor by
R_TAPER = 34 mm.  W = 16 mm = 8 coarse cells, so the taper is resolved at BOTH
meshes; it is a fixed PHYSICAL width, never a cell count.

surface_patch.radius = 50 mm, well beyond R_TAPER
-------------------------------------------------
A taper in h_expr does NOT remove the hard mask: `r2 > r2max -> continue`
(:1745) and `in_patch` (:1579) and `col_in_patch` (:1626) all cut at
surface_patch.radius.  Leaving radius at 30 mm would MOVE the cliff, not
remove it.  50 mm puts the binary cutoff 16 mm into the floor region, where
h is 1e-3 and the face-form/cell-form switch at :2653 is meaningless.

The floor is on h, not on the multiplier
----------------------------------------
:1793 aborts unless h is finite and > 0, so the taper cannot reach zero.  The
packet says "a floor of 1.0e-3" and cites C1, whose 1.0e-3 is a value of h.
Read as a multiplier floor instead, h would land at 0.7 W/m^2K, which against
T_gas = 1600 K equilibrates the unheated far field at ~350 K rather than 293 K
-- not "emulating h = 0".  So F0 = 1e-3/700 and h lands on 1e-3 below the
plane, 2.03e-4 above.  Recorded in PREFLIGHT.md as a reading, not a silent choice.
"""
import argparse
import concurrent.futures as cf
import importlib.util
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)


def load(name, path):
    sp = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


B = load("d2j0b_run", os.path.join(STUD, "d2j0b_plane", "run.py"))
J = B.J

OUT = os.path.join(HERE, "output")
B.OUT = OUT
J.OUT = OUT

H_IN = B.H_IN                     # 700, below the plane, unchanged
H_ABOVE = 142.0                   # the annulus derivation; never tuned
F_BELOW = 0.2                     # frozen stand-in, never scanned
R_FLAT, W, R_PATCH = 0.018, 0.016, 0.050
R_TAPER = R_FLAT + W        # 0.034; W is the primitive so the literal is exact
H_FLOOR = 1.0e-3                  # the same value C1 uses, on h
F0 = H_FLOOR / H_IN

# In-memory only: the frozen d2j0b_plane/run.py is not touched.
B.WALL.setdefault("TA", (H_ABOVE, B.T_IN))

# Numeric first in every min/max (D2o-0b's rule); the denominator is a literal,
# never a max() -- the AMReX parser evaluates a/max(x, num) as 0.
TAPER = (f"({F0!r}+{1.0 - F0!r}*0.5*(1.0+cos(3.141592653589793*"
         f"(min({R_TAPER!r},max({R_FLAT!r},r))-{R_FLAT!r})/{W!r})))")


def taper_of(r):
    """The same raised cosine in numpy.  A FILTER ONLY: the D2o-0 division bug
    lived inside the AMReX parser, so this cannot prove the compiled profile.
    That is what taper_trace.py is for."""
    rc = np.clip(r, R_FLAT, R_TAPER)
    return F0 + (1.0 - F0) * 0.5 * (1.0 + np.cos(np.pi * (rc - R_FLAT) / W))


# name -> (mesh, side_face_factor_above, key OVERRIDES)
# Overrides REPLACE the base key in place rather than being appended, so the
# result never depends on how ParmParse resolves a duplicated name.
LEGS = {
    "P_2mm":  ("2mm", None, {}),
    "P_1mm":  ("1mm", None, {}),
    "PA_2mm": ("2mm", 1.0,  {}),
    "PA_1mm": ("1mm", 1.0,  {}),
    # Item 2, the compiled-taper thermal trace: 10 steps at 2 mm, stopped long
    # before the first spall (the frozen logs put that at 4.58 s).  TR_above
    # drops the nozzle plane below the rock so the SAME expression is read on
    # its other branch; the collision guard is switched off for it because
    # every column is then above the plane by construction.
    "TR_below": ("2mm", None, {"stop_time": "0.16", "amr.plot_int": "10"}),
    "TR_above": ("2mm", None, {"stop_time": "0.16", "amr.plot_int": "10",
                               "surface_patch.nozzle_z0": "0.14",
                               "surface_patch.nozzle_collision_radius": "0.0"}),
}
CASES = list(LEGS)
TRACE = ["TR_below", "TR_above"]


def job(name):
    mesh, above, over = LEGS[name]
    _, cmd, meta = B.job(f"TA_{mesh}")
    pf = os.path.join(OUT, name)
    h_expr = f'surface_patch.h_expr="{TAPER}*if(s>0.0,{H_IN!r},{H_ABOVE!r})"'
    out = []
    for c in cmd:
        if c.startswith("plot_file="):
            c = f"plot_file={pf}"
        elif c.startswith("surface_patch.h_expr="):
            c = h_expr
        elif c.startswith("surface_patch.radius="):
            c = f"surface_patch.radius={R_PATCH!r}"
        out.append(c)
    out += ["surface_patch.side_face_flux=1",
            f"surface_patch.side_face_factor={F_BELOW!r}"]
    if above is not None:
        out += [f"surface_patch.side_face_factor_above={float(above)!r}"]
    for k, v in over.items():
        pre = k + "="
        hit = [n for n, c in enumerate(out) if c.startswith(pre)]
        if len(hit) > 1:
            raise SystemExit(f"{name}: {k} appears {len(hit)} times in the base keys")
        if hit:
            out[hit[0]] = pre + v
        else:
            out.append(pre + v)
    meta = dict(meta, name=name, patch_r=R_PATCH, h_expr=h_expr.split("=", 1)[1],
                r_flat=R_FLAT, r_taper=R_TAPER, taper_w=W, h_floor=H_FLOOR, f0=F0,
                h_in=H_IN, h_above=H_ABOVE, side_face_flux=1,
                side_face_factor=F_BELOW, side_face_factor_above=above,
                trace=name in TRACE)
    if name in TRACE:
        meta = dict(meta, stop=0.16, plot_int=10)
    return pf, out, meta


def run_all(names, jobs, force):
    J.keep_awake()
    js = [job(n) for n in names]
    print(f"{len(js)} runs, {min(jobs, len(js))} in parallel -> {OUT}", flush=True)
    for _, _, m in js:
        print(f"  {m['name']:9s} dz {m['dz'] * 1e3:g} mm, dt {m['dt'] * 1e3:g} ms, "
              f"stop {m['stop']:g} s, patch_r {m['patch_r'] * 1e3:g} mm, "
              f"f_above {m['side_face_factor_above']}", flush=True)
    ok = True
    with cf.ThreadPoolExecutor(max_workers=min(jobs, len(js))) as ex:
        for fu in cf.as_completed([ex.submit(J.run_one, pf, c, m, force) for pf, c, m in js]):
            m, status, wall = fu.result()
            print(f"  {m['name']:9s} {status} at t = {m.get('t_end', 0):g} s "
                  f"({wall / 60:.1f} min)", flush=True)
            ok = ok and not status.startswith("FAILED")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", nargs="*")
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--kill")
    a = ap.parse_args()
    if a.kill:
        return J.kill(a.kill)
    if a.list:
        print(f"binary {J.binary_sha()[:8]}...  (provenance only, NOT an oracle)")
        print(f"taper: flat to {R_FLAT * 1e3:g} mm, floor {H_FLOOR:g} by "
              f"{R_TAPER * 1e3:g} mm, W {W * 1e3:g} mm = {W / 2e-3:g} coarse cells; "
              f"patch_r {R_PATCH * 1e3:g} mm = r_taper + {(R_PATCH - R_TAPER) / 2e-3:g} coarse cells")
        print(f"h_expr = {TAPER}*if(s>0.0,{H_IN!r},{H_ABOVE!r})")
        for n in CASES:
            _, cmd, m = job(n)
            print(f"{n:9s} dz {m['dz'] * 1e3:g} mm {m['nxy']}^2x{m['nz']} "
                  f"dt {m['dt'] * 1e3:g} ms stop {m['stop']:g} s")
        return
    bad = [n for n in (a.cases or []) if n not in CASES]
    if not a.cases or bad:
        sys.exit(f"give --cases from {CASES}" + (f" (unknown: {bad})" if bad else ""))
    sys.exit(0 if run_all(a.cases, a.jobs, a.force) else 1)


if __name__ == "__main__":
    main()
