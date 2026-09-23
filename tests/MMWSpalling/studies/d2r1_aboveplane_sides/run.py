#!/usr/bin/env python3
"""D2r-1: does heating the SIDE faces of above-plane columns close the mesh gap?

  run.py --cases B_2mm B_1mm T_2mm T_1mm S_2mm --jobs 4
  run.py --cases B_2mm --jobs 1        Goal item 1: the baseline-validity proof
  run.py --list
  run.py --kill NAME                   never pkill -f

Imports studies/d2j0b_plane/run.py by path (FROZEN; never edited) and reuses its
C1 arbiter keys verbatim. The above-plane (h, T) pair is a new WALL entry added
IN MEMORY, so the frozen file on disk is untouched.

The one variable
----------------
Every leg carries side_face_flux = 1 and side_face_factor = 0.2 below the plane
(the frozen stand-in that holds the baseline's treatment constant).

  B   above-plane h = 1e-3  (C1 verbatim), key UNSET  -> must reproduce D2l's
      A_f02 / A_f02_1mm byte for byte. That is Goal item 1, and it is what makes
      the frozen 6 / 8 / 14 % a valid baseline.
  T   above-plane h = 142, side_face_factor_above = 1.0  -> the test.
  S   above-plane h = 142, side_face_factor_above = 0.0  -> the separator:
      identical to T except the above-plane SIDE faces are forced to zero, so
      T - S is exactly those faces and nothing else.

142 is the annulus area-expansion derivation's value (anchored to the floor's
605 W/m^2K at r = 45.8 mm). It is NOT a dial and is never varied to get a result.
No absolute number here is physics -- only the differences T - B and T - S.
"""
import argparse
import concurrent.futures as cf
import importlib.util
import os
import sys

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

H_ABOVE = 142.0                  # derived, declared, never tuned
F_BELOW = 0.2                    # frozen stand-in, never scanned
# In-memory only: the frozen d2j0b_plane/run.py is not touched.
B.WALL.setdefault("TA", (H_ABOVE, B.T_IN))

LEGS = {                         # name -> (d2j0b case, side_face_factor_above)
    "B_2mm": ("C1", None), "B_1mm": ("C1", None),
    "T_2mm": ("TA", 1.0),  "T_1mm": ("TA", 1.0),
    "S_2mm": ("TA", 0.0),  "S_1mm": ("TA", 0.0),
}
CASES = list(LEGS)
FROZEN = {"B_2mm": os.path.join(STUD, "d2l_scan", "output", "A_f02"),
          "B_1mm": os.path.join(STUD, "d2l_scan", "output", "A_f02_1mm")}


def job(name):
    case, above = LEGS[name]
    mesh = name.split("_")[1]
    _, cmd, meta = B.job(f"{case}_{mesh}")
    pf = os.path.join(OUT, name)
    cmd = [f"plot_file={pf}" if c.startswith("plot_file=") else c for c in cmd]
    cmd += ["surface_patch.side_face_flux=1",
            f"surface_patch.side_face_factor={F_BELOW!r}"]
    if above is not None:
        cmd += [f"surface_patch.side_face_factor_above={float(above)!r}"]
    meta = dict(meta, name=name, leg=case, side_face_flux=1,
                side_face_factor=F_BELOW, side_face_factor_above=above,
                h_above=B.WALL[case][0], T_above=B.WALL[case][1])
    return pf, cmd, meta


def run_all(names, jobs, force):
    J.keep_awake()
    js = [job(n) for n in names]
    print(f"{len(js)} runs, {min(jobs, len(js))} in parallel -> {OUT}", flush=True)
    for _, _, m in js:
        print(f"  {m['name']:8s} dz {m['dz'] * 1e3:g} mm, dt {m['dt'] * 1e3:g} ms, "
              f"stop {m['stop']:g} s, h_above {m['h_above']:g}, "
              f"f_above {m['side_face_factor_above']}", flush=True)
    ok = True
    with cf.ThreadPoolExecutor(max_workers=min(jobs, len(js))) as ex:
        for fu in cf.as_completed([ex.submit(J.run_one, pf, c, m, force) for pf, c, m in js]):
            m, status, wall = fu.result()
            print(f"  {m['name']:8s} {status} at t = {m.get('t_end', 0):g} s "
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
        for n in CASES:
            _, cmd, m = job(n)
            print(f"{n:8s} dz {m['dz'] * 1e3:g} mm {m['nxy']}^2x{m['nz']} "
                  f"dt {m['dt'] * 1e3:g} ms stop {m['stop']:g} s")
            print("      ", " ".join(c for c in cmd if "side_face" in c or "h_expr" in c)[:150])
        return
    bad = [n for n in (a.cases or []) if n not in CASES]
    if not a.cases or bad:
        sys.exit(f"give --cases from {CASES}" + (f" (unknown: {bad})" if bad else ""))
    sys.exit(0 if run_all(a.cases, a.jobs, a.force) else 1)


if __name__ == "__main__":
    main()
