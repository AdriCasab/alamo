#!/usr/bin/env python3
"""D2q-2c arbiter: the D2j-0b C1 cascade with PER-CELL REMOVAL on.

  run.py --cases C1_2mm C1_1mm C0_2mm C0_1mm [--jobs 4] [--force]
  run.py --cases C1_0p5mm --jobs 1        only to support a pass
  run.py --list
  run.py --kill NAME

Imports studies/d2j0b_plane/run.py by path (frozen; never edited) and reuses its
C0/C1 keys verbatim, adding only:

    surface_patch.side_face_flux   = 1
    surface_patch.side_face_h_max  = 80.0     (the default; stated for provenance)
    spall.per_cell_removal         = 1

The question: does per-cell removal convert the corner drain into removal
without drilling faster? The mechanism is live but has not fired -- at 90 s the
lateral Sp maxes at 0.4533, T <= 656.6 K against a 822.2 K firing point, and the
cap needs ~149 s of exposure to carry a wall there. 250 s is the test.

`side_face_h_max` and `side_face_factor` are DERIVED and DECLARED. Neither is
ever adjusted to make something fire -- an inert campaign is a real result (P4).

The 0.5 mm leg needs a mesh the frozen harness does not define, so MESH/DT/MGS
are extended IN MEMORY here (the file on disk is never touched). The dt ladder
continues 16 -> 8 -> 4 ms and max_grid_size 20 -> 40 -> 80, both by the same
factor the frozen harness uses between its own meshes.
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
J = B.J                                   # the d2j0_hotdisk harness underneath

OUT = os.path.join(HERE, "output")
B.OUT = OUT
J.OUT = OUT                               # run_one's makedirs and kill()'s token

# In-memory extension for the 0.5 mm leg; the frozen file is untouched.
B.MESH.setdefault("0p5mm", 0.5e-3)
B.DT.setdefault("0p5mm", 4e-3)
B.MGS.setdefault("0p5mm", 80)

H_MAX = 80.0                              # W/m^2K, declared, never fitted
PAIR = ["C1_2mm", "C1_1mm", "C0_2mm", "C0_1mm"]
FINE = ["C1_0p5mm"]
CASES = PAIR + FINE


def job(name):
    case, mesh = name.split("_")
    pf_b, cmd, meta = B.job(f"{case}_{mesh}")
    pf = os.path.join(OUT, name)
    cmd = [f"plot_file={pf}" if c.startswith("plot_file=") else c for c in cmd]
    cmd += ["surface_patch.side_face_flux=1",
            f"surface_patch.side_face_h_max={H_MAX!r}",
            "spall.per_cell_removal=1"]
    meta = dict(meta, name=name, case=case, side_face_flux=1,
                side_face_h_max=H_MAX, per_cell_removal=1)
    return pf, cmd, meta


def run_all(names, jobs, force):
    J.keep_awake()
    js = [job(n) for n in names]
    print(f"{len(js)} runs, {min(jobs, len(js))} in parallel -> {OUT}", flush=True)
    for _, _, m in js:
        print(f"  {m['name']:10s} dz {m['dz'] * 1e3:g} mm, dt {m['dt'] * 1e3:g} ms, "
              f"stop {m['stop']:g} s, h_max {m['side_face_h_max']:g}, "
              f"per_cell_removal {m['per_cell_removal']}", flush=True)
    ok = True
    with cf.ThreadPoolExecutor(max_workers=min(jobs, len(js))) as ex:
        for fu in cf.as_completed([ex.submit(J.run_one, pf, c, m, force) for pf, c, m in js]):
            m, status, wall = fu.result()
            print(f"  {m['name']:10s} {status} at t = {m.get('t_end', 0):g} s "
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
        print(f"binary {J.binary_sha()[:8]}…  (provenance only, NOT an oracle)")
        for n in CASES:
            _, cmd, m = job(n)
            print(f"{n:10s} dz {m['dz'] * 1e3:g} mm {m['nxy']}²×{m['nz']} "
                  f"dt {m['dt'] * 1e3:g} ms stop {m['stop']:g} s")
            print("      keys:", " ".join(c for c in cmd
                                          if "side_face" in c or "per_cell" in c))
        return
    bad = [n for n in (a.cases or []) if n not in CASES]
    if not a.cases or bad:
        sys.exit(f"give --cases from {CASES}" + (f" (unknown: {bad})" if bad else ""))
    sys.exit(0 if run_all(a.cases, a.jobs, a.force) else 1)


if __name__ == "__main__":
    main()
