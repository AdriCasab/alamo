#!/usr/bin/env python3
"""D2k arbiter: the D2j-0b C1 cascade with side-face flux on.

  run.py --cases C1_2mm C1_1mm C1f07_2mm [--jobs 3] [--force]
  run.py --cases LI0_2mm LI0_1mm [--jobs 2]     conditional Meier pair
  run.py --list
  run.py --kill NAME

Imports studies/d2j0b_plane/run.py by path (frozen; not edited) and reuses its
C1 keys verbatim, adding only surface_patch.side_face_flux = 1 and
side_face_factor. C1 is the case whose 1 mm leg cascades: disk-mean gap by
block 50-100 / 100-150 / 150-200 s = -13.5 / -28.7 / -45.5 %, r_front walking
26 -> 22 -> 18 mm. The question is whether heating the exposed vertical faces
removes that.

  C1_2mm / C1_1mm   f = 1.0, the pre-registered value
  C1f07_2mm         f = 0.7, the pre-registered sensitivity (a staircase of
                    one-cell steps has area dx + dz per column against a true
                    sloped sqrt(dx^2 + dz^2), so f = 1 over-heats by up to
                    41 % at 45 deg and ~16 % at the observed 80 deg cone).
                    f is NEVER fitted.

The Meier pair (LI0_*) is the D2i keys with the key on, to 450 s, and runs
only if the arbiter passes; it goes through the d2i_gate harness so the keys
are exactly D2i's.
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

FACTOR = {"C1": 1.0, "C1f07": 0.7}
ARBITER = ["C1_2mm", "C1_1mm", "C1f07_2mm"]
MEIER = ["LI0_2mm", "LI0_1mm"]
CASES = ARBITER + MEIER
MEIER_STOP = 450.0


def job(name):
    case, mesh = name.split("_")
    if case == "LI0":
        return meier_job(name)
    pf, cmd, meta = B.job(f"C1_{mesh}")
    pf = os.path.join(OUT, name)
    cmd = [c if not c.startswith("plot_file=") else f"plot_file={pf}" for c in cmd]
    cmd = [c.replace(B.OUT, OUT) if c.startswith("plot_file=") else c for c in cmd]
    cmd += ["surface_patch.side_face_flux=1",
            f"surface_patch.side_face_factor={FACTOR[case]!r}"]
    meta = dict(meta, name=name, case=case, side_face_flux=1, side_face_factor=FACTOR[case])
    return pf, cmd, meta


def meier_job(name):
    """The D2i LI0 keys (Q09 + idle clock off) with the side-face key on, to
    450 s. Loaded through d2i_gate/run.py so the keys are exactly D2i's."""
    I = load("d2i_run", os.path.join(STUD, "d2i_gate", "run.py"))
    mesh = name.split("_")[1]
    pf_i, cmd, meta = I.job(f"LI0_{mesh}")
    pf = os.path.join(OUT, name)
    out = []
    for c in cmd:
        if c.startswith("plot_file="):
            out.append(f"plot_file={pf}")
        elif c.startswith("stop_time="):
            out.append(f"stop_time={MEIER_STOP!r}")
        else:
            out.append(c)
    out += ["surface_patch.side_face_flux=1", "surface_patch.side_face_factor=1.0"]
    meta = dict(meta, name=name, case="LI0", stop=MEIER_STOP,
                side_face_flux=1, side_face_factor=1.0)
    return pf, out, meta


def run_all(names, jobs, force):
    J.keep_awake()
    js = [job(n) for n in names]
    print(f"{len(js)} runs, {min(jobs, len(js))} in parallel -> {OUT}", flush=True)
    for _, _, m in js:
        print(f"  {m['name']}: f = {m['side_face_factor']}; dz {m['dz'] * 1e3:g} mm, "
              f"dt {m['dt'] * 1e3:g} ms, stop {m['stop']:g} s", flush=True)
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
        print(f"binary {J.binary_sha()[:8]}…")
        for n in CASES:
            _, cmd, m = job(n)
            print(f"{n:10s} f {m['side_face_factor']} dz {m['dz'] * 1e3:g} mm dt "
                  f"{m['dt'] * 1e3:g} ms stop {m['stop']:g} s")
            print("      side keys:", " ".join(c for c in cmd if "side_face" in c))
        return
    bad = [n for n in (a.cases or []) if n not in CASES]
    if not a.cases or bad:
        sys.exit(f"give --cases from {CASES}" + (f" (unknown: {bad})" if bad else ""))
    sys.exit(0 if run_all(a.cases, a.jobs, a.force) else 1)


if __name__ == "__main__":
    main()
