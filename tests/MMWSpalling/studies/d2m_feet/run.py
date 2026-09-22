#!/usr/bin/env python3
"""D2m: is the cut diameter set by the feet or by the flame? Key-only, no build.

  run.py --cases R40 R32 R44 R48 S3648 --jobs 5     ~20 min
  run.py --identity      R40 vs the D2l f = 0 control, every shared thermo row
  run.py --list
  run.py --kill NAME     binary path + plot_file match (never pkill -f)

Binary: bin/mmwspalling-3d-g++ = 62dea451... (the D2k campaign build). This
packet changes no source and builds nothing.

Base configuration: the D2i `LI0` Meier keys at 2 mm (D2f Q09 support rule,
pads, quantile 0.9, clearance off, idle clock off), `side_face_flux` OFF
(D2l closed that line: f_min 0.5 > f_rate 0), stop 250 s. Loaded from
studies/d2i_gate/run.py by path; that study is frozen and is not edited.

Cases -- only the foot annulus changes
--------------------------------------
  R40     foot_r_inner 0.028, foot_r_outer 0.040   CONTROL (unchanged)
  R32     0.028, 0.032
  R44     0.028, 0.044
  R48     0.028, 0.048
  S3648   0.036, 0.048   the ring SHIFTED OUTWARD at constant width, to
                         separate "the outer edge sets it" from "the whole
                         ring sets it"

Two couplings the packet does not mention, both found by reading the source
before launch. Declared in PREDICTIONS.md §1 before anything was run.
--------------------------------------------------------------------------
1. `foot_r_inner` SILENTLY sets `nozzle_collision_radius`
   (MMWSpalling.H:5695, `feet ? foot_r_inner : radius`), and the D2i keys never
   set it explicitly. Without care, S3648 would also move the collision guard
   28 -> 36 mm. **Every case therefore pins
   `surface_patch.nozzle_collision_radius = 0.028`**, which is the control's own
   effective value -- so the control is unchanged (the identity guard proves
   it) and the coupling is removed from the ladder.
2. `foot_r_outer` ALSO defines the free-surface boundary
   (MMWSpalling.H:2062: a column is free-surface when r > foot_r_outer and it
   lies inside the aspect cone measured from r - foot_r_outer). There is no
   separate key for that boundary, so this one cannot be pinned away. Its size
   is BOUNDED from D2l instead of by extra runs: `P_fsoff0` removed the free
   surface **entirely** and moved min Ø by +0.3 mm and ROP by -2.3 %. Moving
   foot_r_outer only modulates it partially, so the effect on min Ø is well
   under the 64 -> 96 mm swing H-feet predicts -- but it is material at the
   few-per-cent level for ROP, so `jet_fs_cols` is reported per case.
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


B = load("d2j0b_run", os.path.join(STUD, "d2j0b_plane", "run.py"))   # for J: the runner
J = B.J                                   # d2j0: run_one / kill / keep_awake / binary_sha
I = load("d2i_run", os.path.join(STUD, "d2i_gate", "run.py"))        # the LI0 Meier keys

OUT = os.path.join(HERE, "output")
B.OUT = OUT
J.OUT = OUT                               # run_one's makedirs and kill()'s token
I.OUT = OUT
I.E.OUT = OUT
I.R.OUT = OUT

# The D2l f = 0 control. ACTIVE_STEP calls it "D2l's M_f0", but no such run
# exists: D2l's f = 0 row was read from d2i_gate/output/LI0_2mm (side key off),
# which is the run R40 must reproduce. Named explicitly so the substitution is
# not silent.
CONTROL_REF = os.path.join(STUD, "d2i_gate", "output", "LI0_2mm")

STOP = 250.0
COLLISION_R = 0.028                       # pinned in every case; see the docstring
RING = {                                  # name -> (foot_r_inner, foot_r_outer)
    "R40":   (0.028, 0.040),
    "R32":   (0.028, 0.032),
    "R44":   (0.028, 0.044),
    "R48":   (0.028, 0.048),
    "S3648": (0.036, 0.048),
}
CASES = list(RING)


def job(name):
    ri, ro = RING[name]
    _, cmd, meta = I.job("LI0_2mm")
    pf = os.path.join(OUT, name)
    out = []
    for c in cmd:
        if c.startswith("plot_file="):
            out.append(f"plot_file={pf}")
        elif c.startswith("stop_time="):
            out.append(f"stop_time={STOP!r}")
        elif c.startswith("surface_patch.foot_r_inner="):
            out.append(f"surface_patch.foot_r_inner={ri!r}")
        elif c.startswith("surface_patch.foot_r_outer="):
            out.append(f"surface_patch.foot_r_outer={ro!r}")
        else:
            out.append(c)
    # the base keys must actually have contained both, or the replacement above
    # silently did nothing and every case would be the control
    for k, v in (("foot_r_inner", ri), ("foot_r_outer", ro)):
        hit = [c for c in out if c.startswith(f"surface_patch.{k}=")]
        if len(hit) != 1 or hit[0] != f"surface_patch.{k}={v!r}":
            sys.exit(f"{name}: failed to set surface_patch.{k} (got {hit})")
    if any(c.startswith("surface_patch.nozzle_collision_radius=") for c in out):
        sys.exit(f"{name}: base keys already set nozzle_collision_radius; "
                 "the pin below would be ambiguous")
    if any("side_face_flux" in c for c in out):
        sys.exit(f"{name}: side_face_flux must stay off in this packet")
    out.append(f"surface_patch.nozzle_collision_radius={COLLISION_R!r}")
    meta = dict(meta, name=name, case=name, stop=STOP, foot_r_inner=ri, foot_r_outer=ro,
                ring_mm=(ri * 1e3, ro * 1e3), collision_r=COLLISION_R, side_face_flux=0)
    return pf, out, meta


def run_one(pf, cmd, meta, force):
    p = pf + "_jet_profile.csv"
    if force and os.path.exists(p):
        os.remove(p)                      # the Meier runs append to it
    return J.run_one(pf, cmd, meta, force)


def run_all(names, jobs, force):
    J.keep_awake()
    js = [job(n) for n in names]
    print(f"binary {J.binary_sha()[:8]}…  {len(js)} runs, {min(jobs, len(js))} in parallel "
          f"-> {OUT}", flush=True)
    for _, _, m in js:
        print(f"  {m['name']:6s} ring {m['ring_mm'][0]:.0f}–{m['ring_mm'][1]:.0f} mm "
              f"(Ø {2 * m['ring_mm'][1]:.0f} stance); dz {m['dz'] * 1e3:g} mm, "
              f"dt {m['dt'] * 1e3:g} ms, stop {m['stop']:g} s", flush=True)
    ok = True
    with cf.ThreadPoolExecutor(max_workers=min(jobs, len(js))) as ex:
        for fu in cf.as_completed([ex.submit(run_one, pf, c, m, force) for pf, c, m in js]):
            m, status, wall = fu.result()
            print(f"  {m['name']:6s} {status} at t = {m.get('t_end', 0):g} s "
                  f"({wall / 60:.1f} min)", flush=True)
            ok = ok and not status.startswith("FAILED")
    return ok


def rows_upto(path, t):
    out = []
    for i, line in enumerate(open(path)):
        if i == 0 or not line.strip():
            out.append(line)
            continue
        if float(line.split()[0]) > t:
            break
        out.append(line)
    return out


def identity():
    """R40 must reproduce the D2l f = 0 control byte for byte. If it does not,
    the harness is reaching through the frozen studies incorrectly (or the
    pinned collision radius is not inert) and nothing else here is meaningful."""
    a = os.path.join(OUT, "R40", "thermo.dat")
    b = os.path.join(CONTROL_REF, "thermo.dat")
    if not (os.path.exists(a) and os.path.exists(b)):
        print(f"R40 or the control is missing ({a} / {b})")
        return False
    ra, rb = rows_upto(a, STOP), rows_upto(b, STOP)
    same = ra == rb
    print(f"| run | vs | rows (t <= {STOP:g} s) | verdict |\n|---|---|---|---|")
    print(f"| R40 | d2i_gate/LI0_2mm (D2l's f = 0 control) | {len(ra) - 1} | "
          f"{'byte-identical' if same else 'DIFFER'} |")
    print(f"\nIdentity: {'PASS' if same else 'FAIL'}")
    return same


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", nargs="*")
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--identity", action="store_true")
    ap.add_argument("--kill")
    a = ap.parse_args()
    if a.kill:
        return J.kill(a.kill)
    if a.identity:
        sys.exit(0 if identity() else 1)
    if a.list:
        print(f"binary {J.binary_sha()[:8]}…   control ref {CONTROL_REF}")
        for n in CASES:
            _, cmd, m = job(n)
            print(f"{n:6s} ring {m['ring_mm'][0]:.0f}–{m['ring_mm'][1]:.0f} mm  "
                  f"stance Ø {2 * m['ring_mm'][1]:.0f} mm  stop {m['stop']:g} s")
            print("      foot:", " ".join(c for c in cmd if "foot_r" in c),
                  "|", " ".join(c for c in cmd if "collision" in c))
        return
    bad = [n for n in (a.cases or []) if n not in CASES]
    if not a.cases or bad:
        sys.exit(f"give --cases from {CASES}" + (f" (unknown: {bad})" if bad else ""))
    sys.exit(0 if run_all(a.cases, a.jobs, a.force) else 1)


if __name__ == "__main__":
    main()
