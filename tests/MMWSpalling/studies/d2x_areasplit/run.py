#!/usr/bin/env python3
"""D2x: close `A`'s wall band without losing its rate.

  run.py --cases C0 --jobs 1          the identity witness (200 s)
  run.py --cases C D E --jobs 3
  run.py --list
  run.py --kill NAME                  never pkill -f

`A` -- the confined map alone, no side-face path, no free parameter -- already
gives Meier's shape (2.27 mm at 150 mm descent) at 1.296 m/h. It is the CONTROL
THAT ALREADY PASSES; the bar is not "beat M0". What `A` still gets wrong is a
3-6 mm narrowing of the wall band at 100-125 mm depth, and the question is
whether heating that band closes it WITHOUT losing `A`'s rate, pit lead or
volume.

`f = 0.2` is the obstacle to answering cleanly: it multiplies an AREA inside
h*A*(T_gas - T), where an area factor and an h factor are algebraically
indistinguishable, so one number absorbs two unrelated errors and no output can
falsify it. And it serves two surfaces at once -- the cone flank (which makes
any side term explosive: M1 doubles the rate) and the wall band (which needs
more heat). This 2x2 separates them:

                 f = 0.2                 derived per-edge area
  all faces      M1  (on disk)           D
  r >= 38 mm     C                       E

with `A` (no side path) as the reference, also on disk. `C` tests the
drain-compensation reading; `D`/`E` test the derive-f reading. `A` and `M1` are
imported by path and NEVER re-run.

The map string is copied VERBATIM from d2w_attribution/output/A.log and asserted
byte-identical to M1.log's before launch (D2w takeaway 9): re-deriving it would
make these legs incomparable with `M1`.
"""
import argparse
import concurrent.futures as cf
import importlib.util
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)


def load(name, path):
    sp = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


D2W = load("d2w_run", os.path.join(STUD, "d2w_attribution", "run.py"))
D2V = D2W.D2V
R, E = D2W.R, D2W.E                 # d2c run_one (frozen watchdogs), d2e kill
A_OUT = D2W.OUT                     # d2w output: A, B
M_OUT = D2V.OUT                     # d2v output: M0, M1

OUT = os.path.join(HERE, "output")
for _m in (E, R):
    _m.OUT = OUT

R_MIN = 0.038          # first ring where A's collapse appears (-13.1 %)
F_FUDGE = 0.2          # frozen stand-in, never scanned
C0_STOP = 200.0        # identity witness only

# leg -> (extra keys beyond A's argv, stop_time override)
LEGS = {
    "C0": (["surface_patch.side_face_flux=1",
            f"surface_patch.side_face_factor={F_FUDGE!r}",
            "surface_patch.side_face_factor_above=0.0",
            "surface_patch.side_face_r_min=0.0",
            "surface_patch.side_face_area_exact=0"], C0_STOP),
    "C":  (["surface_patch.side_face_flux=1",
            f"surface_patch.side_face_factor={F_FUDGE!r}",
            "surface_patch.side_face_factor_above=0.0",
            f"surface_patch.side_face_r_min={R_MIN!r}"], None),
    "D":  (["surface_patch.side_face_flux=1",
            "surface_patch.side_face_factor_above=0.0",
            "surface_patch.side_face_area_exact=1"], None),
    "E":  (["surface_patch.side_face_flux=1",
            "surface_patch.side_face_factor_above=0.0",
            f"surface_patch.side_face_r_min={R_MIN!r}",
            "surface_patch.side_face_area_exact=1"], None),
}
CASES = list(LEGS)


def argv_of(leg, out):
    return open(os.path.join(out, leg + ".log")).readline().split()


def map_token():
    """A's h_expr, VERBATIM, asserted byte-identical to M1's."""
    ta = [c for c in argv_of("A", A_OUT) if c.startswith("surface_patch.h_expr=")]
    tm = [c for c in argv_of("M1", M_OUT) if c.startswith("surface_patch.h_expr=")]
    if len(ta) != 1 or len(tm) != 1:
        raise SystemExit("expected exactly one h_expr token in each of A.log / M1.log")
    if ta[0] != tm[0]:
        raise SystemExit("A's map differs from M1's -- the legs would be incomparable")
    return ta[0]


def base_argv():
    """A's resolved argv: the frozen d2i_gate key set plus the confined map."""
    cmd = list(argv_of("A", A_OUT))
    tok = map_token()
    hit = [n for n, c in enumerate(cmd) if c.startswith("surface_patch.h_expr=")]
    assert len(hit) == 1 and cmd[hit[0]] == tok, "A.log's map is not the asserted token"
    return cmd


def _set(cmd, key, val):
    pre = key + "="
    hit = [n for n, c in enumerate(cmd) if c.startswith(pre)]
    if len(hit) > 1:
        raise SystemExit(f"{key} appears {len(hit)} times")
    if hit:
        cmd[hit[0]] = pre + val
    else:
        cmd.append(pre + val)


def job(name):
    extra, stop = LEGS[name]
    cmd = base_argv()
    pf = os.path.join(OUT, name)
    _set(cmd, "plot_file", pf)
    if stop is not None:
        _set(cmd, "stop_time", repr(stop))
    cmd += list(extra)
    # D/E must not carry side_face_factor (the parser aborts, but fail here first)
    if name in ("D", "E") and any(c.startswith("surface_patch.side_face_factor=") for c in cmd):
        raise SystemExit(f"{name} must not set side_face_factor")
    if any(c.startswith(("surface_patch.jet_above_h", "surface_patch.jet_T_above")) for c in cmd):
        raise SystemExit(f"{name}: above-plane heating is out of scope")
    if not any(c == "surface_patch.side_face_factor_above=0.0" for c in cmd):
        raise SystemExit(f"{name} must set side_face_factor_above=0.0 explicitly")
    # Read the labels back off the RESOLVED argv, never off the recipe: an
    # `in` test against the key list silently compares whole elements, not
    # substrings, and would mislabel the .done metadata.
    def val(key, default=None, cast=float):
        hit = [c for c in cmd if c.startswith("surface_patch." + key + "=")]
        return cast(hit[0].split("=", 1)[1]) if hit else default
    meta = dict(name=name, dz=0.002, nxy=60, nz=200, lz=0.4, dt=0.016,
                stop=float([c for c in cmd if c.startswith("stop_time=")][0].split("=")[1]),
                r_min=val("side_face_r_min", 0.0),
                area_exact=val("side_face_area_exact", 0, int),
                f=val("side_face_factor", None),
                f_above=val("side_face_factor_above", None))
    return pf, cmd, meta


def keep_awake():
    try:
        return subprocess.Popen(["caffeinate", "-dis", "-w", str(os.getpid())],
                                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except (OSError, ValueError) as e:
        print(f"  (no caffeinate: {e})", flush=True)
        return None


def run_all(names, jobs, force):
    keep_awake()
    js = [job(n) for n in names]
    print(f"{len(js)} runs, {min(jobs, len(js))} in parallel -> {OUT}", flush=True)
    for _, _, m in js:
        print(f"  {m['name']:2s} r_min {m['r_min']*1e3:5.1f} mm, area_exact {m['area_exact']}, "
              f"f {m['f']}, stop {m['stop']:g} s", flush=True)
    ok = True
    with cf.ThreadPoolExecutor(max_workers=min(jobs, len(js))) as ex:
        for fu in cf.as_completed([ex.submit(R.run_one, pf, c, m, force) for pf, c, m in js]):
            m, status, wall = fu.result()
            print(f"  {m['name']:2s} {status} at t = {m.get('t_end', 0):g} s "
                  f"({wall / 60:.1f} min)", flush=True)
            ok = ok and not status.startswith("FAILED")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", nargs="*", default=CASES)
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--kill")
    a = ap.parse_args()
    if a.kill:
        return E.kill(a.kill)
    if a.list:
        b = base_argv()
        print(f"map token {len(map_token())} chars, verbatim from A.log, "
              f"asserted == M1.log's")
        for n in CASES:
            _, cmd, m = job(n)
            print(f"\n{n}: stop {m['stop']:g} s")
            for c in cmd:
                if c not in b:
                    print(f"    + {c[:100]}")
        return
    bad = [n for n in a.cases if n not in CASES]
    if bad:
        sys.exit(f"give --cases from {CASES} (unknown: {bad})")
    sys.exit(0 if run_all(a.cases, a.jobs, a.force) else 1)


if __name__ == "__main__":
    main()
