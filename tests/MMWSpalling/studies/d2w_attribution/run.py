#!/usr/bin/env python3
"""D2w: complete the 2x2 — is the near-cylinder the MAP or the SIDE FACES?

  run.py --cases A B --jobs 2
  run.py --list
  run.py --kill NAME                never pkill -f

D2v's `M1` changed TWO things at once: the confined h map AND the below-plane
side-face path at f = 0.2.  Its headline "M1 - M0 is the funnel fix" is
therefore confounded, and both mechanisms are large: the side faces add +59 %
to the top-face power (624.8 W on 1063.8 W of tops at t = 400 s), and burner
descent goes 1.50 -> 3.12 m/h, i.e. from INSIDE Meier's 1.3-1.6 band to twice it.

The factorial:

    cell   confined map   side-face path
    M0         no              no          on disk (== frozen d2i_gate/LI0_2mm)
    A         YES              no          NEW
    B          no             YES          NEW
    M1        YES             YES          on disk (D2v)

M0 and M1 are imported by path and NEVER re-run.

Both new legs start from the frozen d2i_gate key set that M0 uses.

  A  adds surface_patch.h_expr ONLY, copied VERBATIM from M1's resolved argv
     (output/M1.log), not re-derived -- re-deriving it would make A
     incomparable with M1 (N5).  A sets NO side key at all: setting
     side_face_factor = 0 without side_face_flux ABORTS (:6099), so the
     correct "off" cell is the key being ABSENT, not zeroed.

  B  adds the three side keys and leaves h_expr at M0's.
     side_face_factor_above = 0.0 is MANDATORY: unset, the sentinel falls back
     to side_face_factor = 0.2, and with h = 0 against 293 K above the plane the
     Newton solve returns q_in < 0, so B would COOL its above-plane faces
     instead of leaving them alone (N4 -- the same trap D2v documented for M1).

No src/ change.
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


D2V = load("d2v_run", os.path.join(STUD, "d2v_hmap", "run.py"))
R = D2V.R                       # d2c: run_one + the frozen stall/bottom watchdogs
E = D2V.E                       # d2e: kill, rows_upto
D2V_OUT = D2V.OUT               # frozen M0 / M1 live here; never written to

OUT = os.path.join(HERE, "output")
for _m in (E, R):
    _m.OUT = OUT

SIDE_KEYS = ("surface_patch.side_face_flux=1",
             "surface_patch.side_face_factor=0.2",
             "surface_patch.side_face_factor_above=0.0")
LEGS = ("A", "B")
FROZEN = ("M0", "M1")


def argv_of(leg):
    """The resolved argv of a frozen D2v leg, exactly as it was run."""
    return open(os.path.join(D2V_OUT, leg + ".log")).readline().split()


def map_token():
    """M1's h_expr token, VERBATIM. Copy, never rebuild (N5)."""
    t = [c for c in argv_of("M1") if c.startswith("surface_patch.h_expr=")]
    if len(t) != 1:
        raise SystemExit(f"expected exactly 1 h_expr token in M1.log, found {len(t)}")
    return t[0]


def base_argv():
    """M0's argv = the frozen d2i_gate key set at stop_time 500."""
    return list(argv_of("M0"))


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
    cmd = base_argv()
    pf = os.path.join(OUT, name)
    _set(cmd, "plot_file", pf)
    if name == "A":
        tok = map_token()
        hit = [n for n, c in enumerate(cmd) if c.startswith("surface_patch.h_expr=")]
        assert len(hit) == 1, "M0 must carry exactly one h_expr"
        cmd[hit[0]] = tok                       # verbatim substitution
    elif name == "B":
        cmd += list(SIDE_KEYS)
    else:
        raise SystemExit(f"unknown leg {name}")
    # Invariants the packet's inert table names.
    n_side = sum(1 for c in cmd if c.startswith("surface_patch.side_face"))
    if name == "A" and n_side:
        raise SystemExit("N1: leg A must carry NO side key")
    if name == "B" and n_side != 3:
        raise SystemExit(f"leg B must carry 3 side keys, has {n_side}")
    h = [c for c in cmd if c.startswith("surface_patch.h_expr=")][0]
    if name == "B" and h != [c for c in base_argv()
                             if c.startswith("surface_patch.h_expr=")][0]:
        raise SystemExit("N2: leg B must keep M0's h_expr")
    if name == "A" and h != map_token():
        raise SystemExit("N5: leg A's map is not M1's verbatim")
    meta = dict(name=name, dz=0.002, nxy=60, nz=200, lz=0.4, dt=0.016, stop=500.0,
                confined_map=(name == "A"), side_path=(name == "B"),
                h_expr=h.split("=", 1)[1], n_side_keys=n_side)
    return pf, cmd, meta


def keep_awake():
    """These runs are not checkpointed; a host sleep loses one outright."""
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
        print(f"  {m['name']}: confined map {m['confined_map']}, side path "
              f"{m['side_path']} ({m['n_side_keys']} side keys), stop {m['stop']:g} s",
              flush=True)
    ok = True
    with cf.ThreadPoolExecutor(max_workers=min(jobs, len(js))) as ex:
        for fu in cf.as_completed([ex.submit(R.run_one, pf, c, m, force) for pf, c, m in js]):
            m, status, wall = fu.result()
            print(f"  {m['name']} {status} at t = {m.get('t_end', 0):g} s "
                  f"({wall / 60:.1f} min)", flush=True)
            ok = ok and not status.startswith("FAILED")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", nargs="*", default=list(LEGS))
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--kill")
    a = ap.parse_args()
    if a.kill:
        return E.kill(a.kill)
    if a.list:
        b = base_argv()
        for n in LEGS:
            _, cmd, m = job(n)
            extra = [c for c in cmd if c not in b]
            drop = [c for c in b if c not in cmd]
            print(f"\n{n}: +{len(extra)} keys, -{len(drop)}")
            for c in extra:
                print(f"    + {c[:110]}")
            for c in drop:
                print(f"    - {c[:110]}")
        print(f"\nmap token length {len(map_token())} chars, taken verbatim from "
              f"{os.path.join(D2V_OUT, 'M1.log')}")
        return
    bad = [n for n in a.cases if n not in LEGS]
    if bad:
        sys.exit(f"give --cases from {list(LEGS)} (unknown: {bad})")
    sys.exit(0 if run_all(a.cases, a.jobs, a.force) else 1)


if __name__ == "__main__":
    main()
