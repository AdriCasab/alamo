#!/usr/bin/env python3
"""D2l Gate 0: why the model drills 2.7x too fast. Key-only, no build.

  run.py --cases A_f10 A_f05 A_f02 A_f01 --jobs 4     0a, 2 mm arbiter, ~4 min
  run.py --cases A_f05_1mm A_f02_1mm --jobs 2         0a, the boundary 1 mm legs
  run.py --cases M_f10 M_f05 M_f02 M_f01 --jobs 4     0b, 2 mm Meier to 250 s
  run.py --cases P_fsoff P_flat --jobs 2              0c, the competing explanation
  run.py --identity        A_f10 vs d2k C1_2mm, M_f10 vs d2k LI0_2mm (t <= 250 s)
  run.py --list
  run.py --kill NAME       binary path + plot_file match (never pkill -f)

Binary: bin/mmwspalling-3d-g++ = 62dea451... (the D2k campaign build). Gate 0
changes no source and builds nothing.

Loads by path, never edits: studies/d2j0b_plane/run.py (the C1 arbiter keys)
and studies/d2i_gate/run.py (the LI0 Meier keys, which pull d2h -> d2f -> d2e
-> d2c). Both are frozen; their OUT is redirected here at import so nothing is
written into them.

The scan dial
-------------
`surface_patch.side_face_factor` (f) scales q on exposed lateral faces
linearly, at exactly the place an h scale would. That is why it maps the
window. **Its physical meaning is still AREA, not h, and no value found here
may be kept as a calibrated setting** (ACTIVE_STEP guardrail). f = 0 is not
run: it is the side key off, i.e. D2i, whose runs already exist.

Cases
-----
  A_f*     the D2j-0b C1 cascade (static disk, descending exclusion plane),
           2 mm unless suffixed _1mm, stop 250 s. f = 1.0 is a re-run of the
           d2k arbiter leg and must come out byte-identical (--identity).
  M_f*     the D2i LI0 Meier keys (Q09 support rule, idle clock off) with the
           side key on, 2 mm, stop 250 s. f = 1.0 must reproduce the first
           250 s of d2k's LI0_2mm byte-for-byte (--identity).
  P_fsoff  M_f10 with surface_patch.jet_free_surface removed (with its aspect
           and exponent keys, which the parser rejects when the flag is off).
           The wellhead finding: phase 1 was confined exhaust, not 293 K air.
  P_flat   M_f10 with the walljet far-field exponent n = 1.0 -> 0.5.
  P_fsoff0 / P_flat0
           the same two probes with the side key OFF (f = 0, i.e. the D2i
           configuration) to 250 s.

Why 0c is run at BOTH ends of the dial
--------------------------------------
ACTIVE_STEP puts 0c at f = 1. Measured on the existing d2k LI0_2mm run,
**both probes are inert there**, for reasons that are properties of the
resulting geometry, not of the keys:

  * the far law is pow(12D/max(12D, s), n) with 12D = 90 mm, and with side
    faces on the hole is self-similar at jet_s_c = 67-69 mm and
    patch_min_standoff = 63 mm, so every column sits in the clamp region
    where the factor is exactly 1 for any n  -> P_flat is byte-identical;
  * jet_fs_cols = 0.0 through the whole 150-250 s window and m_ratio = 1.000,
    so the free-surface dilution is already switched off by the hole depth
    -> P_fsoff differs only through the early transient.

With the side key off (D2i) both are live: jet_s_c runs 88 -> 166 mm, well
past the 90 mm clamp, and jet_fs_cols = 17.9 in the window. That is also the
configuration the competing explanation is ABOUT -- the runaway centre pit
with a skirt. So the probes are run at f = 1 as the packet specifies (and to
settle the inertness by byte-identity rather than by argument) AND at f = 0,
where they can actually move something. Both are pre-registered in
PREDICTIONS.md before launch. No new key is used at either end.
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


B = load("d2j0b_run", os.path.join(STUD, "d2j0b_plane", "run.py"))   # arbiter keys
J = B.J                                                              # d2j0: run_one/kill/keep_awake
I = load("d2i_run", os.path.join(STUD, "d2i_gate", "run.py"))        # Meier LI0 keys

OUT = os.path.join(HERE, "output")
B.OUT = OUT
J.OUT = OUT                      # run_one's makedirs and kill()'s token
I.OUT = OUT
I.E.OUT = OUT
I.R.OUT = OUT

D2K = os.path.join(STUD, "d2k_sideface", "output")                   # read-only, for --identity

# the pre-registered ladder; f = 0 is the side key off (D2i), already run
F = {"f10": 1.0, "f05": 0.5, "f02": 0.2, "f01": 0.1}
MEIER_STOP = 250.0
FS_KEYS = ("surface_patch.jet_free_surface", "surface_patch.jet_fs_aspect",
           "surface_patch.jet_fs_exponent")
FAR_FROM, FAR_TO = "pow(0.09/max(0.09,s),1.0)", "pow(0.09/max(0.09,s),0.5)"

ARBITER = [f"A_{k}" for k in F] + [f"A_{k}_1mm" for k in F]
MEIER = [f"M_{k}" for k in F]
PROBES = ["P_fsoff", "P_flat", "P_fsoff0", "P_flat0"]
CASES = ARBITER + MEIER + PROBES


def side_keys(f):
    """f = 0 is the side key OFF, not the key on at zero: that is the D2i
    configuration exactly, and it stays byte-comparable with D2i's runs."""
    if f == 0.0:
        return []
    return ["surface_patch.side_face_flux=1", f"surface_patch.side_face_factor={float(f)!r}"]


def arbiter_job(name):
    """D2j-0b C1 at 2 mm (or 1 mm), side-face flux on at f."""
    parts = name.split("_")
    f = F[parts[1]]
    mesh = parts[2] if len(parts) > 2 else "2mm"
    _, cmd, meta = B.job(f"C1_{mesh}")
    pf = os.path.join(OUT, name)
    cmd = [f"plot_file={pf}" if c.startswith("plot_file=") else c for c in cmd]
    cmd += side_keys(f)
    meta = dict(meta, name=name, gate="0a", side_face_flux=1, side_face_factor=f)
    return pf, cmd, meta


def meier_job(name, f, probe=None):
    """The D2i LI0 keys verbatim, side key on at f, stop 250 s."""
    _, cmd, meta = I.job("LI0_2mm")
    pf = os.path.join(OUT, name)
    out = []
    for c in cmd:
        if c.startswith("plot_file="):
            out.append(f"plot_file={pf}")
        elif c.startswith("stop_time="):
            out.append(f"stop_time={MEIER_STOP!r}")
        elif probe == "fsoff" and c.split("=")[0] in FS_KEYS:
            continue                       # the parser rejects aspect/exponent with the flag off
        elif probe == "flat" and c.startswith("surface_patch.h_expr="):
            assert FAR_FROM in c, f"far law {FAR_FROM} not found in h_expr"
            out.append(c.replace(FAR_FROM, FAR_TO))
        else:
            out.append(c)
    out += side_keys(f)
    meta = dict(meta, name=name, gate="0b" if probe is None else "0c", stop=MEIER_STOP,
                side_face_flux=1, side_face_factor=f, probe=probe or "none")
    return pf, out, meta


def job(name):
    if name.startswith("A_"):
        return arbiter_job(name)
    if name.startswith("M_"):
        return meier_job(name, F[name.split("_")[1]])
    if name.startswith("P_"):
        probe = "fsoff" if "fsoff" in name else "flat"
        return meier_job(name, 0.0 if name.endswith("0") else 1.0, probe=probe)
    raise KeyError(name)


def run_one(pf, cmd, meta, force):
    """J.run_one, plus the jet profile CSV in the list it clears on --force
    (the Meier runs append to it, so a stale one would corrupt the rerun)."""
    p = pf + "_jet_profile.csv"
    if force and os.path.exists(p):
        os.remove(p)
    return J.run_one(pf, cmd, meta, force)


def run_all(names, jobs, force):
    J.keep_awake()
    js = [job(n) for n in names]
    print(f"binary {J.binary_sha()[:8]}…  {len(js)} runs, {min(jobs, len(js))} in parallel "
          f"-> {OUT}", flush=True)
    for _, _, m in js:
        print(f"  {m['name']:10s} gate {m['gate']} f = {m['side_face_factor']}; "
              f"dz {m['dz'] * 1e3:g} mm, dt {m['dt'] * 1e3:g} ms, stop {m['stop']:g} s"
              + (f", probe {m['probe']}" if m.get("probe", "none") != "none" else ""), flush=True)
    ok = True
    with cf.ThreadPoolExecutor(max_workers=min(jobs, len(js))) as ex:
        for fu in cf.as_completed([ex.submit(run_one, pf, c, m, force) for pf, c, m in js]):
            m, status, wall = fu.result()
            print(f"  {m['name']:10s} {status} at t = {m.get('t_end', 0):g} s "
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
    """The f = 1.0 legs re-run d2k's keys through this harness; they must come
    out byte-identical. This is the guard that the D2l harness has not drifted
    from D2k while reaching through two frozen studies."""
    checks = [("A_f10", os.path.join(D2K, "C1_2mm"), 1e9),
              ("M_f10", os.path.join(D2K, "LI0_2mm"), MEIER_STOP)]
    ok = True
    L = ["| run | vs d2k | rows | verdict |", "|---|---|---|---|"]
    for name, ref, t in checks:
        a = os.path.join(OUT, name, "thermo.dat")
        b = os.path.join(ref, "thermo.dat")
        if not (os.path.exists(a) and os.path.exists(b)):
            L.append(f"| {name} | {os.path.basename(ref)} | — | NOT RUN |")
            ok = False
            continue
        ra, rb = rows_upto(a, t), rows_upto(b, t)
        same = ra == rb
        ok = ok and same
        L.append(f"| {name} | {os.path.basename(ref)} | {len(ra) - 1} (t <= {t:g} s) | "
                 f"{'byte-identical' if same else 'DIFFER'} |")
    L.append(f"\nIdentity: {'PASS' if ok else 'FAIL'}")
    print("\n".join(L))
    return ok


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
        print(f"binary {J.binary_sha()[:8]}…")
        for n in CASES:
            _, cmd, m = job(n)
            print(f"{n:10s} gate {m['gate']} f {m['side_face_factor']:<4g} dz {m['dz'] * 1e3:g} mm "
                  f"dt {m['dt'] * 1e3:g} ms stop {m['stop']:g} s probe {m.get('probe', 'none')}")
            print("      side:", " ".join(c for c in cmd if "side_face" in c))
            fs = [c for c in cmd if "jet_fs" in c or "jet_free_surface" in c]
            print("      fs  :", " ".join(fs) if fs else "(none — free surface off)")
            print("      far :", FAR_TO if any(FAR_TO in c for c in cmd) else FAR_FROM)
        return
    bad = [n for n in (a.cases or []) if n not in CASES]
    if not a.cases or bad:
        sys.exit(f"give --cases from {CASES}" + (f" (unknown: {bad})" if bad else ""))
    sys.exit(0 if run_all(a.cases, a.jobs, a.force) else 1)


if __name__ == "__main__":
    main()
