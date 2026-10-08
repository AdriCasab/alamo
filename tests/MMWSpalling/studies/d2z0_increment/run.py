#!/usr/bin/env python3
"""D2z-0 Stage B: spall.removal_increment = full_cell on the frozen arbiter pair.

  run.py --cases H2 F2 --jobs 2      2 mm identity + full_cell   ~8 min
  run.py --cases F1    --jobs 1      1 mm full_cell              ~28 min
  run.py --list
  run.py --kill NAME                 binary path + plot_file match (never pkill -f)

EVERY leg's argv is LIFTED VERBATIM from the frozen d2l_scan logs
(`output/A_f02.log`, `output/A_f02_1mm.log`) and differs from it by AT MOST the
single key `spall.removal_increment`. The frozen studies are read-only: nothing
here writes into d2l_scan or d2j0b_plane.

  H2  A_f02's argv, key ABSENT     -> must be byte-identical to A_f02 (crit d)
  F2  A_f02's argv + full_cell     -> 2 mm scored leg
  F1  A_f02_1mm's argv + full_cell -> 1 mm scored leg

WHY H2 EXISTS. The 2639-file tree gate cannot stand in for it: no input in that
set turns spalling on with an event log at the arbiter's configuration, so the
frozen A_f02 baseline is not covered by it. H2 re-validates that baseline
against the current binary AND is criterion (d)'s identity leg. It costs 4 min.

THE CONFIGURATION IS THE POINT. input_hotdisk runs
`surface_patch.mode = convective_flame` with
`spall.stress_source = local_thermoelastic`, so sp_weibull_presc_T_face is
FALSE and the depth scan's surface region collapses to T_face = T_top. That is
the branch full_cell acts on; a prescribed_T leg would not reach it.
"""
import argparse
import concurrent.futures as cf
import os
import shlex
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)
ROOT = os.path.abspath(os.path.join(STUD, "..", "..", ".."))
OUT = os.path.join(HERE, "output")
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")

# Frozen sources. READ ONLY.
FROZEN = {
    "2mm": os.path.join(STUD, "d2l_scan", "output", "A_f02.log"),
    "1mm": os.path.join(STUD, "d2l_scan", "output", "A_f02_1mm.log"),
}

CASES = {
    #  name: (mesh, extra keys)
    "H2": ("2mm", []),
    "F2": ("2mm", ["spall.removal_increment=full_cell"]),
    "F1": ("1mm", ["spall.removal_increment=full_cell"]),
}


def frozen_argv(mesh):
    """The resolved command line of the frozen leg, minus mpirun/binary/input.

    Lifted from the log's first line so no key is re-typed. Asserts the
    invariants this packet depends on rather than trusting the file."""
    with open(FROZEN[mesh]) as f:
        line = f.readline().strip()
    tok = shlex.split(line)
    # mpirun --oversubscribe --bind-to none -np 4 <bin> <input> <keys...>
    binary = tok[6]
    inp = tok[7]
    keys = tok[8:]
    assert binary.endswith("mmwspalling-3d-g++"), binary
    assert inp.endswith("input_hotdisk"), inp
    keys = [k for k in keys if not k.startswith("plot_file=")]
    # The invariants criterion (a) and the mechanism rely on.
    joined = " ".join(keys)
    assert "surface_patch.side_face_factor=0.2" in joined, "not the f = 0.2 leg"
    assert "surface_patch.side_face_flux=1" in joined
    assert "stop_time=250.0" in joined
    assert "spall.removal_increment" not in joined, "frozen leg already sets the key"
    assert ("timestep=0.016" if mesh == "2mm" else "timestep=0.008") in joined
    return inp, keys


def meta_of_argv(name, rc=0, wall=0.0):
    """Geometry read back off this leg's own resolved argv -- never re-typed."""
    mesh, extra = CASES[name]
    _, keys = frozen_argv(mesh)
    kv, i = {}, 0
    toks = keys + extra
    while i < len(toks):
        t = toks[i]
        if "=" in t:
            k, v = t.split("=", 1)
            vals = [v]
            j = i + 1
            while j < len(toks) and "=" not in toks[j]:
                vals.append(toks[j]); j += 1
            kv[k] = vals; i = j
        else:
            i += 1
    nxy = int(kv["amr.n_cell"][0])
    nz = int(kv["amr.n_cell"][2])
    lxy = float(kv["geometry.prob_hi"][0])
    return {"name": name, "case": "C1", "mesh": mesh,
            "dz": lxy / nxy, "nxy": nxy, "nz": nz, "lxy": lxy,
            "lz": float(kv["geometry.prob_hi"][2]),
            "dt": float(kv["timestep"][0]), "stop": float(kv["stop_time"][0]),
            "side_face_flux": 1, "side_face_factor": 0.2,
            "removal_increment": "full_cell" if extra else "half_cell (unset)",
            "rc": rc, "seconds": round(wall, 1),
            "status": "ok" if rc == 0 else "FAILED"}


def job(name):
    mesh, extra = CASES[name]
    inp, keys = frozen_argv(mesh)
    pf = os.path.join(OUT, name)
    return (["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4",
             BIN, inp, f"plot_file={pf}"] + keys + extra)


def run_one(name):
    os.makedirs(OUT, exist_ok=True)
    cmd = job(name)
    log = os.path.join(OUT, name + ".log")
    t0 = time.time()
    with open(log, "w") as fh:
        fh.write(" ".join(shlex.quote(c) for c in cmd) + "\n")
        fh.flush()
        rc = subprocess.call(cmd, cwd=ROOT, stdout=fh, stderr=subprocess.STDOUT)
    dt = time.time() - t0
    # The .done must be a Python dict literal: d2j0_hotdisk/run.py:meta_of
    # ast.literal_eval()s it, and analyze.py imports that reader UNCHANGED so
    # the geometry comes from the same place the frozen legs' did.
    with open(os.path.join(OUT, name + ".done"), "w") as fh:
        fh.write(repr(meta_of_argv(name, rc, dt)))
    print(f"  {name}: rc={rc} in {dt/60:.1f} min", flush=True)
    return rc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", nargs="+")
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--kill")
    a = ap.parse_args()

    if a.list:
        for n in CASES:
            print(f"{n:4s} {' '.join(job(n)[7:])}")
        return 0

    if a.kill:
        # Match on BOTH the binary path and this study's plot_file, never
        # `pkill -f` (which kills the shell waiters too).
        pf = os.path.join(OUT, a.kill)
        out = subprocess.run(["pgrep", "-af", "mmwspalling-3d-g\\+\\+"],
                             capture_output=True, text=True).stdout
        n = 0
        for ln in out.splitlines():
            if f"plot_file={pf}" in ln:
                pid = int(ln.split()[0])
                os.kill(pid, 15)
                n += 1
        print(f"killed {n} process(es) for {a.kill}")
        return 0

    if not a.cases:
        ap.error("--cases is required")
    for n in a.cases:
        if n not in CASES:
            ap.error(f"unknown case {n}; have {list(CASES)}")

    # Runs are not checkpointed: a host sleep loses the leg outright. Assert the
    # hold here so it cannot be forgotten at the call site.
    subprocess.Popen(["caffeinate", "-dis", "-w", str(os.getpid())])

    print(f"launching {a.cases} with {a.jobs} job(s) at {time.strftime('%T')}")
    with cf.ThreadPoolExecutor(max_workers=a.jobs) as ex:
        rcs = list(ex.map(run_one, a.cases))
    print(f"all done at {time.strftime('%T')}")
    return 0 if all(r == 0 for r in rcs) else 1


if __name__ == "__main__":
    sys.exit(main())
