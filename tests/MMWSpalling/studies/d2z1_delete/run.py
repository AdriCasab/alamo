#!/usr/bin/env python3
"""D2z-1: whole-cell removal unconditional, `spall.removal_increment` deleted.

  run.py --cases W2 W1 --jobs 2        arbiter at both meshes      ~18 min
  run.py --cases REBASE --jobs 2       Stage C anchors             ~15 min
  run.py --list
  run.py --kill NAME                   binary path + plot_file match (never pkill -f)

EVERY leg's argv is LIFTED VERBATIM from a frozen log and adds NOTHING. The key
is gone, so a whole-cell run is now just the frozen argv against the new binary.
That is what makes criterion (c0) exact: W2/W1 must be byte-identical to D2z-0's
F2/F1, which ran the same physics with `spall.removal_increment=full_cell` set.

  W2  d2l_scan/A_f02.log       -> the 2 mm arbiter, and ALSO Stage C's
  W1  d2l_scan/A_f02_1mm.log   -> the 1 mm arbiter, and ALSO Stage C's

⚠ W2/W1 ARE the re-baselined d2l_scan pair -- the same argv under whole-cell --
so Stage C does not re-run them. REBASE is only the other two anchors:

  RB_C1_2mm  d2j0b_plane/C1_2mm.log   the f = 0 cascade reference
  RB_D       d2x_areasplit/D.log      Z0's byte-identity target, inherited by D2z

Frozen studies are READ-ONLY: every output lands under this study's own
output/, and nothing is written into d2l_scan, d2j0b_plane or d2x_areasplit. A
later packet points at these re-baselined outputs by path.
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

# name -> (frozen log, expected input basename)
CASES = {
    "W2":        (os.path.join(STUD, "d2l_scan", "output", "A_f02.log"), "input_hotdisk"),
    "W1":        (os.path.join(STUD, "d2l_scan", "output", "A_f02_1mm.log"), "input_hotdisk"),
    "RB_C1_2mm": (os.path.join(STUD, "d2j0b_plane", "output", "C1_2mm.log"), "input_hotdisk"),
    "RB_D":      (os.path.join(STUD, "d2x_areasplit", "output", "D.log"), "input_drilling"),
}
GROUPS = {"REBASE": ["RB_C1_2mm", "RB_D"]}


def frozen_argv(name):
    """The resolved command line of the frozen leg, minus mpirun/binary/input."""
    log, want_input = CASES[name]
    with open(log) as f:
        line = f.readline().strip()
    tok = shlex.split(line)
    binary, inp, keys = tok[6], tok[7], tok[8:]
    assert binary.endswith("mmwspalling-3d-g++"), binary
    assert inp.endswith(want_input), f"{name}: expected {want_input}, got {inp}"
    keys = [k for k in keys if not k.startswith("plot_file=")]
    joined = " ".join(keys)
    # The key is DELETED -- a frozen leg that set it would no longer parse.
    assert "spall.removal_increment" not in joined, f"{name}: frozen leg sets the deleted key"
    return inp, keys


def meta_of_argv(name, rc=0, wall=0.0):
    """Geometry read back off this leg's own resolved argv -- never re-typed.
    .done must be a dict literal: d2j0_hotdisk/run.py:meta_of literal_evals it."""
    _, keys = frozen_argv(name)
    kv, i = {}, 0
    while i < len(keys):
        if "=" in keys[i]:
            k, v = keys[i].split("=", 1)
            vals, j = [v], i + 1
            while j < len(keys) and "=" not in keys[j]:
                vals.append(keys[j]); j += 1
            kv[k] = vals; i = j
        else:
            i += 1
    m = {"name": name, "rc": rc, "seconds": round(wall, 1),
         "status": "ok" if rc == 0 else "FAILED", "removal_increment": "whole_cell (no key)"}
    if "amr.n_cell" in kv and "geometry.prob_hi" in kv:
        nxy = int(kv["amr.n_cell"][0]); lxy = float(kv["geometry.prob_hi"][0])
        m.update({"nxy": nxy, "nz": int(kv["amr.n_cell"][2]), "dz": lxy / nxy,
                  "lxy": lxy, "lz": float(kv["geometry.prob_hi"][2])})
    for k in ("timestep", "stop_time"):
        if k in kv:
            m[{"timestep": "dt", "stop_time": "stop"}[k]] = float(kv[k][0])
    return m


def job(name):
    inp, keys = frozen_argv(name)
    return (["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4",
             BIN, inp, f"plot_file={os.path.join(OUT, name)}"] + keys)


def run_one(name):
    os.makedirs(OUT, exist_ok=True)
    cmd = job(name)
    t0 = time.time()
    with open(os.path.join(OUT, name + ".log"), "w") as fh:
        fh.write(" ".join(shlex.quote(c) for c in cmd) + "\n")
        fh.flush()
        rc = subprocess.call(cmd, cwd=ROOT, stdout=fh, stderr=subprocess.STDOUT)
    dt = time.time() - t0
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
            print(f"{n:10s} {' '.join(job(n)[7:])}")
        return 0

    if a.kill:
        # Match on BOTH the binary path and this study's plot_file, never
        # `pkill -f` (which kills the shell waiters too).
        # ⚠ macOS pgrep ignores -a and prints BARE PIDs, so reading the command
        # line off pgrep's output silently matches nothing and the kill reports
        # "killed 0" while the run carries on. Ask `ps` for the command line.
        pf = os.path.join(OUT, a.kill)
        pids = subprocess.run(["pgrep", "-f", "mmwspalling-3d"],
                              capture_output=True, text=True).stdout.split()
        n = 0
        for pid in pids:
            cmd = subprocess.run(["ps", "-o", "command=", "-p", pid],
                                 capture_output=True, text=True).stdout
            if f"plot_file={pf}" in cmd:
                os.kill(int(pid), 15); n += 1
        print(f"killed {n} process(es) for {a.kill}")
        return 0

    if not a.cases:
        ap.error("--cases is required")
    cases = []
    for n in a.cases:
        cases.extend(GROUPS.get(n, [n]))
    for n in cases:
        if n not in CASES:
            ap.error(f"unknown case {n}; have {list(CASES)} and groups {list(GROUPS)}")

    # Runs are not checkpointed: a host sleep loses the leg outright.
    subprocess.Popen(["caffeinate", "-dis", "-w", str(os.getpid())])

    print(f"launching {cases} with {a.jobs} job(s) at {time.strftime('%T')}")
    with cf.ThreadPoolExecutor(max_workers=a.jobs) as ex:
        rcs = list(ex.map(run_one, cases))
    print(f"all done at {time.strftime('%T')}")
    return 0 if all(r == 0 for r in rcs) else 1


if __name__ == "__main__":
    sys.exit(main())
