#!/usr/bin/env python3
"""D2aa: heat delivery under the burner body.

  run.py --cases WL WL81 PD BO S0 --jobs 2     ~76 min
  run.py --list
  run.py --kill NAME    binary path + plot_file match (never pkill -f)

Base is FROZEN `SW` (d2z2_multicell/output/SW.log), one edge per leg. `SW`'s own
argv is lifted verbatim; each leg adds only its literal key lines.

⚠ `WL` runs at BOTH ends of the von Rohr band and neither end is chosen: 142 is
the annulus derivation's value (standing constraint, never varied) and is the
PRIMARY; 81 is the Dittus-Boelter end. Picking 81 alone would land P3 and would
be choosing h to hit a rate.
⚠ `side_face_T = 1463.5` is `jet_T_exhaust` MEASURED on SW over 150-400 s
(1463.46), not chosen.
⚠ `PD` is BOUNDED: the slot width is on record nowhere, so pad_shield_open = 0
is one end of the bracket and frozen `SW` (open = 1) is the other.

⚠ macOS pgrep ignores -a and prints bare PIDs, so --kill reads command lines from
`ps -o command= -p`.
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

Z2 = os.path.join(STUD, "d2z2_multicell", "output")
SWLOG = os.path.join(Z2, "SW.log")
RBLOG = os.path.join(STUD, "d2z1_delete", "output", "RB_D.log")

T_EXH = "1463.5"          # measured jet_T_exhaust on SW, 150-400 s
WALL = ["surface_patch.side_face_h={h}", f"surface_patch.side_face_T={T_EXH}"]
PADS = ["surface_patch.pad_shield_open=0.0"]


def wall(h):
    return [k.format(h=h) for k in WALL]


# name -> (source log, expected input basename, extra keys)
CASES = {
    "WL":   (SWLOG, "input_drilling", wall("142.0")),
    "WL81": (SWLOG, "input_drilling", wall("81.0")),
    "PD":   (SWLOG, "input_drilling", PADS),
    "BO":   (SWLOG, "input_drilling", wall("142.0") + PADS),
    "S0":   (RBLOG, "input_drilling", []),
}


def frozen_argv(name):
    log, want_input, extra = CASES[name]
    with open(log) as f:
        line = f.readline().strip()
    tok = shlex.split(line)
    binary, inp, keys = tok[6], tok[7], tok[8:]
    assert binary.endswith("mmwspalling-3d-g++"), binary
    assert inp.endswith(want_input), f"{name}: expected {want_input}, got {inp}"
    keys = [k for k in keys if not k.startswith("plot_file=")]
    joined = " ".join(keys)
    # The source leg must not already carry this packet's keys.
    for k in ("side_face_h=", "side_face_T=", "pad_shield_open"):
        assert k not in joined, f"{name}: source leg already sets {k}"
    # The pad shield and the region split both key off the foot annulus, which
    # is query_required ONLY under the feet rule -- so assert the base runs it.
    if extra:
        assert "nozzle_descent=feet" in joined, f"{name}: base is not a feet run"
        assert "foot_r_inner" in joined and "foot_r_outer" in joined, \
            f"{name}: base does not set the foot annulus"
    # Anything `extra` sets explicitly overrides rather than appending twice.
    drop = {e.split("=", 1)[0] for e in extra if "=" in e}
    kept = [k for k in keys if not ("=" in k and k.split("=", 1)[0] in drop)]
    return inp, kept + extra


def meta_of_argv(name, rc=0, wall_s=0.0):
    """.done must be a dict literal: d2j0_hotdisk/run.py:meta_of literal_evals it."""
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
    g = lambda k, d: kv.get("surface_patch." + k, [d])[0]
    m = {"name": name, "rc": rc, "seconds": round(wall_s, 1),
         "status": "ok" if rc == 0 else "FAILED",
         "side_face_h": g("side_face_h", "(inherit)"),
         "side_face_T": g("side_face_T", "(inherit)"),
         "pad_shield_open": g("pad_shield_open", "(off)"),
         "n_new_keys": sum(1 for k in kv if k.split(".")[-1] in
                           ("side_face_h", "side_face_T", "pad_shield_open"))}
    if "amr.n_cell" in kv and "geometry.prob_hi" in kv:
        nxy = int(kv["amr.n_cell"][0]); lxy = float(kv["geometry.prob_hi"][0])
        m.update({"nxy": nxy, "nz": int(kv["amr.n_cell"][2]), "dz": lxy / nxy,
                  "lxy": lxy, "lz": float(kv["geometry.prob_hi"][2])})
    for k, nm in (("timestep", "dt"), ("stop_time", "stop")):
        if k in kv:
            m[nm] = float(kv[k][0])
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
            m = meta_of_argv(n)
            print(f"{n:5s} h={m['side_face_h']:>9s} T={m['side_face_T']:>9s} "
                  f"pads={m['pad_shield_open']:>6s}")
        return 0

    if a.kill:
        pf = os.path.join(OUT, a.kill)
        pids = subprocess.run(["pgrep", "-f", "mmwspalling-3d"],
                              capture_output=True, text=True).stdout.split()
        n = 0
        for pid in pids:
            c = subprocess.run(["ps", "-o", "command=", "-p", pid],
                               capture_output=True, text=True).stdout
            if f"plot_file={pf}" in c:
                os.kill(int(pid), 15); n += 1
        print(f"killed {n} process(es) for {a.kill}")
        return 0

    if not a.cases:
        ap.error("--cases is required")
    for n in a.cases:
        if n not in CASES:
            ap.error(f"unknown case {n}; have {list(CASES)}")

    # Runs are not checkpointed: a host sleep loses the leg outright.
    subprocess.Popen(["caffeinate", "-dis", "-w", str(os.getpid())])

    print(f"launching {a.cases} with {a.jobs} job(s) at {time.strftime('%T')}")
    with cf.ThreadPoolExecutor(max_workers=a.jobs) as ex:
        rcs = list(ex.map(run_one, a.cases))
    print(f"all done at {time.strftime('%T')}")
    return 0 if all(r == 0 for r in rcs) else 1


if __name__ == "__main__":
    sys.exit(main())
