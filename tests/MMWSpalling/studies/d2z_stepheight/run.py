#!/usr/bin/env python3
"""D2z: the step-height closure on the arbiter, then Meier.

  run.py --cases S0 S2 S2w S2f --jobs 4     identity + the 2 mm legs   ~8 min
  run.py --cases S1 S1w --jobs 2            the GATE at 1 mm           ~30 min
  run.py --cases SW --jobs 1                Meier, ONLY if the gate passes
  run.py --list
  run.py --kill NAME    binary path + plot_file match (never pkill -f)

Every leg's argv is LIFTED VERBATIM from a frozen/re-baselined log and adds only
the keys in the table below. ⚠ The reference pair is the WHOLE-CELL
d2z1_delete/{W2,W1}; the frozen d2l_scan/A_f02 pair is HALF-CELL and no longer
the reference. Nothing is written into any imported study.

  S2  / S1   W2 / W1 argv + step key, area_law = step      the closure
  S2w / S1w  W2 / W1 argv + step key, area_law = wetted    the other heating law
  S2f        W2 argv + flux, factor = 0.0, above = 0.0     f = 0 physics WITH
             (NO step key)                                 the drain registered
  S0         RB_D argv, step key ABSENT                    identity vs RB_D
  SW         A's argv + the shared block                   Meier

⚠ macOS pgrep ignores -a and prints bare PIDs, so --kill reads command lines
from `ps -o command= -p` instead; a kill that reports success without killing is
worse than one that errors.
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

Z1 = os.path.join(STUD, "d2z1_delete", "output")
AW = os.path.join(STUD, "d2w_attribution", "output", "A.log")

SHARED = ["surface_patch.side_face_flux=1",
          "surface_patch.side_face_step_height=1",
          "surface_patch.side_face_factor_above=0.0"]

# name -> (source log, expected input basename, extra keys)
CASES = {
    "S2":  (os.path.join(Z1, "W2.log"), "input_hotdisk",
            SHARED + ["surface_patch.side_face_area_law=step"]),
    "S1":  (os.path.join(Z1, "W1.log"), "input_hotdisk",
            SHARED + ["surface_patch.side_face_area_law=step"]),
    "S2w": (os.path.join(Z1, "W2.log"), "input_hotdisk",
            SHARED + ["surface_patch.side_face_area_law=wetted"]),
    "S1w": (os.path.join(Z1, "W1.log"), "input_hotdisk",
            SHARED + ["surface_patch.side_face_area_law=wetted"]),
    # f = 0 physics, but with side_face_factor_above set so the drain columns
    # register. Also the FIRST test of "key-on at f = 0 equals key-off".
    "S2f": (os.path.join(Z1, "W2.log"), "input_hotdisk",
            ["surface_patch.side_face_flux=1",
             "surface_patch.side_face_factor=0.0",
             "surface_patch.side_face_factor_above=0.0"]),
    "S0":  (os.path.join(Z1, "RB_D.log"), "input_drilling", []),
    "SW":  (AW, "input_drilling", SHARED),
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
    for k in ("side_face_step_height", "side_face_area_law"):
        assert k not in joined, f"{name}: source leg already sets {k}"
    # ⚠ The source arbiter legs (W2/W1) carry side_face_factor=0.2 -- they ARE
    # the f = 0.2 configuration. The step key REPLACES f, and the parser refuses
    # both, so f is stripped here rather than silently shadowed. Same for the
    # D2x derived weight. This is the one deliberate edit to a lifted argv and
    # it is exactly what "no free parameter" means.
    drop = set()
    if any("side_face_step_height" in e for e in extra):
        assert "side_face_area_exact" not in joined, f"{name}: source sets area_exact"
        drop.add("surface_patch.side_face_factor")
    # Anything `extra` sets explicitly overrides the source rather than being
    # appended twice, so --list shows the argv the binary actually receives.
    for e in extra:
        if "=" in e:
            drop.add(e.split("=", 1)[0])
    kept = [k for k in keys
            if not ("=" in k and k.split("=", 1)[0] in drop)]
    return inp, kept + extra


def meta_of_argv(name, rc=0, wall=0.0):
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
    m = {"name": name, "rc": rc, "seconds": round(wall, 1),
         "status": "ok" if rc == 0 else "FAILED",
         "area_law": kv.get("surface_patch.side_face_area_law", ["(none)"])[0],
         "step_height": kv.get("surface_patch.side_face_step_height", ["0"])[0]}
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
            print(f"{n:4s} {' '.join(job(n)[7:])}")
        return 0

    if a.kill:
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
