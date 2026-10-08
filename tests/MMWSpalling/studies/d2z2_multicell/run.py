#!/usr/bin/env python3
"""D2z-2: repair the three measurement defects, then settle the multi-cell branch.

  run.py --cases R2 R2w R2d S0 --jobs 2      the repairs + identity   ~12 min
  run.py --cases R1d --jobs 1                THE RIM TEST at 1 mm     ~21 min
  run.py --cases SW SWs SW0 --jobs 1         Meier                    ~36 min
  run.py --list
  run.py --kill NAME     binary path + plot_file match (never pkill -f)

Every leg's argv is LIFTED VERBATIM from a frozen/re-baselined log and adds only
the literal key lines in ACTIVE_STEP.md. Nothing is written into any imported
study.

⚠ TWO BINARIES, and which one a leg used is the difference under test.
`StepConductance`'s normaliser is a SOURCE change, not a key -- the packet
forbids leaving a second form alive behind a switch -- so the comparison legs
were run from a binary built with the superseded D2z form and then the source
was restored. BOTH binaries are fingerprinted into the leg's .done record, and
the harness REFUSES to launch if the binary on disk is not the one the leg
declares. That is the only thing standing between this and an unattributable
result, so it is an assertion, not a comment.

  D2z form  (normalised to n = 1 regardless of n) : R2, R2w, SW0
  derived   ((1 - e^{-delta a})/(1 - e^{-n a}))   : R2d, R1d, SW, SWs, S0

⚠ macOS pgrep ignores -a and prints bare PIDs, so --kill reads command lines
from `ps -o command= -p` instead.
"""
import argparse
import concurrent.futures as cf
import hashlib
import os
import shlex
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)
ROOT = os.path.abspath(os.path.join(STUD, "..", "..", ".."))
OUT = os.path.join(HERE, "output")

BIN_DERIVED = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
# Built by temporarily restoring the D2z normaliser, then reverting. Kept out of
# bin/ (a managed build directory) on purpose.
BIN_D2Z = os.environ.get(
    "D2Z_FORM_BIN",
    "/private/tmp/claude-501/-Users-tzetze20-amr-tools-alamo/"
    "ac0a0f5a-4af6-48bf-b61d-9009cf0c8681/scratchpad/mmwspalling-3d-g++-d2zform")

SHA_DERIVED = "2dde94489142ce84597102a171ffcfde1f88d41324783250c49437a97bda75b6"
SHA_D2Z     = "f9085c65ca28959aa650d2248a0f9b2277afedd93bf56032ab2846e42c15ef5b"

Z1 = os.path.join(STUD, "d2z1_delete", "output")
AW = os.path.join(STUD, "d2w_attribution", "output", "A.log")

# The packet's literal key lines, verbatim.
SHARED_ARB = ["surface_patch.side_face_flux=1",
              "surface_patch.side_face_step_height=1",
              "surface_patch.side_face_factor_above=0.0",
              "surface_patch.step_diag_radius=0.025"]
SHARED_MEI = ["surface_patch.side_face_flux=1",
              "surface_patch.side_face_step_height=1",
              "surface_patch.side_face_factor_above=0.0",
              "surface_patch.step_diag_radius=0.050"]

# ⚠ Deliverable 5: side_face_area_law DEFAULTS to `step`, so every Meier leg
# states its law as a literal key line. SW as the D2z harness wrote it would
# have silently run the law the tie rule did NOT pick.
def arb(law):
    return SHARED_ARB + [f"surface_patch.side_face_area_law={law}"]


def mei(law):
    return SHARED_MEI + [f"surface_patch.side_face_area_law={law}"]


# name -> (source log, expected input basename, extra keys, binary, sha)
CASES = {
    "R2":  (os.path.join(Z1, "W2.log"), "input_hotdisk", arb("step"),   BIN_D2Z,     SHA_D2Z),
    "R2w": (os.path.join(Z1, "W2.log"), "input_hotdisk", arb("wetted"), BIN_D2Z,     SHA_D2Z),
    "R2d": (os.path.join(Z1, "W2.log"), "input_hotdisk", arb("step"),   BIN_DERIVED, SHA_DERIVED),
    "R1d": (os.path.join(Z1, "W1.log"), "input_hotdisk", arb("step"),   BIN_DERIVED, SHA_DERIVED),
    "SW":  (AW, "input_drilling", mei("wetted"), BIN_DERIVED, SHA_DERIVED),
    "SWs": (AW, "input_drilling", mei("step"),   BIN_DERIVED, SHA_DERIVED),
    "SW0": (AW, "input_drilling", mei("wetted"), BIN_D2Z,     SHA_D2Z),
    "S0":  (os.path.join(Z1, "RB_D.log"), "input_drilling", [], BIN_DERIVED, SHA_DERIVED),
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def frozen_argv(name):
    log, want_input, extra, _, _ = CASES[name]
    with open(log) as f:
        line = f.readline().strip()
    tok = shlex.split(line)
    binary, inp, keys = tok[6], tok[7], tok[8:]
    assert binary.endswith("mmwspalling-3d-g++"), binary
    assert inp.endswith(want_input), f"{name}: expected {want_input}, got {inp}"
    keys = [k for k in keys if not k.startswith("plot_file=")]
    joined = " ".join(keys)
    for k in ("side_face_step_height", "side_face_area_law", "step_diag_radius"):
        assert k not in joined, f"{name}: source leg already sets {k}"
    # ⚠ W2/W1 carry side_face_factor=0.2 -- they ARE the f = 0.2 configuration.
    # The step key REPLACES f and the parser refuses both, so f is stripped
    # rather than silently shadowed. The one deliberate edit to a lifted argv,
    # and it is exactly what "no free parameter" means.
    drop = set()
    if any("side_face_step_height" in e for e in extra):
        assert "side_face_area_exact" not in joined, f"{name}: source sets area_exact"
        drop.add("surface_patch.side_face_factor")
    for e in extra:
        if "=" in e:
            drop.add(e.split("=", 1)[0])
    kept = [k for k in keys if not ("=" in k and k.split("=", 1)[0] in drop)]
    return inp, kept + extra


def meta_of_argv(name, rc=0, wall=0.0, sha=""):
    """.done must be a dict literal: d2j0_hotdisk/run.py:meta_of literal_evals it."""
    _, keys = frozen_argv(name)
    _, _, _, binpath, want = CASES[name]
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
         # THE attribution record: which conductance form produced this leg.
         "conductance": "derived" if want == SHA_DERIVED else "d2z_form",
         "binary": os.path.basename(binpath), "binary_sha256": sha,
         "area_law": kv.get("surface_patch.side_face_area_law", ["(none)"])[0],
         "step_height": kv.get("surface_patch.side_face_step_height", ["0"])[0],
         "step_diag_radius": kv.get("surface_patch.step_diag_radius", ["0"])[0]}
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
    binpath = CASES[name][3]
    return (["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4",
             binpath, inp, f"plot_file={os.path.join(OUT, name)}"] + keys)


def run_one(name):
    os.makedirs(OUT, exist_ok=True)
    binpath, want = CASES[name][3], CASES[name][4]
    got = sha256(binpath)
    # The whole packet is a comparison between two conductance forms. A leg run
    # against the wrong binary is not a wrong number, it is an unattributable
    # one, so this refuses rather than warns.
    if got != want:
        print(f"  {name}: ABORT -- {os.path.basename(binpath)} is {got[:16]}..., "
              f"leg declares {want[:16]}...", flush=True)
        return 1
    cmd = job(name)
    t0 = time.time()
    with open(os.path.join(OUT, name + ".log"), "w") as fh:
        fh.write(" ".join(shlex.quote(c) for c in cmd) + "\n")
        fh.flush()
        rc = subprocess.call(cmd, cwd=ROOT, stdout=fh, stderr=subprocess.STDOUT)
    dt = time.time() - t0
    with open(os.path.join(OUT, name + ".done"), "w") as fh:
        fh.write(repr(meta_of_argv(name, rc, dt, got)))
    print(f"  {name}: rc={rc} in {dt/60:.1f} min  [{meta_of_argv(name)['conductance']}]",
          flush=True)
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
            print(f"{n:4s} [{m['conductance']:9s}] {' '.join(job(n)[7:])}")
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
