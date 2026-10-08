#!/usr/bin/env python3
"""D2ac: Meier at Meier's kinematics -- the prescribed feed.

  run.py --cases FM FL FH S0 --jobs 2      ~80 min
  run.py --cases FM13 --jobs 1             ~23 min
  run.py --cases SMOKE_FM SMOKE_FL --jobs 2   the pre-campaign smoke, ~2 min
  run.py --list
  run.py --kill NAME    binary path + plot_file match (never pkill -f)

Base is FROZEN `A` (d2w_attribution/output/A.log), lifted verbatim, plus D2ab's
closure block. NO SOURCE CHANGE in this packet.

⚠ THE RATE IS AN INPUT. `nozzle_feed` is set to Meier's kinematics; burner
descent is therefore NOT a result and is never scored. The deliverable is the
FLOOR-FLUX EXCESS -- how much faster the rock recedes than the feed demands.

⚠⚠ FORCED KEY DROPS. `nozzle_descent = prescribed` REFUSES A's feet block:
MMWSpalling.H:7593 aborts on `jet_free_surface = 1` and :7596 on
`foot_body_clearance`. Ten keys are dropped (DROP_FEET below); with the whole
block gone the unused-key list is empty. Measured with a parse probe, not
reasoned. ⇒ `FM` differs from `A` by MORE than the feed, and every comparison
says so. See PREFLIGHT.md (d).

⚠ `lz` 0.4 -> 0.5 m (the block's real height) so `amr.n_cell` z must go
200 -> 250: dz stays 2 mm AND cubic, which side_face_flux enforces.

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

ALOG = os.path.join(STUD, "d2w_attribution", "output", "A.log")
RBLOG = os.path.join(STUD, "d2z1_delete", "output", "RB_D.log")

# PREFLIGHT.md (a). DO NOT EDIT without re-hashing.
NOZZLE_Z0 = "0.550"
LZ, NZ = "0.5", "250"
FEED15 = "0.0004166666666666667"      # 1.5 m/h
FEED13 = "0.0003611111111111111"      # 1.3 m/h

# PREFLIGHT.md (d): forced by the code, not chosen.
DROP_FEET = ["surface_patch.jet_free_surface", "surface_patch.jet_fs_aspect",
             "surface_patch.jet_fs_exponent", "surface_patch.foot_r_inner",
             "surface_patch.foot_r_outer", "surface_patch.foot_standoff",
             "surface_patch.foot_rule", "surface_patch.foot_npads",
             "surface_patch.foot_pad_quantile",
             "surface_patch.foot_body_clearance"]

SHARED = ["surface_patch.nozzle_descent=prescribed",
          "surface_patch.nozzle_z0=" + NOZZLE_Z0,
          "surface_patch.side_face_flux=1",
          "surface_patch.side_face_step_height=1",
          "surface_patch.side_face_area_law=wetted",
          "surface_patch.side_face_factor_above=0.0",
          "surface_patch.step_diag_radius=0.050",
          # ⚠ PREFLIGHT_ADDENDUM.md. The guard's default is
          #     v.feet ? v.foot_r_inner : v.radius      (:7676)
          # so under PRESCRIBED it defaults to surface_patch.radius = 0.2 m and
          # every column in the domain is a collision column. All four legs
          # aborted at t = 120 s on unrecessed rock at r = 59 mm, which is the
          # NORMAL state of a burner inside its own hole (:7670 says so).
          # 0.028 is the code's OWN feet-mode default (foot_r_inner = the Ø56
          # bore radius), not a tuned value, and rc2 appears at exactly two
          # lines -- its definition and the abort -- so it moves no trajectory.
          "surface_patch.nozzle_collision_radius=0.028",
          # ⚠ multi-value keys are SEPARATE tokens, exactly as the lifted argv
          # writes them; the drop below is token-aware for the same reason.
          "amr.n_cell=60", "60", NZ,
          "geometry.prob_hi=0.12", "0.12", LZ]

# D2ab's derived annulus bracket, unchanged.
WALL = lambda h: ["surface_patch.side_face_h=" + h,
                  "surface_patch.side_face_T=live",
                  "surface_patch.side_face_h_r_min=0.040"]


def feed(v):
    return ["surface_patch.nozzle_feed=" + v]


# name -> (source log, expected input basename, extra keys, stop override)
CASES = {
    "FM":   (ALOG, "input_drilling", SHARED + feed(FEED15), None),
    "FL":   (ALOG, "input_drilling", SHARED + feed(FEED15) + WALL("38.4"), None),
    "FH":   (ALOG, "input_drilling", SHARED + feed(FEED15) + WALL("169.3"), None),
    "FM13": (ALOG, "input_drilling", SHARED + feed(FEED13), None),
    "S0":   (RBLOG, "input_drilling", [], None),
    # Smoke: same argv, 40 s, so the crater has real side faces and the
    # prescribed descent is visible in nozzle_z.
    "SMOKE_FM": (ALOG, "input_drilling", SHARED + feed(FEED15), "40.0"),
    "SMOKE_FL": (ALOG, "input_drilling", SHARED + feed(FEED15) + WALL("38.4"), "40.0"),
    "SMOKE_S0": (RBLOG, "input_drilling", [], "40.0"),
}


def frozen_argv(name):
    log, want_input, extra, stop = CASES[name]
    with open(log) as f:
        line = f.readline().strip()
    tok = shlex.split(line)
    binary, inp, keys = tok[6], tok[7], tok[8:]
    assert binary.endswith("mmwspalling-3d-g++"), binary
    assert inp.endswith(want_input), f"{name}: expected {want_input}, got {inp}"
    keys = [k for k in keys if not k.startswith("plot_file=")]
    joined = " ".join(keys)
    for k in ("nozzle_z0=", "nozzle_feed=", "side_face_h=", "side_face_T="):
        assert k not in joined, f"{name}: source leg already sets {k}"
    # ⚠ "has extras" is NOT the same as "is a prescribed leg": the identity
    # legs carry a stop_time override and nothing else, and they MUST keep
    # their feet block. Key on the marker itself.
    presc = any(e == "surface_patch.nozzle_descent=prescribed" for e in extra)
    drop = set()
    if presc:
        assert "nozzle_descent=feet" in joined, f"{name}: base is not a feet run"
        assert "jet_closure=enthalpy" in joined, f"{name}: base has no march"
        # The feet block CANNOT survive the prescribed path -- see PREFLIGHT (d).
        drop |= set(DROP_FEET) | {"surface_patch.nozzle_descent"}
    drop |= {e.split("=", 1)[0] for e in extra if "=" in e}
    if stop:
        drop.add("stop_time")
        extra = extra + ["stop_time=" + stop]
    # ⚠ TOKEN-AWARE DROP. ParmParse multi-value keys arrive as `k=v0 v1 v2`,
    # i.e. one `=` token followed by bare ones. Dropping only the `=` token
    # would leave `60 200` dangling in the argv and silently re-bind them to
    # whatever key came before. Drop the key AND its trailing bare tokens.
    kept, i = [], 0
    while i < len(keys):
        k = keys[i]
        if "=" in k and k.split("=", 1)[0] in drop:
            i += 1
            while i < len(keys) and "=" not in keys[i]:
                i += 1
            continue
        kept.append(k)
        i += 1
    if presc:
        # The drop list must have removed every feet-only key, or the run dies
        # on ParmParse's strict unused-key check partway through init.
        left = [k for k in kept if k.split("=", 1)[0] in set(DROP_FEET)]
        assert not left, f"{name}: feet-only keys survived: {left}"
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
         "nozzle_descent": g("nozzle_descent", "feet"),
         "nozzle_z0": g("nozzle_z0", "(unset)"),
         "nozzle_feed": g("nozzle_feed", "0"),
         "side_face_h": g("side_face_h", "(inherit)"),
         "side_face_T": g("side_face_T", "(inherit)"),
         "side_face_h_r_min": g("side_face_h_r_min", "(unset)")}
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
            print(f"{n:9s} descent={m['nozzle_descent']:<10s} "
                  f"z0={m['nozzle_z0']:>6s} feed={m['nozzle_feed']:>22s} "
                  f"h={m['side_face_h']:>9s} lz={m.get('lz','')} "
                  f"stop={m.get('stop','')}")
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
