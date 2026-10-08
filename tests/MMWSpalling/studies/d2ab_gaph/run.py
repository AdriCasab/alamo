#!/usr/bin/env python3
"""D2ab: the gap-derived wall coefficient, and the face assignment.

  run.py --cases GL GH WLv WLr SWc S0 --jobs 2   ~71 min
  run.py --cases SMOKE_GL SMOKE_S0 --jobs 2    the pre-campaign smoke, ~1 min
  run.py --list
  run.py --kill NAME    binary path + plot_file match (never pkill -f)

Base is FROZEN `SW` (d2z2_multicell/output/SW.log), lifted verbatim, one edge per
leg. `S0` is FROZEN `RB_D` (d2z1_delete/output/RB_D.log) with no new key at all.

⚠ h IS NOT CHOSEN. 38 and 169 W/m^2K are the two ends of the annulus bracket
derived in `preflight.py` / `PREFLIGHT.md` (hashed BEFORE this file ran): the
Meier-MEASURED Ø85-93 hole against the Ø80 skirt, the nozzle mass flow, and the
laminar/Dittus-Boelter pair at Re = 2035-2134. No point value inside the bracket
is taken; both ends are run. A rate outside Meier's band IS the result.

⚠ `side_face_T = live` is the packet's BLOCKER FIX: D2aa fed the wall 1463.5 K
lifted from SW while WL/WL81's own exhaust ran at 1601/1631 K, so their driving
dT was 22-27 % too small. Every D2aa wall number is provisional until this runs.

⚠ THE ATTRIBUTION PAIR is `WL81`(frozen, const T, all faces) -> `WLv`(live T,
all faces) -> `WLr`(live T, wall only). ONE edge each, or the frozen-T fix and
the face assignment are confounded and neither is attributable.

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

SWLOG = os.path.join(STUD, "d2z2_multicell", "output", "SW.log")
RBLOG = os.path.join(STUD, "d2z1_delete", "output", "RB_D.log")

# preflight.py section 3. DO NOT EDIT without re-hashing PREFLIGHT.md.
H_LOW, H_HIGH = "38.4", "169.3"
R_PADS = "0.040"          # == foot_r_outer, the burner radius


def wall(h, rmin):
    return ["surface_patch.side_face_h=" + h,
            "surface_patch.side_face_T=live",
            "surface_patch.side_face_h_r_min=" + rmin]


# name -> (source log, expected input basename, extra keys, stop_time override)
CASES = {
    "GL":  (SWLOG, "input_drilling", wall(H_LOW,  R_PADS), None),
    "GH":  (SWLOG, "input_drilling", wall(H_HIGH, R_PADS), None),
    "WLr": (SWLOG, "input_drilling", wall("81.0", R_PADS), None),
    "WLv": (SWLOG, "input_drilling", wall("81.0", "0.0"),  None),
    "S0":  (RBLOG, "input_drilling", [],                   None),
    # `SWc` -- the BRIDGE LEG. SW's argv verbatim, no new key. It does three
    # jobs the packet needs and has no other source for:
    #   1. /implement step zero: re-validate a frozen baseline whose keys the
    #      identity set does NOT exercise. `S0`/`RB_D` never touches
    #      side_face_flux, so nothing else proves SW still reproduces.
    #   2. criterion (c) needs an INHERITED-MAP leg carrying side_P_{bore,pads,
    #      wall}. Frozen SW predates those columns (checked: 0 of 3 present), so
    #      (c) is literally unscoreable without this leg.
    #   3. it discharges D2aa's open item-0(a) debt, which was exactly "a bridge
    #      leg should have been run".
    "SWc": (SWLOG, "input_drilling", [],                   None),
    # S0R: S0 re-run against the POST-campaign binary. The live-without-march
    # parse guard was deliberately added AFTER the sweep (never rebuild
    # mid-sweep), so this leg is what PROVES that guard is parse-only rather
    # than asserting it: S0R must be byte-identical to S0.
    "S0R": (RBLOG, "input_drilling", [],                   None),
    # Smoke: the SAME argv, 40 s, so the crater has real side faces. The witness (side_T_wall) must be
    # non-zero and tracking jet_T_exhaust on SMOKE_GL, and the column must be
    # ABSENT on SMOKE_S0. A grep proves reachability in principle; this proves
    # it happened.
    "SMOKE_GL": (SWLOG, "input_drilling", wall(H_LOW, R_PADS), "40.0"),
    "SMOKE_WV": (SWLOG, "input_drilling", wall(H_LOW, "0.0"),  "40.0"),
    "SMOKE_S0": (RBLOG, "input_drilling", [],                  "40.0"),
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
    # The source leg must not already carry this packet's keys.
    for k in ("side_face_h=", "side_face_T=", "side_face_h_r_min"):
        assert k not in joined, f"{name}: source leg already sets {k}"
    if extra:
        # The radius gate is foot_r_outer, which is query_required ONLY under
        # the feet rule -- so assert the base actually runs it, or the gate
        # would be comparing against a number the base never read.
        assert "nozzle_descent=feet" in joined, f"{name}: base is not a feet run"
        assert "surface_patch.foot_r_outer=0.040" in joined, \
            f"{name}: base foot_r_outer is not 0.040, so {R_PADS} is not the pad edge"
        # `live` needs the march: without jet_closure = enthalpy there is no
        # jet_T_exhaust and the key would silently sit on its nozzle seed.
        assert "jet_closure=enthalpy" in joined, f"{name}: base has no march"
    drop = {e.split("=", 1)[0] for e in extra if "=" in e}
    if stop:
        drop.add("stop_time")
        extra = extra + ["stop_time=" + stop]
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
         "side_face_h_r_min": g("side_face_h_r_min", "(unset)"),
         "n_new_keys": sum(1 for k in kv if k.split(".")[-1] in
                           ("side_face_h", "side_face_T", "side_face_h_r_min"))}
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
            print(f"{n:9s} h={m['side_face_h']:>9s} T={m['side_face_T']:>9s} "
                  f"r_min={m['side_face_h_r_min']:>7s} stop={m.get('stop', '')}")
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
