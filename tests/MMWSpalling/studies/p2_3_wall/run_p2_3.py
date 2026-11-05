#!/usr/bin/env python3
"""P2-3 -- the floor coefficient to Kant, the wall on, the gas marched up the bore.

Legs, all on the P2-2f CTRL8 argv with the floor peak scaled to Kant and the
wall switched on (`side_face_factor_above = 1`):
  W1    march ON,  side_face_h_scale 1.0
  W074  march ON,  side_face_h_scale 0.743 (the floor's own factor)
  NC    march OFF  -- the NEGATIVE CONTROL: uniform floor-exit gas on every wall
        face, the D2v M2 configuration that ballooned the mouth to Ø 137 mm.

⚠⚠ THE RADIUS GATE, AND WHY LOWERING IT IS SAFE **NOW**. The bore wall sits at
r ~ 39.7 mm (measured `wall_r_bore` 0.0397679 on m13) and the kernel's gate
`side_face_h_r_min` was 0.040 -- the wall was 0.23 mm, about a tenth of a cell,
INSIDE it. Two things were gated on that one number: the kernel's choice of the
explicit (h, T) pair, and the wall law's own HEIGHT DATUM. Lowering it alone
would have shifted the datum across the whole band (P2-1: the datum moves the
answer 1.04-1.74x) and confounded exactly what this step measures.
⇒ `side_face_h_datum = bore` first decouples them: the datum now truncates at
  the per-step MEASURED bore radius and no longer reads the key at all.
⇒ `side_face_h_r_min = 0.036` then admits the bore wall to the explicit pair
  WITHOUT touching the datum. Both halves are recorded, per deliverable 4.
⚠ It also moves the 7a/7b Newton split for the 36-40 mm annulus, which is the
  CORRECT classification -- those faces are hole wall, not floor riser -- and is
  local to that annulus rather than global like the datum.

⚠ Carried argv must ROUND-TRIP (P2-2f rev 1 lost 100 mm of rock to a parser).
Every base key is passed through verbatim, one line = one argv entry, so
multi-value keys (`amr.n_cell=60 60 300`) and the quoted `h_expr` stay intact.
The only edit is the ONE peak literal, and `argv.sh` records what actually ran.

⚠ Waits on a SENTINEL FILE, never on a `pgrep` of its own command line -- that
pattern matches the waiting shell itself and cannot exit (P2-2f takeaway 12).

Usage: run_p2_3.py short [legs...] | full [legs...] | reduce
"""
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
DECK = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "meier",
                    "sp_meier_pilot", "input_drilling")
ARGV_FILE = os.path.join(ROOT, "tests", "MMWSpalling", "studies",
                         "p2_2c_meier", "feet_argv.txt")
OUT = os.path.join(HERE, "legs")

# the scaled peak, from derive_scale.py -- NOT retyped. 1449.0 * 0.742499.
NEW_PEAK = "1075.9"

COMMON = [
    "surface_patch.side_face_factor_above=1.0",   # THE WALL ON
    "surface_patch.side_face_h_datum=bore",       # datum off the typed radius
    "surface_patch.side_face_h_r_min=0.036",      # gate admits the bore wall
    "weibull.m=13",
]
LEGS = {
    "W1":   ["surface_patch.jet_wall_march=1", "surface_patch.side_face_h_scale=1.0"],
    "W074": ["surface_patch.jet_wall_march=1", "surface_patch.side_face_h_scale=0.743"],
    "NC":   [],                                   # march OFF
}


def base_keys():
    """The CTRL8 argv, with the floor peak scaled. One line = one argv entry."""
    keys = [l.strip() for l in open(ARGV_FILE) if l.strip()]
    out, n_sub = [], 0
    for k in keys:
        if "1449.0*" in k:
            k2, n = re.subn(re.escape("1449.0*"), NEW_PEAK + "*", k, count=1)
            n_sub += n
            out.append(k2)
        else:
            out.append(k)
    if n_sub != 1:
        raise SystemExit("expected exactly ONE `1449.0*` in the argv, found %d "
                         "-- the scale would be applied in the wrong number of "
                         "places" % n_sub)
    # the gate value must not be carried twice: COMMON overrides it, and a
    # duplicate key is resolved by ParmParse's last-wins, which is implicit.
    # Drop the carried copy so what ran is unambiguous in argv.sh.
    out = [k for k in out if not k.startswith("surface_patch.side_face_h_r_min=")]
    return out


def run(tag, stop_time, extra, nproc=4):
    d = os.path.join(OUT, tag)
    os.makedirs(d, exist_ok=True)
    done = os.path.join(d, ".done")
    if os.path.exists(done):
        os.remove(done)
    argv = (["caffeinate", "-dis", "mpirun", "--oversubscribe", "--bind-to",
             "none", "-np", str(nproc), BIN, DECK]
            + base_keys() + COMMON + extra
            + ["stop_time=%g" % stop_time,
               "amr.thermo.int=1", "amr.thermo.plot_int=6",
               "surface_patch.jet_profile_interval=2.0",
               "energy_ledger.enabled=1", "spall.removal_events_csv=1",
               "plot_file=%s" % os.path.join(d, "run")])
    with open(os.path.join(d, "argv.sh"), "w") as f:
        f.write(" ".join("'%s'" % a if " " in a else a for a in argv) + "\n")
    env = os.environ.copy()
    env.setdefault("OMPI_MCA_btl", "self,vader,tcp")
    t0 = time.time()
    with open(os.path.join(d, "run.log"), "w") as lg:
        rc = subprocess.call(argv, cwd=ROOT, env=env, stdout=lg,
                             stderr=subprocess.STDOUT)
    wall = time.time() - t0
    with open(done, "w") as f:
        f.write("rc=%d wall=%.1f\n" % (rc, wall))
    return rc, wall, d


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "short"
    want = sys.argv[2:] or list(LEGS)
    print("P2-3 -- %s  legs: %s" % (mode, ", ".join(want)))
    print("  binary %s" % subprocess.run(["shasum", "-a", "256", BIN],
                                         capture_output=True,
                                         text=True).stdout[:24])
    print("  base argv %s (%d keys), floor peak -> %s"
          % (os.path.relpath(ARGV_FILE, ROOT), len(base_keys()), NEW_PEAK))
    print("  common: %s" % " ".join(COMMON))
    if mode == "reduce":
        import reduce_p2_3
        return reduce_p2_3.main(want)
    T = {"short": 60.0, "full": 500.0}.get(mode)
    if T is None:
        print("  usage: run_p2_3.py short | full | reduce [legs...]")
        return 2
    rcs = {}
    for tag in want:
        if tag not in LEGS:
            print("  unknown leg %s" % tag)
            return 2
        rc, wall, d = run(tag + ("" if mode == "full" else "_s"), T, LEGS[tag])
        print("  %-6s rc=%d  wall %.1f s (%.2f wall-s per simulated s)"
              % (tag, rc, wall, wall / T))
        rcs[tag] = rc
    bad = [t for t, r in rcs.items() if r != 0]
    if bad:
        print("  ABORTED: %s -- not a physical result" % ", ".join(bad))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
