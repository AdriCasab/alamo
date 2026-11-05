#!/usr/bin/env python3
"""P2-5 -- the Hoeser lab leg: the floor criterion at a MEASURED flux.

Builds every leg's argv the P2-3 way: BASE = P2-3's W1 argv (the Grimsel/spall
configuration, last-wins) with the feet / map / march / body blocks DROPPED,
then the lab COMMON block, then the leg's own keys. The drop list and the final
argv are written to each leg directory, because a leg's argv is part of its
provenance (standing rule 8) and `metadata` is a RESOLVED dump, i.e.
reconstruction and not replay.

REV 2. `jet_core_length` and `nozzle_z0` are LEG keys (the bracket K30/K20 and
the L 13 standoff leg move them); `h_expr` is read per LEG from MAP.txt, since
the lid constant and the standoff law differ by mesh and by L. The flat-h
closure probe lives in `probe/SMOKE0`, NOT in `legs/`: it is not a result.

LEGS -- the table lives in p2_5_legs.py, shared with derive_lab_map.py.
  probe/SMOKE0  2 s, 2 mm, a PROVISIONAL flat map -- D_e only.
          ⚠ Deliverables 2 and 3 are circular as written: the map is derived
          from "a 2 s SMOKE run's thermo" while the smoke of (3.) already runs
          "shape B", which is the map's output. They are separable because the
          closure's T_stag does NOT depend on h at all -- phi = min(1,
          jet_core_length*D_dec/s_c) and T_stag = T_mix + (T_nozzle-T_mix)*phi
          read only the geometry (MMWSpalling.H:3449-3456). SMOKE0 therefore
          carries a flat h and is used for NOTHING but reading jet_T_stag,
          jet_s_c and jet_D_e. Declared, not hidden.
  SMOKE   2 s, 2 mm, shape B -- deliverable 3 / criterion a2.
  LA      shape A, 30 s, 1 mm -- the pre-registered controlled failure (P1).
  LB      shape B, 30 s, 1 mm -- the headline (P2).
  LB2     shape B, 30 s, 2 mm -- the mesh twin (P4). weibull.V0 x8.
  LB120   shape B, 120 s, 1 mm -- the transient (P5).
  LC      shape C, 30 s, 1 mm -- ONLY on a c2 overshoot.

USAGE
  run_p2_5.py SMOKE0 SMOKE          # ord 1
  run_p2_5.py LA LB LB2 LB120       # ord 3-5 (serial)
  run_p2_5.py --probe LB            # 1-step parse probe: rc + unused keys
  run_p2_5.py --argv LB             # print the argv and exit (no run)
"""
import argparse
import os
import shlex
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from p2_5_legs import LEGS as LEGTAB, ORDER
STUD = os.path.dirname(HERE)
ROOT = os.path.abspath(os.path.join(STUD, "..", "..", ".."))
OUT = os.path.join(HERE, "legs")
PROBE = os.path.join(HERE, "probe")
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
DECK = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "hoeser", "input_lab")
W1 = os.path.join(STUD, "p2_3_wall", "legs", "W1", "argv.sh")
MAPTXT = os.path.join(HERE, "MAP.txt")

# ---------------------------------------------------------------------------
# FORCED DROPS. Three groups, each with its reason. The feet group is forced by
# the code (`nozzle_descent = prescribed` refuses it: MMWSpalling.H:10735 aborts
# on foot_body_clearance, :10733 on jet_free_surface); the rest is the packet's
# scope (no body in the hole, no march, no above-plane wall term, and the Meier
# h map replaced by the lab map).
DROP_FEET = ["surface_patch.nozzle_descent",
             "surface_patch.foot_r_inner", "surface_patch.foot_r_outer",
             "surface_patch.foot_standoff", "surface_patch.foot_rule",
             "surface_patch.foot_npads", "surface_patch.foot_pad_quantile",
             "surface_patch.foot_body_clearance",
             "surface_patch.jet_free_surface", "surface_patch.jet_fs_aspect",
             "surface_patch.jet_fs_exponent"]
DROP_SCOPE = ["surface_patch.h_expr",               # replaced by the lab map
              "surface_patch.jet_wall_march",       # no annulus to march
              "surface_patch.side_face_h_scale",    # a march key (P2-3)
              "surface_patch.side_face_factor_above",
              "surface_patch.jet_body", "surface_patch.jet_body_T",
              "surface_patch.jet_body_h_scale", "surface_patch.jet_body_mdot",
              "surface_patch.step_diag_radius"]     # Meier's pad ring
DROP_OVERRIDDEN = ["amr.blocking_factor", "timestep", "stop_time", "plot_file", "amr.plot_int",
                   "amr.thermo.plot_int", "amr.thermo.int", "amr.n_cell",
                   "amr.max_grid_size", "geometry.prob_lo", "geometry.prob_hi",
                   "weibull.V0", "energy_ledger.enabled",
                   "spall.removal_events_csv",
                   "surface_patch.radius", "surface_patch.x0",
                   "surface_patch.y0", "surface_patch.face",
                   "surface_patch.jet_T_ent_mode", "surface_patch.jet_T_ent",
                   "surface_patch.jet_mdot", "surface_patch.jet_cp",
                   "surface_patch.jet_D", "surface_patch.jet_T_nozzle",
                   "surface_patch.jet_core_length", "surface_patch.jet_De_ref",
                   "surface_patch.nozzle_z0",
                   "surface_patch.jet_profile_interval",
                   "surface_patch.nozzle_collision_radius"]
DROP = DROP_FEET + DROP_SCOPE + DROP_OVERRIDDEN

# ---------------------------------------------------------------------------
# The lab block, literal from ACTIVE_STEP.md goal 1.
#   jet_mdot  = the quarter's share of 13.6 g/s: CH4 0.684 g/s at 55.5 MJ/kg
#               (= 38.0 kW), air x 17.2 x 1.1  =>  13.63 g/s, /4 = 3.4e-3.
#   jet_De_ref= 0.0079 x the Meier deck's De_ref/D ratio 0.4740795578655565.
#   nozzle_z0 = 0.1632 = the 0.10 m rock top + her 63.2 mm (L = 8) standoff,
#               with nozzle_feed = 0 so the standoff only ever OPENS.
COMMON = [
    "surface_patch.enabled=1",
    "surface_patch.face=zhi",
    "surface_patch.x0=0.0", "surface_patch.y0=0.0",
    # ⚠ beyond the domain corner (0.0707 m) so the HARD circular mask at
    # MMWSpalling.H:1398/:1714/:2108/:2579 never bites inside the opening;
    # the flange is the lid in h_expr instead.
    "surface_patch.radius=0.08",
    "surface_patch.nozzle_descent=prescribed",
    # ⚠ REV 2: `nozzle_z0` and `jet_core_length` are LEG keys now, not COMMON
    # -- the K30/K20 bracket and the L 13 leg move them.
    "surface_patch.nozzle_feed=0.0",
    # standing rule 6: defaults to `radius` under prescribed = no guard.
    "surface_patch.nozzle_collision_radius=0.004",
    "surface_patch.jet_closure=enthalpy",
    "surface_patch.jet_stagnation=decay",
    "surface_patch.jet_entrained_mass=1",
    "surface_patch.jet_T_ent_mode=fixed",
    "surface_patch.jet_T_ent=293.15",
    "surface_patch.jet_mdot=3.4e-3",
    "surface_patch.jet_cp=1250.0",
    "surface_patch.jet_D=0.0079",
    "surface_patch.jet_T_nozzle=1866.0",
    "surface_patch.jet_decay_diameter=momentum",
    "surface_patch.jet_De_ref=0.003745",
    "surface_patch.jet_bin_update=exponential",
    "surface_patch.jet_negative_flux=count",
    "surface_patch.jet_profile_interval=1.0",
    "confining.p=0.0",
    "amr.thermo.plot_int=6", "amr.thermo.int=1",
    "energy_ledger.enabled=1", "spall.removal_events_csv=1",
]

MESH1 = ["amr.n_cell=50 50 100", "amr.max_grid_size=25",
         "geometry.prob_lo=0.0 0.0 0.0", "geometry.prob_hi=0.05 0.05 0.10",
         "timestep=0.008", "weibull.V0=1.0e-9", "amr.plot_int=625"]
# ⚠ PACKET DEFECT, fixed here and declared: the 2 mm twin's `amr.n_cell=25 25
# 50` is NOT divisible by the deck's `amr.blocking_factor = 2` (inherited from
# the Meier deck), and AMReX aborts at init -- "domain size in direction 0 is
# 25, blocking_factor is 2". Measured, rc 6 on the first SMOKE0 attempt. The
# twin must keep the SAME 50 x 50 x 100 mm box as the 1 mm leg or the mesh row
# compares two geometries, so the blocking factor moves instead: at
# amr.max_level = 0 with no regrid it only constrains the box decomposition,
# and decomposition-independence of this closure is already a measured property
# (P2-4 (a1) row (h): 0 of 12,499 printed tokens differ between 4 and 5 ranks).
MESH2 = ["amr.n_cell=25 25 50", "amr.max_grid_size=25", "amr.blocking_factor=1",
         "geometry.prob_lo=0.0 0.0 0.0", "geometry.prob_hi=0.05 0.05 0.10",
         "timestep=0.016", "weibull.V0=8.0e-9", "amr.plot_int=313"]

# A flat, deliberately meaningless coefficient for SMOKE0. Any positive value
# gives the same jet_T_stag; 1000 keeps the 2 s run numerically ordinary.
# A flat, deliberately meaningless coefficient for the probe. Any positive value
# gives the same jet_D_e / jet_T_stag; 1000 keeps the 2 s run numerically
# ordinary. ⚠ Quoted for D_e only, never as a result (packet guardrail).
FLAT = "surface_patch.h_expr=1000.0"

LEGS = {"SMOKE0": MESH2 + ["stop_time=2.0", "surface_patch.jet_core_length=30.0",
                           "surface_patch.nozzle_z0=0.1632", FLAT]}
for _n in ORDER:
    _dz, _st, _K, _z0, _sh = LEGTAB[_n]
    LEGS[_n] = ((MESH1 if _dz == 0.001 else MESH2)
                + [f"stop_time={_st}",
                   f"surface_patch.jet_core_length={_K}",
                   f"surface_patch.nozzle_z0={_z0}",
                   f"@MAP:{_n}"])
# The probe is not a leg: it writes to probe/, every real leg to legs/.
PROBE_LEGS = {"SMOKE0"}


def w1_base():
    """W1's argv with the drop list applied, token-aware.

    ⚠ TOKEN-AWARE: ParmParse multi-value keys arrive as `k=v0 v1 v2`, one `=`
    token then bare ones. Dropping only the `=` token leaves `60 300` dangling
    and ParmParse re-binds them to whatever key came before (d2ac_feed/run.py).
    ⚠ LAST-WINS: W1's argv has DUPLICATE keys (side_face_factor_above appears
    as 0.0 then 1.0) and ParmParse takes the FIRST. Rebuilding here with
    last-wins semantics is what P2-4 had to do after reading the wall OFF.
    """
    toks = shlex.split(open(W1).read())
    i = toks.index(DECK.replace("/tests/MMWSpalling/validation/hoeser/input_lab",
                                "")) if False else 0
    # keep everything after the input deck path
    for n, t in enumerate(toks):
        if t.endswith("input_drilling"):
            i = n + 1
            break
    else:
        sys.exit("W1 argv: could not find the input deck path")
    keys, j = [], i
    while j < len(toks):
        k = toks[j]
        if "=" in k and k.split("=", 1)[0] in DROP:
            j += 1
            while j < len(toks) and "=" not in toks[j]:
                j += 1
            continue
        keys.append(k)
        j += 1
    # collapse `k=v0 v1 v2` runs into single entries, then de-duplicate
    # LAST-WINS so a value the packet overrides later actually takes effect.
    merged, j = [], 0
    while j < len(keys):
        e = keys[j]
        j += 1
        while j < len(keys) and "=" not in keys[j]:
            e += " " + keys[j]
            j += 1
        merged.append(e)
    seen = {}
    for e in merged:
        seen[e.split("=", 1)[0]] = e
    left = [e for e in seen.values() if e.split("=", 1)[0] in set(DROP)]
    assert not left, f"dropped keys survived: {left}"
    return list(seen.values())


def read_map(leg):
    """The h_expr string for one LEG, from MAP.txt (derive_lab_map.py).

    ⚠ Per LEG and not per shape: the lid constant depends on the mesh and on
    nozzle_z0, and the standoff exponent on K. Two legs of the same shape
    (G05 / G05M2 / G05K20 / G05L13) carry four different strings."""
    if not os.path.exists(MAPTXT):
        sys.exit(f"{MAPTXT} missing -- run derive_lab_map.py first (ord 1)")
    for ln in open(MAPTXT):
        if ln.startswith(f"H_EXPR_{leg}="):
            return "surface_patch.h_expr=" + ln.split("=", 1)[1].strip()
    sys.exit(f"MAP.txt has no H_EXPR_{leg}= line")


def argv_of(name):
    base_dir = PROBE if name in PROBE_LEGS else OUT
    leg = os.path.join(base_dir, name)
    extra = []
    for e in LEGS[name]:
        extra.append(read_map(e.split(":", 1)[1]) if e.startswith("@MAP:") else e)
    # The leg block wins over the base: drop any base key it names.
    over = {e.split("=", 1)[0] for e in extra + COMMON if "=" in e}
    base = [e for e in w1_base() if e.split("=", 1)[0] not in over]
    keys = base + COMMON + extra + ["plot_file=" + os.path.join(leg, "run")]
    # ⚠ Rank count from the BOX COUNT, not a literal. 1 mm: 50/25 x 50/25 x
    # 100/25 = 16 boxes -> 4 ranks, 4 each. 2 mm: 25 x 25 x 50 at mgs 25 = 2
    # boxes -> 2 ranks (4 would idle two). The packet's "-np 5 if the grid is 5
    # boxes" is a P2-4 fact about the 60^3 Meier deck and does not apply here.
    nc = [e for e in extra if e.startswith("amr.n_cell=")][0]
    n = [int(v) for v in nc.split("=", 1)[1].split()]
    nbox = 1
    for v in n:
        nbox *= max(1, -(-v // 25))
    nb = min(4, nbox)
    cmd = (["caffeinate", "-dis", "mpirun", "--oversubscribe", "--bind-to",
            "none", "-np", str(nb), BIN, DECK] + keys)
    return leg, cmd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("legs", nargs="+")
    ap.add_argument("--probe", action="store_true",
                    help="1-step parse probe: rc + ParmParse unused keys")
    ap.add_argument("--argv", action="store_true", help="print argv, do not run")
    a = ap.parse_args()
    for name in a.legs:
        if name not in LEGS:
            sys.exit(f"unknown leg {name}; have {sorted(LEGS)}")
    for name in a.legs:
        leg, cmd = argv_of(name)
        if a.probe:
            # ⚠ The probe must redirect BOTH its log and its plot_file. The
            # first rev moved only the log directory, so a one-step probe wrote
            # its `run/` into the REAL leg's directory -- stale single-step
            # output sitting where a result belongs, which is exactly the
            # provenance confusion standing rule 8 is about.
            leg = leg + "_probe"
            cmd = [c for c in cmd
                   if not c.startswith(("stop_time=", "plot_file="))]
            cmd += ["stop_time=1e-9", "plot_file=" + os.path.join(leg, "run")]
        os.makedirs(leg, exist_ok=True)
        if a.argv:
            print(" \\\n  ".join(shlex.quote(c) for c in cmd))
            continue
        with open(os.path.join(leg, "argv.sh"), "w") as f:
            f.write(" ".join(shlex.quote(c) for c in cmd) + "\n")
        with open(os.path.join(leg, "DROPPED.txt"), "w") as f:
            f.write("P2-5 forced key drops from P2-3 W1's argv\n")
            for tag, lst in (("FEET (forced by the code)", DROP_FEET),
                             ("SCOPE (the packet)", DROP_SCOPE),
                             ("OVERRIDDEN by the lab block", DROP_OVERRIDDEN)):
                f.write(f"\n-- {tag}\n")
                for k in lst:
                    f.write(f"   {k}\n")
        t0 = time.time()
        log = os.path.join(leg, "run.log")
        with open(log, "w") as f:
            rc = subprocess.call(cmd, stdout=f, stderr=subprocess.STDOUT,
                                 cwd=ROOT)
        wall = time.time() - t0
        print(f"  {name:7s} rc={rc} np={cmd[7]} wall {wall:.1f} s  -> {log}")
        if a.probe:
            txt = open(log).read()
            bad = [l for l in txt.splitlines()
                   if "nused" in l or "Abort" in l or "ERROR" in l]
            print("    " + ("\n    ".join(bad[:25]) if bad
                            else "no unused-key or abort lines"))
        if rc == 0:
            open(os.path.join(leg, ".done"), "w").write(f"{{'rc': 0}}\n")


if __name__ == "__main__":
    main()
