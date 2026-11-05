#!/usr/bin/env python3
"""P2-2e -- the constraint depth: every leg's argv, built from the decks and
harnesses the packet names, each key emitted ONCE (last value wins here, so a
duplicated key can never be resolved by ParmParse's own rule -- P2-4's W1
trap). Every leg writes argv.sh (standing rule 8), run.log and a .done sentinel.

GROUPS
  ord 2  S0s                 S0's argv, stop_time 60, key 0           (c6)
         A3                  decomposition legs at d = 8.5 mm         (a3)
  ord 3  K0_h500             Kant column, Grimsel props, h 500, key 0 (b1)
         LADDER              C<d>_<SOD>  Gotthard thermal             (b2)
         GLADDER             G<d>_<SOD>  Grimsel thermal   (reported, never selected)
         SLADDER             S<d>_<SOD>  Gotthard + robin_form=surface (reported)
  ord 4  CHECKS              V50 V80 V110 V_h500 V_h150 M80_2 M80_1 M80_05
                             GB B2000_0 B2000_d B100_0 B100_d       (b3, c1-c4)
  ord 5  H0 HD               Hoeser G05 on a 200 mm block, key 0 / d_ref (c5)
  ord 6  SD                  S0's argv + d_ref, 500 s, -np 5          (c7-c11)

d_ref for ord 4-6 is READ from D_REF.txt (written by reduce_cref.py CAL); it is
chosen once by the b2 rule and never typed here.

Usage: run_cref.py <group|leg> [...] [--jobs N] [--probe] [--argv]
"""
import argparse
import concurrent.futures as cf
import os
import re
import shlex
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)
ROOT = os.path.abspath(os.path.join(STUD, "..", "..", ".."))
OUT = os.path.join(HERE, "legs")
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
KEY = "spallation.sp.constraint_depth"

sys.path.insert(0, os.path.join(STUD, "p2_2d_kic"))
import run_prediction as RP                                   # noqa: E402
sys.path.insert(0, os.path.join(STUD, "s1_surface_resolution"))
import run_sweep as S1                                        # noqa: E402

KANT_DECK = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "kant",
                         "sp_kant_onset_meiercfg", "input_model")
GB_DECK = os.path.join(ROOT, "tests", "MMWSpalling", "unit", "gradient_blind",
                       "input")
S0_ARGV = os.path.join(STUD, "p2_4_body", "legs", "S0", "argv.sh")
G05_ARGV = os.path.join(STUD, "p2_5_hoeser", "legs", "G05", "argv.sh")

# Kant 2016 Table 4: SOD [mm] -> (h_fl [W/m2K], T_f [K]); tabulated only.
SOD = {50: (766, 1723.15), 80: (1009, 1683.15), 110: (873, 1633.15),
       130: (747, 1523.15), 150: (633, 1373.15)}
ANCH = {"h500": (500, 1673.15), "h150": (150, 1673.15)}
# Calibration rock = Kant's: Gotthard THERMAL properties on the Grimsel law
# shape (lambda(25 C) = 1923.5/(25 + 563) = 3.27 W/mK; same theta0). The
# STRENGTH side stays Grimsel's, so the level is the banked one.
GOTTHARD = ["microstructure.phase0.grimsel.lambda.a=1923.5",
            "microstructure.phase0.grimsel.rho=2594"]
LADDER_D = [round(3.0 + 0.5 * i, 1) for i in range(19)]      # 3.0 .. 12.0 mm
# ⚠ THE PRE-REGISTERED EXTENSION (b2 rule): SOD 130 FIRED at the bottom rung
# (3.0 mm), so the ladder extends DOWN by 2 mm, repeated to 1 mm. Measured, not
# chosen: every rung 3.0-12.0 fired at both standoffs on all three ladders.
LADDER_DN = [1.0, 1.5, 2.0, 2.5]
A3_D = 8.5e-3                                                 # decomposition legs only


# THE CONSTANT-PROPERTY DIAGNOSTIC (not a packet leg; it decides whether the
# calibrated d_ref is physics or a read error). The OFFLINE model's own
# assumptions in the binary: constant conductivity and diffusivity (Gotthard
# 3.27 W/mK, 1.6e-6 m2/s), constant E/beta (law and stress_law = constant) and
# constant toughness (linear law with T_melt -> infinity; T_melt must be > 0,
# so 1e12 K makes K_Ic constant to 1e-9). If the reference is read at the right
# depth, this ladder must reproduce the offline interval for the SAME threshold
# the binary itself reports at key 0.
CONST = ["microstructure.phase0.law=constant", "spall.stress_law=scalar",
         "microstructure.phase0.kappa=3.27", "microstructure.phase0.rho=2587.2",
         "microstructure.phase0.Cp=790.0", "hc.heat.alpha=1.6e-6",
         "spallation.sp.kic.law=linear", "spallation.sp.kic.K_Ic0=0.563",
         "spallation.sp.kic.T_ref=293.15", "spallation.sp.kic.T_melt=1.0e12"]
CONST_D = [round(1.0 + 0.5 * i, 1) for i in range(23)]       # 1.0 .. 12.0 mm


# ------------------------------------------------------------------ helpers
def kv_list(toks):
    """argv tokens after the deck -> [(key, value)] with multi-token values
    re-joined (`amr.n_cell=50 50 100` arrives as three tokens)."""
    out, i = [], 0
    while i < len(toks):
        t = toks[i]
        i += 1
        if "=" not in t:
            continue
        k, v = t.split("=", 1)
        while i < len(toks) and "=" not in toks[i]:
            v += " " + toks[i]
            i += 1
        out.append((k, v))
    return out


def merge(*groups):
    """Ordered, de-duplicated, LAST value wins -- emitted once each."""
    d = {}
    for g in groups:
        for e in g:
            if isinstance(e, tuple):
                k, v = e
            else:
                k, v = e.split("=", 1)
            d.pop(k, None)
            d[k] = v
    return [f"{k}={v}" for k, v in d.items()]


def from_argv(path, deck_suffix):
    """(deck, [(k, v)]) from a recorded argv.sh."""
    toks = shlex.split(open(path).read())
    for n, t in enumerate(toks):
        if t.endswith(deck_suffix):
            return t, kv_list(toks[n + 1:])
    sys.exit(f"{path}: no deck ending in {deck_suffix}")


def d_ref():
    p = os.path.join(HERE, "D_REF.txt")
    if not os.path.exists(p):
        sys.exit("D_REF.txt missing -- run the calibration (ord 3) and "
                 "`reduce_cref.py CAL > D_REF.txt` first")
    m = re.search(r"^D_REF_M=([0-9.eE+-]+)\s*$", open(p).read(), re.M)
    if not m:
        sys.exit("D_REF.txt has no D_REF_M= line")
    return float(m.group(1))


# ------------------------------------------------------------- leg builders
def kant(h, Tf, d_m, thermal="gotthard", extra=(), geom=("4 4 150", "0.150"),
         dt=None, stop="200"):
    """The P2-2d Kant configuration (run_prediction COMMON + P_grimsel) on a
    deep column. ⚠ COMMON does NOT pass surface_patch.robin_form, so the column
    inherits the deck's (the BC the banked level was measured in) -- kept
    verbatim; SLADDER overrides it on purpose and is reported only."""
    p_grimsel = dict(RP.LEGS_DICT)["P_grimsel"] if hasattr(RP, "LEGS_DICT") \
        else [e for (t, x, _) in RP.LEGS if t == "P_grimsel" for e in x]
    nc, depth = geom
    keys = merge(RP.COMMON, p_grimsel,
                 [f"geometry.prob_hi=0.004 0.004 {depth}",
                  f"amr.n_cell={nc}", "amr.max_grid_size=25",
                  f"stop_time={stop}",
                  f"surface_patch.h_conv={h}", f"surface_patch.T_flame={Tf}",
                  "amr.plot_int=4000"],
                 GOTTHARD if thermal == "gotthard" else [],
                 [f"timestep={dt}"] if dt else [],
                 list(extra),
                 [f"{KEY}={d_m!r}"] if d_m > 0 else [])
    return KANT_DECK, keys


def beam(alpha, d_m):
    """s1 column at 1 mm, 5 MW/m2, via run_sweep.overrides() with a synthetic
    beam case, plus the P2-2d criterion keys."""
    case = f"P2E-{int(alpha)}"
    S1.CASES[case] = ({"kind": "beam", "q": 5.0e6, "alpha": float(alpha)},
                      None, None)
    ov = S1.overrides(case, 1.0, 30.0, "@PLOT@")
    ov = [e for e in ov if not e.startswith("plot_file=")]
    crit = ["microstructure.phase0.law=grimsel", "spall.stress_law=grimsel",
            "spallation.sp.kic.law=grimsel_table", "spallation.sp.kic.scale=1.0",
            "spallation.sp.stress_relief=%.6f" % RP.R_DEF,
            "spall.surface_T=reconstructed", "spallation.sp.a0=8.0e-6",
            "weibull.enabled=0", "spall.removal_events_csv=1",
            "energy_ledger.enabled=1", "amr.thermo.int=1",
            "amr.thermo.plot_int=10"]
    return S1.INPUT, merge(ov, crit, [f"{KEY}={d_m!r}"] if d_m > 0 else [])


def hoeser(d_m, block200=True, stop="30.0", extra=()):
    deck, kv = from_argv(G05_ARGV, "input_lab")
    kv = [(k, v) for (k, v) in kv if k != "plot_file"]
    g = []
    if block200:
        # ⚠ the standoff and the lid constant are unchanged: the rock top moves
        # 0.10 -> 0.20 and the nozzle 0.1632 -> 0.2632, s0 = 63.2 mm either way,
        # and h_expr reads s (the standoff), never z.
        g = ["geometry.prob_hi=0.05 0.05 0.20", "amr.n_cell=50 50 200",
             "surface_patch.nozzle_z0=0.2632"]
    return deck, merge(kv, g, [f"stop_time={stop}"], list(extra),
                       [f"{KEY}={d_m!r}"] if d_m > 0 else [])


def meier(d_m, stop="500"):
    deck, kv = from_argv(S0_ARGV, "input_drilling")
    kv = [(k, v) for (k, v) in kv if k != "plot_file"]
    return deck, merge(kv, [f"stop_time={stop}"],
                       [f"{KEY}={d_m!r}"] if d_m > 0 else [])


# ------------------------------------------------------------------ catalogue
def catalogue(names):
    """{leg: (deck, keys, np)} for the requested groups/legs."""
    legs = {}
    for nm in names:
        if nm == "S0s":
            dk, ks = meier(0.0, stop="60")
            # ⚠ -np 4, NOT 5: S0 itself ran on 4 ranks (its argv.sh), and c6
            # is a BYTE comparison -- a rank-count change must not stand in
            # for a binary change.
            legs["S0s"] = (dk, ks, 4)
        elif nm == "A3":
            h, Tf = SOD[130]
            for np_, mgs in ((1, 25), (2, 25), (4, 25), (4, 50)):
                dk, ks = kant(h, Tf, A3_D, extra=[f"amr.max_grid_size={mgs}"],
                              stop="200")
                legs[f"A3C_np{np_}_m{mgs}"] = (dk, ks, np_)
            for np_, mgs in ((1, 25), (2, 25), (4, 25), (4, 50)):
                dk, ks = hoeser(A3_D, block200=False, stop="5.0",
                                extra=[f"amr.max_grid_size={mgs}",
                                       "amr.plot_int=625"])
                legs[f"A3H_np{np_}_m{mgs}"] = (dk, ks, np_)
            # (iii) per-cell removal ON: the Sp_field site's verdicts reach the
            # event log only through PerCellRemoval.
            # ⚠ mgs 8 puts x/y seams at i, j = 8, 16 -- THROUGH the pit walls
            # (the jet axis is the domain corner, the pit ~19 mm wide). At mgs
            # 25 the only x/y seam is at 25, which the walls never reach, so the
            # lateral branch's one-cell horizontal read never crossed a seam.
            for np_, mgs in ((2, 25), (4, 50), (4, 8)):
                dk, ks = hoeser(A3_D, block200=False, stop="5.0",
                                extra=[f"amr.max_grid_size={mgs}",
                                       "amr.plot_int=625",
                                       "spall.per_cell_removal=1"])
                legs[f"A3P_np{np_}_m{mgs}"] = (dk, ks, np_)
            dk, ks = hoeser(A3_D, block200=False, stop="5.0",
                            extra=["amr.max_grid_size=8", "amr.plot_int=625"])
            legs["A3H_np4_m8"] = (dk, ks, 4)
        elif nm == "K0_h500":
            dk, ks = kant(*ANCH["h500"], 0.0, thermal="grimsel")
            legs["K0_h500"] = (dk, ks, 2)
        elif nm.rstrip("_DN") in ("LADDER", "GLADDER", "SLADDER"):
            base = nm[:-3] if nm.endswith("_DN") else nm
            tag = {"LADDER": "C", "GLADDER": "G", "SLADDER": "S"}[base]
            th = "grimsel" if base == "GLADDER" else "gotthard"
            ex = ["surface_patch.robin_form=surface"] if base == "SLADDER" else []
            for d in (LADDER_DN if nm.endswith("_DN") else LADDER_D):
                for sod in (130, 150):
                    dk, ks = kant(*SOD[sod], d * 1e-3, thermal=th, extra=ex)
                    legs[f"{tag}{d:.1f}_{sod}"] = (dk, ks, 1)
        elif re.fullmatch(r"[CGS]\d+(\.\d)?_(130|150)", nm):
            tag, rest = nm[0], nm[1:]
            d, sod = rest.split("_")
            th = "grimsel" if tag == "G" else "gotthard"
            ex = ["surface_patch.robin_form=surface"] if tag == "S" else []
            dk, ks = kant(*SOD[int(sod)], float(d) * 1e-3, thermal=th, extra=ex)
            legs[nm] = (dk, ks, 1)
        elif nm == "CONST":
            for nmk, (h, Tf) in (("K0C_130", SOD[130]), ("K0C_150", SOD[150]),
                                  ("K0C_h500", ANCH["h500"])):
                legs[nmk] = (*kant(h, Tf, 0.0, thermal="grimsel", extra=CONST), 1)
            for d in CONST_D:
                for sod in (130, 150):
                    legs[f"X{d:.1f}_{sod}"] = (*kant(*SOD[sod], d * 1e-3,
                                                      thermal="grimsel",
                                                      extra=CONST), 1)
        elif nm == "CHECKS":
            dr = d_ref()
            for sod in (50, 80, 110):
                legs[f"V{sod}"] = (*kant(*SOD[sod], dr), 1)
            for an in ("h500", "h150"):
                legs[f"V_{an}"] = (*kant(*ANCH[an], dr, thermal="grimsel"), 1)
            for tag, nc in (("2", "2 2 50"), ("1", "4 4 100"), ("05", "8 8 200")):
                legs[f"M80_{tag}"] = (*kant(*SOD[80], dr,
                                             geom=(nc, "0.100"), dt="0.0125"), 1)
            # ⚠ amr.thermo.* set EXPLICITLY: the gradient_blind deck never sets
            # them and amr.thermo.plot_int defaults to -1 = NEVER, so the first
            # run exited rc 0 with no thermo.dat and no cref_n_vert witness.
            legs["GB"] = (GB_DECK, merge([f"{KEY}={dr!r}",
                                           "spall.removal_events_csv=1",
                                           "amr.thermo.int=1",
                                           "amr.thermo.plot_int=1"]), 4)
            for a in (2000, 100):
                legs[f"B{a}_0"] = (*beam(a, 0.0), 1)
                legs[f"B{a}_d"] = (*beam(a, dr), 1)
        elif nm in ("H0", "HD"):
            legs[nm] = (*hoeser(0.0 if nm == "H0" else d_ref()), 4)
        elif nm == "SD":
            legs["SD"] = (*meier(d_ref()), 5)
        else:
            sys.exit(f"unknown group/leg {nm}")
    return legs


def run_leg(name, deck, keys, np_, probe=False, stop=None):
    # ⚠ a probe or a smoke writes to its OWN directory, plot_file included,
    # never into the real leg's (the P2-5 --probe wart).
    sfx = "_probe" if probe else (f"_smoke{stop}" if stop else "")
    d = os.path.join(OUT, name + sfx)
    os.makedirs(d, exist_ok=True)
    ks = list(keys)
    if probe:
        ks = merge(ks, ["stop_time=1e-9"])
    elif stop:
        ks = merge(ks, [f"stop_time={stop}"])
    cmd = (["mpirun", "--oversubscribe", "--bind-to", "none", "-np", str(np_),
            BIN, deck] + ks + ["plot_file=" + os.path.join(d, "run")])
    with open(os.path.join(d, "argv.sh"), "w") as f:
        f.write(" ".join(shlex.quote(c) for c in cmd) + "\n")
    env = os.environ.copy()
    env.setdefault("OMPI_MCA_btl", "self,vader,tcp")
    t0 = time.time()
    with open(os.path.join(d, "run.log"), "w") as lg:
        rc = subprocess.call(cmd, cwd=ROOT, env=env, stdout=lg,
                             stderr=subprocess.STDOUT)
    wall = time.time() - t0
    if rc == 0:
        with open(os.path.join(d, ".done"), "w") as f:
            f.write("{'rc': 0, 'wall_s': %.1f}\n" % wall)
    unused = ""
    if probe:
        txt = open(os.path.join(d, "run.log")).read()
        m = re.findall(r"\[TOP\]::(\S+)\(nvals", txt)
        unused = " unused: " + (",".join(m) if m else "-")
    return name, rc, np_, wall, unused


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("names", nargs="+")
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--probe", action="store_true")
    ap.add_argument("--argv", action="store_true")
    ap.add_argument("--stop", default=None,
                    help="SMOKE: override stop_time, write to <leg>_smoke<stop>")
    ap.add_argument("--force", action="store_true",
                    help="re-run legs that already carry a .done")
    a = ap.parse_args()
    # ⚠ the harness holds the host awake itself (standing rule 11)
    subprocess.Popen(["caffeinate", "-dis", "-w", str(os.getpid())])
    legs = catalogue(a.names)
    if a.argv:
        for n, (dk, ks, np_) in legs.items():
            print(f"--- {n} (-np {np_})\n  {dk}\n  " + "\n  ".join(ks))
        return
    todo = [(n, v) for n, v in legs.items()
            if a.probe or a.force or a.stop
            or not os.path.exists(os.path.join(OUT, n, ".done"))]
    print(f"{len(todo)} of {len(legs)} legs to run, {a.jobs} at a time")
    with cf.ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = [ex.submit(run_leg, n, dk, ks, np_, a.probe, a.stop)
                for n, (dk, ks, np_) in todo]
        for f in cf.as_completed(futs):
            n, rc, np_, wall, un = f.result()
            print(f"  {n:16s} rc={rc} np={np_} wall {wall:7.1f} s{un}",
                  flush=True)


if __name__ == "__main__":
    main()
