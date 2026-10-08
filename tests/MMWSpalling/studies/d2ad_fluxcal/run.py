#!/usr/bin/env python3
"""D2ad: calibrate the FLOOR FLUX on Fig. 8.8's standoff, and band it.

  run.py --cases SMOKE_K13 SMOKE_K16 SMOKE_K19 SMOKE_K22 \
                 SMOKE_Q10 SMOKE_Q5 SMOKE_Q35 SMOKE_G1600 --jobs 2   item 0
  run.py --cases C0 G1600 Q5 Q55T870 Q5T870 Q5T920 Q35 S0 --jobs 2    the ladder
  run.py --list
  run.py --kill NAME    binary path + plot_file match (never pkill -f)

**NO SOURCE CHANGE.** Inputs only; the binary must still hash to ea833f2f... at
the end.

Base is D2ac's `FH` argv, obtained by CALLING d2ac_feed/run.py's own
`frozen_argv('FH')` -- the packet's instruction is to re-use the lifting
machinery, not to retype 56 keys. That carries D2ac's forced fixes for free:
`nozzle_collision_radius=0.028`, the ten dropped feet keys, lz = 0.5 / n_cell z
= 250, the prescribed feed, and the D2ab wall pair at 169.3.

THE TWO LEVERS
  F_q   a uniform multiplicative prefix on surface_patch.h_expr. VERIFIED: the
        lifted expression is 728 chars whose bracket depth returns to 0 ONLY at
        the final character, so `F*(<expr>)` needs no re-bracketing. It is built
        by STRING-PREFIXING the lifted value, never retyped.
        ⚠ h_expr is the FLOOR h under jet_closure=enthalpy (:2059 guards on
        s > 0.0, floor faces have s > 0, assigned to flame_h at :2067). This is
        NOT the D2u inert-lever trap, which was about ABOVE-plane faces.
  T_s   spallation.sp.kic.scale (parse :8194, live only under kic.law =
        nasseri_table which is the default :8043; fires at Removal.H:759).
        Absent from the lifted argv -- verified.

⚠ SMOKE_K10 would be identical to SMOKE_Q5 (F_q = 0.5, scale = 1.0), so it is
NOT run twice; the item-0(b) table reads SMOKE_Q5 as its scale-1.0 row.

⚠ KIC870 / KIC920 are set by ITEM 0(b) and must be filled in before the ladder.
The packet forbids using the brief's hand mapping (+/-10 %).

⚠ macOS pgrep ignores -a and prints bare PIDs, so --kill reads command lines from
`ps -o command= -p`.
"""
import argparse
import concurrent.futures as cf
import importlib.util
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


def _load(name, path):
    sp = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


D2AC = _load("d2ac_run", os.path.join(STUD, "d2ac_feed", "run.py"))

# ---- ITEM 0(b) OUTPUT. Filled from the smokes, before CRITERION.md is hashed.
# MEASURED by item 0(b) on the 60 s smokes: T_s(scale) = 821.3 / 864.6 / 895.5 /
# 919.1 / 950.1 K at scale 1.0 / 1.3 / 1.6 / 1.9 / 2.2, firing window 10.5-18.8 K
# (bar 30 K, PASS). Interpolated, not taken from the brief's hand mapping.
KIC870 = "1.3523"      # kic.scale for T_s = 870 K
KIC920 = "1.9088"      # kic.scale for T_s = 920 K

SMOKE_STOP = "60.0"
LADDER_STOP = "800.0"
C0_STOP = "650.0"      # C3: C0's floor would reach the adiabatic base by 800 s


def base_keys():
    """FH's argv, verbatim, through D2ac's own lifting machinery."""
    inp, keys = D2AC.frozen_argv("FH")
    assert any(k.startswith("surface_patch.h_expr=") for k in keys), "no h_expr"
    assert not any("kic.scale" in k for k in keys), "base already sets kic.scale"
    assert not any("pad_shield_open" in k for k in keys), "base shields the pads"
    return inp, keys


def lifted_h_expr(keys):
    h = [k for k in keys if k.startswith("surface_patch.h_expr=")]
    assert len(h) == 1, h
    e = h[0].split("=", 1)[1]
    # Re-verify the bracket structure on every call rather than trusting the
    # note: if the base argv is ever re-lifted from a different leg, a silent
    # re-bracketing would scale only part of the map.
    d, zeros = 0, []
    for i, c in enumerate(e):
        if c == "(":
            d += 1
        elif c == ")":
            d -= 1
            if d == 0:
                zeros.append(i)
    assert e[0] == "(" and zeros == [len(e) - 1], \
        "h_expr is not a single top-level bracket pair; prefixing is unsafe"
    return e


def leg(fq=None, scale=None, nozzle=None, stop=LADDER_STOP):
    """Return the extra key lines for one leg. `fq` prefixes h_expr."""
    return dict(fq=fq, scale=scale, nozzle=nozzle, stop=stop)


CASES = {
    # ---- the ladder ------------------------------------------------------
    "C0":       leg(stop=C0_STOP),
    "G1600":    leg(nozzle="1600.0"),
    "Q5":       leg(fq="0.5"),
    # ⚠⚠ THE LADDER IS RECOMPUTED -- see PREFLIGHT.md (g). Item 0 measured
    # R90 = 4.237 m/h, not the 3.3 the packet called validated (+28.4 %, over
    # its own 10 % recompute trigger), and the C1 back-reaction was 53.5 K out
    # (bar 40 K). Both recomputes are MANDATED by the packet itself. The
    # corrected window is F in (0.354, 0.393], i.e. s* in (90, 100] mm, and
    # s* = 254.2*F not 198*F.
    #   Q5T870  the packet's PRE-REGISTERED PASS -- kept, to test P1 directly
    #   Q35T870 / Q30  the two legs that now reach the window, by DIFFERENT
    #                  lever pairs (criterion (i)'s degeneracy, measured)
    #   Q30T870 nearest bracket BELOW (contact)
    # DROPPED: Q55T870 (145 mm) and Q5T920 (117 mm) -- both now fail high.
    "Q5T870":   leg(fq="0.5",  scale="KIC870"),
    "Q35":      leg(fq="0.35"),
    "Q35T870":  leg(fq="0.35", scale="KIC870"),
    "Q30":      leg(fq="0.30"),
    "Q30T870":  leg(fq="0.30", scale="KIC870"),
    # ---- item 0 smokes ---------------------------------------------------
    "SMOKE_Q10":   leg(stop=SMOKE_STOP),                       # 0(c) F_q = 1.0
    "SMOKE_Q5":    leg(fq="0.5",  stop=SMOKE_STOP),            # 0(c) + 0(b) s=1.0
    "SMOKE_Q35":   leg(fq="0.35", stop=SMOKE_STOP),            # 0(c) F_q = 0.35
    "SMOKE_G1600": leg(nozzle="1600.0", stop=SMOKE_STOP),      # 0(d)
    "SMOKE_K13":   leg(fq="0.5", scale="1.3", stop=SMOKE_STOP),  # 0(b)
    "SMOKE_K16":   leg(fq="0.5", scale="1.6", stop=SMOKE_STOP),
    "SMOKE_K19":   leg(fq="0.5", scale="1.9", stop=SMOKE_STOP),
    "SMOKE_K22":   leg(fq="0.5", scale="2.2", stop=SMOKE_STOP),
}


def frozen_argv(name):
    if name == "S0":
        return D2AC.frozen_argv("S0")
    c = CASES[name]
    inp, keys = base_keys()
    extra = []
    if c["fq"] is not None:
        e = lifted_h_expr(keys)
        # ⚠ ParmParse is LAST-VALUE-WINS, so appending replaces the lifted
        # token. Built by prefixing the lifted string, never retyped.
        extra.append(f'surface_patch.h_expr={c["fq"]}*{e}')
    if c["scale"] is not None:
        v = c["scale"]
        if v in ("KIC870", "KIC920"):
            v = {"KIC870": KIC870, "KIC920": KIC920}[v]
            assert v is not None, (
                f"{name}: item 0(b) has not set {c['scale']} yet -- run the "
                "SMOKE_K* smokes and fill it in before the ladder")
        extra.append(f"spallation.sp.kic.scale={v}")
    if c["nozzle"] is not None:
        extra.append(f"surface_patch.jet_T_nozzle={c['nozzle']}")
    extra.append(f"stop_time={c['stop']}")
    return inp, keys + extra


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
    c = CASES.get(name, {})
    m = {"name": name, "rc": rc, "seconds": round(wall_s, 1),
         "status": "ok" if rc == 0 else "FAILED",
         "F_q": c.get("fq") or "1.0",
         "kic_scale": kv.get("spallation.sp.kic.scale", ["1.0"])[0],
         "jet_T_nozzle": kv.get("surface_patch.jet_T_nozzle", ["1900.0"])[0],
         "nozzle_feed": kv.get("surface_patch.nozzle_feed", ["0"])[0]}
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
        for n in list(CASES) + ["S0"]:
            try:
                m = meta_of_argv(n)
            except AssertionError as e:
                print(f"{n:13s} (blocked) {e}")
                continue
            print(f"{n:13s} F_q={m['F_q']:>5s} kic.scale={m['kic_scale']:>5s} "
                  f"T_nozzle={m['jet_T_nozzle']:>7s} stop={m.get('stop','')}")
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
        if n not in CASES and n != "S0":
            ap.error(f"unknown case {n}; have {list(CASES) + ['S0']}")

    # Runs are not checkpointed: a host sleep loses the leg outright.
    subprocess.Popen(["caffeinate", "-dis", "-w", str(os.getpid())])

    print(f"launching {a.cases} with {a.jobs} job(s) at {time.strftime('%T')}")
    with cf.ThreadPoolExecutor(max_workers=a.jobs) as ex:
        rcs = list(ex.map(run_one, a.cases))
    print(f"all done at {time.strftime('%T')}")
    return 0 if all(r == 0 for r in rcs) else 1


if __name__ == "__main__":
    sys.exit(main())
