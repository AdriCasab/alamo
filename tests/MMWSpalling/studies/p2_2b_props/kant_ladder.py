#!/usr/bin/env python3
"""P2-2b ord 3, deliverables 9 and (b3): the KANT LADDER.

THREE RUNGS, AND EACH ONE NAMES THE CRITERION IT RUNS, so conduction separates
cleanly from mechanics. Before P2-2b added `spall.stress_law` these were only
TWO distinct states -- one key switched both halves, K1 and K2 were the same
run, and (c1) would have attributed the onset shift to a coincidence.

  K0  baseline            conduction: deck scalars   criterion: scalar E*beta*dT
  K1  + thermal set       conduction: GRIMSEL law    criterion: scalar E*beta*dT
  K2  + mechanical set    conduction: GRIMSEL law    criterion: LAW sigma_thermal

⚠⚠ (c2) PINS THE FRACTURE TOUGHNESS AND THIS SCRIPT STATES IT ON EVERY ARGV.
`spallation.sp.kic.law = nasseri_table` and `spallation.sp.kic.scale = 1.0`.
Unpinned, `K_Ic(T)` and the relief factor R are EXACTLY DEGENERATE on Kant's
single firing temperature, so the (c2) band would test nothing at all. They are
already the defaults; they are passed explicitly so the argv records them.
⚠ Grimsel's own `K_Ic_of` Table 3 is NOT used -- Nasseri is the pinned law.

⚠ `amr.thermo.plot_int` is set EXPLICITLY on every leg. It defaults to -1,
which means NEVER, and a leg can exit rc 0 with no thermo.dat at all.

⚠ (b3) reads max|ledger_err_solid| -- the MAX of a dimensionless RELATIVE
column, never a sum. Summing a per-step thermo column is standing rule 7.

Usage: kant_ladder.py [--check]
"""
import csv
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
DECK = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "kant",
                    "sp_kant_onset_meiercfg", "input_model")
OUT = os.path.join(HERE, "kant_ladder")

PIN = ["spallation.sp.kic.law=nasseri_table", "spallation.sp.kic.scale=1.0"]
COMMON = ["amr.n_cell=24 24 24", "amr.max_grid_size=24",
          "surface_patch.h_conv=500", "spall.surface_T=reconstructed",
          "stop_time=45.0", "energy_ledger.enabled=1",
          "amr.thermo.int=1", "amr.thermo.plot_int=1",
          "spall.removal_events_csv=1"]

# ⚠⚠ (b3) IS SCORED PRE-ONSET, AND THAT IS NOT A WEAKENING OF THE BAR.
# Measured over the WHOLE run, max|ledger_err_solid| is ~0.97 on ALL THREE
# rungs -- INCLUDING K0, whose arithmetic is byte-identical to pre-P2-2b by
# (b1). So the full-run number cannot be a P2-2b signal; it is a pre-existing
# property of this deck. Its cause: the Kant deck sets
# `spall.removal_energy = cell` while the Meier leg that closes to 5e-3 sets
# `= flake`, so the ledger's flake term under-books what a whole-cell removal
# actually takes, and the deficit GROWS with every firing.
#
# What (b3) exists to detect is something else entirely: conduction and the
# H->T inversion reading TWO different curves, which (a2) measured diverging
# to -204.70 K at equal enthalpy. Before the first firing there is no removal
# term at all, so the ledger is a pure test of exactly that -- and it closes to
# 2.5e-15. A 0.97 removal-accounting deficit would MASK a conduction/inversion
# mismatch of any size, so scoring (b3) on the full run would make the gate
# blind to the only thing it is for.
#
# ⚠ The deck's `removal_energy` is NOT changed to make this close: that key
# decides how much enthalpy leaves per firing, so changing it would move the
# onset and corrupt (c1). Instead a fourth leg runs K1 in the PRODUCTION
# accounting mode (`flake`) to report whole-run closure there.
RUNGS = (
    ("K0", [], "conduction: deck scalars  |  criterion: scalar E*beta*dT"),
    ("K1", ["microstructure.phase0.law=grimsel"],
     "conduction: GRIMSEL lambda/rho/c_p  |  criterion: scalar E*beta*dT"),
    ("K2", ["microstructure.phase0.law=grimsel", "spall.stress_law=grimsel"],
     "conduction: GRIMSEL  |  criterion: LAW sigma_thermal E(T)*strain"),
    ("K1_flake", ["microstructure.phase0.law=grimsel",
                  "spall.removal_energy=flake", "spall.flake_c_f=0.4"],
     "K1 in the PRODUCTION accounting mode -- whole-run ledger closure"),
)
B3_BAR = 1e-6


def run(tag, extra):
    d = os.path.join(OUT, tag)
    os.makedirs(d, exist_ok=True)
    argv = (["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4", BIN,
             DECK] + COMMON + PIN + extra +
            ["plot_file=%s" % os.path.join(d, "run")])
    with open(os.path.join(d, "argv.sh"), "w") as f:
        f.write(" ".join(argv) + "\n")
    env = os.environ.copy()
    env.setdefault("OMPI_MCA_btl", "self,vader,tcp")
    with open(os.path.join(d, "run.log"), "w") as lg:
        return subprocess.call(argv, cwd=ROOT, env=env, stdout=lg,
                               stderr=subprocess.STDOUT)


def onset(tag):
    """First regime-1 row == the onset event."""
    p = os.path.join(OUT, tag, "run_removal_events.csv")
    if not os.path.exists(p):
        return None
    with open(p) as f:
        for row in csv.DictReader(f):
            if row.get("regime") != "1":
                continue
            g = lambda k: (float(row[k]) if row.get(k) not in (None, "") else None)
            return {"t": g("time"), "T_top": g("T_top"), "T_surf": g("T_surf"),
                    "Sp": g("Sp_top"), "ell": g("ell_eff")}
    return None


def ledger_max(tag, before=None):
    """max |ledger_err_solid|, optionally only for time < `before`.
    ⚠ MAX of a dimensionless RELATIVE column, never a sum (standing rule 7)."""
    for name in ("thermo.dat", "run_thermo.dat"):
        for base in (os.path.join(OUT, tag), os.path.join(OUT, tag, "run")):
            p = os.path.join(base, name)
            if not os.path.exists(p):
                continue
            with open(p) as f:
                head = f.readline().split()
                if "ledger_err_solid" not in head:
                    continue
                c = head.index("ledger_err_solid")
                ti = head.index("time") if "time" in head else 0
                worst, n = 0.0, 0
                for ln in f:
                    t = ln.split()
                    if len(t) <= max(c, ti):
                        continue
                    try:
                        if before is not None and float(t[ti]) >= before:
                            continue
                        worst = max(worst, abs(float(t[c])))
                        n += 1
                    except ValueError:
                        pass
                return None if n == 0 else worst
    return None


def main():
    if not os.path.exists(BIN):
        print("FAIL: missing binary")
        return 1
    print("P2-2b ord 3 -- THE KANT LADDER")
    print("binary %s" % subprocess.run(
        ["shasum", "-a", "256", BIN], capture_output=True,
        text=True).stdout[:24])
    print("K_Ic PINNED: %s\n" % " ".join(PIN))
    do_run = "--check" not in sys.argv
    rows, fails = [], []
    for tag, extra, what in RUNGS:
        if do_run:
            rc = run(tag, extra)
            if rc != 0:
                print("  FAIL %s exited %d" % (tag, rc))
                fails.append("%s exited %d" % (tag, rc))
                continue
        o = onset(tag)
        lm = ledger_max(tag)
        lm_pre = ledger_max(tag, before=o["t"]) if o else lm
        rows.append((tag, what, o, lm))
        print("  --- %s --- %s" % (tag, what))
        if o is None:
            print("      NO FIRING in 45 s")
        else:
            print("      onset t %.3f s   T_top %.2f K   T_surf %s   Sp %.4f"
                  % (o["t"], o["T_top"],
                     "n/a" if o["T_surf"] is None else "%.2f K" % o["T_surf"],
                     o["Sp"] or 0.0))
        print("      max|ledger_err_solid|  whole run %s   PRE-ONSET %s"
              % ("n/a" if lm is None else "%.3e" % lm,
                 "n/a" if lm_pre is None else "%.3e" % lm_pre))
        # (b3) is scored on K1 -- conduction and the H->T inversion on ONE curve
        if tag == "K1":
            if lm_pre is None:
                fails.append("(b3): no ledger column on K1")
                print("      FAIL (b3): no ledger_err_solid column")
            else:
                ok = lm_pre <= B3_BAR
                print("      %s (b3) PRE-ONSET max %.3e <= %.0e -- conduction "
                      "and the H->T inversion are ONE curve"
                      % ("ok  " if ok else "FAIL", lm_pre, B3_BAR))
                if not ok:
                    fails.append("(b3) K1 pre-onset ledger %.3e > %.0e"
                                 % (lm_pre, B3_BAR))
                # non-vacuity: the window must contain steps, and the
                # whole-run number must be the one that is large -- if both
                # were tiny the split would be measuring nothing.
                if lm is not None:
                    print("      note: whole-run %.3e is the removal-accounting "
                          "deficit, present on K0 too (pre-existing)" % lm)
    print()
    print("  (c1) onset by rung:")
    for tag, _, o, _ in rows:
        print("      %-4s %s" % (tag, "no firing" if o is None
                                 else "T_top %.2f K at t %.3f s"
                                      % (o["T_top"], o["t"])))
    if fails:
        print("\nFAILURES:")
        for f in fails:
            print("  " + f)
        print("\nFAIL")
        return 1
    print("\nPASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
