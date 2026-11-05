#!/usr/bin/env python3
"""P2-2d rev 6, ord 3 and 4 — the re-runs, the slopes, and THE DISCRIMINATOR.

⚠⚠ DELIVERABLE 7: RE-RUN, DO NOT RE-SCORE. Every leg whose size factor exceeded
1 fired LATE and its onset is wrong AS STORED, so no re-reduction of the old
files can fix it. All q = 300 legs plus the capped q = 500 are re-run here.

⚠⚠ THE CAPPED FAMILY'S PAIRING, DECLARED. rev 5 ran it with
`stress_relief = 0.531294` -- the FULL-expansion relief -- so it fired at
1032-1098 K, above the 873 K clamp, and is void. The two defensible pairings
are:
    coherent   : capped expansion WITH capped stiffness, R = 0.4648, 845.1 K
    incoherent : capped expansion with FULL stiffness,  R = 0.8050
⚠ The MODEL CAN ONLY EXPRESS THE INCOHERENT ONE -- `grimsel.alpha_cap_theta`
caps expansion and there is NO key that caps stiffness. So the leg below runs
the incoherent pairing at its OWN relief R = 0.8050 (not 0.531294), labelled as
a diagnostic; the coherent figure stays offline (derive_onset.py's (a7)).

⚠⚠ (c4) IS THE DISCRIMINATOR AND IT NEEDS A CONSTANT TOUGHNESS. `kic.law =
linear` with `kic.T_melt = 1e9` is a DECLARED DEVICE, not a melting
temperature: `KIcLinear` returns `K_Ic0*(1 - dT/T_melt)`, constant to seven
digits there. ⚠⚠ And the level then costs a FITTED parameter -- `K_Ic0` must be
set to land 850 K, which the at-temperature route did not need. Reported as a
0-fitted-parameter model against a 1-fitted-parameter one.

Usage: run_rev6.py
"""
import csv
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
sys.path.insert(0, HERE)
import derive_onset as D  # noqa: E402

BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
DECK = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "kant",
                    "sp_kant_onset_meiercfg", "input_model")
OUT = os.path.join(HERE, "rev6")
T_I = 293.15
FLUXES = (300.0, 500.0, 800.0)
R_FULL = 0.531294          # UCS/sigma(850), full expansion + full stiffness
R_INCOH = 0.8050           # capped expansion, FULL stiffness -- incoherent
CLAMP = 873.15

BASE = ["amr.n_cell=24 24 24", "amr.max_grid_size=24",
        "spall.surface_T=reconstructed", "stop_time=45.0",
        "amr.thermo.int=1", "amr.thermo.plot_int=1",
        "spall.removal_events_csv=1", "spall.flake_length=flux",
        "spallation.sp.kic.scale=1.0"]
WB = ["weibull.enabled=1", "weibull.a0_ig=8.0e-6", "weibull.a0_gb=40.0e-6",
      "weibull.size_volume=flux"]


def leg(tag, extra):
    d = os.path.join(OUT, tag)
    os.makedirs(d, exist_ok=True)
    argv = (["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4", BIN,
             DECK] + BASE + extra + ["plot_file=%s" % os.path.join(d, "run")])
    with open(os.path.join(d, "argv.sh"), "w") as f:
        f.write(" ".join(argv) + "\n")
    env = os.environ.copy()
    env.setdefault("OMPI_MCA_btl", "self,vader,tcp")
    with open(os.path.join(d, "run.log"), "w") as lg:
        rc = subprocess.call(argv, cwd=ROOT, env=env, stdout=lg,
                             stderr=subprocess.STDOUT)
    # ⚠⚠ "ABORTED" AND "DID NOT FIRE" ARE DIFFERENT ANSWERS. Returning None for
    # both made the constant-property ladder print "no fire" at all three
    # fluxes when every leg had in fact ABORTED on a missing required key. The
    # identical mistake produced an impossible bracket in P2-2b's fit_relief.
    if rc != 0:
        return "ABORT"
    p = os.path.join(d, "run_removal_events.csv")
    if not os.path.exists(p):
        return "NOFIRE"
    with open(p) as f:
        for row in csv.DictReader(f):
            if row.get("regime") == "1":
                g = lambda k: (float(row[k]) if row.get(k) not in (None, "")
                               else 0.0)
                return {"Ts": g("T_surf"), "rise": g("T_surf") - T_I,
                        "Sp": g("Sp_top"), "lq": g("ell_q")}
    return "NOFIRE"


def slope(qs, rs):
    pts = [(math.log(q), math.log(r)) for q, r in zip(qs, rs) if r and r > 0]
    if len(pts) < 2:
        return None
    n = len(pts)
    sx = sum(x for x, _ in pts); sy = sum(y for _, y in pts)
    sxx = sum(x * x for x, _ in pts); sxy = sum(x * y for x, y in pts)
    return (n * sxy - sx * sy) / (n * sxx - sx * sx)


def grimsel_family(m, R, extra, label):
    print("\n  --- %s, m = %g, R = %.6f ---" % (label, m, R))
    print("      %-7s %9s %9s %9s %10s" % ("q", "T_surf", "rise", "Sp_top",
                                           "ell_q[mm]"))
    got = []
    for q in FLUXES:
        r = leg("%s_m%g_q%.0f" % (label, m, q),
                ["surface_patch.h_conv=%.0f" % q,
                 "microstructure.phase0.law=grimsel", "spall.stress_law=grimsel",
                 "spallation.sp.kic.law=grimsel_table",
                 "spallation.sp.stress_relief=%.6f" % R,
                 "weibull.m=%g" % m, "weibull.V_ref=1.7188e-08"] + WB + extra)
        got.append(r)
        if isinstance(r, str):
            print("      %-7.0f %9s   <= NOT a physical result" % (q, r))
        else:
            print("      %-7.0f %9.2f %9.2f %9.4f %10.4f"
                  % (q, r["Ts"], r["rise"], r["Sp"], r["lq"] * 1e3))
    ok = [g for g in got if not isinstance(g, str) and g]
    rs = [g["rise"] if (g and not isinstance(g, str)) else None for g in got]
    lqs = [g["lq"] if (g and not isinstance(g, str)) else None for g in got]
    e = slope(FLUXES, rs)
    # the three factors, each measurable
    f1 = slope(FLUXES, lqs)
    print("      three-point slope e = %s" % ("n/a" if e is None
                                              else "%+.4f" % e))
    if f1 is not None:
        print("      factor 1  |dln l_q/dln q| = %.4f" % abs(f1))
        print("      factor 2  3/m             = %.4f" % (3.0 / m))
        print("      factor 3  |dln dT/dln a|  = %.4f (implied)"
              % (abs(e) / (abs(f1) * 3.0 / m) if e and f1 else float("nan")))
    if len(rs) == 3 and all(rs):
        s1 = math.log(rs[1] / rs[0]) / math.log(FLUXES[1] / FLUXES[0])
        s2 = math.log(rs[2] / rs[1]) / math.log(FLUXES[2] / FLUXES[1])
        print("      (c2) slope 300->500 %+.4f  vs  500->800 %+.4f  ratio %.2f"
              % (s1, s2, max(s1, s2) / min(s1, s2) if min(s1, s2) else 0))
    return e, got


def main():
    fails = []

    def need(ok, msg):
        print("  %s  %s" % ("ok  " if ok else "FAIL", msg))
        if not ok:
            fails.append(msg)

    print("P2-2d rev 6 -- ord 3 and 4")
    print("=" * 72)

    # ---- (c1)/(c2)/(c3): the Grimsel ladder, re-run ---------------------
    print("\n[ord 3] (c1)/(c2)/(c3) THE GRIMSEL LADDER, RE-RUN AFTER THE FIXES")
    res = {}
    for m in (8.0, 20.0):
        e, got = grimsel_family(m, R_FULL, [], "grimsel_full")
        res[m] = (e, got)
        sp = [g["Sp"] for g in got if g and not isinstance(g, str)]
        if sp:
            need(all(abs(v - 1.0) <= 0.01 for v in sp),
                 "(c1) m=%g: every Sp_top within 1.00 +/- 0.01 (%s) -- the two "
                 "gates agree and no slope is half-pinned"
                 % (m, ", ".join("%.4f" % v for v in sp)))

    # ---- (c6): the capped family, declared pairing ----------------------
    print("\n[ord 3] (c6) THE CAPPED FAMILY, RE-RUN WITH ITS OWN RELIEF")
    print("      pairing: capped expansion + FULL stiffness = INCOHERENT.")
    print("      R = %.4f (its own), NOT the 0.531294 rev 5 used." % R_INCOH)
    print("      ⚠ the COHERENT pairing (capped expansion + capped stiffness,")
    print("        R = 0.4648, 845.1 K) is NOT expressible: no key caps E.")
    _, capped = grimsel_family(8.0, R_INCOH,
                               ["microstructure.phase0.grimsel"
                                ".alpha_cap_theta=300"], "capped_incoh")
    cts = [g["Ts"] for g in capped if g and not isinstance(g, str)]
    if cts:
        need(max(cts) < CLAMP,
             "(c6) the capped family fires BELOW the 873 K clamp (max %.2f K) "
             "-- rev 5's 1032-1098 K was void" % max(cts))

    # ---- (c4)/(c5): THE DISCRIMINATOR ----------------------------------
    print("\n[ord 4] (c4) THE DISCRIMINATOR -- constant properties, draw ON")
    print("      law = constant: deck scalars E 30.0e9, beta 8.0e-6,")
    print("      kappa 1.5 (input_model:72-75) => sigma(850) 190.9 MPa,")
    print("      R = 0.8852.")
    print("      ⚠⚠ constant toughness via `kic.law = linear` +")
    print("         `kic.T_melt = 1e9` -- A DECLARED DEVICE, not a melting")
    print("         temperature. KIcLinear returns K_Ic0*(1 - dT/T_melt).")
    K_IC0_FIT = 0.951     # ⚠ FITTED to land the level; see (c5)
    print("      ⚠⚠ K_Ic0 = %.3f is FITTED to land ~850 K. The"
          % K_IC0_FIT)
    print("         at-temperature route needed NO fitted parameter, so this")
    print("         is a 1-fitted-parameter model against a 0-fitted one.")
    # ⚠⚠ THE ANCHOR IS PER-LADDER. `V_ref` = 1.7188e-08 was measured under the
    # PER-CELL lambda of the Grimsel ladder. Here lambda is the deck scalar 1.5,
    # so `ell_q` is different and reusing that anchor would put the size factor
    # far from 1 for a reason that has nothing to do with the modulus. Measured
    # on this ladder's OWN deterministic h_conv = 500 leg.
    a0 = leg("const_anchor", ["surface_patch.h_conv=500",
                              "spallation.sp.kic.law=linear",
                              "spallation.sp.kic.K_Ic0=%.6f" % K_IC0_FIT,
                              "spallation.sp.kic.T_ref=293.15",
                              "spallation.sp.kic.T_melt=1e9",
                              "spallation.sp.stress_relief=0.885200"])
    if isinstance(a0, str) or not a0 or not a0["lq"]:
        print("      FAIL: the constant-property anchor leg gave %s" % a0)
        return 1
    V_CONST = a0["lq"] ** 3
    print("\n      anchor (this ladder's own): ell_q = %.4f mm => V_ref = %.4e"
          % (a0["lq"] * 1e3, V_CONST))
    print("      ⚠ the Grimsel ladder's 1.7188e-08 is NOT reused -- that was")
    print("        measured under the per-cell lambda; here lambda = 1.5.")
    print("      (c5) LEVEL, deterministic flaw: %.2f K  %s   (target 850 K)"
          % (a0["Ts"], D.tag(a0["Ts"])))

    for m in (4.0, 8.0):
        print("\n  --- constant properties, m = %g ---" % m)
        print("      %-7s %9s %9s %9s %10s" % ("q", "T_surf", "rise", "Sp_top",
                                               "ell_q[mm]"))
        got = []
        for q in FLUXES:
            r = leg("const_m%g_q%.0f" % (m, q),
                    ["surface_patch.h_conv=%.0f" % q,
                     "spallation.sp.kic.law=linear",
                     "spallation.sp.kic.K_Ic0=%.6f" % K_IC0_FIT,
                     # ⚠ DEVICE: a constant K_Ic, not a melting temperature.
                     # ⚠ `linear` query_requires THREE keys, not two:
                     # K_Ic0, T_ref AND T_melt (:11115-11117). Omitting T_ref
                     # aborts, and the first harness read that as "no fire".
                     "spallation.sp.kic.T_ref=293.15",
                     "spallation.sp.kic.T_melt=1e9",
                     "spallation.sp.stress_relief=0.885200",
                     "weibull.m=%g" % m,
                     "weibull.V_ref=%.9e" % V_CONST] + WB)
            got.append(r)
            if isinstance(r, str):
                print("      %-7.0f %9s   <= NOT a physical result" % (q, r))
            else:
                print("      %-7.0f %9.2f %9.2f %9.4f %10.4f"
                      % (q, r["Ts"], r["rise"], r["Sp"], r["lq"] * 1e3))
        rs = [g["rise"] if (g and not isinstance(g, str)) else None for g in got]
        lqs = [g["lq"] if (g and not isinstance(g, str)) else None for g in got]
        e = slope(FLUXES, rs)
        f1 = slope(FLUXES, lqs)
        print("      slope e = %s" % ("n/a" if e is None else "%+.4f" % e))
        if f1 is not None:
            print("      factor 1 |dln l_q/dln q| = %.4f   factor 2 3/m = %.4f"
                  % (abs(f1), 3.0 / m))
            print("      ⚠ factor 1 was measured under the PER-CELL lambda in")
            print("        the Grimsel ladder; here lambda is the deck scalar")
            print("        1.5, so it does NOT carry over -- hence it is")
            print("        reported per ladder, not reused.")
        if got[1] and not isinstance(got[1], str):
            print("      (c5) level at q = 500: %.2f K  %s"
                  % (got[1]["Ts"], D.tag(got[1]["Ts"])))

    print("\n" + "=" * 72)
    if fails:
        print("FAILURES:")
        for f in fails:
            print("  " + f)
        print("\nFAIL")
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
