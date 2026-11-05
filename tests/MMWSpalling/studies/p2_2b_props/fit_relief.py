#!/usr/bin/env python3
"""P2-2b ord 3 deliverable 10 / (c2) / (c7): FIT R ON K2, FLAW HELD.

THE FIT. One parameter, `spallation.sp.stress_relief` (= R), bisected so K2's
onset temperature RISE lands on Kant Table 2's predicted band centre. dT(R) is
monotonically DECREASING in R -- a larger R means a larger driving stress and
therefore an earlier, cooler firing -- so bisection is well posed.

⚠⚠ THE FLAW IS HELD. `spallation.sp.a0` and the weibull radii are NOT touched.
P1's refit moved the flaw to 80 um against a REDUCED drive and was withdrawn in
P2-2a precisely because a flaw size matched to one drive and to no other is the
calibration-coupling defect this programme has now hit three times.

⚠⚠ K_Ic IS PINNED, AND WITHOUT THAT THE BAND TESTS NOTHING. `K_Ic(T)` and R are
EXACTLY DEGENERATE on Kant's single firing temperature: any onset reachable by
scaling one is reachable by scaling the other. `kic.law = nasseri_table` and
`kic.scale = 1.0` are stated on every argv.

(c7) RUNS THE FIT TWICE -- alpha as given, and alpha held at its 300 degC value
(`grimsel.alpha_cap_theta=300`). ⚠ If R leaves the (c2) band BETWEEN the two
fits, (c2) has no discriminating power: the quartz transition at 573 degC sits
between the Kant and Meier firing temperatures, so which alpha is right is not
settled by this test.

⚠ The cap key binds only because P2-2b parses the per-phase grimsel block
unconditionally. Inside the thermal-law branch it was silently ignored for a
mechanical-only leg, and both fits would have returned the SAME R.

Usage: fit_relief.py
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
OUT = os.path.join(HERE, "fit_relief")

T_I = 293.15
BAND_LO, BAND_HI = 390.0, 560.0          # Kant Table 2 predicted band
TARGET = 0.5 * (BAND_LO + BAND_HI)       # 475.0 K
R_BAND = (0.41, 0.69)                    # (c2)'s bar

PIN = ["spallation.sp.kic.law=nasseri_table", "spallation.sp.kic.scale=1.0"]
COMMON = ["amr.n_cell=24 24 24", "amr.max_grid_size=24",
          "surface_patch.h_conv=500", "spall.surface_T=reconstructed",
          "stop_time=45.0", "amr.thermo.int=1", "amr.thermo.plot_int=1",
          "spall.removal_events_csv=1",
          "microstructure.phase0.law=grimsel", "spall.stress_law=grimsel"]


def leg(tag, R, cap):
    d = os.path.join(OUT, tag)
    os.makedirs(d, exist_ok=True)
    extra = ["spallation.sp.stress_relief=%.6f" % R]
    if cap:
        extra.append("microstructure.phase0.grimsel.alpha_cap_theta=300")
    argv = (["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4", BIN,
             DECK] + COMMON + PIN + extra +
            ["plot_file=%s" % os.path.join(d, "run")])
    with open(os.path.join(d, "argv.sh"), "w") as f:
        f.write(" ".join(argv) + "\n")
    env = os.environ.copy()
    env.setdefault("OMPI_MCA_btl", "self,vader,tcp")
    with open(os.path.join(d, "run.log"), "w") as lg:
        rc = subprocess.call(argv, cwd=ROOT, env=env, stdout=lg,
                             stderr=subprocess.STDOUT)
    # ⚠⚠ "ABORT" and "did not fire" are DIFFERENT ANSWERS and must not share a
    # return value. The first version returned None for both; R = 1.5 aborted
    # (R is constrained to (0,1] -- it RELIEVES), the harness read that as
    # "no fire", and the bracket it then reported contradicted K2's own R = 1
    # result without anything noticing.
    if rc != 0:
        return "ABORT"
    p = os.path.join(d, "run_removal_events.csv")
    if not os.path.exists(p):
        return "NOFIRE"
    with open(p) as f:
        for row in csv.DictReader(f):
            if row.get("regime") == "1":
                return float(row["T_top"]) - T_I
    return "NOFIRE"


def fit(label, cap):
    print("\n  === FIT: %s ===" % label)
    print("      target dT %.1f K (Kant band %.0f-%.0f centre)"
          % (TARGET, BAND_LO, BAND_HI))
    # ⚠ R is constrained to (0, 1] by `ParseSpallSettings` -- it RELIEVES the
    # stress, so a value above 1 would amplify it and the code aborts. The
    # bracket's upper end is therefore the LARGEST LEGAL R, not a free choice.
    lo, hi = 0.10, 1.00
    n = 0
    hist = []
    # bracket check first -- a fit inside an unverified bracket is meaningless
    d_hi = leg("%s_hi" % label, hi, cap)
    d_lo = leg("%s_lo" % label, lo, cap)
    fmt = lambda d: d if isinstance(d, str) else "%.2f" % d
    print("      bracket: R=%.2f -> dT %s   R=%.2f -> dT %s"
          % (hi, fmt(d_hi), lo, fmt(d_lo)))
    if d_hi == "ABORT" or d_lo == "ABORT":
        print("      FAIL bracket: a leg ABORTED -- not a physical result")
        return None, hist
    if not isinstance(d_hi, float) or d_hi >= TARGET:
        print("      FAIL bracket: at the largest legal R the onset is still "
              "at or above the target, so no R in (0,1] can reach it")
        return None, hist
    if isinstance(d_lo, float) and d_lo <= TARGET:
        print("      FAIL bracket: even R=%.2f fires below the target" % lo)
        return None, hist
    while n < 10 and (hi - lo) > 0.005:
        mid = 0.5 * (lo + hi)
        d = leg("%s_%02d" % (label, n), mid, cap)
        print("      R %.4f -> dT %s" % (mid, fmt(d)))
        hist.append((mid, d))
        if d == "ABORT":
            print("      FAIL: leg aborted mid-bisection")
            return None, hist
        # NOFIRE means the stress never reached threshold in 45 s, i.e. the
        # onset is beyond the window -- same side of the bracket as a dT above
        # the target, and labelled so rather than silently merged.
        if d == "NOFIRE" or d > TARGET:
            lo = mid           # too little stress, raise R
        else:
            hi = mid
        n += 1
    R = 0.5 * (lo + hi)
    print("      FITTED R = %.4f" % R)
    return R, hist


def main():
    if not os.path.exists(BIN):
        print("FAIL: missing binary")
        return 1
    print("P2-2b ord 3 deliverable 10 -- R FITTED ON K2, FLAW HELD")
    print("K_Ic PINNED: %s" % " ".join(PIN))
    print("flaw: spallation.sp.a0 and weibull.* UNTOUCHED")
    fails = []
    R_given, _ = fit("alpha_given", False)
    R_cap, _ = fit("alpha_cap300", True)

    print("\n  (c2) / (c7) RESULT")
    for nm, R in (("alpha as given", R_given), ("alpha capped at 300 C", R_cap)):
        if R is None:
            print("      %-24s FIT FAILED" % nm)
            fails.append("%s: fit failed" % nm)
            continue
        ok = R_BAND[0] <= R <= R_BAND[1]
        print("      %-24s R = %.4f   %s (c2) band %.2f-%.2f"
              % (nm, R, "IN" if ok else "OUT OF", *R_BAND))
        if not ok:
            fails.append("%s: R %.4f outside %.2f-%.2f" % (nm, R, *R_BAND))
    if R_given is not None and R_cap is not None:
        span = abs(R_given - R_cap)
        print("      span between the two fits  %.4f" % span)
        inside = all(R_BAND[0] <= r <= R_BAND[1] for r in (R_given, R_cap))
        print("      (c7) verdict: %s"
              % ("both fits inside the band -- (c2) HAS power" if inside else
                 "R leaves the band between the two alpha choices -- (c2) has "
                 "NO power and which alpha is right is not settled here"))
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
