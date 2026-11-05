#!/usr/bin/env python3
"""P2-2d rev 6 — deliverable 3 / (a3) / (a4) / (c5)'s UCS sensitivity.

⚠⚠ THE PROVENANCE GAP THIS CLOSES. The 7.73 um headline was quoted for several
revisions with no script that printed it. Standing rule 8: a number is backed
only if both the file it was read from and the code that reduced it are in the
repository. This reads the ARCHIVED prediction leg and prints it.

Usage: reduce_level.py
"""
import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import derive_onset as D  # noqa: E402

A_HELD = 8.0e-6
T_REF = 293.15


def onset_from(path):
    if not os.path.exists(path):
        return None
    with open(path) as f:
        for row in csv.DictReader(f):
            if row.get("regime") == "1":
                return float(row["T_surf"])
    return None


def main():
    print("P2-2d rev 6 -- THE LEVEL, REDUCED FROM THE ARCHIVED FILES")
    print("=" * 68)
    src = os.path.join(HERE, "prediction", "P_grimsel",
                       "run_removal_events.csv")
    Ts = onset_from(src)
    print("\n[a3] THE IMPLIED-FLAW REDUCER")
    print("  source: %s" % os.path.relpath(src, HERE))
    if Ts is None:
        print("  FAIL: archived prediction leg not found")
        return 1
    print("  onset surface temperature, as stored: %.3f K" % Ts)
    # a = pi*[K_Ic(T)/(W*UCS*s_ratio)]^2 -- the collapsed condition inverted,
    # with R = UCS/sigma(850) so R*sigma(T) = UCS*sigma(T)/sigma(850).
    sr = D.sigma(Ts, False, False) / D.sigma(850.0, False, False)
    k = D.kic_grimsel(Ts)
    a = math.pi * (k * 1e6 / (D.W * D.UCS0 * sr)) ** 2
    print("  K_Ic(T_s)  = %.4f MPa.sqrt(m)   (Grimsel::K_Ic_of, via KIC_CPP)" % k)
    print("  sigma ratio= %.4f               (sigma(T_s)/sigma(850))" % sr)
    print("  W          = %.6f               (DERIVED Tada integral)" % D.W)
    print("  => implied flaw a = %.4f um  against the held %.2f um"
          % (a * 1e6, A_HELD * 1e6))
    print("     ratio to held: %.4f" % (a / A_HELD))
    # ⚠⚠ TWO DIFFERENT QUANTITIES HAVE BOTH BEEN CALLED "the implied flaw",
    # AND THE PACKET'S 7.73 um IS NOT THE ONE THE ARCHIVED RUN IMPLIES.
    #   (i)  the flaw implied by the TABLE at the 850 K target -- a property of
    #        K_Ic(T) and UCS, needing no run at all;
    #   (ii) the flaw implied by WHERE THE MODEL ACTUALLY FIRED (843.757 K),
    #        which is what a reducer over the archived files produces.
    # They differ because K_Ic(T) FALLS with T, so evaluating 6.2 K lower gives
    # a higher toughness (0.9441 vs 0.9347) and hence a larger flaw. Both are
    # meaningful; only (ii) is a reduction of the archived leg, and deliverable
    # 3 asks for a reducer over those files.
    a850 = math.pi * (D.kic_grimsel(850.0) * 1e6 / (D.W * D.UCS0)) ** 2
    print("\n  the two quantities, separated:")
    print("    (i)  table at the 850 K TARGET      : %.3f um   <- the packet's"
          % (a850 * 1e6))
    print("                                            7.73, reproduced")
    print("    (ii) implied by the ARCHIVED FIRING : %.3f um   <- this reducer"
          % (a * 1e6))
    print("                                            at %.3f K" % Ts)
    ok_i = abs(a850 * 1e6 - 7.73) <= 0.02
    print("  %s (a3): the packet's 7.73 um reproduces as %.3f um (table at 850 K)"
          % ("ok  " if ok_i else "FAIL", a850 * 1e6))
    print("  ⚠ the run-implied value is %.3f um, ratio to the held 8.00 of"
          % (a * 1e6))
    print("    %.4f -- CLOSER to the held flaw than the 850 K figure's %.4f,"
          % (a / A_HELD, a850 / A_HELD))
    print("    so the agreement is better, not worse. But it sits OUTSIDE the")
    print("    rev-5 (c2) window of 7.2-7.8 um, which was set from quantity (i)")
    print("    while being described as the implied flaw of the run.")
    ok = ok_i

    print("\n[a4] THE SENSITIVITY AT THE FIRING TEMPERATURE")
    # |dln(dT)/dln a| from the collapsed condition: K_Ic(T) = W*R*sigma(T)*sqrt(a/pi)
    # => dln a = 2*[dln K_Ic - dln sigma]; with s = dln sigma/dln(dT) and
    #    kappa = dln K_Ic/dln(dT):  dln(dT)/dln a = 1/(2*(kappa - s)).
    print("  from the collapsed condition, |dln(dT)/dln a| = 1/(2|kappa - s|)")
    print("  %-14s %9s %9s %12s %s" % ("evaluated at", "s", "kappa",
                                       "sensitivity", "formula e (m=8)"))
    for T, note in ((806.8, "where they fire"), (850.0, "the old error")):
        dT = 0.5
        s = (math.log(D.sigma(T + dT, False, False))
             - math.log(D.sigma(T - dT, False, False))) \
            / (math.log(T + dT - T_REF) - math.log(T - dT - T_REF))
        kap = (math.log(D.kic_grimsel(T + dT)) - math.log(D.kic_grimsel(T - dT))) \
            / (math.log(T + dT - T_REF) - math.log(T - dT - T_REF))
        sens = 1.0 / (2.0 * abs(kap - s))
        e = (1.5 / 8) / (s + 1.5 / 8 - kap)
        print("  %-14s %+9.3f %+9.3f %12.3f %13.4f   (%s)"
              % ("%.1f K" % T, s, kap, sens, e, note))
    print("  ⚠ the planner evaluated this at 850 K for four revisions; the")
    print("    full legs fire at 803-824 K.")

    print("\n[c5] THE CONSTANT-PROPERTY LEVEL'S SENSITIVITY TO UCS")
    # ⚠⚠ the collapsed condition makes the level depend only on UCS and the
    # flaw, so a +/-1 SD band on UCS is the honest error bar on that level --
    # and K_Ic0 was FITTED to land it, so the band is what the fit absorbs.
    print("  the collapsed condition is K_Ic_req = W*UCS*sqrt(a/pi): the level")
    print("  depends on UCS and the flaw ALONE, so the fitted K_Ic0 absorbs")
    print("  whatever UCS is wrong by.")
    print("  %-14s %16s %s" % ("UCS [MPa]", "K_Ic req [MPa.m^.5]", "vs fitted 0.951"))
    for u in (D.UCS0 - D.UCS_SD, D.UCS0, D.UCS0 + D.UCS_SD):
        kq = D.W * u * math.sqrt(A_HELD / math.pi) / 1e6
        print("  %-14.0f %16.4f %+.1f %%" % (u / 1e6, kq, 100 * (kq / 0.951 - 1)))
    print("  => one SD of UCS moves the required toughness by -22 %% / +22 %%,")
    print("     so 'K_Ic0 = 0.951' is only as meaningful as UCS = 169 +/- 37.")
    print("\n" + "=" * 68)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
