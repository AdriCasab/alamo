#!/usr/bin/env python3
"""P2-2d rev 6 — deliverables 1 and 2 / (a1) and (a2), both BLOCKING.

Does a measured onset-vs-flux series exist for this rock at all? Everything
below is derived against the files in this repository, not quoted.

Usage: target_audit.py
"""
import math
import os
import re
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "..", "..", "..", ".."))
KANT = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "kant",
                    "kant_von_rohr_validation_test_for_mmwspalling_plan.txt")
MD = os.path.join(ROOT, "Claude_markdowns", "2026-10-06b.md")


def main():
    fails = []

    def need(ok, msg):
        print("  %s  %s" % ("ok  " if ok else "FAIL", msg))
        if not ok:
            fails.append(msg)

    print("P2-2d rev 6 -- (a1)/(a2) THE TARGET AUDIT")
    print("=" * 70)

    lines = open(KANT).read().splitlines()
    print("\n[a1] IS THERE A MEASURED ONSET-VS-FLUX SERIES FOR THIS ROCK?")
    print("  file: %s" % os.path.relpath(KANT, ROOT))
    print("  length: %d lines" % len(lines))
    need(len(lines) == 151, "the file is 151 lines")
    # the citations five revisions used
    for cite in (251, 253, 258):
        need(cite > len(lines),
             "citation [%d] points PAST the end of a %d-line file"
             % (cite, len(lines)))
    hits = [(i + 1, l.strip()) for i, l in enumerate(lines)
            if re.search(r"0\.229|0\.306|exponen", l, re.I)]
    for n, l in hits:
        print("       ** %d: %s" % (n, l[:70]))
    need(not hits,
         "grep for `0.229`, `0.306`, `exponen` in the Kant file returns "
         "NOTHING (%d hits)" % len(hits))

    print("\n  WHAT THE FILE DOES CARRY (derived):")
    for i, l in enumerate(lines):
        if re.search(r"theta_SP|theta_FP|Minimal required surface", l):
            print("    %3d: %s" % (i + 1, l.strip()[:72]))
    print("  => single MEANS, not a flux series.")

    print("\n  WHERE THE RETIRED NUMBERS ACTUALLY COME FROM:")
    if os.path.exists(MD):
        for i, l in enumerate(open(MD).read().splitlines(), 1):
            if "0.306" in l or "0.229" in l:
                print("    %s:%d  %s" % (os.path.relpath(MD, ROOT), i,
                                         l.strip()[:76]))
    print("  => WILKINSON's, on BARRE and WESTERLY granite, and already the")
    print("     OUTPUT of a Weibull size-effect fit. Scoring our Weibull route")
    print("     against them is partly circular, and the rock is wrong.")
    print("\n  AND THE SECOND MISATTRIBUTION: Table 4's temperature column is")
    print("     `theta_s,sens` -- the surface temperature of the HEAT-FLUX")
    print("     SENSOR PLATE that replaced the rock during burner calibration.")
    print("     Fitting it measures how a metal plate warms. VOID as a target.")

    # the resolvability bound
    q_lo, q_hi = 450.0, 866.0
    unc = 0.055                      # 5-6 % on T, per the file
    bound = unc / math.log(q_hi / q_lo)
    print("\n  THE RESOLVABILITY BOUND, DERIVED:")
    print("    flux range %.0f-%.0f kW/m2 => ratio %.3fx" % (q_lo, q_hi, q_hi / q_lo))
    print("    surface-T uncertainty ~5-6 %%")
    print("    => any exponent fitted across it is bounded by +/- %.3f" % bound)
    need(abs(bound - 0.08) < 0.02,
         "the bound is +/-%.2f, comparable to the whole model-target gap"
         % bound)

    print("\n  ⇒ (a1) ANSWER: NO. There is no measured onset-vs-flux series for")
    print("    this rock on disk. EVERY exponent row in this packet is an")
    print("    INTERNAL CONSISTENCY TEST between two property ladders, NOT a")
    print("    validation, and is labelled so.")

    print("\n[a2] WHAT WOULD SUPPLY A REAL SERIES, AND WHAT IT WOULD LICENSE")
    print("  (1) Kant & von Rohr's 2016/2017 papers, IF they tabulate per-run")
    print("      onset temperature by standoff. Licenses: a same-rock,")
    print("      same-apparatus onset-vs-flux comparison.")
    print("      ⚠ CONFOUND: Kant varies STANDOFF, so q and T_flame co-vary")
    print("        (1450 -> 1100 degC across the table), while the model swept")
    print("        `h_conv` at FIXED gas temperature. The two ladders are not")
    print("        the same experiment and a slope comparison needs the model")
    print("        re-run along Kant's own standoff locus to be like-for-like.")
    print("      ⚠ `q_sens` carries ~17 %% sensor calibration uncertainty and is")
    print("        an incident COLD-FACE flux, not absorbed flux at a hot rock.")
    print("  (2) Wilkinson & Tester's raw Barre/Westerly points, as an")
    print("      explicitly DIFFERENT-ROCK reference. Licenses: an order-of-")
    print("      magnitude sanity check only.")
    print("      ⚠ Those exponents are ALREADY a Weibull-model fit output, so")
    print("        scoring our Weibull route on them is partly circular; and")
    print("        REFERENCES.md states three of its four primary sources are")
    print("        PAYWALLED AND NEVER FETCHED.")
    print("  ⇒ neither is on disk today. Until one is, no exponent row here")
    print("    may be called a validation.")

    print("\n" + "=" * 70)
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
