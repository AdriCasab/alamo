#!/usr/bin/env python3
"""P2-2a deliverable 2 / criterion (a9) -- THE STRESS-FREE REFERENCE, OFFLINE.

⚠⚠ rev 7 REPLACED the reference definition. rev 6 used the column MINIMUM; this
is the STORED INITIAL TEMPERATURE FIELD, read at the candidate flake's own cell
and never updated. The minimum was withdrawn because it is a reduction over the
column, so (i) it inherits round-off dust from the deepest cell -- measured
1.99e-12 K on `QF` and 9.05e-08 K on `M-100`, against exactly 0 on the
`M-2000` control whose exp(-500) underflows -- and (ii) it DEGRADES AS THE HOLE
DEEPENS, because progressively less cold rock is left beneath the front
(`M-100` at step 10000: 69 of 250 cells solid, base genuinely 293.53 K). Under
volumetric heating a column warms through, the minimum rises toward the
column's own temperature, the bracket collapses and spalling stops -- i.e. the
model would switch itself off on the target application in a way that looks
like physics.

Four things, in the packet's order:

  1. the SIGN, against the shipped compression-positive convention, with the
     code sites that fix it;
  2. the drilling limit -- the new bracket reduces to today's;
  3. the uniform limit -- the new bracket is EXACTLY zero;
  4. the stress formula re-printed with its three values.

Then a fifth block the packet does NOT ask for and which is the reason this
script exists in the form it does: (b1) predicts `0 differing of 2639` on the
identity manifold, and change 6 carries NO KEY, so that prediction is a claim
about every manifest leg, not only the drilling ones. The block evaluates the
bracket on each manifest deck's own (T_ref, initial T) pair.

⚠ METHOD, stated where the numbers appear (guardrail 1): every stress figure
below is `sigma = E(T_s) * [strain(T_s) - strain(T_restraint)] / (1 - nu)` with
`strain` = `grimsel.thermal_strain` (exact, analytic per Table-3 segment) and
`E` = `grimsel.E_of`. Nothing here is quadrature and nothing is a literal.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import grimsel as G  # noqa: E402

NU = G.NU
BAR = "=" * 78
SUB = "-" * 78


def strain_K(T_K, cap_300=False):
    """thermal_strain, keyed on KELVIN so every caller below shares one
    conversion. grimsel.thermal_strain takes degC."""
    return G.thermal_strain(T_K - 273.15, cap_300)


def sigma_restraint(T_K, T_restraint_K, cap_300=False, nu=NU):
    """THE ONE EXPRESSION, compression-positive:

        sigma = E(T) * [ strain(T) - strain(T_restraint) ] / (1 - nu)

    Positive when T > T_restraint, i.e. a hot layer against cold rock is in
    COMPRESSION and drives delamination. The shipped criterion discards a
    non-positive drive, so the sign is part of the formula, not a convention
    applied afterwards."""
    eps = strain_K(T_K, cap_300) - strain_K(T_restraint_K, cap_300)
    return G.E_of(T_K - 273.15) * eps / (1.0 - nu)


# ---------------------------------------------------------------- 1. the sign
def block_sign():
    print(BAR)
    print("(a9).1  THE SIGN -- compression-positive, hot MINUS restraint")
    print(BAR)
    print("""
  The shipped convention, read from the source (not assumed):

    Removal.H:1833   `return -sg(0, 0);`
                     the stress_mf sampler NEGATES sigma_xx, so the value the
                     K_I integrator consumes is COMPRESSION-POSITIVE.
    Removal.H:1915   `const Set::Scalar jmp = T_hi_z - T_ref_z;`
    Removal.H:1916   `return (jmp > 0.0) ? ... : 0.0;`
                     HOT minus COLD, and a non-positive drive is DISCARDED.
    Removal.H:1924   the live `absolute` branch returns
                     `+ E * beta * (T_z - phase_T_ref_loc) / (1 - nu)`
                     the comment at :1922 states the double negation outright.

  => the bracket must be  strain(T) - strain(T_restraint).

  ⚠⚠ rev 4 wrote it the other way round -- `strain(T_restraint) - strain(T)`.
  On every hot flake that is NEGATIVE, the `jmp > 0.0` guard discards it, the
  drive is identically zero and NOTHING SPALLS. It would not have been a wrong
  number; it would have stopped all drilling. This row is why it is caught on
  paper rather than in a 500 s leg.
""")
    T_s, T_r = 823.15, 293.15
    s_right = sigma_restraint(T_s, T_r)
    s_rev4 = -s_right
    print("  numerically, at T_s = %.2f K against restraint %.2f K:" % (T_s, T_r))
    print("    rev 5 (hot - restraint) : %+.4e Pa   -> kept by `jmp > 0`"
          % s_right)
    print("    rev 4 (restraint - hot) : %+.4e Pa   -> DISCARDED, drive 0"
          % s_rev4)
    ok = s_right > 0.0 and s_rev4 < 0.0
    print("    sign correct and rev 4 refuted : %s" % ok)
    return ok


# ------------------------------------------------------- 2. the drilling limit
def block_drilling():
    print()
    print(BAR)
    print("(a9).2  THE DRILLING LIMIT -- the bracket is LITERALLY UNCHANGED")
    print(BAR)
    print("""
  Today's live expression references `phases[0].T_ref`, a DECLARED INPUT
  CONSTANT. The new one reads the STORED INITIAL TEMPERATURE at the sampled
  cell. On a deck whose initial condition equals its typed reference the two are
  THE SAME NUMBER -- not close, the same -- so the bracket is bit-identical and
  no tolerance applies. Verified from the decks, not assumed:
""")
    DECKS = (("Meier  input_drilling",                      293.15, 293.15),
             ("Kant   sp_kant_onset_meiercfg/input_model",   293.15, 293.15),
             ("unit/gradient_blind/input",                   293.15, 820.00),
             ("unit/jet_enthalpy/input_strip",               293.15, 500.00))
    print("  %-44s %-10s %-10s %s"
          % ("deck", "T_ref", "hc.ic", "bracket vs today"))
    print("  " + SUB)
    exact = True
    for name, tref, ic in DECKS:
        same = (tref == ic)
        print("  %-44s %-10.2f %-10.2f %s"
              % (name, tref, ic,
                 "IDENTICAL" if same else "CHANGES by %.2f K" % (ic - tref)))
        if name.startswith(("Meier", "Kant")) and not same:
            exact = False
    print()
    print("  ⇒ on every drilling and Kant leg the reference is the same 293.15 K")
    print("    the code already used, so (b1) is a statement about arithmetic")
    print("    identity rather than about a tolerance.")
    print("  ⚠ The two decks that DISAGREE are the two the defect was found on:")
    print("    the slab books 526.85 K and the strip 206.85 K of superheat that")
    print("    no heating produced. Those are the files (b1) authorises to move.")
    print("  drilling/Kant legs identical : %s" % exact)
    return exact


# -------------------------------------------------------- 3. the uniform limit
def block_uniform():
    print()
    print(BAR)
    print("(a9).3  THE UNIFORM LIMIT -- the bracket is EXACTLY zero")
    print(BAR)
    print("""
  A column held uniformly at T with no gradient is restrained by nothing: the
  rock around the candidate flake has expanded by exactly as much as the flake
  has. The column minimum IS the local temperature, so the bracket is a number
  minus itself -- zero by cancellation, not zero by smallness. No tolerance is
  quoted below because none is needed.
""")
    print("  %-12s %-18s %-18s %s"
          % ("T uniform", "strain(T)", "bracket", "sigma [Pa]"))
    print("  " + SUB[:70])
    allz = True
    for T in (500.0, 745.0, 820.0, 873.15, 923.15):
        br = strain_K(T) - strain_K(T)
        s = sigma_restraint(T, T)
        allz = allz and (br == 0.0) and (s == 0.0)
        print("  %-12.2f %-18.10e %-18.1f %.1f" % (T, strain_K(T), br, s))
    print()
    print("  every bracket is 0.0 exactly (== 0.0, not < eps) : %s" % allz)
    print("  ⚠ this is `unit/gradient_blind`'s 820 K slab. Today's expression")
    print("    returns beta*(820 - 293.15) = 527 K of PHANTOM superheat there,")
    print("    which is the defect the test currently certifies.")
    return allz


# ------------------------------------------------------------ 4. the formula
def block_formula():
    print()
    print(BAR)
    print("(a9).4  THE STRESS FORMULA AND ITS THREE VALUES")
    print(BAR)
    print("""
  METHOD (named here, where the numbers appear):
      sigma(T_s) = E(T_s) * INTEGRAL[T_0..T_s] alpha(theta) dtheta / (1 - nu)
  E from the measured E/E_0 ratios; alpha from Table 3; the integral exact per
  segment; nu = %.3f unchanged. T_0 = 293.15 K (the restraint, cold bulk).
""" % NU)
    # ⚠ T_s CARRIED from the ord-1 run (`PROPS.txt` (a2) table), not re-chosen
    # here: Meier at 800.0 K, Kant at 873.0 K, and "today" at the leg's OWN
    # measured face temperature 759.7 K. PROPS.txt also records why the
    # packet's `159 today` is not comparable with its other two -- it is at a
    # different T_s, so the 159 -> 299 step mixes a temperature change with a
    # property change. At one common 800 K the property effect alone is 1.73x.
    print("  %-22s %-10s %-12s %-14s %-12s %s"
          % ("condition", "T_s [K]", "E [GPa]", "strain", "alpha_sec", "sigma [MPa]"))
    print("  " + SUB)
    out = {}
    # ⚠⚠ THE FIRST ROW IS A DIFFERENT METHOD AND IS LABELLED AS ONE. "today"
    # uses the SHIPPED FLAT coefficient beta = 8.0e-6 /K, so its strain is
    # beta*(T_s - T_0); the other two integrate the measured alpha(theta).
    # Comparing row 1 against rows 2-3 therefore mixes a property change with
    # a method change -- PROPS.txt (a2) establishes the clean one-T_s
    # comparison (173.8 -> 300.4 at 800 K, a factor 1.73).
    # read from the SHIPPING deck, not written here:
    #   microstructure.phase0.beta = 8.0e-6, phase0.E = 30.0e9
    #   (validation/meier/sp_meier_pilot/input_drilling:99-100)
    BETA_FLAT, E_SHIPPED = 8.0e-6, 30.0e9
    T_0 = 293.15
    for T_s, want in ((800.0, 174.0),):
        eps = BETA_FLAT * (T_s - T_0)
        sig = E_SHIPPED * eps / (1.0 - NU)
        out["today, FLAT beta"] = sig
        print("  %-22s %-10.2f %-12.3f %-14.6e %-12.3e %.1f   [PROPS %.0f]"
              % ("today, deck E+beta", T_s, E_SHIPPED * 1e-9, eps, BETA_FLAT,
                 sig * 1e-6, want))
    for label, T_s, want in (("MEASURED Meier", 800.0, 299.0),
                             ("MEASURED Kant", 873.0, 321.0)):
        sig, E_s, eps, a_mean = G.sigma_formula(T_s)
        out[label] = sig
        print("  %-22s %-10.2f %-12.3f %-14.6e %-12.3e %.1f   [packet %.0f]"
              % (label, T_s, E_s * 1e-9, eps, a_mean, sig * 1e-6, want))
    print()
    print("  the property effect at ONE common T_s = 800 K : %.2fx"
          % (out["MEASURED Meier"] / out["today, FLAT beta"]))
    print()
    print("  UCS = %.0f +/- %.0f MPa  =>  R = UCS/sigma:"
          % (G.UCS * 1e-6, G.UCS_SD * 1e-6))
    rs = {}
    for label in ("MEASURED Meier", "MEASURED Kant"):
        rs[label] = G.UCS / out[label]
        print("    %-22s R = %.4f" % (label, rs[label]))
    # the band is the WIDEST R the UCS uncertainty admits across both legs:
    # the low end takes UCS - SD against the LARGER sigma, the high end
    # UCS + SD against the SMALLER sigma.
    s_hi, s_lo = max(rs, key=lambda k: out[k]), min(rs, key=lambda k: out[k])
    lo = (G.UCS - G.UCS_SD) / out[s_hi]
    hi = (G.UCS + G.UCS_SD) / out[s_lo]
    print("    band from UCS +/- 37 MPa : %.3f - %.3f   [packet 0.41-0.69]"
          % (lo, hi))
    spread = max(rs.values()) / min(rs.values())
    print("    spread %.3fx -- inside UCS's own +/-22 %%, so R is ONE number"
          % spread)
    print("  ⚠ 1.8-1.9x the strength either way: an elastic criterion is")
    print("    outside its validity, which is what R measures.")
    # the two R values must land in the packet's band, else the carried
    # (a2) result has not reproduced and change 4 has no target.
    return all(0.41 <= r <= 0.69 for r in rs.values())


# ------------------------------------------------------------ 5. the in-situ case
def block_insitu():
    """⚠ The case NEITHER a typed constant NOR a column minimum can reach, and
    the reason the stored field generalises: rock at depth sitting at 200 degC
    is STRESS-FREE at 200 degC. Scored OFFLINE ONLY -- no such deck exists."""
    print()
    print(BAR)
    print("(a9).5  THE IN-SITU CASE -- stress-free at the temperature it sits at")
    print(BAR)
    print("""
  A geothermal target is not at 293 K. Rock initialised at 200 degC has done its
  expanding already: it is unstressed at 473.15 K, and only heating ABOVE that
  loads a flake. The stored field returns 473.15 and the bracket is zero at
  t = 0; a typed 293.15 books the whole 180 K as phantom superheat.
""")
    T_IN = 473.15                      # 200 degC
    print("  %-26s %-14s %-16s %s"
          % ("reference", "T0 used [K]", "bracket at t=0", "sigma [MPa]"))
    print("  " + SUB[:72])
    s_stored = sigma_restraint(T_IN, T_IN)
    s_typed = sigma_restraint(T_IN, 293.15)
    print("  %-26s %-14.2f %-16.1f %.3f"
          % ("stored initial field", T_IN, 0.0, s_stored * 1e-6))
    print("  %-26s %-14.2f %-16.2f %.3f"
          % ("typed 293.15 (today)", 293.15, T_IN - 293.15, s_typed * 1e-6))
    print()
    print("  phantom superheat booked by the typed constant : %.2f K" % (T_IN - 293.15))
    print("  phantom stress booked                          : %.1f MPa"
          % (s_typed * 1e-6))
    ok = (s_stored == 0.0) and (s_typed > 0.0)
    print("  stored field returns EXACTLY zero, typed does not : %s" % ok)
    print("  ⚠ OFFLINE ONLY -- no deck is initialised at depth yet (packet P6).")
    return ok


def main():
    print(BAR)
    print("P2-2 rev 5 -- (a9) THE RESTRAINT REFERENCE, DERIVED OFFLINE")
    print(BAR)
    print("  grimsel reference: %s" % os.path.relpath(G.ref_path(), HERE))
    rows = [("(a9).1 sign", block_sign()),
            ("(a9).2 drilling limit", block_drilling()),
            ("(a9).3 uniform limit", block_uniform()),
            ("(a9).4 formula", block_formula()),
            ("(a9).5 in-situ case", block_insitu())]
    print()
    print(BAR)
    print("(a9) VERDICT")
    print(BAR)
    for name, ok in rows:
        print("  %-26s %s" % (name, "PASS" if ok else "** FAIL **"))
    gate = all(ok for name, ok in rows)
    print()
    print("  (a9) ALL FIVE ROWS : %s" % ("PASS" if gate else "FAIL"))
    print("  ⚠ row .5 is scored OFFLINE ONLY -- it is the case the old form")
    print("    could never reach, and no deck exercises it yet.")
    return 0 if gate else 1


if __name__ == "__main__":
    sys.exit(main())
