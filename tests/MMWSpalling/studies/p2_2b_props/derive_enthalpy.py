#!/usr/bin/env python3
"""P2-2b ord 1 -- OFFLINE, zero compute. Criteria (a1), (a2), (a3).

  (a1)  the hard-coded `rho*c_p` enumeration, re-derived against the CURRENT
        source rather than quoted, and tested for completeness.
  (a2)  H(T) and T(H), `Material::Constant` from the scalars against Grimsel
        from the laws, 293-1200 K, max divergence WITH ITS SIGN.
  (a3)  (d2)'s arithmetic written out: numerator, denominator, the
        quarter-domain factor, and the value on today's baseline.

⚠ METHOD, named where the numbers appear. The enthalpy both ways is
`H(T) = INTEGRAL[T_ref..T] effective_rhocp(theta) dtheta` with
`effective_rhocp = rho(theta) * [c_p(theta) + latent boxcars]`
(`src/Numeric/Material/Material.H:30-43`), integrated the way
`Table::Build` integrates it -- trapezoid per grid interval
(`Table.H:52-54`) -- so this is the same quadrature the code uses and not a
sharper one that would flatter the comparison.
"""

import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "../../../.."))
sys.path.insert(0, os.path.join(HERE, "..", "p2_2_props"))

import grimsel as G  # noqa: E402

BAR = "=" * 78
SUB = "-" * 78

# the shipped deck scalars (validation/meier/sp_meier_pilot/input_drilling:102-106)
DECK_RHO, DECK_CP, DECK_TREF = 2750.0, 790.0, 293.15


# --------------------------------------------------------------- (a1)
# Each row: (label, file, regex, what it feeds). The regexes are matched
# against the CURRENT source, so a moved or deleted site shows up as a count
# change rather than as a stale quotation.
SITES = [
    ("1 pin reference", "src/Integrator/MMWSpalling/Removal.H",
     r"phases\[p_pin\]\.rho \* phases\[p_pin\]\.Cp",
     "rev_rhoCp -> pin_rhoCp, the idle clock's reference"),
    ("2 idle clock", "src/Integrator/MMWSpalling.H",
     r"pin_rhoCp\[c\] \* dz \* \(pin_T\[c\] - T_amb\) / q_pin",
     "t_cell, consuming site 1"),
    ("3 flake energy", "src/Integrator/MMWSpalling/Removal.H",
     r"phases\[0\]\.rho \* phases\[0\]\.Cp",
     "the flake energy estimate"),
    ("4 enthalpy table", "src/Integrator/MMWSpalling.H",
     r"phase_materials\[p\]\.m_rho\s*=\s*phases\[p\]\.rho",
     "phase_tables[p].Build -> H(T) and its inverse, and conduction"),
    ("5 gas side (OUT OF SCOPE)", "src/Integrator/MMWSpalling.H",
     r"St \* wall_rho \* wall_cp \* v",
     "the wall law -- GAS, not rock"),
    # ⚠⚠ THE TWO THE PACKET'S LIST OMITS. Neither is a `rho * Cp` product, so
    # no grep for one can find them: both read the two scalars SEPARATELY.
    ("6 vaporisation budget", "src/Integrator/MMWSpalling/Removal.H",
     r"Set::Scalar Cp_p\s*=\s*phases\[p_top\]\.Cp",
     "L_tot = Cp_p*(T_vap_lo - T0_v) + L_m + L_v, and rop_v = Q/(rho_p*L_tot)"),
    ("7 init fill of rho_mf/Cp_mf", "src/Integrator/MMWSpalling.H",
     r"rh\[p\]\s*=\s*phases\[p\]\.rho",
     "InitializeMicrostructure writes rho_mf/Cp_mf FROM THE SCALARS"),
]


def block_a1():
    print(BAR)
    print("(a1)  THE HARD-CODED rho*c_p ENUMERATION, RE-DERIVED")
    print(BAR)
    print("""
  ⚠⚠ THE PACKET'S LIST IS INCOMPLETE: it names five sites and bars (a1) at
  "the five sites above and no sixth". There are SEVEN. The two it misses read
  `rho` and `Cp` as SEPARATE scalars, so no grep for a `rho * Cp` product can
  find them -- which is exactly how they survived the planner's sweep.
""")
    print("  %-30s %-7s %s" % ("site", "hits", "what it feeds"))
    print("  " + SUB)
    ok = True
    for label, rel, pat, feeds in SITES:
        txt = open(os.path.join(ROOT, rel)).read()
        n = len(re.findall(pat, txt))
        if n != 1:
            ok = False
        print("  %-30s %-7s %s" % (label, n, feeds))
    print()
    print("""  THE TWO NEW ONES, and what they are worth:

  6  VAPORISATION. `rop_v = (1-R)*Q_local / (rho_p * L_tot)` with
     `L_tot = Cp_p*(T_vap_lo - T0_v) + L_m + L_v`. A CONSTANT c_p integrated
     over the whole span up to vaporisation -- the worst place for one, since
     Grimsel's c_p rises ~1.5x across it. ⚠ INERT IN THIS STEP: both decks set
     `vapor.enabled = 0`. Live in general, so it belongs in the list.
  7  INIT FILL. `InitializeMicrostructure` writes `rho_mf`/`Cp_mf` from the
     SCALARS, while `UpdatePropertyFields` writes the same two fields from the
     laws. ⚠ BENIGN BY ORDERING: `UpdatePropertyFields` (`MMWSpalling.H:1501`)
     runs before the conduction dispatch (`:1524`) every step, so conduction
     never reads the init values. But it is a SECOND WRITER of a field this
     step makes temperature-dependent, which is the shape P2-2a's (b9) locked
     against for the stress reference.
""")
    print("  every listed site found exactly once : %s" % ok)
    print("  ⇒ (a1) FAILS AS BARRED (five, no sixth); the enumeration is SEVEN.")
    return False, 7


# --------------------------------------------------------------- (a2)
def rhocp_constant(T):
    """`Material::Constant` + the base boxcars. ⚠ Both decks leave every
    phase-change scalar at its base default (`T_m = 1e30`, `fus_width = 0`,
    `L_m = 0`, `T_vap_hi == T_vap_lo == 1e30`), so the boxcars contribute
    EXACTLY ZERO on every leg this step runs -- verified from the decks."""
    return DECK_RHO * DECK_CP


def rhocp_grimsel(T):
    """Grimsel's laws: rho constant at 2650, c_p(theta) = 1370 - 178000/(theta+271)."""
    return G.RHO * G.c_p(T - 273.15)


def h_table(fn, T_ref, T_hi, n=4096):
    """Trapezoid per interval, from T_ref -- the same quadrature Table::Build
    uses, so this compares implementations and not quadratures."""
    Ts, Hs = [T_ref], [0.0]
    dT = (T_hi - T_ref) / n
    acc = 0.0
    for i in range(n):
        a, b = T_ref + i * dT, T_ref + (i + 1) * dT
        acc += 0.5 * (fn(a) + fn(b)) * dT
        Ts.append(b)
        Hs.append(acc)
    return Ts, Hs


def invert(Ts, Hs, H):
    """T(H) by linear interpolation on the same table -- the inverse the code
    does by bisection on the same curve."""
    if H <= Hs[0]:
        return Ts[0]
    if H >= Hs[-1]:
        return Ts[-1]
    lo, hi = 0, len(Hs) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if Hs[mid] <= H:
            lo = mid
        else:
            hi = mid
    f = (H - Hs[lo]) / (Hs[hi] - Hs[lo])
    return Ts[lo] + f * (Ts[hi] - Ts[lo])


def block_a2():
    print()
    print(BAR)
    print("(a2)  H(T) AND T(H) -- `Constant` FROM THE SCALARS vs GRIMSEL")
    print(BAR)
    print("""
  deck scalars : rho %.1f  c_p %.1f  T_ref %.2f   (input_drilling:102-106)
  Grimsel      : rho %.1f  c_p(theta) = 1370 - 178000/(theta+271)
  ⚠ The latent boxcars are ZERO on both forms here: neither deck sets T_m,
    fus_width, L_m, T_vap_lo/hi or L_v, so they sit at base defaults. ⇒ the
    packet's "latent heat silently vanishes" hazard is REAL IN CODE but has NO
    EFFECT on any leg this step runs, and (b3)'s ledger cannot see it here.
""" % (DECK_RHO, DECK_CP, DECK_TREF, G.RHO))
    Tc, Hc = h_table(rhocp_constant, DECK_TREF, 1200.0)
    Tg, Hg = h_table(rhocp_grimsel, DECK_TREF, 1200.0)
    print("  %-10s %-14s %-14s %-12s %-14s %s"
          % ("T [K]", "rhocp const", "rhocp grim", "ratio", "H const [J/m3]", "H grim / H const"))
    print("  " + SUB)
    worst_h, worst_hT = 0.0, 0.0
    for T in (293.15, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1023.15, 1200.0):
        rc_, rg_ = rhocp_constant(T), rhocp_grimsel(T)
        # H at T on each curve
        hcv = Hc[min(range(len(Tc)), key=lambda i: abs(Tc[i] - T))]
        hgv = Hg[min(range(len(Tg)), key=lambda i: abs(Tg[i] - T))]
        rat = (hgv / hcv) if hcv else float("nan")
        if hcv and abs(rat - 1.0) > abs(worst_h):
            worst_h, worst_hT = rat - 1.0, T
        print("  %-10.2f %-14.4e %-14.4e %-12.4f %-14.4e %.4f"
              % (T, rc_, rg_, rg_ / rc_, hcv, rat))
    print()
    print("  ⚠ THE SIGN MATTERS AND IT CHANGES. At room temperature Grimsel is")
    print("    COLDER-cheaper: rho 2650 < 2750 and c_p(20 C) = %.1f < %.1f, so"
          % (G.c_p(20.0), DECK_CP))
    print("    rho*c_p starts at x%.4f -- a DECREASE -- and crosses 1.0 near"
          % (rhocp_grimsel(293.15) / rhocp_constant(293.15)))
    xo = None
    for i in range(1, len(Tg)):
        if rhocp_grimsel(Tg[i]) >= rhocp_constant(Tg[i]):
            xo = Tg[i]
            break
    print("    T = %.1f K, reaching x%.4f at 1200 K."
          % (xo if xo else float('nan'),
             rhocp_grimsel(1200.0) / rhocp_constant(1200.0)))
    print("  ⇒ the ~1.3x the packet predicts is NOT a uniform scaling; it is the")
    print("    integral over the span the removed rock actually spans, which is")
    print("    why (d3) must report rho*c_p and mean-T-of-removed-rock apart.")
    print()
    print("  max |H_grim/H_const - 1| over 293-1200 K : %+.4f at T = %.1f K"
          % (worst_h, worst_hT))

    print()
    print("  T(H) -- the inverse, at equal enthalpy:")
    print("  %-16s %-14s %-14s %s" % ("H [J/m3]", "T const [K]", "T grim [K]", "diff [K]"))
    print("  " + SUB[:64])
    worst_t, worst_tH = 0.0, 0.0
    for frac in (0.1, 0.25, 0.5, 0.75, 1.0):
        H = frac * Hc[-1]
        tc_ = invert(Tc, Hc, H)
        tg_ = invert(Tg, Hg, H)
        if abs(tg_ - tc_) > abs(worst_t):
            worst_t, worst_tH = tg_ - tc_, H
        print("  %-16.4e %-14.2f %-14.2f %+.2f" % (H, tc_, tg_, tg_ - tc_))
    print()
    print("  max |T_grim - T_const| at equal H : %+.2f K at H = %.3e J/m3"
          % (worst_t, worst_tH))
    print("  ⚠⚠ THIS IS THE NUMBER (b3) IS ABOUT. If conduction integrates one")
    print("     curve and the H->T inversion reads the other, the books drift by")
    print("     this much per cell, and `ledger_err_solid` is the only row that")
    print("     can see it.")
    return True


# --------------------------------------------------------------- (a3)
def block_a3():
    print()
    print(BAR)
    print("(a3)  (d2)'s ARITHMETIC, WRITTEN OUT")
    print(BAR)
    print("""
  (d2) is ENERGY LEAVING AS ROCK, as a ratio to the leg's OWN baseline in the
  SAME descent mode. Per unit of charged energy:

      E_rock  =  price [J/m3]  x  volume removed [m3]
      (d2)    =  E_rock(property leg) / E_rock(its own baseline)

  ⚠⚠ IT IS A RATIO OF PRODUCTS, AND THAT IS ITS LIMIT. rev 5 called this
  "feedback-immune"; the packet corrects that, and the correction is the point:
  the product is immune to the five known absorption paths (firing rate, energy
  budget, wall surface T, wall coefficient, standoff) but NOT to a price error
  and a volume error that cancel. ⇒ (d3) reports the two factors separately,
  and (d2) alone cannot be read as a validation.

  THE FACTORS, and where each number must come from:

      price     = charged energy / charged volume, BOTH PATHS:
                  floor `h_applied` AND wall `lat_bank_hf x 6.0226`
                  ⚠ omitting the wall produced two withdrawn planner figures
      volume    = from the void MASK, five depth bands
      ⚠⚠ QUARTER DOMAIN: the deck is a quarter disk, so every volume is
         multiplied by 4. The ratio (d2) is dimensionless and the factor
         CANCELS -- but it must not cancel by being forgotten on one side,
         which is why it is written here.
""")
    print("  today's baseline, carried (method named where the number appears):")
    print("  %-34s %s" % ("realised price, CTRL5", "0.930 GJ/m3"))
    print("  %-34s %s" % ("predicted price, P2-2b", "~1.25 GJ/m3  (P1)"))
    print("  %-34s %s" % ("  from rho*c_p", "~1.3x   (d3), reported apart"))
    print("  %-34s %s" % ("  from mean T of removed rock", "~1.1-1.2x (d3)"))
    pred = 1.3 * 1.15
    print()
    print("  the two factors multiply to x%.3f, so the price goes" % pred)
    print("  0.930 -> %.3f GJ/m3, against P1's ~1.25." % (0.930 * pred))
    print()
    print("  ⚠⚠ P1'S OWN ARITHMETIC DOES NOT CLOSE, AND THE GAP IS NOT SMALL.")
    print("     P1 states the price goes 0.930 -> ~1.25 GJ/m3, which needs a")
    print("     combined factor of x%.3f. But its two stated factors give:"
          % (1.25 / 0.930))
    for mt in (1.10, 1.15, 1.20):
        print("       rho*c_p 1.3 x mean-T %.2f = x%.3f  ->  %.3f GJ/m3"
              % (mt, 1.3 * mt, 0.930 * 1.3 * mt))
    print("     ⇒ even the LOWEST corner of the stated range lands at 1.330,")
    print("       above the quoted 1.25. To reach 1.25 the mean-T factor would")
    print("       have to be x%.3f, not 1.1-1.2." % ((1.25 / 0.930) / 1.3))
    print("     ⚠ The ~1.3x itself IS supported: (a2) gives the ENTHALPY ratio")
    print("       at the Meier firing temperature (1023 K) as x1.298. So the")
    print("       factor is right and the ENDPOINT is wrong, or the mean-T")
    print("       factor is. The planner must reconcile them before (d3) is")
    print("       scored, or (d3) will 'confirm' a prediction that is")
    print("       inconsistent with itself.")
    print()
    print("  ⇒ (d2) = price_ratio x volume_ratio. With the price up x%.2f, (d2)"
          % pred)
    print("    lands in 0.95-1.10 only if the volume falls by x%.3f-%.3f."
          % (0.95 / pred, 1.10 / pred))
    print("    ⚠ THAT IS THE PREDICTION (d2) ACTUALLY MAKES, and it is a")
    print("      statement about VOLUME. Writing it this way first is what stops")
    print("      a volume change being read as a confirmation of the price.")
    print()
    print("  ⚠ The baseline's own E_rock cannot be computed offline: it needs the")
    print("    charged energy and the mask volume FROM THE LEG UNDER TEST. What is")
    print("    fixed here is the FORM, the quarter-domain factor, and the band's")
    print("    implication for volume -- which is what (a3) bars.")
    return True


def main():
    print(BAR)
    print("P2-2b ord 1 -- OFFLINE: (a1) enumeration, (a2) enthalpy, (a3) (d2)")
    print(BAR)
    print("  grimsel reference : %s" % os.path.relpath(G.ref_path(), ROOT))
    a1_ok, n_sites = block_a1()
    a2_ok = block_a2()
    a3_ok = block_a3()
    print()
    print(BAR)
    print("ORD 1 VERDICT")
    print(BAR)
    print("  (a1) enumeration complete at five  %s  -- it is %d sites"
          % ("** FAIL **", n_sites))
    print("  (a2) H(T)/T(H) divergence reported  %s" % ("PASS" if a2_ok else "FAIL"))
    print("  (a3) (d2) arithmetic written out    %s" % ("PASS" if a3_ok else "FAIL"))
    print()
    print("  ⚠⚠ (a1) is a GATE and it FAILS as written. The enumeration is not")
    print("     wrong in a way that threatens the step -- site 6 is inert here")
    print("     (`vapor.enabled = 0` on both decks) and site 7 is superseded")
    print("     before conduction reads it -- but the bar says 'no sixth' and")
    print("     there are two. The planner must extend the list to seven.")
    # ⚠ (a1) is a GATE, so its failure must set the exit code. Returning 0
    # on a failed gate is the 'fails open' shape that let a broken scorer
    # report 288/288 as 0 in P2-2a.
    return 0 if (a1_ok and a2_ok and a3_ok) else 1


if __name__ == "__main__":
    sys.exit(main())
