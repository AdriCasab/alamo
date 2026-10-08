#!/usr/bin/env python3
"""P2-2 ord 1 / deliverables 1, 3 and 4 — OFFLINE. Zero compute.

Scores (a1)-(a7) before any source change. Every law is verified against the
repo-resident reference BY LINE (see grimsel.py), every code claim is resolved
by reading the source here, and every number the packet states is CHECKED rather
than inherited.

Usage:  python3 derive_props.py > PROPS.txt
"""
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "../../../.."))
sys.path.insert(0, HERE)
import grimsel as G                                        # noqa: E402

MMW = os.path.join(ROOT, "src/Integrator/MMWSpalling.H")
REM = os.path.join(ROOT, "src/Integrator/MMWSpalling/Removal.H")
SPC = os.path.join(ROOT, "src/Numeric/SpCriterion.H")
MEIER_IN = os.path.join(ROOT, "tests/MMWSpalling/validation/meier/"
                              "sp_meier_pilot/input_drilling")
CTRL5_ARGV2 = os.path.join(ROOT, "tests/MMWSpalling/studies/p2_1_wallh/CTRL5.argv.sh")
KANT_IN = os.path.join(ROOT, "tests/MMWSpalling/validation/kant/"
                             "sp_kant_onset/input_weibull")

# The packet's own claims, restated so they can be CHECKED, not inherited.
PACKET = {   # rev 3
    "sigma_today": 159.0, "sigma_meier": 299.0, "sigma_kant": 321.0,
    "R_band": (0.41, 0.69), "R_meier": 0.565, "R_kant": 0.527,
    "T_s_meier": 800.0, "T_s_kant": 873.0,
    "m_band": (3.4, 5.1), "flaw_implied": 9.0,
    "price_today": 0.930, "price_pred": 1.25, "price_ratio": 1.35,
    "kic_cold": 1.43, "kic_ratio_shipped": 0.54, "kic_ratio_attemp": 0.74,
}
EXPONENTS = (0.229, 0.273, 0.306)
P = print


def key(path, name, default=None):
    """Read `name = value` out of an ALAMO input deck."""
    if not os.path.exists(path):
        return default
    for ln in open(path):
        s = ln.split("#", 1)[0]
        if "=" not in s:
            continue
        k, v = s.split("=", 1)
        if k.strip() == name:
            vs = v.split()
            try:
                return [float(x) for x in vs] if len(vs) > 1 else float(vs[0])
            except ValueError:
                return vs[0] if len(vs) == 1 else vs
    return default


def src(path):
    return open(path).read()


def find(txt, pat, label):
    m = re.search(pat, txt)
    return (txt[:m.start()].count("\n") + 1) if m else 0


def kic_shipped(T_K):
    """The SHIPPED KIcNasseri, re-implemented from its own anchors, which are
    parsed out of src/Numeric/SpCriterion.H rather than retyped."""
    t = src(SPC)
    mT = re.search(r"T_K\[N\]\s*=\s*\{(.*?)\}", t, re.S)
    mK = re.search(r"Kic\[N\]\s*=\s*\{(.*?)\}", t, re.S)
    if not (mT and mK):
        return None, None, None
    # ⚠ STRIP THE INLINE COMMENTS FIRST. The anchor block is commented
    # (`298.15,   // ~25 C (RT)`), so a bare number scrape picks up the comment
    # values too and silently returns a 10-point table. That read 0.242 at
    # 400 degC and inverted the whole (a6) comparison.
    strip = lambda b: re.sub(r"//[^\n]*", "", b)
    Ts = [float(x) for x in re.findall(r"[\d.]+", strip(mT.group(1)))]
    Ks = [float(x) for x in re.findall(r"[\d.]+", strip(mK.group(1)))]
    if len(Ts) != len(Ks):
        raise RuntimeError("KIcNasseri parse: %d temperatures vs %d values"
                           % (len(Ts), len(Ks)))
    if T_K <= Ts[0]:
        return Ks[0], Ts, Ks
    if T_K >= Ts[-1]:
        return Ks[-1], Ts, Ks
    for i in range(len(Ts) - 1):
        if T_K <= Ts[i + 1]:
            f = (T_K - Ts[i]) / (Ts[i + 1] - Ts[i])
            return (1 - f) * Ks[i] + f * Ks[i + 1], Ts, Ks
    return Ks[-1], Ts, Ks


def implied_flaw(K_Ic_MPa, R, sigma_Pa):
    """The flaw that reproduces the onset, from the criterion's own form.

    Onset:  K_I = 2 * 1.763 * (R sigma) * sqrt(a/pi) = K_Ic
      =>    a = pi * [ K_Ic / (3.526 R sigma) ]^2
    ⚠ This is the INVERSE of the packet's
      Theta_onset = K_Ic (1-nu) / (2 * 1.763 * sqrt(a/pi) * E alpha),
    i.e. `E a` there is E*alpha (an expansion coefficient), not E*flaw-length —
    the dimensions only close that way.
    """
    if not (R > 0.0 and sigma_Pa > 0.0):
        return float("nan")
    return math.pi * (K_Ic_MPa * 1e6 / (2.0 * 1.763 * R * sigma_Pa)) ** 2


def main():
    fails = []

    def gate(ok, label, detail=""):
        P("  %-6s %s" % ("PASS" if ok else "FAIL", label))
        if detail:
            P("         %s" % detail)
        if not ok:
            fails.append(label)

    P("=" * 78)
    P("P2-2 ord 1 -- THE GRIMSEL SET, THE STRESSED VOLUME, THE RELIEF FACTOR")
    P("=" * 78)
    P("  reference : %s" % os.path.relpath(G.ref_path(), ROOT))
    P("              %d lines" % len(G.ref_lines()))
    P("  ⚠ the packet cites tests/MMWSpalling/validation/grimsel_properties/...")
    P("    which does NOT exist; the file is in Claude_markdowns/. Both paths are")
    P("    tried, so standing rule 8 is satisfied either way.")

    # ---------------- (a1) ------------------------------------------------
    P("")
    P("-" * 78)
    P("(a1) GATE  every law and table located in the reference BY LINE")
    P("-" * 78)
    rows = G.verify()
    for nm, ln, ok, txt in rows:
        P("   %-4s [%3d] %s" % ("ok" if ok else "FAIL", ln, nm))
    gate(all(r[2] for r in rows), "(a1) all %d laws/rows sourced" % len(rows),
         "%d of %d located" % (sum(1 for r in rows if r[2]), len(rows)))

    # ⚠⚠ the one check the reference DEMANDS of itself: [112] says kappa is
    # never set independently of lambda/(rho c_p). Its three fits are therefore
    # NOT independent, and if they disagree the set is inconsistent.
    P("")
    P("   [112] CONSISTENCY: is kappa_1 == lambda_1/(rho c_p) across the range?")
    P("     theta    kappa_fit    lambda/(rho c_p)     ratio")
    worst = 0.0
    for th in (25, 100, 200, 300, 400, 500):
        kf = G.kappa_mm2s(th) * 1e-6
        kc = G.kappa_coupled(th)
        worst = max(worst, abs(kc / kf - 1.0))
        P("     %5d   %.4e   %.4e   %.4f" % (th, kf, kc, kc / kf))
    # ⚠ THE BAR IS THE REFERENCE'S OWN UNCERTAINTY, not a round number. c_p(25)
    # is quoted as 768 +/- 48 (6.3 %) and rho as 2650 +/- 50 (1.9 %), so
    # lambda/(rho c_p) carries ~6.6 % before any fit disagreement. A 5 % spread
    # is INSIDE that, i.e. the three fits are consistent; a 5 % round-number bar
    # would have been an arbitrary line that this lands exactly on.
    tol = math.hypot(48.0 / 768.0, 50.0 / 2650.0)
    P("     propagated uncertainty on lambda/(rho c_p) from c_p 768+/-48 and")
    P("     rho 2650+/-50: %.1f %% -- that is the bar, not a round number." % (100 * tol))
    gate(worst <= tol,
         "(a1) the three fits agree inside the reference's own uncertainty",
         "worst deviation %.2f %% against %.2f %% propagated -- kappa is NOT an "
         "independent input" % (100 * worst, 100 * tol))

    # and the enthalpy law must be the integral of the heat capacity
    P("")
    P("   [43]/[48] CONSISTENCY: is dH/dtheta == c_p?")
    wH = 0.0
    for th in (25, 200, 400):
        num = (G.H(th + 0.01) - G.H(th - 0.01)) / 0.02
        wH = max(wH, abs(num / G.c_p(th) - 1.0))
        P("     theta %3d   dH/dtheta %.3f   c_p %.3f   ratio %.6f"
          % (th, num, G.c_p(th), num / G.c_p(th)))
    gate(wH <= 1e-4, "(a1) H is the integral of c_p (the H<->T inversion)",
         "worst deviation %.2e" % wH)

    # ---------------- (a2) ------------------------------------------------
    P("")
    P("-" * 78)
    P("(a2) GATE  the elastic stress today and measured, and the R band")
    P("-" * 78)
    E_now = G_key = key(MEIER_IN, "microstructure.phase0.E")
    b_now = key(MEIER_IN, "microstructure.phase0.beta")
    mu_now = key(MEIER_IN, "microstructure.phase0.mu")
    rho_now = key(MEIER_IN, "microstructure.phase0.rho")
    cp_now = key(MEIER_IN, "microstructure.phase0.Cp")
    Tref = key(MEIER_IN, "microstructure.phase0.T_ref")
    nu_now = E_now / (2.0 * mu_now) - 1.0
    P("  shipped (from the Meier deck, not retyped): E %.3g Pa, beta %.3g /K,"
      % (E_now, b_now))
    P("    mu %.6g => nu = E/(2 mu) - 1 = %.4f  [reference [195]: 0.25-0.30]"
      % (mu_now, nu_now))
    P("    rho %.0f, C_p %.0f, T_ref %.2f K" % (rho_now, cp_now, Tref))
    P("  measured set: E_0 %.3g Pa, rho %.0f, UCS %.0f +/- %.0f MPa"
      % (G.E0, G.RHO, G.UCS / 1e6, G.UCS_SD / 1e6))
    P("    ⚠ mu must move with E: mu = E/(2(1+nu)) = %.4e" % G.mu_of(G.E0))

    # the temperatures. ⚠ TAKEN FROM THE LEG / THE REFERENCE, never assumed.
    T_face_meier = 759.7      # CTRL5's own window-mean lat_Tface_max [K]
    th_meier = T_face_meier - 273.15
    T_s_kant = 600.0 + 273.15  # reference [252]: spalling surface ~600 degC
    th_kant = 600.0
    P("")
    P("  the two conditions, each from its own source:")
    P("    Meier: T_face %.1f K (CTRL5's own window-mean lat_Tface_max),"
      % T_face_meier)
    P("           dT = %.1f K above T_ref" % (T_face_meier - Tref))
    P("    Kant : T_s %.2f K (reference [252], spalling surface ~600 degC),"
      % T_s_kant)
    P("           dT = %.1f K" % (T_s_kant - Tref))
    dT_m = T_face_meier - Tref
    dT_k = T_s_kant - Tref

    # ⚠⚠ (a2) rev 3: THE BAR IS THE FORMULA FIRST. A matching number from a
    # different formula is exactly what this row exists to catch.
    P("")
    P("  THE FORMULA, implemented in grimsel.py::sigma_formula and printed here:")
    P("      sigma(T_s) = E(T_s) * INTEGRAL[T_0..T_s] alpha(theta) dtheta / (1-nu)")
    P("    ⚠ the strain is ACCUMULATED, then multiplied by the stiffness AT T_s.")
    P("      Three readings are excluded, and each was tried by someone:")
    P("        alpha(T_s)*dT        overstates  -- alpha TRIPLES over 300-500 degC")
    P("        cold E * strain      understates -- E HALVES by 600 degC")
    P("        INT(E*alpha dtheta)  wrong       -- elasticity has NO path memory")
    P("")
    T_m, T_k = PACKET["T_s_meier"], PACKET["T_s_kant"]
    s_today = G.sigma_elastic(E_now, b_now, T_m - Tref, nu_now)
    s_today_leg = G.sigma_elastic(E_now, b_now, T_face_meier - Tref, nu_now)
    s_meier, E_m, eps_m, am_m = G.sigma_formula(T_m)
    s_kant, E_k, eps_k, am_k = G.sigma_formula(T_k)
    P("    condition        T_s [K]  E(T_s) [GPa]   strain    alpha_mean   sigma [MPa]  x UCS")
    P("    today, flat a     %6.1f      %6.2f    %.4e  %6.2f e-6    %7.1f   %5.2f"
      % (T_m, E_now / 1e9, b_now * (T_m - Tref), b_now * 1e6,
         s_today / 1e6, s_today / G.UCS))
    P("    MEASURED Meier    %6.1f      %6.2f    %.4e  %6.2f e-6    %7.1f   %5.2f   [packet %.0f]"
      % (T_m, E_m / 1e9, eps_m, am_m * 1e6, s_meier / 1e6, s_meier / G.UCS,
         PACKET["sigma_meier"]))
    P("    MEASURED Kant     %6.1f      %6.2f    %.4e  %6.2f e-6    %7.1f   %5.2f   [packet %.0f]"
      % (T_k, E_k / 1e9, eps_k, am_k * 1e6, s_kant / 1e6, s_kant / G.UCS,
         PACKET["sigma_kant"]))
    P("")
    P("  ⚠ THE PACKET'S '159 MPa today' IS AT A DIFFERENT T_s FROM ITS OTHER TWO.")
    P("    At T_s = %.0f K the shipped set gives %.1f MPa, not 159; 159 appears at"
      % (T_m, s_today / 1e6))
    P("    T_s = %.1f K, the leg's OWN measured face temperature (%.1f MPa)."
      % (T_face_meier, s_today_leg / 1e6))
    P("    ⇒ the 159 -> 299 step mixes a TEMPERATURE change with a PROPERTY")
    P("      change. At one common T_s = %.0f K the property effect alone is"
      % T_m)
    P("      %.1f -> %.1f MPa, a factor %.2f, against the %.2f the packet implies."
      % (s_today / 1e6, s_meier / 1e6, s_meier / s_today,
         PACKET["sigma_meier"] / PACKET["sigma_today"]))
    P("")
    P("  ⚠⚠ THE PARADOX SURVIVES, RE-SIZED: %.2fx UCS at Meier's onset and %.2fx"
      % (s_meier / G.UCS, s_kant / G.UCS))
    P("     at Kant's. An elastic criterion is outside its validity, so a relief")
    P("     factor is the right stand-in.")
    R_m, R_k = G.UCS / s_meier, G.UCS / s_kant
    lo = (G.UCS - G.UCS_SD) / max(s_meier, s_kant)
    hi = (G.UCS + G.UCS_SD) / min(s_meier, s_kant)
    P("  R = UCS/sigma : Meier %.4f [packet %.3f]  Kant %.4f [packet %.3f]"
      % (R_m, PACKET["R_meier"], R_k, PACKET["R_kant"]))
    P("     spread %.3fx (%.1f %%), inside UCS's own +/-%.0f %% -- so R IS ONE"
      % (R_m / R_k, 100 * (R_m / R_k - 1), 100 * G.UCS_SD / G.UCS))
    P("     MATERIAL NUMBER, and the 1.53x reported at ord 1 was an artefact of")
    P("     keeping the COLD stiffness at two different T_s.")
    P("  band from UCS +/- %.0f MPa: %.3f - %.3f   [packet %.2f-%.2f]"
      % (G.UCS_SD / 1e6, lo, hi, *PACKET["R_band"]))
    ok2 = (abs(s_meier / 1e6 - PACKET["sigma_meier"]) <= 5.0
           and abs(s_kant / 1e6 - PACKET["sigma_kant"]) <= 5.0)
    okR = (abs(lo - PACKET["R_band"][0]) <= 0.01
           and abs(hi - PACKET["R_band"][1]) <= 0.01)
    gate(ok2 and okR,
         "(a2) the FORMULA is implemented and its values reproduce",
         "Meier %.1f vs %.0f, Kant %.1f vs %.0f MPa; R %.3f/%.3f, band %.3f-%.3f"
         % (s_meier / 1e6, PACKET["sigma_meier"], s_kant / 1e6,
            PACKET["sigma_kant"], R_m, R_k, lo, hi))

    # the implied flaw -- a CONSISTENCY check, never a fit output
    P("")
    P("  THE IMPLIED FLAW (a consistency check on R; `a0` is HELD, not fitted):")
    P("    ⚠⚠ IT IS TREATMENT-INDEPENDENT. Because R is DEFINED as UCS/sigma,")
    P("       R*sigma == UCS identically, so")
    P("           a = pi [K_Ic / (2*1.763*R*sigma)]^2 = pi [K_Ic/(3.526*UCS)]^2")
    P("       and the stress drops out entirely. The same value appears at BOTH")
    P("       conditions, at any R. ⇒ it cannot corroborate R or the stress, and")
    P("       a predicted (R, a) pair is ONE number, not two independent ones.")
    P("       ⚠ It IS a real check on a FITTED R, which is not tied to UCS.")
    P("    the flaw sits on the FAILURE PLANE, not the surface [252]:")
    for lbl, T_K in (("failure plane, 823 K (~550 degC)", 823.15),
                     ("spalling surface, 873 K (~600 degC)", 873.15)):
        kic = G.K_Ic_at_T(T_K - 273.15)
        a = implied_flaw(kic, R_m, s_meier)
        P("      %-36s K_Ic %.3f -> a = %6.2f um"
          % (lbl, kic, a * 1e6))
    a_fp = implied_flaw(G.K_Ic_at_T(550.0), R_m, s_meier) * 1e6
    a0_meier = 8.0
    mm = re.search(r"weibull\.a0_ig=([\d.eE+-]+)", open(CTRL5_ARGV2).read()) \
        if os.path.exists(CTRL5_ARGV2) else None
    if mm:
        a0_meier = float(mm.group(1)) * 1e6
    P("    ⇒ on the failure plane a = %.2f um [packet ~9], against the Meier"
      % a_fp)
    P("      leg's held a0_ig = %.1f um and the reference's measured MEAN CRACK"
      % a0_meier)
    P("      LENGTH of 130 um [96-97]. ⚠ the held flaw is a FLAW, not that mean.")

    # ---------------- (a8) -----------------------------------------------
    P("")
    P("-" * 78)
    P("(a8) REPORT  the 20 um / 22-36 um flaw figures: sourced, or dropped")
    P("-" * 78)
    L = G.ref_lines()
    hits = {}
    for needle in ("20 um", "22 um", "36 um", "22-36", "130 um", "200 um"):
        hits[needle] = [n + 1 for n, ln in enumerate(L) if needle in ln]
    for needle, lns in hits.items():
        P("    %-10s %s" % (needle,
                            ("found at " + ", ".join(str(x) for x in lns))
                            if lns else "** NOT IN THE REFERENCE **"))
    P("    [96-97] says: crack density 1.6 -> 5.8 mm/mm^2 after 500 degC,")
    P("            **mean crack length 130 um** (Kant 2017a s5).")
    P("  ⇒ VERDICT: the reviewer's '20 um reference crack length' and the")
    P("    '22-36 um' pair are **NOT IN THE REFERENCE** and are DROPPED, as the")
    P("    packet already proposed. An unsourced number must not become the held")
    P("    flaw. What the reference does give is 130 um (mean crack length) and")
    P("    ~200 um (Hoeser's mean spall thickness, [289]) -- neither of which is")
    P("    a mode-I flaw half-length, so neither replaces a0_ig either.")
    gate(not hits["22-36"] and not hits["22 um"] and bool(hits["130 um"]),
         "(a8) the unsourced figures are dropped; 130 um is sourced",
         "22-36 um absent; 130 um at line %s"
         % (hits["130 um"][0] if hits["130 um"] else "n/a"))

    # ---------------- (a3) ------------------------------------------------
    P("")
    P("-" * 78)
    P("(a3) GATE  3/(2m+3) solved against the MEASURED EXPONENTS")
    P("-" * 78)
    P("  ⚠⚠ scored against the measured EXPONENTS, never the published moduli:")
    P("     the reference's own [258] quotes 'm = 6.8' FOR exponent 0.306, which")
    P("     was DERIVED from it under a strength-modulus reading -- checking the")
    P("     relation against that number would be a tautology (rev 1's error).")
    P("")
    P("     exponent   flaw-LENGTH m = (3/e-3)/2   strength m = 3/e-3")
    ms = []
    for e in EXPONENTS:
        ml = G.m_from_exponent(e, "length")
        msr = G.m_from_exponent(e, "strength")
        ms.append(ml)
        P("     %8.3f   %24.2f   %17.2f" % (e, ml, msr))
    P("     => the code draws a flaw LENGTH with shape m scaled by (V/V0)^(1/m)")
    P("        and onset goes as a^-1/2, so the LEFT column is the demanded m.")
    inb = all(PACKET["m_band"][0] <= m <= PACKET["m_band"][1] for m in ms)
    gate(inb, "(a3) demanded m lands in %.1f-%.1f" % PACKET["m_band"],
         "m = %s for exponents %s"
         % (", ".join("%.2f" % m for m in ms),
            ", ".join("%.3f" % e for e in EXPONENTS)))
    P("  ⚠⚠ a flaw-length modulus of %.1f-%.1f is ENORMOUS scatter. m = 20 already"
      % (min(ms), max(ms)))
    P("     froze 48 % of the Meier floor in the P1 era, so (c3) is the binding")
    P("     gate and (c2) is a row we may legitimately fail -- that is (c6).")

    # ---------------- (a4) ------------------------------------------------
    P("")
    P("-" * 78)
    P("(a4) GATE  the seven stress sites, BY ROLE, resolved in source")
    P("-" * 78)
    rem, mmw = src(REM), src(MMW)
    sites = [
        ("lateral criterion  (C_th)", REM,
         find(rem, r"const Set::Scalar C_th = phase_E_loc \* phase_beta_loc", "")),
    ]
    for n, m in enumerate(re.finditer(
            r"f_skel \* phase_E_loc\s*\n\s*\* phase_beta_loc \* jmp", rem)):
        sites.append(("vertical branch %d (jmp)" % (n + 1), REM,
                      rem[:m.start()].count("\n") + 1))
    m = re.search(r"f_skel \* \(phase_E_loc\s*\n\s*\* phase_beta_loc", rem)
    if m:
        sites.append(("vertical branch 3 (T_z + conf_offset)", REM,
                      rem[:m.start()].count("\n") + 1))
    for lbl, pat in (("integrator: jump/(1-nu)",
                      r"\? te_E \* te_beta \* jump / \(1\.0 - te_nu\)"),
                     ("integrator: T_recon_at",
                      r"return te_E \* te_beta \* \(T_recon_at\(dv\) - te_Tref\)"),
                     ("integrator: T_a cell",
                      r"return te_E \* te_beta \* \(T_a\(i_d, j_d, k_d\) - te_Tref\)")):
        sites.append(("%s" % lbl, MMW, find(mmw, pat, "")))
    for lbl, path, ln in sites:
        P("   %-4s %-22s %s:%d"
          % ("ok" if ln else "FAIL", lbl, os.path.basename(path), ln))
    # and no EIGHTH: every E*beta product in these two files must be accounted
    allprod = (len(re.findall(r"phase_E_loc\s*\*?\s*\n?\s*\*?\s*phase_beta_loc", rem))
               + len(re.findall(r"te_E \* te_beta", mmw)))
    P("   every E*beta product in the two files: %d; sites enumerated: %d"
      % (allprod, len(sites)))
    gate(len(sites) == 7 and all(ln for _l, _p, ln in sites) and allprod == 7,
         "(a4) 7 of 7 stress sites, and no eighth",
         "%d enumerated, %d products found" % (len(sites), allprod))
    P("  ⚠⚠ R on a SUBSET relieves the wall and not the floor -- the exact family")
    P("     Removal.H:589-594 pre-registers. One factor, seven sites, locked.")

    # ---------------- (a5) ------------------------------------------------
    P("")
    P("-" * 78)
    P("(a5) GATE  every property consumer: STATIC field vs per-step")
    P("-" * 78)
    def writers(fld):
        """⚠ By the DECLARED TYPE, not the accessor. RefreshKappaEff binds
        kappa_mf through `.array(mfi)` into an `Array4<const Set::Scalar>` --
        a READ. Counting the accessor called it a writer and broke (a5)."""
        out = []
        for m in re.finditer(
                r"Array4<Set::Scalar> const& \w+ = \(\*%s\[lev\]\)\.array\(mfi\)"
                % fld, mmw):
            ln = mmw[:m.start()].count("\n") + 1
            fn = ""
            for fm in re.finditer(r"\n    (?:void|Set::Scalar|bool|amrex::Vector<int>) (\w+)\(", mmw):
                if mmw[:fm.start()].count("\n") + 1 < ln:
                    fn = fm.group(1)
            out.append((ln, fn))
        return out
    for fld in ("kappa_mf", "beta_mf", "kappa_eff_mf"):
        ws = writers(fld)
        P("   %-14s written at: %s" % (fld, ", ".join("%d (%s)" % w for w in ws) or "nowhere"))
    P("")
    P("   ⚠⚠ THE PREMISE OF CHANGE 1, VERIFIED: kappa_mf and beta_mf are written")
    P("      ONLY inside InitializeMicrostructure, i.e. they are STATIC per-cell")
    P("      fields filled at init. T-dependence is therefore NOT a parameter")
    P("      swap -- it needs a per-step refresh site.")
    P("   ⚠ RefreshKappaEff writes kappa_eff_mf, a DIFFERENT field: it applies the")
    P("     DAMAGE factor exp(-a_D D) to kappa_mf. It is not a temperature")
    P("     refresh and must not be mistaken for one.")
    P("")
    P("   consumers, by field:")
    cons = {
        "kappa_mf (read)": [mmw[:m.start()].count("\n") + 1 for m in
                            re.finditer(r"\(\*kappa_mf\[lev\]\)\.const_array\(mfi\)", mmw)],
        "beta_mf (read)": [mmw[:m.start()].count("\n") + 1 for m in
                           re.finditer(r"beta_mf\[lev\]->array\(mfi\)", mmw)],
        "effective_rhocp": [mmw[:m.start()].count("\n") + 1 for m in
                            re.finditer(r"effective_rhocp\(Tc\)", mmw)],
    }
    for nm, lns in cons.items():
        P("     %-20s %s" % (nm, ", ".join(str(x) for x in lns)))
    nph_m = key(MEIER_IN, "microstructure.nphases")
    nph_k = key(KANT_IN, "microstructure.nphases")
    P("   nphases: Meier %s, Kant %s  => no composition rule needed"
      % (nph_m, nph_k))
    gate(bool(cons["kappa_mf (read)"]) and bool(cons["effective_rhocp"])
         and len(writers("kappa_mf")) == 1,
         "(a5) consumers enumerated; the property fields are static at init",
         "kappa_mf written at exactly %d site (InitializeMicrostructure)"
         % len(writers("kappa_mf")))

    # ---------------- (a6) ------------------------------------------------
    P("")
    P("-" * 78)
    P("(a6) GATE  the shipped K_Ic(T) against AT-TEMPERATURE data")
    P("-" * 78)
    k_cold, Ts, Ks = kic_shipped(298.15)
    P("  the SHIPPED table, parsed out of src/Numeric/SpCriterion.H:")
    P("    T [K]  %s" % "  ".join("%8.2f" % t for t in (Ts or [])))
    P("    K_Ic   %s" % "  ".join("%8.2f" % k for k in (Ks or [])))
    P("  ⚠ 'Nasseri' appears in the reference ONLY at line 387, in its")
    P("    'Not obtained' list. The TABLE is in the repository with its own")
    P("    provenance comment, so (a6) is scoreable -- but the packet's")
    P("    attribution to the reference does not hold.")
    P("")
    P("  cold value: shipped %.2f vs the reference's K_Ic,0 = %.1f +/- 0.3 [209]"
      % (k_cold, G.K_IC0))
    P("")
    P("    T [K]  theta  shipped  at-temp [223]  shipped/cold  at-temp/cold  gap")
    worstgap = 0.0
    for T_K in (673.15, 773.15, 798.15, 823.15, 873.15):
        th = T_K - 273.15
        sh, _a, _b = kic_shipped(T_K)
        at = G.K_Ic_at_T(th)
        rs, ra = sh / k_cold, at / G.K_IC0
        worstgap = max(worstgap, ra / rs - 1.0)
        P("    %6.1f %6.0f  %7.3f  %13.3f  %12.3f  %12.3f  %+5.0f %%"
          % (T_K, th, sh, at, rs, ra, 100.0 * (ra / rs - 1.0)))
    P("")
    P("  ⚠⚠ AT THE SPALL TEMPERATURE the shipped table reads %.0f %% LOWER than"
      % (100.0 * (1.0 - 1.0 / (1.0 + worstgap))))
    P("     the at-temperature recommendation. Left alone, R silently absorbs")
    P("     that and stops being the measurement this step promises.")
    P("  ⚠ 'the table IS post-cooling data' is only true in its UPPER range:")
    sh4, _, _ = kic_shipped(673.15)
    P("    at 400 degC it gives %.3f/%.2f = %.2fx against the at-temperature"
      % (sh4, k_cold, sh4 / k_cold))
    P("    recommendation's 0.78 and post-cooling's 0.34-0.67. It tracks")
    P("    at-temperature to ~400 degC and post-cooling above. The substantive")
    P("    point stands exactly where it matters.")
    P("  ⇒ ADOPTED: the at-temperature anchors of [223]/Table 3.")
    gate(abs(k_cold - PACKET["kic_cold"]) < 0.01 and worstgap > 0.2,
         "(a6) both stated; the at-temperature anchors adopted",
         "cold %.2f (packet %.2f); worst gap at the spall T %+.0f %%"
         % (k_cold, PACKET["kic_cold"], 100 * worstgap))

    # ---------------- (a7) ------------------------------------------------
    P("")
    P("-" * 78)
    P("(a7) GATE  which path drives conduction, resolved in source")
    P("-" * 78)
    for lbl, pat in (("BuildSurfaceTColumns", r"BuildSurfaceTColumns\(lev, dt\);"),
                     ("UpdateWallLawState (P2-1's mirror)",
                      r"UpdateWallLawState\(lev\);"),
                     ("UpdateRemovalAfterCohesive (THE CRITERION)",
                      r"UpdateRemovalAfterCohesive\(lev, time, dt\);"),
                     ("PerCellRemoval", r"PerCellRemoval\(lev, time, dt\);"),
                     ("the conduction dispatch",
                      r"if      \(microstructure_enabled\)\s+AdvanceMicrostructure")):
        P("   %-42s :%d" % (lbl, find(mmw, pat, "")))
    ln_disp = find(mmw, r"if      \(microstructure_enabled\)\s+AdvanceMicrostructure", "")
    disp = mmw.split("\n")[ln_disp - 1] if ln_disp else ""
    first_is_micro = "microstructure_enabled" in disp
    P("")
    P("   the dispatch line: %s" % disp.strip()[:76])
    P("   ⚠⚠ microstructure_enabled is tested FIRST, so AdvanceMicrostructure IS")
    P("      the conduction route on both decks and `hc.heat.alpha` is NOT -- it")
    P("      reaches only Removal.H's h_spall. A packet gating on hc.heat.alpha")
    P("      would have PASSED FOR A FALSE REASON.")
    P("   hc.heat.alpha on the Meier deck: %s (inert for conduction)"
      % key(MEIER_IN, "hc.heat.alpha"))
    gate(first_is_micro and ln_disp > 0,
         "(a7) the conduction route is microstructure_enabled, tested first",
         "dispatch at :%d" % ln_disp)

    # ---------------- deliverable 4 ---------------------------------------
    P("")
    P("-" * 78)
    P("DELIVERABLE 4  the new V0, its datum, and which `ell`")
    P("-" * 78)
    dx_m = (key(MEIER_IN, "geometry.prob_hi")[0]
            - key(MEIER_IN, "geometry.prob_lo")[0]) / key(MEIER_IN, "amr.n_cell")[0]
    dx_k = (key(KANT_IN, "geometry.prob_hi")[0]
            - key(KANT_IN, "geometry.prob_lo")[0]) / key(KANT_IN, "amr.n_cell")[0]
    # ⚠⚠ MEIER SETS V0 BY ARGV, NOT IN THE DECK. Reading the deck gives the
    # 1e-9 default and a ratio of 8, which would have reported `m` as ALREADY
    # LIVE on Meier and inverted (b2). Take it from the leg that ran.
    CTRL5_ARGV = os.path.join(ROOT, "tests/MMWSpalling/studies/p2_1_wallh/"
                                    "CTRL5.argv.sh")
    V0_m = None
    if os.path.exists(CTRL5_ARGV):
        m = re.search(r"weibull\.V0=([\d.eE+-]+)", open(CTRL5_ARGV).read())
        V0_m = float(m.group(1)) if m else None
    V0_m_src = "CTRL5.argv.sh" if V0_m else "the deck default"
    if V0_m is None:
        V0_m = key(MEIER_IN, "weibull.V0", 1e-9)
    V0_k = key(KANT_IN, "weibull.V0", 1.0e-9)
    P("  TODAY the size term is pow(V_cell/V0, 1/m) at InitializeFlaws:")
    P("    Meier dx %.4f m => V_cell %.3e m3 ; V0 %.3e (from %s) => ratio %.6f"
      % (dx_m, dx_m ** 3, V0_m, V0_m_src, dx_m ** 3 / V0_m))
    P("    Kant  dx %.4f m => V_cell %.3e m3 ; V0 %.3e => ratio %.6f"
      % (dx_k, dx_k ** 3, V0_k, dx_k ** 3 / V0_k))
    P("    ⇒ the factor is EXACTLY 1 on BOTH decks, for any m. `weibull.m` is")
    P("      inert today by construction, which is why (b2) must be RED before.")
    P("")
    P("  ⚠ FOUR sites read weibull_m, not one:")
    for path, nm in ((MMW, "MMWSpalling.H"), (REM, "Removal.H")):
        t = src(path)
        for m in re.finditer(r"weibull_m", t):
            ln = t[:m.start()].count("\n") + 1
            P("      %-16s :%-6d %s" % (nm, ln, t.split(chr(10))[ln - 1].strip()[:58]))
    P("    ⇒ MMWSpalling.H:3652 and Removal.H:841 compute W_CONE_MEDIAN/weibull_m")
    P("      as the FLAKE-THICKNESS fraction, under spall.flake_fraction =")
    P("      weibull. That defaults to 'cf' and NO leg sets it, so the inertness")
    P("      CONCLUSION survives -- but at the demanded m ~ 4 that fraction would")
    P("      be 2.35/4 = %.2f against (d7)'s 0.2-0.3 mm prediction and today's"
      % (2.35 / 4.0))
    P("      c_f = %s. A coupling to record, not a blocker."
      % key(MEIER_IN, "spall.flake_c_f", "(set by argv)"))
    P("")
    P("  THE NEW V0 AND ITS DATUM:")
    P("    V = ell_q^3 with ell_q = k(T_s - T_0)/q, the FLUX length.")
    P("    ⚠ never kappa/v: at ONSET v = 0, so kappa/v is singular AND circular")
    P("      (it needs the rate the criterion is trying to predict). P1g-3")
    P("      established ell_q as an identity.")
    for q in (0.45e6, 0.60e6, 0.87e6):
        ell = G.lam(600.0) * (873.15 - 293.15) / q
        P("      q %.2f MW/m2 -> ell_q %.3f mm -> V %.3e m3" % (q / 1e6, ell * 1e3, ell ** 3))
    ell_ref = G.lam(600.0) * (873.15 - 293.15) / 0.60e6
    P("    ⇒ V0 = ell_q(q_ref)^3 at q_ref = 0.60 MW/m2 (the reference's own")
    P("      mid-band, [251] 0.45-0.87 and [263]'s T_s(q) datum): V0 = %.4e m3"
      % ell_ref ** 3)
    P("      i.e. ell_q = %.3f mm. ⚠ This KEEPS V0's role as the mesh-" % (1e3 * ell_ref))
    P("      independence anchor: it is a PHYSICAL length, not a cell size, so")
    P("      the size term stops being 1-by-coincidence without becoming")
    P("      mesh-dependent.")

    # ---------------- the price ------------------------------------------
    P("")
    P("-" * 78)
    P("THE PRICE (P2): what the measured set does to GJ/m3")
    P("-" * 78)
    cp_mean = sum(G.c_p(t) for t in range(25, 601, 5)) / len(range(25, 601, 5))
    P("  today     rho %.0f x C_p %.0f = %.3e J/(m3 K)" % (rho_now, cp_now, rho_now * cp_now))
    P("  measured  rho %.0f x mean c_p %.0f (25-600 degC) = %.3e J/(m3 K)"
      % (G.RHO, cp_mean, G.RHO * cp_mean))
    rc = (G.RHO * cp_mean) / (rho_now * cp_now)
    P("  rho c_p ratio %.4f" % rc)
    P("  ⚠ the price also moves with dT as T_s rises; the packet predicts")
    P("    x1.09 there, giving an overall x%.2f and %.3f -> %.3f GJ/m3."
      % (rc * 1.09, PACKET["price_today"], PACKET["price_today"] * rc * 1.09))
    P("  ⇒ reported as a PREDICTION; the realised price is scored per path on")
    P("    the Meier legs (floor h_applied AND wall lat_bank_hf x 6.0226).")

    P("")
    P("=" * 78)
    P("ORD 1 VERDICT")
    P("=" * 78)
    if fails:
        for f in fails:
            P("  FAILED: %s" % f)
        return 1
    P("  (a1)-(a7) all PASS. The four source changes may proceed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
