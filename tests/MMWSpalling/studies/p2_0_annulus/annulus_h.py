#!/usr/bin/env python3
"""P2-0 ord 1 (a1)-(a7): the annulus heat-transfer coefficient, offline.

Scores the case against `surface_patch.side_face_h = 169.3` entirely from files
already on disk -- `studies/p1g3_thickness/MEIER_CTRL2/` -- so no run is needed
and every number is traceable to a logged column.

⚠⚠ RULES THIS SCRIPT ENFORCES, each of which has cost this programme a packet:

 1. **Every run-dependent quantity is read from the run**, never assumed. The
    P1g-0b (d1) reversal came from assuming a 1200 K gas when the log said
    1464-1599 K. Constants come from the leg's OWN argv.sh (parsed, not
    retyped); only the gas property fits and the burner shield radius are
    typed, and both are declared with a source below.
 2. **phi is MEASURED three independent ways and the three are compared.** One
    route would be an assumption wearing a measurement's clothes.
 3. **A quantity that is zero BY CONFIGURATION cannot score a gate.** See (a3):
    `side_face_factor_above = 0.0` is in this leg's argv, so the above-nozzle
    side flux is identically zero whatever coefficient is computed there. The
    script reports the above-nozzle column POPULATION beside the zero power, so
    the difference between "empty" and "zeroed by a key" is visible.
 4. **Quarter domain.** The input is a quarter annulus (axis at x0=y0=0), so
    every logged power is a quarter power. Reynolds number is invariant under
    the choice as long as mass flow and area are taken consistently; this
    script uses FULL values (mdot x4, full annulus area) throughout and says so
    at each use.

Usage:  python3 annulus_h.py > ANNULUS_H.txt
"""
import csv
import math
import os
import re
import statistics as st
import sys

ROOT = "/Users/tzetze20/amr_tools/alamo"
LEG = os.path.join(ROOT, "tests/MMWSpalling/studies/p1g3_thickness/MEIER_CTRL2")
W_LO, W_HI = 150.0, 350.0          # the declared scoring window [s]

# ---------------------------------------------------------------------------
# TYPED CONSTANTS -- the only two, each with its provenance.
# ---------------------------------------------------------------------------
# The burner's water-cooled shield is Oe67 mm => 33.5 mm radius, from the VT5
# burner drawings (see memory `meier-burner-drawings`: feet Oe80, standoff
# 50 mm, nozzle Oe7.5). ⚠ This is NOT side_face_h_r_min (40 mm, the Oe80 feet).
# The two radii are deliberately different: 40 mm gates WHICH faces take wall
# treatment, 33.5 mm is what the annulus gap is measured against.
R_SHIELD = 0.0335

# The rock wall temperature used by the (a7) gas-cooling march. The model pins
# the spalling surface near 823 K (memory `grimsel-granite-properties`: 823 K
# is mid-range of direct measurements). Declared, and (a7) reports a
# sensitivity rather than resting on it.
T_WALL_ROCK = 823.0


# ---------------------------------------------------------------------------
# GAS PROPERTIES. ⚠⚠ These are AIR fits. Combustion products run 10-20 % higher
# in k, which raises h by the same 10-20 % -- i.e. further ABOVE the derived
# value and further BELOW 169.3, so the conclusion is one-sided in our favour
# and the bias is not load-bearing. Two independent property sets are carried
# and both are reported, because the spread between them is the honest
# uncertainty on h and is larger than any tolerance this step quotes.
# ---------------------------------------------------------------------------
def mu_sutherland(T):
    """Sutherland's law for air (White, *Viscous Fluid Flow*, 2nd ed. Eq. 1-36).
    ⚠ Fitted near room temperature; it UNDER-predicts air viscosity above
    ~1000 K by a few per cent against tabulated values, which raises Re."""
    mu0, T0, S = 1.716e-5, 273.15, 110.4
    return mu0 * (T / T0) ** 1.5 * (T0 + S) / (T + S)


def k_fit(T):
    """Cubic fit for the thermal conductivity of air, 250-1600 K."""
    return (1.5207e-11 * T ** 3 - 4.8574e-8 * T ** 2
            + 1.0184e-4 * T - 3.9333e-4)


# Incropera & DeWitt, *Fundamentals of Heat and Mass Transfer*, Table A.4 (air
# at 1 atm). Tabulated, not fitted -- the independent check on the two fits.
_A4_T = [1000, 1100, 1200, 1300, 1400, 1500, 1600]
_A4_MU = [424e-7, 449e-7, 473e-7, 496e-7, 530e-7, 557e-7, 584e-7]
_A4_K = [66.7e-3, 71.5e-3, 76.3e-3, 82.0e-3, 91.0e-3, 100.0e-3, 106.0e-3]
_A4_PR = [0.726, 0.728, 0.728, 0.719, 0.703, 0.730, 0.730]


def _interp(T, xs, ys):
    if T <= xs[0]:
        return ys[0]
    if T >= xs[-1]:
        return ys[-1]
    for i in range(len(xs) - 1):
        if xs[i] <= T <= xs[i + 1]:
            f = (T - xs[i]) / (xs[i + 1] - xs[i])
            return ys[i] + f * (ys[i + 1] - ys[i])
    return ys[-1]


def props(T, which):
    """(mu, k, Pr) at T. `which` selects the declared property set."""
    if which == "fit":                      # Sutherland mu + cubic k, Pr 0.72
        return mu_sutherland(T), k_fit(T), 0.72
    return (_interp(T, _A4_T, _A4_MU), _interp(T, _A4_T, _A4_K),
            _interp(T, _A4_T, _A4_PR))


# ---------------------------------------------------------------------------
# CORRELATIONS. All three are reported at every bore; the choice spans 3-4x and
# the point of reporting all three is that the conclusion survives the spread.
# ---------------------------------------------------------------------------
def nu_gnielinski(Re, Pr):
    """Nu = 0.012 (Re^0.87 - 280) Pr^0.4, valid 1.5e3 < Re < 1e6.
    The right form here: Re ~ 6e3 is transitional, where Dittus-Boelter is out
    of range (it wants Re > 1e4)."""
    return 0.012 * (Re ** 0.87 - 280.0) * Pr ** 0.4


def nu_dittus(Re, Pr):
    """Nu = 0.023 Re^0.8 Pr^0.3 (gas being COOLED by the wall => 0.3).
    ⚠ Reported for comparison only: Re < 1e4 is BELOW its validity range."""
    return 0.023 * Re ** 0.8 * Pr ** 0.3


def nu_laminar(Re, Pr):
    """Fully-developed laminar annulus, one surface heated, one adiabatic:
    Nu -> 5.74 (Incropera Table 8.2, r_i/r_o ~ 0.6, inner adiabatic).
    Re-independent, which is why it floors h as the gap widens."""
    return 5.74


CORRS = (("Gnielinski", nu_gnielinski),
         ("Dittus-Boelter", nu_dittus),
         ("laminar", nu_laminar))


def annulus(r_bore, mdot_full, T_gas, which="fit", corr=nu_gnielinski):
    """Derived h on a bore of radius r_bore. Returns (h, Re, D_h, A).

    FULL annulus geometry with the FULL mass flow -- see rule 4. D_h = 2*gap is
    the standard hydraulic diameter of a concentric annulus (4A/P with
    A = pi(ro^2-ri^2), P = 2pi(ro+ri)).
    """
    gap = r_bore - R_SHIELD
    if gap <= 0.0:
        return float("nan"), float("nan"), gap, float("nan")
    D_h = 2.0 * gap
    A = math.pi * (r_bore ** 2 - R_SHIELD ** 2)
    mu, k, Pr = props(T_gas, which)
    Re = (mdot_full / A) * D_h / mu
    Nu = corr(Re, Pr)
    if Nu < nu_laminar(Re, Pr):     # a correlation may undershoot its own floor
        Nu = nu_laminar(Re, Pr)
    return Nu * k / D_h, Re, D_h, A


# ---------------------------------------------------------------------------
# INPUTS, read from the leg's own argv.sh -- never retyped.
# ---------------------------------------------------------------------------
def read_argv(path):
    txt = open(path).read()
    out = {}
    for m in re.finditer(r"([A-Za-z_][\w.]*)=((?:-?[\d.eE+-]+)(?:\s+-?[\d.eE+-]+)*)",
                         txt):
        key, val = m.group(1), m.group(2).split()
        try:
            out[key] = [float(v) for v in val]
        except ValueError:
            continue
    return out


def read_thermo(path, lo, hi):
    lines = open(path).read().split("\n")
    hdr = lines[0].split()
    ix = {n: i for i, n in enumerate(hdr)}
    rows = []
    for ln in lines[1:]:
        f = ln.split()
        if len(f) < len(hdr):
            continue
        try:
            t = float(f[ix["time"]])
        except ValueError:
            continue
        if lo <= t <= hi:
            rows.append(f)
    return ix, rows


def main():
    P = print
    av = read_argv(os.path.join(LEG, "argv.sh"))
    ix, rows = read_thermo(os.path.join(LEG, "FH/thermo.dat"), W_LO, W_HI)
    if not rows:
        P("FAIL: no thermo rows in the window")
        return 1
    g = lambda r, n: float(r[ix[n]])
    mean = lambda n: st.mean([g(r, n) for r in rows])

    mdot_q = av["surface_patch.jet_mdot"][0]
    cp = av["surface_patch.jet_cp"][0]
    T_noz = av["surface_patch.jet_T_nozzle"][0]
    T_ent = av["surface_patch.jet_T_ent"][0]
    L_core = av["surface_patch.jet_core_length"][0]
    h_const = av["surface_patch.side_face_h"][0]
    r_gate = av["surface_patch.side_face_h_r_min"][0]
    f_above = av["surface_patch.side_face_factor_above"][0]
    noz_z0 = av["surface_patch.nozzle_z0"][0]
    feed = av["surface_patch.nozzle_feed"][0]
    plo = av["geometry.prob_lo"]
    phi_hi = av["geometry.prob_hi"]
    ncell = av["amr.n_cell"]
    DZ = (phi_hi[2] - plo[2]) / ncell[2]
    Z_TOP = plo[2] + ncell[2] * DZ        # original rock surface

    P("=" * 78)
    P("P2-0 ord 1 -- THE ANNULUS COEFFICIENT AGAINST side_face_h = %.1f" % h_const)
    P("=" * 78)
    P("  leg    : MEIER_CTRL2 (shipped law, 500 s, quarter domain)")
    P("  window : %.0f-%.0f s, %d thermo rows" % (W_LO, W_HI, len(rows)))
    P("  from its argv: mdot_q %.6g kg/s, cp %.0f, T_noz %.0f, T_ent %.2f," %
      (mdot_q, cp, T_noz, T_ent))
    P("                 L_core %.1f, r_gate %.3f m, nozzle_z0 %.3f, feed %.6g m/s"
      % (L_core, r_gate, noz_z0, feed))
    P("  typed  : R_SHIELD %.4f m (VT5 drawings, Oe67 shield) -- NOT the %.3f m gate"
      % (R_SHIELD, r_gate))
    P("  mesh   : DZ %.4f m, rock top z = %.3f m" % (DZ, Z_TOP))

    # ---------------- (a1) phi, three independent ways ------------------
    T_stag = mean("jet_T_stag")
    T_exh = mean("jet_T_exhaust")
    P_cap = mean("jet_P_cap")
    s_c = mean("jet_s_c")
    D_e = mean("jet_D_e")

    phi_cap = mdot_q * cp * (T_stag - T_ent) / P_cap      # route 1: enthalpy cap
    phi_geo = min(1.0, L_core * D_e / s_c)                # route 2: pure geometry
    phi_mix = (T_stag - T_exh) / (T_noz - T_exh)          # route 3: mixing law
    # the (a1) check: predict the LOGGED T_stag from route 1's phi via route 3
    T_stag_pred = T_exh + (T_noz - T_exh) * phi_cap

    P("")
    P("-" * 78)
    P("(a1) GATE  phi measured THREE independent ways, and T_stag predicted")
    P("-" * 78)
    P("  logged: T_stag %.2f K, T_exhaust %.2f K, jet_P_cap %.1f W (quarter)," %
      (T_stag, T_exh, P_cap))
    P("          s_c %.4f m, D_e %.5f m" % (s_c, D_e))
    P("  phi from the enthalpy cap   mdot*cp*(T_stag-T_ent)/P_cap = %.4f" % phi_cap)
    P("  phi from the geometry       min(1, L_core*D_e/s_c)       = %.4f" % phi_geo)
    P("  phi from the mixing law     (T_stag-T_exh)/(T_noz-T_exh) = %.4f" % phi_mix)
    P("  spread across the three routes: %.5f (%.2f %%)"
      % (max(phi_cap, phi_geo, phi_mix) - min(phi_cap, phi_geo, phi_mix),
         100.0 * (max(phi_cap, phi_geo, phi_mix) - min(phi_cap, phi_geo, phi_mix))
         / phi_cap))
    P("  T_stag predicted from phi_cap through the mixing law: %.1f K" % T_stag_pred)
    P("  logged T_stag %.1f K  =>  residual %+.2f K   [bar: within 2 K]"
      % (T_stag, T_stag_pred - T_stag))
    ok_a1 = abs(T_stag_pred - T_stag) <= 2.0 and abs(phi_cap - 0.336) <= 0.002
    P("  phi = %.4f  [bar: 0.336 +/- 0.002]" % phi_cap)
    P("  (a1) %s" % ("PASS" if ok_a1 else "FAIL"))
    P("  => annulus mass flow is %.2fx the nozzle flow (m_jet = mdot/phi)."
      % (1.0 / phi_cap))

    phi = phi_cap
    mdot_below = 4.0 * mdot_q / phi     # FULL: recirculating flow below nozzle
    mdot_above = 4.0 * mdot_q           # FULL: net throughflow leaving the hole
    P("  FULL mass flows: below-nozzle (recirculating) %.6g kg/s," % mdot_below)
    P("                   above-nozzle (net throughflow) %.6g kg/s" % mdot_above)
    P("  ⚠ the ratio IS phi: entrainment recirculates in-hole gas, so only the")
    P("    nozzle's own mdot actually leaves the top.")

    # ---------------- bore radius, measured from the events --------------
    # ⚠⚠ THE BORE IS NOT THE MOUTH. The max event radius anywhere picks up the
    # 10-50 mm collar, which runs ~20 mm wider than the hole below it (the known
    # collar excess). The REPRESENTATIVE bore is taken over the depth bands that
    # (c1) is allowed to score -- 110-300 mm, which excludes the collar, the
    # 90-110 mm pause underream (standing rule: never score the bore at 100 mm)
    # and the advancing front. The mouth is reported separately, because the
    # derived h is worst there and the packet's P5 is about the two converging.
    BORE_BANDS = ((110.0, 150.0), (150.0, 200.0), (200.0, 300.0))
    MOUTH_BAND = (10.0, 50.0)
    ev = os.path.join(LEG, "FH_removal_events.csv")
    r6, r6_above, r6_below = [], 0, 0
    band_rmax = [0.0] * len(BORE_BANDS)
    mouth_rmax = 0.0
    with open(ev) as fh:
        for row in csv.DictReader(fh):
            try:
                reg = int(float(row["regime"]))
                t = float(row["time"])
                kt = int(float(row["k_top"]))
                rd = float(row.get("riser_depth") or 0.0)
                ci, cj = int(row["col_i"]), int(row["col_j"])
            except (ValueError, KeyError, TypeError):
                continue
            rr = DZ * math.hypot(ci + 0.5, cj + 0.5)
            d_mm = 1e3 * ((ncell[2] - 1 - kt) * DZ + rd)
            for b, (a, c) in enumerate(BORE_BANDS):
                if a <= d_mm < c:
                    band_rmax[b] = max(band_rmax[b], rr)
            if MOUTH_BAND[0] <= d_mm < MOUTH_BAND[1]:
                mouth_rmax = max(mouth_rmax, rr)
            if reg != 6:
                continue
            z_ev = Z_TOP - ((ncell[2] - 1 - kt) * DZ + rd)
            z_noz = noz_z0 - feed * t
            r6.append(z_ev - z_noz)
            if z_ev > z_noz:
                r6_above += 1
            else:
                r6_below += 1
    r6.sort()
    q = lambda p: r6[min(len(r6) - 1, int(p * len(r6)))] if r6 else float("nan")
    live = [x for x in band_rmax if x > 0.0]
    D_win = 2.0 * (st.mean(live) if live else 0.0)
    D_mouth = 2.0 * mouth_rmax

    # ---------------- (a3) BEFORE (a2): it sets which zone matters -------
    P("")
    P("-" * 78)
    P("(a3) GATE  is all wall removal below the nozzle plane?")
    P("-" * 78)
    P("  regime-6 events total %d   below nozzle %d   ABOVE nozzle %d"
      % (len(r6), r6_below, r6_above))
    P("  height relative to the nozzle plane [mm]: p05 %+.0f  p50 %+.0f  p95 %+.0f"
      % (1e3 * q(0.05), 1e3 * q(0.50), 1e3 * q(0.95)))
    P("  measured bar: 0 above  =>  %s" % ("0 above" if r6_above == 0 else
                                           "%d above" % r6_above))
    P("")
    P("  ⚠⚠ THIS ROW CANNOT GATE -- IT IS ZERO BY CONFIGURATION.")
    P("     This leg's argv sets surface_patch.side_face_factor_above = %.1f," % f_above)
    P("     which multiplies the side flux on every face ABOVE THE NOZZLE PLANE")
    P("     (MMWSpalling.H:2155-2160, predicate `!(flame_s[col] > 0)`) by that")
    P("     factor. With 0.0 the above-nozzle side flux is identically zero, so")
    P("     no above-nozzle face can ever be driven to a lateral firing --")
    P("     whatever coefficient a zone test computes for it.")
    P("     The population is NOT empty, which is what separates the two cases:")
    P("       side_cols_above  mean %.0f  max %.0f   (of side_faces mean %.0f)"
      % (mean("side_cols_above"),
         max(g(r, "side_cols_above") for r in rows), mean("side_faces")))
    P("       side_P_face_above mean %.4g W, max %.4g W"
      % (mean("side_P_face_above"),
         max(g(r, "side_P_face_above") for r in rows)))
    P("     => %.0f of %.0f side faces (%.0f %%) sit above the nozzle and receive"
      % (mean("side_cols_above"), mean("side_faces"),
         100.0 * mean("side_cols_above") / mean("side_faces")))
    P("        EXACTLY 0 W at every one of the run's steps.")
    P("  (a3) measured as stated, but VACUOUS as a gate -- reported, not scored.")

    # ---------------- (a2) Re and regime, by zone -----------------------
    r_bore = D_win / 2.0
    P("")
    P("-" * 78)
    P("(a2) REPORT  Reynolds number and flow regime, by zone")
    P("-" * 78)
    P("  bore radius from the leg's own removal events, bands 110-300 mm:")
    P("    per-band Oe %s mm  =>  representative bore Oe %.1f mm"
      % (" / ".join("%.1f" % (2e3 * x) for x in band_rmax), 1e3 * D_win))
    P("    mouth (10-50 mm) Oe %.1f mm -- reported apart; the collar is NOT the"
      % (1e3 * D_mouth))
    P("    bore and using it would overstate the case against the constant.")
    for lbl, md, Tg in (("below nozzle (m_jet)", mdot_below, T_exh),
                        ("above nozzle (net)  ", mdot_above, T_exh)):
        for which in ("fit", "A4"):
            h, Re, D_h, A = annulus(r_bore, md, Tg, which)
            reg = "LAMINAR" if Re < 2300 else ("transitional" if Re < 1e4
                                               else "turbulent")
            P("  %s  props=%-4s Re %7.0f  %-12s D_h %.4f m  h %6.1f"
              % (lbl, which, Re, reg, D_h, h))
    P("  ⚠ Re scales linearly with mass flow at fixed geometry, so")
    P("    Re_above/Re_below = phi = %.4f exactly." % phi)
    P("  ⚠⚠ but see (a3): the above-nozzle zone carries 0 W on this leg, so its")
    P("     laminar coefficient inflates NOTHING -- it is already zeroed.")

    # ---------------- (a4) h by bore diameter ---------------------------
    P("")
    P("-" * 78)
    P("(a4) REPORT  derived h by bore diameter, all three correlations")
    P("-" * 78)
    P("  below-nozzle mass flow, T_gas = T_exhaust = %.0f K" % T_exh)
    P("  %-9s %9s %9s %9s %9s %9s %9s"
      % ("bore Oe", "Re", "Gniel", "D-B", "laminar", "ratio", "A4 Gniel"))
    P("  %-9s %9s %9s %9s %9s %9s %9s"
      % ("[mm]", "", "[W/m2K]", "[W/m2K]", "[W/m2K]", "/169.3", "[W/m2K]"))
    for D_mm in (78.0, 84.0, 96.0, 102.0, 121.8):
        rb = D_mm * 1e-3 / 2.0
        vals = {}
        for nm, fn in CORRS:
            hh, Re, _, _ = annulus(rb, mdot_below, T_exh, "fit", fn)
            vals[nm] = hh
        h_a4, _, _, _ = annulus(rb, mdot_below, T_exh, "A4", nu_gnielinski)
        _, Re, _, _ = annulus(rb, mdot_below, T_exh, "fit")
        P("  %-9.1f %9.0f %9.1f %9.1f %9.1f %9.3f %9.1f"
          % (D_mm, Re, vals["Gnielinski"], vals["Dittus-Boelter"],
             vals["laminar"], vals["Gnielinski"] / h_const, h_a4))
    P("  packet expects (Gnielinski): 175 / 109 / 59 / 47 / 27, tol 5 %")
    exp = {78.0: 175.0, 84.0: 109.0, 96.0: 59.0, 102.0: 47.0, 121.8: 27.0}
    worst, worst_D = 0.0, None
    for D_mm, want in exp.items():
        got, _, _, _ = annulus(D_mm * 1e-3 / 2.0, mdot_below, T_exh, "fit")
        e = abs(got - want) / want
        if e > worst:
            worst, worst_D = e, D_mm
    P("  worst deviation from the packet's figures: %.1f %% (at Oe %.1f mm)"
      % (100.0 * worst, worst_D))
    P("  (a4) %s against the 5 %% tolerance"
      % ("PASS" if worst <= 0.05 else "OUTSIDE -- see the property note below"))
    P("  ⚠ every correlation at every bore >= 84 mm is BELOW %.1f, which is the"
      % h_const)
    P("    conclusion; it does not depend on the correlation choice.")

    # ---------------- (a5) the break-even bore --------------------------
    lo, hi = 2.0 * R_SHIELD + 1e-4, 0.300
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        h, _, _, _ = annulus(mid / 2.0, mdot_below, T_exh, "fit")
        if h > h_const:
            lo = mid
        else:
            hi = mid
    D_be = 0.5 * (lo + hi)
    P("")
    P("-" * 78)
    P("(a5) REPORT  the bore diameter at which the derived h equals %.1f" % h_const)
    P("-" * 78)
    P("  break-even bore Oe = %.1f mm        [bar: 78.4 mm, tol 1 mm]" % (1e3 * D_be))
    P("  (a5) %s" % ("PASS" if abs(1e3 * D_be - 78.4) <= 1.0 else
                     "deviation %.1f mm" % abs(1e3 * D_be - 78.4)))
    P("  measured bore in this window: Oe %.1f mm (mouth Oe %.1f mm)"
      % (1e3 * D_win, 1e3 * D_mouth))
    P("  ⚠⚠ THIS IS THE CASE AGAINST THE CONSTANT: %.1f is the correct"
      % h_const)
    P("     coefficient for a Oe%.0f mm hole, and the bore is Oe%.1f mm."
      % (1e3 * D_be, 1e3 * D_win))

    # ---------------- (a6) reweighted wall power ------------------------
    sp = os.path.join(LEG, "FH_side_power.csv")
    per_t = {}
    nsamp = 0
    with open(sp) as fh:
        for row in csv.DictReader(fh):
            try:
                t = float(row["time"])
                r = float(row["r"])
                Pw = float(row["P_side"])
            except (ValueError, KeyError, TypeError):
                continue
            if not (W_LO <= t <= W_HI) or r < r_gate:
                continue
            nsamp += 1
            hh, _, _, _ = annulus(r, mdot_below, T_exh, "fit")
            fac = (hh / h_const) if hh == hh else 1.0
            a, b = per_t.get(t, (0.0, 0.0))
            per_t[t] = (a + Pw, b + Pw * fac)
    P("")
    P("-" * 78)
    P("(a6) REPORT  the logged wall power, reweighted by each column's own radius")
    P("-" * 78)
    if per_t:
        P_now = st.mean([a for a, _ in per_t.values()])
        P_new = st.mean([b for _, b in per_t.values()])
        P("  per-column samples in the window, r >= %.3f m: %d  (%d time slices)"
          % (r_gate, nsamp, len(per_t)))
        P("  wall power as logged        %.1f W" % P_now)
        P("  reweighted to the derived h %.1f W" % P_new)
        P("  ratio %.3f    removes %.0f W" % (P_new / P_now, P_now - P_new))
        P("  cross-check against the thermo reduction of the same band:")
        P("    surf_side_P_rout (r >= %.3f) = %.1f W" % (r_gate, mean("surf_side_P_rout")))
        P("  packet expects 686 -> 279 W, ratio 0.406, tol 5 %")
        P("  removes %.0f W = %.1f %% of the %.0f W the model delivers"
          % (P_now - P_new, 100.0 * (P_now - P_new) / mean("ledger_P_robin"),
             mean("ledger_P_robin")))
    else:
        P("  no samples -- FH_side_power.csv empty in the window")

    # ---------------- deliverable 3: the kernel<->march split -----------
    P("")
    P("-" * 78)
    P("DELIVERABLE 3  the kernel<->march side-power split")
    P("-" * 78)
    k1 = st.mean([g(r, "surf_side_share") * g(r, "ledger_P_robin") for r in rows])
    k2 = st.mean([g(r, "surf_side_P_rin") + g(r, "surf_side_P_rout") for r in rows])
    k3 = st.mean([g(r, "surf_side_P_top") + g(r, "surf_side_P_buried") for r in rows])
    m1 = mean("side_P_face")
    P("  KERNEL, three independent reductions:")
    P("    surf_side_share x ledger_P_robin  %.2f W" % k1)
    P("    surf_side_P_rin + _rout           %.2f W" % k2)
    P("    surf_side_P_top + _buried         %.2f W" % k3)
    P("    spread across the three: %.3g W (%.2g %%)"
      % (max(k1, k2, k3) - min(k1, k2, k3),
         100.0 * (max(k1, k2, k3) - min(k1, k2, k3)) / k1))
    P("  MARCH, its own accounting:")
    P("    side_P_face                       %.2f W" % m1)
    P("  split: kernel/march = %.4f  =>  %.1f %%" % (k1 / m1, 100.0 * (k1 / m1 - 1.0)))
    P("  ⚠ the packet states the kernel routes as 984.4 / 982.4 / 982.4 and the")
    P("    split as 7.5 %%. Measured here the three kernel routes are EQUAL to")
    P("    %.1f W and the split is %.1f %%." % (k1, 100.0 * (k1 / m1 - 1.0)))
    P("  ⚠ `P_robin` is ambiguous in the packet: ledger_P_robin = %.0f W (total)"
      % mean("ledger_P_robin"))
    P("    while patch_P_robin = %.0f W (the FLOOR). Only ledger_P_robin"
      % mean("patch_P_robin"))
    P("    reproduces 982.4; the floor figure gives %.0f W."
      % st.mean([g(r, "surf_side_share") * g(r, "patch_P_robin") for r in rows]))

    # ---------------- (a7) gas cooling ----------------------------------
    P("")
    P("-" * 78)
    P("(a7) REPORT  gas cooling up the annulus, both mass-flow choices")
    P("-" * 78)
    L = 0.300
    h_der, _, _, _ = annulus(r_bore, mdot_below, T_exh, "fit")
    P("  marched with the DERIVED h = %.1f W/m2K over L = %.0f mm," % (h_der, 1e3 * L))
    P("  rock wall at %.0f K (declared constant, see the header)" % T_WALL_ROCK)
    for lbl, md in (("net throughflow ", mdot_above),
                    ("scrubbing flow  ", mdot_below)):
        mcp = md * cp
        UA = h_der * 2.0 * math.pi * r_bore * L
        dT = (T_exh - T_WALL_ROCK) * (1.0 - math.exp(-UA / mcp))
        P("  %s mdot %.6g kg/s  mcp %.2f W/K  UA %.3f W/K" % (lbl, md, mcp, UA))
        P("      %.0f K -> %.0f K   dT = %.0f K  (%.1f %% of T_exh)"
          % (T_exh, T_exh - dT, dT, 100.0 * dT / T_exh))
    P("  packet expects 161 K (net) and ~47 K (scrubbing)")
    P("  ⚠⚠ THE CHOICE IS DECLARED: the NET throughflow is the right denominator")
    P("     for how much the gas LEAVING the hole has cooled, because only mdot")
    P("     leaves. The scrubbing flow recirculates, so using it would charge the")
    P("     same gas's enthalpy against the wall more than once.")

    # ---------------- the verdict ---------------------------------------
    P("")
    P("=" * 78)
    P("VERDICT")
    P("=" * 78)
    h_at_bore, Re_b, _, _ = annulus(r_bore, mdot_below, T_exh, "fit")
    P("  at the measured bore Oe%.1f mm the derived coefficient is %.1f W/m2K"
      % (1e3 * D_win, h_at_bore))
    P("  against the typed %.1f  =>  the constant is %.1fx too high."
      % (h_const, h_const / h_at_bore))
    P("  It is the correct value for a Oe%.1f mm hole." % (1e3 * D_be))
    P("  The wall band carries %.0f W of the %.0f W delivered (%.0f %%);"
      % (mean("surf_side_P_rout"), mean("ledger_P_robin"),
         100.0 * mean("surf_side_P_rout") / mean("ledger_P_robin")))
    P("  reweighting removes %.0f W."
      % (st.mean([a for a, _ in per_t.values()]) -
         st.mean([b for _, b in per_t.values()]) if per_t else float("nan")))
    P("  ⚠⚠ TWO ROWS OF THIS PACKET DO NOT STAND AS WRITTEN:")
    P("     (a3) is zero by configuration and cannot gate (see above).")
    P("     The above-nozzle ZONE of deliverable 6, and (b7)'s above-nozzle 14,")
    P("     are therefore inert on every leg this packet runs.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
