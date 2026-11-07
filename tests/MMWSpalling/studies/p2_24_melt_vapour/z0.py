#!/usr/bin/env python3
"""P2-24 z0 -- the offline numbers, each with its formula. Usage: z0.py > Z0.txt

The Grimsel law is transcribed from src/Numeric/Material/Grimsel.H (rho 2650; c_p = 1370 - 178000/(theta + 271);
lambda = 2000/(theta + 563); theta in C). The phase table is transcribed from src/Numeric/Material/Table.H (Build:
the trapezoid of rho c_p,eff on a 0.5 K grid from 200 K, the boxcars of Material.H effective_rhocp INCLUSIVE at both
band edges, shifted so H(T_ref) = 0). The melt set is the rock file's (Zhang, LIT).
"""
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
T = os.path.join(ROOT, "tests", "MMWSpalling")
DOE = os.path.join(T, "studies", "p2_21_loss_table", "doe")
W23 = os.path.join(T, "studies", "tree_sweep_0922", "witness", "face_p2_23")
DECK = os.path.join(T, "validation", "oglesby", "input")
RHO, SIG = 2650.0, 5.670374419e-8
T_M, FUS_W, L_M, T_VLO, T_VHI, L_V = 1510.5, 200.0, 3.4e5, 3233.0, 3503.0, 4.8e6
EPS_R, TAN_D, LAM = 5.02, 0.01065, 1.070687e-2
W0 = 0.01288
TAB_T = [295.15, 474.15, 573.15, 674.15, 724.15, 773.15, 823.15, 848.15, 873.15, 922.15, 972.15, 1021.15, 1070.15,
         1120.15, 1169.15, 1218.15, 1268.15]
TAB_R = [1.00, 0.98, 1.15, 1.39, 1.62, 2.01, 2.42, 2.73, 3.15, 4.27, 5.58, 5.72, 4.74, 3.79, 3.86, 6.02, 8.88]


def cp(T):
    return 1370.0 - 178000.0 / ((T - 273.15) + 271.0)


def lam(T):
    return 2000.0 / ((T - 273.15) + 563.0)


def cp_int(T1, T2):          # the analytic integral of c_p [J/kg]
    a, b = T1 - 273.15, T2 - 273.15
    return 1370.0 * (b - a) - 178000.0 * math.log((b + 271.0) / (a + 271.0))


def rhocp_eff(T):
    c = cp(T)
    if T_M - FUS_W / 2 <= T <= T_M + FUS_W / 2:
        c += L_M / FUS_W
    if T_VLO <= T <= T_VHI:
        c += L_V / (T_VHI - T_VLO)
    return RHO * c


def table(T_ref):
    """Table::Build: H at the 0.5 K nodes from 200 K to 5000 K, H(T_ref) = 0 (linear interpolation, as H_of_T)"""
    Tn = 200.0 + 0.5 * np.arange(int((5000.0 - 200.0) / 0.5) + 1)
    f = np.array([rhocp_eff(t) for t in Tn])
    H = np.concatenate([[0.0], np.cumsum(0.5 * (f[:-1] + f[1:]) * 0.5)])
    H -= np.interp(T_ref, Tn, H)
    return lambda t: float(np.interp(t, Tn, H))


def alpha_cold(lam_=LAM):
    return 2.0 * math.pi * math.sqrt(EPS_R) * (TAN_D / lam_)


def alpha_tab(T):
    return alpha_cold() * float(np.interp(T, TAB_T, TAB_R))


def f_L(T):
    return min(max((T - (T_M - FUS_W / 2)) / FUS_W, 0.0), 1.0)


def fig34_full():
    rows, sched = {}, None
    for ln in open(os.path.join(DOE, "FIG34_FULL.txt")):
        p = ln.split()
        if len(p) == 4 and re.match(r"^\d+$", p[0]):
            rows[int(p[0])] = float(p[1])
        if ln.startswith("the schedule (P_use"):
            sched = [float(x) for x in ln.split(":", 1)[1].split()]
    return rows, sched


def fig34_top():
    out = {}
    for ln in open(os.path.join(DOE, "FIG34_TOP.txt")):
        p = ln.split()
        if len(p) >= 2 and re.match(r"^\d+$", p[0]):
            try:
                out[int(p[0])] = float(p[1])
            except ValueError:
                pass
    return out


def deck_schedule():
    kv = {}
    for ln in open(DECK):
        if "=" in ln and not ln.lstrip().startswith("#"):
            k, v = ln.split("=", 1)
            kv[k.strip()] = v.split("#")[0].split()
    return [float(x) for x in kv["beam.schedule.t_end"]], [float(x) for x in kv["beam.schedule.P"]]


def hottest(run):
    import yt
    yt.set_log_level(50)
    best = (-1.0, None)
    for p in sorted(os.listdir(run)):
        if p.endswith("cell") and os.path.isdir(os.path.join(run, p)):
            ds = yt.load(os.path.join(run, p))
            tm = float(ds.all_data()["Temp"].to_value().max())
            if tm > best[0]:
                best = (tm, "%s (t %.4g s)" % (p, float(ds.current_time)))
    return best


def main():
    print("P2-24 z0 -- the offline numbers (units: s, W, K, J/m3, 1/m, m printed mm)")
    print("=" * 110)

    # 1 the power trace
    P, sched = fig34_full()
    lit = sorted(t for t, p in P.items() if math.isfinite(p) and 210 <= t <= 1790)
    t_end, Pd = deck_schedule()
    print("\n1. THE POWER TRACE (doe/FIG34_FULL.txt, fig34_digitise.py --full; Fig. 34's clock)")
    print("   lit windows (10 s, centres) %d .. %d s: %d; beam-on 210 s (P2-21's convention: 5 s before the first window "
          "centre, the window [210, 220]) => beam-off at the last window's end 1790 s; heating 1790 - 210 = 1580 s "
          "(DOE Table 1: 26 m 20 s = 1580 s)" % (lit[0], lit[-1], len(lit)))
    print("   ⚠ the packet wrote the shift as -215 s (beam-on at the first window's CENTRE) and beam-off 'about 1795 s'; the "
          "deck and this study use -210 s (the window's start), the digitiser's own convention since P2-21. Every DOE "
          "time below maps model = DOE - 210 s. The 5 s is half a window.")
    n_lit = sum(1 for p in Pd if p > 0)
    same = sched is not None and [round(x) for x in sched[:158]] == [round(x) for x in Pd[:158]]
    print("   the deck's schedule: %d intervals, %d with P > 0, the last lit ending %.0f s, then P = %.0f W to %.0f s; "
          "equal to the digitiser's printed schedule: %s" % (len(Pd), n_lit, t_end[n_lit - 1], Pd[-1], t_end[-1], same))
    print("   rows: " + " ".join("%d:%.0f" % (t_end[i], Pd[i]) for i in range(0, 158, 10)) + " ... 1580:%.0f" % Pd[157])
    for a, b in ((600, 800), (800, 1200), (1200, 1500), (1500, 1800)):
        w = [P[t] for t in lit if a <= t <= b]
        print("   DOE window %4d-%4d s (model %4d-%4d s): mean P %.3f kW over %d lit windows" % (
            a, b, a - 210, min(b, 1790) - 210, sum(w) / len(w), len(w)))
    TOP = fig34_top()
    print("   the radiometer (doe/FIG34_TOP.txt): the jump 272 -> 1195 -> 1649 C at 495 / 515 / 595 s (model onset "
          "convention: DOE onset t_on = 495 s, i.e. %d s after beam-on); beam-off from the POWER trace 1790 s, so "
          "the cooling readings 1805 / 1855 / 2005 s are at +15 / +65 / +215 s (⚠ the packet's +5 / +55 / +205 took "
          "beam-off at 1800 s); TOP rows present: %s" % (495 - 210, ", ".join("%d:%.0f" % (t, TOP[t]) for t in
                                                                              (495, 505, 515, 595, 1805, 1855, 2005)
                                                                              if t in TOP)))

    # 2 the hottest cell of each face-gate leg
    print("\n2. THE HOTTEST CELL OF EACH FACE-GATE LEG (face_p2_23, every written plotfile; a0/a1's identity needs every "
          "flame leg < 1410.5 K, where the table is unchanged)")
    legs = []
    for dp, ds, fs in os.walk(W23):
        if os.path.basename(dp) == "run" and any(x.endswith("cell") for x in ds):
            legs.append(dp)
    worst = 0.0
    for run in sorted(legs):
        tm, where = hottest(run)
        worst = max(worst, tm) if "BEAM_S" not in run else worst
        print("   %-12s max T %8.2f K at %s %s" % (os.path.relpath(os.path.dirname(run), W23), tm, where,
                                                  "" if tm < 1410.5 else "  ⚠ ABOVE the solidus"))
    print("   the hottest flame leg: %.2f K (< 1410.5 K: %s)" % (worst, worst < 1410.5))

    # 3 the enthalpy jump
    print("\n3. THE PHASE TABLE (Table::Build transcribed; H [J/m3], 0 at T_ref)")
    inv = RHO * L_M + RHO * cp_int(T_M - FUS_W / 2, T_M + FUS_W / 2)
    for T_ref, who in ((297.5, "M0 (the deck's T_ref)"), (300.0, "unit/face_vapor")):
        H = table(T_ref)
        print("   %s: H(1410.5) %.10e, H(1610.5) %.10e, H(3233) %.10e, H(3503) %.10e" % (
            who, H(1410.5), H(1610.5), H(3233.0), H(3503.0)))
    H = table(297.5)
    jump = H(1610.5) - H(1410.5)
    print("   the jump H(1610.5) - H(1410.5) = %.10e (table) vs rho L_m + rho int c_p dT = %.0f x %.0f + %.0f x %.4f = "
          "%.10e (rel %.1e; the boxcar is exact on the nodes, the c_p part is the trapezoid's)" % (
              jump, RHO, L_M, RHO, cp_int(1410.5, 1610.5), inv, jump / inv - 1.0))
    print("   the latent share rho L_m / jump = %.3f; the vapour jump H(3503) - H(3233) = %.4e vs rho L_v + rho int c_p "
          "= %.4e" % (RHO * L_M / jump, H(3503.0) - H(3233.0), RHO * L_V + RHO * cp_int(3233.0, 3503.0)))

    # 4 U1 / U1b
    m, a, L, Rm = 4.0, 1.0, 1.0, 0.3
    print("\n4. unit/face_melt U1 (alpha 1/m, L 1 m, liquid top half, m %.0f, R_melt %.1f): top (1-R_melt)(1-e^-(m a L/2)) "
          "%.9f, bottom (1-R_melt) e^-2 (1-e^-0.5) %.9f, transmitted (1-R_melt) e^-2.5 %.9f; U2 (flat) %.9f / %.9f / "
          "%.9f (⚠ the packet printed the top as 0.605298)" % (
              m, Rm, (1 - Rm) * (1 - math.exp(-m * a * L / 2)), (1 - Rm) * math.exp(-2) * (1 - math.exp(-0.5)),
              (1 - Rm) * math.exp(-2.5), 1 - math.exp(-0.5), math.exp(-0.5) * (1 - math.exp(-0.5)), math.exp(-1)))

    # 5 U3
    H3 = table(300.0)
    q = 3141.5926535897934 / (math.pi * 0.01 ** 2)
    t_v = 1e-3 * (H3(3503.0) - H3(300.0)) / (q * (1 - math.exp(-2)))
    print("\n5. unit/face_vapor U3: t_v = dz [H(T_vap_hi) - H(T0)] / (q (1 - e^-2)) = 1e-3 x %.6e / (%.4g x %.6f) = "
          "%.6f s => the crossing on step %d (dt 0.05), the void on step %d (the rule reads H at the step's start)" % (
              H3(3503.0) - H3(300.0), q, 1 - math.exp(-2), t_v, math.ceil(t_v / 0.05), math.ceil(t_v / 0.05) + 1))

    # 6 kappa and dt
    print("\n6. THE EXPLICIT LIMIT dx^2 / (6 kappa), kappa = lambda / (rho c_p,eff)")
    for Tk in (297.5, 1000.0, 1410.5, 1510.5, 1610.5, 2500.0, 3300.0):
        ce = rhocp_eff(Tk) / RHO
        k = lam(Tk) / (RHO * ce)
        print("   T %7.1f K: c_p,eff %7.1f J/kgK, kappa %.4e m2/s, limit %.4f s at 1.5625 mm, %.4f s at 0.78125 mm" % (
            Tk, ce, k, 0.0015625 ** 2 / (6 * k), 0.00078125 ** 2 / (6 * k)))
    print("   the room-temperature limit binds: dt 0.2 s (M0) / 0.05 s (M1) against %.4f / %.4f s" % (
        0.0015625 ** 2 / (6 * lam(297.5) / (RHO * cp(297.5))), 0.00078125 ** 2 / (6 * lam(297.5) / (RHO * cp(297.5)))))

    # 7 radiation-only surface T
    print("\n7. RADIATION ONLY (the surface sheds the absorbed flux by eps sigma (T^4 - T_amb^4), T_amb 297.5 K; "
          "absorbed = (1 - R_melt) = 0.7 of the incident)")
    for Pk in (2.0, 3.3):
        qc = 0.7 * 2 * Pk * 1e3 / (math.pi * W0 ** 2)
        r = 0.019
        qa = 0.7 * Pk * 1e3 * (1 - math.exp(-2 * r * r / W0 ** 2)) / (math.pi * r * r)
        for e in (0.8, 0.5):
            Tc = (qc / (e * SIG) + 297.5 ** 4) ** 0.25
            Ta = (qa / (e * SIG) + 297.5 ** 4) ** 0.25
            print("   P %.1f kW, eps %.1f: centre 0.7 x 2P/(pi w0^2) = %.3f MW/m2 => %.0f K (%.0f C); over the 38 mm "
                  "spot (area mean) %.3f MW/m2 => %.0f K (%.0f C)%s" % (
                      Pk, e, qc / 1e6, Tc, Tc - 273.15, qa / 1e6, Ta, Ta - 273.15,
                      "  ⚠ >= T_vap_lo 3233 K" if Tc >= T_VLO else ""))
    print("   eps 0.8 -> 0.5 raises a radiation-limited T by (0.8/0.5)^(1/4) = %.3f" % ((0.8 / 0.5) ** 0.25))

    # 8 transmission at full melt
    ac = alpha_cold()
    at = alpha_tab(1510.0)
    print("\n8. TRANSMISSION THROUGH A FULLY MOLTEN LAYER: alpha_cold = 2 pi sqrt(eps_r) tan_delta / lambda = %.3f 1/m; "
          "alpha_table(T >= 1268.15 K) = %.3f x 8.88 = %.2f 1/m (the code multiplies the TABLE value by the melt factor)"
          % (ac, ac, at))
    for mf in (4.0, 100.0):
        print("   x%-3.0f: m alpha_table = %.0f 1/m; e^-tau over 1 / 2 / 4 mm = %.3e / %.3e / %.3e; tau = 2 at %.2f mm"
              "  (on the COLD alpha: %.0f 1/m, tau = 2 at %.2f mm)" % (
                  mf, mf * at, math.exp(-mf * at * 1e-3), math.exp(-mf * at * 2e-3), math.exp(-mf * at * 4e-3),
                  2e3 / (mf * at), mf * ac, 2e3 / (mf * ac)))
    print("   ⚠ P1's 'x100: tau = 2 at 1.4 mm' is the COLD alpha (14 x 100); on the table (124 x 100) it is 0.16 mm, a "
          "tenth of one 1.5625 mm cell: at x100 the first molten cell takes the whole beam. P1's x4 number (497 1/m, "
          "4.0 mm) is the table's.")
    print("=" * 110)
    return 0


if __name__ == "__main__":
    sys.exit(main())
