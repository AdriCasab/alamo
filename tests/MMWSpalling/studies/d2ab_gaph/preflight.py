#!/usr/bin/env python3
"""D2ab item 0 -- the annulus coefficient, derived on paper BEFORE any leg runs.

Nothing in this file reads a D2ab output. Every input is either a Meier MEASURED
dimension, a key from the frozen SW/WL argv, or a tabulated gas property with a
named source. Run it, hash it, then launch.

  python3 preflight.py            the derivation + the pre-registered predictions
"""
import numpy as np

# ---------------------------------------------------------------- geometry
# ITEM 0.0 -- WHICH HOLE DIAMETER SETS THE GAP.
# Meier's MEASURED hole, 85-93 mm (Fig. 8.8, a visible width = a LOWER BOUND on
# diameter). NOT the leg's own diameter: that is a model OUTPUT, and using it
# would make h depend on the result it produces -- the same act as choosing
# f = 0.2 because it landed the rate.
D_HOLE = (0.085, 0.093)          # m, Meier measured
D_BURNER = 0.080                 # m, the Ø80 slotted skirt (= 2*foot_r_outer)
D_SHIELD = 0.067                 # m, the water-cooled shield above the skirt

# ------------------------------------------------------------- mass flow
# The QUARTER-hole key, x4. jet_P_cap = mdot*cp*(T_nozzle-T_ent) + P_recirc
# reproduces to 0.0 W on WL, which is what proves the quarter convention.
MDOT_NOZZLE = 4 * 0.0037826388888888888     # kg/s, full hole = 15.13 g/s
# The march's OWN marched mass flow, m_jet = jet_mdot/phi under
# jet_entrained_mass = 1, recovered from phi = (T_stag-T_rec)/(T_nozzle-T_rec)
# on the frozen legs over 150-400 s. SW is the argv base for GL/GH.
MDOT_MARCH = {"SW": 0.020874, "WL": 0.030597, "WL81": 0.032266, "A": 0.034479}

# ------------------------------------------------------- gas properties
# Air at 1 atm, Incropera & DeWitt, "Fundamentals of Heat and Mass Transfer",
# Table A.4. ⚠ The repo holds NO property table for combustion products; air is
# a stated approximation. The check that it is not a bad one: tabulated cp at
# 1600 K is 1248 J/kgK against the model's own jet_cp = 1250.
AIR_T  = np.array([1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000], float)
AIR_MU = np.array([473.0, 496.0, 530.0, 557.0, 584.0, 611.0, 637.0, 663.0, 689.0]) * 1e-7
AIR_K  = np.array([.0763, .0820, .0891, .0946, .1000, .1050, .1110, .1170, .1240])
AIR_PR = np.array([.728, .719, .703, .685, .688, .685, .683, .677, .672])
gas = lambda T: (np.interp(T, AIR_T, AIR_MU), np.interp(T, AIR_T, AIR_K),
                 np.interp(T, AIR_T, AIR_PR))

# Exhaust temperature, 150-400 s, from the frozen legs' own jet_T_exhaust.
T_EXH = {"SW": 1463.5, "WL": 1601.3, "WL81": 1631.0, "A": 1678.8}


def annulus(mdot, D_o, D_i, T):
    """Re and the two Nusselt ends for a concentric annulus.

    D_h = 4A/P = 2(R_o-R_i) = D_o - D_i, so
        Re = G D_h / mu = 2 mdot / (pi (R_o+R_i) mu)
    which is very nearly independent of which hole diameter is chosen.
    """
    mu, k, Pr = gas(T)
    D_h = D_o - D_i
    A = np.pi / 4 * (D_o**2 - D_i**2)
    Re = mdot * D_h / (A * mu)
    # LAMINAR, fully developed, one wall heated and the other insulated.
    # Incropera Table 8.2, annulus at D_i/D_o -> 1 (these run 0.86-0.94, i.e.
    # effectively parallel plates): Nu_o = 5.385.
    Nu_lam = 5.385
    # TURBULENT: Dittus-Boelter, Nu = 0.023 Re^0.8 Pr^0.4 (heating), with D_h
    # for the annulus. Incropera Eq. 8.60.
    Nu_turb = 0.023 * Re**0.8 * Pr**0.4
    return dict(D_h=D_h, A=A, Re=Re, mu=mu, k=k, Pr=Pr,
                h_lam=Nu_lam * k / D_h, h_turb=Nu_turb * k / D_h,
                Nu_lam=Nu_lam, Nu_turb=Nu_turb)


def main():
    L = []
    P = L.append
    P("=" * 78)
    P("D2ab ITEM 0 -- the gap-derived annulus coefficient")
    P("=" * 78)
    P("")
    P("0. WHICH DIAMETER SETS THE GAP")
    P("   Meier's MEASURED 85-93 mm, NOT the leg's own diameter. The leg's")
    P("   diameter is a model OUTPUT; using it would make h depend on the")
    P("   result it produces. Burner Ø80 mm = 2*foot_r_outer = the skirt that")
    P("   rests on the rock.")
    P(f"   gap = (D_hole - 80)/2 = {(D_HOLE[0]-D_BURNER)/2*1e3:.1f} - "
      f"{(D_HOLE[1]-D_BURNER)/2*1e3:.1f} mm")
    P("   ⚠ The packet says 'a gap of 2.5-6 mm'; (93-80)/2 = 6.5, so the top of")
    P("     the packet's range is 6.5 not 6. Reported, not silently replaced.")
    P("")
    P("1. HYDRAULIC DIAMETER, REYNOLDS NUMBER, REGIME")
    P("   D_h = 4A/P = 2(R_o - R_i) = D_o - D_i for a concentric annulus.")
    P("")
    P(f"   {'mdot basis':<22} {'D_hole':>7} {'D_h':>7} {'A':>9} {'T_gas':>7} "
      f"{'Re':>7}  regime")
    rows = []
    for basis, mdot in (("nozzle feed (15.13 g/s)", MDOT_NOZZLE),
                        ("march, SW (20.87 g/s)", MDOT_MARCH["SW"]),
                        ("march, A (34.48 g/s)", MDOT_MARCH["A"])):
        for D_o in D_HOLE:
            T = T_EXH["SW"] if "nozzle" in basis or "SW" in basis else T_EXH["A"]
            a = annulus(mdot, D_o, D_BURNER, T)
            # ⚠ The classical 2300 is the PIPE threshold. Annular transition
            # is reported from ~2000 upward, so 2000-2300 is called what it is:
            # the laminar edge of the transitional band, not settled laminar.
            reg = ("laminar" if a["Re"] < 2000 else
                   "transitional" if a["Re"] < 4000 else "turbulent")
            P(f"   {basis:<22} {D_o*1e3:6.0f}  {a['D_h']*1e3:6.1f}  "
              f"{a['A']*1e4:8.2f}  {T:6.0f}  {a['Re']:6.0f}  {reg}")
            rows.append((basis, D_o, a))
    P("")
    P("   ⚠ Re = 2 mdot / (pi (R_o+R_i) mu) -- the hole diameter cancels almost")
    P("     exactly, so the 85-93 mm uncertainty moves Re by < 5 %. What moves")
    P("     Re is the MASS FLOW BASIS, by a factor of 2.3.")
    P("")
    P("   WHICH MASS FLOW. At steady state only the nozzle feed crosses the")
    P("   exit plane: jet_T_ent_mode = exhaust means the march's entrained gas")
    P("   IS the hole's own recirculated exhaust, so counting it in the annulus")
    P("   double-counts mass that never leaves. PRIMARY = the nozzle feed.")
    P("   The marched value is carried as the upper sensitivity.")
    P("")
    P("2. NUSSELT CORRELATIONS AND THEIR SOURCES")
    P("   laminar   Nu = 5.385   Incropera Table 8.2, annulus D_i/D_o -> 1,")
    P("                          outer wall heated, inner insulated. These run")
    P("                          D_i/D_o = 0.86-0.94, i.e. parallel plates.")
    P("   turbulent Nu = 0.023 Re^0.8 Pr^0.4   Incropera Eq. 8.60 (Dittus-")
    P("                          Boelter, heating), on the annulus D_h.")
    P("   ⚠ Re = 2035-2134 on the primary basis is TRANSITIONAL. The choice")
    P("     between these two IS a factor of two and is NOT derivable. The")
    P("     bracket is their UNION, exactly as D2aa ran von Rohr at both ends.")
    P("   ⚠ Entry length: L_e ~ 0.05 Re D_h = 0.5 m laminar, far longer than")
    P("     the 50-200 mm annulus, so the real flow is in the ENTRANCE region")
    P("     where Nu exceeds the fully-developed value. This biases the true h")
    P("     toward the UPPER end of the bracket. Stated, not corrected for:")
    P("     an entry-length correction would be a third free choice.")
    P("")
    P("3. h_low AND h_high")
    P("")
    P(f"   {'basis':<22} {'D_hole':>7} {'h_lam':>8} {'h_turb':>8}")
    hs = []
    for basis, D_o, a in rows:
        P(f"   {basis:<22} {D_o*1e3:6.0f}  {a['h_lam']:7.1f}  {a['h_turb']:7.1f}")
        if "nozzle" in basis:
            hs += [a["h_lam"], a["h_turb"]]
    h_low, h_high = min(hs), max(hs)
    P("")
    P(f"   PRIMARY BRACKET (nozzle feed, Ø85-93, laminar U turbulent):")
    P(f"     h_low  = {h_low:6.1f} W/m^2K")
    P(f"     h_high = {h_high:6.1f} W/m^2K")
    P(f"   >>> GL uses side_face_h = {h_low:.0f}, GH uses side_face_h = {h_high:.0f} <<<")
    P("")
    allh = [a["h_lam"] for _, _, a in rows] + [a["h_turb"] for _, _, a in rows]
    P(f"   FULL union including the marched mass flow: {min(allh):.0f} - {max(allh):.0f}")
    P("   ⚠ NO point value inside the bracket is taken. Both ends are run.")
    P("")
    P("   SENSITIVITIES, reported and NOT folded into the bracket:")
    for D_o in D_HOLE:
        a = annulus(MDOT_NOZZLE, D_o, D_SHIELD, T_EXH["SW"])
        P(f"     Ø67 water-cooled shield instead of the Ø80 skirt, Ø{D_o*1e3:.0f}: "
          f"D_h {a['D_h']*1e3:.0f} mm, Re {a['Re']:.0f}, h {a['h_lam']:.0f}-{a['h_turb']:.0f}")
    mu_f, k_f, _ = gas(0.5 * (T_EXH["SW"] + 850.0))
    P(f"     film properties at (T_gas+T_wall)/2 = {0.5*(T_EXH['SW']+850):.0f} K: "
      f"k {k_f:.4f} vs {gas(T_EXH['SW'])[1]:.4f} -> h lower by "
      f"{(1-k_f/gas(T_EXH['SW'])[1])*100:.0f} %")
    P("")
    P("4. WHAT EACH END PREDICTS")
    P("=" * 78)
    # -- the power decomposition, measured on the frozen legs (D2aa analysis.txt)
    # P_rock = P_top + P_side, all FULL hole [kW].
    P_top = {"WL": 5.48 - 1.3216, "WL81": 4.94 - 0.8255, "SW": 8.05 - 4.0905}
    P("   TOP-FACE power is invariant to the side treatment:")
    for n, v in P_top.items():
        P(f"     {n:<5} P_top = {v:.2f} kW")
    Ptop = float(np.mean(list(P_top.values())))
    P(f"     mean {Ptop:.2f} kW  <- used below")
    P("")
    P("   BORE + PADS revert to the INHERITED map, so their power is SW's.")
    P("   ⚠ SW was frozen WITHOUT the region columns (the D2aa item-0(a) debt),")
    P("     so SW's own split is not on disk. Bracketed by WL's wall share:")
    P("       WL   wall = 468.6/1321.6 = 35.5 % of side power")
    P("       WL81 wall = 275.6/ 825.5 = 33.4 %")
    P("     so bore+pads = 60-70 % of SW's 4.09 kW side total.")
    bp_lo, bp_hi = 0.60 * 4.0905, 0.70 * 4.0905
    P(f"       bore+pads = {bp_lo:.2f} - {bp_hi:.2f} kW, central {0.5*(bp_lo+bp_hi):.2f}")
    P("")
    P("   WALL power scales with h, measured on the WL pair:")
    P(f"     h=142 -> 468.6 W  ({468.6/142:.2f} W per W/m^2K)")
    P(f"     h= 81 -> 275.6 W  ({275.6/81:.2f} W per W/m^2K)")
    slope = 0.5 * (468.6 / 142 + 275.6 / 81)
    P(f"     slope {slope:.2f} W per (W/m^2K), used linearly")
    P("")
    P("   RATE via the packet's rate ~ P^1.59, anchored on SW (8.05 kW,")
    P("   3.694 m/h). ⚠ On the frozen legs this anchor UNDER-predicts by")
    P("   3-4 % (D2aa: WL -3.3 %, WL81 -4.4 %), which is inside the 15 % bar.")
    P("")
    P(f"   {'leg':<5} {'h':>6} {'P_wall':>8} {'P_rock kW':>11} {'rate m/h':>10}")
    pred = {}
    for nm, h in (("GL", h_low), ("GH", h_high)):
        Pw = slope * h / 1000.0
        lo = Ptop + bp_lo + Pw
        hi = Ptop + bp_hi + Pw
        mid = 0.5 * (lo + hi)
        r = lambda p: 3.694 * (p / 8.05) ** 1.59
        P(f"   {nm:<5} {h:6.1f} {slope*h:7.0f} W {mid:8.2f}     {r(mid):9.2f}"
          f"   [{r(lo):.2f} - {r(hi):.2f}]")
        pred[nm] = (h, mid, r(mid))
    P("")
    P("   PRE-REGISTERED, criterion (a) scores RATE against these:")
    for nm, (h, p, r) in pred.items():
        P(f"     {nm}: h = {h:.0f} W/m^2K -> rate {r:.2f} m/h  "
          f"(15 % bar: {r*0.85:.2f} - {r*1.15:.2f})")
    P("")
    P("   ⚠ DISAGREEMENT WITH THE PACKET, reported not replaced: P2 pre-registers")
    P("     'GL/GH land 2.3-2.9 m/h'. This derivation gives "
      f"{pred['GL'][2]:.2f}-{pred['GH'][2]:.2f}, above that band. The")
    P("     difference is the bore+pads estimate: the packet's 'roughly 2.5 kW'")
    P("     side total implies P_rock 6.6 kW and 2.72 m/h with the same P^1.59,")
    P("     whereas SW's measured 4.09 kW side total at 60-70 % gives 2.45-2.86")
    P("     for bore+pads ALONE, before the wall is added. Criterion (a) is")
    P("     scored against ITEM 0 (this file), as the packet specifies; P2 is")
    P("     reported alongside.")
    P("")
    P("   SHAPE, REPORTED not scored -- by interpolation in rate on the frozen")
    P("   legs (A 1.251, WL81 1.780, WL 2.073, SW 3.694):")
    rr = np.array([1.251, 1.780, 2.073, 3.694])
    for key, vals in (("D@100 [mm]", [83.8, 85.2, 86.7, 87.6]),
                      ("D@125 [mm]", [81.3, 82.7, 83.6, 85.8]),
                      ("pit lead [mm]", [128.0, 98.0, 88.0, 38.0]),
                      ("volume [cm3]", [1.27, 1.22, 1.23, 1.08])):
        v = [np.interp(pred[n][2], rr, vals) for n in ("GL", "GH")]
        P(f"     {key:<14} GL {v[0]:6.1f}   GH {v[1]:6.1f}")
    P("   ⚠ Both predicted rates are ABOVE WL's 2.073 and below SW's 3.694, so")
    P("     these are interpolations, not extrapolations.")
    P("")
    P("=" * 78)
    P("D2l FACTOR-OF-FOUR AUDIT -- the debt owed by three packets. DISCHARGED.")
    P("=" * 78)
    P("   d2l_scan/analyze.py:")
    P("     vdot()      returns 4.0 * sum(h_applied) * dz^2 / dt   -> FULL hole")
    P("     tse_rock    = 4.0 * p_face / vmm3    -> 4 on BOTH sides, cancels")
    P("     tse_hhv     = P_HHV(38 kW, full burner) / vmm3(full hole)  consistent")
    P("     plot :397   4 * p_face, labelled 'absorbed (full hole)'")
    P("   RESULTS.md:190 'absorbed, full hole' f=1 -> 12.82 kW.")
    P("   Reproduced from d2l_scan/output/M_f10/thermo.dat over 150-250 s:")
    P("     jet_P_face (quarter) 3204.2 W  x4 = 12.82 kW   EXACT MATCH")
    P("   The f = 0 row (6.56 kW) is the linear-in-f extrapolation of")
    P("     f=0.1 7.35, f=0.2 7.92, f=0.5 10.10, f=1.0 12.82 kW (full hole).")
    P("   VERDICT: D2l DID apply the factor of four. '17.3 % of 38 kW HHV' and")
    P("   'needs 9.1 %' both stand as written. Nothing to restate. Debt closed.")
    P("")
    print("\n".join(L))


if __name__ == "__main__":
    main()
