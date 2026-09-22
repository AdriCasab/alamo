#!/usr/bin/env python3
"""D2c deliverable 7: isentropic Laval state of Meier's burner nozzle and the
momentum-equivalent decay diameter.

Nozzle: throat Ø 6 mm (thesis), exit Ø 7.5 mm (drawing sheet 14), measured
mdot = 0.01513 kg/s (52 + 2.47 kg/h), products R = 290 J/kg K, gamma = 1.3
(1.25 and 1.35 reported), ambient 1.013 bar. Choked throat:
    mdot = p0 A_t sqrt(gamma/(R T0)) (2/(gamma+1))^((gamma+1)/(2(gamma-1)))
exit Mach from the supersonic root of A_e/A_t, then
    T_e = T0/(1 + (gamma-1)/2 M^2),  p_e = p0 (T_e/T0)^(gamma/(gamma-1)),
    u_e = M sqrt(gamma R T_e),       J = mdot u_e + (p_e - p_amb) A_e.
Momentum diameter in surroundings at T_surr (air, 287 J/kg K):
    D_e = 2 mdot / sqrt(pi rho_surr J),  rho_surr = p_amb/(287 T_surr).
At fixed mdot and J, D_e ~ sqrt(T_surr), which is what ALAMO applies:
D_e = jet_De_ref sqrt(T_mix/jet_T_ent). jet_De_ref is the FULL jet's D_e at
jet_T_ent (293.15 K); the quarter domain's jet_mdot/4 share must not be used.

`python nozzle.py` prints the tables (and checks the ACTIVE_STEP §Sources
table at 1900 K).
"""
import math

MDOT = 0.01513                 # kg/s, full jet (thesis)
D_THROAT = 6.0e-3              # m
D_EXIT = 7.5e-3                # m
R_GAS = 290.0                  # J/kg K, products
R_AIR = 287.0                  # J/kg K, surroundings
P_AMB = 1.013e5                # Pa
T_ENT = 293.15                 # K, jet_T_ent
GAMMA = 1.3
GAMMAS = (1.25, 1.3, 1.35)
T0S = (1600.0, 1750.0, 1900.0)
T_SURRS = (293.0, 800.0, 1500.0)
CORES = (5.0, 8.0)


def area(d):
    return math.pi * 0.25 * d * d


def area_ratio(M, g):
    """A/A* at Mach M."""
    return (1.0 / M) * ((2.0 / (g + 1.0)) * (1.0 + 0.5 * (g - 1.0) * M * M)) ** ((g + 1.0) / (2.0 * (g - 1.0)))


def mach_supersonic(ar, g):
    lo, hi = 1.0, 10.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if area_ratio(mid, g) < ar:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def laval(T0, g=GAMMA, mdot=MDOT, R=R_GAS, p_amb=P_AMB):
    A_t, A_e = area(D_THROAT), area(D_EXIT)
    flux = math.sqrt(g / (R * T0)) * (2.0 / (g + 1.0)) ** ((g + 1.0) / (2.0 * (g - 1.0)))
    p0 = mdot / (A_t * flux)
    M = mach_supersonic(A_e / A_t, g)
    T_e = T0 / (1.0 + 0.5 * (g - 1.0) * M * M)
    p_e = p0 * (T_e / T0) ** (g / (g - 1.0))
    u_e = M * math.sqrt(g * R * T_e)
    J = mdot * u_e + (p_e - p_amb) * A_e
    return dict(T0=T0, gamma=g, p0=p0, M_e=M, p_e=p_e, T_e=T_e, u_e=u_e, J=J,
                J_press=(p_e - p_amb) * A_e)


def rho_air(T_surr, p_amb=P_AMB):
    return p_amb / (R_AIR * T_surr)


def d_e(J, T_surr, mdot=MDOT):
    return 2.0 * mdot / math.sqrt(math.pi * rho_air(T_surr) * J)


def de_ref(T0, g=GAMMA):
    """jet_De_ref [m] for ALAMO: full-jet D_e at jet_T_ent."""
    return d_e(laval(T0, g)["J"], T_ENT)


def check_sources():
    """The ACTIVE_STEP §Sources table (planner and user reproduced), 1900 K."""
    s = laval(1900.0)
    rows = [("p0 [bar]", s["p0"] / 1e5, 5.95, 0.01),
            ("M_e", s["M_e"], 1.86, 0.005),
            ("p_e [bar]", s["p_e"] / 1e5, 0.97, 0.01),
            ("T_e [K]", s["T_e"], 1250.0, 10.0),
            ("u_e [m/s]", s["u_e"], 1277.0, 5.0),
            # the packet's 19.3 N is the momentum thrust mdot*u_e; J used here
            # adds the pressure thrust (p_e - p_amb)*A_e = -0.18 N (-0.9 %)
            ("mdot*u_e [N]", MDOT * s["u_e"], 19.3, 0.1)]
    de = {T: d_e(s["J"], T) * 1e3 for T in T_SURRS}
    # the packet's 800 K value is rounded (5.85 mm with J = mdot*u_e, 5.87 with
    # the pressure thrust); checked to 0.1 mm
    rows += [("D_e 293 K [mm]", de[293.0], 3.55, 0.06),
             ("D_e 800 K [mm]", de[800.0], 5.8, 0.1),
             ("D_e 1500 K [mm]", de[1500.0], 8.0, 0.06),
             ("8 D_e 293 K [mm]", 8.0 * de[293.0], 28.5, 0.6),
             ("8 D_e 1500 K [mm]", 8.0 * de[1500.0], 64.0, 0.6)]
    return [(n, v, ref, tol, abs(v - ref) <= tol) for n, v, ref, tol in rows]


def main():
    print("Laval state (throat Ø 6, exit Ø 7.5 mm, mdot 0.01513 kg/s, R 290, p_amb 1.013 bar)")
    print(f"{'T0 [K]':>7} {'gamma':>5} {'p0 [bar]':>8} {'M_e':>6} {'p_e [bar]':>9} {'T_e [K]':>7} "
          f"{'u_e [m/s]':>9} {'J [N]':>7} {'(p_e-p_a)A_e [N]':>16}")
    for T0 in T0S:
        for g in GAMMAS:
            s = laval(T0, g)
            print(f"{T0:7.0f} {g:5.2f} {s['p0'] / 1e5:8.3f} {s['M_e']:6.3f} {s['p_e'] / 1e5:9.3f} "
                  f"{s['T_e']:7.0f} {s['u_e']:9.1f} {s['J']:7.3f} {s['J_press']:16.4f}")
    print()
    print("Momentum diameter D_e = 2 mdot/sqrt(pi rho J) (gamma = 1.3) and core lengths")
    print(f"{'T0 [K]':>7} {'T_surr [K]':>10} {'D_e [mm]':>8} {'5 D_e [mm]':>10} {'8 D_e [mm]':>10}")
    for T0 in T0S:
        J = laval(T0)["J"]
        for T in (T_ENT,) + T_SURRS[1:]:
            de = d_e(J, T) * 1e3
            print(f"{T0:7.0f} {T:10.2f} {de:8.3f} {5 * de:10.2f} {8 * de:10.2f}")
    print()
    print("jet_De_ref (full jet, at jet_T_ent = 293.15 K) for ALAMO:")
    for T0 in T0S:
        print(f"  T_nozzle {T0:.0f} K: jet_De_ref = {de_ref(T0)!r} m "
              f"(gamma 1.25: {de_ref(T0, 1.25) * 1e3:.3f} mm, 1.35: {de_ref(T0, 1.35) * 1e3:.3f} mm)")
    print()
    print("Check against the ACTIVE_STEP §Sources table (1900 K, gamma 1.3):")
    ok = True
    for n, v, ref, tol, good in check_sources():
        ok = ok and good
        print(f"  {n:20s} {v:9.3f}  (packet {ref:g} ± {tol:g})  {'ok' if good else 'MISMATCH'}")
    s = laval(1900.0)
    print(f"  J used (incl. pressure thrust {s['J_press']:+.3f} N) = {s['J']:.3f} N, "
          f"{100 * (s['J'] / (MDOT * s['u_e']) - 1):+.2f} % vs mdot*u_e -> D_e "
          f"{100 * ((MDOT * s['u_e'] / s['J']) ** 0.5 - 1):+.2f} %")
    print("  all reproduced" if ok else "  MISMATCH: see above")
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
