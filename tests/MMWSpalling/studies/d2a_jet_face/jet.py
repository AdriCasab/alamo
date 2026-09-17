#!/usr/bin/env python3
"""D2a jet face source: Martin (1977) single round nozzle correlation, the local
coefficient derived from it, the hot-jet T_gas(s) decay, and the parser
expressions used by run.py. `python jet.py` prints the README arithmetic and
the h_loc verification.

Martin (1977), Adv. Heat Transfer 13:1-60, single round nozzle, area-averaged
over a disk of radius r (as quoted by Zuckerman & Lior 2006, Adv. Heat Transfer
39, and Incropera et al. section 7.7):
    Nu_avg / Pr^0.42 = G(r/D, H/D) * F(Re),   Nu_avg = h_avg D / k
    G = (D/r) (1 - 1.1 D/r) / (1 + 0.1 (H/D - 6) D/r)
    F = 2 Re^0.5 (1 + 0.005 Re^0.55)^0.5
    2000 <= Re <= 4e5, 2.5 <= r/D <= 7.5, 2 <= H/D <= 12.
Local coefficient: h_loc(r) = (1/2r) d(r^2 h_avg)/dr; with c = 0.1 (H/D - 6)
    h_loc/(k/D Pr^0.42 F) = (D/2r) (r^2 + 2cDr - 1.1cD^2) / (r + cD)^2.
"""
import math

import numpy as np

# Meier 2017 pilot (validation/meier/meier-2017-pilot-660kg-validation-reference.md)
D = 7.1e-3                                  # Laval nozzle outlet [m]
MDOT = (52.0 + 2.47) / 3600.0               # air + CH4 [kg/s] = 0.01513
T_AMB = 293.15
T_FIRE = 821.0                              # C1, V0 = V_cell (2 mm)

# Air, Sutherland's law (White, Viscous Fluid Flow, Table 1-2 / 1-3 constants):
#   mu = mu0 (T/T0)^1.5 (T0 + S_mu)/(T + S_mu),  mu0 = 1.716e-5 Pa s, S_mu = 110.4 K
#   k  = k0  (T/T0)^1.5 (T0 + S_k )/(T + S_k ),  k0 = 0.0241 W/mK,   S_k = 194 K
#   T0 = 273.15 K; cp = 1175 J/kgK (air near 1200 K); Pr = cp mu / k.
# The jet is combustion products (lean CH4/air, ~73 % N2 by mass); air
# properties are the stated approximation.
CP_AIR = 1175.0


def air(T):
    mu = 1.716e-5 * (T / 273.15) ** 1.5 * (273.15 + 110.4) / (T + 110.4)
    k = 0.0241 * (T / 273.15) ** 1.5 * (273.15 + 194.0) / (T + 194.0)
    return mu, k, CP_AIR * mu / k


def G_martin(r, H):
    return (D / r) * (1.0 - 1.1 * D / r) / (1.0 + 0.1 * (H / D - 6.0) * D / r)


def F_martin(Re):
    return 2.0 * Re ** 0.5 * (1.0 + 0.005 * Re ** 0.55) ** 0.5


def hloc_shape(r, H):
    """h_loc / (k/D Pr^0.42 F): (D/2r)(r^2 + 2cDr - 1.1cD^2)/(r + cD)^2."""
    c = 0.1 * (H / D - 6.0)
    return (D / (2.0 * r)) * (r * r + 2.0 * c * D * r - 1.1 * c * D * D) / (r + c * D) ** 2


def verify_hloc(H=7.0 * D, n=200001):
    """Area average of h_loc over the annulus [2.5D, r] plus the Martin disk
    average inside 2.5D must reproduce G(r): r^2 G(r) = (2.5D)^2 G(2.5D) +
    int_{2.5D}^{r} 2 r' h_loc(r') dr'. Returns max |rel error| over 2.5-7.5 D
    (trapezoid on n points) for H/D in (2, 7, 12)."""
    worst = 0.0
    for HD in (2.0, 7.0, 12.0):
        Hh = HD * D
        r = np.linspace(2.5 * D, 7.5 * D, n)
        f = 2.0 * r * hloc_shape(r, Hh)
        integ = np.concatenate(([0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(r))))
        G_rec = ((2.5 * D) ** 2 * G_martin(2.5 * D, Hh) + integ) / r ** 2
        worst = max(worst, float(np.max(np.abs(G_rec / G_martin(r, Hh) - 1.0))))
    return worst


def martin_reference(T_ref=1500.0, T_surf=T_FIRE):
    """h_avg and h_loc at r = 2.5 D, H = 7 D with film properties."""
    T_film = 0.5 * (T_ref + T_surf)
    mu, k, Pr = air(T_film)
    Re = 4.0 * MDOT / (math.pi * D * mu)
    F = F_martin(Re)
    Nu_avg = Pr ** 0.42 * G_martin(2.5 * D, 7.0 * D) * F
    h_avg = Nu_avg * k / D
    h_loc = Pr ** 0.42 * hloc_shape(2.5 * D, 7.0 * D) * F * k / D
    return dict(T_film=T_film, mu=mu, k=k, Pr=Pr, Re=Re, F=F, G=G_martin(2.5 * D, 7.0 * D),
                Nu_avg=Nu_avg, h_avg=h_avg, h_loc=h_loc)


# ---------------------------------------------------------------- expressions
# AMReX 25.12 parser gotcha (ext/amrex/Src/Base/Parser/AMReX_Parser_Y.cpp,
# parser_ast_optimize, PARSER_DIV case "f(.) / pow(x,n) => f(.) * pow(x,-n)"):
# the rewrite is applied to ANY two-argument function whose SECOND argument is a
# number, so a/max(r, 0.5) silently becomes a*max(r, -0.5). Every max/min below
# therefore takes its numeric argument FIRST (max(0.5, r)); run.py --parser-check
# verifies the compiled expressions inside ALAMO against the numpy forms.
def _hl(R, S):
    """Parser text for hloc_shape(R, S) with R, S parser sub-expressions."""
    c = f"(0.1*(({S})/{D!r}-6.0))"
    return (f"(({D!r}/(2.0*({R})))*(({R})*({R})+2.0*{c}*{D!r}*({R})-1.1*{c}*{D!r}*{D!r})"
            f"/((({R})+{c}*{D!r})*(({R})+{c}*{D!r})))")


def h_expr(h_ref, stagnation="flat"):
    """h(r, s) = h_ref * h_loc(max(r, 2.5D), clamp(s, 2D, 12D)) / h_loc(2.5D, 7D).
    stagnation = 'rise': inside 2.5 D the factor rises linearly to 1.5 at r = 0."""
    R = f"max({2.5 * D!r},r)"
    S = f"min({12.0 * D!r},max({2.0 * D!r},s))"
    norm = hloc_shape(2.5 * D, 7.0 * D)
    e = f"{h_ref!r}*{_hl(R, S)}/{norm!r}"
    if stagnation == "rise":
        e = f"({e})*(1.0+0.5*max(0.0,1.0-r/{2.5 * D!r}))"   # number-first already
    return e


def h_py(h_ref, r, s, stagnation="flat"):
    R = np.maximum(r, 2.5 * D)
    S = np.clip(s, 2.0 * D, 12.0 * D)
    h = h_ref * hloc_shape(R, S) / hloc_shape(2.5 * D, 7.0 * D)
    if stagnation == "rise":
        h = h * (1.0 + 0.5 * np.maximum(0.0, 1.0 - r / (2.5 * D)))
    return h


def Tg_expr(T_ref, T_ent=T_AMB):
    """T_gas(s) = T_ent + (T_ref - T_ent) min(7/5, 7D/max(s, 2D)); flat in the
    5 D core, 1/s beyond, = T_ref at s = 7 D. The max(s, 2D) only matters for
    s <= 0 columns outside the collision radius (then T_gas = core value)."""
    return (f"{T_ent!r}+({T_ref!r}-{T_ent!r})*min(1.4,{7.0 * D!r}/max({2.0 * D!r},s))")


def Tg_py(T_ref, s, T_ent=T_AMB):
    return T_ent + (T_ref - T_ent) * np.minimum(1.4, 7.0 * D / np.maximum(s, 2.0 * D))


if __name__ == "__main__":
    m = martin_reference()
    print(f"D = {D} m, mdot = {MDOT:.5f} kg/s, T_film = {m['T_film']:.1f} K")
    print(f"air (Sutherland): mu = {m['mu']:.4e} Pa s, k = {m['k']:.4f} W/mK, Pr = {m['Pr']:.3f}")
    print(f"Re = {m['Re']:.4e}, F = {m['F']:.2f}, G(2.5D, 7D) = {m['G']:.5f}")
    print(f"Nu_avg = {m['Nu_avg']:.1f}, h_avg(2.5D) = {m['h_avg']:.1f} W/m2K, h_loc(2.5D) = {m['h_loc']:.1f} W/m2K")
    print(f"h_loc verification (area average reproduces G over 2.5-7.5 D, H/D = 2, 7, 12): "
          f"max rel err {verify_hloc():.2e}")
