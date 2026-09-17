#!/usr/bin/env python3
"""D2a2 wall-jet helpers. Reuses D2a's jet.py (Martin 1977, h_loc, Sutherland
air, the parser-bug gotcha) at the burner-drawing nozzle D = 7.5 mm and the
drawing stand-off SOD = 50 mm (6.67 D). `python walljet.py` prints the README
arithmetic and the pre-registration numbers.

Closure (surface_patch.jet_closure = enthalpy): march outward from the
stagnation zone with T_gas(0) = T_stag = T_ent + (T_nozzle - T_ent) phi,
phi = min(1, 5 D / s_c) (free-jet centreline decay; phi = 1 for
jet_stagnation = nozzle), mass m = mdot (or mdot/phi with entrained mass),
    mdot cp dT_gas/dr = -q(r) 2 pi r,  q = h(r, s) (T_gas - T_s).
With h = h_ref r_core / r beyond r_core and uniform T_s the excess decays as
exp(-(r - r_core)/L), L = mdot cp / (2 pi h_ref r_core), and the flat
stagnation disk (h = h_ref) spends f_core = r_core/(2 L) of the excess
(to first order).
"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(os.path.dirname(HERE), "d2a_jet_face"))
import jet  # noqa: E402

jet.D = 7.5e-3                     # burner drawing (sheet 14/17), not the thesis 7.1 mm
D = jet.D
SOD = 0.050                        # nozzle exit to foot tip (sheet 1/17)
R_CORE = 2.5 * D                   # 18.75 mm
MDOT = jet.MDOT                    # 0.01513 kg/s, thesis measured (52 + 2.47 kg/h)
CP = 1250.0                        # T-averaged cp of the products (planner's value)
T_FIRE = jet.T_FIRE                # 821 K (V0 = V_cell)
TA = jet.T_AMB
RHOCP = 2750.0 * 790.0
SIG, EPS, HC = 5.670374e-8, 0.8, 10.0
DT_FIRE = T_FIRE - TA
T_NOZZLES = {"19": 1900.0, "14": 1436.0}   # adiabatic flame; chamber TC (uncorrected)
P_MEIER_ROCK = 2.47e-6 * RHOCP * 528.0     # W, rock-side removal power from the data


def martin(T_nozzle, H=SOD):
    """Martin at r = 2.5 D, stand-off H, film T = (T_nozzle + T_fire)/2."""
    T_film = 0.5 * (T_nozzle + T_FIRE)
    mu, k, Pr = jet.air(T_film)
    Re = 4.0 * MDOT / (math.pi * D * mu)
    F = jet.F_martin(Re)
    G = jet.G_martin(2.5 * D, H)
    return dict(T_film=T_film, mu=mu, k=k, Pr=Pr, Re=Re, F=F, G=G, Nu_avg=Pr ** 0.42 * G * F,
                h_avg=Pr ** 0.42 * G * F * k / D,
                h_loc=Pr ** 0.42 * jet.hloc_shape(2.5 * D, H) * F * k / D)


def h_anchor(key, T_nozzle):
    return {"M": round(martin(T_nozzle)["h_loc"]), "5": 5.0e3, "10": 1.0e4}[key]


NORM = jet.hloc_shape(2.5 * D, SOD)


def h_expr(h_ref):
    """h(r, s) = h_ref h_loc(max(r, 2.5D), clamp(s, 2D, 12D)) / h_loc(2.5D, SOD).
    Numeric arguments first in every min/max (AMReX 25.12 rewrite trap)."""
    R = f"max({2.5 * D!r},r)"
    S = f"min({12.0 * D!r},max({2.0 * D!r},s))"
    return f"{float(h_ref)!r}*{jet._hl(R, S)}/{NORM!r}"


def h_py(h_ref, r, s):
    return h_ref * jet.hloc_shape(np.maximum(r, 2.5 * D), np.clip(s, 2.0 * D, 12.0 * D)) / NORM


CORE_LEN = 5.0


def phi_of(s, stag="decay"):
    return 1.0 if stag == "nozzle" else min(1.0, CORE_LEN * D / s)


def t_stag(T_nozzle, s, stag="decay", T_ent=TA):
    return T_ent + (T_nozzle - T_ent) * phi_of(s, stag)


def loss(T):
    return EPS * SIG * (T ** 4 - TA ** 4) + HC * (T - TA)


def rop_closed(h, T_gas, T_f=T_FIRE):
    """Pinned closed-form ROP [m/h] for a face at T_fire."""
    return np.maximum(0.0, h * (T_gas - T_f) - loss(T_f)) / (RHOCP * (T_f - TA)) * 3600.0


def hand(h_ref, T_nozzle, cp=CP, r_out=0.170, mdot=MDOT):
    """The planner's closed forms (1/r beyond r_core, uniform T_s = T_fire)."""
    L = mdot * cp / (2.0 * math.pi * h_ref * R_CORE)
    f_core = R_CORE / (2.0 * L)
    ex0 = (T_nozzle - T_FIRE) * (1.0 - f_core)
    cap = mdot * cp * (T_nozzle - T_FIRE)
    P = cap * (1.0 - (1.0 - f_core) * math.exp(-(r_out - R_CORE) / L))
    # the planner's "hole" radius: excess falls to dT_fire (= 528 K), which is
    # what reproduces the planner's table (see README); None if already below
    r_hole = R_CORE + L * math.log(ex0 / DT_FIRE) if ex0 > DT_FIRE else None
    return dict(L=L, f_core=f_core, cap=cap, P170=P, r_hole=r_hole,
                cap_amb=mdot * cp * (T_nozzle - TA))


def march(h_ref, T_nozzle, cp=CP, s=SOD, r_max=0.170, n=17001, T_s=T_FIRE, mdot=MDOT,
          stag="nozzle", entrained=False, losses=False):
    """Continuous closure with the study's h(r, s) (flat core, Martin decay) and
    uniform T_s: returns r, T_gas(r), absorbed P(r), pinned closed-form ROP(r).
    stag = decay applies the centreline decay at s; entrained carries mdot/phi.
    losses = True takes the face losses at T_s out of the absorbed flux too
    (the reviewer's variant); default False (the planner's)."""
    phi = phi_of(s, stag)
    if entrained:
        mdot = mdot / phi
    r = np.linspace(0.0, r_max, n)
    h = h_py(h_ref, r, s)
    T = np.empty_like(r)
    T[0] = t_stag(T_nozzle, s, stag)
    P = np.zeros_like(r)
    dr = r[1] - r[0]
    for i in range(1, n):
        rm = r[i - 1] + 0.5 * dr
        hm = 0.5 * (h[i - 1] + h[i])
        # exact step for dT/dr = -a (T - T_s), a = 2 pi rm hm / (mdot cp)
        a = 2.0 * math.pi * rm * hm / (mdot * cp)
        ex = max(0.0, T[i - 1] - T_s)
        if losses:
            # net flux h (T - T_s) - loss(T_s), taken out of the gas while > 0
            ql = loss(T_s)
            q = max(0.0, hm * ex - ql)
            T[i] = max(T_s, T[i - 1] - q * 2.0 * math.pi * rm * dr / (mdot * cp))
        else:
            T[i] = T_s + ex * math.exp(-a * dr)
        P[i] = P[i - 1] + mdot * cp * (T[i - 1] - T[i])
    return r, T, P, rop_closed(h, T)


def ref_table(h_ref, T_nozzle, stag, entrained=False, losses=True, feed_mh=1.5):
    """The packet's reference: flat face at s = SOD; ring ROP at 40 mm; hole
    radius = largest r with closed-form ROP >= 0.9 feed."""
    r, T, P, rop = march(h_ref, T_nozzle, s=SOD, stag=stag, entrained=entrained, losses=losses)
    i40 = int(np.searchsorted(r, 0.040))
    ok = np.nonzero(rop >= 0.9 * feed_mh)[0]
    r_hole = float(r[ok[-1]]) if ok.size else 0.0
    return dict(rop40=float(rop[i40]), hole=2.0 * r_hole, P=float(P[-1]),
                T_stag=t_stag(T_nozzle, SOD, stag), T40=float(T[i40]))


def centre_rop(h_ref, T_nozzle, s, stag="decay"):
    return float(rop_closed(h_py(h_ref, 0.0, s), t_stag(T_nozzle, s, stag)))


def centre_rop_max(h_ref, T_nozzle, stag="nozzle", s_min=2.0 * D):
    """Pinned closed-form ROP on axis at the largest value over s in
    [s_min, 12D] (s beyond 12D only lowers it under decay)."""
    s = np.linspace(s_min, 12.0 * D, 201)
    return float(max(centre_rop(h_ref, T_nozzle, x, stag) for x in s))


def s_equilibrium(h_ref, T_nozzle, feed_mh, stag="decay"):
    """Stand-off where the centre's closed-form ROP equals the feed (decay);
    None if the centre outruns the feed at every s (nozzle treatment)."""
    s = np.linspace(SOD, 1.0, 20001)
    rop = np.array([centre_rop(h_ref, T_nozzle, x, stag) for x in s])
    below = np.nonzero(rop <= feed_mh)[0]
    return float(s[below[0]]) if below.size else None


def q_stag_max(h_ref, T_nozzle, stag="nozzle", s_min=2.0 * D):
    s = np.linspace(s_min, 12.0 * D, 201)
    return float(max(h_py(h_ref, 0.0, x) * (t_stag(T_nozzle, x, stag) - T_FIRE) - loss(T_FIRE)
                     for x in s))


if __name__ == "__main__":
    print(f"D = {D} m, SOD = {SOD} m = {SOD / D:.3f} D, r_core = {R_CORE * 1e3:.2f} mm, "
          f"mdot = {MDOT:.5f} kg/s, cp = {CP:g}")
    print(f"h_loc verification at D = 7.5 mm (G reproduced over 2.5-7.5 D, H/D = 2, 7, 12): "
          f"max rel err {jet.verify_hloc():.2e}")
    for Tn in (1900.0, 1436.0, 1500.0):
        m = martin(Tn)
        print(f"T_nozzle {Tn:g}: T_film {m['T_film']:.1f} K, mu {m['mu']:.4e}, k {m['k']:.4f}, "
              f"Pr {m['Pr']:.3f}, Re {m['Re']:.4e}, F {m['F']:.1f}, G(2.5D, SOD) {m['G']:.4f}, "
              f"Nu_avg {m['Nu_avg']:.1f}, h_avg {m['h_avg']:.0f}, h_loc {m['h_loc']:.0f} W/m2K")
    print(f"Meier rock-side removal power {P_MEIER_ROCK:.0f} W")
    print("\nplanner closed forms (cp 1250, T_s = T_fire):")
    for h in (1400.0, 3000.0, 5000.0, 10000.0):
        a, b = hand(h, 1900.0), hand(h, 1436.0)
        f = lambda x: "none" if x is None else f"{2e3 * x:.0f}"
        print(f"  h {h:6.0f}: L {a['L'] * 1e3:5.1f} mm f_core {a['f_core']:.2f} "
              f"P170@1900 {a['P170'] / 1e3:.1f} kW, 'hole' {f(a['r_hole'])} / {f(b['r_hole'])} mm; "
              f"cap {a['cap'] / 1e3:.1f} / {b['cap'] / 1e3:.1f} kW")
    print("\nreference table replica (flat face at SOD; ring ROP @40 mm / hole Ø [ROP >= 1.35 m/h]),"
          " with losses (reviewer) | without (planner):")
    for Tn in (1900.0, 1436.0):
        for stag, ent in (("decay", False), ("nozzle", False), ("decay", True)):
            row = []
            for hk in ("M", "5", "10"):
                h = h_anchor(hk, Tn)
                a, b = ref_table(h, Tn, stag, ent, True), ref_table(h, Tn, stag, ent, False)
                row.append(f"J{hk} {a['rop40']:.2f}/{a['hole'] * 1e3:.0f} | {b['rop40']:.2f}/{b['hole'] * 1e3:.0f}")
            print(f"  {Tn:g} {stag}{' +mass' if ent else ''}: T_stag {t_stag(Tn, SOD, stag):.0f} K; "
                  + "; ".join(row))
    print("\ncentre equilibrium stand-off at 1.5 m/h (decay):")
    for hk in ("M", "5", "10"):
        for tk, Tn in T_NOZZLES.items():
            h = h_anchor(hk, Tn)
            se = s_equilibrium(h, Tn, 1.5)
            print(f"  J{hk}_{tk}: s_eq {se * 1e3 if se else float('nan'):.0f} mm, "
                  f"centre ROP at SOD {centre_rop(h, Tn, SOD):.2f} m/h, "
                  f"q_stag max (s >= SOD - 5 mm) {q_stag_max(h, Tn, 'decay', SOD - 0.005) / 1e6:.2f} MW/m2")
    print("\nstudy anchors: centre ROP bound vs feed 1.5 m/h; continuous march at SOD:")
    for hk in ("M", "5", "10"):
        for tk, Tn in T_NOZZLES.items():
            h = h_anchor(hk, Tn)
            r, T, P, rop = march(h, Tn)
            i40 = np.searchsorted(r, 0.040)
            edge = r[np.nonzero(rop > 0)[0][-1]] if np.any(rop > 0) else 0.0
            print(f"  J{hk}_{tk} h {h:6.0f}: centre ROP <= {centre_rop_max(h, Tn):5.1f} m/h, "
                  f"T_gas(r_core) {T[np.searchsorted(r, R_CORE)]:.0f} K, T_gas(40 mm) {T[i40]:.0f} K, "
                  f"ROP(40 mm) {rop[i40]:.2f} m/h, q>loss edge r {edge * 1e3:.0f} mm, "
                  f"P(170 mm) {P[-1] / 1e3:.1f} kW")
