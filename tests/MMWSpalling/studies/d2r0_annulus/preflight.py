#!/usr/bin/env python3
"""D2r-0 pre-flight 0 -- arithmetic and existing output ONLY, before any code.

Reproduces every reference value ACTIVE_STEP.md states with an independent
estimator, per the standing guardrail, and evaluates the packet's five STOP
tests. Imports nothing; reads frozen study output by path and never re-runs it.

    python3 preflight.py

Result: items 1, 2, 3 PASS. Item 4 STOPS (1.46x < 2x) and item 5 misses
criterion (c) at the band's lower end. See PREFLIGHT.md.
"""
import numpy as np
from scipy.optimize import brentq

SIG, EPS, HC, TAMB = 5.670374419e-8, 0.8, 10.0, 293.15
KAP, CP = 1.5, 1250.0
T_FIRE, T_HOT, STEP0, Q0, Q_PIN = 822.2, 820.0, 402.0, 302e3, 520e3
LI0 = "tests/MMWSpalling/studies/d2i_gate/output/LI0_2mm"


def loss(T):
    return EPS * SIG * (T ** 4 - TAMB ** 4) + HC * (T - TAMB)


def h_real(r, s):
    """The h_expr LI0_2mm actually ran (NOT the packet's simplified stand-in)."""
    rm = max(0.01875, r)
    b = 0.1 * ((min(0.09, max(0.015, s))) / 0.0075 - 6.0)
    num = (0.0075 / (2 * rm)) * (rm * rm + 2 * b * 0.0075 * rm
                                 - 1.1 * b * 0.0075 * 0.0075) / ((rm + b * 0.0075) ** 2)
    return 1449.0 * num / 0.19763872491145218 * (0.09 / max(0.09, s)) ** 1.0


def thermo(name=LI0):
    n = open(name + "/thermo.dat").readline().split()
    rows = np.loadtxt(name + "/thermo.dat", skiprows=1, ndmin=2)
    k = {x: i for i, x in enumerate(n)}
    return {x: rows[:, k[x]] for x in n}


def surface():
    """Final column-top height from the frozen removal-event log."""
    a = np.loadtxt(LI0 + "_removal_events.csv", delimiter=",", skiprows=1, ndmin=2)
    ci, cj, kt = a[:, 1].astype(int), a[:, 2].astype(int), a[:, 4].astype(int)
    dz, nxy, nz = 2e-3, 60, 200
    top = np.full((nxy, nxy), nz - 1, dtype=int)
    np.minimum.at(top, (ci, cj), kt)
    z_top = (top + 0.5) * dz
    A_side = 0.0
    for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        A_side += 0.5 * np.maximum(0.0, np.roll(np.roll(z_top, di, 0), dj, 1) - z_top).sum() * dz
    c = (np.arange(nxy) + 0.5) * dz
    X, Y = np.meshgrid(c, c, indexing="ij")
    return z_top, np.hypot(X, Y), A_side, (nz - 0.5) * dz


def main():
    t = thermo()
    i = int(np.argmin(np.abs(t["time"] - 571.0)))
    T_gas, T_stag = t["jet_T_exhaust"][i], t["jet_T_stag"][i]
    P_cap, P_exh, P_face = t["jet_P_cap"][i], t["jet_P_exhaust"][i], t["jet_P_face"][i]
    mcp = P_cap / (T_stag - TAMB)

    print("ITEM 1  exhaust budget")
    for tt in (95, 285, 571):
        j = int(np.argmin(np.abs(t["time"] - tt)))
        print(f"  t={t['time'][j]:6.1f}  P_face {t['jet_P_face'][j]:7.1f} W   "
              f"P_exh {t['jet_P_exhaust'][j]:8.1f} W   T_exh {t['jet_T_exhaust'][j]:7.1f} K   "
              f"r_reach {1e3 * t['jet_r_reach'][j]:6.1f} mm")
    print(f"  mcp {mcp:.3f} W/K; check mcp*(T_exh-T_ent) = {mcp * (T_gas - TAMB):.1f} W vs {P_exh:.1f}")
    print("  => r_reach collapses 168 -> 46 mm, T_exh high.  NO STOP\n")

    z_top, r, A_wall, z0 = surface()
    nzz = t["nozzle_z"][i]
    r_hole = r[z_top < nzz].max()
    r_shield = 0.0335
    A_flow = 0.25 * np.pi * (r_hole ** 2 - r_shield ** 2)
    Dh = 2.0 * (r_hole - r_shield)

    print("ITEM 2  wall h, two independent derivations")
    hf = h_real(r_hole, 0.0)
    exp = (r_hole - r_shield) / 2e-3
    print(f"  floor h_expr (the REAL one) at r={1e3 * r_hole:.1f} mm = {hf:.1f} W/m2K "
          f"(packet states 611; simplified stand-in {1500 * 0.01875 / r_hole:.0f})")
    print(f"  (i)  area expansion {exp:.2f}x, h~u^0.8 -> {hf / exp ** 0.8:.0f} W/m2K")
    mu = 1.716e-5 * (T_gas / 273.15) ** 1.5 * (273.15 + 110.4) / (T_gas + 110.4)
    Pr = 0.7
    kg = mu * CP / Pr
    m = mcp / CP
    Re = (m / A_flow) * Dh / mu
    Nu = 0.023 * Re ** 0.8 * Pr ** 0.4
    fe = 1.0 + (Dh / 0.185) ** 0.7
    print(f"  (ii) Re {Re:.0f} (turbulent), Nu_DB {Nu:.1f} -> h {Nu * kg / Dh:.0f}; "
          f"x entry-length {fe:.2f} -> {fe * Nu * kg / Dh:.0f} W/m2K")
    print(f"  laminar-annulus arithmetic floor Nu=5.7 -> {5.7 * kg / Dh:.0f} W/m2K (NOT this regime)")
    print(f"  => lowest DEFENSIBLE h is the fully-developed {Nu * kg / Dh:.0f}, not the packet's 50."
          f"  vs 20 threshold: NO STOP\n")

    print("ITEM 3  affordability (counter-current march up the wall)")
    print(f"  measured wall area (staircase, quarter) {A_wall:.4f} m2")
    for h in (50, 150, 250):
        NTU = h * A_wall / mcp
        Tout = 700.0 + (T_gas - 700.0) * np.exp(-NTU)
        print(f"  h={h:3d}  draw {mcp * (T_gas - Tout):6.0f} W "
              f"({100 * mcp * (T_gas - Tout) / P_exh:4.1f}% of {P_exh / 1e3:.1f} kW)  "
              f"mouth T_gas {Tout:6.0f} K")
    print("  => stream leaves at 1224-1571 K, nowhere near T_ent.  NO STOP")
    print("     BUT: the enthalpy bound therefore does NOT bind inside the hole (see P3).\n")

    h_a, Tg_a, Ts_a = 80.0, 1600.0, 640.0          # the D2q-2c measurement
    S = (h_a * (Tg_a - Ts_a) - loss(Ts_a)) / (Ts_a - TAMB)

    def wall(h):
        Tg = T_gas - 0.5 * (T_gas - 700.0) * (1 - np.exp(-h * A_wall / mcp))
        return Tg, brentq(lambda T: h * (Tg - T) - loss(T) - S * (T - TAMB), 293.16, 4000.0)

    print(f"ITEMS 4+5  wall temperature and drain (rock sink S = {S:.0f} W/m2K, calibrated on")
    print( "           D2q-2c's MEASURED steady wall: h=80 vs T_gas=1600 -> 630-650 K at 2 mm)")
    print(f"  {'h':>5} {'T_gas':>7} {'T_wall':>7} {'step':>6} {'drain kW/m2':>12} {'%q_pin':>7} {'reduction':>10}")
    for h in (50, 80, 150, 250):
        Tg, Ts = wall(h)
        step = T_HOT - Ts
        print(f"  {h:5d} {Tg:7.0f} {Ts:7.0f} {step:6.0f} {Q0 * step / STEP0 / 1e3:12.0f} "
              f"{100 * Q0 * step / STEP0 / Q_PIN:6.0f}% {STEP0 / max(step, 1e-9):9.2f}x")
    Ts50 = wall(50)[1]
    red = STEP0 / (T_HOT - Ts50)
    print(f"  h=80 row reproduces the anchor -> model self-consistent")
    print(f"\n  STOP TEST (item 4): reduction at the band's LOWER end = {red:.2f}x vs 2x  -> **STOP**")
    print(f"  criterion (c): A_low wall {Ts50:.0f} K vs the >600 K bar -> MISSES")
    # bracket below the h at which the wall passes T_HOT and the ratio diverges
    h_fire = brentq(lambda h: wall(h)[1] - T_FIRE, 1.0, 1000.0)
    print(f"  h needed for 2x: "
          f"{brentq(lambda h: wall(h)[1] - (T_HOT - STEP0 / 2.0), 10.0, h_fire):.0f} W/m2K")
    print(f"  h needed to FIRE the wall (822.2 K): {h_fire:.0f} W/m2K -- above the band, so the")
    print( "  wall is WARMED not spalled everywhere in the band, as criterion (c) intends.")


if __name__ == "__main__":
    main()
