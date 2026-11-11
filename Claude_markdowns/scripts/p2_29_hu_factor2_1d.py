#!/usr/bin/env python3
"""P2-29 planner pre-registration: how the prescribed-patch factor 2 moves the Hu conduction probes (1-D).

Hu granite 2 / sandstone 2 constant properties (deck values: kappa, rho, cp), Hu's surface ramps
T_s = 638.22 + 3.74 t (granite), 684.64 + 2.46 t (sandstone), T0 = 280.15 K, 40 cells over 0.1 m (dz 2.5 mm, the
deck's mesh), insulated bottom, explicit FV. Factor 1 = AdvanceMaterial's kc*(Ts-Tc)/dz^2 (the held T a full cell
away); factor 2 = the corrected half-cell conductance (production since P2-28). 'exact' = the half-space solution
T0 + (A-T0) erfc(eta) + B t [(1+2 eta^2) erfc(eta) - (2/sqrt(pi)) eta exp(-eta^2)] at the cell-centre depth.
Probes: hu_conduction L1 at 5 mm and L2 at 10 mm depth (both on a cell face: the two neighbouring cells are shown).
Output: p2_29_hu_factor2_1d.txt. 1-D: the 30 mm patch's lateral loss is not in it.
"""
import numpy as np
from math import erfc, exp, sqrt, pi


def run(k, rho, cp, A, B, T0, fac, L=0.1, n=40, t_end=30.0):
    dz = L / n; a = k / (rho * cp); dt = 0.2 * dz * dz / a
    nt = int(np.ceil(t_end / dt)); dt = t_end / nt
    T = np.full(n, T0)
    for s in range(nt):
        Ts = A + B * s * dt
        F = np.zeros(n + 1)
        F[1:n] = k * (T[:-1] - T[1:]) / dz
        F[0] = fac * k * (Ts - T[0]) / dz
        T = T + dt * (F[:-1] - F[1:]) / (rho * cp * dz)
    return T


def exact(k, rho, cp, A, B, T0, z, t):
    a = k / (rho * cp); eta = z / (2 * sqrt(a * t))
    return T0 + (A - T0) * erfc(eta) + B * t * ((1 + 2 * eta ** 2) * erfc(eta) - (2 / sqrt(pi)) * eta * exp(-eta ** 2))


G = (2.5, 2640, 827.2, 638.22, 3.74, 280.15)
S = (4.96, 2480, 829.3, 684.64, 2.46, 280.15)
runs = {f: (run(*G, f), run(*S, f)) for f in (1, 2)}
print("t = 30 s; pct = 100 (T_sandstone - T_granite) / T_granite (hu_conduction's statistic; Hu L1 12.58, L2 18.90)")
for depth in (3.75e-3, 6.25e-3, 8.75e-3, 11.25e-3):
    i = int(round(depth / 2.5e-3 - 0.5)); out = []
    for f in (1, 2):
        Tg, Ts = runs[f][0][i], runs[f][1][i]
        out.append(f"factor {f}: Tg {Tg:6.1f} Ts {Ts:6.1f} pct {100 * (Ts - Tg) / Tg:5.2f}")
    eg, es = exact(*G, depth, 30.0), exact(*S, depth, 30.0)
    print(f"depth {depth * 1e3:5.2f} mm | " + " | ".join(out) + f" | exact Tg {eg:6.1f} Ts {es:6.1f} pct {100 * (es - eg) / eg:5.2f}")
