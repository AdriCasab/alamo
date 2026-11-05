#!/usr/bin/env python3
"""P2-7 (a3) common-mode check: the bracket for HG/HA step-1 patch_P_robin.

A parse fault shared by the shipped and the equivalent HG strings would pass
the pair. The HG/HA power ratio at step 1 must instead sit between
  (i)  sum h_G / sum h_A over the in-opening column centres (P proportional to h), and
  (ii) the same with a half-cell Robin conductance h k'/(k'+h),
       k' = lambda_Grimsel(293.15 K) / (dz/2), which flattens the ratio toward 1.
A dropped radial factor would read 7272.3958/7100 = 1.0243 or above.
"""
import math

D, DX, NX = 0.0079, 0.001, 50
LAM = 2000.0 / (20.0 + 563.0)          # Grimsel.H lambda at 20 degC
KP = LAM / (0.5 * DX)
g = lambda r: 533.0 / (533.0 + 44.0 * (r / D) ** 1.285)
sa = sg = ra = rg = 0.0
n = 0
for i in range(NX):
    for j in range(NX):
        r = math.hypot((i + 0.5) * DX, (j + 0.5) * DX)
        if r > 0.02:
            continue
        n += 1
        ha, hg = 7100.0, 7272.3958 * g(r)
        sa += ha; sg += hg
        ra += ha * KP / (KP + ha); rg += hg * KP / (KP + hg)
print(f"columns in the opening {n}; k' = {KP:.1f} W/m2K")
print(f"(i)  P ~ h          : {sg/sa:.5f}")
print(f"(ii) half-cell Robin: {rg/ra:.5f}")
print(f"dropped radial factor would read >= {7272.3958/7100:.5f}")
