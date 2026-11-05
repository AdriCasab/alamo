#!/usr/bin/env python3
"""P2-2e cross-check (advisor, 2026-10-19): delamination energy release of a surface-parallel crack.

Layer of thickness t above a crack; misfit = layer mean minus substrate mean over [t, t+W] (W = flake width).
sigma_m = R*E*beta*(Tbar_layer - Tbar_sub)/(1-nu);  G(t) = (1-nu^2)*sigma_m^2*t/(2E)  (Hutchinson & Suo,
steady-state delamination of a misfitting layer). Fires when max_t G >= G_c. G_c is CALIBRATED so that a wide
W (100 mm, the thin-skin limit) fires Kant 2016 SOD 80 at T_s = 843.76 K (the banked level), then W is laddered
against the SOD 130 (spall) / 150 (no spall) bracket. Exact Robin half-space erfc profile, constant k, kappa.
"""
from math import erfc, exp, sqrt
k, kap, T0 = 2.5, 1.1e-6, 293.15
E, beta, nu, R = 50e9, 8e-6, 0.25, 0.531294
SOD = {50: (766, 1723.15), 80: (1009, 1683.15), 110: (873, 1633.15), 130: (747, 1523.15), 150: (633, 1373.15),
       "h500": (500, 1673.15), "h150": (150, 1673.15)}
def th(z, t, h):
    s = sqrt(kap*t); b = h*s/k; u = z/(2*s); return erfc(u) - exp(h*z/k + b*b)*erfc(u + b)
def mean(f, a, b, n=40):
    dz = (b-a)/n; return sum(f(a+(i+0.5)*dz) for i in range(n))/n
TS = [x*1e-4 for x in range(1, 31)]          # crack depth 0.1-3.0 mm
def Gmax(time, h, Tf, W):
    f = lambda z: (Tf-T0)*th(z, time, h); best = 0.0
    for t in TS:
        sm = R*E*beta*(mean(f, 0, t, 10) - mean(f, t, t+W))/(1-nu)
        if sm > 0: best = max(best, (1-nu*nu)*sm*sm*t/(2*E))
    return best
def fire(h, Tf, W, Gc, tmax=200):
    t = 0.0; g = 0.0
    while t < tmax:
        t += 0.05 if t < 5 else 0.5
        g = max(g, Gmax(t, h, Tf, W))
        if g >= Gc: return t, T0 + (Tf-T0)*th(0, t, h)
    return None, g/Gc
# calibrate Gc at W = 100 mm: SOD 80 fires when T_s reaches 843.76
h, Tf = SOD[80]; t = 0.0
while T0 + (Tf-T0)*th(0, t+1e-3, h) < 843.76: t += 0.005
Gc = Gmax(t, h, Tf, 0.1)
print(f"calibrated G_c = {Gc:.1f} J/m^2 at SOD 80, t = {t:.3f} s (thin-skin limit, W = 100 mm)")
print("W [mm] | onset T_s [K] or max G/G_c")
for W in (3, 4, 5, 6, 8, 10, 12, 15, 20, 30):
    row = []
    for s in (80, 110, 130, 150, "h500", "h150"):
        tt, v = fire(*SOD[s], W*1e-3, Gc)
        row.append(f"{s}:{'%.0fK' % v if tt else 'no(%.2f)' % v}")
    print(f"  {W:4.0f}  " + "  ".join(row))
