#!/usr/bin/env python3
"""P2-21 step zero (2026-11-06): the deposition closure of the packet's tau kernel, offline, on the #19 grid.

The packet's kernel (Goal 3): each cell receives (1 - R) Pi_k (e^{-tau_top,k} - e^{-tau_bot,k}) / dz, with Pi_k =
BeamIncidentIntensity at the cell's OWN centre depth (untouched, guardrail 1), tau summed down the column. The closure
counter (BeamEnergyBalance) compares P_dep + P_trans with P_inc = sum of (1 - R) Pi at each column's TOP cell.

Per column the tau pair telescopes only if Pi is the same down the column. A Gaussian's Pi changes with depth (the
waist grows), by a relative (d/zR)^2 (2 r^2/w0^2 - 1): down on the axis, up beyond w0/sqrt2. For a laterally UNIFORM
alpha the lateral sum of Pi at every depth is P, so the changes cancel (P2-20 b1: 6.5e-7). Once the centre runs
hotter than the edge (the table's whole point) they no longer cancel. This script measures the size of that term
with the packet's kernel ("per-cell Pi") and with Pi taken once per column at the top ("column Pi"), on the T1 grid
(128 x 128 x 64 over 0.1 x 0.1 x 0.05 m), the Hartlieb table on alpha_cold = 14.003 1/m, P = 1 W, R = 0.

The temperature field is SYNTHETIC (a Gaussian hot zone exp(-2 r^2/w_T^2) exp(-d/ell) of peak rise dT on 297.5 K):
it sets the lateral contrast of alpha, which is what the term needs; the run's own field would give other digits,
not another order. Usage: closure_step0.py > CLOSURE_STEP0.txt
"""
import math

import numpy as np

NX, NZ, LX, LZ = 128, 64, 0.10, 0.05
DX, DZ = LX / NX, LZ / NZ
LAM, A0 = 1.070687e-2, 14.003
# Hartlieb granite, 2467 MHz heating branch, tan d / tan d(22 C); T in K (packet Sources)
TC = [22, 201, 300, 401, 451, 500, 550, 575, 600, 649, 699, 748, 797, 847, 896, 945, 995]
TK = np.array([t + 273.15 for t in TC])
RR = np.array([1.00, 0.98, 1.15, 1.39, 1.62, 2.01, 2.42, 2.73, 3.15, 4.27, 5.58, 5.72, 4.74, 3.79, 3.86, 6.02, 8.88])

xc = (np.arange(NX) + 0.5) * DX - 0.5 * LX
X, Y = np.meshgrid(xc, xc, indexing="ij")
R2 = X ** 2 + Y ** 2
D = (np.arange(NZ) + 0.5) * DZ                    # cell-centre depth below the top face


def closure(w0, T, column_pi):
    zR = math.pi * w0 ** 2 / LAM
    alpha = A0 * np.interp(T, TK, RR)
    tau_bot = np.cumsum(alpha * DZ, axis=2)
    tau_top = tau_bot - alpha * DZ
    w = w0 * np.sqrt(1.0 + (D / zR) ** 2)
    Pi = (2.0 / (math.pi * w0 ** 2)) * (w0 / w) ** 2 * np.exp(-2.0 * R2[..., None] / w ** 2)
    if column_pi:
        Pi = np.repeat(Pi[..., :1], NZ, axis=2)
    A = DX * DX
    inc = Pi[..., 0].sum() * A
    dep = (Pi * (np.exp(-tau_top) - np.exp(-tau_bot))).sum() * A
    tr = (Pi[..., -1] * np.exp(-tau_bot[..., -1])).sum() * A
    return (dep + tr - inc) / inc, dep / inc, tr / inc


def table(w0, wT, ell, rises):
    zR = math.pi * w0 ** 2 / LAM
    print("w0 %.1f mm (zR %.1f mm), hot zone w_T %.1f mm, ell %.0f mm" % (1e3 * w0, 1e3 * zR, 1e3 * wT, 1e3 * ell))
    print("   %8s %9s | %-34s | %s" % ("rise [K]", "T_top [K]", "per-cell Pi (the packet's kernel)",
                                      "column Pi (closure only)"))
    for dT in rises:
        T = 297.5 + dT * np.exp(-2.0 * R2[..., None] / wT ** 2) * np.exp(-D / ell)[None, None, :]
        c, fd, ft = closure(w0, T, False)
        c2, _, _ = closure(w0, T, True)
        print("   %8d %9.1f | closure %+.2e  dep %.4f trans %.4f | closure %+.1e" % (
            dT, 297.5 + dT * math.exp(-D[0] / ell), c, fd, ft, c2))


if __name__ == "__main__":
    print(__doc__.split("\n\n")[0])
    print("b1's bar: |beam_closure_err| <= 1e-6 on every T1b step row\n")
    table(0.0165, 0.0165, 0.006, (0, 100, 190, 250, 400, 700, 1213))      # T1b's spot
    print()
    table(0.020, 0.020 / math.sqrt(2.0), 0.006, (0, 500, 1213))           # T1c's spot, narrower hot zone
    table(0.020, 0.020, 0.010, (500, 1213))                               # T1c's spot, wider and deeper
