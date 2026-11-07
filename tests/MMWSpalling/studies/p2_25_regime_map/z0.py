#!/usr/bin/env python3
"""P2-25 z0 -- the offline numbers, each with its formula, and the gate that fixes the legs' argv.
Usage: z0.py > Z0.txt   (writes STOPS.json; exit 1 on a failed bar)

Constants are READ: beam.eps_r and beam.tan_delta from validation/shared/rock_grimsel.inputs, beam.lambda and
beam.omega0 from the regime deck, the centre-flux relation from BeamIncidentIntensity (Pi = 2P/(pi w0^2) (w0/w)^2
e^{-2r^2/w^2}: at the centre of the surface w = w0, so the absorbed centre flux is q = (1 - R) 2P/(pi w0^2)).
The Grimsel law is transcribed from src/Numeric/Material/Grimsel.H (rho 2650; c_p = 1370 - 178000/(theta + 271);
lambda = 2000/(theta + 563); theta in C).

THE 1-D COLUMN (L11's t_fire_1D): 100 cells of 0.5 mm (the deck's 50 mm), explicit, the Beer-Lambert source as each
cell's exact integral q (e^{-alpha z_top} - e^{-alpha z_bot}) / dz (the power past 50 mm leaves), Grimsel lambda(T),
c_p(T), rho; the bottom adiabatic; the surface temperature T_s from the face balance lambda (T_1 - T_s)/(dz/2) =
eps sigma (T_s^4 - T_amb^4) + h (T_s - T_amb) with eps 0.8, h 56, T_amb 297.5 K (the #19 deck's, robin_form = surface:
the face's own T). Run until T_s >= 843.757 K (the onset level, P2-2d) or 600 s; dt = min(0.4 dz^2/kappa_max, 2 K over
the top cell's heating rate). The stop rule: stop = clamp(3 t_fire_1D, floor, 600 s), floor 60 s (q <= 2) / 20 s
(q >= 4); never by 600 s => 600 s and pre-registered C. A 1-D centre column ignores lateral conduction and the beam's
spread: it is the stop rule's clock, not a prediction of the map.
"""
import hashlib
import json
import math
import os
import re
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
ROCK = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "shared", "rock_grimsel.inputs")
DECK = os.path.join(HERE, "input")
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
BIN24 = os.path.join(ROOT, "tests", "MMWSpalling", "studies", "p2_24_melt_vapour", "bin_p2_24", "mmwspalling-3d-g++")
SIG, EPS, H_CONV, T_AMB, T0, T_ON = 5.670374419e-8, 0.8, 56.0, 297.5, 297.5, 843.757
RHO = 2650.0
DEPTHS = [("0p5", 0.5e-3), ("1", 1e-3), ("2", 2e-3), ("4", 4e-3), ("8", 8e-3), ("16", 16e-3), ("71", None)]
FLUXES = [("0p3", 0.3), ("0p5", 0.5), ("1", 1.0), ("2", 2.0), ("4", 4.0), ("8", 8.0)]
# the packet's tables (ACTIVE_STEP.md "The ladders"), to be reproduced to the printed digit
PK_TD = {"0p5": "1.52111", "1": "0.760555", "2": "0.380278", "4": "0.190139", "8": "0.0950694", "16": "0.0475347",
         "71": None}
PK_P = {"0p3": "91.6", "0p5": "152.7", "1": "305.3", "2": "610.7", "4": "1221.3", "8": "2442.7"}
PK_DT = {"0p5": "0.06 0.03 0.01 0.009 0.004 0.002", "1": "0.08 0.04 0.02 0.01 0.006 0.003",
         "2": "0.1 0.07 0.03 0.01 0.008 0.004", "4": "0.1 0.1 0.05 0.02 0.01 0.007", "8": "0.2 0.2 0.1 0.05 0.02 0.01",
         "16": "0.2 0.2 0.2 0.1 0.05 0.02", "71": "0.2 0.2 0.2 0.2 0.2 0.1"}
FAIL = []


def need(ok, msg):
    print(("  PASS  " if ok else "  FAIL  ") + msg)
    if not ok:
        FAIL.append(msg)


def read_keys(path):
    kv = {}
    for ln in open(path):
        ln = ln.split("#")[0].strip()
        if "=" in ln and not ln.startswith("FILE"):
            k, v = ln.split("=", 1)
            kv[k.strip()] = v.strip()
    return kv


def cp(T):
    return 1370.0 - 178000.0 / ((T - 273.15) + 271.0)


def lam(T):
    return 2000.0 / ((T - 273.15) + 563.0)


def dt_rule(q, alpha, dz, rhocp):
    r = q * 1e6 * (1.0 - math.exp(-alpha * dz)) / (rhocp * dz)
    dt = min(0.2, 6.0 / r)
    e = math.floor(math.log10(dt))
    return float("%.1g" % (math.floor(dt / 10 ** e) * 10 ** e)), r


def surface_T(T1, dz):
    """the face balance lambda(T1)(T1 - Ts)/(dz/2) = eps sigma (Ts^4 - Ta^4) + h (Ts - Ta), Newton from T1"""
    g = lam(T1) / (0.5 * dz)
    Ts = T1.copy() if isinstance(T1, np.ndarray) else T1
    for _ in range(30):
        f = g * (T1 - Ts) - EPS * SIG * (Ts ** 4 - T_AMB ** 4) - H_CONV * (Ts - T_AMB)
        df = -g - 4.0 * EPS * SIG * Ts ** 3 - H_CONV
        Ts = Ts - f / df
    return Ts


def column(q, alpha, t_max=600.0, dz=0.5e-3, n=100):
    """t_fire_1D: the first time the face T reaches T_ON, or None by t_max"""
    z = dz * np.arange(n + 1)
    S = q * 1e6 * (np.exp(-alpha * z[:-1]) - np.exp(-alpha * z[1:])) / dz        # W/m3
    T = np.full(n, T0)
    kmax = lam(T0) / (RHO * cp(T0))
    dt_c = 0.4 * dz * dz / kmax
    r_top = S[0] / (RHO * cp(T0))
    dt = min(dt_c, 2.0 / r_top) if r_top > 0 else dt_c
    t = 0.0
    while t < t_max:
        Ts = surface_T(T[0], dz)
        if Ts >= T_ON:
            return t, Ts
        L = lam(T)
        Lf = 2.0 * L[:-1] * L[1:] / (L[:-1] + L[1:])                              # harmonic face conductivity
        F = Lf * (T[1:] - T[:-1]) / dz                                              # flux INTO cell k from k+1 [W/m2]
        div = np.zeros(n)
        div[:-1] += F
        div[1:] -= F
        div[0] -= lam(T[0]) * (T[0] - Ts) / (0.5 * dz)                               # the face loss through the top
        T = T + dt * (div / dz + S) / (RHO * cp(T))
        t += dt
    return None, surface_T(T[0], dz)


def main():
    rk, dk = read_keys(ROCK), read_keys(DECK)
    er, td_rock = float(rk["beam.eps_r"]), float(rk["beam.tan_delta"])
    lam_b, w0 = float(dk["beam.lambda"]), float(dk["beam.omega0"])
    n_idx = math.sqrt(er)
    R = ((n_idx - 1) / (n_idx + 1)) ** 2
    Pper = math.pi * w0 * w0 / (2.0 * (1.0 - R))
    nz = int(dk["amr.n_cell"].split()[2])
    Lz = float(dk["geometry.prob_hi"].split()[2]) - float(dk["geometry.prob_lo"].split()[2])
    dz = Lz / nz
    rhocp = RHO * cp(T0)
    a_rock = 2.0 * math.pi * math.sqrt(er) * (td_rock / lam_b)
    print("P2-25 z0 -- the ladders, the boundary arithmetic, the stop rule (units: W/m2 printed MW/m2, W, 1/m, mm, s, K)")
    print("=" * 112)
    print("constants READ: rock eps_r %g (n %.6f, R %.6f, 1-R %.5f), tan_delta %g; deck lambda %g m, omega0 %g m, "
          "dz %.6g m (%d cells over %g m); rho c_p(%.1f K) = %.4g J/m3K" % (er, n_idx, R, 1 - R, td_rock, lam_b, w0, dz,
                                                                       nz, Lz, T0, rhocp))
    print("q = (1 - R) 2P/(pi w0^2) => P = q x pi w0^2 / (2 (1 - R)) = q x %.5e W per W/m2 (packet: 3.0533e-4)" % Pper)
    need("%.4e" % Pper == "3.0533e-04", "the P-per-flux factor %.4e = the packet's 3.0533e-4" % Pper)

    # 1 the ladders
    print("\n1. THE LADDERS (alpha = 2 pi sqrt(eps_r) tan_delta / lambda, the code's association, MMWBeam.H)")
    alpha = {}
    for tok, d in DEPTHS:
        if d is None:
            alpha[tok] = a_rock
            print("   1/alpha %7.3f mm: alpha %.4f 1/m = the rock file's tan_delta %g (no token)" % (1e3 / a_rock, a_rock,
                                                                                               td_rock))
            continue
        a = 1.0 / d
        td = a * lam_b / (2.0 * math.pi * n_idx)
        a_tok = 2.0 * math.pi * math.sqrt(er) * (float(PK_TD[tok]) / lam_b)
        alpha[tok] = a_tok
        need("%g" % float("%.6g" % td) == PK_TD[tok] and "%g" % a_tok == "%g" % a,
             "1/alpha %g mm: alpha %g, tan_delta %.9g -> token %s (the packet's %s); the token gives alpha %.6f (%%g %s)"
             % (d * 1e3, a, td, "%.6g" % td, PK_TD[tok], a_tok, "%g" % a_tok))
    P = {}
    for tok, q in FLUXES:
        P[tok] = round(q * 1e6 * Pper, 1)
        need("%.1f" % P[tok] == PK_P[tok], "q %g MW/m2 -> P %.1f W (the packet's %s)" % (q, P[tok], PK_P[tok]))
    DT = {}
    for dtok, _ in DEPTHS:
        row = []
        for qtok, q in FLUXES:
            dt, r = dt_rule(q, alpha[dtok], dz, rhocp)
            DT[(dtok, qtok)] = dt
            row.append("%g" % dt)
        need(" ".join(row) == PK_DT[dtok], "dt row 1/alpha %s: %s (the packet's %s)" % (dtok, " ".join(row), PK_DT[dtok]))

    # 2 the threshold in absorbed terms
    Ts17, q17, q17n = 819.81, 0.4316e6, 0.4368e6
    rad = EPS * SIG * (Ts17 ** 4 - T_AMB ** 4)
    conv = H_CONV * (Ts17 - T_AMB)
    print("\n2. THE MINIMUM FLUX IN ABSORBED TERMS. P2-17 (studies/p2_17_buckle_on/offline.py:15) defines q_net = h (T_g - "
          "T_s) - eps sigma (T_s^4 - T_amb^4) - h_conv (T_s - T_amb): radiation and the CAL deck's h_conv 10 are "
          "SUBTRACTED. At the 12 mm span q_net,min 0.4316 MW/m2 at T_s 819.81 K (0.4368 without h_conv) is the net INTO "
          "the solid. Under the beam the same net needs an absorbed q = q_net + eps sigma (T_s^4 - 297.5^4) + 56 (T_s - "
          "297.5) = %.4f + %.4f + %.4f = %.4f MW/m2 (on the 0.4368 reading: %.4f)" % (
              q17 / 1e6, rad / 1e6, conv / 1e6, (q17 + rad + conv) / 1e6, (q17n + rad + conv) / 1e6))
    print("   => the beam's threshold is %.2f MW/m2 ABSORBED (the packet's range 0.43-0.49: the upper end; the q 0.5 column "
          "sits %.2fx above it, the 0.3 column below)" % ((q17 + rad + conv) / 1e6, 0.5e6 / (q17 + rad + conv)))

    # 3 the predicted boundary (P1)
    kap, lamb, dT = 1.0e-6, 2.2, 550.0
    print("\n3. THE PREDICTED BOUNDARY (P1): t_on = (pi / 4 kappa)(lambda dT / q)^2 for surface heating (lambda 2.2, dT 550 "
          "K, kappa 1e-6), L_on = max(1/alpha, 1.1 sqrt(kappa t_on)) against B's layer limit 2.35-2.4 mm at ~850 K")
    for qtok, q in FLUXES:
        ton = math.pi / (4 * kap) * (lamb * dT / (q * 1e6)) ** 2
        lc = 1.1 * math.sqrt(kap * ton)
        cells = " ".join("%s:%.2f" % (dtok, 1e3 * max(1.0 / alpha[dtok], lc)) for dtok, _ in DEPTHS)
        print("   q %-4g t_on %7.2f s, 1.1 sqrt(kappa t_on) %.2f mm; L_on [mm] per depth %s" % (q, ton, 1e3 * lc, cells))

    # 4 the stop rule
    print("\n4. t_fire_1D AND THE STOP (the 1-D column above; stop = clamp(3 t_fire_1D, floor, 600 s))")
    stops, tfire = {}, {}
    total = 0
    for dtok, _ in DEPTHS:
        line = []
        for qtok, q in FLUXES:
            tf, Ts = column(q, alpha[dtok])
            floor = 60.0 if q <= 2.0 else 20.0
            st = 600.0 if tf is None else min(600.0, max(floor, math.ceil(3.0 * tf)))
            stops["D%s_Q%s" % (dtok, qtok)] = st
            tfire["D%s_Q%s" % (dtok, qtok)] = tf
            steps = int(round(st / DT[(dtok, qtok)]))
            total += steps
            line.append("%s:%s->%g s/%d" % (qtok, "%.2f" % tf if tf is not None else "never(Ts %.0f)" % Ts, st, steps))
        print("   1/alpha %-3s %s" % (dtok, "  ".join(line)))
    need(all(60.0 <= stops["D%s_Q%s" % (d, qt)] <= 600.0 if q <= 2 else 20.0 <= stops["D%s_Q%s" % (d, qt)] <= 600.0
             for d, _ in DEPTHS for qt, q in FLUXES),
         "every stop >= its floor and <= 600 s; %d steps over the 42 legs" % total)
    nC = sum(1 for v in tfire.values() if v is None)
    print("   pre-registered C (the 1-D column never reaches onset by 600 s): %d cells: %s" % (
        nC, " ".join(k for k, v in tfire.items() if v is None) or "none"))

    # 5 the binary
    h = hashlib.sha256(open(BIN, "rb").read()).hexdigest()
    newer = subprocess.run(["find", "src", "-type", "f", "-newer", BIN24], cwd=ROOT, capture_output=True,
                           text=True).stdout.split()
    print("\n5. THE BINARY: bin/mmwspalling-3d-g++ sha256 %s; sources newer than the P2-24 binary: %s" % (h, newer))
    need(h.startswith("c052c348") and newer == ["src/Numeric/Material/Grimsel.H"],
         "sha c052c348..., only Grimsel.H newer")
    json.dump({"stops": stops, "t_fire_1D": tfire, "dt": {"D%s_Q%s" % k: v for k, v in DT.items()}, "P": P,
               "td": PK_TD}, open(os.path.join(HERE, "STOPS.json"), "w"), indent=1)
    print("=" * 112)
    print("z0 %s (%d failure(s)); STOPS.json written" % ("PASS" if not FAIL else "FAIL", len(FAIL)))
    return 0 if not FAIL else 1


if __name__ == "__main__":
    sys.exit(main())
