#!/usr/bin/env python3
"""P2-7 deliverable 1 -- how much heat could Hoeser's rock have taken in?

OFFLINE, ord 1 (criteria a1, a1', a2). No simulation is read.

(a) An UPPER BOUND on what her sample absorbed in 30 s at a face held at
    T_face, as the sum of three terms, each itself an upper bound:
      flakes      V rho c_p (T_face - T0)            every flake leaves at T_face
      conduction  2 k (T_face - T0) sqrt(t/(pi kappa)) A_c
                                                     a stationary face held at
                                                     T_face from t = 0 (semi-
                                                     infinite step response)
      radiation   eps sigma (T_face^4 - T0^4) A_c t
    with kappa = k/(rho c_p). The flake term is repeated with the Grimsel.H
    c_p law (rho 2650, its analytic enthalpy from 20 degC).
(b) What her calorimeter's h delivers over the Ø40 opening for two shapes, each
    normalised to a Ø9 sensing-disc mean of 7,000 W/m2K, at T_f 1,866 K.
(c) The ratio (b)/(a) and the face temperature at which it reaches 1.
(d) Her bowl's depth ratio at r = 17.5 mm read as a removal-rate ratio, turned
    into a flux ratio through her own removal exponents.

Every input is from ACTIVE_STEP.md `## Sources`; every number the packet states
about them is RE-DERIVED below with this file's own estimator and printed next
to the packet's value, so a planner slip shows up here, not in a verdict.
"""
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
SERIES = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "hoeser",
                      "hoeser_2016_drilling_series.csv")

# ---- inputs, exactly as in the packet's Sources --------------------------
T0, T_EXP, EPS = 293.15, 30.0, 0.8
SIG_SB = 5.670374419e-8
RHO, CP, K = 2706.0, 1025.0, 3.22             # Hoeser et al. 2018 sec. 3.1
KAPPA = K / (RHO * CP)
RHO_G = 2650.0                                 # Grimsel.H rho
CP_A, CP_B, CP_TH0 = 1370.0, 178000.0, 271.0   # Grimsel.H c_p law [43]
T_F = 1866.0                                   # her measured nozzle T
H_DISC = 7000.0                                # Fig 4.26, L 8 -> 7.0 kW/m2K
H_DISC_LEGS = 7100.0                           # the legs' MAP.txt value
D = 0.0079
R_OPEN = 0.020                                 # the Ø40 opening
R_DISC = 0.0045                                # the Ø9 sensing disc
R_EDGE = 0.0175
TFACES = (850.0, 1000.0, 1200.0, 1400.0)
LOADS = (("0 MPa", 12.9e-6, (2.0e-3, 2.5e-3)),
         ("47 MPa", 33.0e-6, (4.0e-3,)))
BOWL = ((5.0, 23.4), (10.0, 22.0), (15.0, 19.4), (20.0, 16.5),
        (25.0, 12.0), (28.0, 9.5))             # Fig 5.10 (depth, half-width) mm
BOWL_AXIS = 30.0


def H_grimsel(T):
    """Grimsel.H::H_of -- the exact integral of c_p from 20 degC, J/kg."""
    th = T - 273.15
    return CP_A * (th - 20.0) - CP_B * math.log((th + CP_TH0) / (20.0 + CP_TH0))


def goldstein(r):
    """Goldstein, Behbahani & Heppelmann 1986, the radial factor at fixed L."""
    return 533.0 / (533.0 + 44.0 * (r / D) ** 1.285)


def gauss05(r):
    s = 0.5 * D
    return math.exp(-r * r / (2.0 * s * s))


def disc_mean(f, R, n=200000):
    """(1/pi R^2) int_0^R f(r) 2 pi r dr, midpoint rule on n rings."""
    dr = R / n
    acc = 0.0
    for i in range(n):
        r = (i + 0.5) * dr
        acc += f(r) * r
    return 2.0 * acc * dr / (R * R)


def E_terms(T, V, A, flake="const"):
    dT = T - T0
    if flake == "const":
        e_fl = V * RHO * CP * dT
    else:
        e_fl = V * RHO_G * H_grimsel(T)
    e_cd = 2.0 * K * dT * math.sqrt(T_EXP / (math.pi * KAPPA)) * A
    e_rd = EPS * SIG_SB * (T ** 4 - T0 ** 4) * A * T_EXP
    return e_fl, e_cd, e_rd


def P_deliv(T, shape, h_disc=H_DISC):
    """Power over the Ø40 opening, W, at a uniform face T."""
    if shape == "uniform":
        return h_disc * (T_F - T) * math.pi * R_OPEN ** 2
    h0 = h_disc / disc_mean(goldstein, R_DISC)
    return h0 * disc_mean(goldstein, R_OPEN) * (T_F - T) * math.pi * R_OPEN ** 2


def closure(shape, V, A, flake):
    """T_face where delivered energy over t equals E_max; bisection."""
    f = lambda T: P_deliv(T, shape) * T_EXP - sum(E_terms(T, V, A, flake))
    lo, hi = T0 + 1.0, T_F - 1e-6
    if f(lo) < 0.0:
        return float("nan")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def bowl():
    """Volume, wetted surface and the r = 17.5 mm depth from the digitised
    Fig 5.10 profile. The profile is extended to the original surface at the
    first digitised half-width (23.4 mm at z = 0, a vertical lip -- declared)
    and closed to r = 0 at the 30 mm axis; between points r^2 is linear in z
    for the volume and r is linear in z for the surface."""
    pts = [(0.0, BOWL[0][1])] + list(BOWL) + [(BOWL_AXIS, 0.0)]
    V = S = 0.0
    for (z1, r1), (z2, r2) in zip(pts[:-1], pts[1:]):
        dz = z2 - z1
        V += math.pi * 0.5 * (r1 * r1 + r2 * r2) * dz
        S += math.pi * (r1 + r2) * math.hypot(dz, r2 - r1)
    # r = 17.5 between the digitised points (half-width falls with depth)
    z175 = None
    for (z1, r1), (z2, r2) in zip(BOWL[:-1], BOWL[1:]):
        if r1 >= R_EDGE * 1e3 >= r2:
            z175 = z1 + (r1 - R_EDGE * 1e3) / (r1 - r2) * (z2 - z1)
    return V * 1e-3, S * 1e-6, z175        # cm3, m2, mm


def n_q():
    rows = [l.strip().split(",") for l in open(SERIES)
            if l.strip() and not l.startswith(("#", "series"))]
    out = {}
    for ser in ("firing_rate", "standoff"):
        pts = sorted((float(r[6]), float(r[7])) for r in rows
                     if r[0] == ser and r[4] == "47" and r[5] == "30")
        (q1, v1), (q2, v2) = pts[0], pts[-1]
        end = math.log(v2 / v1) / math.log(q2 / q1)
        xs = [math.log(q) for q, _ in pts]
        ys = [math.log(v) for _, v in pts]
        mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
        fit = (sum((x - mx) * (y - my) for x, y in zip(xs, ys))
               / sum((x - mx) ** 2 for x in xs))
        out[ser] = (end, fit, pts)
    return out


def main():
    P = print
    L = "=" * 96
    P("P2-7 (a1, a1', a2) -- Hoeser's absorbable energy vs her calorimeter's delivery, OFFLINE")
    P(L)
    P(f"  inputs  T0 {T0} K   t {T_EXP:g} s   eps {EPS}   rho {RHO:g}  c_p {CP:g}  k {K}"
      f"  => kappa {KAPPA:.6e} m2/s")
    P(f"          T_f {T_F:g} K   h (Ø9 disc mean) {H_DISC:g} W/m2K   D {D*1e3:g} mm"
      f"   opening Ø{2*R_OPEN*1e3:g} mm")
    P(f"          Grimsel.H flake variant: rho {RHO_G:g}, H(T) = {CP_A:g}(th-20)"
      f" - {CP_B:g} ln((th+{CP_TH0:g})/{20+CP_TH0:g})   (analytic, unclamped)")
    P(f"          conduction depth sqrt(t/(pi kappa)) = {math.sqrt(T_EXP/(math.pi*KAPPA)):.2f} s^0.5/m"
      f"  (diffusion length 2 sqrt(kappa t) = {2*math.sqrt(KAPPA*T_EXP)*1e3:.2f} mm)")

    # ---- re-derived packet numbers ---------------------------------------
    P("\n  RE-DERIVED SOURCE NUMBERS (this file's estimator | the packet's value)")
    for r_mm, want in ((10, 0.899), (15, 0.842), (17.5, 0.813), (20, 0.786)):
        g = goldstein(r_mm * 1e-3)
        P(f"    Goldstein g({r_mm:>4} mm) = {g:.4f}  | {want}   {'ok' if abs(g-want) < 6e-4 else 'DIFFERS'}")
    gd = disc_mean(goldstein, R_DISC)
    P(f"    Goldstein Ø9 disc mean   = {gd:.6f}  | 0.976294   {'ok' if abs(gd-0.976294) < 1e-6 else 'DIFFERS'}")
    H0G = H_DISC_LEGS / gd
    P(f"    H0G = 7100 / disc mean   = {H0G:.4f}  | 7272.3958   {'ok' if abs(H0G-7272.3958) < 5e-5 else 'DIFFERS'}")
    go = disc_mean(goldstein, R_OPEN)
    P(f"    Goldstein Ø40 area mean  = {go:.6f}  => H0G-shape mean over the opening"
      f" {H0G*go:.1f} W/m2K = {H0G*go/H_DISC_LEGS:.4f} x the uniform 7100")
    Vb, Sb, z175 = bowl()
    P(f"    her 47 MPa bowl volume   = {Vb:.2f} cm3 | 33.0 (= 65.7 cm3/min x 0.5 min)"
      f"   {'ok' if abs(Vb-33.0)/33.0 < 0.05 else 'DIFFERS'}")
    P(f"    her 47 MPa bowl surface  = {Sb:.3e} m2 (wall, revolved) | 4.0e-3"
      f"   {'ok' if abs(Sb-4.0e-3)/4.0e-3 < 0.10 else 'DIFFERS'}")
    P(f"    depth at r = 17.5 mm     = {z175:.2f} mm | 18.3 ; / axis {BOWL_AXIS:g} mm"
      f" = {z175/BOWL_AXIS:.4f} | 0.61")
    nq = n_q()
    for ser, want in (("standoff", 1.44), ("firing_rate", 1.83)):
        end, fit, pts = nq[ser]
        P(f"    n_q {ser:<11s} endpoints = {end:.4f} | {want}   (log-log LSQ over"
          f" {len(pts)} points {fit:.4f})   {'ok' if abs(end-want) < 6e-3 else 'DIFFERS'}")
    P(f"    her 0 MPa 30 s volume    = 25.8 cm3/min x 0.5 = {25.8*0.5:.1f} cm3 | 12.9")

    # ---- (a) absorbable --------------------------------------------------
    P("\n  (a) ABSORBABLE ENERGY over 30 s -- an UPPER BOUND, every term")
    P(f"    {'load':>6s} {'A_c [m2]':>9s} {'T_face':>7s} | {'flake':>8s} {'flakeG':>8s}"
      f" {'cond':>8s} {'rad':>7s} | {'E_max':>8s} {'E_maxG':>8s} [kJ] | {'P_max':>6s} {'P_maxG':>6s} [kW]")
    for name, V, As in LOADS:
        for A in As:
            for T in TFACES:
                fl, cd, rd = E_terms(T, V, A, "const")
                flg = E_terms(T, V, A, "grimsel")[0]
                E, EG = fl + cd + rd, flg + cd + rd
                P(f"    {name:>6s} {A:9.1e} {T:7.0f} | {fl/1e3:8.2f} {flg/1e3:8.2f} {cd/1e3:8.2f}"
                  f" {rd/1e3:7.2f} | {E/1e3:8.2f} {EG/1e3:8.2f}      | {E/T_EXP/1e3:6.3f}"
                  f" {EG/T_EXP/1e3:6.3f}")
    P("    flakeG / E_maxG: the flake term on the Grimsel.H c_p law (rho 2650).")

    # ---- (b) delivered ---------------------------------------------------
    P("\n  (b) CALORIMETER-IMPLIED DELIVERY over the Ø40 opening (Ø9 disc mean"
      f" {H_DISC:g} W/m2K, T_f {T_F:g} K)")
    P(f"    {'T_face':>7s} | {'uniform [kW]':>12s} {'Goldstein [kW]':>15s} |"
      f" {'uniform [kJ/30s]':>17s} {'Goldstein [kJ/30s]':>19s}")
    for T in TFACES:
        pu, pg = P_deliv(T, "uniform"), P_deliv(T, "goldstein")
        P(f"    {T:7.0f} | {pu/1e3:12.3f} {pg/1e3:15.3f} | {pu*T_EXP/1e3:17.2f}"
          f" {pg*T_EXP/1e3:19.2f}")
    P(f"    (at the legs' 7,100: x{H_DISC_LEGS/H_DISC:.4f} on every entry --"
      f" uniform {P_deliv(850,'uniform',H_DISC_LEGS)/1e3:.3f} kW at 850 K)")

    # ---- (c) ratio and closure ------------------------------------------
    P("\n  (c) RATIO delivered / absorbable over 30 s, and the CLOSURE T_face (ratio = 1)")
    P(f"    {'load':>6s} {'A_c':>8s} {'flake law':>10s} {'shape':>10s} | "
      + " ".join(f"{'@'+str(int(T)):>7s}" for T in TFACES) + f" | {'closure T [K]':>13s}")
    summ = {}
    for name, V, As in LOADS:
        for A in As:
            for fl in ("const", "grimsel"):
                for sh in ("uniform", "goldstein"):
                    rs = [P_deliv(T, sh) * T_EXP / sum(E_terms(T, V, A, fl)) for T in TFACES]
                    tc = closure(sh, V, A, fl)
                    summ[(name, A, fl, sh)] = (rs[0], tc)
                    P(f"    {name:>6s} {A:8.1e} {fl:>10s} {sh:>10s} | "
                      + " ".join(f"{x:7.3f}" for x in rs) + f" | {tc:13.1f}")

    # ratio at 870 K (the a1 bar's T_face <= 870 K)
    P("\n    at T_face = 870 K (the a1 bar's upper face T):")
    worst0 = 1e30
    for name, V, As in LOADS:
        for A in As:
            for fl in ("const", "grimsel"):
                for sh in ("uniform", "goldstein"):
                    r870 = P_deliv(870.0, sh) * T_EXP / sum(E_terms(870.0, V, A, fl))
                    if name == "0 MPa":
                        worst0 = min(worst0, r870)
                    P(f"      {name:>6s} A_c {A:.1e} {fl:>8s} {sh:>10s}  ratio {r870:6.3f}")

    # ---- (d) the edge flux ratio ---------------------------------------
    P("\n  (d) HER IMPLIED EDGE FLUX RATIO at r = 17.5 mm")
    dr = z175 / BOWL_AXIS
    P(f"    depth ratio (bowl, 47 MPa)   {dr:.4f}  (first order: a depth ratio read as a"
      f" removal-RATE ratio at the original radius)")
    imp = {}
    for ser in ("standoff", "firing_rate"):
        n = nq[ser][0]
        imp[ser] = dr ** (1.0 / n)
        P(f"    flux ratio = {dr:.4f}^(1/{n:.4f}) [{ser}] = {imp[ser]:.4f}")
    gg, ga = goldstein(R_EDGE), gauss05(R_EDGE)
    P(f"    Goldstein g(17.5 mm)         {gg:.4f}")
    P(f"    Gaussian sigma = 0.5 D       {ga:.3e}   (implied / Gaussian = "
      f"{min(imp.values())/ga:.3e} .. {max(imp.values())/ga:.3e})")

    # ---- verdict rows ---------------------------------------------------
    P("\n" + L)
    P("  CRITERIA (ord 1)")
    r0 = [summ[k][0] for k in summ if k[0] == "0 MPa"]
    t0 = [summ[k][1] for k in summ if k[0] == "0 MPa"]
    r47 = [summ[k][0] for k in summ if k[0] == "47 MPa"]
    t47 = [summ[k][1] for k in summ if k[0] == "47 MPa"]
    a1 = worst0 >= 3.0 and min(t0) >= 1300.0
    P(f"  a1   0 MPa @850 K ratio {min(r0):.2f}-{max(r0):.2f} (all A_c, both flake laws, both"
      f" shapes); min @870 K {worst0:.2f} (bar >= 3);")
    P(f"       closure T_face {min(t0):.0f}-{max(t0):.0f} K (bar >= 1,300)   -> "
      f"{'PASS' if a1 else 'MISS'}")
    P(f"       P1 predicted 5.0-6.5x and closure 1,400-1,450 K")
    a1p = max(r47) < min(r0)
    P(f"  a1'  47 MPa @850 K ratio {min(r47):.2f}-{max(r47):.2f}; closure {min(t47):.0f}-{max(t47):.0f} K;"
      f" below the 0 MPa ratio {min(r0):.2f}-{max(r0):.2f}  -> {'PASS' if a1p else 'MISS'}")
    P(f"       P1' predicted ~3.2 kW absorbable => 2.5-2.8x")
    lo, hi = min(imp.values()), max(imp.values())
    a2 = 0.55 <= lo and hi <= 0.85 and lo >= 1e3 * 5.5e-5
    P(f"  a2   implied edge flux ratio {lo:.3f}-{hi:.3f} (bar 0.55-0.85 and >= 1e3 x 5.5e-5 ="
      f" 0.055)  -> {'PASS' if a2 else 'MISS'}")
    P(f"       Goldstein {gg:.3f}; Gaussian 0.5 D {ga:.2e}; P2 predicted 0.71-0.76")
    P(L)


if __name__ == "__main__":
    main()
