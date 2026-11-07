#!/usr/bin/env python3
"""P2-20 z0 (offline, report) -- every number the packet's predictions rest on, each printed with its formula.

  (1) the dielectric description: lambda(28 GHz) = c/f; n = 2.24 (DOE p. 52) => eps' = n^2; the POWER coefficient
      from the 30 % single-pass absorption in 25 mm (DOE p. 44/52): alpha = -ln(0.7)/0.025; tan_delta = alpha lambda /
      (2 pi n); R = ((n-1)/(n+1))^2; and the values the CODE computes from the declared pair eps_r 5.02, tan_delta
      0.01065 (MMWBeam.H: alpha = 2 pi sqrt(eps_r) tan_delta / lambda, R = ((sqrt eps_r - 1)/(sqrt eps_r + 1))^2)
  (2) absorbed fraction of 25 mm, 1 - exp(-alpha L); centre intensity into the rock (1 - R) 2P/(pi w0^2); the Rayleigh
      range z_R = pi w0^2 / lambda; the waist growth over the block; the Gaussian closure order (dz/z_R)^2
  (3) THE FIXTURE STRINGS: for each (alpha_old, lambda string) a beam-ON fixture carried, a tan_delta string whose
      double reproduces alpha_old BIT-FOR-BIT through the code's own expression with eps_r = 1 (sqrt(1.0) == 1.0, so
      R == 0 exactly): fl(fl(2 pi * 1.0) * fl(td / lambda)) == alpha_old (the code's order), searched +-64 ulp around alpha lambda / (2 pi).
      Python floats are IEEE doubles evaluated left to right like the C++ expression (no add, so no FMA contraction);
      lambda is the strtod of the deck's string with a unit factor of exactly 1.0 (Unit::Parse), as float() gives.
  (4) the 1-D column (Zhang's laws, Zhang.H, with effective_rhocp's smeared latent heat; adiabatic bottom; the top
      cell loses eps sigma (T^4 - T_amb^4) + h (T - T_amb), eps 0.7 below T_m 1510.5 K and 0.3 above; backward Euler,
      properties at T_old, harmonic face conductivity; the source is the code's analytic cell integral of the axial
      Gaussian (waist growth included) on the deck's dz): Z1 / Z1w / Z2 / Z3a-c -- top-cell T, peak T and its depth,
      the time the PEAK first reaches T_m.
Writes Z0.txt. Imported by apply.py for the fixture strings (FIXTURES, fixture_string).
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
C0 = 299792458.0
SIG = 5.670374419e-8

# every (alpha_old, lambda string) a beam-ON fixture or the s1 sweep carries, plus Zhang's recorded input (Z2)
FIXTURES = [(10.0, "1.0e-3"), (2.0, "1.0e-3"), (0.0, "1.0e-3"), (2000.0, "1.0e-3"),
            (2000.0, "2.19e-3"), (300.0, "2.19e-3"),
            (32000.0, "1.0e-3"), (100.0, "1.0e-3"),
            (14.0, "1.070687e-2")]


def code_alpha(eps_r, td, lam):
    """MMWBeam.H's expression in ITS evaluation order: 2 pi sqrt(eps_r) * (tan_delta / lambda)"""
    return 2.0 * math.pi * math.sqrt(eps_r) * (td / lam)


def code_R(eps_r):
    n = math.sqrt(eps_r)
    return ((n - 1.0) / (n + 1.0)) * ((n - 1.0) / (n + 1.0))


def fixture_string(alpha_old, lam_s):
    """(string, exact?, ulp offset, alpha reproduced): the shortest round-trip repr of the tan_delta double"""
    lam = float(lam_s)
    if alpha_old == 0.0:
        return "0.0", True, 0, 0.0
    td0 = alpha_old * lam / (2.0 * math.pi)
    best = None
    for k in range(0, 65):
        for sgn in ((1,) if k == 0 else (1, -1)):
            td = td0
            for _ in range(k):
                td = math.nextafter(td, math.inf if sgn > 0 else -math.inf)
            a = code_alpha(1.0, td, lam)
            if a == alpha_old:
                return repr(td), True, sgn * k, a
            err = abs(a - alpha_old)
            if best is None or err < best[0]:
                best = (err, td, sgn * k, a)
    return repr(best[1]), False, best[2], best[3]


# ------------------------------------------------------------------ Zhang's laws (Zhang.H defaults)
def rho(T):
    return 2300.0 if T >= 1273.0 else -0.2877 * (T - 273.15) + 2884.0


def kappa(T):
    Tc = T if T < 1273.0 else 1273.0
    return 1.0 / (0.3154 + 3.8e-4 * (Tc - 273.15))


def cp(T):
    th = T - 273.15
    return (-2.73e-4 * th * th + 0.894 * th + 788.0) if T <= 1510.5 else (0.184 * th + 1249.0)


def rhocp_eff(T, Tm=1510.5, fw=200.0, Lm=3.4e5):
    c = cp(T)
    if Tm - 0.5 * fw <= T <= Tm + 0.5 * fw:
        c += Lm / fw
    if 3233.0 <= T <= 3503.0:
        c += 4.8e6 / (3503.0 - 3233.0)
    return rho(T) * c


def column(P_of_t, w0, alpha, R, lam, L, N, t_end, dt=0.25, T0=297.5, Tamb=297.5, h=56.0, Tm=1510.5,
           report=(25, 50, 100, 150)):
    dz = L / N
    zR = math.pi * w0 * w0 / lam
    T = [T0] * N              # index 0 = the TOP cell
    out, onset, t = {}, None, 0.0
    rep = sorted(report)
    while t < t_end - 1e-9:
        P = P_of_t(t)
        Q = []
        for i in range(N):
            zc = (i + 0.5) * dz
            w = w0 * math.sqrt(1.0 + (zc / zR) ** 2)
            I = 2.0 * P / (math.pi * w0 * w0) * (w0 / w) ** 2
            s_top, s_bot = max(zc - 0.5 * dz, 0.0), zc + 0.5 * dz
            Q.append((1.0 - R) * I * (math.exp(-alpha * s_top) - math.exp(-alpha * s_bot)) / dz)
        eps = 0.7 if T[0] <= Tm else 0.3
        qloss = eps * SIG * (T[0] ** 4 - Tamb ** 4) + h * (T[0] - Tamb)
        a_d = [rhocp_eff(x) / dt for x in T]
        k = [kappa(x) for x in T]
        kf = [2 * k[i] * k[i + 1] / (k[i] + k[i + 1]) for i in range(N - 1)]
        A, B, Cc, D = [0.0] * N, [0.0] * N, [0.0] * N, [0.0] * N
        for i in range(N):
            B[i] = a_d[i]
            D[i] = a_d[i] * T[i] + Q[i]
            if i > 0:
                g = kf[i - 1] / dz ** 2
                A[i] = -g; B[i] += g
            if i < N - 1:
                g = kf[i] / dz ** 2
                Cc[i] = -g; B[i] += g
        D[0] -= qloss / dz
        for i in range(1, N):            # Thomas
            m = A[i] / B[i - 1]
            B[i] -= m * Cc[i - 1]
            D[i] -= m * D[i - 1]
        Tn = [0.0] * N
        Tn[-1] = D[-1] / B[-1]
        for i in range(N - 2, -1, -1):
            Tn[i] = (D[i] - Cc[i] * Tn[i + 1]) / B[i]
        T = Tn
        t += dt
        ip = max(range(N), key=lambda i: T[i])
        if onset is None and T[ip] >= Tm:
            onset = (t, (ip + 0.5) * dz, T[0])
        while rep and t >= rep[0] - 1e-9:
            out[rep.pop(0)] = (T[0], T[ip], (ip + 0.5) * dz)
    return out, onset


def main():
    L = ["z0 -- P2-20 offline derivation (every number with its formula)", "=" * 100]
    lam = C0 / 28.0e9
    n = 2.24
    a_meas = -math.log(0.7) / 0.025
    L.append("  (1) lambda = c / 28 GHz = %.7e m (deck: 1.070687e-2); n = 2.24 (DOE p. 52) => eps' = n^2 = %.4f "
             "(declared 5.02)" % (lam, n * n))
    L.append("      POWER alpha from 30 %% single-pass absorption in 25 mm: -ln(0.7)/0.025 = %.3f 1/m (DOE 0.14 Np/cm "
             "rounds it to 14); tan_delta = alpha lambda / (2 pi n) at alpha 14: %.5f (declared 0.01065); at the "
             "measured %.2f: %.5f" % (a_meas, 14.0 * lam / (2 * math.pi * n), a_meas, a_meas * lam / (2 * math.pi * n)))
    a = code_alpha(5.02, 0.01065, 1.070687e-2)
    R = code_R(5.02)
    L.append("      the CODE from the declared pair (eps_r 5.02, tan_delta 0.01065, lambda 1.070687e-2): alpha = 2 pi "
             "sqrt(eps_r) tan_delta / lambda = %.4f 1/m; R = ((sqrt eps_r - 1)/(sqrt eps_r + 1))^2 = %.4f (DOE's "
             "Fresnel 0.14); main.tex's 'low-loss' pi tan_delta sqrt(eps')/lambda = %.3f is the FIELD coefficient "
             "(x 1/2)" % (a, R, a / 2))
    frac = 1 - math.exp(-a * 0.025)
    L.append("  (2) absorbed fraction of 25 mm: 1 - exp(-alpha 0.025) = %.4f (b2 bar 0.2953 +- 0.002); transmitted "
             "%.4f" % (frac, 1 - frac))
    for w0 in (0.02, 0.01414):
        zR = math.pi * w0 * w0 / lam
        L.append("      w0 %.5f m: (1 - R) 2P/(pi w0^2) at 800 W = %.4f MW/m2; z_R = pi w0^2 / lambda = %.4f m; waist "
                 "growth over 25 mm %.2f %%; closure order (dz/z_R)^2 at dz 0.025/32 = %.2e" % (
                     w0, 1e-6 * (1 - R) * 2 * 800 / (math.pi * w0 * w0), zR,
                     100 * (math.sqrt(1 + (0.025 / zR) ** 2) - 1), (0.025 / 32 / zR) ** 2))
    L.append("  (3) the fixture strings (eps_r = 1.0 => R == 0 exactly; alpha reproduced through the code's expression):")
    for ao, ls in FIXTURES:
        s, ex, k, got = fixture_string(ao, ls)
        L.append("      alpha %-8g lambda %-11s tan_delta = %-24s %s (ulp offset %+d; alpha %r)" % (
            ao, ls, s, "BIT-EXACT" if ex else "NEAREST (not exact)", k, got))
    L.append("  (4) the 1-D column (Zhang's laws, top-cell T, peak T and depth, peak-melt onset):")
    sched28 = lambda t: 800.0 if t < 150.0 else 1800.0
    a14 = code_alpha(1.0, float(fixture_string(14.0, "1.070687e-2")[0]), 1.070687e-2)
    legs = [("Z1", sched28, 0.02, a, R, 0.025, 32, 150.0),
            ("Z1w", sched28, 0.01414, a, R, 0.025, 32, 150.0),
            ("Z2", sched28, 0.02, a14, 0.0, 0.025, 32, 150.0)]
    for w0, tag in ((0.015, "Z3a"), (0.020, "Z3b"), (0.025, "Z3c")):
        legs.append((tag, lambda t: 2000.0, w0, a, R, 0.05, 64, 400.0))
    for tag, sch, w0, al, Rr, Lz, N, tend in legs:
        rep = (25, 50, 100, 150) if tend <= 150 else (100, 200, 300, 400)
        out, onset = column(sch, w0, al, Rr, 1.070687e-2, Lz, N, tend, report=rep)
        row = "  ".join("t %3d: T_top %6.1f K, peak %6.1f K at %5.2f mm" % (tt, v[0], v[1], 1e3 * v[2])
                        for tt, v in sorted(out.items()))
        L.append("      %-4s w0 %.4f alpha %.3f R %.4f L %.3f m: %s;  peak reaches T_m: %s" % (
            tag, w0, al, Rr, Lz, row,
            "t %.1f s at %.2f mm (T_top %.1f K)" % (onset[0], 1e3 * onset[1], onset[2]) if onset else "never"))
    L.append("z0 (report)")
    t = "\n".join(L)
    print(t)
    open(os.path.join(HERE, "Z0.txt"), "w").write(t + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
