#!/usr/bin/env python3
"""P2-2e deliverable 1(a) -- the OFFLINE constraint-depth ladder, pre-registration.

THE MODEL. Exact Robin half-space solution with constant properties:
    T(z,t) - T0 = (T_f - T0) * [ erfc(u) - exp(h z/k + b^2) * erfc(u + b) ]
    u = z / (2 sqrt(kappa t)),   b = h sqrt(kappa t) / k
The constraint-depth criterion fires when the DRIVE reaches the banked level:
    max_t [ T(0,t) - T(d,t) ]  >=  dT_c = 843.757 - 293.15 = 550.607 K
d = 0 recovers today's absolute criterion (T(0) - T0 >= dT_c) because T(inf) = T0.
In steady state T(0) - T(d) = q d / k, so firing needs q >= k dT_c / d: the
MINIMUM FLUX the law exists to produce, and d_ref ~ k dT_c / q_min -- which is
why the CONDUCTIVITY of the calibration rock decides the answer.

⚠ This is the pre-registration the packet's table states; it is NOT the binary.
  The binary runs the full depth scan (K_I over the flake with a Weibull/
  deterministic flaw), so the offline numbers bound and order the ladder, they
  do not predict the run to the kelvin. The discriminating rows are h150, the
  beam pair, Hoeser invariance and Meier -- not landing inside the band.
"""
import math
import sys

T0 = 293.15
DT_C = 843.757 - T0            # 550.607 K, the P2-2d banked level

# Kant 2016 Table 4: SOD [mm] -> (h_fl [W/m2K], T_f [K]). Tabulated only, no
# interpolation. The bracket is SOD 130 spalls / SOD 150 does not.
SOD = {50: (766, 1723.15), 80: (1009, 1683.15), 110: (873, 1633.15),
       130: (747, 1523.15), 150: (633, 1373.15)}
ANCHORS = {"h500": (500, 1673.15), "h150": (150, 1673.15)}

# (label, k [W/mK], kappa [m2/s]) -- the four property rows of the packet table.
PROPS = [
    ("lambda 2.50, kappa 1.10e-6 (rev 1's offline)", 2.50, 1.10e-6,
     (5.5, 7.5)),
    ("GOTTHARD lambda 3.27, kappa 1.60e-6", 3.27, 1.60e-6, (7.0, 10.0)),
    ("Grimsel at 25 C lambda 3.40, kappa 1.67e-6", 3.40, 1.67e-6, (7.0, 10.5)),
    ("Grimsel at 550 C lambda 1.80, kappa 5.90e-7", 1.80, 5.90e-7, (4.0, 5.5)),
]
T_MAX = 200.0                  # the calibration deck's stop_time


def theta(z, t, h, k, kap):
    """(T(z,t) - T0)/(T_f - T0), Robin half-space. 0 at t <= 0."""
    if t <= 0.0:
        return 0.0
    s = math.sqrt(kap * t)
    b = h * s / k
    u = z / (2.0 * s)
    # ⚠ exp(h z/k + b^2) overflows for large b; the product with erfc(u+b)
    # tends to 0. Use the scaled form exp(x^2)erfc(x) via a stable asymptote.
    arg = h * z / k + b * b
    if arg > 650.0:
        # exp(arg)*erfc(u+b) with erfc(x) ~ exp(-x^2)/(x sqrt(pi)):
        #   = exp(h z/k + b^2 - (u+b)^2) / ((u+b) sqrt(pi))
        e = h * z / k + b * b - (u + b) ** 2
        return math.erfc(u) - math.exp(e) / ((u + b) * math.sqrt(math.pi))
    return math.erfc(u) - math.exp(arg) * math.erfc(u + b)


def drive_max(h, Tf, k, kap, d, tmax=T_MAX):
    """max over t of [T(0,t) - T(d,t)], and the time and T_s at which the
    criterion is first met (None if never)."""
    best = 0.0
    t = 0.0
    t_fire = Ts_fire = None
    while t < tmax:
        t += 0.002 if t < 2.0 else (0.05 if t < 20.0 else 0.5)
        dr = (Tf - T0) * (theta(0.0, t, h, k, kap) - theta(d, t, h, k, kap))
        if dr > best:
            best = dr
        if t_fire is None and dr >= DT_C:
            t_fire = t
            Ts_fire = T0 + (Tf - T0) * theta(0.0, t, h, k, kap)
    return best, t_fire, Ts_fire


def ladder(k, kap, ds):
    """{SOD: set of d that fire}"""
    out = {}
    for sod, (h, Tf) in sorted(SOD.items()):
        out[sod] = [d for d in ds if drive_max(h, Tf, k, kap, d)[1] is not None]
    return out


def main():
    P = print
    ds = [x * 0.5e-3 for x in range(1, 41)]        # 0.5 .. 20.0 mm
    P("P2-2e (a1 a) -- the OFFLINE constraint-depth ladder")
    P("=" * 94)
    P(f"  drive = max_t [T(0,t) - T(d,t)];  fires at >= dT_c ="
      f" 843.757 - 293.15 = {DT_C:.3f} K;  t <= {T_MAX:.0f} s")
    P(f"  Kant 2016 Table 4 standoffs (h [W/m2K], T_f [K]):"
      f" {', '.join(f'{s}:({h},{Tf})' for s,(h,Tf) in sorted(SOD.items()))}")
    P(f"  bracket: SOD 130 SPALLS, SOD 150 does NOT (Gotthard blocks)")
    P("")
    ok = True
    for label, k, kap, want in PROPS:
        L = ladder(k, kap, ds)
        f130 = L[130]
        n150 = [d for d in ds if d not in L[150]]
        d_lo = min(f130) * 1e3 if f130 else None
        d_hi = max(n150) * 1e3 if n150 else None
        mid = (d_lo + d_hi) / 2.0 if (d_lo and d_hi) else None
        P(f"  {label}")
        P(f"    SOD  fires for d >= [mm]        (none = never fires up to"
          f" {ds[-1]*1e3:.1f} mm)")
        for sod in sorted(SOD):
            fs = L[sod]
            lab = ("%.1f" % (min(fs) * 1e3)) if fs else "never"
            P(f"    {sod:4d}  {lab:>8s}   (fires at {len(fs)} of {len(ds)} rungs)")
        got = (d_lo, d_hi)
        agree = (got == want)
        ok = ok and agree
        P(f"    separating interval  {d_lo} - {d_hi} mm   midpoint"
          f" {mid}   packet says {want[0]} - {want[1]}, midpoint"
          f" {(want[0]+want[1])/2}   {'ok' if agree else 'MISMATCH'}")
        P("")
    P("  THE MODE SWITCH (the packet's warning that small d is NOT nearly today)")
    k, kap = PROPS[1][1], PROPS[1][2]
    P(f"    on Gotthard properties, the smallest d at which each standoff fires:")
    for sod in sorted(SOD):
        h, Tf = SOD[sod]
        fs = [d for d in ds if drive_max(h, Tf, k, kap, d)[1] is not None]
        lab = ("%.1f mm" % (min(fs) * 1e3)) if fs else "never"
        P(f"      SOD {sod:3d}: {lab:>9s}")
    for d_mm in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0):
        anyf = [s for s in sorted(SOD)
                if drive_max(SOD[s][0], SOD[s][1], k, kap, d_mm * 1e-3)[1]
                is not None]
        P(f"      d = {d_mm:4.1f} mm: fires at {anyf if anyf else 'NO standoff'}")
    P("")
    P("  THE OLD ROBIN ANCHORS at the Gotthard midpoint d = 8.5 mm")
    for nm, (h, Tf) in ANCHORS.items():
        for label, k2, kap2, _ in PROPS:
            best, tf, Ts = drive_max(h, Tf, k2, kap2, 8.5e-3)
            if tf is not None:
                what = "fires t=%.2f s at T_s=%.1f K" % (tf, Ts)
            else:
                what = ("NEVER (max drive %.1f K = %.0f %% of dT_c)"
                        % (best, 100.0 * best / DT_C))
            P(f"    {nm:5s} {label[:34]:34s} -> {what}")
    P("")
    P("  THE MINIMUM FLUX the law implies (steady state: q_min = k dT_c / d)")
    for label, k2, _, _ in PROPS:
        P(f"    {label[:34]:34s}  d 8.5 mm -> q_min ="
          f" {k2*DT_C/8.5e-3/1e6:.3f} MW/m2"
          f"   (Kant measured 0.3-0.5)")
    P("=" * 94)
    P("PASS (a1 a): all four property rows reproduce" if ok else
      "FAIL (a1 a): a property row does not reproduce")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
