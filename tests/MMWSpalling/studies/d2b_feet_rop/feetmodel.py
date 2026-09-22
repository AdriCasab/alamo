#!/usr/bin/env python3
"""D2b hand model: burner-on-feet equilibrium with the D2a2 wall-jet closure.

Written and run BEFORE any D2b simulation (deliverable 5). Uses
d2a2_walljet/walljet.py (D = 7.5 mm; Martin h shape; cp = 1250).

Model (full jet, axisymmetric, steady descent at rate v):
  * Ring. The pad rule takes the 0.9 quantile of the annulus face heights.
    Recession falls with r, so the highest 10 % are the outermost columns;
    the ring is represented by the 0.9 area-quantile radius of the annulus,
    r_q = sqrt(ri^2 + q (ro^2 - ri^2)) = 38.97 mm, at stand-off s = foot_standoff.
  * Centre. The axis column sits at s_c, set by centre ROP(s_c) = v (the centre
    self-regulates, D2a2).
  * Stand-off between them: s = s_c for r <= r_core (the flat stagnation disk),
    linear in r from s_c at r_core to foot_standoff at r_q.
  * Stagnation: phi = min(1, 5 D / s_c), T_stag = T_mix + (T_nozzle - T_mix) phi,
    m = mdot / phi (entrained mass); T_mix = T_rec (exhaust recirculation) or
    293.15 K (fixed). Nozzle-mass variant: m = mdot.
  * March: m cp dT/dr = -h(r, s(r)) (T - T_s) 2 pi r with T_s = T_fire = 821 K
    (the firing face; losses are not taken from the gas, as in ALAMO's
    max(0, h (T_gas - T_s))), from r = 0 to r_q. Beyond r_q every column falls
    behind the burner and freezes, so in the steady state the march ends at r_q
    and T_exhaust = T_gas(r_q); T_rec = T_exhaust is iterated to its fixed
    point.
  * Ring ROP = pinned closed form rop_closed(h(r_q, foot_standoff), T_gas(r_q))
    (with face losses). The equilibrium is v = ring ROP(v).
  * Depth profile. A column at r > r_q recedes at v(r) (the march continued at
    s = foot_standoff from T_gas(r_q): its rate while it is still well inside
    the jet, an upper estimate) until the burner is foot_standoff ahead of it:
        own depth at freeze  D(r) = foot_standoff * v(r) / (v - v(r)).
    The packet's d(r) = foot_standoff * v / (v - v(r)) is the burner's (feet's)
    depth at that moment, exactly foot_standoff deeper than the rock; the
    wall is at D(r). r_wall(z) = max{r : D(r) > z}, with r_wall >= r_q.
    Depth-mean diameter = mean over z of 2 r_wall(z).
`python feetmodel.py` prints everything RESULTS.md section 0 quotes.

D2c extension (studies/d2c_steady): every model function takes an optional
cfg = Cfg(far, n, core, De_ref). The default Cfg() is the D2b model exactly
(clamped h, nozzle diameter D, core 5): h = walljet.h_py(..., far, n),
phi = min(1, core D_dec/s_c) with D_dec = D, or (De_ref given) the momentum
diameter De_ref sqrt(T_mix/T_ent).
"""
import math
import os
import sys
from dataclasses import dataclass
from typing import Optional

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "d2a2_walljet"))
import walljet as wj  # noqa: E402

D = wj.D
RI, RO, STANDOFF, QUANT = 0.028, 0.040, 0.050, 0.9
R_Q = math.sqrt(RI * RI + QUANT * (RO * RO - RI * RI))
T_ENT = wj.TA
BLOCK = 0.50


@dataclass(frozen=True)
class Cfg:
    far: str = "clamp"               # walljet far law: clamp | power
    n: float = 1.0                   # power-law exponent beyond 12 D
    core: float = 5.0                # core length in decay-diameter units
    De_ref: Optional[float] = None   # None: nozzle D; else D_e at T_ENT [m]


DEF = Cfg()


def h_of(cfg, h_ref, r, s):
    return wj.h_py(h_ref, r, s, cfg.far, cfg.n)


def d_dec(cfg, T_mix):
    return D if cfg.De_ref is None else cfg.De_ref * math.sqrt(T_mix / T_ENT)


def phi_of(cfg, s_c, T_mix):
    return min(1.0, cfg.core * d_dec(cfg, T_mix) / s_c)


def s_of_r(r, s_c):
    if r <= wj.R_CORE:
        return s_c
    if r >= R_Q:
        return STANDOFF
    return s_c + (STANDOFF - s_c) * (r - wj.R_CORE) / (R_Q - wj.R_CORE)


def s_profile(r, s_c, s_beyond=STANDOFF):
    s = np.where(r <= wj.R_CORE, s_c, s_c + (STANDOFF - s_c) * (r - wj.R_CORE) / (R_Q - wj.R_CORE))
    return np.where(r >= R_Q, s_beyond, s)


def march(h_ref, T_stag, m, s_c, r_end, n=2001, s_beyond=STANDOFF, r_start=0.0, T0=None, cfg=DEF):
    """T_gas(r) on [r_start, r_end]: the excess over T_fire decays as
    exp(-int 2 pi r h / (m cp) dr) (midpoint rule)."""
    r = np.linspace(r_start, r_end, n)
    dr = r[1] - r[0]
    rm = r[:-1] + 0.5 * dr
    h = h_of(cfg, h_ref, rm, s_profile(rm, s_c, s_beyond))
    a = 2.0 * math.pi * rm * h / (m * wj.CP) * dr
    ex0 = max(0.0, (T_stag if T0 is None else T0) - wj.T_FIRE)
    T = wj.T_FIRE + ex0 * np.exp(-np.concatenate(([0.0], np.cumsum(a))))
    return r, T


def stag(T_nozzle, s_c, T_mix, entrained=True, cfg=DEF):
    phi = phi_of(cfg, s_c, T_mix)
    T_stag = T_mix + (T_nozzle - T_mix) * phi
    m = wj.MDOT / phi if entrained else wj.MDOT
    return T_stag, m


def centre_rop_at(h_ref, T_nozzle, s_c, T_mix, entrained, cfg=DEF):
    T_stag, _ = stag(T_nozzle, s_c, T_mix, entrained, cfg)
    return float(wj.rop_closed(h_of(cfg, h_ref, 0.0, s_c), T_stag))


S_C_MAX = 2.0   # deep-pit limit: s_c at the bisection bound means "no centre equilibrium"


def s_c_for(h_ref, T_nozzle, v, T_mix, entrained, cfg=DEF):
    """Centre stand-off with centre ROP = v (bisection on s_c in [2D, 2 m]);
    None if the centre cannot reach v even at s_c = 2D."""
    lo, hi = 2.0 * D, S_C_MAX
    if centre_rop_at(h_ref, T_nozzle, lo, T_mix, entrained, cfg) < v:
        return None
    if centre_rop_at(h_ref, T_nozzle, hi, T_mix, entrained, cfg) > v:
        return hi
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if centre_rop_at(h_ref, T_nozzle, mid, T_mix, entrained, cfg) > v:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def state(h_ref, T_nozzle, v, mode="exhaust", entrained=True, cfg=DEF):
    """Fixed point in T_rec for a given v. Returns dict or None."""
    T_rec = T_ENT
    for _ in range(400):
        T_mix = T_rec if mode == "exhaust" else T_ENT
        s_c = s_c_for(h_ref, T_nozzle, v, T_mix, entrained, cfg)
        if s_c is None:
            return None
        T_stag, m = stag(T_nozzle, s_c, T_mix, entrained, cfg)
        r, T = march(h_ref, T_stag, m, s_c, R_Q, cfg=cfg)
        T_new = T[-1]
        if mode != "exhaust" or abs(T_new - T_rec) < 1e-6:
            T_rec = T_new
            break
        T_rec = 0.5 * T_rec + 0.5 * T_new
    ring = float(wj.rop_closed(h_of(cfg, h_ref, R_Q, STANDOFF), T[-1]))
    P_face = m * wj.CP * (T_stag - T[-1])
    T_mix = T_rec if mode == "exhaust" else T_ENT
    return dict(v=v, s_c=s_c, T_stag=T_stag, T_rec=T_rec if mode == "exhaust" else T_ENT,
                T_ring=T[-1], ring=ring, m=m, P_face=P_face, phi=phi_of(cfg, s_c, T_mix),
                D_e=d_dec(cfg, T_mix))


def equilibria(h_ref, T_nozzle, mode="exhaust", entrained=True, vmin=0.05, vmax=30.0, n=80, cfg=DEF):
    """Roots of g(v) = ring(v) - v on a log grid, refined by bisection.
    Stability: d ring/dv < 1 at the root (a faster burner slows the ring)."""
    vs = np.geomspace(vmin, vmax, n)
    g, st = [], []
    for v in vs:
        s = state(h_ref, T_nozzle, v, mode, entrained, cfg)
        st.append(s)
        g.append(np.nan if s is None else s["ring"] - v)
    g = np.array(g)
    roots = []
    for i in range(n - 1):
        if np.isfinite(g[i]) and np.isfinite(g[i + 1]) and g[i] * g[i + 1] < 0.0:
            lo, hi, glo = vs[i], vs[i + 1], g[i]
            for _ in range(30):
                mid = math.sqrt(lo * hi)
                sm = state(h_ref, T_nozzle, mid, mode, entrained, cfg)
                gm = sm["ring"] - mid
                if gm * glo > 0.0:
                    lo, glo = mid, gm
                else:
                    hi = mid
            s = state(h_ref, T_nozzle, math.sqrt(lo * hi), mode, entrained, cfg)
            s["stable"] = bool(g[i] > 0.0 > g[i + 1])
            roots.append(s)
    return roots, vs, g


def depth_profile(h_ref, eq, r_max=0.20, n=2001, z_max=BLOCK, cfg=DEF):
    """v(r) beyond r_q, freeze depth D(r) and r_wall(z)."""
    r, T = march(h_ref, eq["T_stag"], eq["m"], eq["s_c"], r_max, n=n, r_start=R_Q, T0=eq["T_ring"],
                 cfg=cfg)
    vr = wj.rop_closed(h_of(cfg, h_ref, r, STANDOFF), T)
    v = eq["v"]
    with np.errstate(divide="ignore"):
        Dz = np.where(vr < v, STANDOFF * vr / np.maximum(v - vr, 1e-30), np.inf)
    d_pkt = np.where(vr < v, STANDOFF * v / np.maximum(v - vr, 1e-30), np.inf)
    z = np.linspace(0.0, z_max, 501)
    rw = np.array([max(R_Q, float(r[Dz > zz].max())) if np.any(Dz > zz) else R_Q for zz in z])
    return dict(r=r, vr=vr, D=Dz, d_pkt=d_pkt, z=z, r_wall=rw)


def depth_mean(prof, z_lo, z_hi):
    sel = (prof["z"] >= z_lo) & (prof["z"] <= z_hi)
    return float(np.mean(2.0 * prof["r_wall"][sel]))


CASES = {
    "J-M 1900 exhaust (scored)": ("M", 1900.0, "exhaust", True),
    "J-5 1900 exhaust": ("5", 1900.0, "exhaust", True),
    "J-10 1900 exhaust": ("10", 1900.0, "exhaust", True),
    "J-M 1900 fixed 293 K": ("M", 1900.0, "fixed", True),
    "J-M 1900 fixed, nozzle mass": ("M", 1900.0, "fixed", False),
    "J-M 1436 exhaust": ("M", 1436.0, "exhaust", True),
}
T_STAGE_B = 600.0


def report():
    out = []
    out.append(f"r_q = {R_Q * 1e3:.2f} mm (0.9 area quantile of [28, 40] mm); r_core = {wj.R_CORE * 1e3:.2f} mm; "
               f"foot_standoff = {STANDOFF * 1e3:.0f} mm; D = {D * 1e3} mm; cp = {wj.CP:g}")
    rows = []
    for name, (hk, Tn, mode, ent) in CASES.items():
        h = float(wj.h_anchor(hk, Tn))
        roots, vs, g = equilibria(h, Tn, mode, ent)
        if not roots:
            gmax = np.nanmax(g)
            rows.append((name, h, None, None, f"no equilibrium on [0.05, 30] m/h (max ring - v = {gmax:.2f})"))
            continue
        for eq in roots:
            prof = depth_profile(h, eq)
            d_drill = min(BLOCK, eq["v"] / 3600.0 * T_STAGE_B)
            rows.append((name, h, eq, prof,
                         dict(mean_block=depth_mean(prof, 0.0, BLOCK),
                              mean_drill=depth_mean(prof, 0.0, d_drill), d_drill=d_drill,
                              min_d=2.0 * float(prof["r_wall"].min()),
                              r_edge=float(prof["r"][np.nonzero(prof["vr"] > 0)[0][-1]])
                              if np.any(prof["vr"] > 0) else R_Q)))
    return out, rows


def d2a2_check():
    """Feed D2a2's measured early-window centre stand-off (and ring stand-off)
    into the march with D2a2's treatments (T_mix = 293.15 K) and compare the
    ring ROP at r = 40 mm with D2a2's early-window values."""
    meas = {  # name: (h key, T_nozzle, entrained, s_c E mean [m], ring s at 75 s [m], D2a2 ring ROP)
        "JM_19": ("M", 1900.0, False, 0.062, None, 0.54),
        "J5_19": ("5", 1900.0, False, 0.084, None, 0.57),
        "J10_19": ("10", 1900.0, False, 0.083, None, 0.71),
        "J5_19_ent": ("5", 1900.0, True, 0.087, None, 0.62),
        "J10_19_ent": ("10", 1900.0, True, 0.089, None, 1.27),
    }
    return meas


def march_d2a2(h_ref, T_nozzle, s_c, s_ring, entrained):
    """D2a2 geometry: stand-off linear from s_c (r <= r_core) to s_ring at 40 mm,
    T_mix = 293.15 K; ring ROP at r = 40 mm."""
    phi = min(1.0, 5.0 * D / s_c)
    T_stag = T_ENT + (T_nozzle - T_ENT) * phi
    m = wj.MDOT / phi if entrained else wj.MDOT
    r = np.linspace(0.0, 0.040, 2001)
    dr = r[1] - r[0]
    rm = r[:-1] + 0.5 * dr
    s = np.where(rm <= wj.R_CORE, s_c, s_c + (s_ring - s_c) * (rm - wj.R_CORE) / (0.040 - wj.R_CORE))
    a = 2.0 * math.pi * rm * wj.h_py(h_ref, rm, s) / (m * wj.CP) * dr
    T40 = wj.T_FIRE + (T_stag - wj.T_FIRE) * math.exp(-a.sum())
    return float(wj.rop_closed(wj.h_py(h_ref, 0.040, s_ring), T40)), T_stag, T40


def main():
    head, rows = report()
    print("\n".join(head))
    for name, h, eq, prof, extra in rows:
        if eq is None:
            print(f"{name}: h_ref {h:g}: {extra} -> the ring cannot keep up with the centre-set "
                  f"stand-off at any v: STALL predicted")
            hk, Tn, mode, ent = CASES[name]
            s_stall = 5.0 * D * (Tn - T_ENT) / (wj.T_FIRE - T_ENT)
            parts = []
            for sc in (0.050, 0.062, 0.080, 0.100):
                rop, Ts, T40 = march_d2a2(h, Tn, sc, STANDOFF, ent)
                parts.append(f"s_c {sc * 1e3:.0f} mm -> ring {rop:.2f} m/h (T_stag {Ts:.0f} K)")
            print(f"   centre self-stall stand-off (T_stag = T_fire) {s_stall * 1e3:.0f} mm; "
                  f"ring ROP at r = 40 mm if the centre were held: " + "; ".join(parts))
            continue
        deep = eq["s_c"] >= S_C_MAX
        print(f"{name}: h_ref {h:g}: v = {eq['v']:.3f} m/h ({'stable' if eq['stable'] else 'UNSTABLE'}), "
              + ("NO centre equilibrium (centre outruns the ring at every s_c; deep-pit limit "
                 "s_c -> inf, evaluated at 2 m)," if deep else f"s_c = {eq['s_c'] * 1e3:.1f} mm,")
              + f" phi = {eq['phi']:.3f}, T_stag = {eq['T_stag']:.0f} K, "
              f"T_rec = {eq['T_rec']:.0f} K, T_gas(r_q) = {eq['T_ring']:.0f} K, "
              f"P_face = {eq['P_face'] / 1e3:.2f} kW; v(r) > 0 out to {extra['r_edge'] * 1e3:.0f} mm; "
              f"depth-mean Ø 0-0.5 m {extra['mean_block'] * 1e3:.1f} mm, "
              f"0-{extra['d_drill'] * 1e3:.0f} mm (600 s) {extra['mean_drill'] * 1e3:.1f} mm, "
              f"min Ø {extra['min_d'] * 1e3:.1f} mm")
        zz = [0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5]
        print("   r_wall(z): " + ", ".join(
            f"{z * 1e3:.0f} mm -> {np.interp(z, prof['z'], prof['r_wall']) * 1e3:.1f}" for z in zz))
    print("D2a2 cross-check (T_mix = 293.15 K, s linear from s_c to s_ring at 40 mm):")
    for n, (hk, Tn, ent, s_c, _, rop_meas) in d2a2_check().items():
        h = float(wj.h_anchor(hk, Tn))
        for s_ring in (0.030, 0.020):
            rop, Ts, T40 = march_d2a2(h, Tn, s_c, s_ring, ent)
            print(f"  {n}: s_c {s_c * 1e3:.0f} mm, s_ring {s_ring * 1e3:.0f} mm: T_stag {Ts:.0f} K, "
                  f"T_gas(40) {T40:.0f} K, ring ROP {rop:.2f} m/h (D2a2 measured {rop_meas:.2f})")


if __name__ == "__main__":
    main()
