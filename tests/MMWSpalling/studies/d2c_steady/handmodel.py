#!/usr/bin/env python3
"""D2c hand model: steady state (d2b_feet_rop/feetmodel with the D2c
closures) and a transient burner-on-feet model for the time to steady state,
the mouth (free-surface dilution) and the whole excavation.

Steady state: feetmodel.equilibria(cfg), unchanged D2b geometry (ring at the
0.9 area quantile r_q = 38.97 mm, stand-off linear from s_c at r_core to 50 mm
at r_q, T_s = T_fire, march ends at r_q, T_rec = T_gas(r_q)). In a deep hole
no free-surface column marches (the nozzle plane is below the original
surface once the feet are 50 mm deep), so dilution does not enter it.

Transient (explicit in time, dt = 0.25 s, axisymmetric, full jet):
  state: feet (ring) depth d, centre depth z_c, outer column depths z(r) on the
         mesh's 2 mm grid from the annulus edge 40 mm (centres 41, 43, ... mm)
         to R_MAX, plus the annulus band r_q < r <= 40 mm (one column at
         39.5 mm), T_rec; depths below the original surface, all zero at
         t = 0, T_rec(0) = T_ent.
  step:  s_c = 0.05 + z_c - d;  T_mix = T_rec;  D_e and phi from cfg;
         T_stag, m = mdot/phi;  inner march 0 -> r_q (feetmodel.march);
         outer columns in order of r, each with s = 0.05 + z - d:
           s <= 0: frozen (above the nozzle plane), gas passes unchanged;
           s > 0:  T_gas -> T_fire + (T_gas - T_fire) exp(-2 pi r h dr/(m cp)),
                   v(r) = rop_closed(h(r, s), inlet T_gas);
                   free-surface (r > 40 mm and z <= aspect (r - 40 mm)):
                   then m -> m ((r + dr/2)/(r - dr/2))^e and
                   T_gas -> T_ent + (T_gas - T_ent) m_old/m (ALAMO's per-bin rule);
         T_rec(next) = inlet T_gas of the first free-surface column (the
         mouth), else the exhaust T; v_ring = rop_closed(h(r_q, 0.05),
         T_gas(r_q)), v_c = rop_closed(h(0, s_c), T_stag);
         d += v_ring dt, z_c += v_c dt, z += v dt (marching columns);
         clear = True (foot_body_clearance): the annulus band is cut to d
         whenever it stands proud (z < d), the cut volume is logged.
  The ring stands in for the pad rule (feet on r_q, as the D2b model). Face
  temperature T_fire everywhere (cold columns would draw more). Volume =
  integral of 2 pi r z(r) dr to R_MAX (inside r_q: z = d - 0.05 + s(r)).
Wall profile: r_wall(z) = max{r : z(r) > z}, floor r_q (Ø = 2 r_wall).
"""
import math
import os
import sys
from dataclasses import dataclass

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "d2b_feet_rop"))
sys.path.insert(0, HERE)
import feetmodel as fm  # noqa: E402
import nozzle as nz  # noqa: E402

wj = fm.wj
T_ENT = fm.T_ENT
R_MAX = 0.12            # quarter domain half-width (inscribed circle)
DR = 2.0e-3            # the mesh
DT = 0.25
STEADY_DSC = 0.02e-3 / 1.0        # m/s  (0.02 mm/s)
STEADY_LEN = 150.0                # s
SUB = 50.0                        # s, sub-fit length
STEADY_TOL = 0.10


@dataclass(frozen=True)
class FS:
    on: bool = False
    aspect: float = 1.0
    exponent: float = 1.0


def cfg_for(T_nozzle, far="power", n=1.0, core=8.0, momentum=True):
    return fm.Cfg(far, n, core, nz.de_ref(T_nozzle) if momentum else None)


def h_ref_for(T_nozzle):
    """Martin anchor at T_nozzle (the D2b J-M rule, walljet.h_anchor)."""
    return float(wj.h_anchor("M", T_nozzle))


def steady(T_nozzle, cfg, h_ref=None):
    h = h_ref_for(T_nozzle) if h_ref is None else h_ref
    roots, _, _ = fm.equilibria(h, T_nozzle, "exhaust", True, cfg=cfg)
    stable = [e for e in roots if e.get("stable")]
    return (stable or roots or [None])[0]


def transient(T_nozzle, cfg, fs=FS(), h_ref=None, t_end=900.0, dt=DT, z_stop=None, n_inner=401,
              clear=True):
    """Returns a dict of time series and the outer grid state at the end.
    z_stop: stop when the centre depth reaches it (domain bottom rule)."""
    h_ref = h_ref_for(T_nozzle) if h_ref is None else h_ref
    RQ, RO = fm.R_Q, fm.RO
    r_o = np.concatenate(([0.5 * (RQ + RO)], RO + DR * (np.arange(int(round((R_MAX - RO) / DR))) + 0.5)))
    w_o = np.concatenate(([RO - RQ], np.full(r_o.size - 1, DR)))
    band = r_o <= RO
    h_ring = float(fm.h_of(cfg, h_ref, RQ, fm.STANDOFF))
    z = np.zeros_like(r_o)
    d = z_c = 0.0
    T_rec = T_ENT
    out = {k: [] for k in ("t", "d", "z_c", "s_c", "T_rec", "T_stag", "T_ring", "T_exh", "v_ring",
                           "v_c", "phi", "D_e", "m_ratio", "n_fs", "V", "P_face", "V_mech")}
    V_mech = 0.0
    t = 0.0
    nsteps = int(round(t_end / dt))
    rc = np.linspace(0.0, RQ, n_inner)
    for _ in range(nsteps):
        s_c = max(1.0e-4, fm.STANDOFF + z_c - d)
        T_stag, m = fm.stag(T_nozzle, s_c, T_rec, True, cfg)
        phi = fm.phi_of(cfg, s_c, T_rec)
        _, Tin = fm.march(h_ref, T_stag, m, s_c, RQ, n=n_inner, cfg=cfg)
        T = float(Tin[-1])
        T_ring = T
        v_ring = float(wj.rop_closed(h_ring, T_ring)) / 3600.0
        v_c = float(wj.rop_closed(fm.h_of(cfg, h_ref, 0.0, s_c), T_stag)) / 3600.0
        s = fm.STANDOFF + z - d
        v = np.zeros_like(z)
        m_cur = m
        T_mouth = None
        nfs = 0
        for i in range(r_o.size):
            if s[i] <= 0.0:
                continue
            r = r_o[i]
            T_in = T
            hi = float(fm.h_of(cfg, h_ref, r, s[i]))
            v[i] = float(wj.rop_closed(hi, T)) / 3600.0
            if T > wj.T_FIRE:
                T = wj.T_FIRE + (T - wj.T_FIRE) * math.exp(-2.0 * math.pi * r * hi * w_o[i] / (m_cur * wj.CP))
            if fs.on and r > RO and z[i] <= fs.aspect * (r - RO):
                if T_mouth is None:
                    T_mouth = T_in
                nfs += 1
                m_new = m_cur * ((r + 0.5 * w_o[i]) / (r - 0.5 * w_o[i])) ** fs.exponent
                T = T_ENT + (T - T_ENT) * m_cur / m_new
                m_cur = m_new
        T_exh = T
        T_rec_next = T_mouth if T_mouth is not None else T_exh
        # volume
        rin = rc
        zin = d - fm.STANDOFF + fm.s_profile(rin, s_c)
        V = float(np.trapezoid(2.0 * math.pi * rin * zin, rin)) + float(np.sum(2.0 * math.pi * r_o * z * w_o))
        for k, val in (("t", t), ("d", d), ("z_c", z_c), ("s_c", s_c), ("T_rec", T_rec),
                       ("T_stag", T_stag), ("T_ring", T_ring), ("T_exh", T_exh), ("v_ring", v_ring),
                       ("v_c", v_c), ("phi", phi), ("D_e", fm.d_dec(cfg, T_rec)),
                       ("m_ratio", m_cur / m), ("n_fs", nfs), ("V", V),
                       ("P_face", m * wj.CP * (T_stag - T_ring)), ("V_mech", V_mech)):
            out[k].append(val)
        d += v_ring * dt
        z_c += v_c * dt
        z = z + v * dt
        if clear:
            cut = band & (z < d)
            V_mech += float(np.sum(2.0 * math.pi * r_o[cut] * w_o[cut] * (d - z[cut])))
            z[cut] = d
        T_rec = T_rec_next
        t += dt
        if z_stop is not None and z_c >= z_stop:
            break
    res = {k: np.array(v) for k, v in out.items()}
    res.update(r_o=r_o, z=z, h_ref=h_ref, cfg=cfg, fs=fs, T_nozzle=T_nozzle)
    return res


def steady_window(res):
    """Earliest [t0, t0 + 150 s] with |ds_c/dt| < 0.02 mm/s throughout and the
    50 s sub-fits of d(t) within ±10 % of the window fit (the score's rule on
    the smooth hand series). Returns (t0, rate m/s) or (None, None)."""
    t, d, sc = res["t"], res["d"], res["s_c"]
    if t.size < 3:
        return None, None
    dsc = np.gradient(sc, t)
    step = t[1] - t[0]
    nwin = int(round(STEADY_LEN / step))
    nsub = int(round(SUB / step))
    ok_sc = np.abs(dsc) < STEADY_DSC
    for i0 in range(0, t.size - nwin, int(round(5.0 / step))):
        sl = slice(i0, i0 + nwin + 1)
        if not np.all(ok_sc[sl]):
            continue
        rate = np.polyfit(t[sl], d[sl], 1)[0]
        good = True
        for j in range(i0, i0 + nwin - nsub + 1, nsub):
            sub = slice(j, j + nsub + 1)
            rs = np.polyfit(t[sub], d[sub], 1)[0]
            if abs(rs / rate - 1.0) > STEADY_TOL:
                good = False
                break
        if good:
            return float(t[i0]), float(rate)
    return None, None


def wall_profile(res, z_grid):
    r_o, z = res["r_o"], res["z"]
    return np.array([max(fm.R_Q, float(r_o[z > zz].max())) if np.any(z > zz) else fm.R_Q
                     for zz in z_grid])


def profile_metrics(res, drilled):
    zg = np.linspace(0.0, max(drilled, 1e-3), 401)
    rw = wall_profile(res, zg)
    return dict(d20=2.0 * float(wall_profile(res, [0.020])[0]),
                d50=2.0 * float(wall_profile(res, [0.050])[0]),
                mean=float(np.mean(2.0 * rw)), min=float(np.min(2.0 * rw)), zg=zg, rw=rw)


def volume_rate(res, t0, t1):
    sel = (res["t"] >= t0) & (res["t"] <= t1)
    return float(np.polyfit(res["t"][sel], res["V"][sel], 1)[0])
