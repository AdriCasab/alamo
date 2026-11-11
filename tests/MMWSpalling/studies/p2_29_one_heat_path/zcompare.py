#!/usr/bin/env python3
"""P2-29 stage 1 -- THE ZHANG LEGS, compared (S/legs/<leg>/run, written by chain.sh 1).

  a3  ZF (the quarter deck's FULL-block twin, production kernel, bin_1) vs ZO (the OLD deck on bin_pre, solver=explicit,
      timestep 0.05): Temp at every cell of every plotfile, max |dT| <= 1e-6 K (prediction ~1e-9); the beam thermo
      columns. ⚠ thermo.dat prints 6 significant digits (Integrator.cpp, the stream default), so the planner's
      "rel <= 1e-9" is applied to the PRINTED tokens: equal, or one unit apart in the 6th digit (a rounding flip of a
      sub-print difference) -- counted and reported. The FIRST difference above 1e-9 K is written whatever the bar.
  a4  ZQ (the quarter block) vs ZF's quadrant, cell (i, j, k) <-> (64 + i, 64 + j, k): max |dT| <= 1e-6 K; ZF's own
      mirror asymmetry (cell 63 - i vs 64 + i, and in y) as the noise floor; thermo beam totals x 4 vs ZF within the
      print's half-units (4 x ZQ's + ZF's); the central-column bc_* tokens equal; beam_P_beam = 800 on every row.
  a6  kappa witness: kappa_phase vs Zhang's kappa(T) (Zhang.H:73-77, a CONDUCTIVITY) at 20 cells spanning the T range,
      rel <= 1e-12 -- at the t = 0 plotfile, and at the LAST plotfile against the T the refresh read: the kernel
      refreshes kappa from the step's START temperature, so the last plotfile's kappa is checked against ZQ1's last
      Temp (ZQ1 = the same run one step shorter, bit-identical up to its end). The dt message's alpha vs
      kappa/(rho c_p) at its T_cold.
  m3  ZQ vs the implicit record (leg ZI: the old deck as written, implicit, 2 s, on bin_pre): c1 surface-centre T
      at 25/50/100/150 s, c5 onset (first finite bc_z_melt) -- the zhang test's own reductions.
  m8  ZL (ZQ on the face list, evaluate-only) vs ZQ: surface-centre T at 25/50/100/150 s, the bc_z_melt onset, the
      total losses int ledger_P_loss dt (thermo on every step).
Usage: zcompare.py a3|a4|a6|m3|m8   (exit 1 on a failed gate; m3 / m8 report)
"""
import glob
import importlib.machinery
import math
import os
import re
import sys

import numpy as np
import yt

yt.set_log_level(50)
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
L = lambda leg: os.path.join(HERE, "legs", leg, "run")
ZT = importlib.machinery.SourceFileLoader(
    "zt", os.path.join(ROOT, "tests/MMWSpalling/validation/zhang/zhang_oglesby/test")).load_module()
# the implicit record: the old deck AS WRITTEN (solver = implicit, timestep 2 s) re-run on bin_pre as leg ZI -- the
# deck's own output/ holds no plotfile and no thermo.dat (only metadata), so there is no on-disk record to read
RECORD = os.path.join(HERE, "legs", "ZI", "run")
BEAM = ["beam_P_inc", "beam_P_dep", "beam_P_trans", "beam_E_inc", "beam_E_dep", "beam_E_trans"]
FAIL = []


def need(ok, msg):
    print(("  PASS  " if ok else "  FAIL  ") + msg)
    if not ok:
        FAIL.append(msg)


def pfs(run):
    return sorted(glob.glob(os.path.join(run, "*cell")))


def grid(pf, field="Temp"):
    ds = yt.load(pf)
    cg = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    return float(ds.current_time), np.array(cg[field].to_value(), dtype=np.float64)


def thermo(run):
    T = open(os.path.join(run, "thermo.dat")).read().splitlines()
    h = [c.strip() for c in T[0].split("\t")]
    return h, [[c.strip() for c in r.split("\t")] for r in T[1:] if r.strip()]


def half_unit(tok):
    """half a unit in the last printed significant digit of a %g token"""
    if tok in ("nan", "-nan", "inf", "-inf") or float(tok) == 0.0:
        return 0.0
    m = re.match(r"-?(\d+)(?:\.(\d*))?(?:e([+-]\d+))?$", tok)
    frac = len(m.group(2) or "")
    exp = int(m.group(3) or 0)
    return 0.5 * 10.0 ** (exp - frac)


def tok_eq(a, b):
    return a == b or (a.lower().endswith("nan") and b.lower().endswith("nan"))


def first_diff(times, diffs):
    for t, (mx, idx, va, vb) in zip(times, diffs):
        if mx > 1e-9:
            return "first difference > 1e-9 K: t = %g s, cell %s, Temp %.12f vs %.12f K (|d| %.3e)" % (
                t, idx, va, vb, mx)
    return "no difference above 1e-9 K at any plotfile"


def a3():
    print("a3 -- ZF (production kernel, full block, bin_1) vs ZO (the old explicit material kernel, bin_pre)")
    A, B = pfs(L("ZF")), pfs(L("ZO"))
    na, nb = [os.path.basename(p) for p in A], [os.path.basename(p) for p in B]
    need(na == nb and len(na) >= 2, "the same %d plotfiles on both legs (%s .. %s)" % (len(na), na[0], na[-1]))
    worst, times, diffs = 0.0, [], []
    for pa, pb in zip(A, B):
        ta, Ta = grid(pa)
        tb, Tb = grid(pb)
        d = np.abs(Ta - Tb)
        k = np.unravel_index(int(np.argmax(d)), d.shape)
        times.append(ta)
        diffs.append((float(d.max()), tuple(int(x) for x in k), float(Ta[k]), float(Tb[k])))
        worst = max(worst, float(d.max()))
        print("    t %6.1f s  max|dT| %.3e K at %s  (T %.3f K)" % (ta, d.max(), tuple(int(x) for x in k), Ta[k]))
    print("  " + first_diff(times, diffs))
    need(worst <= 1e-6, "Temp: max |ZF - ZO| over every cell of every plotfile %.3e K (bar 1e-6; prediction ~1e-9)" % worst)
    ha, ra = thermo(L("ZF"))
    hb, rb = thermo(L("ZO"))
    need(len(ra) == len(rb), "thermo rows %d vs %d" % (len(ra), len(rb)))
    eq = flip = bad = 0
    worst_rel, where = 0.0, ""
    for c in BEAM:
        ia, ib = ha.index(c), hb.index(c)
        for x, y in zip(ra, rb):
            a, b = x[ia], y[ib]
            if tok_eq(a, b):
                eq += 1
                continue
            r = abs(float(a) - float(b)) / max(abs(float(a)), abs(float(b)))
            if abs(float(a) - float(b)) <= half_unit(a) + half_unit(b):
                flip += 1
            else:
                bad += 1
            if r > worst_rel:
                worst_rel, where = r, "%s at t %s: %s vs %s" % (c, x[0], a, b)
    need(bad == 0, "thermo beam columns (6 printed digits): %d tokens equal, %d one-unit rounding flips, %d beyond "
         "print resolution; worst printed rel %.2e%s" % (eq, flip, bad, worst_rel, (" (" + where + ")") if where else ""))


def a4():
    print("a4 -- ZQ (the quarter block, axis on its corner) vs ZF's quadrant")
    Q, F = pfs(L("ZQ")), pfs(L("ZF"))
    nq, nf = [os.path.basename(p) for p in Q], [os.path.basename(p) for p in F]
    need(nq == nf, "the same %d plotfiles on both legs" % len(nq))
    worst, asym, times, diffs = 0.0, 0.0, [], []
    for pq, pf in zip(Q, F):
        tq, Tq = grid(pq)
        tf, Tf = grid(pf)
        nx, ny, nz = Tq.shape
        Fq = Tf[nx:2 * nx, ny:2 * ny, :]
        d = np.abs(Tq - Fq)
        k = np.unravel_index(int(np.argmax(d)), d.shape)
        mx = np.abs(Tf[nx - 1::-1, :, :] - Tf[nx:2 * nx, :, :]).max()          # x mirror: 63 - i vs 64 + i
        my = np.abs(Tf[:, ny - 1::-1, :] - Tf[:, ny:2 * ny, :]).max()          # y mirror
        times.append(tq)
        diffs.append((float(d.max()), tuple(int(x) for x in k), float(Tq[k]), float(Fq[k])))
        worst, asym = max(worst, float(d.max())), max(asym, float(mx), float(my))
        print("    t %6.1f s  max|ZQ - ZF quadrant| %.3e K at %s;  ZF mirror asymmetry x %.3e, y %.3e K" % (
            tq, d.max(), tuple(int(x) for x in k), mx, my))
    print("  " + first_diff(times, diffs))
    need(worst <= 1e-6, "Temp: max |ZQ - ZF's quadrant| %.3e K (bar 1e-6; ZF's own mirror asymmetry %.3e K, the noise "
         "floor)" % (worst, asym))
    hq, rq = thermo(L("ZQ"))
    hf, rf = thermo(L("ZF"))
    need(len(rq) == len(rf), "thermo rows %d vs %d" % (len(rq), len(rf)))
    bad, worst_rel, where = 0, 0.0, ""
    for c in BEAM:
        iq, i_f = hq.index(c), hf.index(c)
        for x, y in zip(rq, rf):
            a, b = 4.0 * float(x[iq]), float(y[i_f])
            tol = 4.0 * half_unit(x[iq]) + half_unit(y[i_f])
            if abs(a - b) > tol:
                bad += 1
            r = abs(a - b) / max(abs(a), abs(b), 1e-300)
            if r > worst_rel:
                worst_rel, where = r, "%s at t %s: 4 x %s vs %s" % (c, x[0], x[iq], y[i_f])
    need(bad == 0, "beam totals x 4 vs the full block within the print's half-units on every row (%d outside); worst "
         "printed rel %.2e (%s)" % (bad, worst_rel, where))
    bcs = [c for c in hq if c.startswith("bc_")]
    neq = [(c, x[0], x[hq.index(c)], y[hf.index(c)]) for c in bcs for x, y in zip(rq, rf)
           if not tok_eq(x[hq.index(c)], y[hf.index(c)])]
    need(not neq, "central column: %d bc_* columns, every printed token equal on every row%s" % (
        len(bcs), ("; differ: %s" % neq[:3]) if neq else ""))
    pb = {x[hq.index("beam_P_beam")] for x in rq[1:]}
    need(pb == {"800"}, "beam_P_beam on every step row: %s (the schedule's full 800 W, not hand-quartered)" % sorted(pb))
    i0 = hq.index("beam_P_inc")
    print("    first step: ZQ beam_P_inc %s W x 4 = %.6g vs ZF %s W" % (rq[1][i0], 4 * float(rq[1][i0]),
                                                                       rf[1][hf.index("beam_P_inc")]))


def zk(T):
    """Zhang.H:73-77 defaults -- a CONDUCTIVITY [W/(m K)], capped at 1273 K"""
    Tc = np.minimum(T, 1273.0)
    return 1.0 / (0.3154 + 3.8e-4 * (Tc - 273.15))


def zrho(T):
    return np.where(T >= 1273.0, 2300.0, -0.2877 * (T - 273.15) + 2884.0)


def zcp(T):
    th = T - 273.15
    return np.where(T <= 1510.5, -2.73e-4 * th * th + 0.894 * th + 788.0, 0.184 * th + 1249.0)


def a6():
    print("a6 -- the zhang law SELECTED (kappa_phase = Zhang's kappa(T), a conductivity)")
    Q, Q1 = pfs(L("ZQ")), pfs(L("ZQ1"))
    t0, T0 = grid(Q[0])
    _, K0 = grid(Q[0], "kappa_phase")
    r0 = float(np.max(np.abs(K0 - zk(T0)) / zk(T0)))
    need(t0 == 0.0 and r0 <= 1e-12, "t = 0: kappa_phase vs kappa(Temp) over every cell, max rel %.2e (bar 1e-12); "
         "kappa %.6f W/(m K) at %.2f K" % (r0, float(K0.mean()), float(T0.mean())))
    tl, Kl = grid(Q[-1], "kappa_phase")
    ts, Ts = grid(Q1[-1])
    dt = tl - ts
    print("    last plotfile %s (t %g s) vs ZQ1's last %s (t %g s): one step apart (%g s)" % (
        os.path.basename(Q[-1]), tl, os.path.basename(Q1[-1]), ts, dt))
    need(abs(dt - 0.05) < 1e-9, "ZQ1 ends exactly one step (0.05 s) before ZQ")
    flat = np.argsort(Ts, axis=None)
    pick = [flat[int(q)] for q in np.linspace(0, flat.size - 1, 20)]
    rows, worst = [], 0.0
    for p in pick:
        idx = np.unravel_index(p, Ts.shape)
        want = float(zk(np.array(Ts[idx])))
        got = float(Kl[idx])
        worst = max(worst, abs(got - want) / want)
        rows.append((tuple(int(x) for x in idx), float(Ts[idx]), got, want))
    for r in rows[::4] + [rows[-1]]:
        print("    cell %-14s T %9.3f K  kappa_phase %.12f  Zhang %.12f" % r)
    need(worst <= 1e-12, "last plotfile: kappa_phase at 20 cells spanning T %.1f .. %.1f K vs Zhang's kappa at the "
         "step's start T (ZQ1), max rel %.2e (bar 1e-12); kappa spans %.4f .. %.4f W/(m K) (non-vacuous)" % (
             rows[0][1], rows[-1][1], worst, min(r[2] for r in rows), max(r[2] for r in rows)))
    log = open(os.path.join(HERE, "legs", "ZQ", "run.log")).read()
    m = re.search(r"alpha (\S+) m\^2/s at T_cold (\S+) K", log)
    a_pr, Tc = float(m.group(1)), float(m.group(2))
    a_want = float(zk(np.array(Tc)) / (zrho(np.array(Tc)) * zcp(np.array(Tc))))
    need(abs(a_pr / a_want - 1) <= 5e-6, "the dt message's alpha %s m^2/s at T_cold %s K = kappa/(rho c_p) %.6g (a "
         "DIFFUSIVITY; at 297.5 K it would be %.5g -- the packet's number; T_cold = min(T_ref 273.15, T_amb))" % (
             m.group(1), m.group(2), a_want,
             float(zk(np.array(297.5)) / (zrho(np.array(297.5)) * zcp(np.array(297.5))))))


def surf_T(run):
    ax = ZT.axial(run, ZT.meta(run))
    return [(t, Tc[0]) for t, z, Tc in ax]


def onset_t(run):
    h, rows = thermo(run)
    i = h.index("bc_z_melt")
    for r in rows[1:]:
        if math.isfinite(float(r[i])):
            return float(r[0]), float(r[i])
    return None, None


def m3():
    print("m3 (report) -- ZQ (explicit, dt 0.05 s, production kernel) vs the implicit record (dt 2 s, %s)" %
          os.path.relpath(RECORD, ROOT))
    a, b = surf_T(L("ZQ")), surf_T(RECORD)
    for t in (25.0, 50.0, 100.0, 150.0):
        ta, tb = ZT.at(a, t), ZT.at(b, t)
        print("    c1 surface-centre T at %5.1f s: ZQ %8.2f K, implicit %8.2f K  (%+.2f %%)" % (t, ta, tb, 100 * (ta / tb - 1)))
    (oa, za), (ob, zb) = onset_t(L("ZQ")), onset_t(RECORD)
    print("    c5 onset (first finite bc_z_melt): ZQ %s, implicit %s%s" % (
        "%.2f s (z %.4g m)" % (oa, za) if oa else "never", "%.2f s (z %.4g m)" % (ob, zb) if ob else "never",
        ("  (%+.1f %%; prediction 5-10 %% earlier)" % (100 * (oa / ob - 1))) if oa and ob else ""))
    print("    prediction: T within 2 %, onset 5-10 % earlier (P2-20..22's explicit-vs-implicit record)")


def m8():
    print("m8 (report) -- ZL (the face list, evaluate-only: losses on the surface T_s) vs ZQ (the cell form, top-cell T)")
    a, b = surf_T(L("ZL")), surf_T(L("ZQ"))
    for t in (25.0, 50.0, 100.0, 150.0):
        ta, tb = ZT.at(a, t), ZT.at(b, t)
        print("    surface-centre (top-cell) T at %5.1f s: ZL %8.2f K, ZQ %8.2f K  (%+.3f %%)" % (t, ta, tb,
                                                                                         100 * (ta / tb - 1)))
    (oa, za), (ob, zb) = onset_t(L("ZL")), onset_t(L("ZQ"))
    print("    bc_z_melt onset: ZL %s, ZQ %s%s" % ("%.2f s" % oa if oa else "never", "%.2f s" % ob if ob else "never",
                                                    ("  (%+.1f %%)" % (100 * (oa / ob - 1))) if oa and ob else ""))
    for leg in ("ZL", "ZQ"):
        h, rows = thermo(L(leg))
        i = h.index("ledger_P_loss")
        ts = [float(r[0]) for r in rows]
        steps = sorted(set(round(t2 - t1, 9) for t1, t2 in zip(ts, ts[1:])))
        E = sum(float(r[i]) * 0.05 for r in rows[1:])
        print("    %s total losses int ledger_P_loss dt = %.4f kJ over %d step rows (thermo dt %s); final P_loss %s W "
              "(one QUARTER block)" % (leg, E / 1e3, len(rows) - 1, steps, rows[-1][i]))
    print("    prediction: ZL cooler, |dT|/T <= 0.5 % at 50 s and <= 3 % at 150 s, onset later by <= 10 %")


if __name__ == "__main__":
    {"a3": a3, "a4": a4, "a6": a6, "m3": m3, "m8": m8}[sys.argv[1]]()
    print("RESULT: %s" % ("FAIL (%d)" % len(FAIL) if FAIL else "PASS"))
    sys.exit(1 if FAIL else 0)
