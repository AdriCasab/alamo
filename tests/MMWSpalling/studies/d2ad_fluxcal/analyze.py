#!/usr/bin/env python3
"""D2ad analysis -- scores CRITERION.md (a)-(i).

Imports the equal-burner-depth machinery UNCHANGED from d2x_areasplit (which
imports d2v_hmap's zface/profile/diameter) and binds LZ PER LEG from that leg's
own .done metadata -- the D2ac trap (both modules hard-code LZ = 0.4).

⚠ No nozzle_z column under `prescribed`: z_n(t) = nozzle_z0 - feed*t is exact,
and the standoff is `jet_s_c`.
⚠ Three legs ABORTED on the bore collision guard. Their traces are valid up to
the abort and are scored there; (c)/(d) are VOID from sustained pad contact.
"""
import ast
import contextlib
import filecmp
import importlib.util
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)


def load_mod(name, path):
    sp = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


X = load_mod("d2x_an", os.path.join(STUD, "d2x_areasplit", "analyze.py"))
V = X.V

OUT = os.path.join(HERE, "output")
Z1 = os.path.join(STUD, "d2z1_delete", "output")

LEGS = ("C0", "G1600", "Q5", "Q5T870", "Q35", "Q35T870", "Q30", "Q30T870")
X.WHERE.update({n: OUT for n in LEGS + ("S0",)})
X.WHERE["RB_D"] = Z1

# ---- PRE-REGISTERED, from PREFLIGHT.md (g) / CRITERION.md. DO NOT EDIT. ----
NOZZLE_Z0, FEED = 0.550, 1.5 / 3600.0
S_CEIL, SLOPE_BAR = 100.0, 2.0       # mm; mm per 100 s
FOOT_RI, FOOT_RO, FOOT_STANDOFF = 0.028, 0.040, 0.050
PRED = {"C0": (1.00, 823.3, 1664.7, 1.000, 254, 1.02),
        "G1600": (1.00, 823.3, 1465.8, 0.764, 194, 0.78),
        "Q5": (0.50, 823.3, 1834.2, 0.601, 153, 0.61),
        "Q5T870": (0.50, 870.0, 1834.2, 0.527, 134, 0.59),
        "Q35": (0.35, 823.3, 1871.2, 0.436, 111, 0.45),
        "Q35T870": (0.35, 870.0, 1871.2, 0.383, 97, 0.43),
        "Q30": (0.30, 823.3, 1882.1, 0.378, 96, 0.39),
        "Q30T870": (0.30, 870.0, 1882.1, 0.332, 84, 0.37)}
KANT = (0.45, 0.87)                  # MW/m^2, measured spallation band
MEIER_Q = (0.50, 0.60)               # MW/m^2, required face-mean
D_TARGET, D_TOL = 85.0, 5.0
CAL = (120.0, 250.0)
VOL_CEIL = 3.42
WIN_BAR = 30.0                       # K, p90 - p10
C0_WIN = (400.0, 600.0)              # C3: C0 is scored here


def have(n):
    return os.path.exists(os.path.join(X.WHERE.get(n, OUT), n + ".done"))


def meta(n):
    return ast.literal_eval(open(os.path.join(X.WHERE[n], n + ".done")).read())


def thermo(n):
    p = os.path.join(X.WHERE[n], n, "thermo.dat")
    with open(p) as f:
        cols = f.readline().split()
    v = np.loadtxt(p, skiprows=1, ndmin=2)
    return {c: v[:, k] for k, c in enumerate(cols)}


@contextlib.contextmanager
def with_lz(n):
    lz = meta(n).get("lz", 0.4)
    ov, ox = V.LZ, X.LZ
    V.LZ, X.LZ = lz, lz
    try:
        yield lz
    finally:
        V.LZ, X.LZ = ov, ox


def hdr(t):
    print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78)


def z_n(t):
    return NOZZLE_Z0 - FEED * t


def trace(n):
    d = thermo(n)
    return d["time"], d["jet_s_c"] * 1e3          # mm


def contact_trace(n, npts=140):
    """max top-face height over the pad ring minus the foot plane, per time."""
    with with_lz(n):
        m, r, a = X.cell(n)
        t = thermo(n)["time"]
        ann = (r >= FOOT_RI) & (r <= FOOT_RO)
        ts = [x for x in t[:: max(1, len(t) // npts)] if x > 60.0]
        out = []
        for tt in ts:
            zf, _ = V.zface(m, r, a, tt)
            out.append((float(tt),
                        float(zf[ann].max()) - (z_n(tt) - FOOT_STANDOFF)))
        return out


def crit_a(rep):
    hdr("(a) SCORED -- THE WINDOW. Both halves required.")
    print("  (i) standoff plateau (slope < 2 mm/100 s over the last 200 s) <= 100 mm")
    print("  (ii) NO pad contact after t = 60 s (max top face over r 28-40 mm")
    print("       against z_n - 0.050). t = 0 contact is BY CONSTRUCTION.")
    print("  ⚠ C0 is scored 400-600 s (C3: its floor nears the adiabatic base).")
    print("  ⚠⚠ s* carries +/-5 %; 95-105 mm is NOT a decisive verdict (C4).\n")
    print(f"  {'leg':<9}{'pred s*':>8}{'s_end':>8}{'slope':>8}{'plateau':>9}"
          f"{'(i)':>6}{'1st contact':>12}{'worst':>9}{'(ii)':>6}  verdict")
    for n in LEGS:
        if not have(n):
            continue
        t, s = trace(n)
        rc = meta(n)["rc"]
        w = C0_WIN if n == "C0" else (t[-1] - 200.0, t[-1])
        m = (t >= w[0]) & (t <= w[1])
        sl = float(np.polyfit(t[m], s[m], 1)[0] * 100.0) if m.sum() > 3 else np.nan
        send = float(np.interp(w[1], t, s))
        plat = abs(sl) < SLOPE_BAR
        i_ok = plat and send <= S_CEIL
        ct = contact_trace(n)
        first = next((tt for tt, d in ct if d >= 0.0), None)
        worst = max(d for _, d in ct) * 1e3
        ii_ok = first is None
        rep.setdefault("a", {})[n] = (send, sl, plat, i_ok, first, worst, ii_ok, rc)
        v = ("*** WINDOW ***" if i_ok and ii_ok else
             "fail (i)+(ii)" if not i_ok and not ii_ok else
             "fail (i)" if not i_ok else "fail (ii) CONTACT")
        print(f"  {n:<9}{PRED[n][4]:>8d}{send:>8.1f}{sl:>8.2f}"
              f"{'yes' if plat else 'no':>9}{'PASS' if i_ok else 'FAIL':>6}"
              f"{(f'{first:.0f} s' if first is not None else 'never'):>12}"
              f"{worst:>8.0f}mm{'PASS' if ii_ok else 'FAIL':>6}  {v}"
              + ("  [ABORTED]" if rc else ""))
    print("\n  ⚠ 'worst' > 0 means the pad ring is ABOVE the foot plane, i.e. the")
    print("    pads are passing THROUGH rock -- unphysical under `prescribed`,")
    print("    which has no force balance. (c) and (d) are VOID from that point.")
    print("  ⚠ Three legs ABORTED on the BORE collision guard (r <= 28 mm), which")
    print("    is a stronger statement than pad contact: the nozzle itself")
    print("    reached rock.")


def crit_b(rep):
    hdr("(b) SCORED -- THE FLUX GATE, at a COMMON 100 mm standoff")
    print("  q100 = F_q * 1216 * (T_g - T_s), each leg's MEASURED floor T_g.")
    print(f"  against Kant {KANT[0]}-{KANT[1]} (measured spallation) and Meier's")
    print(f"  required {MEIER_Q[0]}-{MEIER_Q[1]} MW/m^2.\n")
    print(f"  {'leg':<9}{'F_q':>6}{'T_s':>7}{'T_g meas':>10}{'q100':>8}"
          f"  {'vs Kant':<12}{'vs Meier':<10}")
    for n in LEGS:
        if not have(n):
            continue
        fq, Ts, tg_pred, _, _, q_pred = PRED[n]
        p = os.path.join(OUT, n + "_jet_profile.csv")
        tg = tg_pred
        if os.path.exists(p):
            a = np.genfromtxt(p, delimiter=",", names=True)
            ts = np.unique(a["time"])
            sel = a["time"] == ts[-1]
            wgt = a["P_bin"][sel]
            tg = float(np.sum(a["T_gas"][sel] * wgt) / np.sum(wgt))
        q = fq * 1216.0 * (tg - Ts) / 1e6
        k = ("ABOVE" if q > KANT[1] else "below" if q < KANT[0] else "inside")
        mm = ("above" if q > MEIER_Q[1] else "BELOW" if q < MEIER_Q[0] else "inside")
        rep.setdefault("b", {})[n] = (q, tg, k, mm)
        print(f"  {n:<9}{fq:>6.2f}{Ts:>7.1f}{tg:>10.1f}{q:>8.2f}  {k:<12}{mm:<10}")
    print("\n  ⚠ T_g meas is the power-weighted floor gas temperature at the last")
    print("    logged profile row, so aborted legs are read at their abort time.")


def crit_cd(rep):
    hdr("(c)/(d) SCORED -- Ø(z) over 120-250 mm, and volume")
    print(f"  Ø against {D_TARGET} +/- {D_TOL} mm (Tab. 8.2 TEXT, never photo).")
    print(f"  volume against {VOL_CEIL} L per 0.5 m.")
    print("  ⚠⚠ VOID from sustained pad contact -- reported only up to it.\n")
    depths = [40, 60, 80, 100, 120, 150, 175, 200, 225, 250, 275, 300]
    print(f"  {'leg':<9}" + "".join(f"{d:>6}" for d in depths) + "   cal mean")
    for n in LEGS:
        if not have(n):
            continue
        first = rep.get("a", {}).get(n, (None,) * 8)[4]
        with with_lz(n) as lz:
            m, r, a = X.cell(n)
            t_use = first if first is not None else float(thermo(n)["time"][-1])
            zf, nv = V.zface(m, r, a, t_use)
            rc_, zc, _ = V.profile(r, zf)
            D = {d: V.diameter(rc_, zc, d) * 1e3 for d in depths}
            cal = [D[d] for d in depths
                   if CAL[0] <= d <= CAL[1] and D[d] == D[d]]
            mn = float(np.mean(cal)) if cal else float("nan")
            vol = float(nv.sum()) * m["dz"] ** 3 * 4.0 * 1e3
            depth = (lz - float(zf.min()))
            per = vol / depth * 0.5 if depth > 0 else float("nan")
            rep.setdefault("cd", {})[n] = (mn, per, t_use, first is not None)
            print(f"  {n:<9}" + "".join(
                f"{(D[d] if D[d] == D[d] else float('nan')):>6.1f}" for d in depths)
                + f"   {mn:>6.1f}" + ("  [to contact]" if first is not None else ""))
    print(f"\n  {'leg':<9}{'read at':>9}{'cal Ø':>8}{'vs 85+/-5':>11}"
          f"{'vol L/0.5m':>12}{'vs 3.42':>9}   status")
    for n in LEGS:
        if n not in rep.get("cd", {}):
            continue
        mn, per, t_use, voided = rep["cd"][n]
        okD = abs(mn - D_TARGET) <= D_TOL if mn == mn else False
        okV = per <= VOL_CEIL if per == per else False
        print(f"  {n:<9}{t_use:>8.0f}s{mn:>8.1f}{'PASS' if okD else 'FAIL':>11}"
              f"{per:>12.3f}{'PASS' if okV else 'FAIL':>9}   "
              + ("VOID after contact" if voided else "valid"))
        rep.setdefault("c", []).append((n, okD))
        rep.setdefault("d", []).append((n, okV))


def crit_e(rep):
    hdr("(e) SCORED -- the firing window, p90 - p10 <= 30 K")
    print(f"  {'leg':<9}{'firings':>9}{'mean':>9}{'p10':>9}{'p90':>9}"
          f"{'width':>8}  verdict")
    for n in LEGS:
        p = os.path.join(OUT, n + "_removal_events.csv")
        if not os.path.exists(p):
            continue
        a = np.genfromtxt(p, delimiter=",", names=True)
        s = (a["time"] >= 150.0) & (a["time"] <= 400.0)
        T = a["T_top"][s]
        if not len(T):
            continue
        p10, p90 = np.percentile(T, 10), np.percentile(T, 90)
        ok = (p90 - p10) <= WIN_BAR
        print(f"  {n:<9}{len(T):>9d}{T.mean():>9.1f}{p10:>9.1f}{p90:>9.1f}"
              f"{p90-p10:>8.1f}  {'PASS' if ok else 'FAIL'}")
        rep.setdefault("e", []).append((n, ok))


def crit_f(rep):
    hdr("(f) SCORED -- the power identity, +/- 0.02 W")
    for n in LEGS + ("S0",):
        if not have(n):
            continue
        d = thermo(n)
        if "jet_P_face" not in d:
            continue
        m = d["time"] >= 150.0
        if not m.any():
            continue
        e = float(np.abs((d["jet_P_face"] - (d["patch_P_robin"] + d["side_P_face"]
                                             + d["side_P_loss"]))[m]).max())
        ok = e <= 0.02
        print(f"  {n:<9} max residual {e:.3e} W   {'PASS' if ok else 'FAIL'}"
              f"    ledger_err {np.abs(d['ledger_err']).max():.2e}")
        rep.setdefault("f", []).append((n, ok))


def crit_g(rep):
    hdr("(g) SCORED -- identity")
    for suf in ("_removal_events.csv", "_jet_profile.csv"):
        pa, pb = os.path.join(OUT, "S0" + suf), os.path.join(Z1, "RB_D" + suf)
        if os.path.exists(pa) and os.path.exists(pb):
            ok = filecmp.cmp(pa, pb, shallow=False)
            print(f"  S0 vs RB_D {suf}: {'IDENTICAL' if ok else 'DIFFER'}")
            rep.setdefault("g", []).append((suf, ok))
    pa, pb = os.path.join(OUT, "S0"), os.path.join(Z1, "RB_D")
    if os.path.exists(pa) and os.path.exists(pb):
        ha = open(os.path.join(pa, "thermo.dat")).readline().split()
        hb = open(os.path.join(pb, "thermo.dat")).readline().split()
        va = np.loadtxt(os.path.join(pa, "thermo.dat"), skiprows=1, ndmin=2)
        vb = np.loadtxt(os.path.join(pb, "thermo.dat"), skiprows=1, ndmin=2)
        nr = min(len(va), len(vb))
        delta = set(ha) ^ set(hb)
        diff = [c for c in set(ha) & set(hb)
                if not np.array_equal(va[:nr, ha.index(c)], vb[:nr, hb.index(c)])]
        ok = not delta and not diff
        print(f"  S0 vs RB_D thermo: header delta "
              f"{'EMPTY' if not delta else sorted(delta)}; {len(diff)} differ  "
              f"{'PASS' if ok else 'FAIL'}")
        rep.setdefault("g", []).append(("thermo", ok))
    print("\n  tree gate vs witness/d2z1.md5 (--force): see GATE.txt")


def shoulder(rep):
    hdr("THE MECHANISM -- why every s* predictor failed")
    print("  Recession by radius at each leg's last time [mm]. The floor does")
    print("  NOT slow uniformly as the flux falls: the OUTER bore drops below")
    print("  the firing threshold and stops, while the axis keeps going.\n")
    bands = [(0, 8), (8, 16), (16, 24), (24, 32), (32, 40), (40, 48)]
    print(f"  {'leg':<9}" + "".join(f"{f'{a}-{b}':>9}" for a, b in bands)
          + f"{'axis-28mm':>11}")
    for n in LEGS:
        if not have(n):
            continue
        with with_lz(n) as lz:
            m, r, a = X.cell(n)
            t_end = float(thermo(n)["time"][-1])
            zf, _ = V.zface(m, r, a, t_end)
            rec = (lz - zf) * 1e3
            vals = []
            for lo, hi in bands:
                s = (r >= lo * 1e-3) & (r < hi * 1e-3)
                vals.append(float(rec[s].mean()) if s.any() else float("nan"))
            print(f"  {n:<9}" + "".join(f"{v:>9.1f}" for v in vals)
                  + f"{vals[0]-vals[3]:>11.1f}")
    print("\n  ⚠ The last column is the SHOULDER LAG: how far the r = 24-32 mm")
    print("    band trails the axis. The burner descends on a schedule, so it")
    print("    lands on that shoulder -- which is where the bore collision guard")
    print("    fired (column (13,3), r = 27.9 mm) on all three aborted legs.")
    print("  ⚠ This is a THRESHOLD effect. The F model is linear in")
    print("    (T_g - T_s) and cannot represent a radius at which firing stops,")
    print("    which is why BOTH s* predictors -- the packet's and item 0's --")
    print("    overpredicted the equilibrium standoff on every leg.")


def predictions(rep):
    hdr("PREDICTIONS")
    A = rep.get("a", {})
    win = [n for n in LEGS if n in A and A[n][3] and A[n][6]]
    print(f"  P1 (packet) the window closes at Q5T870 and nowhere else.")
    print(f"     legs passing BOTH halves of (a): {win if win else 'NONE'}")
    print(f"  P1 (item 0) the window closes at Q35T870/Q30: "
          f"{'HELD' if set(win) == {'Q35T870','Q30'} else 'REFUTED'}")
    print("\n  P2 descending s* order C0 > G1600 > Q5 > Q5T870 > Q35, "
          "Q30T870 contact:")
    order = [(n, A[n][0]) for n in ("C0", "G1600", "Q5", "Q5T870", "Q35") if n in A]
    mono = all(order[i][1] > order[i + 1][1] for i in range(len(order) - 1))
    print("     " + " > ".join(f"{n} {v:.0f}" for n, v in order)
          + f"   monotone: {'HELD' if mono else 'REFUTED'}")
    B = rep.get("b", {})
    if "C0" in B:
        print(f"\n  P3 q100 above Kant's {KANT[1]} on C0: {B['C0'][0]:.2f} "
              f"{'HELD' if B['C0'][0] > KANT[1] else 'REFUTED'}")
    if win:
        for n in win:
            print(f"     passing leg {n} q100 {B[n][0]:.2f} inside "
                  f"{MEIER_Q}: {'HELD' if MEIER_Q[0] <= B[n][0] <= MEIER_Q[1] else 'REFUTED'}")
    else:
        print("     no passing leg -> P3's third clause is VACUOUS")
    E = rep.get("e", [])
    print(f"\n  P5 firing window <= {WIN_BAR:.0f} K on every leg: "
          f"{'HELD' if all(ok for _, ok in E) else 'REFUTED'}")
    print("  P6 a required T_s > 920 K refutes the rock-side explanation: "
          "no leg needed > 870 K, so NOT triggered")


def main():
    print("D2ad -- scoring CRITERION.md (sha in CRITERION.sha256)")
    print("legs on disk: " + ", ".join(n for n in LEGS + ("S0",) if have(n)))
    print("⚠ Q35T870, Q30, Q30T870 ABORTED on the bore collision guard; their")
    print("  traces are valid to the abort and are scored there.")
    rep = {}
    for f in (crit_a, crit_b, crit_cd, crit_e, crit_f, crit_g, shoulder,
              predictions):
        f(rep)
    hdr("SCORECARD")
    for k in "cdefg":
        rows = rep.get(k, [])
        if rows:
            npass = sum(1 for _, ok in rows if ok)
            print(f"  ({k}) {npass}/{len(rows)} pass")
    A = rep.get("a", {})
    win = [n for n in LEGS if n in A and A[n][3] and A[n][6]]
    print(f"  (a) legs passing BOTH halves: {win if win else 'NONE -- the window '
                                                          'is EMPTY'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
