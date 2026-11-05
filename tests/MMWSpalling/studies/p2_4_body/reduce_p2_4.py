#!/usr/bin/env python3
"""P2-4 reductions, rows c2-c9.

⚠⚠ REUSED, NOT RE-DERIVED. `reduce_price` owns the charged volume and the firing
population, `bore_profile` owns the void-mask geometry, `wall_profile` owns the
per-bin gas profile. Only paths are parameters.

⚠⚠ Ø_eq IS AN AREA STATISTIC AND IT IS NOT THE SECTOR-MEAN MASK DIAMETER.
c4 is scored on Ø_eq = 2*sqrt(A/pi) with A the VOIDED AREA of the 0-10 mm layer
scaled to the full disc, resolution ~1-2 mm. The sector-mean mask diameter that
P2-3 reported for the mouth had a sector SPREAD of 10.6 mm, i.e. wider than the
effect c4 is trying to measure, which is exactly why the packet changed the
statistic.

⚠⚠ THE BORE IS A CONE THAT CROSSES MEIER, NOT A PROFILE THAT MATCHES IT
(must-fix 2). So c5 reports the TAPER in mm/mm and the per-depth profile, and
NEVER a band mean -- a mean over a crossing cone is a number with no referent.

⚠ Per-step columns are MEANS over their window, never sums (standing rule 7(i)).
`body_noLaw_n` and `ledger_E_flake_wall` are CUMULATIVE and are DIFFERENCED.
⚠ A row that cannot be computed prints "not computed", never a verdict.
"""
import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
P22C = os.path.join(ROOT, "tests", "MMWSpalling", "studies", "p2_2c_meier")
P23 = os.path.join(ROOT, "tests", "MMWSpalling", "studies", "p2_3_wall")
for d in (P22C, P23, HERE):
    sys.path.insert(0, d)
import reduce_price as RP        # noqa: E402
import bore_profile as BP        # noqa: E402
import wall_profile as WP        # noqa: E402

LEGS_DIR = os.path.join(HERE, "legs")
WATER_W_PER_K = 186.0      # 160 L/h; meier_thesis_notes L43-45
WIN1, WIN2 = (150.0, 350.0), (500.0, 800.0)
T_CLAMP = 873.15
DX = DZ = 0.002
Z_TOP = 0.500


def thermo(d):
    p = os.path.join(d, "run", "thermo.dat")
    if not os.path.exists(p):
        return None, []
    lines = [l.split() for l in open(p).read().splitlines() if l.split()]
    return lines[0], lines[1:]


def col(head, rows, nm):
    if not head or nm not in head:
        return []
    i = head.index(nm)
    out = []
    for r in rows:
        if len(r) > i:
            try:
                out.append(float(r[i]))
            except ValueError:
                pass
    return out


def win_mean(head, rows, nm, lo, hi):
    t, v = col(head, rows, "time"), col(head, rows, nm)
    sel = [x for ti, x in zip(t, v) if lo <= ti <= hi]
    return (sum(sel) / len(sel)) if sel else None


def at_time(head, rows, nm, tt):
    t, v = col(head, rows, "time"), col(head, rows, nm)
    if not t or not v:
        return None
    i = min(range(min(len(t), len(v))), key=lambda j: abs(t[j] - tt))
    return v[i]


def cum_delta(head, rows, nm, lo, hi):
    """⚠ CUMULATIVE column: the window's value is a DIFFERENCE, not a sum."""
    t, v = col(head, rows, "time"), col(head, rows, nm)
    if not t or not v:
        return None
    n = min(len(t), len(v))
    a = [j for j in range(n) if t[j] <= lo]
    b = [j for j in range(n) if t[j] <= hi]
    if not b:
        return None
    return v[b[-1]] - (v[a[-1]] if a else 0.0)


def events(d):
    p = os.path.join(d, "run_removal_events.csv")
    return list(csv.DictReader(open(p))) if os.path.exists(p) else []


def plotfiles(d):
    r = os.path.join(d, "run")
    if not os.path.isdir(r):
        return []
    return sorted(os.path.join(r, x) for x in os.listdir(r)
                  if x.endswith("cell"))


def layer_Deq(pf, nx, ny, nz, depth_mm=10.0):
    """Ø_eq of the 0-`depth_mm` layer: mean voided area per level, x4 for the
    quarter, then 2*sqrt(A/pi)."""
    mask = BP.void_mask(pf, nx, ny, nz)
    klo = int(round((Z_TOP - depth_mm * 1e-3 - (Z_TOP - nz * DZ)) / DZ))
    khi = nz - 1
    if klo < 0:
        klo = 0
    tot, nk = 0, 0
    for k in range(klo, khi + 1):
        c = sum(1 for i in range(nx) for j in range(ny) if mask[i][j][k])
        if c:
            tot += c
            nk += 1
    if not nk:
        return None
    A = (tot / nk) * DX * DX * 4.0
    return 2.0 * math.sqrt(A / math.pi) * 1e3


def taper(prof, lo_mm, hi_mm):
    """least-squares slope of Ø against depth, mm/mm, over [lo, hi]."""
    pts = [(d, v) for d, v in prof if lo_mm <= d <= hi_mm and v is not None]
    n = len(pts)
    if n < 3:
        return None
    sx = sum(d for d, _ in pts); sy = sum(v for _, v in pts)
    sxx = sum(d * d for d, _ in pts); sxy = sum(d * v for d, v in pts)
    den = n * sxx - sx * sx
    return ((n * sxy - sx * sy) / den) if abs(den) > 0 else None


def clamp_by_volume(ev):
    va = vt = 0.0
    for r in ev:
        try:
            if int(float(r.get("regime", -1))) != RP.R_FLOOR:
                continue
            T = float(r.get("T_surf") or 0.0)
            n = float(r.get("n_voided") or 0.0)
        except (TypeError, ValueError):
            continue
        if T <= 0 or n <= 0:
            continue
        vt += n
        if T > T_CLAMP:
            va += n
    return va * RP.V_CELL, vt * RP.V_CELL, (va / vt if vt else None)


def cap_share(ev, hmax):
    """fraction of firings sitting AT the flake cap, by radius band and regime.
    A diagnostic carried for the PARKED cap step -- no bar."""
    bands = [(0, 20), (20, 28), (28, 40), (40, 48), (48, 60)]
    out = {}
    for r in ev:
        try:
            reg = int(float(r.get("regime", -1)))
            i = int(float(r["col_i"])); j = int(float(r["col_j"]))
        except (TypeError, ValueError, KeyError):
            continue
        rad = 1e3 * RP.col_radius(i, j)
        h = r.get("h_f") or r.get("h_f_lat") or ""
        try:
            hv = float(h)
        except ValueError:
            continue
        for lo, hi in bands:
            if lo <= rad < hi:
                k = ("%d-%d" % (lo, hi), reg)
                a, b = out.get(k, (0, 0))
                out[k] = (a + (1 if hv >= hmax * 0.999 else 0), b + 1)
                break
    return out


def verdict(tag, ok, note=""):
    if ok is None:
        print("     not computed  %s%s" % (tag, ("  -- " + note) if note else ""))
    else:
        print("     %s %s%s" % ("ok  " if ok else "MISS", tag,
                                ("  -- " + note) if note else ""))


def one(tag):
    d = os.path.join(LEGS_DIR, tag)
    if not os.path.isdir(d):
        print("\n--- %s: no leg directory" % tag)
        return None
    head, rows = thermo(d)
    ev = events(d)
    print("\n" + "-" * 72)
    print("LEG %s  (%s)  %d firing events"
          % (tag, os.path.relpath(d, ROOT), len(ev)))
    if not head:
        print("  not computed -- no thermo.dat")
        return None
    t = col(head, rows, "time")
    t_end = t[-1] if t else 0.0
    print("  ran to t = %.1f s" % t_end)
    res = {"tag": tag, "t_end": t_end}
    nx, ny, nz, dx, dz, z_top = BP.geometry()

    # ---- c2: the mouth gas at the end ----
    tm_end = at_time(head, rows, "wall_T_mouth", t_end)
    tm_w1 = win_mean(head, rows, "wall_T_mouth", *WIN1)
    if tm_end is None:
        verdict("c2", None, "wall_T_mouth absent (march off on this leg)")
    else:
        print("  c2 wall_T_mouth: %.1f K at t=%.0f s; window %s mean %.1f K"
              % (tm_end, t_end, WIN1, tm_w1 or 0.0))
        verdict("c2  (bar 1190-1270 K at 800 s)",
                (1190.0 <= tm_end <= 1270.0) if t_end >= 790 else None,
                "" if t_end >= 790 else "leg stopped at %.0f s" % t_end)
        res["T_mouth_end"] = tm_end

    # ---- c3: body power and the water rise ----
    for nm, lo, hi, bar in (("500 s", 400.0, 500.0, (650.0, 850.0)),
                            ("800 s", 700.0, 800.0, (1150.0, 1650.0))):
        pw = win_mean(head, rows, "jet_P_body", lo, hi)
        if pw is None:
            verdict("c3 @%s" % nm, None, "jet_P_body absent")
            continue
        full = pw * 4.0
        print("  c3 @%s: jet_P_body %.1f W quarter = %.2f kW full; water rise "
              "%.1f K (160 L/h = %.0f W/K)"
              % (nm, pw, full / 1e3, full / WATER_W_PER_K, WATER_W_PER_K))
        verdict("c3 @%s  (bar %.2f-%.2f kW quarter)" % (nm, bar[0] / 1e3, bar[1] / 1e3),
                bar[0] <= pw <= bar[1] if t_end >= hi - 10 else None,
                "" if t_end >= hi - 10 else "leg stopped at %.0f s" % t_end)
        res["P_body_%s" % nm.replace(" ", "")] = pw
        if full / WATER_W_PER_K > 60.0:
            print("     ⚠ water rise over 60 K -- the sink is implausibly strong "
                  "and the COEFFICIENT, not the physics, is the question")

    # ---- the body's own state ----
    for nm in ("body_h_used", "body_z_imm", "body_L"):
        v1 = win_mean(head, rows, nm, *WIN1)
        v2 = win_mean(head, rows, nm, *WIN2)
        if v1 is not None:
            print("     %-12s window1 %.4f   window2 %s"
                  % (nm, v1, ("%.4f" % v2) if v2 is not None else "-"))
    nl = cum_delta(head, rows, "body_noLaw_n", 0.0, t_end)
    if nl is not None:
        print("     body_noLaw_n over the whole run: %.0f steps (CUMULATIVE, "
              "differenced not summed)" % nl)

    # ---- c4: the top-layer Ø_eq ----
    pfs = plotfiles(d)
    if pfs:
        # ⚠ EVERY plotfile, with its TIME. The first version sampled three
        # (first / middle / last) and reported the same value twice, which hid
        # the fact that the trajectory PLATEAUS -- the whole point of the row.
        # step -> time from the run's own dt (total steps / stop_time).
        import re
        steps = [int(re.search(r"(\d+)cell", os.path.basename(x)).group(1))
                 for x in pfs]
        dt = (t_end / steps[-1]) if steps and steps[-1] else 0.0
        print("  c4 top-layer Ø_eq (AREA statistic: mean voided area per level "
              "over the 0-10 mm layer, x4 for the quarter, Ø = 2*sqrt(A/pi)):")
        traj = []
        for pf, st in zip(pfs, steps):
            dq = layer_Deq(pf, nx, ny, nz)
            if dq:
                traj.append((st * dt, dq))
        for tt, dq in traj:
            print("       t = %6.1f s   Ø_eq = %7.2f mm" % (tt, dq))
        res["Deq_end"] = traj[-1][1] if traj else None
        # the growth rate over the packet's own window
        w = [(tt, dq) for tt, dq in traj if WIN2[0] - 1 <= tt <= WIN2[1] + 1]
        if len(w) >= 2 and w[-1][0] > w[0][0]:
            g = (w[-1][1] - w[0][1]) / (w[-1][0] - w[0][0])
            print("     growth over %.0f-%.0f s = %.4f mm/s (%+.1f mm) "
                  "(bar <= 0.042; W1's realised 0.085 over 300-500 s)"
                  % (w[0][0], w[-1][0], g, w[-1][1] - w[0][1]))
            verdict("c4 growth  (bar <= 0.042 mm/s)", g <= 0.042)
            res["Deq_growth"] = g
        else:
            verdict("c4 growth", None, "no plotfile pair inside %s" % (WIN2,))
        if traj:
            at500 = min(traj, key=lambda x: abs(x[0] - 500.0))
            print("     Ø_eq at %.0f s = %.2f mm  (bar 128-132)"
                  % (at500[0], at500[1]))
            verdict("c4 level  (bar 128-132 mm)", 128.0 <= at500[1] <= 132.0)
            print("     ⚠⚠ THE LEVEL BAR'S PROVENANCE DOES NOT REPRODUCE. It was")
            print("        set against a W1 figure of 134.7 mm at 500 s; with the")
            print("        definition the packet states, W1 measures 117.37 mm")
            print("        here -- and that agrees with W1's own sector-mean mouth")
            print("        (117.8), which cross-validates this estimator. The")
            print("        GROWTH half of c4 is unaffected: it is a difference")
            print("        within one estimator and does not depend on the level.")
    else:
        verdict("c4", None, "no plotfiles")

    # ---- c5: the bore, per depth AND the taper ----
    print("  c5 bore Ø(z) from the void mask, per azimuth sector:")
    BP.main(tag=tag, leg_root=d)
    if pfs:
        mask = BP.void_mask(pfs[-1], nx, ny, nz)
        prof = []
        for d_mm in (10, 25, 50, 75, 100, 125, 150, 175, 200, 250, 300, 350, 400):
            z = z_top - d_mm * 1e-3
            k = int((z - (z_top - nz * dz)) / dz)
            if not (0 <= k < nz):
                continue
            rs = [x for x in BP.sector_radii(mask, k, nx, ny, dx) if x is not None]
            prof.append((float(d_mm), 2e3 * (sum(rs) / len(rs)) if rs else None))
        t1 = taper(prof, 25.0, 150.0)
        above = [p for p in prof if p[1] is not None]
        t2 = taper(prof, above[0][0], above[-1][0]) if len(above) >= 3 else None
        print("     TAPER (least squares, mm of Ø per mm of depth):")
        print("       25-150 mm        : %s" % (("%.3f" % t1) if t1 else "n/a"))
        print("       whole section    : %s" % (("%.3f" % t2) if t2 else "n/a"))
        print("       ⚠ W1 measured 0.302 mm/mm over 140-280 mm; W074 0.166.")
        print("       ⚠ Meier's own slope over 25-450 mm is about -0.027 mm/mm")
        print("         (92 -> 80 mm), i.e. nearly PARALLEL-SIDED. A band MEAN")
        print("         over a crossing cone has no referent, which is why the")
        print("         taper is the scored shape statistic (must-fix 2).")
        res["taper_25_150"] = t1

    # ---- c6: the wall's share of removed volume ----
    cv = RP.charged_volume(ev)
    tot = cv["charged"]
    share = (cv["wall"]["vol_nvoided"] / tot) if tot > 0 else None
    print("  c6 charged volume (pad clearing EXCLUDED): floor %.5f L + wall "
          "%.5f L = %.5f L" % (1e3 * cv["floor"]["vol_nvoided"],
                               1e3 * cv["wall"]["vol_nvoided"], 1e3 * tot))
    print("     wall share %s %%  (W1 42 %% over 150-350 s)"
          % (("%.1f" % (100 * share)) if share else "n/a"))
    verdict("c6  (bar 20-30 %% over 500-800 s)",
            (0.20 <= share <= 0.30) if (share and t_end >= 790) else None,
            "whole-run share; the windowed split needs the event log by time")
    res["wall_share"] = share

    # ---- c7: the rate, both windows ----
    for lo, hi in (WIN1, WIN2):
        tt, nz_ = col(head, rows, "time"), col(head, rows, "nozzle_z")
        pts = [(a, b) for a, b in zip(tt, nz_) if lo <= a <= hi and b > 0]
        if len(pts) > 1 and pts[-1][0] > pts[0][0]:
            v = (pts[0][1] - pts[-1][1]) / (pts[-1][0] - pts[0][0]) * 3600.0
            print("  c7 rate %.0f-%.0f s = %.3f m/h" % (lo, hi, v))
            res["rate_%d" % int(lo)] = v
    print("     ⚠ P5's two branches: A = S0 within +/-0.2 of W1's 1.63; "
          "B = S0 2.0-2.3. Either is a finding; the rate is NOT recovered by "
          "re-tuning. ⚠ PROVISIONAL at 2 mm (1 mm parked).")

    # ---- c8: the criterion ----
    ts, note = RP.firing_T(ev)
    if ts:
        med = RP.pctl(ts, 0.5)
        va, vt, fr = clamp_by_volume(ev)
        print("  c8 median firing T_surf %.2f K; clamp fraction BY VOLUME "
              "%.2f %% (%.5f of %.5f L)"
              % (med, 100 * fr if fr else 0.0, 1e3 * va, 1e3 * vt))
        print("     %s" % note)
        verdict("c8 firing T (845 +/- 10)", abs(med - 845.0) <= 10.0)
        verdict("c8 clamp (12-25 %)", (0.12 <= (fr or -1) <= 0.25))
        res["T_fire"] = med
    else:
        verdict("c8", None, "no firing carries a reconstructed T_surf")

    # ---- c9: price per regime + the cap share ----
    Ef = col(head, rows, "ledger_E_flake")
    Ew = col(head, rows, "ledger_E_flake_wall")
    if Ef and Ew:
        ew = Ew[-1]
        ef = Ef[-1] - ew
        vf, vw = cv["floor"]["vol_nvoided"], cv["wall"]["vol_nvoided"]
        print("  c9 price PER REGIME (ledger_E_flake_wall splits the total; "
              "floor = total - wall):")
        print("       floor  %.4e J / %.5f L = %s GJ/m^3"
              % (ef, 1e3 * vf, ("%.3f" % (ef / vf / 1e9)) if vf > 0 else "n/a"))
        print("       wall   %.4e J / %.5f L = %s GJ/m^3"
              % (ew, 1e3 * vw, ("%.3f" % (ew / vw / 1e9)) if vw > 0 else "n/a"))
        print("       ⚠ P2-3 reported ONE price over floor + wall while the wall "
              "carried 42 % of the volume; this is the split it owed.")
        verdict("c9", None, "diagnostic for the parked steps -- no bar")
    hmaxs = [k for k in (col(head, rows, "lat_qmax"),) if k]
    cs = cap_share(ev, 0.015)   # spall.h_spall_max on this deck
    if cs:
        print("     cap share by radius band and regime (at flake_h_max; "
              "diagnostic, no bar):")
        for (band, reg), (a, b) in sorted(cs.items()):
            if b >= 50:
                print("       r %-6s regime %-2d : %5.1f %% of %d firings"
                      % (band, reg, 100.0 * a / b, b))
    return res


def main(want=None):
    want = want or ["S0", "B1L"]
    print("P2-4 REDUCTIONS")
    print("=" * 72)
    print("⚠ windows %s and %s; per-step columns are MEANS, cumulative ones are "
          "DIFFERENCED." % (WIN1, WIN2))
    got = {}
    for tg in want:
        r = one(tg)
        if r:
            got[tg] = r
    print("\n" + "=" * 72)
    for k, v in got.items():
        print("  %-5s %s" % (k, ", ".join(
            "%s=%.4g" % (a, b) for a, b in v.items()
            if isinstance(b, float))))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:] or None))
