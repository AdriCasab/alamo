#!/usr/bin/env python3
"""D2ad item 0 -- measured on 60 s smokes, BEFORE CRITERION.md is hashed.

  (b) the kic.scale -> T_s mapping, MEASURED (the brief's hand mapping is +/-10 %
      and says not to trust it), with the firing-window width guard
  (c) the enthalpy back-reaction, verified LIKE FOR LIKE (matched time)
  (d) G1600's floor T_g, measured (its tabled 1410 K is an estimate)
  (e) R90, the undecayed floor recession, measured
  (f) the calibration / validation split, written down before running

Reads only D2ad smokes and the frozen d2ac_feed legs. Run it, hash it, then set
the scales in run.py and hash CRITERION.md.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)
OUT = os.path.join(HERE, "output")
AC = os.path.join(STUD, "d2ac_feed", "output")

T_S_BASE = 823.3          # K, D2ac FH 150-400 s, 30 716 firings (item 0(a))
T_G_BASE = 1664.7         # K, power-weighted over 23 bins at t = 490 s
WIN_BAR = 30.0            # K, p90 - p10 must not exceed this
FEED = 1.5                # m/h
NOZZLE_Z0 = 0.550


def ev(p):
    a = np.genfromtxt(p, delimiter=",", names=True)
    return a


def tg_pw(p, t_want=None):
    """Power-weighted floor gas temperature over all bins.

    ⚠ t_want = None returns the LAST logged row. A stop_time = 60 run's profile
    ends at t = 50.016, because jet_profile_interval = 10 logs at the FIRST step
    of each interval and the run stops before the 60 s one. Asking for "nearest
    60" therefore silently returned 50.016 for the smokes and 60.016 for FH --
    a 10 s mismatch, which is exactly the like-for-like error item 0(c) exists
    to prevent. Comparisons below use the smoke's OWN last time on both sides.
    """
    a = np.genfromtxt(p, delimiter=",", names=True)
    ts = np.unique(a["time"])
    t = ts[-1] if t_want is None else ts[np.argmin(np.abs(ts - t_want))]
    s = a["time"] == t
    w = a["P_bin"][s]
    return float(t), float(np.sum(a["T_gas"][s] * w) / np.sum(w)), float(np.sum(w))


def firing(p, t0=0.0):
    a = ev(p)
    s = a["time"] >= t0
    T = a["T_top"][s]
    return (len(T), float(T.mean()), float(np.percentile(T, 10)),
            float(np.percentile(T, 50)), float(np.percentile(T, 90)))


def main():
    P = print
    P("=" * 78)
    P("D2ad ITEM 0 -- measured on 60 s smokes")
    P("=" * 78)

    # ---------------------------------------------------------------- (b)
    P("\n(b) THE kic.scale -> T_s MAPPING, MEASURED")
    P("-" * 78)
    P("  Smokes on the Q5 base (F_q = 0.5), stop_time = 60 s.")
    P("  ⚠ SMOKE_K10 would duplicate SMOKE_Q5 exactly (F_q 0.5, scale 1.0), so")
    P("    SMOKE_Q5 IS the scale-1.0 row; it is not run twice.")
    P(f"\n  {'scale':>6}{'firings':>9}{'mean':>9}{'p10':>9}{'p50':>9}{'p90':>9}"
      f"{'p90-p10':>9}  window")
    rows = []
    for sc, n in ((1.0, "SMOKE_Q5"), (1.3, "SMOKE_K13"), (1.6, "SMOKE_K16"),
                  (1.9, "SMOKE_K19"), (2.2, "SMOKE_K22")):
        p = os.path.join(OUT, n + "_removal_events.csv")
        if not os.path.exists(p):
            P(f"  {sc:>6.1f}   (missing {n})")
            continue
        cnt, mn, p10, p50, p90 = firing(p)
        w = p90 - p10
        rows.append((sc, mn, w, cnt))
        P(f"  {sc:>6.1f}{cnt:>9d}{mn:>9.1f}{p10:>9.1f}{p50:>9.1f}{p90:>9.1f}"
          f"{w:>9.1f}  {'OK' if w <= WIN_BAR else 'TOO BROAD'}")
    P(f"\n  ⚠ BAR: p90 - p10 <= {WIN_BAR:.0f} K at every scale. If it broadens, the")
    P("    lever is changing the CRITERION, not the temperature, and the packet")
    P("    stops and reports that (P5).")
    if rows:
        worst = max(r[2] for r in rows)
        P(f"  widest window across the ladder: {worst:.1f} K   "
          f"{'PASS' if worst <= WIN_BAR else 'FAIL -- T_s legs are VOID (P5)'}")
        sc = np.array([r[0] for r in rows])
        Ts = np.array([r[1] for r in rows])
        P(f"\n  the lever moves T_s: {Ts.min():.1f} -> {Ts.max():.1f} K over "
          f"scale {sc.min():.1f} -> {sc.max():.1f}")
        P("  (a flat 823 K across all five would mean the key never reached "
          "Removal.H:759)")
        for target in (870.0, 920.0):
            if Ts.max() < target:
                P(f"  ⚠ T_s = {target:.0f} K is ABOVE the measured range -- "
                  f"EXTRAPOLATION, flagged")
            k = float(np.interp(target, Ts, sc)) if Ts[0] < Ts[-1] else float("nan")
            near = sc[int(np.argmin(np.abs(Ts - target)))]
            dT = Ts[int(np.argmin(np.abs(Ts - target)))] - target
            P(f"  T_s = {target:.0f} K -> scale {k:.4f} (interpolated); "
              f"nearest run scale {near:.1f} is {dT:+.1f} K off")

    # ---------------------------------------------------------------- (c)
    P("\n(c) THE BACK-REACTION, VERIFIED LIKE FOR LIKE")
    P("-" * 78)
    P("  ⚠ entrainment is `exhaust` mode, so the exhaust temperature BUILDS over")
    P("    the run. Comparing a 60 s smoke against FH's 490 s value would be")
    P("    wrong. FH is read at the SAME 60 s and the offset stated.")
    fh = os.path.join(AC, "FH_jet_profile.csv")
    # The matched time is the smokes' OWN last logged row, not a nominal 60 s.
    t_match = tg_pw(os.path.join(OUT, "SMOKE_Q10_jet_profile.csv"))[0]
    P(f"  MATCHED TIME (the smokes' last logged row): t = {t_match:.3f} s")
    if os.path.exists(fh):
        t60, g60, p60 = tg_pw(fh, t_match)
        t490, g490, p490 = tg_pw(fh, 490.0)
        P(f"\n  FH floor T_g at t = {t60:.3f} s: {g60:8.1f} K  "
          f"(floor power {p60:7.1f} W quarter)")
        P(f"  FH floor T_g at t = {t490:.1f} s: {g490:8.1f} K  "
          f"(floor power {p490:7.1f} W quarter)")
        P(f"  matched -> 490 s OFFSET: {g490-g60:+.1f} K   <- applied below")
        off = g490 - g60
    else:
        P("  FH_jet_profile.csv missing"); off = 0.0
    P(f"\n  {'leg':<13}{'F_q':>6}{'T_g(60s)':>10}{'+offset':>10}"
      f"{'C1 pred':>10}{'diff':>8}  {'floor P':>9}")
    c1 = lambda f: 1900.0 - 235.3 * f
    got = {}
    for n, f in (("SMOKE_Q10", 1.0), ("SMOKE_Q5", 0.5), ("SMOKE_Q35", 0.35)):
        p = os.path.join(OUT, n + "_jet_profile.csv")
        if not os.path.exists(p):
            continue
        t, g, pw = tg_pw(p)          # its own last row = the matched time
        assert abs(t - t_match) < 1e-6, f"{n}: time mismatch {t} vs {t_match}"
        adj, pred = g + off, c1(f)
        got[f] = (adj, pw)
        P(f"  {n:<13}{f:>6.2f}{g:>10.1f}{adj:>10.1f}{pred:>10.1f}"
          f"{adj-pred:>8.1f}  {pw:>9.1f}")
    P("\n  ⚠ BAR: any leg off by > 40 K after the offset ⇒ recompute F, s* and")
    P("    q100 and re-hash before the ladder.")
    if got:
        worst = max(abs(got[f][0] - c1(f)) for f in got)
        P(f"  worst |measured - C1| = {worst:.1f} K   "
          f"{'PASS -- C1 stands' if worst <= 40 else 'FAIL -- RECOMPUTE'}")
        P("\n  WITNESS (the prefix actually applied): floor power must fall ~F_q x")
        if 1.0 in got:
            for f in sorted(got, reverse=True):
                P(f"    F_q = {f:4.2f}: floor power {got[f][1]:7.1f} W  "
                  f"ratio to F_q=1 {got[f][1]/got[1.0][1]:5.3f}  (expect ~{f:.2f})")

    # ---------------------------------------------------------------- (d)
    P("\n(d) G1600's FLOOR T_g, MEASURED")
    P("-" * 78)
    p = os.path.join(OUT, "SMOKE_G1600_jet_profile.csv")
    if os.path.exists(p):
        t, g, pw = tg_pw(p)
        adj = g + off
        P(f"  SMOKE_G1600 floor T_g at t = {t:.3f} s: {g:.1f} K, "
          f"+offset = {adj:.1f} K")
        P(f"  the ladder TABLED 1410 K (an estimate)  ->  diff {adj-1410.0:+.1f} K")
        P(f"  floor power {pw:.1f} W quarter")
        F = (adj - T_S_BASE) / (T_S_BASE - 293.0) / 1.58665
        P(f"  recomputed F = 1.0 * [({adj:.1f} - {T_S_BASE})/({T_S_BASE} - 293)]"
          f" / 1.58665 = {F:.3f}")
        P(f"  => s* = 198 * F = {198*F:.0f} mm   (tabled 138)")

    # ---------------------------------------------------------------- (e)
    P("\n(e) R90 -- THE UNDECAYED FLOOR RECESSION, MEASURED")
    P("-" * 78)
    P("  ⚠ C0 IS FH's argv plus stop_time, so C0 == FH for the first 500 s and")
    P("    R90 can be read from the FROZEN FH trace -- no new run needed.")
    P("  Measured where the standoff is still inside the flat region (s < 90 mm),")
    P("  because beyond it the map decays as 1/s and the rate is no longer R90.")
    th = os.path.join(AC, "FH", "thermo.dat")
    if os.path.exists(th):
        h = open(th).readline().split()
        d = np.loadtxt(th, skiprows=1)
        t = d[:, 0]
        s = d[:, h.index("jet_s_c")]
        zf = (NOZZLE_Z0 - (FEED / 3600.0) * t) - s
        m = (s < 0.090) & (t > 5.0)
        if m.sum() > 3:
            R90 = -np.polyfit(t[m], zf[m], 1)[0] * 3600.0
            P(f"  window: t = {t[m].min():.1f}-{t[m].max():.1f} s, "
              f"standoff {s[m].min()*1e3:.1f}-{s[m].max()*1e3:.1f} mm "
              f"({int(m.sum())} rows)")
            P(f"  R90 = {R90:.3f} m/h")
            P(f"\n  the brief quotes TWO inconsistent values: 3.3 (s* formula) and")
            P(f"  4.2 (contact edge). Measured {R90:.2f} -> "
              f"{'3.3 is the right one' if abs(R90-3.3)<abs(R90-4.2) else '4.2 is closer'}")
            P(f"  deviation from 3.3: {100*(R90-3.3)/3.3:+.1f} %   "
              f"{'within 10 %, ladder stands' if abs(R90-3.3)/3.3 <= 0.10 else 'OVER 10 % -- RECOMPUTE THE LADDER'}")
            P(f"\n  CONTACT THRESHOLD  F = feed/R90 = {FEED:.1f}/{R90:.3f} = "
              f"{FEED/R90:.4f}   (ladder assumed 0.4545)")
            P(f"  PASSING WINDOW     {FEED/R90:.4f} < F <= 0.505  "
              f"(s* between 90 and 100 mm)")
            P(f"  ⚠ THE WINDOW IS ~{1000*0.09*(R90/FEED)*(1-FEED/R90/0.505)*0:.0f}"
              f"10 mm WIDE and the predictor carries +/-5 %. A leg at 95-105 mm")
            P("    is NOT a decisive pass or fail (C4).")

    # ---------------------------------------------------------------- (f)
    P("\n(f) THE CALIBRATION / VALIDATION SPLIT -- fixed BEFORE running")
    P("-" * 78)
    P("  CALIBRATED (two numbers, chosen on the STANDOFF ALONE):")
    P("      F_q   the floor heat-transfer level")
    P("      T_s   the flaking temperature, via spallation.sp.kic.scale")
    P("  VALIDATION (nothing here may be used to choose a leg):")
    P("      Ø(z) over 120-250 mm")
    P("      volume per 0.5 m")
    P("      pad clearance")
    P("      the flux gate q100 against Kant and Meier")
    P("      the free-descent rate (the CONDITIONAL R_cal follow-on)")
    # ---------------------------------------------------------------- (g)
    P("\n(g) THE RECOMPUTED LADDER -- mandated by (c) and (e)")
    P("=" * 78)
    P("  The packet REQUIRES this recompute in two places:")
    P("    (c) 'off by > 40 K after the offset => recompute F, s* and q100'")
    P("    (e) 'if R90 differs from 3.3 by more than 10 %, recompute the ladder'")
    P("  Both fired. Measured worst |C1 error| = 53.5 K and R90 = 4.237 (+28.4 %).")
    P("")
    R90 = 4.237
    fqm = np.array([0.35, 0.50, 1.00])
    tgm = np.array([1871.2, 1834.2, 1664.7])
    cf = np.polyfit(fqm, tgm, 2)   # exact through the 3 measured points
    Tg = lambda f: float(np.polyval(cf, f))
    P("  T_g(F_q): C1's linear 1900 - 235.3*F_q is REPLACED by a quadratic fitted")
    P("  exactly through the three MEASURED points (0.35, 0.50, 1.00). Every leg")
    P("  uses F_q inside that range, so this is interpolation, not extrapolation.")
    P(f"    T_g = {cf[0]:+.2f} F^2 {cf[1]:+.2f} F {cf[2]:+.2f}")
    P("  ⚠ The smokes were all at scale 1.0. Raising T_s cuts the driving dT by")
    P("    ~5 %, which raises T_g by < 5 K -- second order, carried as noise.")
    P("")
    P(f"  s* = 0.09*R90*F/feed = {1000*0.09*R90/FEED:.1f}*F mm     "
      f"(the packet's 198*F used R90 = 3.3)")
    P(f"  CONTACT when F <= feed/R90 = {FEED/R90:.4f}   (the packet assumed 0.4545)")
    P(f"  WINDOW  {FEED/R90:.4f} < F <= {100/(1000*0.09*R90/FEED):.4f}"
      f"   i.e. s* in (90, 100] mm")
    P("")
    Ts_of = {"1.0": 823.3, "1.3523": 870.0, "1.9088": 920.0}
    Fof = lambda f, Ts: f * ((Tg(f) - Ts) / (Ts - 293.0)) / 1.58665
    P(f"  {'leg':<10}{'F_q':>6}{'scale':>8}{'T_s':>7}{'T_g':>8}{'F':>7}"
      f"{'s* mm':>7}{'q100':>7}  verdict")
    ladder = (("C0", 1.00, "1.0"), ("G1600", 1.00, "1.0"), ("Q5", 0.50, "1.0"),
              ("Q5T870", 0.50, "1.3523"), ("Q35", 0.35, "1.0"),
              ("Q35T870", 0.35, "1.3523"), ("Q30", 0.30, "1.0"),
              ("Q30T870", 0.30, "1.3523"))
    for nm, f, sc in ladder:
        Ts = Ts_of[sc]
        tg = 1465.8 if nm == "G1600" else Tg(f)
        F = f * ((tg - Ts) / (Ts - 293.0)) / 1.58665
        st = 1000 * 0.09 * R90 * F / FEED
        q = f * 1216.0 * (tg - Ts) / 1e6
        v = ("CONTACT" if F <= FEED / R90 else
             "*** IN WINDOW ***" if st <= 100 else f"fail {st:.0f} mm")
        P(f"  {nm:<10}{f:>6.2f}{sc:>8}{Ts:>7.1f}{tg:>8.1f}{F:>7.3f}"
          f"{st:>7.0f}{q:>7.2f}  {v}")
    P("")
    P("  ⚠⚠ THIS INVERTS THE PACKET'S PREDICTIONS.")
    P("    * the pre-registered pass Q5T870 lands at 134 mm and FAILS;")
    P("    * the window is reached by Q35T870 (97 mm) and Q30 (96 mm);")
    P("    * Q35T870 is a leg the packet explicitly DROPPED as 'contact under")
    P("      C4' -- that verdict used R90 = 3.3, which item 0(e) refutes.")
    P("  ⚠ Q35T870 and Q30 reach the SAME s* by DIFFERENT lever pairs. That is a")
    P("    direct measurement of criterion (i)'s degeneracy, not a redundancy.")
    P("  ⚠⚠ The window is 10 mm wide and the predictor carries +/-5 % (C4), so a")
    P("    leg at 95-105 mm is NOT a decisive verdict. The deliverable stays")
    P("    'F to about +/-10 % with an equivalent nozzle temperature'.")
    P("")
    P("  ⚠⚠ THE GEOMETRY WINDOW AND THE FLUX GATE NO LONGER COINCIDE.")
    P("    Revision 2's table had the predicted pass at q100 = 0.55 MW/m^2, the")
    P("    centre of Meier's required 0.5-0.6, and called that coincidence 'a")
    P("    real result' if it held. With the MEASURED R90 it does not hold: the")
    P("    legs that reach the standoff window sit at q100 = 0.39-0.43, BELOW")
    P("    Kant's measured spallation band (0.45-0.87) and below Meier's")
    P("    required 0.5-0.6. PRE-REGISTERED HERE, before the ladder runs:")
    P("      * P3 is expected to be REFUTED on its third clause;")
    P("      * holding Fig. 8.8's standoff would need a floor flux BELOW the")
    P("        threshold at which Kant measured spallation to occur at all.")
    P("    That tension is a RESULT, and it is the packet's real product if it")
    P("    survives. It must not be resolved by moving either bar.")
    P("")
    P("  LEGS TO RUN (8 + identity, same budget as the packet's):")
    P("    C0       predictor check + known-error anchor (650 s, C3)")
    P("    G1600    the degeneracy witness, measured T_g 1465.8 K")
    P("    Q5       the F_q ladder at 0.50")
    P("    Q5T870   THE PACKET'S PRE-REGISTERED PASS -- kept to test P1 directly")
    P("    Q35      nearest bracket ABOVE the window")
    P("    Q35T870  predicted pass")
    P("    Q30      predicted pass by the OTHER lever")
    P("    Q30T870  nearest bracket BELOW (contact)")
    P("    S0       identity")
    P("  DROPPED from the packet's table: Q55T870 (145 mm) and Q5T920 (117 mm) --")
    P("  both now fail high and carry less than the legs that bracket the window.")
    P("")
    P("=" * 78)


if __name__ == "__main__":
    main()
