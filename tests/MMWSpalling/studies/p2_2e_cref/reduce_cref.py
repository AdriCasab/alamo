#!/usr/bin/env python3
"""P2-2e reductions. Modes:
  CAL     b1 + b2: the level witness, the three ladders, the selection rule
          -> printed; redirect to D_REF.txt (the D_REF_M= line is what every
          later leg reads)
  CHECKS  b3, c1-c4 (+ c12 counters on those legs)
  PROD    c5, c7-c12 on H0/HD and SD
Per-step thermo columns are averaged or maxed, never summed (SR 7); every
volume is x4 (SR 3); every number carries its leg directory.
"""
import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
L = os.path.join(HERE, "legs")
sys.path.insert(0, HERE)
import run_cref as RC                                         # noqa: E402

LEVEL = 843.757


def onset(leg):
    p = os.path.join(L, leg, "run_removal_events.csv")
    if not os.path.exists(p):
        return "MISSING"
    for r in csv.DictReader(open(p)):
        if r["regime"] == "1":
            return (float(r["time"]), float(r["T_surf"]))
    return None


def thermo(leg):
    T = open(os.path.join(L, leg, "run", "thermo.dat")).read().split("\n")
    h = T[0].split()
    rows = [r.split() for r in T[1:] if r.strip()]
    return {n: [float(r[h.index(n)]) for r in rows] for n in h}


def ladder(tag):
    """{(d, sod): onset-or-None} over every rung that exists on disk."""
    out = {}
    for d in sorted(set(RC.LADDER_D) | set(RC.LADDER_DN)):
        for sod in (130, 150):
            leg = f"{tag}{d:.1f}_{sod}"
            if os.path.isdir(os.path.join(L, leg)):
                out[(d, sod)] = onset(leg)
    return out


def interval(lad):
    ds = sorted({d for d, _ in lad})
    fires = lambda d, s: isinstance(lad.get((d, s)), tuple)
    f130 = [d for d in ds if fires(d, 130)]
    n150 = [d for d in ds if not fires(d, 150)]
    d_lo = min(f130) if f130 else None
    d_hi = max(n150) if n150 else None
    bottom_ok = not fires(ds[0], 130)
    top_ok = fires(ds[-1], 150)
    mono = all(fires(d, s) <= fires(d2, s)
               for s in (130, 150) for d, d2 in zip(ds, ds[1:]))
    return ds, d_lo, d_hi, bottom_ok, top_ok, mono


def round_half_up_05(x):
    return math.floor(x / 0.5 + 0.5) * 0.5


def cal():
    P = print
    P("P2-2e (b1, b2) -- the calibration of the constraint depth on Kant 2016")
    P("=" * 92)
    o = onset("K0_h500")
    ok1 = isinstance(o, tuple) and abs(o[1] - LEVEL) <= 2.0
    P(f"  b1  K0_h500 (key 0, GRIMSEL thermal, h 500, T_f 1673.15): onset "
      f"T_surf {o[1]:.3f} K at t = {o[0]:.3f} s   bar {LEVEL:.2f} +- 2   "
      f"{'PASS' if ok1 else 'FAIL'}")
    P("      -> the column deck IS the anchor's configuration: the calibration"
      " sits on the banked level")
    P("")
    res = {}
    for tag, lab in (("C", "GOTTHARD thermal -- THE SELECTION LADDER"),
                     ("G", "GRIMSEL thermal -- reported, never selected from"),
                     ("S", "GOTTHARD + robin_form=surface -- reported, never"
                           " selected from")):
        lad = ladder(tag)
        ds, d_lo, d_hi, bok, tok, mono = interval(lad)
        res[tag] = (d_lo, d_hi, bok, tok, mono, lad)
        P(f"  {tag}-ladder: {lab}")
        P(f"    {'d [mm]':>7s} | {'SOD 130 (spalls in Kant)':>26s} |"
          f" {'SOD 150 (does NOT spall)':>26s}")
        for d in ds:
            c = []
            for s in (130, 150):
                v = lad[(d, s)]
                c.append(f"fires t={v[0]:6.2f}s {v[1]:7.2f}K" if isinstance(v, tuple)
                         else ("no fire in 200 s" if v is None else "MISSING"))
            P(f"    {d:7.1f} | {c[0]:>26s} | {c[1]:>26s}")
        P(f"    d_lo (smallest d where SOD 130 fires)       = {d_lo} mm")
        P(f"    d_hi (largest d where SOD 150 does NOT fire) = {d_hi} mm")
        P(f"    interior: SOD 130 silent at the bottom rung ({ds[0]}) {bok};"
          f" SOD 150 fires at the top rung ({ds[-1]}) {tok};"
          f" fire/no-fire monotone in d {mono}")
        P("")
    d_lo, d_hi, bok, tok, mono, lad = res["C"]
    ok2 = (d_lo is not None and d_hi is not None and d_lo <= d_hi
           and bok and tok and mono)
    dref = round_half_up_05((d_lo + d_hi) / 2.0) if ok2 else None
    P("  THE SELECTION RULE (pre-registered, deterministic, Gotthard ladder)")
    P("    ladder 3.0 .. 12.0 mm in 0.5 mm: SOD 130 FIRED at the bottom rung,")
    P("    so it was extended DOWN by 2 mm to 1.0 mm as the rule requires")
    P("    (1.0 / 1.5 / 2.0 / 2.5 added); the top condition held at 12.0 mm.")
    P(f"    separating interval [{d_lo}, {d_hi}] mm, interior on both ends:"
      f" {bok and tok}")
    if ok2:
        P(f"    d_ref = round_half_up_0.5mm(({d_lo} + {d_hi})/2 ="
          f" {(d_lo + d_hi) / 2:.3f}) = {dref:.1f} mm")
    for tag in ("G", "S"):
        g = res[tag]
        same = (g[0] == d_lo and g[1] == d_hi)
        P(f"    {tag}-ladder interval [{g[0]}, {g[1]}] mm -> "
          f"{'IDENTICAL to the selected one' if same else 'DIFFERS'}"
          + (" (S: within one rung, so the verdict flag does not trigger)"
             if tag == "S" and same else ""))
    P("")
    P("  ⚠⚠ FAR BELOW THE PRE-REGISTRATION (7.0-10.0 mm, offline midpoint 8.5).")
    P("    Even the hottest-property offline row (Grimsel at 550 C, 4.0-5.5 mm)")
    P("    does not reach it, so the falling conductivity is not the whole")
    P("    cause. The offline model fires at a CONSTANT drive dT_c = 550.6 K;")
    P("    the binary fires on the full K_I scan against K_Ic(T) with E(T) and")
    P("    alpha(T), so when a short d forces the face hotter before firing")
    P("    (948 K at 2.5 mm on SOD 150) the drive that suffices is a different")
    P("    one. The guardrail governs: d_ref is chosen ONCE by this rule and is")
    P("    NOT adjusted; a disagreement is a finding about the law.")
    P("")
    P("=" * 92)
    if ok1 and ok2:
        P("PASS (b1, b2)")
        P(f"D_REF_M={dref * 1e-3:.4f}")
        P(f"D_REF_MM={dref:.1f}")
        P(f"D_LO_MM={d_lo}")
        P(f"D_HI_MM={d_hi}")
    else:
        P("FAIL (b1/b2) -- no D_REF_M is written, so no later leg can run")




def counters(leg):
    """c12: cref_n_clamp / cref_n_neg as a fraction of reads, ledger max."""
    th = thermo(leg)
    out = {}
    if "cref_n_vert" in th:
        reads = sum(th["cref_n_vert"]) + sum(th["cref_n_lat"])
        # ⚠ a RATIO of two per-step columns sampled at the same cadence: the
        # thermo plot_int factor cancels, so summing both sides is legitimate
        # here (SR 7 forbids reading one summed column as a total, not this).
        out["reads"] = reads
        out["clamp"] = sum(th["cref_n_clamp"]) / reads if reads else 0.0
        out["neg"] = sum(th["cref_n_neg"]) / reads if reads else 0.0
        out["vert_max"] = max(th["cref_n_vert"])
    if "ledger_err_solid" in th:
        out["ledger"] = max(abs(x) for x in th["ledger_err_solid"])
    return out


def n_fire(leg):
    p = os.path.join(L, leg, "run_removal_events.csv")
    return sum(1 for r in csv.DictReader(open(p)) if r["regime"] == "1") \
        if os.path.exists(p) else 0


def checks():
    P = print
    dref = float([l for l in open(os.path.join(HERE, "D_REF.txt"))
                  if l.startswith("D_REF_MM=")][0].split("=")[1])
    P(f"P2-2e (b3, c1-c4, c12) -- the checks at d_ref = {dref:.1f} mm")
    P("=" * 92)
    # ---- b3
    g = n_fire("GB")
    c = counters("GB")
    ok = (g == 0) and c.get("vert_max", 0) > 0
    P(f"  b3  GB (uniform 820 K slab, no net flux, d_ref): {g} regime-1 firings"
      f" (bar 0); reference READ {c.get('vert_max', 0):.0f} times in the max"
      f" step  -> {'PASS' if ok else 'FAIL'}")
    P("      ⚠ weak evidence by the unit test's own account: the key-0 path")
    P("      (T_init = 820 K) also reads 0 here. What it proves is that the")
    P("      constraint reference REACHES the criterion and cancels a uniform state.")
    P("")
    # ---- c1
    P("  c1  the LEVEL at Kant's fast standoffs (Gotthard thermal, d_ref)")
    P("      bar: within +-15 K of 843.76 -- REPORT, a consistency check whose")
    P("      premise is that T(z + d_ref) is still ~293 K at onset")
    for sod in (50, 80, 110):
        o = onset(f"V{sod}")
        k0 = None
        P(f"      V{sod:<3d} h {RC.SOD[sod][0]:4d}, T_f {RC.SOD[sod][1]:.2f}:"
          + (f" onset {o[1]:8.2f} K at {o[0]:6.3f} s   shift {o[1]-LEVEL:+7.2f} K"
             f"   {'within' if abs(o[1]-LEVEL) <= 15 else 'OUTSIDE'} +-15"
             if isinstance(o, tuple) else " NO FIRE in 200 s"))
    P("")
    # ---- c2
    P("  c2  THE OLD ANCHORS (Grimsel thermal, d_ref) -- h150 is the discriminator")
    for an in ("h150", "h500"):
        o = onset(f"V_{an}")
        tm = thermo(f"V_{an}")["time"][-1]
        P(f"      V_{an}: " + (f"FIRES at {o[1]:.2f} K, t = {o[0]:.3f} s"
                               if isinstance(o, tuple)
                               else f"NO FIRE to t = {tm:.1f} s"))
    o150 = onset("V_h150")
    P(f"      h150 never fires in 200 s: {'YES' if o150 is None else 'NO'}"
      f"  (offline: never, at every property set)")
    P("")
    # ---- c3
    P("  c3  MESH at SOD 80 (Gotthard thermal, d_ref, dt 0.0125 on all)")
    on = {}
    for tag, dz in (("2", 2.0), ("1", 1.0), ("05", 0.5)):
        on[tag] = onset(f"M80_{tag}")
        v = on[tag]
        P(f"      M80_{tag:<3s} dz {dz:3.1f} mm: " + (f"onset {v[1]:8.3f} K at {v[0]:6.3f} s"
                                                 if isinstance(v, tuple) else "NO FIRE"))
    vals = [v[1] for v in on.values() if isinstance(v, tuple)]
    if len(vals) == 3:
        rise = sum(vals) / 3 - 293.15
        spr = max(vals) - min(vals)
        P(f"      spread {spr:.3f} K over a rise of {rise:.1f} K = "
          f"{100 * spr / rise:.3f} %   bar <= 1 %   "
          f"{'PASS' if spr / rise <= 0.01 else 'MISS'}")
    P("")
    # ---- c4
    P("  c4  THE BEAM COLUMNS (s1, 1 mm, 5 MW/m2, P2-2d criterion keys, 30 s)")
    for a, lab in ((2000, "surface-absorbed, 1/alpha 0.5 mm"),
                   (100, "volumetric, 1/alpha 10 mm")):
        o0, od = onset(f"B{a}_0"), onset(f"B{a}_d")
        f0 = f"{o0[1]:8.2f} K at {o0[0]:6.3f} s" if isinstance(o0, tuple) else "no fire"
        fd = f"{od[1]:8.2f} K at {od[0]:6.3f} s" if isinstance(od, tuple) else "no fire"
        P(f"      B{a:<5d} ({lab})")
        P(f"         key 0     : {f0}")
        P(f"         key d_ref : {fd}")
        if isinstance(o0, tuple) and isinstance(od, tuple):
            P(f"         shift {od[1] - o0[1]:+.2f} K")
        n0, nd = n_fire(f"B{a}_0"), n_fire(f"B{a}_d")
        P(f"         regime-1 firings in 30 s: key 0 {n0}, key d_ref {nd}")
    P("")
    # ---- c12 on every check leg
    P("  c12 COUNTERS on every check leg (clamp and negative-bracket fractions of"
      " reads; ledger max)")
    for leg in ("V50", "V80", "V110", "V_h500", "V_h150", "M80_2", "M80_1",
                "M80_05", "GB", "B2000_d", "B100_d"):
        c = counters(leg)
        if "reads" in c:
            P(f"      {leg:8s} reads {c['reads']:12.0f}  clamp {100*c['clamp']:7.4f} %"
              f"  neg {100*c['neg']:7.4f} %"
              + (f"  ledger {c['ledger']:.2e}" if "ledger" in c else ""))
    P("=" * 92)


# ============================================================== PROD (c5, c7-c12)
STUD = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(STUD, "p2_2c_meier"))
sys.path.insert(0, os.path.join(STUD, "p2_4_body"))
sys.path.insert(0, os.path.join(STUD, "p2_2_props"))
S0_DIR = os.path.join(STUD, "p2_4_body", "legs", "S0")
WIN = (150.0, 350.0)


def ev_rows(d):
    return list(csv.DictReader(open(os.path.join(d, "run_removal_events.csv"))))


def in_win(r, lo=WIN[0], hi=WIN[1]):
    t = float(r["time"])
    return lo <= t <= hi


def windowed(d, label):
    """The packet's WINDOWED rows, computed identically for S0 and SD with
    reduce_price's own charged_volume / firing_T / col_radius (Meier, 2 mm)."""
    import reduce_price as RPm
    ev = ev_rows(d)
    w = [r for r in ev if in_win(r)]
    cv = RPm.charged_volume(w)
    share = cv["wall"]["vol_nvoided"] / cv["charged"] if cv["charged"] else None
    cvall = RPm.charged_volume(ev)
    share_all = (cvall["wall"]["vol_nvoided"] / cvall["charged"]
                 if cvall["charged"] else None)
    near = [r for r in w if r["regime"] == "1"
            and RPm.col_radius(int(float(r["col_i"])), int(float(r["col_j"]))) < 0.020
            and float(r["T_surf"] or 0) > 0]
    ts = sorted(float(r["T_surf"]) for r in near)
    med = ts[len(ts) // 2] if ts else None
    allts = sorted(float(r["T_surf"]) for r in ev
                   if r["regime"] == "1" and float(r["T_surf"] or 0) > 0)
    return dict(label=label, share=share, share_all=share_all,
                wall_cells=cv["wall"]["cells"], floor_cells=cv["floor"]["cells"],
                Tfire=med, nTfire=len(ts),
                Tfire_all=allts[len(allts) // 2] if allts else None)


def rate(d, lo=WIN[0], hi=WIN[1]):
    T = open(os.path.join(d, "run", "thermo.dat")).read().split("\n")
    h = T[0].split()
    pts = [(float(r.split()[h.index("time")]), float(r.split()[h.index("nozzle_z")]))
           for r in T[1:] if r.strip()]
    pts = [(a, b) for a, b in pts if lo <= a <= hi and b > 0]
    return (pts[0][1] - pts[-1][1]) / (pts[-1][0] - pts[0][0]) * 3600.0


def c7(d, d_ref_m, n_want=20):
    """At fired LATERAL faces: the profile model's T_ref = T0 + dT_s e^{-d/ell}
    vs the FIELD read d behind the face in the plotfile written just before
    the firing. Ratio of the two brackets (T_face - T_ref)."""
    import fabcmp
    T0 = 293.15
    run = os.path.join(d, "run")
    pfs = sorted(x for x in os.listdir(run) if x.endswith("cell"))
    dt = 0.016
    pf_t = {int(x[:5]) * dt: x for x in pfs}         # plot_int is in steps
    ev = [r for r in ev_rows(d) if r["regime"] == "6" and in_win(r)]
    # events fired within 2 s AFTER a plotfile time: the field is <= 2 s stale
    picks = []
    for tp in sorted(t for t in pf_t if WIN[0] <= t <= WIN[1]):
        cand = [r for r in ev if 0.0 < float(r["time"]) - tp <= 2.0
                and float(r["T_face_lat"] or 0) > T0 and float(r["ell_lat"] or 0) > 0]
        step = max(1, len(cand) // 4)
        picks += [(tp, r) for r in cand[::step][:4]]
    picks = picks[:n_want]
    cache = {}

    def fld(tp):
        if tp in cache:
            return cache[tp]
        lvl = os.path.join(run, pf_t[tp], "Level_0")
        files = sorted(f for f in os.listdir(lvl) if f.startswith("Cell_D_"))
        first = fabcmp.read_fabs_boxed(os.path.join(lvl, files[0]))
        names = fabcmp.var_names(os.path.join(run, pf_t[tp], "Header"), first[0]["ncomp"])
        it, ir = names.index("Temp"), names.index("removed")
        T, R = {}, {}
        for fn in files:
            for fab in fabcmp.read_fabs_boxed(os.path.join(lvl, fn)):
                lo, hi = fab["lo"], fab["hi"]
                nx, ny = hi[0] - lo[0] + 1, hi[1] - lo[1] + 1
                for o, (tv, rv) in enumerate(zip(fab["comp"][it], fab["comp"][ir])):
                    key = (lo[0] + o % nx, lo[1] + (o // nx) % ny, lo[2] + o // (nx * ny))
                    T[key] = tv; R[key] = rv
        cache[tp] = (T, R)
        return cache[tp]

    rows = []
    dx = 0.002
    for tp, r in picks:
        i, j = int(float(r["col_i"])), int(float(r["col_j"]))
        k = int(float(r["k_top"])) - int(round(float(r["riser_depth"] or 0) / dx))
        T, R = fld(tp)
        if R.get((i, j, k), 1.0) > 0.5:
            continue                                  # already void at tp
        # the void lateral neighbour gives the face; inward is opposite
        inward = None
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            if R.get((i - di, j - dj, k), 0.0) > 0.5 and (i + di, j + dj, k) in T:
                inward = (di, dj); break
        if inward is None:
            continue
        # field T at d behind the face: cell centres sit at dx/2, 3dx/2, ...
        mf = d_ref_m / dx - 0.5
        m0 = int(math.floor(mf)); w = mf - m0
        a = T.get((i + inward[0] * m0, j + inward[1] * m0, k))
        b = T.get((i + inward[0] * (m0 + 1), j + inward[1] * (m0 + 1), k))
        if a is None or b is None:
            continue
        Tf_field = (1 - w) * a + w * b
        Tface = float(r["T_face_lat"]); ell = float(r["ell_lat"])
        Tref_model = T0 + (Tface - T0) * math.exp(-d_ref_m / ell)
        den = Tface - Tf_field
        if abs(den) < 1e-9:
            continue
        rows.append(((Tface - Tref_model) / den, Tface, Tref_model, Tf_field, ell, tp))
    return rows


def prod():
    P = print
    dref = float([l for l in open(os.path.join(HERE, "D_REF.txt"))
                  if l.startswith("D_REF_M=")][0].split("=")[1])
    P(f"P2-2e (c5, c7-c12) -- the production legs at d_ref = {dref*1e3:.1f} mm")
    P("=" * 94)
    # ---- c5 Hoeser invariance
    for leg in ("H0", "HD"):
        if not os.path.isdir(os.path.join(L, leg)):
            P(f"  c5  {leg}: not run"); return
    out = {}
    for leg in ("H0", "HD"):
        ev = ev_rows(os.path.join(L, leg))
        nv = sum(int(float(r["n_voided"] or 0)) for r in ev)
        ts = sorted(float(r["T_surf"]) for r in ev
                    if r["regime"] == "1" and float(r["T_surf"] or 0) > 0)
        out[leg] = (nv * 1e-9 * 4 * 1e6, ts[len(ts) // 2] if ts else None, len(ev))
    v0, vd = out["H0"][0], out["HD"][0]
    t0, td = out["H0"][1], out["HD"][1]
    P(f"  c5  HOESER at 7 MW/m2 (G05 on a 200 mm block, 30 s, 1 mm)")
    P(f"      H0 key 0     : V = {v0:8.3f} cm3 (x4), firing T median {t0:.2f} K")
    P(f"      HD key d_ref : V = {vd:8.3f} cm3 (x4), firing T median {td:.2f} K")
    P(f"      volume {100*(vd-v0)/v0:+.2f} % (bar within 10 %)   firing T "
      f"{td-t0:+.2f} K (bar within 10 K)   "
      f"{'PASS' if abs(vd-v0)/v0 <= 0.10 and abs(td-t0) <= 10 else 'MISS'}")
    P("")
    if not os.path.isdir(os.path.join(L, "SD")) or \
       not os.path.exists(os.path.join(L, "SD", ".done")):
        P("  SD not finished -- c7-c11 pending"); return
    sd = os.path.join(L, "SD")
    # ---- c8, c10 windowed, for S0 AND SD with the same code
    w0, wd = windowed(S0_DIR, "S0"), windowed(sd, "SD")
    P("  ⚠ THE PACKET'S S0 REFERENCES ARE WHOLE-RUN NUMBERS: 12.14 % and 851.5 K come")
    P("    from reduce_p2_4.one(), which reports the share over 0-500 s ('whole-run")
    P("    share; the windowed split needs the event log by time') and the median")
    P("    firing T over ALL radii and the whole run. c8 and c10 are WINDOWED bars,")
    P("    so S0 is recomputed below in exactly each bar's window with the same code.")
    P("")
    P(f"  c8  WALL SHARE (regime-6 charged / charged; pad clearing excluded)")
    for w in (w0, wd):
        P(f"      {w['label']}: 150-350 s {100*w['share']:6.2f} %  "
          f"({w['wall_cells']} wall cells / {w['wall_cells']+w['floor_cells']}) | "
          f"whole run {100*w['share_all']:6.2f} %")
    P(f"      bar <= 4 % : SD {100*wd['share']:.2f} %  -> "
      f"{'PASS' if wd['share'] <= 0.04 else 'MISS'}   (like-for-like S0 "
      f"{100*w0['share']:.2f} %; the packet quoted S0's whole-run 12.1 %)")
    P("")
    r0, rd = rate(S0_DIR), rate(sd)
    P(f"  c9  RATE 150-350 s (nozzle_z descent): S0 {r0:.3f}  SD {rd:.3f} m/h"
      f"   bar 1.6-2.1  -> {'PASS' if 1.6 <= rd <= 2.1 else 'MISS'}")
    import reduce_p2_4 as R4
    R4.LEGS_DIR = L
    import io, contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        res = R4.one("SD")
    deq = res.get("Deq_end") if res else None
    P(f"      Oeq at 10 mm (mask, end of run, reduce_p2_4.one): SD "
      f"{deq if deq is None else round(deq, 2)} mm vs S0 108.2   bar <= S0  -> "
      f"{'PASS' if deq is not None and deq <= 108.2 else 'MISS'}")
    P("")
    P(f"  c10 FIRING T, regime 1, r < 20 mm, 150-350 s (median)")
    P(f"      S0 {w0['Tfire']:.2f} K (n {w0['nTfire']})   SD {wd['Tfire']:.2f} K "
      f"(n {wd['nTfire']})   bar 860-910  -> "
      f"{'PASS' if 860 <= wd['Tfire'] <= 910 else 'MISS'}")
    P(f"      (whole-run all-radius medians, the packet's 851.5 K convention: "
      f"S0 {w0['Tfire_all']:.2f}, SD {wd['Tfire_all']:.2f})")
    P("")
    # ---- c11
    import s0_census as SC
    P("  c11 FROZEN-COLUMN CENSUS (P1's definition), SD vs S0")
    l0, t0m = SC.last_fired(S0_DIR); ld, tdm = SC.last_fired(sd)
    dx = SC.dx_of(sd)
    worst = 0.0
    for t in SC.TIMES:
        a = SC.census(l0, dx, min(t, t0m)); b = SC.census(ld, dx, min(t, tdm))
        pa, pb = 100*a[0]/max(a[1],1), 100*b[0]/max(b[1],1)
        worst = max(worst, pb - pa)
        P(f"      t {t:5.0f} s: S0 {a[0]:3d}/{a[1]} = {pa:5.2f} %   SD {b[0]:3d}/{b[1]}"
          f" = {pb:5.2f} %   delta {pb-pa:+.2f} points")
    P(f"      bar: SD <= S0 + 5 points at every time  -> "
      f"{'PASS' if worst <= 5.0 else 'MISS'} (worst {worst:+.2f})")
    P("")
    # ---- c7
    rows = c7(sd, dref)
    P(f"  c7  LATERAL WITNESS on SD: {len(rows)} fired lateral faces sampled within 2 s")
    P(f"      after a plotfile (cadence 50 s); ratio = (T_face - T_ref,model) /")
    P(f"      (T_face - T_ref,field), field read {dref*1e3:.1f} mm behind the face")
    if rows:
        rs = sorted(x[0] for x in rows)
        med = rs[len(rs) // 2]
        for x in rows[:8]:
            P(f"        t_pf {x[5]:5.0f}  T_face {x[1]:7.2f}  model {x[2]:7.2f}  field"
              f" {x[3]:7.2f}  ell {x[4]*1e3:5.2f} mm  ratio {x[0]:6.3f}")
        P(f"      median ratio {med:.3f} over {len(rows)}   bar 0.7-1.3  -> "
          f"{'PASS' if 0.7 <= med <= 1.3 else 'MISS'}")
    else:
        P("      NO usable lateral firings in the window -- c7 unscoreable on SD")
    P("")
    # ---- c12
    P("  c12 COUNTERS and LEDGER on the production legs")
    for leg in ("HD", "SD", "H0"):
        c = counters(leg)
        P(f"      {leg:3s} " + (f"reads {c['reads']:.0f} clamp {100*c['clamp']:.4f} %"
                               f" neg {100*c['neg']:.4f} %  " if "reads" in c else "key 0  ")
          + (f"ledger_err_solid max {c['ledger']:.3e}" if "ledger" in c else ""))
    P("=" * 94)


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "CAL"
    if mode == "CAL":
        cal()
    elif mode == "CHECKS":
        checks()
    elif mode == "PROD":
        prod()
    else:
        sys.exit(f"mode {mode}: not yet")


if __name__ == "__main__":
    main()
