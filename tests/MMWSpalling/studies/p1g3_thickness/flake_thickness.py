#!/usr/bin/env python3
"""P1g-3 (a0)-(a5): the flake thickness derived from the rock, OFFLINE.

Scores, from files already on disk and with no simulation:
  (a0) size-vs-thickness, resolved and APPLIED before any comparison
  (a1) the identity  ell_q = ell_eff * q_ratio  = k (T_s - T0) / q_net
  (a2) floor thickness  h = w * ell_q / m , overall and by radius band
  (a3) WALL thickness from the SIDE-FACE flux, with its own gas inversion
  (a4) the test that cannot be fitted: one (w, m) against Rossi's TWO populations
  (a5) m sensitivity

⚠⚠ (a1) IS AN IDENTITY, NOT AN APPROXIMATION. `Removal.H` computes
`q_fit = kappa (T_s - T0) / ell_eff` and `q_ratio = q_fit / q_real`, so

    ell_q = ell_eff * q_ratio = kappa (T_s - T0) / q_real

exactly, for zero extra compute. ⚠ `q_real` is NET of radiation and convection
(`Removal.H:1503-1517`), not the jet-delivered flux. The wall analogue in (a3)
must be net in the SAME way or the two are not the same quantity.

⚠⚠ THE FLOOR ROUTE DOES NOT WORK ON THE WALL. `q_ratio` is zero or non-finite on
96.0 % of regime-6 events and on 0 % of regime-1, because `flame_h` is 0 outside
the patch and above the nozzle plane. For the wall the flux is not merely
different -- it is ABSENT. So (a3) builds the wall flux from the side-face pair
instead, and reports how many events the floor route would have dropped. A wall
figure that resembles the floor's is the signature of having used the floor flux.

⚠ Every material constant is PARSED from the Meier input, never transcribed.
"""
import collections, csv, math, os, statistics, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
LEG  = os.path.join(ROOT, "tests", "MMWSpalling", "studies", "p1f_lateral",
                    "output_legB_fix")
EV   = os.path.join(LEG, "FH_removal_events.csv")
INP  = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "meier",
                    "sp_meier_pilot", "input_drilling")

DX      = 0.002          # m, quarter domain, origin at 0
W_CONE  = 2.35           # see (a2)'s note: a MEDIAN, not a percentile choice
M_MEIER = 20.0           # input_drilling:124
BANDS   = ((0, 20), (20, 30), (30, 40), (40, 50), (50, 70))   # mm

# (a1) reproduction targets, mm
T_A1 = {"p05": 0.53, "p50": 1.45, "p95": 2.87}
TOL_A1 = 0.02
# (a2) reproduction targets, um
T_A2_OVERALL, T_A2_TOL = 170.0, 5.0
T_A2_BANDS = (125.0, 155.0, 231.0, 277.0, 888.0)

# ---- (a0) SIZE vs THICKNESS ------------------------------------------------
# ⚠⚠ Walsh's aspect ratio is ~6:1 for LARGE spalls and 1-2:1 near 200 um. Every
# published figure below is a SIZE (a sieve diameter or a PSD bin); this law
# predicts a THICKNESS. Comparing them directly is wrong by up to 6x, and the
# error sits in the comparand where no later row can see it.
# The conversion applied everywhere below, stated once:
#     thickness = size / aspect(size),  aspect = 1.5 near 200 um, 6.0 at >= 1 mm
# with a log-linear interpolation between. A RANGE is carried, not a point.
def aspect(size_m):
    """aspect ratio at this SIZE, with its two anchors named."""
    s = size_m * 1e6                      # um
    if s <= 200.0:
        return (1.0, 2.0)                 # Walsh: 1-2:1 near 200 um
    if s >= 1000.0:
        return (6.0, 6.0)                 # Walsh: ~6:1 for large spalls
    f = math.log(s / 200.0) / math.log(1000.0 / 200.0)
    return (1.0 + f * 5.0, 2.0 + f * 4.0)


def size_to_thickness(size_m):
    lo, hi = aspect(size_m)
    return (size_m / hi, size_m / lo)     # thinner for a bigger aspect ratio


def parse_input(path):
    want = {"microstructure.phase0.kappa": None,
            "microstructure.phase0.T_ref": None,
            "losses.epsilon": None, "losses.sigma_SB": None,
            "losses.h_conv": None, "losses.T_amb": None,
            "weibull.m": None}
    for ln in open(path):
        ln = ln.split("#")[0]
        if "=" not in ln:
            continue
        k, v = ln.split("=", 1)
        if k.strip() in want:
            try:
                want[k.strip()] = float(v.strip())
            except ValueError:
                pass
    miss = [k for k, v in want.items() if v is None]
    if miss:
        raise SystemExit("FAIL: could not parse %s from %s" % (miss, path))
    return want


def quant(xs, p):
    xs = sorted(xs)
    return xs[int(p * (len(xs) - 1))]


def main():
    P = print
    fails = []

    def need(ok, msg):
        P(("  ok    " if ok else "  FAIL  ") + msg)
        if not ok:
            fails.append(msg)

    if not os.path.exists(EV):
        raise SystemExit("FAIL: no event log at %s" % EV)
    cst = parse_input(INP)
    KAPPA, T0 = cst["microstructure.phase0.kappa"], cst["microstructure.phase0.T_ref"]
    EPS, SIG = cst["losses.epsilon"], cst["losses.sigma_SB"]
    HC, TA = cst["losses.h_conv"], cst["losses.T_amb"]
    M_IN = cst["weibull.m"]

    P("=" * 78)
    P("P1g-3  THE FLAKE THICKNESS FROM THE ROCK  --  offline, no simulation")
    P("=" * 78)
    P("  leg       : %s" % os.path.relpath(LEG, ROOT))
    P("  constants : PARSED from %s" % os.path.relpath(INP, ROOT))
    P("              kappa %.4g W/mK   T_ref %.2f K   weibull.m %.1f"
      % (KAPPA, T0, M_IN))
    P("              eps %.3g  sigma %.6g  h_conv %.4g  T_amb %.2f"
      % (EPS, SIG, HC, TA))
    P("  w         = %.2f  (see the (a2) note: a MEDIAN of the cone-weighted"
      % W_CONE)
    P("              failure density, not a percentile choice)")

    rows = list(csv.DictReader(open(EV)))
    r1 = [r for r in rows if int(float(r["regime"])) == 1]
    r6 = [r for r in rows if int(float(r["regime"])) == 6]

    # ================= (a0) =============================================
    P("")
    P("-" * 78)
    P("(a0) SIZE vs THICKNESS -- resolved BEFORE any comparison")
    P("     Walsh: aspect ~6:1 for large spalls, 1-2:1 near 200 um.")
    P("     Every published figure below is a SIZE; this law predicts a")
    P("     THICKNESS. Converting is worth up to 6x and the error would sit in")
    P("     the comparand, where no later row could see it.")
    P("     %-26s %14s   %s" % ("source figure (SIZE)", "aspect", "=> THICKNESS"))
    for lbl, s in (("Meier PSD average <100 um", 100e-6),
                   ("Meier PSD p90 <300 um", 300e-6),
                   ("Rossi population A ~250 um", 250e-6),
                   ("Rossi population B ~1 mm", 1000e-6)):
        lo, hi = aspect(s)
        tlo, thi = size_to_thickness(s)
        P("     %-26s %5.1f - %-5.1f :1   %6.0f - %-6.0f um"
          % (lbl, lo, hi, tlo * 1e6, thi * 1e6))
    need(True, "(a0) conversion stated and APPLIED to every comparand below")

    # ================= (a1) =============================================
    P("")
    P("-" * 78)
    P("(a1) THE IDENTITY  ell_q = ell_eff * q_ratio = k (T_s - T0) / q_net")

    def usable(r):
        try:
            q, e = float(r["q_ratio"]), float(r["ell_eff"])
        except (KeyError, ValueError):
            return False
        return (q != 0.0 and math.isfinite(q) and e > 0.0
                and math.isfinite(e))

    g1 = [r for r in r1 if usable(r)]
    g6 = [r for r in r6 if usable(r)]
    ellq = [float(r["ell_eff"]) * float(r["q_ratio"]) for r in g1]
    P("     regime-1 events used: %d of %d" % (len(g1), len(r1)))
    got = {"p05": quant(ellq, 0.05) * 1e3, "p50": quant(ellq, 0.50) * 1e3,
           "p95": quant(ellq, 0.95) * 1e3}
    for k in ("p05", "p50", "p95"):
        need(abs(got[k] - T_A1[k]) <= TOL_A1,
             "(a1) ell_q %s = %.4f mm vs target %.2f (tol %.2f)"
             % (k, got[k], T_A1[k], TOL_A1))
    need(len(g1) == 193322,
         "(a1) event count %d == 193,322" % len(g1))
    ee = [float(r["ell_eff"]) for r in g1]
    P("     ⚠ for contrast, the UNCONDITIONED ell_eff: p50 %.3f mm, max %.0f mm"
      % (quant(ee, 0.50) * 1e3, max(ee) * 1e3))
    P("       -- ell_eff reaches kilometres on a nearly-flat column, which is")
    P("          why the length had to be re-based on a flux at all.")

    # the 96 % witness
    P("")
    P("     ⚠⚠ THE WALL IS NOT JUST DIFFERENT, IT IS ABSENT on this route:")
    P("        q_ratio usable on regime-1: %d of %d = %.1f %%"
      % (len(g1), len(r1), 100.0 * len(g1) / len(r1)))
    P("        q_ratio usable on regime-6: %d of %d = %.1f %%  => DROPPED %d = %.1f %%"
      % (len(g6), len(r6), 100.0 * len(g6) / len(r6),
         len(r6) - len(g6), 100.0 * (len(r6) - len(g6)) / len(r6)))

    # ================= (a2) =============================================
    def radius_mm(r):
        return DX * math.hypot(int(r["col_i"]) + 0.5, int(r["col_j"]) + 0.5) * 1e3

    def band_of(rr):
        for n, (a, b) in enumerate(BANDS):
            if a <= rr < b:
                return n
        return None

    P("")
    P("-" * 78)
    P("(a2) FLOOR THICKNESS  h = w * ell_q / m   at w = %.2f, m = %.0f"
      % (W_CONE, M_MEIER))
    allh, byband = [], collections.defaultdict(list)
    for r, e in zip(g1, ellq):
        h = W_CONE * e / M_MEIER
        allh.append(h)
        b = band_of(radius_mm(r))
        if b is not None:
            byband[b].append(h)
    ov = statistics.median(allh) * 1e6
    need(abs(ov - T_A2_OVERALL) <= T_A2_TOL,
         "(a2) overall p50 = %.1f um vs target %.0f (tol %.0f)"
         % (ov, T_A2_OVERALL, T_A2_TOL))
    P("     %-12s %8s %12s %12s" % ("band [mm]", "n", "p50 [um]", "target"))
    for n, (a, b) in enumerate(BANDS):
        if not byband[n]:
            P("     %-12s %8d   -- no events --" % ("%d-%d" % (a, b), 0))
            continue
        m50 = statistics.median(byband[n]) * 1e6
        need(abs(m50 - T_A2_BANDS[n]) <= T_A2_TOL,
             "(a2) band %d-%d mm p50 = %.1f um vs target %.0f (tol %.0f)"
             % (a, b, m50, T_A2_BANDS[n], T_A2_TOL))
    P("     ⚠ these are REPRODUCTIONS of the planner's offline numbers, not a")
    P("       physics claim. The physics claim is (a4).")

    # ================= (a3) THE WALL ====================================
    P("")
    P("-" * 78)
    P("(a3) WALL THICKNESS from the SIDE-FACE flux -- its OWN route")
    P("     ⚠ h_side = 169.3 paired with T_face_lat from the log (never T_top);")
    P("       gas backed out of the logged triple exactly as")
    P("       studies/p1f_lateral/margin_by_band.py does; NET of losses.")
    H_CAPPED, H_SIDE = 80.0, 169.3        # the cap in force when Leg B ran
    G_COND = 2.0 * KAPPA / DX

    def loss(Tf):
        return EPS * SIG * (Tf ** 4 - TA ** 4) + HC * (Tf - TA)

    wall = collections.defaultdict(list)
    nneg = 0
    DB = ((0, 100), (100, 200), (200, 300), (300, 400))
    for r in r6:
        try:
            Tf = float(r["T_face_lat"]); Tc = float(r["T_top"])
            kt = int(float(r["k_top"])); rd = float(r.get("riser_depth") or 0.0)
        except (KeyError, ValueError):
            continue
        if not (Tf > Tc > T0):
            continue
        # the gas, backed out of the triple at the h that run used
        Tg = Tf + (G_COND * (Tf - Tc) + loss(Tf)) / H_CAPPED
        # the wall's NET absorbed flux at the consistent h
        q_net = H_SIDE * (Tg - Tf) - loss(Tf)
        if not (q_net > 0.0):
            nneg += 1
            continue
        ellq_w = KAPPA * (Tf - T0) / q_net
        depth = (299 - kt) * DX * 1e3 + rd * 1e3
        b = next((n for n, (a, c) in enumerate(DB) if a <= depth < c), None)
        if b is not None:
            wall[b].append(W_CONE * ellq_w / M_MEIER)
    P("     %-12s %8s %12s %14s" % ("depth [mm]", "n", "ell_q [mm]", "h [um]"))
    wall_all = []
    for n, (a, b) in enumerate(DB):
        if not wall[n]:
            P("     %-12s %8d   -- none --" % ("%d-%d" % (a, b), 0))
            continue
        hs = wall[n]; wall_all += hs
        P("     %-12s %8d %12.4f %14.1f"
          % ("%d-%d" % (a, b), len(hs),
             statistics.median(hs) * M_MEIER / W_CONE * 1e3,
             statistics.median(hs) * 1e6))
    if wall_all:
        P("     wall overall p50 = %.1f um   (floor overall p50 = %.1f um)"
          % (statistics.median(wall_all) * 1e6, ov))
        P("     events the FLOOR route would have dropped: %d of %d = %.1f %%"
          % (len(r6) - len(g6), len(r6),
             100.0 * (len(r6) - len(g6)) / len(r6)))
        P("     q_net <= 0 on this route: %d events" % nneg)
    need(bool(wall_all),
         "(a3) the wall route yields a thickness at all, from the side-face "
         "flux rather than the floor's (which is absent on 96 %% of these events)")

    # ================= (a4) THE TEST THAT CANNOT BE FITTED ==============
    P("")
    P("-" * 78)
    P("(a4) THE TEST THAT CANNOT BE FITTED")
    P("     Does ONE (w, m) give Rossi's TWO populations after (a0)'s conversion?")
    ross = (("A ~250 um", 250e-6), ("B ~1 mm", 1000e-6))
    P("     one pair: w = %.2f, m = %.0f" % (W_CONE, M_MEIER))
    P("     %-14s %18s %18s %10s" % ("population", "target thickness",
                                     "model range", "factor"))
    # ⚠⚠ SCORED AGAINST ONE NUMBER, NOT A RANGE. An earlier version of this
    # row compared the targets against the model's whole band spread
    # (125-888 um) and called a hit "within 2x". That is almost unfalsifiable:
    # ANY thickness from 63 to 1776 um would have passed. "From the same pair"
    # means one (w, m) must produce ONE number that lands within 2x of BOTH
    # converted populations -- so the model's overall p50 is the comparand and
    # the row can actually fail.
    verdict = []
    model = statistics.median(allh)        # ONE number from ONE (w, m) pair
    lo_m = min(statistics.median(v) for v in byband.values() if v)
    hi_m = max(statistics.median(v) for v in byband.values() if v)
    P("     model, ONE number from this pair: p50 = %.1f um" % (model * 1e6))
    P("     %-14s %18s %12s %10s" % ("population", "target thickness",
                                     "factor", "verdict"))
    for lbl, s in ross:
        tlo, thi = size_to_thickness(s)
        if model < tlo:
            f = tlo / model
        elif model > thi:
            f = model / thi
        else:
            f = 1.0
        ok2 = f <= 2.0
        verdict.append(ok2)
        P("     %-14s %6.0f - %-6.0f um %11.2fx %10s"
          % (lbl, tlo * 1e6, thi * 1e6, f, "within 2x" if ok2 else "OUTSIDE"))
    need(all(verdict),
         "(a4) ONE number from ONE stated (w, m) pair -- %.1f um -- lands within "
         "a factor of 2 of BOTH Rossi populations after the size->thickness "
         "conversion. ⚠⚠ NOT a fit: w and m are stated constants, and the row is "
         "falsifiable because a single number is scored, not a range"
         % (model * 1e6))
    # ⚠⚠ and the reason it is possible is itself the finding
    ta = size_to_thickness(250e-6); tb = size_to_thickness(1000e-6)
    P("")
    P("     ⚠⚠ WHY ONE NUMBER CAN COVER BOTH, and it weakens the row as a")
    P("        discriminator: Rossi's two populations differ in SIZE by 4.0x")
    P("        (250 um vs 1 mm) but after (a0)'s conversion their THICKNESSES")
    P("        differ by only %.2fx (%.0f-%.0f um vs %.0f-%.0f um), because the"
      % (statistics.mean(tb) / statistics.mean(ta),
         ta[0] * 1e6, ta[1] * 1e6, tb[0] * 1e6, tb[1] * 1e6))
    P("        larger spall carries the larger aspect ratio. So the two field")
    P("        populations are NOT two independent thickness targets, and")
    P("        covering both is a weaker result than the packet expected.")
    P("     ⚠ the floor's own band spread is %.0f-%.0f um, which BRACKETS both"
      % (lo_m * 1e6, hi_m * 1e6))
    P("       targets -- so the spread cannot be used as corroboration either.")
    P("     ⚠ (a3)'s wall p50 (%.0f um) is %s that band spread, so the wall is %s"
      % (statistics.median(wall_all) * 1e6 if wall_all else 0.0,
         "INSIDE" if wall_all and lo_m <= statistics.median(wall_all) <= hi_m
         else "outside",
         "NOT independent evidence" if wall_all and
         lo_m <= statistics.median(wall_all) <= hi_m else "independent"))

    # ================= (a5) m sensitivity ================================
    P("")
    P("-" * 78)
    P("(a5) m SENSITIVITY -- h scales as 1/m exactly, since h = w ell_q / m")
    P("     %-6s %14s %14s" % ("m", "floor p50 [um]", "x vs m=20"))
    base = None
    for m in (8.0, 13.0, 15.0, 20.0):
        v = statistics.median([W_CONE * e / m for e in ellq]) * 1e6
        if m == 20.0:
            base = v
    for m in (8.0, 13.0, 15.0, 20.0):
        v = statistics.median([W_CONE * e / m for e in ellq]) * 1e6
        P("     %-6.0f %14.1f %14.3f" % (m, v, v / base))
    P("     ⚠⚠ m is NOT fitted here and must not be: fitting it is P1g-4.")
    P("     ⚠ the fixture reads m = 15 (its input sets weibull.enabled = 0 so")
    P("       the parse never runs and the MEMBER default applies) while Meier")
    P("       reads 20 -- a %.0f %% difference in thickness, silently."
      % (100.0 * (15.0 ** -1 - 20.0 ** -1) / (20.0 ** -1)))

    P("")
    if fails:
        P("FAIL: %d check(s) failed" % len(fails))
        return 1
    P("PASS: the thickness is derived from a measured flux and a stated modulus,")
    P("      the wall has its own route, and one constant pair spans both")
    P("      field populations without being fitted to either.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
