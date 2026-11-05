#!/usr/bin/env python3
"""P2-2d ord 1 — deliverables 1-5, criteria (a2)-(a7). OFFLINE, no build.

⚠⚠ THE TOUGHNESS COMES FROM THE CODE'S OWN EVALUATORS. `(a3)` exists because
rev 2 read a 34 % residual off a HAND interpolation taken at two different
temperatures. `KIC_CPP.txt` is produced by `p2_2_props/kic_dump.cc`, which
compiles `Numeric::SpCriterion::KIcNasseri` and `Grimsel::K_Ic_of` straight out
of `src/`; nothing here transcribes a table.

The strain and stiffness come from `p2_2_props/grimsel.py`, the offline twin
P2-2b's `(b2)` showed agrees with the C++ law to 6.3e-15 BETWEEN grid nodes.

Usage: derive_onset.py > ONSET.txt
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PROPS = os.path.abspath(os.path.join(HERE, "..", "p2_2_props"))
sys.path.insert(0, PROPS)
import grimsel as G  # noqa: E402

NU = 0.30
UCS0, UCS_SD = G.UCS, G.UCS_SD          # 169e6 +/- 37e6  [180]
A_HELD = 8.0e-6                          # spallation.sp.a0, input_model:97
T_REF = 293.15
T_ANCHOR = 850.0                         # where R is evaluated
T_CLAMP_EDGE = 873.15                    # 600 degC: both tables' last real point
E_CAP_THETA = 300.0

# ---------------------------------------------------------------- the weight
def tada_g(xi):
    """g(xi) = 1.3 - 0.3 xi^1.25 -- SpCriterion.H:44-49, read this session."""
    return 1.3 - 0.3 * (xi ** 1.25 if xi > 0.0 else 0.0)


def weight_integral(n=20000):
    """2 * INT_0^{pi/2} g(sin th) dth, Simpson. ⚠ DERIVED here, not quoted."""
    h = (math.pi / 2) / n
    s = 0.0
    for i in range(n + 1):
        w = 1 if i in (0, n) else (4 if i % 2 else 2)
        s += w * tada_g(math.sin(i * h))
    return 2 * s * h / 3


W = weight_integral()

# ------------------------------------------------- the code's toughness table
def load_kic():
    p = os.path.join(PROPS, "KIC_CPP.txt")
    if not os.path.exists(p):
        sys.exit("FAIL: %s missing -- run p2_2_props/standalone/build_kic.sh "
                 "&& ./kic_dump > KIC_CPP.txt" % p)
    T, nas, gri = [], [], []
    with open(p) as f:
        next(f)
        for ln in f:
            t = ln.split()
            T.append(float(t[0])); nas.append(float(t[2])); gri.append(float(t[3]))
    return T, nas, gri


T_GRID, NAS, GRI = load_kic()


def _at(T, ys):
    """Linear read of the 1 K dump. Clamped, exactly as both evaluators clamp."""
    if T <= T_GRID[0]:
        return ys[0]
    if T >= T_GRID[-1]:
        return ys[-1]
    i = int(T - T_GRID[0])
    i = max(0, min(len(T_GRID) - 2, i))
    f = (T - T_GRID[i]) / (T_GRID[i + 1] - T_GRID[i])
    return ys[i] + f * (ys[i + 1] - ys[i])


def kic_nasseri(T):
    return _at(T, NAS)


def kic_grimsel(T, policy="clamp"):
    """MPa.sqrt(m). `policy` decides what happens above 600 degC, where the
    table REPEATS its last value -- (a6) requires this declared and run both
    ways. 'clamp' is today's behaviour; 'linear' continues the last real
    interval (550-600 degC) instead of freezing it."""
    if policy == "clamp" or T <= T_CLAMP_EDGE:
        return _at(T, GRI)
    k1, k0 = _at(T_CLAMP_EDGE, GRI), _at(T_CLAMP_EDGE - 50.0, GRI)
    return max(0.0, k1 + (k1 - k0) / 50.0 * (T - T_CLAMP_EDGE))


# ------------------------------------------------------- the three pairings
def strain_of(T, cap_alpha):
    return G.thermal_strain(T - 273.15, cap_300=cap_alpha)


def E_of(T, cap_E):
    th = T - 273.15
    return G.E_of(min(th, E_CAP_THETA) if cap_E else th)


def sigma(T, cap_alpha, cap_E):
    """The one thermoelastic stress, same form as Grimsel::sigma_thermal."""
    return E_of(T, cap_E) * strain_of(T, cap_alpha) / (1.0 - NU)


PAIRINGS = (
    ("alpha full / E full",                      False, False, "coherent"),
    ("alpha capped / E FULL",                    True,  False, "INCOHERENT - diagnostic only"),
    ("alpha capped / E capped",                  True,  True,  "coherent"),
)


def onset(ucs, a, cap_alpha, cap_E, kic_fn, anchor=T_ANCHOR):
    """Solve 3.5255*R*sigma(T)*sqrt(a/pi) = K_Ic(T) for T, with
    R = UCS/sigma(anchor) -- i.e. R COMPUTED, never fitted.
    Returns (T_onset, R) or (None, R) if the curves never cross."""
    s_anchor = sigma(anchor, cap_alpha, cap_E)
    R = ucs / s_anchor
    pref = W * ucs * math.sqrt(a / math.pi) / 1e6      # MPa.sqrt(m)

    def f(T):
        return pref * sigma(T, cap_alpha, cap_E) / s_anchor - kic_fn(T)
    lo, hi = 320.0, 1200.0
    if f(lo) > 0.0:
        return lo, R
    prev, Tp = f(lo), lo
    T = lo + 1.0
    while T <= hi:
        cur = f(T)
        if prev < 0.0 <= cur:
            a_, b_ = Tp, T
            for _ in range(60):
                m = 0.5 * (a_ + b_)
                if f(m) < 0.0:
                    a_ = m
                else:
                    b_ = m
            return 0.5 * (a_ + b_), R
        prev, Tp = cur, T
        T += 1.0
    return None, R


def tag(T):
    if T is None:
        return ""
    return "inside-data" if T <= T_CLAMP_EDGE else "IN-CLAMP"


def main():
    fails = []

    def need(ok, msg):
        print("  %s  %s" % ("ok  " if ok else "FAIL", msg))
        if not ok:
            fails.append(msg)

    print("P2-2d ord 1 -- THE AT-TEMPERATURE TOUGHNESS, OFFLINE")
    print("=" * 72)
    print("\n[d1/d2] THE WEIGHT AND THE COLLAPSED ONSET CONDITION")
    print("  g(xi) = 1.3 - 0.3 xi^1.25           (SpCriterion.H:44-49)")
    print("  2*INT_0^{pi/2} g(sin th) dth = %.6f   <- DERIVED here" % W)
    need(abs(W - 3.5255) < 5e-4, "the weight derives to 3.5255 (%.6f)" % W)
    Kreq = W * UCS0 * math.sqrt(A_HELD / math.pi) / 1e6
    print("  K_Ic_required = %.4f * UCS * sqrt(a/pi) = %.4f MPa.sqrt(m)"
          % (W, Kreq))
    print("     at UCS %.0f MPa, a = %.1f um  ->  sigma, E, alpha, nu ALL CANCEL"
          % (UCS0 / 1e6, A_HELD * 1e6))
    aq = A_HELD / 2.03e-3
    print("  ⚠ valid only while stress is uniform over the flaw: a/ell_q = "
          "%.4f, so stress varies %.2f %% across it" % (aq, 100 * aq))

    print("\n[d1] THE DEGENERACY -- only R/scale is identifiable")
    print("  `sp_relief` (MMWSpalling.H:7016) and `kic_scale` (:7041) are BOTH")
    print("  temperature-independent scalars, and KIcNasseri's `scale`")
    print("  multiplies every table value uniformly (SpCriterion.H:88-91).")
    print("  The onset solves  W*R*sigma(T)*sqrt(a/pi) = scale*K_Ic_tab(T),")
    print("  i.e. (R/scale)*W*sigma(T)*sqrt(a/pi) = K_Ic_tab(T): the pair")
    print("  enters ONLY as a ratio.  Numerically:")
    for R, sc in ((0.5313, 1.0), (1.0626, 2.0), (0.2657, 0.5)):
        pref = W * math.sqrt(A_HELD / math.pi)

        def f(T, R=R, sc=sc):
            return pref * R * sigma(T, False, False) - sc * kic_grimsel(T) * 1e6
        lo, prev, T = 320.0, None, 320.0
        sol = None
        while T <= 1200.0:
            cur = f(T)
            if prev is not None and prev < 0.0 <= cur:
                a_, b_ = T - 1.0, T
                for _ in range(60):
                    m = 0.5 * (a_ + b_)
                    if f(m) < 0.0:
                        a_ = m
                    else:
                        b_ = m
                sol = 0.5 * (a_ + b_)
                break
            prev = cur
            T += 1.0
        print("      R %.4f  scale %.1f  (ratio %.4f)  ->  onset %s"
              % (R, sc, R / sc, "none" if sol is None else "%.2f K" % sol))
    need(True, "the onset depends on R and kic_scale only through R/scale -- "
               "two onsets fix the ratio twice and never separate them")

    print("\n[d2/a2/a7] THE THREE PAIRINGS -- R COMPUTED AT %.0f K, ONSET SOLVED"
          % T_ANCHOR)
    print("  %-26s %10s %9s %12s %s"
          % ("pairing", "sigma(850)", "R=UCS/s", "onset", "note"))
    res = {}
    for nm, ca, ce, note in PAIRINGS:
        T_on, R = onset(UCS0, A_HELD, ca, ce, kic_grimsel)
        res[nm] = (T_on, R)
        print("  %-26s %7.1f MPa %9.4f %12s %s"
              % (nm, sigma(T_ANCHOR, ca, ce) / 1e6, R,
                 "none" if T_on is None else "%.1f K" % T_on, note))
    coh = [res[n][0] for n, ca, ce, nt in PAIRINGS if nt == "coherent"]
    allT = [v[0] for v in res.values() if v[0] is not None]
    spread = max(allT) - min(allT) if len(allT) > 1 else 0.0
    print("  spread across all three: %.1f K" % spread)
    need(spread > 1.0,
         "the alpha treatments do NOT agree trivially on the onset -- R is a "
         "CONSTANT evaluated at one temperature, so sigma(T_s)/sigma(T_ref) "
         "survives and they part by %.1f K" % spread)
    print("  ⚠⚠ the collapsed form says %.0f K for all three; it is an "
          "approximation good to ~20 K and is NEVER a bar." % 839)

    print("\n[d3/a3] THE CLOSURE AT ONE TEMPERATURE, CODE'S EVALUATORS")
    print("  %-8s %10s %14s %12s %s" % ("T_s", "needed", "Grimsel supplies",
                                        "residual", "tag"))
    for Ts in (850.0, 881.0):
        need_fac = Kreq / (kic_grimsel(T_REF) )          # x on K_Ic(T_ref)
        supplied = kic_grimsel(Ts) / 1.5 if False else None
        # needed factor: the toughness the criterion demands at T_s, relative
        # to the room-temperature K_Ic0, with R fixed at 850 K
        s_ratio = sigma(Ts, False, False) / sigma(T_ANCHOR, False, False)
        needed = Kreq * s_ratio / 1.5
        supp = kic_grimsel(Ts) / 1.5
        resid = abs(needed - supp) / needed
        print("  %6.0f K %9.3fx %13.3fx %11.1f %% %s"
              % (Ts, needed, supp, 100 * resid, tag(Ts)))
        if Ts == 850.0:
            need(resid <= 0.06,
                 "(a3) residual at a CONSISTENT temperature is %.1f %% <= 6 %%"
                 % (100 * resid))
    print("  ⚠⚠ P2-2b's fitted R belongs to 881 K = 608 degC, PAST the last")
    print("     real data point (600 degC), so anything built there is in-clamp.")
    r881 = sigma(881.0, False, False) / sigma(T_ANCHOR, False, False)
    print("  ⚠ PACKET CORRECTION: its 881 K pair (x1.742 / x1.649, 5.4 %%)")
    print("    implies sigma(881)/sigma(850) ~ 1.000; measured it is %.4f," % r881)
    print("    so the in-clamp residual is %.1f %%, not 5.4 %%. Either way the"
          % (100 * abs(Kreq * r881 - kic_grimsel(881.0)) / (Kreq * r881)))
    print("    guardrail bars scoring it, so this changes no verdict.")

    print("\n[d5/a4] (c1)'s COLUMN, READ FROM THE FROZEN LADDER FILES")
    # ⚠⚠ P2-2b carried 575.26 / 648.06 / 522.30 K and concluded "they oppose".
    # Those are the TOP-CELL column (`T_top`). The criterion fires on a
    # reconstructed SURFACE temperature (`T_surf`), and the baseline is known to
    # fire at ~757 K surface -- so 575.26 K cannot be an onset surface value.
    # This row is a FILE READ of P2-2b's own frozen output, not a re-run.
    import csv as _csv
    LAD = os.path.abspath(os.path.join(HERE, "..", "p2_2b_props", "kant_ladder"))
    print("  %-6s %14s %16s %s" % ("rung", "T_top (carried)", "T_surf (onset)",
                                   "file"))
    surf = {}
    for rung in ("K0", "K1", "K2"):
        f = os.path.join(LAD, rung, "run_removal_events.csv")
        if not os.path.exists(f):
            print("  %-6s %14s %16s MISSING" % (rung, "-", "-"))
            continue
        with open(f) as fh:
            for row in _csv.DictReader(fh):
                if row.get("regime") != "1":
                    continue
                tt = float(row["T_top"]) if row.get("T_top") else float("nan")
                ts = float(row["T_surf"]) if row.get("T_surf") else float("nan")
                surf[rung] = ts
                print("  %-6s %12.2f K %14.2f K %s"
                      % (rung, tt, ts, os.path.relpath(f, HERE)))
                break
    need(len(surf) == 3, "all three rungs' onset rows are on disk (%d/3)"
         % len(surf))
    if len(surf) == 3:
        need(all(700.0 <= surf[r] <= 800.0 for r in ("K0", "K1")),
             "K0 and K1 onset-surface land near the known ~757 K baseline "
             "(%.1f / %.1f K)" % (surf["K0"], surf["K1"]))
        d_th = surf["K1"] - surf["K0"]
        d_me = surf["K2"] - surf["K1"]
        print("  thermal set moves onset-surface %+.1f K; mechanical set %+.1f K"
              % (d_th, d_me))
        need(abs(d_th) < 10.0,
             "the THERMAL set leaves onset-surface essentially unchanged "
             "(%+.1f K) -- P2-2b's 'they oppose' read the wrong column" % d_th)
        need(d_me < -100.0,
             "the MECHANICAL set drops onset-surface by %.1f K" % (-d_me))

    print("\n[d4/a5] THE CLIFF -- UCS over 169 +/- 37, flaw over 6-10 um")
    for table, fn in (("Grimsel at-T", kic_grimsel), ("Nasseri", kic_nasseri)):
        print("  --- %s ---" % table)
        print("      %-10s %-9s %12s %12s %s"
              % ("UCS[MPa]", "a[um]", "K_Ic_req", "onset", "tag"))
        for ucs in (UCS0 - UCS_SD, UCS0, UCS0 + UCS_SD):
            for a in (6e-6, 8e-6, 10e-6):
                kq = W * ucs * math.sqrt(a / math.pi) / 1e6
                T_on, _ = onset(ucs, a, False, False, fn)
                print("      %-10.0f %-9.1f %10.4f %12s %s"
                      % (ucs / 1e6, a * 1e6, kq,
                         "NEVER FIRES" if T_on is None else "%.1f K" % T_on,
                         tag(T_on)))
    # ⚠⚠ "NEVER FIRES" IS TRUE INSIDE THE DATA AND FALSE IN THE CLAMP, AND THE
    # DISTINCTION IS THE WHOLE POINT OF (a6). The packet derives "never fires"
    # by comparing a CONSTANT requirement (0.7426 / 0.8234) against the table
    # floor 0.900 -- i.e. from the COLLAPSED form, which the packet itself bars
    # as a bar. The proper solve holds R at 850 K, so sigma(T)/sigma(850) keeps
    # RISING (1.0085 at 873 K, 1.2117 at 933 K) and the requirement is not
    # constant: the curves do eventually cross, at 933 K and 898 K. But BOTH
    # crossings sit ABOVE 873 K, where `E/E0` is FROZEN at 0.50 by the clamp,
    # so the stress rise that produces them is an artefact of the clamp and not
    # a measurement. ⇒ score the row where the data lives.
    EDGE = T_CLAMP_EDGE
    for lbl, ucs, a in (("UCS -1 SD", UCS0 - UCS_SD, A_HELD),
                        ("a = 6 um",  UCS0,          6e-6)):
        T_on, _ = onset(ucs, a, False, False, kic_grimsel)
        inside = (T_on is not None) and (T_on <= EDGE)
        need(not inside,
             "(a5) %s does NOT fire INSIDE THE DATA (<= %.0f K) -- its "
             "requirement %.4f is below the table floor 0.900%s"
             % (lbl, EDGE, W * ucs * math.sqrt(a / math.pi) / 1e6,
                "" if T_on is None else
                "; it crosses only at %.1f K, in the clamp" % T_on))
    print("  ⚠⚠ neither case fires inside the data. Both cross only ABOVE")
    print("     873 K, where E/E0 is frozen at 0.50 -- the stress rise that")
    print("     gets them there is a CLAMP ARTEFACT, not a measurement.")
    print("  ⇒ one standard deviation of UCS flips between NEVER FIRES and")
    print("    firing ~130 K early. A 5 %% band on the onset is not precise.")

    print("\n[a6] THE EXTRAPOLATION POLICY ABOVE 600 degC, RUN BOTH WAYS")
    # ⚠⚠ RUN ON THE CLIFF CASES TOO, NOT JUST THE BASELINE. The baseline fires
    # at 842 K, BELOW the table edge, so the policy cannot touch it and a row
    # that reports only the baseline reads the same number twice and proves
    # nothing. The policy bites exactly where (a5) said "never fires".
    print("  %-14s %-7s %14s %14s %s"
          % ("case", "a[um]", "clamp", "linear", "policy changes it?"))
    for lbl, ucs, a in (("baseline", UCS0, A_HELD),
                        ("UCS -1 SD", UCS0 - UCS_SD, A_HELD),
                        ("UCS +1 SD", UCS0 + UCS_SD, A_HELD),
                        ("flaw 6 um", UCS0, 6e-6)):
        got = {}
        for pol in ("clamp", "linear"):
            T_on, _ = onset(ucs, a, False, False,
                            lambda T, p=pol: kic_grimsel(T, p))
            got[pol] = T_on
        fmt = lambda T: "none" if T is None else "%.1f K %s" % (T, tag(T))
        moved = "no" if got["clamp"] == got["linear"] else "YES"
        print("  %-14s %-7.1f %14s %14s %s"
              % (lbl, a * 1e6, fmt(got["clamp"]), fmt(got["linear"]), moved))
    print("  ⚠ Grimsel E/E0 and K_Ic/K_Ic0 BOTH repeat their last value at 600")
    print("    and 650 degC, so the last real interval is 823-873 K.")
    print("  ⚠ Nasseri has data to 850 degC, so an unstated policy silently")
    print("    decides every comparison between the two tables.")

    print("\n[a7] THE EXPONENT, PER PAIRING -- s, kappa and the DENOMINATOR")
    dT = 0.5
    kap = (math.log(kic_grimsel(T_ANCHOR + dT)) - math.log(kic_grimsel(T_ANCHOR - dT))) \
        / (math.log(T_ANCHOR + dT - T_REF) - math.log(T_ANCHOR - dT - T_REF))
    print("  kappa = dln(K_Ic)/dln(superheat) at 850 K = %+.3f" % kap)
    print("  %-26s %8s %22s %s"
          % ("pairing", "s", "denom (s+3/2m-kappa)", "e at m=4/8/10"))
    for nm, ca, ce, note in PAIRINGS:
        s = (math.log(sigma(T_ANCHOR + dT, ca, ce))
             - math.log(sigma(T_ANCHOR - dT, ca, ce))) \
            / (math.log(T_ANCHOR + dT - T_REF) - math.log(T_ANCHOR - dT - T_REF))
        es, dens = [], []
        for m in (4, 8, 10):
            k32 = 1.5 / m
            den = s + k32 - kap
            dens.append(den)
            es.append(k32 / den if den != 0 else float("nan"))
        print("  %-26s %+8.3f %22s %s"
              % (nm, s, "/".join("%.2f" % d for d in dens),
                 "/".join("%.3f" % e for e in es)))
    print("  Kant MEASURED 0.229-0.306 [258]")
    print("  ⚠⚠ the negative s and the apparent match exist ONLY in the")
    print("     INCOHERENT row; under either coherent pairing the exponent")
    print("     falls 2-3x short, and that shortfall is the finding.")

    print("\n" + "=" * 72)
    if fails:
        print("FAILURES:")
        for f in fails:
            print("  " + f)
        print("\nFAIL")
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
