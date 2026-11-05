#!/usr/bin/env python3
"""P2-3 -- the reductions, rows a5 and c1-c8.

⚠⚠ REUSED, NOT RE-DERIVED (deliverable 11). `reduce_price.py` owns the price,
the frozen floor and the firing-temperature population; `bore_profile.py` owns
Ø(z) from the void mask with the geometry parsed from the argv. Both are
imported from `studies/p2_2c_meier` with only a PATH parameter added, so P2-3
cannot disagree with P2-2f about what a price or a bore is.

⚠⚠ c6 IS BY REMOVED VOLUME, NOT BY FIRING COUNT, and that is a different
number. P2-2f reported 22.59 % of FIRINGS above the 873.15 K table ceiling; the
same leg's volume-weighted fraction is what this packet asks for, because a
firing that voids three cells carries three times the rock. Each firing's
`T_surf` is weighted by its `n_voided`.

⚠ Every new thermo column is PER-STEP STATE: averaged over the window, never
summed (standing rule 7(i)).

⚠ No artefact prints a verdict it cannot compute: a row with no value prints
"not computed" (P2-2f deliverable 6).
"""
import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
P22C = os.path.join(ROOT, "tests", "MMWSpalling", "studies", "p2_2c_meier")
sys.path.insert(0, P22C)
sys.path.insert(0, HERE)
import reduce_price as RP          # noqa: E402
import bore_profile as BP          # noqa: E402
import wall_profile as WP          # noqa: E402

LEGS_DIR = os.path.join(HERE, "legs")
T_CLAMP = 873.15
WIN = (150.0, 350.0)


def thermo(d):
    p = os.path.join(d, "run", "thermo.dat")
    if not os.path.exists(p):
        return None, []
    head = open(p).readline().split()
    rows = [l.split() for l in open(p).read().splitlines()[1:] if l.split()]
    return head, rows


def col(head, rows, name):
    if not head or name not in head:
        return []
    i = head.index(name)
    out = []
    for r in rows:
        if len(r) > i:
            try:
                out.append(float(r[i]))
            except ValueError:
                pass
    return out


def events(d):
    p = os.path.join(d, "run_removal_events.csv")
    return list(csv.DictReader(open(p))) if os.path.exists(p) else []


def win_mean(head, rows, name, lo, hi):
    """⚠ MEAN over the window. These columns are per-step STATE."""
    t = col(head, rows, "time")
    v = col(head, rows, name)
    sel = [x for ti, x in zip(t, v) if lo <= ti <= hi]
    return (sum(sel) / len(sel)) if sel else None


def clamp_by_volume(ev):
    """(vol_above, vol_tot, fraction) above the table ceiling, VOLUME-weighted."""
    va = vt = 0.0
    for r in ev:
        try:
            if int(float(r.get("regime", -1))) != RP.R_FLOOR:
                continue
            T = float(r.get("T_surf") or 0.0)
            n = float(r.get("n_voided") or 0.0)
        except (TypeError, ValueError):
            continue
        if T <= 0.0 or n <= 0.0:
            continue
        vt += n
        if T > T_CLAMP:
            va += n
    return va * RP.V_CELL, vt * RP.V_CELL, (va / vt if vt else None)


def verdict(tag, ok, note=""):
    if ok is None:
        print("     not computed  %s%s" % (tag, ("  -- " + note) if note else ""))
    else:
        print("     %s %s%s" % ("ok  " if ok else "MISS", tag,
                                ("  -- " + note) if note else ""))


def one_leg(tag, short=False):
    d = os.path.join(LEGS_DIR, tag + ("_s" if short else ""))
    if not os.path.isdir(d):
        print("\n--- %s: no leg directory (%s)" % (tag, d))
        return None
    head, rows = thermo(d)
    ev = events(d)
    print("\n" + "-" * 72)
    print("LEG %s   (%s)  %d firing events"
          % (tag, os.path.relpath(d, ROOT), len(ev)))
    if not head:
        print("  not computed -- no thermo.dat")
        return None
    lo, hi = (0.0, 60.0) if short else WIN
    res = {"tag": tag}

    # ---- a5: did the derived law REACH the bore wall? ----
    na = col(head, rows, "side_h_used_n_above")
    nall = col(head, rows, "side_h_used_n")
    dr = col(head, rows, "wall_datum_r")
    rb = col(head, rows, "wall_r_bore")
    print("  a5 the law's reach:")
    if not na:
        verdict("a5", None, "side_h_used_n_above column absent (derived law off?)")
    else:
        print("     side_h_used_n_above max %d (band-wide side_h_used_n max %d)"
              % (max(na), max(nall) if nall else 0))
        print("     wall_datum_r final %.6f m vs measured wall_r_bore %.6f m"
              % (dr[-1] if dr else 0.0, rb[-1] if rb else 0.0))
        print("     ⚠ the BAND-WIDE count is NOT the witness: it reads 8968 on")
        print("       the P2-2f m13 control, whose above-plane power is exactly")
        print("       0. Only the above-plane count discriminates.")
        verdict("a5  (bar > 0)", max(na) > 0)
        res["n_above"] = max(na)

    # ---- the ledger (b2 / c5) ----
    les = col(head, rows, "ledger_err_solid")
    worst = max((abs(x) for x in les[1:]), default=0.0)
    print("  ledger: max |ledger_err_solid| = %.3e  (bar <= 5e-3)" % worst)
    verdict("b2/c5", worst <= 5e-3)
    res["ledger"] = worst

    # ---- c1: the wall, from the march's own profile ----
    print("  c1 the gas up the bore:")
    wpr = WP.report(d, lo, hi, label=tag)
    pw = win_mean(head, rows, "jet_P_wall", lo, hi)
    Lf = win_mean(head, rows, "wall_L_fit", lo, hi)
    Tm = win_mean(head, rows, "wall_T_mouth", lo, hi)
    if pw is None:
        verdict("c1", None, "jet_P_wall absent -- the march is OFF on this leg")
    else:
        print("     window means: jet_P_wall %.1f W, wall_L_fit %.4f m, "
              "wall_T_mouth %.1f K" % (pw, Lf or 0.0, Tm or 0.0))
        print("     (bars: L 0.25-0.55 m, mouth 850-1150 K, P_wall 0.8-2.5 kW)")
        verdict("c1 L_fit", (0.25 <= (Lf or -1) <= 0.55))
        verdict("c1 mouth T", (850.0 <= (Tm or -1) <= 1150.0))
        verdict("c1 P_wall", (800.0 <= pw <= 2500.0))
        res.update({"P_wall": pw, "L_fit": Lf, "T_mouth": Tm})

    # ---- c2: the bore, with the isotropy witness ----
    print("  c2 bore Ø(z) -- from the VOID MASK, per azimuth sector:")
    BP.main(tag=tag, leg_root=d)

    # ---- c3: rate ----
    tt, nz = col(head, rows, "time"), col(head, rows, "nozzle_z")
    pts = [(t, z) for t, z in zip(tt, nz) if lo <= t <= hi and z > 0]
    if len(pts) > 1 and pts[-1][0] > pts[0][0]:
        rate = (pts[0][1] - pts[-1][1]) / (pts[-1][0] - pts[0][0]) * 3600.0
        print("  c3 rate over %.0f-%.0f s = %.3f m/h" % (pts[0][0], pts[-1][0], rate))
        print("     (prediction at x0.743: 2.3-2.9; Meier 1.2-1.8 quoted beside)")
        print("     ⚠ PROVISIONAL: 1 mm is parked (rule 12)")
        verdict("c3", 2.3 <= rate <= 2.9,
                "inside Meier's 1.2-1.8 would be a FINDING, not a pass")
        res["rate"] = rate
    else:
        verdict("c3", None, "no descent window")

    # ---- c4 / c7: volume and price ----
    cv = RP.charged_volume(ev)
    Ef = col(head, rows, "ledger_E_flake")
    print("  c4/c7 charged volume from n_voided (quarter, unscaled):")
    for nm in ("floor", "wall", "pad"):
        print("       %-6s %7d events %7d cells %9.5f L%s"
              % (nm, cv[nm]["events"], cv[nm]["cells"],
                 1e3 * cv[nm]["vol_nvoided"],
                 "   <- EXCLUDED (mechanical)" if nm == "pad" else ""))
    print("     charged = floor + wall = %.5f L; x4 = %.4f L (bar 1.7-2.1, "
          "Meier 1.270)" % (1e3 * cv["charged"], 1e3 * 4.0 * cv["charged"]))
    verdict("c4", 1.7 <= 1e3 * 4.0 * cv["charged"] <= 2.1)
    res["vol4"] = 1e3 * 4.0 * cv["charged"]
    if Ef and cv["charged"] > 0 and worst <= 5e-3:
        price = Ef[-1] / cv["charged"] / 1e9
        print("     price = E_flake %.4e J / %.5f L = %.3f GJ/m^3 (bar 1.20-1.30)"
              % (Ef[-1], 1e3 * cv["charged"], price))
        verdict("c7", 1.20 <= price <= 1.30)
        res["price"] = price
    else:
        verdict("c7", None, "ledger gates it" if worst > 5e-3 else "E_flake absent")

    # ---- c6: firing T and the clamp fraction BY VOLUME ----
    ts, note = RP.firing_T(ev)
    if ts:
        med = RP.pctl(ts, 0.5)
        va, vt, fr = clamp_by_volume(ev)
        print("  c6 median firing T_surf = %.2f K  (bar 845 +/- 10)" % med)
        print("     %s" % note)
        print("     clamp fraction BY REMOVED VOLUME: %.5f of %.5f L = %.2f %% "
              "(bar 12-25 %%)" % (1e3 * va, 1e3 * vt,
                                  100.0 * fr if fr else 0.0))
        print("     ⚠ by FIRING COUNT the same leg reads %.2f %% -- the two are"
              % (100.0 * RP.clamp_fraction(ev)[2]))
        print("       different numbers and the packet asks for the VOLUME one.")
        verdict("c6 firing T", abs(med - 845.0) <= 10.0)
        verdict("c6 clamp", (0.12 <= (fr or -1) <= 0.25))
        res.update({"T_fire": med, "clamp_vol": fr})
    else:
        verdict("c6", None, "no firing carries a reconstructed T_surf")
    return res


def main(want=None):
    want = want or ["W1", "W074", "NC"]
    print("P2-3 REDUCTIONS")
    print("=" * 72)
    print("⚠ window %0.0f-%0.0f s for the 500 s legs; per-step columns are "
          "MEANS over it, never sums." % WIN)
    got = {}
    for tag in want:
        for short in (False, True):
            r = one_leg(tag, short=short)
            if r:
                got[tag + ("_s" if short else "")] = r
                break
    # ---- b3: the negative control must balloon ----
    print("\n" + "=" * 72)
    print("b3 THE NEGATIVE CONTROL -- NC must widen the mouth, W1 must not")
    print("  (bar: Ø at 10 mm, NC - W1 >= 10 mm; W1 <= Meier 96 + 10 = 106)")
    print("  read the mouth rows from each leg's c2 table above: a single")
    print("  number here would hide which azimuth sector it came from.")
    print("\n" + "=" * 72)
    for k, v in got.items():
        print("  %-8s %s" % (k, ", ".join("%s=%.4g" % (a, b)
                                          for a, b in v.items()
                                          if isinstance(b, float))))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:] or None))
