#!/usr/bin/env python3
"""P1g-3 (c1): does the derived flake thickness narrow the hole?

Scores MEIER_LAT (both sites derived) against MEIER_CTRL2 (the shipped law), both
run on the SAME binary. ⚠ MEIER_BOTH is also reported where available: it ran
before the LATERAL site was wired, so its wall thickness stayed pinned at the
1.500 mm cap and it is a VERTICAL-ONLY leg. Keeping it in the table is what
separates the two halves of the law.

⚠⚠ THREE RULES THIS SCRIPT ENFORCES, because each has burned a packet here:

 1. **Never score the bore at 100 mm depth** (standing rule, Meier p. 223 pause
    underream). The 90-110 mm band is excluded by construction below.
 2. **A band near either leg's advancing FRONT is not comparable.** The two legs
    reach different depths, so a band that one has drilled and the other has not
    reports front position, not width. Every band is tagged with the margin to
    the shallower leg's deepest removal, and bands without clearance are marked
    FRONT and excluded from the verdict.
 3. **Quarter domain.** Radii are measured in the quarter, so diameter = 2*r_max
    and any volume needs x4.
"""
import collections, csv, math, os, statistics, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DX, KTOP0 = 0.002, 299
# ⚠ 90-110 mm is deliberately absent: never score the bore at 100 mm depth.
BANDS = ((10, 50), (50, 90), (110, 150), (150, 200), (200, 300), (300, 400))
FRONT_MARGIN_MM = 20.0      # a band must end this far above the shallower front
MEIER_BAND = (84.0, 96.0)   # Meier's measured bore diameter


def load(leg):
    f = os.path.join(HERE, leg, "FH_removal_events.csv")
    if not os.path.exists(f):
        return None
    rows = []
    for r in csv.DictReader(open(f)):
        try:
            int(float(r["regime"]))
            rows.append(r)
        except (ValueError, KeyError, TypeError):
            continue
    if not rows:
        return None
    rmax = {b: 0.0 for b in range(len(BANDS))}
    deepest = 0.0
    tmax = 0.0
    for r in rows:
        try:
            kt = int(float(r["k_top"]))
            rd = float(r.get("riser_depth") or 0.0)
            d = (KTOP0 - kt) * DX * 1e3 + rd * 1e3
            rr = DX * math.hypot(int(r["col_i"]) + 0.5,
                                 int(r["col_j"]) + 0.5) * 1e3
            t = float(r["time"])
        except (ValueError, KeyError, TypeError):
            continue
        deepest = max(deepest, d)
        tmax = max(tmax, t)
        for b, (a, c) in enumerate(BANDS):
            if a <= d < c and rr > rmax[b]:
                rmax[b] = rr
    r1 = [r for r in rows if int(float(r["regime"])) == 1 and r.get("h_f")]
    r6 = [r for r in rows if int(float(r["regime"])) == 6
          and r.get("h_f_lat") and float(r["h_f_lat"]) > 0]
    hl = [float(r["h_f_lat"]) for r in r6]
    return dict(
        n1=len(r1), n6=len(r6), deepest=deepest, tmax=tmax, rmax=rmax,
        hf=statistics.median([float(r["h_f"]) for r in r1]) * 1e3 if r1 else float("nan"),
        hl=statistics.median(hl) * 1e3 if hl else float("nan"),
        cap=100.0 * sum(1 for x in hl if abs(x - 0.0015) < 1e-12) / len(hl)
            if hl else float("nan"))


def main():
    P = print
    legs = {}
    for leg in ("MEIER_CTRL", "MEIER_CTRL2", "MEIER_BOTH", "MEIER_LAT"):
        d = load(leg)
        if d:
            legs[leg] = d
    if "MEIER_CTRL2" not in legs or "MEIER_LAT" not in legs:
        P("FAIL: need MEIER_CTRL2 and MEIER_LAT")
        return 1

    P("=" * 78)
    P("P1g-3 (c1)  DOES THE DERIVED THICKNESS NARROW THE HOLE?")
    P("=" * 78)
    P("  scored pair : MEIER_LAT (both sites derived) vs MEIER_CTRL2 (shipped law)")
    P("  ⚠ same binary for the scored pair -- the toolchain is not")
    P("    byte-reproducible, so a cross-build comparison would be unbacked.")
    P("")
    P("  %-26s %12s %12s %12s %12s"
      % ("", "CTRL", "CTRL2", "BOTH(vert)", "LAT(both)"))
    order = ("MEIER_CTRL", "MEIER_CTRL2", "MEIER_BOTH", "MEIER_LAT")
    for lbl, k, f in (("regime-1 firings", "n1", "%d"),
                      ("regime-6 firings", "n6", "%d"),
                      ("vertical h_f p50 [mm]", "hf", "%.4f"),
                      ("LATERAL h_f_lat p50 [mm]", "hl", "%.4f"),
                      ("lateral AT 1.5mm cap [%]", "cap", "%.1f"),
                      ("deepest removal [mm]", "deepest", "%.0f"),
                      ("last event [s]", "tmax", "%.1f")):
        P("  %-26s %12s %12s %12s %12s"
          % ((lbl,) + tuple((f % legs[n][k]) if n in legs else "-"
                            for n in order)))

    # ---- the bore, with the front rule applied -------------------------
    front = min(legs["MEIER_CTRL2"]["deepest"], legs["MEIER_LAT"]["deepest"])
    P("")
    P("  shallower leg's deepest removal = %.0f mm; a band must end %.0f mm"
      % (front, FRONT_MARGIN_MM))
    P("  above that to be comparable.")
    P("")
    P("  %-12s %10s %10s %9s %9s  %s"
      % ("depth [mm]", "CTRL2 O", "LAT O", "diff", "vs Meier", "status"))
    scored = []
    for b, (a, c) in enumerate(BANDS):
        dc = 2 * legs["MEIER_CTRL2"]["rmax"][b]
        db = 2 * legs["MEIER_LAT"]["rmax"][b]
        if dc == 0.0 and db == 0.0:
            continue
        comparable = c <= front - FRONT_MARGIN_MM
        # excess over the TOP of Meier's measured band
        exc = db - MEIER_BAND[1]
        tag = "scored" if comparable else "FRONT -- not comparable"
        if a <= 100 < c:
            tag = "EXCLUDED (100 mm rule)"
            comparable = False
        P("  %-12s %10.1f %10.1f %+9.1f %+9.1f  %s"
          % ("%d-%d" % (a, c), dc, db, db - dc, exc, tag))
        if comparable:
            scored.append((a, c, dc, db))

    P("")
    if not scored:
        P("  ⚠⚠ NO BAND IS COMPARABLE -- the verdict is UNSCOREABLE, not null.")
        return 1
    diffs = [db - dc for _, _, dc, db in scored]
    P("  SCORED BANDS: %d   mean diff %+.2f mm   worst %+.2f mm"
      % (len(scored), statistics.mean(diffs), min(diffs)))
    exc_c = statistics.mean([dc - MEIER_BAND[1] for _, _, dc, _ in scored])
    exc_l = statistics.mean([db - MEIER_BAND[1] for _, _, _, db in scored])
    P("  excess over Meier's %.0f mm: CTRL2 %+.1f mm -> LAT %+.1f mm"
      % (MEIER_BAND[1], exc_c, exc_l))
    if exc_c > 0:
        P("  => the derived thickness closes %.0f %% of the control's excess."
          % (100.0 * (exc_c - exc_l) / exc_c))
    P("")
    P("  ⚠⚠ THE COMPENSATION, which is why the effect is small:")
    for leg, ref in (("MEIER_BOTH", "vertical"), ("MEIER_LAT", "lateral")):
        if leg not in legs:
            continue
        c2 = legs["MEIER_CTRL2"]
        if ref == "vertical":
            tr = legs[leg]["hf"] / c2["hf"]
            fr = legs[leg]["n1"] / c2["n1"]
        else:
            tr = legs[leg]["hl"] / c2["hl"]
            fr = legs[leg]["n6"] / c2["n6"]
        P("     %-8s thickness x%.3f, firings x%.2f, product %.2f"
          % (ref, tr, fr, tr * fr))
    P("     A product near 1 means the volume is set by the ENERGY BUDGET and")
    P("     the thickness only decides how many events it is divided into.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
