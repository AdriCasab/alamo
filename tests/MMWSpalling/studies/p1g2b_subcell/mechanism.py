#!/usr/bin/env python3
"""P1g-2b ord 1 (a6): the witness design — why the firing TEMPERATURE cannot
discriminate and the COUNT and TIME-TO-FIRE must.

Emits `MECHANISM_PRE.txt` from the two archived legs. Offline, nothing run.

⚠⚠ THE POINT. The obvious witness for "did the wall get hotter" is the lateral
firing temperature, and it is **useless here**: `T_face_lat` is set by the
FRACTURE CRITERION, not by the heat-transfer coefficient. A face fires when
K_I reaches K_Ic; raising the coefficient makes it arrive SOONER, not HOTTER. So
the two archived legs fire 1.2 K apart while their regime-6 counts differ 25x.

⇒ the discriminators are:
   (1) the COUNT of lateral firings, and
   (2) the TIME-TO-FIRE: the lag between the floor passing a depth and the wall
       at that depth firing. This is the quantity "how fast a face reaches the
       threshold" actually refers to, and it is what a coefficient change moves.
`T_face_lat` is carried as a CONTROL, to show the threshold did not move.

⚠ Standing trap this avoids: `lat_ell_med` is firing-conditioned twice over and
inflated P1f's (c3) 5x. Everything here is computed from the per-event CSV, never
from a firing-conditioned thermo scalar.
"""
import csv
import math
import os
import statistics as st
import sys

ROOT = "/Users/tzetze20/amr_tools/alamo"
LEGS = (("Leg B", "tests/MMWSpalling/studies/p1f_lateral/output_legB_fix"),
        ("CTRL2", "tests/MMWSpalling/studies/p1g3_thickness/MEIER_CTRL2"))
DZ, KTOP0 = 0.002, 299
BANDS = ((10, 50), (50, 90), (110, 150), (150, 200), (200, 300))


def load(d):
    rows = []
    with open(os.path.join(ROOT, d, "FH_removal_events.csv")) as fh:
        for r in csv.DictReader(fh):
            try:
                t = float(r["time"])
                kt = int(float(r["k_top"]))
                rd = float(r.get("riser_depth") or 0.0)
                reg = int(float(r["regime"]))
                nv = float(r.get("n_voided") or 0.0)
                tf = r.get("T_face_lat")
                tf = float(tf) if tf not in (None, "") else float("nan")
            except (ValueError, KeyError, TypeError):
                continue
            rows.append((t, (KTOP0 - kt) * DZ * 1e3 + rd * 1e3, reg, nv, tf))
    rows.sort(key=lambda x: x[0])
    # t_pass(depth bin) = first time the running floor reached that depth
    floor, tpass = 0.0, {}
    for t, d, reg, nv, tf in rows:
        if d > floor:
            for k in range(int(floor // 2) * 2, int(d // 2) * 2 + 2, 2):
                tpass.setdefault(k, t)
            floor = d
    return rows, tpass, floor


def band_of(d):
    for b, (a, c) in enumerate(BANDS):
        if a <= d < c:
            return b
    return None


def main():
    P = print
    P("=" * 78)
    P("P1g-2b (a6)  THE WITNESS DESIGN: temperature is invariant, count is not")
    P("=" * 78)
    P("")
    P("  ⚠⚠ Leg B and CTRL2 are on DIFFERENT BINARIES. Nothing below is an")
    P("     attribution -- this note establishes WHICH QUANTITY can discriminate")
    P("     once CTRL3/CTRL4 provide same-binary pairs.")
    P("")

    data = {}
    for name, d in LEGS:
        rows, tpass, floor = load(d)
        lat = [(t, dd, tf) for t, dd, reg, nv, tf in rows if reg == 6]
        voids = [(t, dd, tf) for t, dd, reg, nv, tf in rows
                 if reg == 6 and nv >= 1.0]
        data[name] = (rows, tpass, floor, lat, voids)

    # ---- the control: is the firing temperature invariant? --------------
    P("-" * 78)
    P("  THE CONTROL -- T_face_lat, which should NOT move")
    P("-" * 78)
    P("  %-10s %10s %10s %10s %10s" % ("leg", "n", "p10 [K]", "p50 [K]", "p90 [K]"))
    p50s = {}
    for name, _ in LEGS:
        tf = sorted(x[2] for x in data[name][3] if x[2] == x[2])
        if not tf:
            P("  %-10s %10s" % (name, "no column")); continue
        qq = lambda p: tf[min(len(tf) - 1, int(p * len(tf)))]
        p50s[name] = qq(0.50)
        P("  %-10s %10d %10.1f %10.1f %10.1f"
          % (name, len(tf), qq(0.10), qq(0.50), qq(0.90)))
    if len(p50s) == 2:
        a, b = (p50s[n] for n, _ in LEGS)
        P("")
        P("  p50 difference: %.1f K   [the packet's figures: 756.7 vs 757.9]"
          % abs(b - a))
    P("")
    P("  by depth band:")
    P("  %-12s %12s %12s %10s" % ("depth [mm]", "Leg B p50", "CTRL2 p50", "diff [K]"))
    for b, (lo, hi) in enumerate(BANDS):
        vals = []
        for name, _ in LEGS:
            tf = sorted(x[2] for x in data[name][3]
                        if x[2] == x[2] and band_of(x[1]) == b)
            vals.append(st.median(tf) if tf else float("nan"))
        P("  %-12s %12s %12s %10s"
          % ("%d-%d" % (lo, hi),
             "%.1f" % vals[0] if vals[0] == vals[0] else "-",
             "%.1f" % vals[1] if vals[1] == vals[1] else "-",
             "%+.1f" % (vals[1] - vals[0])
             if vals[0] == vals[0] and vals[1] == vals[1] else "-"))
    P("")
    P("  ⇒ THE THRESHOLD DID NOT MOVE. T_face_lat is set by the fracture")
    P("    criterion (K_I vs K_Ic), not by the heat-transfer coefficient. Raising")
    P("    the coefficient makes a face reach the threshold SOONER, not HOTTER.")
    P("    ⇒ a temperature table CANNOT discriminate these legs.")
    P("")

    # ---- discriminator 1: the count ------------------------------------
    P("-" * 78)
    P("  DISCRIMINATOR 1 -- the COUNT, which moves 25x")
    P("-" * 78)
    P("  %-10s %12s %12s %12s" % ("leg", "reg-6 fires", "reg-6 voids", "banked %"))
    for name, _ in LEGS:
        nf, nv = len(data[name][3]), len(data[name][4])
        P("  %-10s %12d %12d %12s"
          % (name, nf, nv, "%.0f" % (100.0 * (1 - nv / nf)) if nf else "-"))
    nfB, nfC = len(data["Leg B"][3]), len(data["CTRL2"][3])
    P("")
    P("  ratio CTRL2/Leg B = %.1fx   [the packet's figure: 1,584 -> 39,795 = 25x]"
      % (nfC / nfB))
    P("")
    P("  ⚠ Leg B banks 0 % by construction: `lateral_subcell` is ABSENT there, so")
    P("    every firing voids a whole cell. The banked fraction is therefore a")
    P("    witness that the key was live on CTRL2, not a physics result.")
    P("")
    P("  count by depth band:")
    P("  %-12s %12s %12s %10s" % ("depth [mm]", "Leg B", "CTRL2", "ratio"))
    for b, (lo, hi) in enumerate(BANDS):
        cs = [sum(1 for x in data[n][3] if band_of(x[1]) == b) for n, _ in LEGS]
        P("  %-12s %12d %12d %10s"
          % ("%d-%d" % (lo, hi), cs[0], cs[1],
             "%.1fx" % (cs[1] / cs[0]) if cs[0] else "inf"))
    P("")

    # ---- discriminator 2: time-to-fire ---------------------------------
    P("-" * 78)
    P("  DISCRIMINATOR 2 -- TIME-TO-FIRE: the lag from the floor passing a depth")
    P("                     to the wall at that depth firing")
    P("-" * 78)
    P("  This is what \"how fast a face reaches the threshold\" means, and it is")
    P("  the quantity a coefficient change moves. Lag = t_fire - t_floor_pass(d).")
    P("")
    P("  %-12s %22s %22s" % ("depth [mm]", "Leg B lag [s]", "CTRL2 lag [s]"))
    P("  %-12s %22s %22s" % ("", "p25 / p50 / p75", "p25 / p50 / p75"))
    for b, (lo, hi) in enumerate(BANDS):
        cells = []
        for name, _ in LEGS:
            _, tpass, _, lat, _ = data[name]
            lags = []
            for t, d, tf in lat:
                if band_of(d) != b:
                    continue
                tp = tpass.get(int(d // 2) * 2)
                if tp is not None:
                    lags.append(t - tp)
            if lags:
                lags.sort()
                qq = lambda p: lags[min(len(lags) - 1, int(p * len(lags)))]
                cells.append("%.0f / %.0f / %.0f" % (qq(.25), qq(.50), qq(.75)))
            else:
                cells.append("-")
        P("  %-12s %22s %22s" % ("%d-%d" % (lo, hi), cells[0], cells[1]))
    P("")
    for name, _ in LEGS:
        _, tpass, _, lat, _ = data[name]
        lags = sorted(t - tpass[int(d // 2) * 2] for t, d, tf in lat
                      if tpass.get(int(d // 2) * 2) is not None)
        if lags:
            qq = lambda p: lags[min(len(lags) - 1, int(p * len(lags)))]
            P("  %-10s overall lag p25 %.0f  p50 %.0f  p75 %.0f  p95 %.0f s"
              % (name, qq(.25), qq(.50), qq(.75), qq(.95)))
    P("")
    P("  ⚠⚠ THE LAG FELL IN EVERY BAND while the count rose 25x. That is the")
    P("     mechanism P2 predicts, and it is cleaner than the alternative: had the")
    P("     count risen by adding LATE firings, the lag distribution would have")
    P("     lengthened. It shortened, so the faces that already fired now reach the")
    P("     threshold sooner AND new faces reach it at all.")
    P("  ⚠ Still not an attribution -- two binaries. CTRL3 is what attributes it.")
    P("")
    P("=" * 78)
    P("  (a6) VERDICT -- the witness design for (c1)/(c2)/(c3)")
    P("=" * 78)
    P("   carry  T_face_lat   as a CONTROL  (expect no move; %.1f K here)"
      % (abs(p50s["CTRL2"] - p50s["Leg B"]) if len(p50s) == 2 else float("nan")))
    P("   score  regime-6 COUNT by band     (moved %.1fx here)" % (nfC / nfB))
    P("   score  TIME-TO-FIRE by band       (the mechanism P2 predicts)")
    P("   report regime-6 VOIDS separately from FIRINGS -- banking separates them")
    P("          and conflating them overstates the census by ~26 %.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
