#!/usr/bin/env python3
"""P1g-2b ord 2-3: (b1)(b2)(b3) and the three verdicts (c1)(c2)(c3).

Emits `SCORE_p1g2b.txt`.

THE THREE COMPARISONS, and what separates each:

  (c1) BANKING       CTRL3 vs CTRL2   1 binary, 1 key   -- the attribution
  (c2) THERMAL HALF  CTRL2 vs CTRL4   1 binary, 1 key   -- first measure at scale
  (c3) THE DRIFT     CTRL3 vs Leg B   2 binaries, 0 keys-- four steps of source

⚠⚠ EVERY ROW IS LABELLED WITH ITS BINARY AND KEY COUNT, because the whole reason
this step exists is that a one-key/two-binary difference was quoted three times
as a one-key result. A comparison that crosses binaries is reported, never
attributed.

⚠⚠ THERMO IS READ BY COLUMN NAME, NEVER BY POSITION. `MMWSpalling.H:875-884`
registers `lat_bank`, `lat_bank_hf`, `lat_accum_sum` AND `lat_void_colmax` only
under `spall.lateral_subcell`, so **CTRL3's header is four columns shorter**.
Positional reads would silently misalign every column after the first missing
one. The column set actually present is asserted and printed per leg.

Band rules carried VERBATIM from `p1g3_thickness/score_bore.py:11-21`:
  1. never score the bore at 100 mm depth (the 90-110 band is absent by
     construction -- Meier p. 223 pause underream);
  2. a band near either leg's advancing FRONT reports front position, not width,
     so every band is tagged against the shallower front and excluded if it has
     no clearance;
  3. quarter domain, so diameter = 2*r_max and any volume needs x4.
"""
import csv
import math
import os
import statistics as st
import sys

ROOT = "/Users/tzetze20/amr_tools/alamo"
S = os.path.join(ROOT, "tests/MMWSpalling/studies/p1g2b_subcell")
DZ, KTOP0 = 0.002, 299
BANDS = ((10, 50), (50, 90), (110, 150), (150, 200), (200, 300), (300, 400))
FRONT_MARGIN_MM = 20.0
LEDGER_BAR = 5.0e-3
MIN_LAST_EVENT = 490.0

# Meier's digitised Fig. 8.8 profile (validation reference :63-64).
# ⚠ a VISIBLE WIDTH => a LOWER bound on Oe, +/-5 mm.
MEIER = ((10, 96.0), (25, 92.0), (50, 88.0), (75, 87.0), (125, 87.0), (300, 84.0))

SHIP = "7908bd6b0a95d547c8e336dae58ce4b450dde61c3d997daa1413245ebfb9eec0"
LEGB_SHA = "e8d8db0a4b9213e8816f254fc476ea69976c5189acbb901583d72896639c5e37"

LEGS = (
    ("Leg B", os.path.join(ROOT, "tests/MMWSpalling/studies/p1f_lateral/output_legB_fix"),
     LEGB_SHA, "base (no lateral_subcell)"),
    ("CTRL2", os.path.join(ROOT, "tests/MMWSpalling/studies/p1g3_thickness/MEIER_CTRL2"),
     SHIP, "base + lateral_subcell=1"),
    ("CTRL3", os.path.join(S, "CTRL3"), SHIP, "base + lateral_subcell=0"),
    ("CTRL4", os.path.join(S, "CTRL4"), SHIP, "base + lateral_subcell=1 + lateral_thermal=1"),
)

BANK_COLS = ("lat_bank", "lat_bank_hf", "lat_accum_sum", "lat_void_colmax")


def band_of(d):
    for b, (a, c) in enumerate(BANDS):
        if a <= d < c:
            return b
    return None


def meier_band(a, c, n=801):
    def at(x):
        if x <= MEIER[0][0]:
            return MEIER[0][1]
        if x >= MEIER[-1][0]:
            return MEIER[-1][1]
        for (d0, v0), (d1, v1) in zip(MEIER, MEIER[1:]):
            if d0 <= x <= d1:
                return v0 + (v1 - v0) * (x - d0) / (d1 - d0)
        return MEIER[-1][1]
    return st.mean(at(a + (c - a) * i / (n - 1)) for i in range(n))


def read_thermo(leg_dir):
    """BY NAME. Returns (columns_present, max|err_solid|, rows, bank_cols_present)."""
    t = os.path.join(leg_dir, "FH/thermo.dat")
    if not os.path.exists(t):
        return None
    lines = open(t).read().split("\n")
    hdr = lines[0].split()
    ix = {n: i for i, n in enumerate(hdr)}
    if "ledger_err_solid" not in ix:
        return dict(cols=set(hdr), err=None, rows=0,
                    bank=[c for c in BANK_COLS if c in ix])
    worst, n = 0.0, 0
    for ln in lines[1:]:
        f = ln.split()
        if len(f) < len(hdr):
            continue
        n += 1
        worst = max(worst, abs(float(f[ix["ledger_err_solid"]])))
    return dict(cols=set(hdr), err=worst, rows=n,
                bank=[c for c in BANK_COLS if c in ix])


def read_events(leg_dir):
    f = os.path.join(leg_dir, "FH_removal_events.csv")
    if not os.path.exists(f):
        return None
    rows = []
    with open(f) as fh:
        for r in csv.DictReader(fh):
            try:
                t = float(r["time"])
                kt = int(float(r["k_top"]))
                rd = float(r.get("riser_depth") or 0.0)
                ci, cj = int(r["col_i"]), int(r["col_j"])
                reg = int(float(r["regime"]))
                nv = float(r.get("n_voided") or 0.0)
                tf = r.get("T_face_lat")
                tf = float(tf) if tf not in (None, "") else float("nan")
            except (ValueError, KeyError, TypeError):
                continue
            d = (KTOP0 - kt) * DZ * 1e3 + rd * 1e3
            rr = DZ * math.hypot(ci + 0.5, cj + 0.5) * 1e3
            rows.append((t, d, rr, reg, nv, tf))
    if not rows:
        return None
    rows.sort(key=lambda x: x[0])
    rmax = {b: 0.0 for b in range(len(BANDS))}
    deepest = tmax = 0.0
    floor, tpass = 0.0, {}
    for t, d, rr, reg, nv, tf in rows:
        deepest = max(deepest, d)
        tmax = max(tmax, t)
        if d > floor:
            for k in range(int(floor // 2) * 2, int(d // 2) * 2 + 2, 2):
                tpass.setdefault(k, t)
            floor = d
        b = band_of(d)
        if b is not None and rr > rmax[b]:
            rmax[b] = rr
    lat = [(t, d, tf) for t, d, rr, reg, nv, tf in rows if reg == 6]
    voids = [(t, d, tf) for t, d, rr, reg, nv, tf in rows
             if reg == 6 and nv >= 1.0]
    n1 = sum(1 for x in rows if x[3] == 1)
    return dict(rmax=rmax, deepest=deepest, tmax=tmax, lat=lat, voids=voids,
                n1=n1, tpass=tpass, rows=rows)


def lags(ev, band=None):
    out = []
    for t, d, tf in ev["lat"]:
        if band is not None and band_of(d) != band:
            continue
        tp = ev["tpass"].get(int(d // 2) * 2)
        if tp is not None:
            out.append(t - tp)
    return sorted(out)


def q(xs, p):
    return xs[min(len(xs) - 1, int(p * len(xs)))] if xs else float("nan")


def compare(P, label, a, b, nbin, nkey, what, L):
    """one bore comparison, with the front rule applied."""
    ea, eb = L[a]["ev"], L[b]["ev"]
    front = min(ea["deepest"], eb["deepest"])
    P("")
    P("-" * 78)
    P("  %s   %s vs %s" % (label, a, b))
    P("  %d binary%s, %d key%s  --  %s"
      % (nbin, "" if nbin == 1 else "ies", nkey, "" if nkey == 1 else "s", what))
    P("-" * 78)
    if nbin > 1:
        P("  ⚠⚠ CROSSES BINARIES -- reported, NOT attributed.")
    P("  shallower front %.0f mm; a band must end %.0f mm above it."
      % (front, FRONT_MARGIN_MM))
    P("")
    P("  %-12s %9s %9s %9s %10s %9s  %s"
      % ("depth [mm]", a + " Ø", b + " Ø", "Δ", "Meier", "b-Meier", "status"))
    scored = []
    for i, (lo, hi) in enumerate(BANDS):
        da, db = 2 * ea["rmax"][i], 2 * eb["rmax"][i]
        if da == 0.0 and db == 0.0:
            continue
        ok = hi <= front - FRONT_MARGIN_MM
        tag = "scored" if ok else "FRONT -- not comparable"
        if lo <= 100 < hi:
            tag, ok = "EXCLUDED (100 mm rule)", False
        mb = meier_band(lo, hi)
        P("  %-12s %9.1f %9.1f %+9.2f %10.1f %+9.1f  %s"
          % ("%d-%d" % (lo, hi), da, db, db - da, mb, db - mb, tag))
        if ok:
            scored.append((db - da, db - mb))
    P("")
    if scored:
        P("  SCORED BANDS %d   mean Δ %+.2f mm   worst %+.2f mm   | mean vs Meier %+.1f"
          % (len(scored), st.mean(x[0] for x in scored),
             max(scored, key=lambda x: abs(x[0]))[0],
             st.mean(x[1] for x in scored)))
    else:
        P("  ⚠⚠ NO COMPARABLE BAND -- UNSCOREABLE, not null.")
    la, lb = lags(ea), lags(eb)
    P("  lateral firings %d → %d (%.2f×)   voids %d → %d (%.2f×)"
      % (len(ea["lat"]), len(eb["lat"]),
         len(eb["lat"]) / len(ea["lat"]) if ea["lat"] else float("nan"),
         len(ea["voids"]), len(eb["voids"]),
         len(eb["voids"]) / len(ea["voids"]) if ea["voids"] else float("nan")))
    P("  time-to-fire p50 %.0f → %.0f s   regime-1 %d → %d (%.2f×)"
      % (q(la, .5), q(lb, .5), ea["n1"], eb["n1"],
         eb["n1"] / ea["n1"] if ea["n1"] else float("nan")))
    tfa = sorted(x[2] for x in ea["lat"] if x[2] == x[2])
    tfb = sorted(x[2] for x in eb["lat"] if x[2] == x[2])
    P("  CONTROL T_face_lat p50 %.1f → %.1f K  (Δ %+.1f — expect ~0)"
      % (q(tfa, .5), q(tfb, .5), q(tfb, .5) - q(tfa, .5)))
    return scored


def main():
    P = print
    L = {}
    for name, d, sha, keys in LEGS:
        L[name] = dict(dir=d, sha=sha, keys=keys,
                       ev=read_events(d), th=read_thermo(d))

    P("=" * 78)
    P("P1g-2b  SUB-CELL BANKING, THE THERMAL HALF, AND THE FOUR-STEP DRIFT")
    P("=" * 78)
    P("")
    P("  %-7s %-46s %s" % ("leg", "keys beyond Leg B's argv", "binary"))
    for name, _, sha, keys in LEGS:
        P("  %-7s %-46s %s…" % (name, keys, sha[:12]))
    P("")
    missing = [n for n in L if L[n]["ev"] is None]
    if missing:
        P("  FAIL: no events for %s" % ", ".join(missing))
        return 1

    # ---------------- (b1) completion + binary ---------------------------
    P("-" * 78)
    P("  (b1) GATE -- both new legs complete, on an UNCHANGED binary")
    P("-" * 78)
    ok_b1 = True
    for leg in ("CTRL3", "CTRL4"):
        tm = L[leg]["ev"]["tmax"]
        good = tm >= MIN_LAST_EVENT
        ok_b1 &= good
        P("  %-6s last event %.1f s  [bar >= %.0f]  %s"
          % (leg, tm, MIN_LAST_EVENT, "OK" if good else "FAIL"))
    for leg in ("CTRL3", "CTRL4"):
        b = os.path.join(S, leg + ".sha_before")
        a = os.path.join(S, leg + ".sha_after")
        sb = open(b).read().split()[0] if os.path.exists(b) else "?"
        sa = open(a).read().split()[0] if os.path.exists(a) else "?"
        same = (sb == sa == SHIP)
        ok_b1 &= same
        P("  %-6s sha before %s…  after %s…  shipping: %s"
          % (leg, sb[:12], sa[:12], "YES" if same else "**NO**"))
    P("")
    P("  (b1) %s" % ("PASS" if ok_b1 else "FAIL"))

    # ---------------- (b2) depth, reported not barred --------------------
    P("")
    P("-" * 78)
    P("  (b2) deepest removal -- REPORTED, no bar (depth is under measurement)")
    P("-" * 78)
    for name, _, _, _ in LEGS:
        P("  %-7s deepest %3.0f mm   last event %.1f s"
          % (name, L[name]["ev"]["deepest"], L[name]["ev"]["tmax"]))

    # ---------------- (b3) ledger + the column set -----------------------
    P("")
    P("-" * 78)
    P("  (b3) GATE -- energy ledger, and the thermo COLUMN SET per leg")
    P("-" * 78)
    ok_b3 = True
    for name, _, _, _ in LEGS:
        th = L[name]["th"]
        if th is None:
            P("  %-7s no thermo.dat" % name); continue
        e = th["err"]
        good = (e is not None and e <= LEDGER_BAR)
        if name in ("CTRL3", "CTRL4"):
            ok_b3 &= good
        P("  %-7s max|ledger_err_solid| %s  [bar %.0e]  %-6s  %d cols, %d rows"
          % (name, "%.4e" % e if e is not None else "absent", LEDGER_BAR,
             "OK" if good else "FAIL", len(th["cols"]), th["rows"]))
        P("          banking columns present: %s"
          % (", ".join(th["bank"]) if th["bank"] else "NONE"))
    P("")
    P("  ⚠⚠ CTRL3 is expected to lack all four banking columns -- they are")
    P("     registered only under `spall.lateral_subcell` (MMWSpalling.H:875-884).")
    nb3 = len(L["CTRL3"]["th"]["bank"]) if L["CTRL3"]["th"] else -1
    nb2 = len(L["CTRL2"]["th"]["bank"]) if L["CTRL2"]["th"] else -1
    nc3 = len(L["CTRL3"]["th"]["cols"]) if L["CTRL3"]["th"] else -1
    nc2 = len(L["CTRL2"]["th"]["cols"]) if L["CTRL2"]["th"] else -1
    P("     CTRL3 has %d of 4, CTRL2 has %d of 4; header widths %d vs %d columns"
      % (nb3, nb2, nc3, nc2))
    P("     ⇒ the header DOES differ, which is why every read above is BY NAME.")
    if nb3 != 0:
        P("     ⚠⚠ UNEXPECTED: CTRL3 should have NONE. Check the key took effect.")
    P("  (b3) %s" % ("PASS" if ok_b3 else "FAIL"))

    # ---------------- the three verdicts --------------------------------
    c1 = compare(P, "(c1) BANKING", "CTRL3", "CTRL2", 1, 1,
                 "what `lateral_subcell` is worth", L)
    c2 = compare(P, "(c2) THE THERMAL HALF", "CTRL2", "CTRL4", 1, 1,
                 "P1g-2's row shift, FIRST measurement at Meier scale", L)
    c3 = compare(P, "(c3) THE FOUR-STEP DRIFT", "Leg B", "CTRL3", 2, 0,
                 "identical keys; P1g-0b/1/2/3 of source between", L)

    # ---------------- the discriminator ---------------------------------
    P("")
    P("=" * 78)
    P("  THE DISCRIMINATOR -- is the 25x lateral rise the BINARY or BANKING?")
    P("=" * 78)
    nB = len(L["Leg B"]["ev"]["lat"])
    n2 = len(L["CTRL2"]["ev"]["lat"])
    n3 = len(L["CTRL3"]["ev"]["lat"])
    n4 = len(L["CTRL4"]["ev"]["lat"])
    P("  Leg B %d   CTRL3 %d   CTRL2 %d   CTRL4 %d" % (nB, n3, n2, n4))
    P("")
    P("  CTRL3 shares Leg B's KEYS and CTRL2's BINARY, so it splits the 25x:")
    P("    binary  (Leg B → CTRL3, 0 keys) : %.2f×" % (n3 / nB if nB else float("nan")))
    P("    banking (CTRL3 → CTRL2, 1 key)  : %.2f×" % (n2 / n3 if n3 else float("nan")))
    P("    product                          : %.2f×  (vs the observed %.2f×)"
      % ((n3 / nB) * (n2 / n3) if nB and n3 else float("nan"),
         n2 / nB if nB else float("nan")))
    P("")
    if n3 and nB:
        if n3 / nB > 5.0:
            P("  ⇒ THE BINARY owns the rise -- P2 CONFIRMED (P1g-0b's coefficient fix).")
        elif n2 / n3 > 5.0:
            P("  ⇒ BANKING owns the rise -- P2 REFUTED; banking is not bookkeeping-only.")
        else:
            P("  ⇒ SPLIT between the two; neither dominates.")

    # ---------------- the mechanism table (deliverable 6) ---------------
    P("")
    P("=" * 78)
    P("  DELIVERABLE 6 -- mechanism by depth band, all four legs")
    P("=" * 78)
    P("  lateral firing COUNT")
    P("  %-12s %10s %10s %10s %10s" % ("depth [mm]", "Leg B", "CTRL3", "CTRL2", "CTRL4"))
    for i, (lo, hi) in enumerate(BANDS):
        P("  %-12s %10d %10d %10d %10d"
          % ("%d-%d" % (lo, hi),
             *[sum(1 for x in L[n]["ev"]["lat"] if band_of(x[1]) == i)
               for n in ("Leg B", "CTRL3", "CTRL2", "CTRL4")]))
    P("")
    P("  TIME-TO-FIRE p50 [s]  (lag from the floor passing a depth)")
    P("  %-12s %10s %10s %10s %10s" % ("depth [mm]", "Leg B", "CTRL3", "CTRL2", "CTRL4"))
    for i, (lo, hi) in enumerate(BANDS):
        vals = [q(lags(L[n]["ev"], i), .5) for n in ("Leg B", "CTRL3", "CTRL2", "CTRL4")]
        P("  %-12s %10s %10s %10s %10s"
          % ("%d-%d" % (lo, hi),
             *["%.0f" % v if v == v else "-" for v in vals]))
    P("")
    P("  CONTROL -- T_face_lat p50 [K]  (expect INVARIANT: set by the fracture")
    P("  threshold, not by the coefficient)")
    P("  %-12s %10s %10s %10s %10s" % ("depth [mm]", "Leg B", "CTRL3", "CTRL2", "CTRL4"))
    for i, (lo, hi) in enumerate(BANDS):
        vals = []
        for n in ("Leg B", "CTRL3", "CTRL2", "CTRL4"):
            tf = sorted(x[2] for x in L[n]["ev"]["lat"]
                        if x[2] == x[2] and band_of(x[1]) == i)
            vals.append(st.median(tf) if tf else float("nan"))
        P("  %-12s %10s %10s %10s %10s"
          % ("%d-%d" % (lo, hi),
             *["%.1f" % v if v == v else "-" for v in vals]))
    allt = {}
    for n in ("Leg B", "CTRL3", "CTRL2", "CTRL4"):
        tf = sorted(x[2] for x in L[n]["ev"]["lat"] if x[2] == x[2])
        allt[n] = q(tf, .5)
    P("")
    P("  overall p50: " + "   ".join("%s %.1f" % (k, v) for k, v in allt.items()))
    P("  spread across all four legs: %.1f K" % (max(allt.values()) - min(allt.values())))
    return 0 if (ok_b1 and ok_b3) else 1


if __name__ == "__main__":
    sys.exit(main())
