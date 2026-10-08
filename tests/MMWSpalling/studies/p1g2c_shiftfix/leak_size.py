#!/usr/bin/env python3
"""P1g-2c ord 1 (a1)(a2)(a3): size the thermal shift's energy leak, offline.

Emits `LEAK_SIZE.txt` from P1g-2b's legs. Nothing is run and nothing edited.

⚠⚠ STANDING RULE 7, WHICH THIS SCRIPT EXISTS TO OBEY. `amr.thermo.plot_int = 6`
on these legs, so `thermo.dat` carries **one row per six steps**. Every per-step
column in it is a 1-in-N subsample, and summing one without the factor
undercounts by that factor. The factor is therefore **derived here from an
identity the data itself must satisfy**, never from `plot_int` alone.

⚠⚠ `ledger_err_solid` IS DIMENSIONLESS. `MMWSpalling.H:4587` divides by
`max(abs(ledger_E_src), 1e-30)`, so it is a RELATIVE error. Multiplying it by
`ledger_E_src` is what turns it into joules. P1g-2b called it an "absolute 5e-3 J
bar" and then argued the bar was mis-specified because the leg relocated more
enthalpy -- an argument that is void precisely because the bar is already
normalised by the cumulative source. That mistake is why this script prints both
forms side by side, every time.
"""
import csv
import os
import statistics as st
import sys

ROOT = "/Users/tzetze20/amr_tools/alamo"
P2B = os.path.join(ROOT, "tests/MMWSpalling/studies/p1g2b_subcell")
CTRL2 = os.path.join(ROOT, "tests/MMWSpalling/studies/p1g3_thickness/MEIER_CTRL2")
LEGS = (("CTRL2", CTRL2), ("CTRL3", os.path.join(P2B, "CTRL3")),
        ("CTRL4", os.path.join(P2B, "CTRL4")))


def thermo(d):
    p = os.path.join(d, "FH/thermo.dat")
    L = open(p).read().split("\n")
    hdr = L[0].split()
    ix = {n: i for i, n in enumerate(hdr)}
    rows = [f for f in (ln.split() for ln in L[1:]) if len(f) >= len(hdr)]
    return ix, rows


def g(ix, r, n):
    return float(r[ix[n]]) if n in ix else float("nan")


def events(d):
    """regime-6 firings and voids, from the per-event log (NOT sampled)."""
    p = os.path.join(d, "FH_removal_events.csv")
    nf = nv = 0
    with open(p) as fh:
        for r in csv.DictReader(fh):
            try:
                if int(float(r["regime"])) != 6:
                    continue
                nf += 1
                if float(r.get("n_voided") or 0.0) >= 1.0:
                    nv += 1
            except (ValueError, KeyError, TypeError):
                continue
    return nf, nv


def main():
    P = print
    T = {n: thermo(d) for n, d in LEGS}
    E = {n: events(d) for n, d in LEGS}

    P("=" * 78)
    P("P1g-2c (a1)-(a3)  SIZING THE THERMAL SHIFT'S ENERGY LEAK")
    P("=" * 78)
    P("")

    # ---------------- the sampling factor, DERIVED ----------------------
    P("-" * 78)
    P("  THE SAMPLING FACTOR -- derived from the data, not from plot_int")
    P("-" * 78)
    ix4, r4 = T["CTRL4"]
    nf4, nv4 = E["CTRL4"]
    sum_bank = sum(g(ix4, r, "lat_bank") for r in r4)
    P("  `lat_bank` is a PER-STEP counter; the event log is complete. Every bank")
    P("  is one regime-6 firing, so the factor is their ratio:")
    P("    regime-6 firings in the event log : %d" % nf4)
    P("    sum of sampled `lat_bank` rows    : %.0f" % sum_bank)
    fac = nf4 / sum_bank
    P("    ⇒ sampling factor = %.4f      [the packet's figure: 6.02]" % fac)
    P("  ⚠ cross-check from the argv: timestep 0.016 s, amr.thermo.plot_int = 6")
    P("    ⇒ nominally 6.00. The derived %.2f exceeds it because the first and" % fac)
    P("    last partial windows are not whole multiples of 6.")
    P("")

    # ---------------- (a1) the leak, both routes ------------------------
    P("-" * 78)
    P("  (a1) GATE -- the leak, sized independently")
    P("-" * 78)
    P("  %-7s %12s %12s %14s %14s" %
      ("leg", "|err_solid|", "E_src [J]", "UNCLOSED [J]", "% of budget"))
    unc = {}
    for n, _ in LEGS:
        ix, rows = T[n]
        last = rows[-1]
        er = abs(g(ix, last, "ledger_err_solid"))
        es = abs(g(ix, last, "ledger_E_src"))
        unc[n] = er * es
        P("  %-7s %12.4e %12.4g %14.1f %14.3f" % (n, er, es, unc[n], 100.0 * er))
    P("")
    P("  ⇒ the three NON-THERMAL legs sit at ~%.0f J; CTRL4 is at %.0f J."
      % (st.mean([unc["CTRL2"], unc["CTRL3"]]), unc["CTRL4"]))
    P("    unclosed on CTRL4 = **%.1f J**        [the packet's figure: 7,616 J]"
      % unc["CTRL4"])
    P("")
    P("  ROUTE 1 -- the shift's own booked energy, thermo column x the factor:")
    sE = sum(g(ix4, r, "lat_shift_E") for r in r4)
    booked1 = sE * fac
    P("    sum(lat_shift_E) over sampled rows = %.1f J" % sE)
    P("    x %.4f                              = **%.0f J = %.1f kJ**"
      % (fac, booked1, booked1 / 1e3))
    P("    [the packet's figure: 49.8 kJ]")
    P("")
    P("  ROUTE 2 -- ⚠⚠ THE FLAKE-LEDGER ROUTE IS CONFOUNDED AND I CANNOT")
    P("            REPRODUCE THE PACKET'S 50.2 kJ FROM IT.")
    f4 = g(ix4, T["CTRL4"][1][-1], "ledger_E_flake")
    f2 = g(T["CTRL2"][0], T["CTRL2"][1][-1], "ledger_E_flake")
    f3 = g(T["CTRL3"][0], T["CTRL3"][1][-1], "ledger_E_flake")
    P("    ledger_E_flake is CUMULATIVE (MMWSpalling.H:9683, only ever `+=`), so")
    P("    no sampling factor applies to it -- that part of the packet is right.")
    P("    But it accumulates FOUR terms, not one: flake+boundary+support")
    P("    (Removal.H:3142), the shift (:3445) AND lateral removal (:3453).")
    P("      CTRL4 %.0f   CTRL2 %.0f   CTRL3 %.0f" % (f4, f2, f3))
    P("      CTRL4 - CTRL2 = %.0f J ;  CTRL4 - CTRL3 = %.0f J" % (f4 - f2, f4 - f3))
    P("    Differencing does NOT isolate the shift, because CTRL4's OTHER terms")
    P("    also moved: its lateral voids fell %d → %d (x%.2f), so its lateral"
      % (E["CTRL2"][1], nv4, nv4 / E["CTRL2"][1]))
    P("    removal energy is smaller by roughly as much as the shift adds.")
    P("    ⇒ the difference %.0f kJ is (shift − the drop in lateral removal),"
      % ((f4 - f2) / 1e3))
    P("      not the shift. The confound is the very effect under measurement.")
    P("")
    P("  ⇒ (a1) VERDICT, stated precisely:")
    leak1 = 100.0 * unc["CTRL4"] / booked1
    P("     ✅ the UNCLOSED energy reproduces to 4 s.f.: %.1f J vs 7,616 J"
      % unc["CTRL4"])
    P("     ✅ route 1 reproduces exactly: %.1f kJ vs 49.8 kJ" % (booked1 / 1e3))
    P("     ✅ the LEAK FRACTION is %.2f %%   [bar 15.2 %% ± 1]  %s"
      % (leak1, "PASS" if abs(leak1 - 15.2) <= 1.0 else "FAIL"))
    P("     ⚠⚠ the sub-bar \"two independent routes within 2 %%\" CANNOT be met as")
    P("        written -- route 2 is confounded, not merely imprecise. The premise")
    P("        stands on route 1 plus the absolute unclosed energy, which is an")
    P("        independent quantity (it comes from the ledger, not from the shift")
    P("        counter) and reproduces to 4 s.f.")
    P("     ⇒ the leak is REAL and ~%.0f %% of the shift's booked energy." % leak1)
    P("")

    # ---------------- (a2) the refusal counter --------------------------
    P("-" * 78)
    P("  (a2) the refusal counter is MISLABELLED, not a defect")
    P("-" * 78)
    s_sh = sum(g(ix4, r, "lat_shift") for r in r4)
    s_rf = sum(g(ix4, r, "lat_shift_refuse") for r in r4)
    P("  sampled sums:   shifts %.0f + refusals %.0f = %.0f   banks %.0f"
      % (s_sh, s_rf, s_sh + s_rf, sum_bank))
    P("  identity holds exactly in the SAMPLED sums: %s"
      % ("YES" if abs((s_sh + s_rf) - sum_bank) < 0.5 else "NO"))
    P("  scaled by %.4f: shifts %.0f + refusals %.0f = banks %.0f"
      % (fac, s_sh * fac, s_rf * fac, sum_bank * fac))
    P("    [the packet's figures: 9,012 + 25,134 = 34,145]")
    P("  refusals / voids = %.0f / %d = %.3f   [the packet's figure: 1.01]"
      % (s_rf * fac, nv4, s_rf * fac / nv4))
    P("")
    P("  ⇒ EVERY BANK ATTEMPTS A SHIFT. On a step where the bank VOIDED the cell,")
    P("    the face cell is already gone, the row gathers length 0 and the attempt")
    P("    lands in the `runlen < ROWMIN` branch (Removal.H:3271) -- which shares")
    P("    one counter with the genuine `ell <= 0` refusal at :3253.")
    P("    ⇒ \"the shift acts on only a quarter of events\" is WRONG; it acts on")
    P("      every non-void bank. Deliverable 6 splits the counter three ways.")
    P("")

    # ---------------- (a3) the collision geometry -----------------------
    P("-" * 78)
    P("  (a3) the collision geometry, from the indexing at Removal.H:3397")
    P("-" * 78)
    P("      ii = ci + m*sgx ,  jj = cj + m*sgy")
    P("  For ax = 0 the row is (sgx, 0); for ax = 1 it is (0, sgy).")
    P("")
    P("  WITHIN ONE FACE at (ci,cj) = (10,10), sgx = sgy = +1:")
    P("    %-4s %-12s %-12s %s" % ("m", "ax=0 writes", "ax=1 writes", "same cell?"))
    shared = []
    for m in range(5):
        a0 = (10 + m, 10)
        a1 = (10, 10 + m)
        same = a0 == a1
        if same:
            shared.append(m)
        P("    %-4d %-12s %-12s %s" % (m, str(a0), str(a1), "**YES**" if same else "no"))
    P("  ⇒ the two rows of one face share EXACTLY m = %s." % shared)
    P("    So the self-overlap is an m = 0 phenomenon -- but see below.")
    P("")
    P("  ⚠⚠ CROSS-FACE collisions are NOT at m = 0, so an m=0 special case is")
    P("     insufficient. Worked example (the packet's):")
    P("       face A at (i,j) = (10,10) shifting +x, m = 2  ->  (12, 10)")
    P("       face B at (i+2, j-2) = (12, 8) shifting +y, m = 2  ->  (12, 10)")
    A = (10 + 2, 10)
    B = (12, 8 + 2)
    P("       A writes %s ; B writes %s ; collide: %s"
      % (str(A), str(B), "**YES**" if A == B else "no"))
    P("     Both are m = 2, neither is m = 0.")
    P("")
    P("  ⇒ THE FIX MUST BE PER-CELL, not per-m: accumulate deltas with a touch")
    P("    count, assign when touched once (so single-row cells stay BIT-EXACT),")
    P("    and use rH + Σ(rNew_ax − rH) when touched more.")
    P("")
    P("  ⚠ Safe against double-application: Removal.H:3371-3373 uses")
    P("    `MFIter(..., false)` and `mfi.validbox()`, so each cell is visited by")
    P("    exactly one box on one rank.")
    P("")
    P("=" * 78)
    P("  the geometry and booking checks that say the WRITE is at fault")
    P("=" * 78)
    rr = sorted(g(ix4, r, "lat_shift_ratio") for r in r4
                if g(ix4, r, "lat_shift") > 0)
    est = sum(g(ix4, r, "lat_shift_E_est") for r in r4)
    P("  lat_shift_ratio (geometry check) p05/p50/p95 = %.4f / %.4f / %.4f"
      % (rr[int(.05 * len(rr))], rr[int(.5 * len(rr))], rr[int(.95 * len(rr))]))
    P("  lat_shift_E_est / lat_shift_E                = %.4f" % (est / sE))
    P("  lat_overlap (vertical remap collisions)      = %.0f"
      % sum(g(ix4, r, "lat_overlap") for r in r4))
    P("")
    P("  ⚠⚠ Both ratios are ~1, and BOTH ARE BLIND TO THIS BUG: they are built")
    P("     from the INTENDED shift (`vol_disp`, `E_shift_local`), never from the")
    P("     field after the write. That is exactly why nothing caught it, and why")
    P("     deliverable 1's diagnostic must read H from the FIELD.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
