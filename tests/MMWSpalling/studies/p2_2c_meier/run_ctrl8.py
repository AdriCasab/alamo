#!/usr/bin/env python3
"""P2-2f (rev 2) — CTRL8 on Meier. The ladder is CTRL8 alone (user, 2026-10-12).

Run order: short 60 s (b1, b2) -> m = 13 (b3, b4, c1-c4) -> m = 8 (c5).
⚠ m = 13 FIRST: b3 halts the ladder on a frozen floor, and ordering m = 8
first would mean a freeze there teaches nothing about 13.

⚠⚠ CELL MODE THROUGHOUT -- no `size_volume`, no `flake_length`, no `V_ref`.
`weibull.V0` stays per deck. Setting V_ref here would give the Meier legs a
flux dependence nobody intended.

⚠⚠ THE PRICE CHARGES BOTH HEAT-DRIVEN PATHS -- floor (regime 1) AND wall
(regime 6) -- COUNTED FROM `n_voided`, with pad clearing (regime 4, mechanical)
EXCLUDED from both sides. Floor-only produced the withdrawn 1.36x and the
withdrawn 1.27 GJ/m^3, and the second outlived the retraction of the first.
⚠ The `lat_bank_hf x 6.0226` route named in rev 1 is RETIRED -- see
`reduce_price.py`, which owns every reduction here.

Usage: run_ctrl8.py short | m13 | m8
"""
import csv
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
DECK = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "meier",
                    "sp_meier_pilot", "input_drilling")
OUT = os.path.join(HERE, "ctrl8")
ARGV_FILE = os.path.join(HERE, "feet_argv.txt")
# ⚠⚠ LAT_FACTOR IS DELETED, NOT MOVED. ARCHIVE_DONE L18670 states 6.0226 is
# "a sampling factor derived from the event log" OF P1g-2c'S RUN. It is not a
# constant, it does not transfer, and it was the route by which a per-step
# column got read as a total. The wall path is counted from `n_voided`.
sys.path.insert(0, HERE)
import reduce_price as RP  # noqa: E402


def base_keys():
    return [l.strip() for l in open(ARGV_FILE) if l.strip()]


def run(tag, stop_time, extra):
    d = os.path.join(OUT, tag)
    os.makedirs(d, exist_ok=True)
    argv = (["caffeinate", "-dis", "mpirun", "--oversubscribe", "--bind-to",
             "none", "-np", "4", BIN, DECK]
            + base_keys() + extra
            + ["stop_time=%g" % stop_time,
               "amr.thermo.int=1", "amr.thermo.plot_int=6",
               "energy_ledger.enabled=1", "spall.removal_events_csv=1",
               "plot_file=%s" % os.path.join(d, "run")])
    with open(os.path.join(d, "argv.sh"), "w") as f:
        f.write(" ".join(argv) + "\n")
    env = os.environ.copy()
    env.setdefault("OMPI_MCA_btl", "self,vader,tcp")
    import time as _t
    t0 = _t.time()
    with open(os.path.join(d, "run.log"), "w") as lg:
        rc = subprocess.call(argv, cwd=ROOT, env=env, stdout=lg,
                             stderr=subprocess.STDOUT)
    return rc, _t.time() - t0, d


def events(d):
    p = os.path.join(d, "run_removal_events.csv")
    return list(csv.DictReader(open(p))) if os.path.exists(p) else []


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


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "short"
    print("P2-2f (rev 2) CTRL8 -- %s" % which)
    print("  binary %s" % subprocess.run(["shasum", "-a", "256", BIN],
                                         capture_output=True,
                                         text=True).stdout[:24])
    print("  keys from %s (%d)" % (os.path.relpath(ARGV_FILE, ROOT),
                                   len(base_keys())))

    if which == "short":
        rc, wall, d = run("short60", 60.0, [])
        print("  rc=%d  wall %.1f s (%.2f wall-s per simulated s)"
              % (rc, wall, wall / 60.0))
        ev = events(d)
        print("  firing events: %d" % len(ev))
        head, rows = thermo(d)
        # b1: the descent rate.
        # ⚠⚠ SKIP THE t = 0 ROW. `nozzle_z` is UNSET there and logs as 0.0, so
        # including it reads the drop as 0.0 -> 0.52 and gives v = -31.2 m/h --
        # a negative rate from a column that simply had not been written yet.
        # The real descent is 0.550 -> 0.520 over 59.8 s = 1.806 m/h.
        tt = col(head, rows, "time")
        for cand in ("nozzle_z", "feet_z", "surf_z_min"):
            raw = col(head, rows, cand)
            if not raw:
                continue
            pts = [(t, z) for t, z in zip(tt, raw) if t > 0.0 and z > 0.0]
            if len(pts) < 2:
                continue
            drop = pts[0][1] - pts[-1][1]
            dt = pts[-1][0] - pts[0][0]
            print("  %-12s %.6f -> %.6f over t=%.3f-%.3f  (drop %.4f m)"
                  % (cand, pts[0][1], pts[-1][1], pts[0][0], pts[-1][0], drop))
            vrate = drop / dt * 3600.0 if dt > 0 else 0.0
            print("  b1 descent rate v = %.3f m/h   %s (bar v > 0, rc 0)"
                  % (vrate, "ok  " if vrate > 0 and rc == 0 else "MISS"))
            if vrate <= 0:
                print("     no t_stop can be projected from a non-positive rate")
                break
            t_stop = min(500.0, 1620.0 / vrate)
            print("  => t_stop = min(500, 1620/v) = %.1f s" % t_stop)
            # b2: the margin, PROJECTED from v -- this is a pre-run gate, taken
            # before a 31-minute leg is committed.
            # datum: rock runs from the top (prob_hi[2] = 0.5) to the domain
            # bottom (prob_lo[2] = -0.1) = 600 mm, less a 40 mm bottom guard
            # => 560 mm usable.
            ROCK, GUARD = 0.600, 0.040
            usable = ROCK - GUARD
            proj = vrate / 3600.0 * t_stop
            margin = usable - proj
            print("  b2 projected depth at t_stop = %.4f m; usable %.3f m "
                  "(600 mm rock - 40 mm guard)" % (proj, usable))
            print("     margin %.0f mm   %s (bar >= 60 mm)"
                  % (1e3 * margin, "ok  " if margin >= 0.060 else "MISS"))
            print("     ⚠ THIS IS A PROJECTION AT A CONSTANT RATE AND THE RATE")
            print("       ACCELERATES. Rev 1 projected 251 mm and the leg")
            print("       actually reached 489 mm, leaving 71 mm -- still over")
            print("       the bar, but 4x tighter than the projection says. The")
            print("       m13 reduction re-checks it against the REALISED floor.")
            break
        else:
            print("  ⚠ no usable descent column in thermo.dat; columns are:")
            print("    %s" % ", ".join(head or []))
        return 0 if rc == 0 else 1

    # `reduce <tag>` re-scores an EXISTING leg's output with no run. It exists
    # so the reduction can be exercised on a log before a 31-minute leg is
    # spent on it, and so (a3) -- a property of the LOG, not of the binary --
    # can be gated in the packet's run order, ahead of every long leg.
    # ⚠ it re-reduces output written by whatever binary produced it: it is NOT
    # a substitute for re-running after a source change (standing rule 9).
    if which == "reduce":
        tag = sys.argv[2] if len(sys.argv) > 2 else "m13"
        d = os.path.join(OUT, tag)
        print("  REDUCE-ONLY on %s (no run; the binary above is NOT the one"
              % os.path.relpath(d, ROOT))
        print("  that wrote this output unless you just ran it)")
        head, rows = thermo(d)
        ev = events(d)
        print("  firing events: %d" % len(ev))
        return _reduce(head, rows, ev)

    # m = 13 FIRST, then m = 8 (b3 halts the ladder on a frozen floor).
    m = {"m13": 13.0, "m8": 8.0}.get(which)
    if m is None:
        print("  usage: run_ctrl8.py short | m13 | m8")
        return 2
    # ⚠ t_stop from the 60 s leg: v = 1.806 m/h => min(500, 1620/v) = 500 s.
    #   b2 margin 309 mm above the 40 mm bottom guard.
    T_STOP = 500.0
    rc, wall, d = run(which, T_STOP, ["weibull.m=%g" % m])
    print("  rc=%d  wall %.1f s (%.2f wall-s per simulated s)"
          % (rc, wall, wall / T_STOP))
    if rc != 0:
        print("  ABORT -- not a physical result")
        return 1
    head, rows = thermo(d)
    ev = events(d)
    print("  firing events: %d" % len(ev))
    return _reduce(head, rows, ev)


def _reduce(head, rows, ev):
    # ------------------------------------------------------------------
    # P2-2f: the reduction lives in `reduce_price`, imported at the top, so
    # the standalone (a3) gate and this harness cannot drift apart.
    #
    # ⚠⚠ DELIVERABLE 6: NO ROW MAY PRINT A VERDICT ITS PROSE RETRACTS. Rev 1's
    # CTRL8_m13.txt printed "MISS c1" against a reduction the notes withdrew in
    # the same breath, and the withdrawn 1.27 GJ/m^3 then OUTLIVED the
    # retraction of the claim it came from. `verdict()` is the only way a row
    # prints a judgement here, and it refuses to print one without a value.
    # ------------------------------------------------------------------
    def verdict(tag, ok, note=""):
        """The ONLY judgement printer. `ok` None => "not computed"."""
        if ok is None:
            print("     not computed  %s%s" % (tag, ("  -- " + note) if note else ""))
        else:
            print("     %s %s%s" % ("ok  " if ok else "MISS", tag,
                                    ("  -- " + note) if note else ""))

    # b3: THE FROZEN FLOOR, per deliverable 4. ⚠ NOT `foot_carry_cols`.
    fz = RP.frozen_floor(ev)
    if fz is None or fz["frac"] is None:
        print("  b3 frozen floor: no floor firing in the log")
        verdict("b3", None, "no population to score")
        b3_ok = None
    else:
        print("  b3 FROZEN floor (fired once, then silent for the final 30 s):")
        print("     active radius  %.4f m  (%s)" % (fz["r_active"], fz["r_from"]))
        print("     population     %d columns inside that radius that ever fired"
              % fz["n_pop"])
        print("     frozen         %d  =>  %.2f %%" % (fz["n_frozen"],
                                                       100.0 * fz["frac"]))
        print("     excluded: %d columns fired but lie OUTSIDE the radius (hole"
              % (fz["n_fired_total"] - fz["n_pop"]))
        print("       wall -- silent by design once the floor passes them), and")
        print("       every column that never fired (never part of the floor).")
        print("     ⚠ rev 1's 13.6 % is WITHDRAWN: it scored `foot_carry_cols`,")
        print("       which MMWSpalling.H:2884 defines as columns at or above")
        print("       their pad height -- load-bearing contact, max 162 = 100 %.")
        b3_ok = fz["frac"] < 0.10
        verdict("b3  (bar < 10 %)", b3_ok)

    # b4: THE LEDGER. ⚠ `ledger_err_solid`, NOT `ledger_err`.
    les = col(head, rows, "ledger_err_solid")
    le = col(head, rows, "ledger_err")
    fme = col(head, rows, "foot_mech_E")
    if not les:
        print("  b4 ledger_err_solid: column ABSENT -- the leg needs "
              "spall.surface_T=reconstructed and spall.removal_energy=flake")
        verdict("b4", None, "column absent")
        b4_ok = None
    else:
        worst = max(abs(x) for x in les[1:]) if len(les) > 1 else abs(les[0])
        print("  b4 max |ledger_err_solid| = %.3e   (bar <= 5e-3)" % worst)
        if le:
            print("     ⚠ for contrast max |ledger_err| = %.3e, and that is"
                  % max(abs(x) for x in le[1:] or le))
            print("       CORRECT, not a second leak: the TOTAL row never")
            print("       credits enthalpy leaving with removed rock, so it")
            print("       cannot be small on a flake leg. The accepted CTRL5")
            print("       reads -1.7258e-01 on it and 7.78e-04 on the row above.")
        if fme:
            print("     foot_mech_E (P2-2f's pad-clearance sink) = %.4e J"
                  % fme[-1])
        b4_ok = worst <= 5e-3
        verdict("b4", b4_ok)

    # ------------------------------------------------------------------
    # c1: THE PRICE. ⚠⚠ Gated on b4 -- no price may be quoted from a leg whose
    # energy does not close. That gate is the packet's, and deliverable 6 means
    # the row prints "not computed" rather than a MISS it would then retract.
    # ------------------------------------------------------------------
    Ef = col(head, rows, "ledger_E_flake")
    cv = RP.charged_volume(ev)
    print("  c1 charged volume from n_voided, ONE QUADRANT, unscaled:")
    for nm in ("floor", "wall", "pad"):
        print("       %-6s regime %-2s %7d events %7d cells %9.5f L%s"
              % (nm, {"floor": 1, "wall": 6, "pad": 4}[nm],
                 cv[nm]["events"], cv[nm]["cells"], 1e3 * cv[nm]["vol_nvoided"],
                 "   <- EXCLUDED (mechanical)" if nm == "pad" else ""))
    print("       charged = floor + wall = %.5f L" % (1e3 * cv["charged"]))
    ok3, rel3, det3 = RP.a3_cross_check(ev)
    print("     a3 cross-check: %s" % det3)
    verdict("a3  (bar < 0.1 %)", ok3)
    if not Ef:
        verdict("c1", None, "ledger_E_flake absent")
    elif b4_ok is not True:
        print("     ⚠⚠ b4 did not pass, so the price is NOT QUOTED.")
        verdict("c1", None, "b4 gates it")
    elif cv["charged"] <= 0.0:
        verdict("c1", None, "charged volume is zero")
    else:
        # ⚠⚠ NEITHER SIDE IS SCALED. Both the event log and ledger_E_flake are
        # one quadrant, so the x4 cancels in the ratio; rev 1 applied it to the
        # volume only and read 4x low. A PRICE IS A RATIO.
        price = Ef[-1] / cv["charged"] / 1e9
        print("     E_flake %.4e J / %.5f L  =>  price %.3f GJ/m^3  (bar 1.21 "
              "+/- 0.05)" % (Ef[-1], 1e3 * cv["charged"], price))
        print("     ⚠ NUMERATOR CAVEAT: ledger_E_flake is cumulative over FOUR")
        print("       += terms (Removal.H:3102 flake+boundary+support, :3535")
        print("       shift, :3543 lateral) -- ARCHIVE_DONE L18675. A miss of a")
        print("       per cent or two is a numerator question before it is")
        print("       physics. P2-2f's foot_mech_E is a FIFTH sink and is")
        print("       deliberately NOT in this total.")
        print("     ⚠ counting pad clearing too would give %.3f GJ/m^3 -- INSIDE"
              % (Ef[-1] / (cv["charged"] + cv["pad"]["vol_nvoided"]) / 1e9))
        print("       the band. The wrong method yields the better number.")
        verdict("c1  (report, headline)", abs(price - 1.21) <= 0.05)

    # c2: THE PHYSICS ROW -- median reconstructed T_surf AND the clamp fraction,
    # on ONE stated population.
    ts, note = RP.firing_T(ev)
    nc, ntot, fclamp = RP.clamp_fraction(ev)
    if not ts:
        verdict("c2", None, "no firing carries a reconstructed T_surf")
    else:
        med = RP.pctl(ts, 0.5)
        print("  c2 median reconstructed T_surf = %.2f K  (bar 845 +/- 10)" % med)
        print("     %s" % note)
        print("     spread: p75 %.2f  p90 %.2f  p99 %.2f  max %.2f"
              % (RP.pctl(ts, 0.75), RP.pctl(ts, 0.90), RP.pctl(ts, 0.99),
                 max(ts)))
        print("     ⚠⚠ ABOVE THE 873.15 K TABLE CLAMP: %d of %d = %.2f %%."
              % (nc, ntot, 100.0 * fclamp))
        print("       Both Grimsel tables repeat their last value at 600 and")
        print("       650 C, so the last real interval is 823-873 K and NOTHING")
        print("       above 873.15 K is scored by either property law. The")
        print("       median is referred to Kant's 843.757 K with nothing")
        print("       refitted, but %.0f %% of the population sits outside the"
              % (100.0 * fclamp))
        print("       tables' validity -- both numbers or the row overclaims.")
        verdict("c2  (report, THE PHYSICS ROW)", abs(med - 845.0) <= 10.0,
                "clamp fraction %.2f %% stated alongside" % (100.0 * fclamp))

    # c3: BARE REPORT -- no bar. ⚠ standing rule 3: a volume reported ON ITS
    # OWN carries the x4, unlike c1's ratio where it cancels.
    pm = col(head, rows, "patch_min_standoff")
    print("  c3 charged volume x4 (full domain) = %.4f L vs Meier 1.270 L"
          % (1e3 * 4.0 * cv["charged"]))
    print("     ⚠ the x4 IS applied here and is NOT applied in c1: c1 is a")
    print("       ratio of two quarter-domain quantities, where it cancels.")
    if pm:
        print("     standoff median %.4f m (feet pin it; the archived 2.27x ran"
              % RP.pctl([x for x in pm if x > 0], 0.5))
        print("       out to ~0.200 m, h(axis) 1449 vs 589.8 = 2.46x), so a")
        print("       ratio to 2.27x would be CROSS-MODE.")
    verdict("c3", None, "bare report by specification -- no bar to score")

    # c4: rate, depth-matched 150-350 s per the packet.
    tt, nz = col(head, rows, "time"), col(head, rows, "nozzle_z")
    win = [(t, z) for t, z in zip(tt, nz) if 150.0 <= t <= 350.0 and z > 0]
    if len(win) > 1 and win[-1][0] > win[0][0]:
        rate = (win[0][1] - win[-1][1]) / (win[-1][0] - win[0][0]) * 3600.0
        print("  c4 rate over %.0f-%.0f s = %.3f m/h  (Meier 1.5 +/- 20 %% = "
              "1.2-1.8)" % (win[0][0], win[-1][0], rate))
        print("     ⚠ PROVISIONAL: 1 mm is parked (rule 12), so no ABSOLUTE 2 mm")
        print("       rate is quotable; the firing skin is ~0.35 cells at 2 mm")
        print("       under Grimsel conductivity.")
        verdict("c4  (report, PROVISIONAL)", 1.2 <= rate <= 1.8)
    else:
        verdict("c4", None, "no 150-350 s descent window in thermo.dat")
    return 0


if __name__ == "__main__":
    sys.exit(main())
