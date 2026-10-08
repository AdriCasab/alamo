#!/usr/bin/env python3
"""P1g-1 (d1): the UNCAPPED lateral flake thickness, c_f * ell_lat.

⚠⚠ WHY THIS IS COMPUTED AND NOT READ. The new `h_f_unc_lat` event-log column
was added BY THIS STEP, so the archived leg it must reproduce
(`output_legB_fix`, written before the column existed) does not carry it. The
quantity is therefore recomputed from the columns the archive DOES carry --
`ell_lat` -- times the leg's own `spall.flake_c_f`. That is the same arithmetic
the kernel performs at `Removal.H:750` immediately before the cap clips it:

    h_f      = c_f * ell            <-- what this script reproduces
    h_f_lat  = min(h_f, flake_h_max) <-- what the archive logged

⚠ `c_f` and the cap are read from the leg's OWN `argv.sh`, never transcribed
here. Standing rule 8: a number is backed only if both the file it was read from
and the code that reduced it are in the repository.

⚠⚠ WHY IT MATTERS. The cap is what makes the lateral flake sub-cell at all: the
uncapped median is ABOVE dx, so without `spall.flake_h_max` the "sub-cell"
recession this step implements would not be sub-cell. A later step re-derives
that bound from evidence, and this column is the evidence. Do NOT touch the cap
in this step.
"""
import csv, math, os, statistics, sys

HERE = os.path.dirname(os.path.abspath(__file__))
LEG  = os.path.join(HERE, "output_legB_fix")
EV   = os.path.join(LEG, "FH_removal_events.csv")
ARGV = os.path.join(LEG, "argv.sh")

# the packet's (d1) reproduction target
TGT_P50_MM, TOL_MM = 2.465, 0.05


def read_argv(path, key):
    """the leg's own argv, so the constant is the one that run used."""
    if not os.path.exists(path):
        raise SystemExit("FAIL: no argv.sh at %s -- the constant cannot be "
                         "backed to the run that produced the data" % path)
    for tok in open(path).read().split():
        if tok.startswith(key + "="):
            return float(tok.split("=", 1)[1])
    raise SystemExit("FAIL: %s not found in %s" % (key, path))


def main():
    P = print
    if not os.path.exists(EV):
        raise SystemExit("FAIL: no event log at %s" % EV)
    c_f   = read_argv(ARGV, "spall.flake_c_f")
    h_max = read_argv(ARGV, "spall.flake_h_max")

    ell, capped = [], 0
    for r in csv.DictReader(open(EV)):
        try:
            if int(float(r["regime"])) != 6:
                continue
            e  = float(r["ell_lat"])
            hf = float(r["h_f_lat"])
        except (KeyError, ValueError):
            continue
        if not (e > 0.0):
            continue
        ell.append(e)
        if hf >= h_max * (1.0 - 1e-9):
            capped += 1

    if not ell:
        raise SystemExit("FAIL: no regime-6 events with ell_lat > 0 in %s" % EV)

    unc = [c_f * e for e in ell]
    n   = len(unc)
    q   = lambda p: float(statistics.quantiles(unc, n=100)[p - 1]) if n > 1 \
                    else unc[0]
    p50 = float(statistics.median(unc))

    P("=" * 78)
    P("P1g-1 (d1) THE UNCAPPED LATERAL FLAKE THICKNESS  c_f * ell_lat")
    P("=" * 78)
    P("  leg        : %s" % os.path.relpath(LEG, HERE))
    P("  constants  : READ from that leg's own argv.sh --")
    P("               spall.flake_c_f   = %g" % c_f)
    P("               spall.flake_h_max = %g m  (= %.3f mm)" % (h_max, h_max * 1e3))
    P("  n regime-6 events with ell_lat > 0 : %d" % n)
    P("  at the cap (h_f_lat == flake_h_max): %d of %d = %.2f %%"
      % (capped, n, 100.0 * capped / n))
    P("")
    P("  UNCAPPED c_f * ell_lat [mm]:")
    for lbl, v in (("p05", q(5)), ("p50", p50), ("p95", q(95)),
                   ("max", max(unc))):
        P("    %-4s %10.4f" % (lbl, v * 1e3))
    P("")
    # dx is a property of the leg's mesh, stated by the packet's Sources table
    DX_MM = 2.0
    P("  uncapped p50 / dx = %.4f mm / %.1f mm = %.3f cells"
      % (p50 * 1e3, DX_MM, p50 * 1e3 / DX_MM))
    P("  capped   p50 / dx = %.4f mm / %.1f mm = %.3f cells"
      % (h_max * 1e3, DX_MM, h_max * 1e3 / DX_MM))
    P("  ⚠⚠ the UNCAPPED median exceeds dx, so the cap is what makes the lateral")
    P("     flake sub-cell. It is NOT touched by this step.")
    P("")
    P("-" * 78)
    ok = abs(p50 * 1e3 - TGT_P50_MM) <= TOL_MM
    P("  (d1) reproduction target p50 = %.3f mm +- %.3f" % (TGT_P50_MM, TOL_MM))
    P("       measured                 = %.4f mm   -> %s"
      % (p50 * 1e3, "PASS" if ok else "FAIL"))
    if not ok:
        P("  ⚠⚠ a disagreement here is a defect in THIS SCRIPT or a change in the")
        P("     archived leg, not a finding about the physics.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
