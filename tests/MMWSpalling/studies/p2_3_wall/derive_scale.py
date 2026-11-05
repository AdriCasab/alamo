#!/usr/bin/env python3
"""P2-3 (a1) -- the floor coefficient scaled to Kant, derived and round-tripped.

THE ARGUMENT. P2-2f left the rock half closed and the residual in DELIVERY. The
energy partition said the delivery is COEFFICIENT-limited, and the model's floor
map is high: its axis value at the m13 leg's own standoff is 1,326.5 W/m^2K
against 985.0 interpolated from Kant & von Rohr 2016 Table 4 -- the only
measured coefficients on a same-lab methane-air burner. So the map is scaled by
985.0/1326.5.

⚠⚠ THE SCALE IS TAKEN AT THE LEG'S OWN STANDOFF, NOT AT THE LITERAL. The map's
peak literal is `1449.0`, which is its value at s = 50 mm, NOT at the 85.3 mm
the feet actually hold. Dividing Kant's 985 by 1449 would be comparing two
different standoffs and would over-scale by 8.5 %. This script prints the map at
BOTH standoffs so that cannot be fudged.

⚠⚠ ONE LITERAL CHANGES, AND THE PARSE-BACK PROVES IT. The map is a 1,400-
character expression with `142.0` appearing twice and radial blends either side
of the peak; an edit in two places would be invisible by eye. The scaled string
is diffed against the original token by token.

⚠ THE +/- 0.01 IN THE BAR IS AN ARITHMETIC CHECK, NOT A PRECISION CLAIM. Table 4's
plate carries ~17 % sensor-calibration uncertainty, so the PHYSICAL band on this
scale is 0.62-0.87. Nothing here is fitted to a Meier number: if the pad ring
starves (b1), the finding is the map's SHAPE, not a different scale.

⚠ Rule 17: the plate's TEMPERATURE was retired as a target because the plate
replaced the rock. Its COEFFICIENT is a property of the FLOW, which is why this
use is legitimate -- a confined flame-jet calorimeter study found h unchanged
over target temperatures 170-530 C.

Usage: derive_scale.py > SCALE.txt
"""
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
ARGV = os.path.join(ROOT, "tests", "MMWSpalling", "studies", "p2_2c_meier",
                    "feet_argv.txt")
KANT = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "kant",
                    "kant_von_rohr_validation_test_for_mmwspalling_plan.txt")

# Kant & von Rohr 2016 Table 4, PARSED from the reference file rather than
# retyped -- a hand-copied validation constant is how a wrong number survives.
SOD_KEYS = (50, 80, 110, 130, 150)


def kant_table():
    rows = {}
    for ln in open(KANT):
        m = re.match(r"\s*(50|80|110|130|150)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*$", ln)
        if m:
            sod = int(m.group(1))
            rows[sod] = {"q_sens": float(m.group(2)),
                         "theta_s": float(m.group(3)),
                         "T_flame_C": float(m.group(4)),
                         "h_fl": float(m.group(5))}
    return rows


def kant_h(rows, s_m):
    """h_fl linearly interpolated at standoff s_m [m]."""
    xs = [k * 1e-3 for k in SOD_KEYS]
    ys = [rows[k]["h_fl"] for k in SOD_KEYS]
    if s_m <= xs[0] or s_m >= xs[-1]:
        return None
    for i in range(len(xs) - 1):
        if xs[i] <= s_m <= xs[i + 1]:
            t = (s_m - xs[i]) / (xs[i + 1] - xs[i])
            return ys[i] + t * (ys[i + 1] - ys[i])
    return None


def h_expr_string():
    txt = open(ARGV).read()
    m = re.search(r"surface_patch\.h_expr=(\S+)", txt)
    if not m:
        raise SystemExit("no surface_patch.h_expr in %s" % ARGV)
    return m.group(1)


def ev(e, r, s):
    """Evaluate ALAMO's h_expr at radius r [m], standoff s [m].
    ⚠ `if(c,a,b)` is ALAMO's ternary, not Python's; both branches are safe to
    evaluate here because every divisor is guarded by a max()."""
    def _if(c, a, b):
        return a if c else b
    py = e.replace("if(", "_if(").replace("pow(", "_pow(")
    return eval(py, {"_if": _if, "_pow": math.pow, "max": max, "min": min,
                     "cos": math.cos, "r": r, "s": s})


def median_standoff():
    """The m13 leg's own median standoff, from its thermo, not retyped."""
    p = os.path.join(ROOT, "tests", "MMWSpalling", "studies", "p2_2c_meier",
                     "ctrl8", "m13", "run", "thermo.dat")
    if not os.path.exists(p):
        return None
    lines = open(p).read().splitlines()
    head = lines[0].split()
    if "patch_min_standoff" not in head:
        return None
    i = head.index("patch_min_standoff")
    vals = []
    for ln in lines[1:]:
        f = ln.split()
        if len(f) > i:
            try:
                v = float(f[i])
            except ValueError:
                continue
            if v > 0.0:
                vals.append(v)
    if not vals:
        return None
    vals.sort()
    return vals[len(vals) // 2]


def main():
    fails = []

    def need(ok, msg):
        print("  %s  %s" % ("ok  " if ok else "FAIL", msg))
        if not ok:
            fails.append(msg)

    print("P2-3 (a1) -- THE FLOOR SCALE, DERIVED")
    print("=" * 72)
    rows = kant_table()
    print("\nKant & von Rohr 2016 Table 4, parsed from")
    print("  %s" % os.path.relpath(KANT, ROOT))
    print("   SOD[mm]  h_fl[W/m2K]  T_flame[K]  q_sens[kW/m2]")
    for k in SOD_KEYS:
        print("   %6d  %11.1f  %10.2f  %13.1f"
              % (k, rows[k]["h_fl"], rows[k]["T_flame_C"] + 273.15,
                 rows[k]["q_sens"]))
    need(len(rows) == 5, "all five Table 4 rows parsed")

    s_med = median_standoff()
    if s_med is None:
        print("\n  not computed -- no m13 thermo to read the standoff from")
        return 1
    print("\nThe leg's own standoff (median patch_min_standoff, m13): %.6f m"
          % s_med)

    hk = kant_h(rows, s_med)
    expr = h_expr_string()
    axis_at_leg = ev(expr, 0.0, s_med)
    axis_at_50 = ev(expr, 0.0, 0.050)
    print("\nKant interpolated at that standoff : %10.3f W/m2K" % hk)
    print("map axis value at that standoff    : %10.3f W/m2K" % axis_at_leg)
    print("map axis value at s = 50 mm        : %10.3f W/m2K  <- the LITERAL"
          % axis_at_50)
    need(abs(axis_at_50 - 1449.0) < 1e-6,
         "the peak literal IS the map's 50 mm axis value (1449.0)")

    scale = hk / axis_at_leg
    print("\nscale = %.3f / %.3f = %.6f" % (hk, axis_at_leg, scale))
    need(abs(scale - 0.743) <= 0.01,
         "scale = 0.743 +/- 0.01 (arithmetic check; PHYSICAL band 0.62-0.87 "
         "from Table 4's ~17 %% plate uncertainty)")
    # ⚠ had the scale been taken on the literal instead:
    print("  ⚠ for contrast, Kant/1449.0 = %.6f -- taking the scale on the"
          % (hk / axis_at_50))
    print("    LITERAL rather than at the leg's standoff would over-scale by")
    print("    %.1f %%, because 1449.0 is the map at 50 mm, not at %.1f mm."
          % (100.0 * (1.0 - (hk / axis_at_50) / scale), 1e3 * s_med))

    new_peak = 1449.0 * scale
    print("\nnew peak literal = 1449.0 x %.6f = %.3f" % (scale, new_peak))

    # ---- the scaled string, ONE literal changed, and the proof ----
    new_expr, n_sub = re.subn(re.escape("1449.0*"), "%.1f*" % new_peak, expr, count=1)
    need(n_sub == 1, "exactly one substitution was made")
    need(expr.count("1449.0") == 1,
         "`1449.0` occurs ONCE in the original, so there is no second peak to miss")
    # token-by-token diff
    tok_a = re.findall(r"[0-9]*\.?[0-9]+(?:[eE][+-]?[0-9]+)?|[A-Za-z_]+|.", expr)
    tok_b = re.findall(r"[0-9]*\.?[0-9]+(?:[eE][+-]?[0-9]+)?|[A-Za-z_]+|.", new_expr)
    diffs = [(i, a, b) for i, (a, b) in enumerate(zip(tok_a, tok_b)) if a != b]
    print("\nparse-back: %d tokens before, %d after, %d differing"
          % (len(tok_a), len(tok_b), len(diffs)))
    for i, a, b in diffs:
        print("   token %4d : %r -> %r" % (i, a, b))
    need(len(tok_a) == len(tok_b) and len(diffs) == 1
         and diffs[0][1] == "1449.0" and diffs[0][2] == "%.1f" % new_peak,
         "the scaled map differs from the original in EXACTLY ONE TOKEN, and it "
         "is the peak literal")
    need(abs(ev(new_expr, 0.0, s_med) - hk) < 0.5,
         "the scaled map reads Kant's %.1f at the leg's standoff (%.2f)"
         % (hk, ev(new_expr, 0.0, s_med)))
    # the ring and the far field must NOT have moved
    need(abs(ev(expr, 0.050, s_med) - ev(new_expr, 0.050, s_med)) > 0.0,
         "the ring blend region DOES change (it mixes the peak) -- stated, not hidden")
    need(abs(ev(expr, 0.080, s_med) - ev(new_expr, 0.080, s_med)) < 1e-12,
         "beyond 76 mm the map is UNCHANGED (the 0.001 far field is untouched)")

    print("\nradial profile at the leg's standoff, before -> after:")
    print("    r[mm]     before      after     ratio")
    for rr in (0.0, 0.010, 0.020, 0.028, 0.035, 0.040, 0.050, 0.060, 0.080):
        b = ev(expr, rr, s_med)
        a = ev(new_expr, rr, s_med)
        print("   %6.1f  %9.2f  %9.2f    %.4f"
              % (1e3 * rr, b, a, (a / b) if b else float("nan")))
    print("  ⚠ the peak scales by %.3f and the 142.0 ring plateau does NOT, so"
          % scale)
    print("    the profile FLATTENS -- the direction D2ad's measured floor")
    print("    (flat 0.48-0.60 MW/m^2 out to 35 mm) asks for. Scaling the whole")
    print("    string instead would have preserved the 2:1 peak.")

    # ---- the pad ring, which is what b1 risks ----
    T_G, T_S = 1504.0, 850.0
    print("\nPAD-RING FLUX (the b1 risk, D2ad). T_gas %.0f K, T_s %.0f K:" % (T_G, T_S))
    print("    r[mm]   h[W/m2K]   q[MW/m2]")
    for rr in (0.028, 0.032, 0.035, 0.040):
        a = ev(new_expr, rr, s_med)
        print("   %6.1f   %8.2f   %8.4f" % (1e3 * rr, a, a * (T_G - T_S) / 1e6))
    q35 = ev(new_expr, 0.035, s_med) * (T_G - T_S)
    print("  => at r = 35 mm, q = %.4f MW/m^2 against Kant's threshold class"
          % (q35 / 1e6))
    print("     0.45-0.60 MW/m^2  ==>  %s"
          % ("MARGINAL -- BELOW the class" if q35 / 1e6 < 0.45 else "inside"))
    print("  ⚠ PRE-REGISTERED (P5): if b1's pads stall, the finding is the map's")
    print("    SHAPE and it is queued as P2-4. THE SCALE IS NOT REVERTED.")

    print("\nTHE SCALED h_expr (for the leg argv, quoting preserved):")
    print("surface_patch.h_expr=%s" % new_expr)

    print("\n" + "=" * 72)
    if fails:
        print("FAIL %d: %s" % (len(fails), "; ".join(fails)))
        return 1
    print("PASS -- scale %.6f, peak literal 1449.0 -> %.1f" % (scale, new_peak))
    return 0


if __name__ == "__main__":
    sys.exit(main())
