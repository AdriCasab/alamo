#!/usr/bin/env python3
"""P1g-0b (d1): the firing margin the coefficient fix opens, PER DEPTH BAND.

⚠⚠ WHY THIS FILE IS IN THE REPOSITORY AND NOT A SCRATCHPAD. Revision 1's version
lived in a session scratchpad, so no reviewer could open it -- and it HARD-CODED
`T_gas = 1200 K`, which is the UNIT FIXTURE's constant
(`unit/lateral_drive/input`: `surface_patch.side_face_T = 1200.0`). Leg B is a
Meier run with `side_face_T = live`, where the gas temperature is the jet exhaust
and runs 1464-1599 K. That one constant moved the requirement ~45 K and REVERSED
the conclusion. Standing rule 8: a number is backed only if both the file it was
read from and the code that reduced it are in the repository.

⚠⚠ NOTHING THE RUN LOGS IS HARD-CODED HERE. The balance

    G (Tf - Tc) = h (Tg - Tf) - eps sig (Tf^4 - Ta^4) - h_c (Tf - Ta)

is exact and has ONE unknown once the logged TRIPLE is in hand. Each regime-6
event logs `T_face_lat` (= Tf, the face the criterion judged) and `T_top` (= Tc,
the cell it read); `h` was the 80.0 W/m^2K cap in force when that run was made.
So the GAS TEMPERATURE IS BACKED OUT, not assumed:

    Tg = Tf + [ G (Tf - Tc) + loss(Tf) ] / h

Then the same balance is re-solved at the CONSISTENT h = `side_face_h` for the
cell temperature that reaches the same threshold face. ⚠ The material constants
(kappa, the loss pair, T_ref) are PARSED from the Meier input, not transcribed.

⚠⚠ POOLING IS FORBIDDEN. The gas runs ~1464 K at the mouth to ~1599 K deep and
the vertical spread is broad shallow and narrow deep, so a single pooled number
would hide that the MOUTH is roughly a coin-flip race while the DEEP HOLE is a
walkover -- and the deep hole is exactly where P1g's prerequisite bites.
"""
import collections, csv, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
LEG = os.path.join(HERE, "output_legB_fix")
INP = os.path.join(ROOT, "MMWSpalling", "validation", "meier",
                   "sp_meier_pilot", "input_drilling")
if not os.path.exists(INP):           # tests/ layout differs by checkout depth
    INP = os.path.join(os.path.dirname(os.path.dirname(HERE)),
                       "validation", "meier", "sp_meier_pilot", "input_drilling")

H_CAPPED = 80.0        # surface_patch.side_face_h_max default -- in force for
                       # output_legB_fix, which predates the cap fix
DBANDS = ((0, 100), (100, 200), (200, 300), (300, 400))
RING = (40.0, 60.0)    # mm, the wall ring P1g-0 scored
KTOP0, DZ, DX = 299, 0.002, 0.002


def parse_input(path):
    """the material constants, PARSED -- not transcribed into this file."""
    want = {"microstructure.phase0.kappa": None,
            "microstructure.phase0.T_ref": None,
            "losses.epsilon": None, "losses.epsilon_high": None,
            "losses.T_eps_break": None, "losses.sigma_SB": None,
            "losses.h_conv": None, "losses.T_amb": None}
    for ln in open(path):
        ln = ln.split("#")[0]
        if "=" not in ln:
            continue
        k, v = ln.split("=", 1)
        k = k.strip()
        if k in want:
            try:
                want[k] = float(v.strip())
            except ValueError:
                pass
    missing = [k for k, v in want.items() if v is None]
    if missing:
        raise SystemExit("FAIL: could not parse %s from %s" % (missing, path))
    return want


def main():
    P = print
    cst = parse_input(INP)
    KAPPA = cst["microstructure.phase0.kappa"]
    T_REF = cst["microstructure.phase0.T_ref"]
    EPS = cst["losses.epsilon"]
    SIG = cst["losses.sigma_SB"]
    HC = cst["losses.h_conv"]
    TA = cst["losses.T_amb"]
    # ⚠ `side_face_h` is NOT in the input file -- it is set on the argv. Read it
    # from the run's OWN committed argv, so the number is the one that leg used
    # rather than one transcribed here.
    H_TRUE = None
    av = os.path.join(LEG, "argv.sh")
    if os.path.exists(av):
        for tok in open(av).read().split():
            if tok.startswith("surface_patch.side_face_h="):
                H_TRUE = float(tok.split("=", 1)[1])
    if H_TRUE is None:
        raise SystemExit("FAIL: no surface_patch.side_face_h in %s -- the "
                         "consistent h cannot be read from the run itself" % av)
    G = 2.0 * KAPPA / DX

    P("=" * 78)
    P("P1g-0b (d1) THE FIRING MARGIN, PER DEPTH BAND -- from the LOGGED TRIPLE")
    P("=" * 78)
    P("  leg        : %s" % os.path.relpath(LEG, ROOT))
    P("  constants  : PARSED from %s" % os.path.relpath(INP, ROOT))
    P("               kappa %.4g W/mK  => G = 2k/dx = %.1f W/m2K at dx = %g mm"
      % (KAPPA, G, DX * 1e3))
    P("               eps %.3g  sigma %.6g  h_conv %.4g  T_amb %.2f  T_ref %.2f"
      % (EPS, SIG, HC, TA, T_REF))
    P("               side_face_h (CONSISTENT h) = %.4g, from the run's own argv.sh" % H_TRUE)
    P("  h in force when this leg ran = %.1f (the side_face_h_max cap)" % H_CAPPED)
    P("  ⚠⚠ the GAS TEMPERATURE IS BACKED OUT of the triple, never assumed.")

    def loss(Tf):
        return EPS * SIG * (Tf ** 4 - TA ** 4) + HC * (Tf - TA)

    def gas_from_triple(Tf, Tc, h):
        return Tf + (G * (Tf - Tc) + loss(Tf)) / h

    def cell_for_face(Tf, Tg, h):
        """invert for the CELL temperature that drives the face to Tf at this h."""
        return Tf - (h * (Tg - Tf) - loss(Tf)) / G

    ev = os.path.join(LEG, "FH_removal_events.csv")
    if not os.path.exists(ev):
        raise SystemExit("FAIL: no event log at %s" % ev)
    rows = list(csv.DictReader(open(ev)))

    lat = collections.defaultdict(list)   # band -> (Tf, Tc, ell, a_f)
    ver = collections.defaultdict(list)   # band -> T_top
    for r in rows:
        try:
            reg = int(float(r["regime"]))
            i, j = int(r["col_i"]), int(r["col_j"])
            kt = int(float(r["k_top"]))
            rd = float(r.get("riser_depth") or 0.0)
        except (KeyError, ValueError):
            continue
        rr = DX * math.hypot(i + 0.5, j + 0.5) * 1e3
        if not (RING[0] <= rr < RING[1]):
            continue
        depth = (KTOP0 - kt) * DZ * 1e3 + rd * 1e3
        b = next((n for n, (a, c) in enumerate(DBANDS) if a <= depth < c), None)
        if b is None:
            continue
        if reg == 6:
            try:
                Tf = float(r["T_face_lat"]); Tc = float(r["T_top"])
                el = float(r["ell_lat"]); af = float(r["a_f"])
            except (KeyError, ValueError):
                continue
            if Tf > Tc > T_REF and el > 0:
                lat[b].append((Tf, Tc, el, af))
        elif reg == 1:
            ver[b].append(float(r["T_top"]))

    P("")
    P("-" * 78)
    P("  ring %g-%g mm. n6 = regime-6 events, n1 = regime-1 (vertical) events."
      % RING)
    P("  %-10s %5s %9s %9s %9s %9s %7s %9s %9s %8s"
      % ("band[mm]", "n6", "Tface[K]", "Tcell[K]", "gas[K]", "need[K]",
         "n1", "vert p10", "vert p50", "lead/p10"))
    out = []
    for b, (d0, d1) in enumerate(DBANDS):
        L, V = lat[b], ver[b]
        if not L:
            P("  %-10s %5d   -- no regime-6 events in this band --"
              % ("%d-%d" % (d0, d1), 0))
            continue
        Tf = float(np.median([x[0] for x in L]))
        Tc = float(np.median([x[1] for x in L]))
        # ⚠ the gas is backed out PER EVENT and then reduced, not reduced first:
        # the balance is nonlinear in Tf, so median-then-invert != invert-then-median.
        Tg = float(np.median([gas_from_triple(x[0], x[1], H_CAPPED) for x in L]))
        need = cell_for_face(Tf, Tg, H_TRUE)
        p10 = float(np.percentile(V, 10)) if V else float("nan")
        p50 = float(np.percentile(V, 50)) if V else float("nan")
        P("  %-10s %5d %9.1f %9.1f %9.0f %9.1f %7d %9.0f %9.0f %+8.1f"
          % ("%d-%d" % (d0, d1), len(L), Tf, Tc, Tg, need, len(V), p10, p50,
             p10 - need))
        out.append((d0, d1, Tg, need, p10, p50))

    P("")
    P("  'need' = the CELL temperature that reaches the SAME threshold face once")
    P("  the criterion uses the kernel's own h = %.4g instead of the %.0f cap."
      % (H_TRUE, H_CAPPED))
    P("  'lead/p10' = vertical p10 MINUS that requirement: POSITIVE means the")
    P("  lateral mode reaches threshold BELOW the coolest decile at which the")
    P("  vertical mode removes, i.e. it wins essentially every race there.")
    P("")
    P("-" * 78)
    P("  ⚠⚠ REPRODUCTION TARGETS (packet P3 -- a disagreement is a defect in")
    P("     THIS SCRIPT, not a finding about the physics):")
    P("       gas  1464 / 1562 / 1590 / 1599 K")
    P("       need 689.5 / 678.7 / 674.7 / 674.3 K")
    P("       vert p10 605 / 689 / 703 / 706 K")
    P("  measured:")
    P("       gas  %s" % " / ".join("%.0f" % o[2] for o in out))
    P("       need %s" % " / ".join("%.1f" % o[3] for o in out))
    P("       vert p10 %s" % " / ".join("%.0f" % o[4] for o in out))
    tgt_gas = (1464.0, 1562.0, 1590.0, 1599.0)
    tgt_need = (689.5, 678.7, 674.7, 674.3)
    tgt_p10 = (605.0, 689.0, 703.0, 706.0)
    bad = []
    for n, o in enumerate(out):
        if n >= len(tgt_gas):
            break
        if abs(o[2] - tgt_gas[n]) > 25.0:
            bad.append("gas band %d: %.0f vs %.0f" % (n, o[2], tgt_gas[n]))
        if abs(o[3] - tgt_need[n]) > 5.0:
            bad.append("need band %d: %.1f vs %.1f" % (n, o[3], tgt_need[n]))
        if abs(o[4] - tgt_p10[n]) > 8.0:
            bad.append("vert p10 band %d: %.0f vs %.0f" % (n, o[4], tgt_p10[n]))
    if bad:
        P("  ⚠⚠ DISAGREES with the reproduction targets:")
        for x in bad:
            P("       %s" % x)
    else:
        P("  ✅ reproduces every target within tolerance (gas 25 K, need 5 K, "
          "p10 8 K)")
    P("")
    P("  ⚠⚠ POOLED NUMBERS ARE NOT PRINTED, DELIBERATELY. The mouth and the deep")
    P("     hole are different races; one number hides which.")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
