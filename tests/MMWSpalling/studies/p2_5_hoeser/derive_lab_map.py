#!/usr/bin/env python3
"""P2-5 rev 2, deliverable 2 -- the lab h map, derived OFFLINE.

Nothing is read from a run except `jet_D_e` from `probe/SMOKE0` (the flat-h
closure probe, which is NOT a result and is quoted for that one number only).

  (i)   T_stag(L) for K = 20 and 30, from the closure formula
  (ii)  h0 = her 7,100 W/m2K and the fitted standoff exponent n per K, with the
        flux residual at each of L 10/11/13/14
  (iii) the lid constant per leg, s_lid = (nozzle_z0 - 0.10) + dz/2
  (iv)  per shape: f_disc, h0_shape, q0 at 850 K, P over the Ø 40 opening,
        the ratio to the 0.55-0.68 kW flake bound, r_0.5
  (v)   one h_expr string per leg + the BRANCH WITNESS
  (vi)  (run_p2_5.py --probe does the parse-back)

⚠ THE REFERENCE FRAME (rev 2's decision). Her h is SENSOR-REFERENCED:
q = h*(T_f - T_s) with T_f = 1,866 K her nozzle-exit temperature. The model's
Robin form is q = h*(T_stag - T_face). The two coincide at the calibration point
only if phi = 1 there, i.e. jet_core_length*D_e >= s0 -- which is what K30 and
K20 buy. Rev 1 failed because jet_core_length = 8 put the core at 3.8 D and the
closure had already cooled the gas below her sensing face beyond L 11.
⚠ No number here is adjusted to her removed VOLUME (packet guardrail).
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from p2_5_legs import (D_NOZ, S0, T_NOZ, T_ENT, H0, R_DISC, R_OPEN, T_FACE_REF,
                       LID_H, FIT_L, SHAPES, LEGS, ORDER, s_lid, s_leg)

ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
PROBE = os.path.join(HERE, "probe", "SMOKE0", "run", "thermo.dat")
HT = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "hoeser",
                  "hoeser_2016_heat_transfer.csv")
FLAKE_LO, FLAKE_HI = 550.0, 680.0      # W of flakes, 26 cm3/min at 1.27-1.56 GJ/m3

# The planner's pre-computed table (ACTIVE_STEP.md "The calibration,
# pre-computed"), reproduced here so a1 can assert against it rather than
# against prose.
WANT_N = {30: 1.319, 20: 0.429}
WANT_TSTAG = {30: {10: 1866, 11: 1866, 13: 1866, 14: 1866},
              20: {10: 1784, 11: 1649, 13: 1440, 14: 1358}}
WANT_RES = {30: {10: -7.2, 11: -6.7, 13: 1.2, 14: 6.1},
            20: {10: 4.6, 11: -1.0, 13: -1.3, 14: -0.1}}
WANT_P = {"A": 9060.0, "G10": 2940.0, "G075": 1830.0, "G05": 960.0,
          "G035": 630.0}
WANT_R05 = {"A": None, "G10": 18.5, "G075": 14.0, "G05": 9.6, "G035": 7.1}


def read_De():
    """jet_D_e as the BINARY printed it under jet_T_ent_mode = fixed."""
    L = open(PROBE).read().split("\n")
    hdr = L[0].split()
    rows = [l.split() for l in L[1:] if l.strip()]
    return float(rows[-1][hdr.index("jet_D_e")])


def her_table():
    """L -> (q [W/m2], h [W/m2K], T_sens [K]) at 38 kW, from her CSV."""
    out = {}
    for ln in open(HT):
        if ln.startswith("#") or ln.startswith("L,"):
            continue
        p = ln.strip().split(",")
        if len(p) < 4 or not p[3] or int(p[1]) != 38:
            continue
        L, q, h = int(p[0]), float(p[2]) * 1e6, float(p[3]) * 1e3
        out[L] = (q, h, T_NOZ - q / h)
    return out


def T_stag(s, K, De):
    phi = min(1.0, K * De / s)
    return T_ENT + (T_NOZ - T_ENT) * phi, phi


def fit_n(K, De, her):
    """n = sum(x*y)/sum(x^2), x = ln(s0/s), y = ln(h_req/H0), over FIT_L."""
    sxy = sxx = 0.0
    rows = []
    for L in FIT_L:
        s = L * D_NOZ
        q, _, Ts = her[L]
        Tg, phi = T_stag(s, K, De)
        hreq = q / (Tg - Ts)
        x = math.log(S0 / s)
        y = math.log(hreq / H0)
        sxy += x * y
        sxx += x * x
        rows.append((L, s, phi, Tg, Ts, q, hreq))
    return sxy / sxx, rows


def sigma_m(shape):
    """⚠⚠ SHAPES holds sigma in units of D_NOZ. Converting it in ONE place is
    deliberate: rev 2's first draft passed the unit-of-D value straight into
    f_disc and 2*sigma^2, which silently emitted h0_shape = 7,100 for EVERY
    shape (f_disc ~ 1 at sigma = 0.5 "m") and a Gaussian 160x too wide. The
    table half had converted correctly, so the two disagreed -- which is the
    only reason it was caught. Never let a unit conversion happen at two sites."""
    sig = SHAPES[shape]
    return None if sig is None else sig * D_NOZ


def h0_of(shape):
    """The shape's axis coefficient: her disc-MEAN 7,100 W/m2K deconvolved to a
    peak, since her Ø 9 mm sensor reports the mean over its own face."""
    return H0 / f_disc(sigma_m(shape))


def f_disc(sig):
    """Mean/peak of a Gaussian over her Ø 9 mm sensing disc:
    f = (2 sig^2/R^2)(1 - exp(-R^2/2 sig^2)).  Plateau (A) -> 1."""
    if sig is None:
        return 1.0
    a = R_DISC * R_DISC / (2.0 * sig * sig)
    return (1.0 / a) * (1.0 - math.exp(-a))


def P_open(sig, h0s):
    """Delivered power over the Ø 40 opening at T_face = 850 K, full circle."""
    dT = T_NOZ - T_FACE_REF
    if sig is None:
        return h0s * dT * math.pi * R_OPEN * R_OPEN
    a = R_OPEN * R_OPEN / (2.0 * sig * sig)
    return h0s * dT * 2.0 * math.pi * sig * sig * (1.0 - math.exp(-a))


def r_half(sig, h0s):
    """r where q(850 K) falls below 0.5 MW/m2 (the c4 periphery boundary)."""
    h_t = 0.5e6 / (T_NOZ - T_FACE_REF)
    if sig is None:
        return None
    if h0s <= h_t:
        return 0.0
    return sig * math.sqrt(2.0 * math.log(h0s / h_t))


def g(v):
    """A decimal literal the parser will read, with no exponent and no space."""
    return f"{v:.10f}".rstrip("0").rstrip(".") if abs(v) < 1e7 else f"{v:.6f}"


def h_expr(leg, n_of_K):
    """⚠⚠ NO WHITESPACE ANYWHERE. ParmParse re-tokenises a value on whitespace
    and a std::string query keeps only the FIRST token, so `... and ...` would
    truncate the expression into a syntax error (the P2-2f h_expr trap). `and`
    is safe written as `)and(`: the lexer's "and" rule (amrex_parser.l:68)
    precedes the SYMBOL rule, so it wins flex's longest-match tie."""
    dz, _, K, z0, shape = LEGS[leg]
    sig = sigma_m(shape)
    h0s = h0_of(shape)
    n = n_of_K[K]
    # ⚠⚠ ARGUMENT ORDER IS LOAD-BEARING. `max(CONST, var)`, never
    # `max(var, CONST)`. Measured on this binary: with a max/min call in a
    # DIVISOR, the variable-first order is mis-evaluated (the call becomes a
    # MULTIPLIER), while the constant-first order is exact. Calibration-free
    # proof: `0.09/max(0.09,s)` and `0.09/max(s,0.09)` are the same number
    # mathematically and returned 381.871 W vs 2.236 W one step apart. Rev 2's
    # first draft wrote `max(s,0.0632)` and delivered 0.26 W where the design
    # asked ~240 W -- a factor 250. The Meier h map is const-first in all 23 of
    # its max/min calls, so it is NOT affected.
    decay = f"pow({g(S0)}/max({g(S0)},s),{g(n)})"
    core = (f"{g(h0s)}*{decay}" if sig is None
            else f"{g(h0s)}*exp(-(r*r)/({g(2.0*sig*sig)}))*{decay}")
    e = f"if((r>{g(R_OPEN)})and(s<{g(s_lid(leg))}),{g(LID_H)},{core})"
    assert " " not in e and "&&" not in e, e
    for bad in ("max(s", "max(r", "max(t", "max(x", "max(y",
                "min(s", "min(r", "min(t", "min(x", "min(y"):
        assert bad not in e, f"{leg}: variable-first {bad}...) -- see the note " \
                             f"in h_expr(); use max(CONST,var)"
    return e


def mirror(leg, r, s, n_of_K):
    """My own evaluation of the same formula -- the branch witness compares
    BRANCHES, not sizes (the shape branch is far below the lid floor for the
    narrow shapes, which is physical and not a defect)."""
    dz, _, K, z0, shape = LEGS[leg]
    sig = sigma_m(shape)
    h0s = h0_of(shape)
    n = n_of_K[K]
    if r > R_OPEN and s < s_lid(leg):
        return LID_H, "lid"
    d = (S0 / max(S0, s)) ** n
    v = h0s * d if sig is None else h0s * math.exp(-(r * r) / (2 * sig * sig)) * d
    return v, "shape"


def main():
    P = print
    De = read_De()
    her = her_table()
    n_of_K = {}
    fails = []

    P("P2-5 rev 2 (a1) -- the lab h map, derived")
    P("=" * 88)
    P(f"  jet_D_e  {De:.6f} m   READ from probe/SMOKE0 (jet_T_ent_mode = fixed"
      f" => D_dec = jet_De_ref)")
    P(f"  her table at 38 kW:  L -> q [MW/m2], h [kW/m2K], T_sens = 1866 - q/h")
    for L in sorted(her):
        q, h, Ts = her[L]
        P(f"    L {L:2d}   s {L*D_NOZ*1e3:6.1f} mm   q {q/1e6:5.2f}"
          f"   h {h/1e3:5.2f}   T_sens {Ts:7.2f} K")
    P("")

    # ---- (i) + (ii) -------------------------------------------------------
    for K in (30, 20):
        n, rows = fit_n(K, De, her)
        n_of_K[K] = n
        core_m = K * De
        P(f"  (i)+(ii) K{K}:  jet_core_length = {K:g} x D_e = {core_m*1e3:.1f}"
          f" mm = {core_m/D_NOZ:.1f} D   (phi = 1 out to s = {core_m*1e3:.1f} mm)")
        P(f"      {'L':>3s} {'s[mm]':>7s} {'phi':>7s} {'T_stag':>9s}"
          f" {'want':>6s} {'T_sens':>8s} {'h_req':>9s} {'h_law':>9s}"
          f" {'q_law':>9s} {'q_her':>8s} {'resid':>8s} {'want':>7s}")
        for (L, s, phi, Tg, Ts, q, hreq) in rows:
            hlaw = H0 * (S0 / s) ** n
            qlaw = hlaw * (Tg - Ts)
            res = 100.0 * (qlaw - q) / q
            wt, wr = WANT_TSTAG[K][L], WANT_RES[K][L]
            P(f"      {L:3d} {s*1e3:7.1f} {phi:7.4f} {Tg:9.2f} {wt:6d}"
              f" {Ts:8.2f} {hreq:9.1f} {hlaw:9.1f} {qlaw/1e6:8.3f}M"
              f" {q/1e6:7.2f}M {res:+7.2f}% {wr:+6.1f}%")
            if abs(Tg - wt) > 1.0:
                fails.append(f"a1(i) K{K} L{L}: T_stag {Tg:.2f} vs {wt}")
            if abs(res - wr) > 0.3:
                fails.append(f"a1(ii) K{K} L{L}: residual {res:+.2f} vs {wr:+.1f}")
            if abs(res) > 15.0:
                fails.append(f"a1(ii) K{K} L{L}: |residual| {abs(res):.1f} > 15 %")
        P(f"      n = {n:.6f}   (planner {WANT_N[K]})   h0 = {H0:.0f} W/m2K"
          f"   worst |residual| {max(abs(100.0*(H0*(S0/s)**n*(Tg-Ts)-q)/q) for (L,s,phi,Tg,Ts,q,_) in rows):.2f} %")
        if abs(n - WANT_N[K]) > 0.001:
            fails.append(f"a1(ii) K{K}: n {n:.6f} vs {WANT_N[K]}")
        P("")

    # ---- (iv) the shapes --------------------------------------------------
    P("  (iv) THE SHAPE LADDER   (sigma in units of D = 7.9 mm; R_disc = 4.5 mm)")
    P(f"      {'shape':>6s} {'sigma[mm]':>10s} {'f_disc':>8s} {'want':>6s}"
      f" {'h0_shape':>9s} {'q0@850K':>9s} {'P(Ø40)':>9s} {'want':>8s}"
      f" {'/flake':>14s} {'r_0.5[mm]':>10s} {'want':>6s}")
    for shape in ("A", "G10", "G075", "G05", "G035"):
        sig = SHAPES[shape]
        sigm = sigma_m(shape)
        f = f_disc(sigm)
        h0s = h0_of(shape)
        q0 = h0s * (T_NOZ - T_FACE_REF)
        P_ = P_open(sigm, h0s)
        rh = r_half(sigm, h0s)
        wp, wr = WANT_P[shape], WANT_R05[shape]
        P(f"      {shape:>6s} {('plateau' if sigm is None else f'{sigm*1e3:.3f}'):>10s}"
          f" {f:8.4f} {'1.000' if sig is None else f'{f:.3f}':>6s}"
          f" {h0s:9.1f} {q0/1e6:8.2f}M {P_/1e3:8.2f}k {wp/1e3:7.2f}k"
          f" {P_/FLAKE_HI:6.1f}-{P_/FLAKE_LO:.1f}x"
          f" {('n/a' if rh is None else f'{rh*1e3:.2f}'):>10s}"
          f" {('n/a' if wr is None else f'{wr:.1f}'):>6s}")
        if abs(P_ - wp) / wp > 0.02:
            fails.append(f"a1(iv) {shape}: P {P_:.0f} W vs {wp:.0f} (>2 %)")
        if wr is not None and abs(rh * 1e3 - wr) > 0.1:
            fails.append(f"a1(iv) {shape}: r_0.5 {rh*1e3:.2f} vs {wr}")
    P(f"      flake bound {FLAKE_LO:.0f}-{FLAKE_HI:.0f} W (full circle, so both"
      f" sides of the ratio are full-scale -- SR 3)")
    P("")

    # ---- (iii) + (v) the per-leg strings and the branch witness ----------
    P("  (iii)+(v) THE PER-LEG MAP, with the BRANCH WITNESS")
    P(f"      {'leg':>9s} {'dz':>5s} {'s_leg':>8s} {'s_lid':>8s} {'K':>4s}"
      f" {'shape':>6s}   {'below s_lid':>22s}   {'above s_lid':>22s}")
    for leg in ORDER:
        dz, st, K, z0, shape = LEGS[leg]
        sl = s_lid(leg)
        vb, bb = mirror(leg, 0.021, sl - dz / 2.0, n_of_K)
        va, ba = mirror(leg, 0.021, sl + dz, n_of_K)
        ok = (bb == "lid" and ba == "shape")
        P(f"      {leg:>9s} {dz*1e3:4.0f}mm {s_leg(leg):8.5f} {sl:8.5f}"
          f" {K:4.0f} {shape:>6s}   {bb:>5s} h={vb:<12.6g}   {ba:>5s}"
          f" h={va:<12.6g} {'ok' if ok else 'FAIL'}")
        if not ok:
            fails.append(f"a1(v) {leg}: branch witness {bb}/{ba}")
    P("      (r = 0.021 m, just outside the Ø 40 opening; the test is WHICH")
    P("       branch, not its size -- the shape branch is far below the 0.001")
    P("       lid floor out here for the narrow shapes, which is the physics.)")
    P("")
    P("  (v2) STRING-vs-TABLE CROSS-CHECK -- the emitted constants must BE the")
    P("       table's. Rev 2's first draft passed sigma in units of D into the")
    P("       string builder and emitted h0 = 7,100 for every shape while the")
    P("       table read 9,651; this row is what made the two disagree aloud.")
    for leg in ORDER:
        shape = LEGS[leg][4]
        e = h_expr(leg, n_of_K)
        h0s = h0_of(shape)
        if g(h0s) not in e:
            fails.append(f"a1(v) {leg}: h0_shape {h0s:.1f} absent from the string")
        sg = sigma_m(shape)
        if sg is not None and g(2.0 * sg * sg) not in e:
            fails.append(f"a1(v) {leg}: 2*sigma^2 {2*sg*sg:.3e} absent")
    P(f"       {'ok' if not fails else 'FAIL'}: every leg's string carries its"
      f" own h0_shape and 2*sigma^2")
    P("")
    for leg in ORDER:
        P(f"H_EXPR_{leg}={h_expr(leg, n_of_K)}")
    P("")
    P("=" * 88)
    if fails:
        P(f"FAIL (a1): {len(fails)} assertion(s)")
        for f in fails:
            P(f"   {f}")
    else:
        P("PASS (a1 offline half): T_stag tables, h0, both exponents, every flux")
        P("  residual <= 15 %, per-shape P within 2 %, r_0.5, the lid constants")
        P("  and the branch witness all reproduce the pre-computed table.")
        P("  ⚠ a1 also requires `run_p2_5.py --probe` on every leg (the")
        P("    parse-back); that half is scored in PROBE.txt.")


if __name__ == "__main__":
    main()
