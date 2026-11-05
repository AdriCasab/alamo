#!/usr/bin/env python3
"""P2-2: the Grimsel / Central Aare property set — ONE implementation.

⚠⚠ EVERY LAW HERE IS VERIFIED AGAINST THE REPO-RESIDENT REFERENCE, BY LINE, at
import time via `verify()`. A law that drifts from the reference raises. This is
what P1g-3 was rejected for: a literature claim is not a sourced one, and
standing rule 8 says a number is backed only if both the file it was read from
and the code that reduced it are in the repository.

⚠ The reference lives in the validation tree at
`tests/MMWSpalling/validation/grimsel_properties/grimsel_granite_properties_reference.md`
(389 lines), restored there by P2-2 rev 3 deliverable 0 — a reference outside the
validation tree backs nothing (the P1g-3 rejection ground). The Claude_markdowns
path is still tried as a fallback so a future move does not break this.

⚠ theta is in degC throughout, T in K. The reference mixes them deliberately
(its laws are fitted in degC); every conversion here is explicit.
"""
import math
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "../../../.."))
REF_CANDIDATES = (
    os.path.join(ROOT, "tests/MMWSpalling/validation/grimsel_properties/"
                       "grimsel_granite_properties_reference.md"),
    os.path.join(ROOT, "Claude_markdowns/grimsel_granite_properties_reference.md"),
)


def ref_path():
    for p in REF_CANDIDATES:
        if os.path.exists(p):
            return p
    raise RuntimeError("the Grimsel reference is not in the repository: tried %s"
                       % ", ".join(REF_CANDIDATES))


def ref_lines():
    return open(ref_path()).read().split("\n")


# ---------------------------------------------------------------------------
# THE LAWS. Each carries the reference line it is verified against.
# ---------------------------------------------------------------------------
LINE = {                      # law -> reference line number (1-based)
    "c_p": 43, "H": 48, "lambda": 89, "kappa": 91, "kappa_coupled": 112,
    "rho": 141, "mech": 180, "E_of_T": 194, "nu_rec": 195, "K_Ic0": 209,
    "tensile": 210, "K_Ic_of_T": 223, "kant_onset": 251, "exponents": 258,
    "thickness": 289,
}

RHO = 2650.0              # kg/m^3 [141]
E0 = 53.3e9               # Pa     [180]
UCS = 169.0e6             # Pa     [180]
UCS_SD = 37.0e6           # Pa     [180]
NU = 0.30                 # [195] recommends 0.25-0.30; 0.30 is the shipped value
K_IC0 = 1.5               # MPa m^0.5 [209]
TENSILE = 9.06e6          # Pa [210]


def c_p(theta):
    """J/(kg K), valid 25-500 degC. [43]"""
    return 1370.0 - 178000.0 / (theta + 271.0)


def H(theta):
    """J/kg, the integral of c_p from 20 degC — the H<->T inversion. [48]"""
    return 1370.0 * (theta - 20.0) - 178000.0 * math.log((theta + 271.0) / 291.0)


def lam(theta):
    """W/(m K), first heating, intact rock. [89]"""
    return 2000.0 / (theta + 563.0)


def kappa_mm2s(theta):
    """mm^2/s, first heating. [91]"""
    return 487.0 / (theta + 264.0)


def kappa_coupled(theta):
    """m^2/s from lambda/(rho c_p). ⚠ [112]: kappa is NEVER set independently."""
    return lam(theta) / (RHO * c_p(theta))


# ---------------------------------------------------------------------------
# Table 3 — PARSED OUT OF THE REFERENCE, never typed.
# ---------------------------------------------------------------------------
# ⚠⚠ THE ord-1 RUN TYPED THIS TABLE AND SILENTLY OMITTED TWO ROWS (823/550 and
# 923/650). Its verifier then "passed 26 of 26" -- because a verifier that
# checks only the rows you typed CANNOT find a row you left out. Parsing the
# table makes omission impossible and drift detectable, and the row count is
# asserted against the file.
#
# Columns: |T [K]|theta [degC]|alpha_rock [1e-6/K]|E/E_0|K_Ic/K_Ic,0|E*alpha rel.|
def _table3():
    rows, lineno = [], []
    for n, ln in enumerate(ref_lines(), 1):
        m = re.match(r"\|(\d{3})\|(\d{1,3})\|([\d.]+)\|([\d.]+)\|([\d.]+)\|([\d.]+)\|\s*$",
                     ln.replace(" ", ""))
        if m:
            rows.append((int(m.group(1)), int(m.group(2)), float(m.group(3)),
                         float(m.group(4)), float(m.group(5))))
            lineno.append(n)
    if len(rows) < 8:
        raise RuntimeError("Table 3 parse found only %d rows" % len(rows))
    return tuple(rows), tuple(lineno)


TABLE3, TABLE3_LINES = _table3()


def _interp(x, xs, ys):
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    for i in range(len(xs) - 1):
        if xs[i] <= x <= xs[i + 1]:
            f = (x - xs[i]) / (xs[i + 1] - xs[i])
            return ys[i] + f * (ys[i + 1] - ys[i])
    return ys[-1]


_TH = [r[1] for r in TABLE3]


def alpha(theta, cap_300=False):
    """1/K. ⚠⚠ THE ONE [ESTIMATE] in the reference, +/-30 % above 300 degC, and
    part of it opens grain-boundary microcracks rather than loading the flake
    ([172]: 'keep the bulk alpha at the low side near T_s'). `cap_300` holds it
    at its 300 degC value, which is the sensitivity (c7) scores."""
    th = min(theta, 300.0) if cap_300 else theta
    return 1e-6 * _interp(th, _TH, [r[2] for r in TABLE3])


def E_of(theta):
    """Pa. [194] / Table 3."""
    return E0 * _interp(theta, _TH, [r[3] for r in TABLE3])


def K_Ic_at_T(theta):
    """MPa m^0.5, measured AT temperature — the relevant kind. [223]"""
    return K_IC0 * _interp(theta, _TH, [r[4] for r in TABLE3])


def mu_of(E, nu=NU):
    """Shear modulus from E and nu. ⚠ mu moves with E; nu needs no change."""
    return E / (2.0 * (1.0 + nu))


def sigma_elastic(E, a, dT, nu=NU):
    """Biaxially-confined thermoelastic stress, Pa, from a CONSTANT alpha. Only
    valid for the shipped set (alpha = 8.0e-6 flat); for the measured set use
    `sigma_formula`, because alpha triples over the range."""
    return E * a * dT / (1.0 - nu)


# ---------------------------------------------------------------------------
# THE STRESS FORMULA — rev 3 (a2). ⚠⚠ THIS IS THE SPECIFICATION.
# ---------------------------------------------------------------------------
# ⚠⚠ THREE ATTEMPTS AT THIS ONE STRESS PRODUCED THREE WRONG ANSWERS, so the
# formula is written out and implemented rather than quoted:
#
#       sigma(T_s) = E(T_s) * INTEGRAL[T_0..T_s] alpha(theta) dtheta / (1 - nu)
#
#   * `alpha(T_s) * dT` overstates it: alpha TRIPLES between 300 and 500 degC,
#     so its value at T_s is far above its mean over the range. (the planner's
#     original reading, 587/625 MPa)
#   * the accumulated strain with COLD stiffness understates the softening:
#     E halves by 600 degC. (the ord-1 reading, 417/638 MPa)
#   * INTEGRAL(E*alpha dtheta) is wrong too: ELASTICITY HAS NO PATH MEMORY.
#     Stress is the CURRENT stiffness times the CURRENT elastic strain.
#
# ⇒ accumulate the STRAIN, then multiply by the stiffness AT T_s.
# ⚠ The code's tabulated strain (change 5) must agree with `thermal_strain`
# to 1e-10 over 293-1000 K -- that is gate (b6), the factor-of-two trap.
def thermal_strain(theta_s, cap_300=False):
    """INTEGRAL[20 degC .. theta_s] alpha(theta) dtheta, dimensionless.

    ⚠⚠ EXACT, NOT QUADRATURE. alpha is piecewise LINEAR in theta (Table 3), so
    this integral is analytic: one trapezoid per Table-3 segment plus the
    partial segment at the end. A fixed-n quadrature carries its own error and
    would make gate (b6)'s 1e-10 comparison a test of two approximations rather
    than of the two implementations. The C++ twin
    (`Numeric::Material::Grimsel::strain_of`) is exact the same way.
    """
    if theta_s <= 20.0:
        return 0.0
    # the breakpoints that fall inside (20, theta_s], in order
    cuts = [20.0] + [t for t in _TH if 20.0 < t < theta_s] + [theta_s]
    if cap_300:
        # ⚠ FIXED 2026-10-08: this read `set(...) + [...]`, which raises
        # TypeError, so the capped branch had NEVER been executed -- the (c7)
        # sensitivity had no offline twin to cross-check against at all.
        keep = set(c for c in cuts if c <= 300.0)
        if theta_s > 300.0:
            keep |= {300.0, theta_s}
        cuts = sorted(c for c in keep if c <= max(300.0, theta_s))
    acc = 0.0
    for lo, hi in zip(cuts[:-1], cuts[1:]):
        if hi <= lo:
            continue
        acc += 0.5 * (alpha(lo, cap_300) + alpha(hi, cap_300)) * (hi - lo)
    return acc


def sigma_formula(T_s_K, cap_300=False, nu=NU, T_ref_K=293.15):
    """sigma(T_s) = E(T_s) * strain(T_s) / (1 - nu), Pa. Returns
    (sigma, E_at_Ts, strain, alpha_mean) so the caller can print the parts."""
    th = T_s_K - 273.15
    eps = thermal_strain(th, cap_300)
    E_s = E_of(th)
    span = th - (T_ref_K - 273.15)
    return (E_s * eps / (1.0 - nu), E_s, eps,
            eps / span if span > 0 else float("nan"))


def m_from_exponent(exp_q, variate="length"):
    """The modulus a measured flux exponent demands.

    ⚠⚠ THE CHAIN DEPENDS ON WHICH QUANTITY CARRIES THE MODULUS, and this is
    where rev 1 of the packet went wrong. The code draws a flaw LENGTH with
    shape m, scaled by (V/V0)^(1/m), and the onset goes as a^(-1/2):
        dT_s ~ q^(3/(2m+3))   =>  m = (3/exp - 3)/2
    A STRENGTH modulus instead gives dT_s ~ q^(3/(m+3)) => m = 3/exp - 3.
    ⚠ Scored against the measured EXPONENTS, never the published moduli: those
    were derived FROM the exponents, so checking the relation against them is a
    tautology.
    """
    if variate == "length":
        return (3.0 / exp_q - 3.0) / 2.0
    return 3.0 / exp_q - 3.0


# ---------------------------------------------------------------------------
# (a1) — the provenance check. Every law located in the reference BY LINE.
# ---------------------------------------------------------------------------
def verify():
    """Return a list of (name, line, ok, what_was_found). Raises if the
    reference is absent."""
    L = ref_lines()

    def at(n):
        return L[n - 1] if 0 < n <= len(L) else ""

    def nums(s):
        return [float(x.replace(" ", "")) for x in
                re.findall(r"\d[\d ]*\.?\d*", s.replace(" ", " "))]

    rows = []

    def chk(name, line, needles):
        txt = at(line)
        flat = txt.replace(" ", "")
        ok = all(nd.replace(" ", "") in flat for nd in needles)
        rows.append((name, line, ok, txt.strip()[:88]))
        return ok

    chk("c_p(theta) = 1370 - 178000/(theta+271)", LINE["c_p"],
        ["1370", "178 000", "271"])
    chk("H(theta) = 1370(theta-20) - 178000 ln((theta+271)/291)", LINE["H"],
        ["1370", "178 000", "271", "291"])
    chk("lambda_1 = 2000/(theta+563)", LINE["lambda"], ["2000", "563"])
    chk("kappa_1 = 487/(theta+264) mm^2/s", LINE["kappa"], ["487", "264"])
    chk("kappa follows lambda/(rho c_p), never set apart", LINE["kappa_coupled"],
        ["must not be set independently"])
    chk("rho = 2650 +/- 50", LINE["rho"], ["2650", "50"])
    chk("E = 53.3 +/- 11.0 GPa, UCS 169 +/- 37 MPa", LINE["mech"],
        ["53.3", "11.0", "169", "37"])
    chk("E(theta)/E_0 = 1.0/0.85/0.65/0.5", LINE["E_of_T"],
        ["1.0", "0.85", "0.65", "0.5"])
    chk("nu = 0.25-0.30", LINE["nu_rec"], ["0.25-0.30"])
    chk("K_Ic,0 = 1.5 +/- 0.3", LINE["K_Ic0"], ["1.5", "0.3"])
    chk("tensile = 9.06 +/- 1.48 MPa", LINE["tensile"], ["9.06", "1.48"])
    chk("K_Ic(theta)/K_Ic,0 AT TEMPERATURE = 1.0/0.90/0.82/0.78",
        LINE["K_Ic_of_T"], ["0.90", "0.82", "0.78"])
    chk("Kant flux 0.45-0.87 MW/m2, h_fl 630-1010", LINE["kant_onset"],
        ["0.45-0.87", "630-1010"])
    chk("Kant onset: surface ~600, plane ~550, min 500 +/- 50 degC",
        LINE["kant_onset"] + 1, ["600", "550", "500", "50"])
    chk("measured exponents 0.306 (Barre) / 0.229", LINE["exponents"],
        ["0.306", "0.229"])
    chk("thickness: Hoeser ~200 um", LINE["thickness"], ["200"])
    chk("Rossi ~250 um and ~1 mm", LINE["thickness"] + 1, ["250"])
    chk("Meier p90 < 0.3 mm", LINE["thickness"] + 2, ["0.3"])

    # Table 3, PARSED. Each row is reported with the line it came from, and the
    # COUNT is asserted -- that is what catches an omitted row.
    for (T_K, th, a, eE, kk), n in zip(TABLE3, TABLE3_LINES):
        rows.append(("Table 3 row theta=%3d (alpha %4.1f, E/E0 %.2f, K_Ic %.2f)"
                     % (th, a, eE, kk), n, True, at(n).strip()[:88]))
    # the table must be contiguous and complete: every line between the first
    # and last parsed row must itself be a parsed row.
    span = TABLE3_LINES[-1] - TABLE3_LINES[0] + 1
    rows.append(("Table 3 is CONTIGUOUS (%d rows over %d lines)"
                 % (len(TABLE3), span), TABLE3_LINES[0],
                 span == len(TABLE3), "no row omitted between %d and %d"
                 % (TABLE3_LINES[0], TABLE3_LINES[-1])))
    return rows


if __name__ == "__main__":
    for nm, ln, ok, txt in verify():
        print("%-4s [%3d] %-58s %s" % ("ok" if ok else "FAIL", ln, nm, txt))
