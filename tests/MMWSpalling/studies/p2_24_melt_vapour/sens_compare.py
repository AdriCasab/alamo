#!/usr/bin/env python3
"""P2-24 ADDENDUM (2026-11-07) -- the melt heat-transport sensitivity: M0 (the record) beside LAM (lambda held at its
1000 degC value above it) and EPS (emissivity 0.5), the same observables as validation/oglesby/test --onset (m2-m6),
labelled correctly, plus the identity witness M0R (M0 re-run on the addendum binary: thermo.dat byte for byte) and the
numbers the review asked for: the vapour as MASS, the hottest cell, the liquid depth, the closure residual.
Reads only thermo.dat / metadata / run.log / the event CSV under legs/<leg>/run. Usage: sens_compare.py > RESULTS_SENS.txt"""
import hashlib
import importlib.machinery
import importlib.util
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
TEST = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "oglesby", "test")
loader = importlib.machinery.SourceFileLoader("oglesby_test", TEST)
spec = importlib.util.spec_from_loader("oglesby_test", loader)
ot = importlib.util.module_from_spec(spec)
loader.exec_module(ot)

RHO = 2650.0                      # the Grimsel law's rho (rock_grimsel.inputs grimsel.rho)
DX = 0.1 / 64                     # the deck's cell (1.5625 mm)
V_CELL = DX ** 3
LEGS = [("M0", "the record (P2-24 binary 98a261d3)"),
        ("LAM", "lambda held at 2000/1563 = 1.280 W/mK above 1000 degC (grimsel.lambda.cap_theta = 1000)"),
        ("EPS", "losses.epsilon = epsilon_high = 0.5 (the rock file's 0.8)")]


TIMERS = {"fr_t_removal", "fr_t_step"}   # wall-clock seconds per step: the only columns a re-run cannot reproduce


def thermo_identical(pa, pb):
    """byte-for-byte on every non-timer column; returns (identical, the set of moved columns)"""
    def load(p):
        L = open(p).read().splitlines()
        return [c.strip() for c in L[0].split("\t")], [r.split("\t") for r in L[1:] if r.strip()]
    ha, A = load(pa)
    hb, B = load(pb)
    if ha != hb or len(A) != len(B):
        return False, {"<header or row count>"}
    moved = set()
    for ra, rb in zip(A, B):
        for j, (x, y) in enumerate(zip(ra, rb)):
            if x != y and ha[j] not in TIMERS:
                moved.add(ha[j])
    return not moved, moved


def run_dir(n):
    return os.path.join(HERE, "legs", n, "run")


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()[:16]


def leg_binary(n):
    return open(os.path.join(HERE, "legs", n, "BINARY.txt")).read().split()[0][:16]


def main():
    P, V, TT = ot.doe_tables()
    print("P2-24 ADDENDUM: the melt heat-transport sensitivity (units: s, mm, C, W, g; the DOE clock = model + 210 s)")
    print("=" * 118)
    print("0. THE BINARIES AND THE IDENTITY WITNESS")
    for n in ("M0", "M0R", "LAM", "EPS"):
        print("   %-4s binary %s  (%s)" % (n, leg_binary(n), open(os.path.join(HERE, "legs", n, "argv.sh")).read()
                                            .split("oglesby/input")[1].split("plot_file")[0].strip() or "no extra key"))
    same, moved = thermo_identical(os.path.join(run_dir("M0"), "thermo.dat"), os.path.join(run_dir("M0R"), "thermo.dat"))
    print("   M0R vs M0 thermo.dat, every column but the two wall-clock timers (fr_t_removal, fr_t_step): %s%s" % (
        "IDENTICAL" if same else "DIFFERENT", "" if same else " -- moved: " + " ".join(sorted(moved))))
    print("   => the addendum binary with the key at its default (-1) reproduces P2-24's M0: %s" % ("YES" if same else "NO"))
    legs = [(n, ot.MLeg(run_dir(n))) for n, _ in LEGS]
    print("\n1. THE LEGS")
    for (n, desc), (_, L) in zip(LEGS, legs):
        print("   %-4s %s; onset %.1f s, %d step rows, beam lit to %.0f s, dt %.3g s" % (n, desc, L.t_on, len(L.rows),
                                                                                        L.t_off, L.dt))
    print("\n2. THE OBSERVABLES (every post-melt row aligned on the leg's own onset)")
    doe_post = [t for t in sorted(TT) if t >= ot.DOE_ON + 100 and t in P and math.isfinite(P[t])]
    d19 = ot.binned(doe_post, lambda t: P[t], 1.9, 2.1, lambda t: TT[t])
    d33 = ot.binned(doe_post, lambda t: P[t], 3.2, 3.4, lambda t: TT[t])
    pre_doe = [TT[t] for t in sorted(TT) if ot.SHIFT <= t <= ot.DOE_ON]
    hdr = "   %-44s %10s %10s %10s %12s" % ("observable", "M0", "LAM", "EPS", "DOE #19")
    print(hdr)
    print("   " + "-" * (len(hdr) - 3))
    rows = []

    def obs(name, f, doe):
        vals = []
        for n, L in legs:
            try:
                vals.append(f(L))
            except Exception:
                vals.append(math.nan)
        rows.append((name, vals, doe))

    def fmt(v):
        if isinstance(v, str):
            return v
        if v is None or (isinstance(v, float) and math.isnan(v)):
            return "-"
        if abs(v) >= 100:
            return "%.0f" % v
        if abs(v) >= 1:
            return "%.2f" % v
        return "%.4f" % v

    def at70(L, k):
        r = L.at(L.t_on + 70.0)
        return r[k] / r["beam_P_beam"]

    def plateau(L, lo, hi):
        post = [r for r in L.lit if r["time"] >= L.t_on + 100.0]
        return ot.binned(post, lambda r: r["beam_P_beam"], lo, hi, lambda r: ot.C(r["bc_T_rad"]))[0]

    def rad_at(L, d):
        r = L.at(L.t_on + d)
        return ot.C(r["bc_T_rad"])

    def cool(L, d):
        r = L.at(L.t_off + d)
        return ot.C(r["bc_T_rad"])

    def spot(L):
        return 2e3 * math.sqrt(L.at(L.t_off)["bc_A_melt"] / math.pi)

    def before(L):
        b = [r for r in L.rows if r["time"] < L.t_on]
        return ot.C(b[-1]["bc_T_rad"])

    def first_vap(L):
        fv = next((r for r in L.rows if r.get("fr_n_vapor", 0) > 0), None)
        return fv["time"] if fv else "none"

    def first_buried(L):
        fv = next((r for r in L.rows if r.get("fr_n_vapor_int", 0) > 0), None)
        return fv["time"] if fv else "none"

    def vap_mass(L):
        return L.rows[-1].get("bc_n_vapor", 0) * V_CELL * RHO * 1e3

    def removals_post(L):
        return sum(ot.removals(r) for r in L.rows if r["time"] >= L.t_on)

    def tpeak(L):
        return max(r["bc_T_peak"] for r in L.rows if math.isfinite(r["bc_T_peak"]))

    def zliq_end(L):
        return 1e3 * L.at(L.t_off)["bc_z_liq"]

    def zliq_max(L):
        return 1e3 * max(r["bc_z_liq"] for r in L.rows if math.isfinite(r["bc_z_liq"]))

    def closure_post(L):
        return max(r["beam_closure_err"] for r in L.lit if r["time"] >= L.t_on)

    def collapse(L):
        tc, _ = L.collapse(0.148)
        return tc if math.isfinite(tc) else "never"

    def plat_at_off(L):
        return ot.C(L.at(L.t_off)["bc_T_rad"])

    obs("onset: first finite bc_z_melt [s]", lambda L: L.t_on, "285 (495 DOE)")
    obs("absorbed P_dep/P_beam at onset + 70 s", lambda L: at70(L, "beam_P_dep"), "0.70 (model)")
    obs("absorbed P_dep/P_beam, 1.9-2.1 kW rows", lambda L: ot.binned(
        [r for r in L.lit if r["time"] >= L.t_on], lambda r: r["beam_P_beam"], 1900, 2100,
        lambda r: r["beam_P_dep"] / r["beam_P_beam"])[0], "0.70 (model)")
    obs("transmission to x0.148 after onset [s]", collapse, "70 (pick-off)")
    obs("radiometer just before onset [C]", before, "153 mean / 272 max")
    obs("radiometer at onset + 20 s [C]", lambda L: rad_at(L, 20), "%.0f" % TT.get(515, math.nan))
    obs("radiometer at onset + 100 s [C]", lambda L: rad_at(L, 100), "%.0f" % TT.get(595, math.nan))
    obs("plateau at 1.9-2.1 kW (>= onset + 100 s) [C]", lambda L: plateau(L, 1900, 2100), "%.0f" % d19[0])
    obs("plateau at 3.2-3.4 kW [C]", lambda L: plateau(L, 3200, 3400), "%.0f" % d33[0])
    obs("plateau rise 2.0 -> 3.3 kW [%]", lambda L: 100 * (plateau(L, 3200, 3400) / plateau(L, 1900, 2100) - 1),
        "%+.1f" % (100 * (d33[0] / d19[0] - 1)))
    obs("radiometer at beam-off [C]", plat_at_off, "~2281 (last)")
    obs("cooling: beam-off + 15 s [C]", lambda L: cool(L, 15), "%.0f" % TT.get(1805, math.nan))
    obs("cooling: beam-off + 65 s [C]", lambda L: cool(L, 65), "%.0f" % TT.get(1855, math.nan))
    obs("cooling: beam-off + 215 s [C]", lambda L: cool(L, 215), "%.0f" % TT.get(2005, math.nan))
    obs("melt spot at the end of heating [mm]", spot, "~38")
    obs("first surface vapour void [s]", first_vap, "none seen (#28: 1.3 g / 4)")
    obs("first buried vapour count [s]", first_buried, "-")
    obs("vapour voided by 1800 s [g] (cells x V x rho)", vap_mass, "0 (#19); 1.3/4 exp. (#28)")
    obs("vapour cells voided (mesh-dependent)", lambda L: L.rows[-1].get("bc_n_vapor", 0), "-")
    obs("hottest cell over the run [K]", tpeak, "< 3233 if no vapour")
    obs("liquid depth on the axis at beam-off [mm]", zliq_end, "-")
    obs("liquid depth on the axis, max [mm]", zliq_max, "-")
    obs("post-melt removals (single-cell wall collapses)", removals_post, "0 flakes")
    obs("beam closure residual, post-onset max", closure_post, "<= 3e-2 (P2-24 b2r)")
    for name, vals, doe in rows:
        print("   %-44s %10s %10s %10s %12s" % ((name,) + tuple(fmt(v) for v in vals) + (doe,)))
    print("\n3. THE MELT SPOT AGAINST TIME / POWER (2 sqrt(A_melt / pi), mm)")
    for n, L in legs:
        grow = " ".join("%.0fs/%.2fkW/%.1fmm" % (r["time"], r["beam_P_beam"] / 1e3,
                                                 2e3 * math.sqrt(r["bc_A_melt"] / math.pi))
                        for r in L.lit if abs(r["time"] / 200.0 - round(r["time"] / 200.0)) < 1e-9)
        print("   %-4s %s" % (n, grow))
    print("\n4. THE VAPOUR AGAINST TIME (bc_n_vapor cells / bc_n_vapor_int cell-steps)")
    for n, L in legs:
        prog = " ".join("%.0fs:%d/%d" % (r["time"], r.get("bc_n_vapor", 0), r.get("bc_n_vapor_int", 0))
                        for r in L.rows if r["time"] in (300.0, 600.0, 900.0, 1200.0, 1580.0, 1800.0))
        print("   %-4s %s; ledger_E_vapor %.4g J" % (n, prog, L.rows[-1].get("ledger_E_vapor", math.nan)))
    print("\n5. THE RADIOMETER AGAINST TIME (bc_T_rad, C; the DOE's 137 GHz top panel where a reading exists)")
    print("   model t [s]  " + " ".join("%6d" % t for t in (50, 100, 150, 200, 300, 400, 600, 800, 1000, 1200, 1400, 1580)))
    for n, L in legs:
        print("   %-12s " % n + " ".join("%6.0f" % ot.C(L.at(float(t))["bc_T_rad"]) for t in
                                        (50, 100, 150, 200, 300, 400, 600, 800, 1000, 1200, 1400, 1580)))
    doe_row = []
    for t in (50, 100, 150, 200, 300, 400, 600, 800, 1000, 1200, 1400, 1580):
        td = t + ot.SHIFT
        near = [k for k in TT if abs(k - td) <= 5]
        doe_row.append("%6.0f" % TT[near[0]] if near else "     -")
    print("   %-12s " % "DOE" + " ".join(doe_row))
    print("=" * 118)
    print(READING)
    return 0


READING = """6. THE READING (advisor, 2026-11-07)
   identity   M0R reproduces M0 on every physics column (only the two wall-clock timers differ): the key at its default
              moves nothing. The addendum binary c052c348 is otherwise UNVERIFIED (no tree / face gate re-run); the P2-24
              binary 98a261d3 is kept in bin_p2_24/.
   LAM        doubling lambda above 1000 degC (1.28 vs the fit's 0.60 W/mK at 3000 K) leaves the surface where it was
              (plateau 1445 / 1630 vs 1449 / 1623 C; onset + 20 s 509 vs 513 C) and halves the vapour (3.28 -> 1.74 g,
              first void 92 -> 113 s) -- but the pool deepens (axis liquid 3.9 -> 7.0 mm at beam-off), widens (34.9 -> 36.7
              mm) and the cooling match M0 HAD is lost (+65 s 812 vs 699 C, DOE 627; +215 s 427 vs 327, DOE 349).
              Conduction is isotropic: more of it moves the buried heat down and sideways as much as up. A conductivity
              value, whatever it is, cannot warm the surface without breaking the cooling curve. Its 136 removals are
              all h_f = 0 wall collapses, every one AFTER beam-off (1585-1605 s, the lid solidifying).
   EPS        emissivity 0.5 warms the view mean by 2-6 % (1475 / 1721 C; the rise 16.6 % vs the DOE's 34.7 %), still
              300-670 C short, ADDS vapour (3.96 g, +21 %; the hottest cell 4183 K) and worsens the cooling (1058 / 798 /
              378 C vs 1011 / 627 / 349). As Z0 item 7 predicted: (0.8/0.5)^(1/4) = 1.125 on a radiation-limited T.
   common     the onset, the absorbed fraction (0.66-0.70) and the closure residual (4.5 %) do not move: none of them is a
              thermal-transport quantity, and the residual is the beam's spreading, not the melt's. Every non-vapour
              removal on all three legs has h_f = 0 (the no-thin-flake null stands through melting, every variant).
   verdict    BOTH recorded alternatives are EXCLUDED as the explanation. A cold surface over a vaporising interior is not
              a factor of two in lambda and not the emissivity; it needs a transport that carries heat PREFERENTIALLY to the
              free surface -- convection in the pool, or radiative transfer through the semi-transparent melt. That is a
              physics packet with its own evidence, not a knob. The alternative is to carry the two misses as the stated
              limit of a conduction-only melt and move on to the regime check.
"""


if __name__ == "__main__":
    sys.exit(main())
