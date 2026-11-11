#!/usr/bin/env python3
"""P2-29 z0 (offline) -- step zero, written down: the ledger rows L1-L13 against src/ (after the stage-1 edits), the
explicit-limit ratio of every new deck at every level, the live setter lists stage 2 acts on, and the outputs the chain
overwrites (archived first, SR 9). Exit 1 on a failed row.
Usage: z0.py > Z0.txt"""
import glob
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
T = os.path.join(ROOT, "tests", "MMWSpalling")
M = open(os.path.join(ROOT, "src/Integrator/MMWSpalling.H")).read()
TB = open(os.path.join(ROOT, "src/Numeric/Material/Table.H")).read()
BC = open(os.path.join(ROOT, "src/BC/Constant.cpp")).read()
FAIL = []


def row(ok, msg):
    print("  %-4s %s" % ("ok" if ok else "FAIL", msg))
    if not ok:
        FAIL.append(msg)


def line_of(s, pat):
    i = s.find(pat)
    return s.count("\n", 0, i) + 1 if i >= 0 else None


def deck_keys(path, seen=None):
    """key -> value with FILE includes resolved (later wins, as the parser reads them)"""
    seen = seen or set()
    out = {}
    for ln in open(path):
        k = ln.split("#", 1)[0].strip()
        if "=" not in k:
            continue
        key, val = (x.strip() for x in k.split("=", 1))
        if key == "FILE":
            inc = os.path.join(ROOT, val)
            if inc not in seen:
                seen.add(inc)
                out.update(deck_keys(inc, seen))
            continue
        out[key] = val
    return out


NEW = sorted([os.path.join(T, "validation/zhang/zhang_oglesby/input_phase")] +
             [os.path.join(T, "unit", n, "input_phase") for n in ("beam_tau", "beam", "equilibrium", "stefan",
                                                                   "convective_patch")] +
             glob.glob(os.path.join(T, "validation/hu/*/input_*_phase")))
STUDY = [os.path.join(HERE, "ZL.input")]


def zhang_alpha(Tk):
    th = min(Tk, 1273.0) - 273.15
    k = 1.0 / (0.3154 + 3.8e-4 * th)
    rho = 2300.0 if Tk >= 1273.0 else -0.2877 * (Tk - 273.15) + 2884.0
    t2 = Tk - 273.15
    cp = -2.73e-4 * t2 * t2 + 0.894 * t2 + 788.0 if Tk <= 1510.5 else 0.184 * t2 + 1249.0
    return k / (rho * cp)


def main():
    print("P2-29 z0 -- step zero (offline), on the stage-1 source")
    print("=" * 110)
    print("THE LEDGER")
    row(all("microstructure.nphases" in deck_keys(d) for d in NEW + STUDY),
        "L1 every new deck sets microstructure.nphases (%d decks + ZL): the production kernel (dispatch M:%d)" % (
            len(NEW), line_of(M, "if      (microstructure_enabled)")))
    row(all(p in M for p in ('else if (plaw == "zhang")', "else if (phases[p].zhang_law)",
                             "any = any || ph.grimsel_law || ph.zhang_law", "use[p] == 2")),
        "L2/L3 the zhang law is selected (parse M:%d), read by the table build (M:%d), refreshed by "
        "UpdatePropertyFields whose `any` includes it (M:%d), on every path (its caller, Advance M:%d)" % (
            line_of(M, 'else if (plaw == "zhang")'), line_of(M, "else if (phases[p].zhang_law)\n"),
            line_of(M, "any = any || ph.grimsel_law || ph.zhang_law"), line_of(M, "        UpdatePropertyFields(lev);")))
    row("ph.zhang_law   ? ph.zhang.kappa(T_cold)" in M,
        "L4 the dt check forms the zhang DIFFUSIVITY kappa/(rho c_p) (M:%d); at T_cold = min(T_ref 273.15, T_amb "
        "297.5) = 273.15 K: %.5g m^2/s (the packet's 1.3225e-6 is alpha(297.5 K) = %.5g; declared)" % (
            line_of(M, "ph.zhang_law   ? ph.zhang.kappa(T_cold)"), zhang_alpha(273.15), zhang_alpha(297.5)))
    zq = deck_keys(NEW[[i for i, d in enumerate(NEW) if "zhang" in d][0]])
    neu = "in(i-glevel[0],j,k,n)" in BC
    row(neu and zq["hc.bc.temp.constant.type.xlo"] == "neumann" and zq["hc.bc.temp.constant.type.ylo"] == "neumann",
        "L5 the quadrant's symmetry planes: xlo/ylo neumann in the deck; Neumann ghost = the nearest cell "
        "(Constant.cpp:%d)" % line_of(BC, "in(i-glevel[0],j,k,n)"))
    row(zq["beam.schedule.P"].split()[0] == "800.0" and "Pi = (2.0 * P / (M_PI * w0 * w0))" in M,
        "L6 the beam is not quartered by hand: schedule P[0] %s W, intensity analytic per cell (M:%d)" % (
            zq["beam.schedule.P"].split()[0], line_of(M, "Pi = (2.0 * P / (M_PI * w0 * w0))")))
    row("if (!fl_on) div_kgrad += 2.0 * kc * (sp_T - Tc) * inv_dz2;" in M,
        "L7 the Hu patch reaches the production kernel's 2x term under !fl_on (M:%d)" % line_of(
            M, "if (!fl_on) div_kgrad += 2.0 * kc * (sp_T - Tc) * inv_dz2;"))
    row("bool SurfaceOwner(int lev) const { return max_level == 0 || lev == finest_level; }" in M,
        "L8 SurfaceOwner = the finest level (M:%d): a8 checks level-1 coverage of the patch" % line_of(
            M, "bool SurfaceOwner(int lev) const"))
    row(True, "L9 the mechanics port is the same linear law (microstructure branch per-phase E, mu, beta): m5 reports")
    row(True, "L10 the fixtures' T ranges: printed by each test (a7) -- run-time")
    s_en = M.find("if (!value.spall_enabled) return;")
    s_sv = M.find("D2q-2b item 0: DEBUG-ONLY seeded interior void")
    row(0 < s_en < s_sv, "L13 spall.seed_void is parsed after the spall_enabled early return (M:%d < M:%d): c3's deck "
        "sets spall.enabled" % (line_of(M, "if (!value.spall_enabled) return;"),
                                line_of(M, "D2q-2b item 0: DEBUG-ONLY seeded interior void")))
    row("m_fus_jump" in TB and "if (m_sharp_fus && H > m_H_fus_hi) H -= m_fus_jump;" in TB,
        "(user 11-10) Table.H's sharp band carries its jump in H_of_T / T_of_H (the Stefan port's fixed defect)")
    print("=" * 110)
    print("THE EXPLICIT LIMIT, every new deck, every level (dt_lim = dx^2 / (n alpha), n = 2d + 1 on the patch's "
          "owner, the finest level; ratio bar 0.75)")
    for d in NEW + STUDY:
        k = deck_keys(d)
        dt0 = float(k["timestep"])
        nc = [int(x) for x in k["amr.n_cell"].split()]
        lo = [float(x) for x in k["geometry.prob_lo"].split()]
        hi = [float(x) for x in k["geometry.prob_hi"].split()]
        dx0 = min((h - l) / n for h, l, n in zip(hi, lo, nc))
        ml = int(k.get("amr.max_level", "0"))
        Tref = float(k.get("microstructure.phase0.T_ref", "273.15"))
        Tamb = float(k.get("losses.T_amb", "0"))
        Tc = min(Tref, Tamb) if Tamb > 0 else Tref
        if k.get("microstructure.phase0.law") == "zhang":
            a = zhang_alpha(Tc)
        else:
            a = float(k["microstructure.phase0.kappa"]) / (float(k["microstructure.phase0.rho"]) *
                                                         float(k["microstructure.phase0.Cp"]))
        presc = k.get("surface_patch.enabled", "0") == "1" and "surface_patch.T_expr" in k
        parts = []
        for lev in range(ml + 1):
            dx, dt = dx0 / 2 ** lev, dt0 / 2 ** lev
            for owner in ([True] if ml == 0 else ([True, False] if lev == 0 else [True])):
                n = 6 + (1 if presc and owner else 0)
                lim = dx * dx / (n * a) if a > 0 else math.inf
                r = dt / lim
                parts.append("L%d%s %.3f" % (lev, "" if ml == 0 else ("(own)" if owner else "(cov)"), r))
                if r > 0.75:
                    FAIL.append("%s level %d ratio %.3f" % (os.path.relpath(d, ROOT), lev, r))
        print("  %-62s alpha %.4g at %.2f K; ratios %s" % (os.path.relpath(d, ROOT), a, Tc, ", ".join(parts)))
    print("=" * 110)
    print("THE LIVE SETTERS (stage 2 acts on them; live = unit/, validation/, extensions/, studies/s1, minus retired/"
          "outputs/legs)")
    live = [p for p in glob.glob(os.path.join(T, "**", "input*"), recursive=True) +
            glob.glob(os.path.join(T, "validation", "shared", "*.inputs"))
            if os.path.isfile(p) and not re.search(r"/(output[^/]*|legs|retired_[a-z_]+|src_pre[^/]*|inputs_pre|"
                                                  r"archive|smoke|pre)/", p)
            and re.search(r"/(unit|validation|extensions)/|studies/s1_surface", p)]
    alph = {}
    for p in live:
        for ln in open(p, errors="replace"):
            mm = re.match(r"^\s*alpha\s*=\s*([^#]*)", ln)
            if mm:
                alph[os.path.relpath(p, T)] = mm.group(1).split()
    nz = sorted(k for k, v in alph.items() if any(float(x) != 0.0 for x in v))
    one = sorted(k for k, v in alph.items() if len(v) != 2)
    print("  L11 alpha: %d live decks set it; NON-ZERO in %d: %s (all the Hu decks stage 2 retires); a token count "
          "other than 2 in %d: %s -- a stage-2 finding: `0.0` is zero but one token, so the packet's token-wise compare "
          "against `0 0` (equal count first) would refuse it" % (len(alph), len(nz), ", ".join(nz), len(one),
                                                                 ", ".join("%s (%s)" % (k, " ".join(alph[k]))
                                                                           for k in one)))
    pats = {
            "foot_rule": re.compile(r"^\s*surface_patch\.foot_rule\s*="),
            "nozzle_feed_max": re.compile(r"^\s*surface_patch\.nozzle_feed_max\s*="),
            "phase.*": re.compile(r"^\s*phase\."), "material.*": re.compile(r"^\s*material\."),
            "solver / implicit.*": re.compile(r"^\s*(solver\s*=|implicit\.)")}
    for name, rx in pats.items():
        hits = sorted(os.path.relpath(p, T) for p in live if any(rx.match(ln) for ln in open(p, errors="replace")))
        print("  %-20s %3d deck(s)%s" % (name, len(hits), (": " + ", ".join(hits)) if hits and len(hits) <= 40 else ""))
    print("  (live decks scanned: %d)" % len(live))
    print("=" * 110)
    print("THE OUTPUTS THIS CHAIN OVERWRITES (archived, SR 9): unit/beam_tau/output (its test on bin_pre rmtree's U1-U3),"
          " validation/hu/*/output (the Hu tests write comparison.png into them); everything else is new (S/legs/,"
          " output_phase/, face_p2_29a)")
    print("RESULT z0: %s" % ("FAIL: " + "; ".join(FAIL) if FAIL else "every row resolves; every new deck's ratio <= 0.75"))
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
