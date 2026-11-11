#!/usr/bin/env python3
"""P2-27 z0 (offline gate) -- everything step zero must establish before the first edit.

  1. the packet's Sources re-resolved to TODAY's lines: label, symbol, file:line (cited -> found); a missing anchor
     is a FAIL (the code drifted from the inventory)
  2. the dt table PER LEVEL (L3): every live deck (P2-18's deck walk, FILE = includes resolved, last definition wins)
     and every gate leg as it last ran (its `metadata`, last occurrence wins: the tree gate's outputs and
     witness/face_p2_26/) -- dt(l) = timestep / prod(nsubsteps[1..l]), dx(l) = dx0 / prod(ref_ratio), alpha_max over
     the phases at T_cold = min(phase T_ref, losses.T_amb if > 0) (Grimsel: lambda/(rho c_p) with the leg's own
     coefficients; constant law: kappa/(rho Cp)); dt_lim(l) = dx_min(l)^2/(6 alpha_max). The user's decision
     2026-11-09 (after this script's first run found 18 decks over the limit, 7 of them live unit tests with a uniform
     start temperature): the check ABORTS on face-list decks and WARNS elsewhere, so the bar is: no FACE-LIST deck or
     leg over the limit (STOP); every other deck over it is listed as WARN.
  3. the schedule proof (L4): the oglesby deck's two arrays -> 159 pairs, monotone; the first 158 P equal
     FIG34_FULL.txt's schedule; the 12 live schedule decks monotone; fig34_schedule.txt (once written) equal as text
  4. L6: robin_face/output/form_* Header ncomp (23, no Sp_field) vs ledger (36)
  5. the deck-edit list with DEFAULT/LIVE (L13): A_D lines; every live setter of a key the packet retires; the grimsel
     scalars; argv in live test scripts
  6. the mover lists: stage 1 only-in-new (unit/jet_d2c/output under hash.sh's patterns), stage 2 only-in-ref
     (spall_event's 8), stage 3 (every Cell_D in p2_26.md5 + jet_d2c's; face_p2_24.md5's), stage 4 (the t = 0
     plotfile of every grimsel face leg: Cell_D + Cell_H)
  7. units (L12); bin_pre/ = the P2-26 binary; src_pre_1/ complete and equal to src/
Usage: z0.py > Z0.txt   (exit 1 on a failed bar)
"""
import glob
import hashlib
import importlib.util
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
T = os.path.join(ROOT, "tests", "MMWSpalling")
W = os.path.join(T, "studies", "tree_sweep_0922", "witness")
SRC = {"M": "src/Integrator/MMWSpalling.H", "R": "src/Integrator/MMWSpalling/Removal.H",
       "F": "src/Integrator/MMWSpalling/FaceRemoval.H", "S": "src/Integrator/MMWSpalling/Shared.H",
       "B": "src/Numeric/MMWBeam.H", "G": "src/Numeric/Material/Grimsel.H",
       "I": "src/Integrator/Integrator.cpp"}
_spec = importlib.util.spec_from_file_location("a2", os.path.join(T, "studies", "p2_18_prereq", "a2_apply.py"))
a2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a2)
LINE = a2.LINE
BAD = []


def need(ok, msg):
    print("  %s  %s" % ("ok  " if ok else "FAIL", msg))
    if not ok:
        BAD.append(msg)


def lines(k):
    """the step-zero source: src_pre_1/ once it exists (the anchors are resolved against the code BEFORE any edit),
    else src/"""
    pre = os.path.join(HERE, "src_pre_1", SRC[k])
    return open(pre if os.path.exists(pre) else os.path.join(ROOT, SRC[k])).read().split("\n")


# ------------------------------------------------------------------ 1. sources ----
# (label, file key, regex, cited line). Regexes name a SYMBOL or a distinctive phrase, never a line number.
ANCHORS = [
    # stage 1
    ("retired-at-default table", "M", r"static const RetiredAtDefault retired_at_default\[\]", 158),
    ("retired-at-default compare (string)", "M", r'std::string\(r\.value\) == "\*" \|\| got == r\.value', 187),
    ("weibull.V0 parse", "M", r'"weibull\.V0"', 9294),
    ("weibull_V0 member", "M", r"Set::Scalar\s+weibull_V0\b", 8865),
    ("vol_factor = pow(V_cell / weibull_V0, inv_m)", "M", r"pow\(V_cell / weibull_V0, inv_m\)", 7629),
    ("phase scalar kappa REQUIRED", "M", r'query_required\(pfx \+ "\.kappa"', 411),
    ("phase scalar A_D REQUIRED", "M", r'query_required\(pfx \+ "\.A_D"', 417),
    ("phase .law", "M", r'query_default\(pfx \+ "\.law"', 424),
    ("grimsel_law = true", "M", r"value\.phases\[i\]\.grimsel_law = true", 442),
    ("t = 0 scalar -> field (kap[p] = phases[p].kappa)", "M", r"kap\[p\] = phases\[p\]\.kappa;", 7387),
    ("t = 0 rho/Cp law branch", "M", r"rh\[p\]  = phases\[p\]\.grimsel\.rho\(phases\[p\]\.T_ref\)", 7412),
    ("t = 0 field writes (kapa(i,j,k) = d_kap[p])", "M", r"kapa \(i, j, k\) = d_kap\[p\];", 7453),
    ("UpdatePropertyFields", "M", r"void UpdatePropertyFields\(", 4549),
    ("hc.heat.alpha required", "M", r'query_required\("hc\.heat\.alpha"', 247),
    ("alpha queryarr (D's path)", "M", r'pp_queryarr\("alpha", value\.alpha\)', 342),
    ("matrix start (damage requires microstructure)", "M", r"value\.damage_enabled && !value\.microstructure_enabled", 693),
    ("matrix end: vapor.A_beam derivation", "M", r"if \(value\.vapor_enabled && !value\.vapor_A_beam_override\)", 1195),
    ("sub-parse: flake_h_max vs h_spall_max", "M", r"value\.sp_flake_h_max >= value\.spall_h_spall_max", 5871),
    ("sub-parse: removal_energy flake <-> surface_T", "M", r"value\.sp_removal_energy_flake && !value\.sp_surface_T_apply", 5905),
    ("sub-parse: surface_T/flake <-> local_thermoelastic", "M",
     r"\(value\.sp_surface_T_recon \|\| value\.sp_removal_energy_flake\) &&$", 5920),
    ("sub-parse: bit <-> spall|vapor", "M", r"if \(!value\.spall_enabled && !value\.vapor_enabled\)", 6105),
    ("sub-parse: sp_weibull refuses damage (dies in C1a)", "M", r"incompatible with \"\s*$|is incompatible with", 9148),
    ("runtime cubic-cell abort (BuildFaceList)", "M", r'"surface_patch\.face_list = 1 requires cubic cells; got "', 3707),
    ("SurfacePatch::Parse", "M", r"static void Parse\(SurfacePatch\s*&", 8019),
    ("dt: kernel H = H_old + dt*(...)", "M", r"H_old\w*\(i, j, k\) \+ dt \*", 2254),
    ("dt: Grimsel kappa(T) = lambda/(rho cp)", "G", r"return lambda\(T\) / \(rho\(T\) \* cp\(T\)\);", 136),
    ("dt: timestep (Integrator.cpp)", "I", r'pp_query_required\("timestep", timestep', 30),
    ("dt: dt[i] = dt[i-1]/nsubsteps[i]", "I", r"dt\[i\] = dt\[i - 1\] / \(amrex::Real\)nsubsteps\[i\];", 197),
    ("schedule queryarr t_end", "B", r'pp\.queryarr\("schedule\.t_end", v\.sched_t\)', 146),
    ("schedule equal-length abort", "B", r"schedule\.t_end and \.P must have the same length", 149),
    ("power_at", "B", r"Set::Scalar power_at\(Set::Scalar t\) const", 156),
    # stage 2 (C1a)
    ("ParseDamageSettings", "M", r"static void ParseDamageSettings\(", 5370),
    ("ParseDamageSettings call", "M", r"^\s*ParseDamageSettings\(value, pp\);", 289),
    ("RefreshKappaEff", "M", r"void RefreshKappaEff\(int lev\)", 5469),
    ("DamageGrowthRHS", "M", r"static Set::Scalar DamageGrowthRHS\(", 5487),
    ("IntegrateDamageRK4", "M", r"static Set::Scalar IntegrateDamageRK4\(", 5501),
    ("UpdateDamageAfterMechanics", "M", r"void UpdateDamageAfterMechanics\(", 5519),
    ("UpdateDamageAfterMechanics call", "M", r"^\s*UpdateDamageAfterMechanics\(lev", 1571),
    ("ParseCohesiveSettings", "M", r"static void ParseCohesiveSettings\(", 5630),
    ("ParseCohesiveSettings call", "M", r"^\s*ParseCohesiveSettings\(value, pp\);", 296),
    ("UpdateCohesiveAfterMechanics", "M", r"void UpdateCohesiveAfterMechanics\(", 5683),
    ("UpdateCohesiveAfterMechanics call", "M", r"^\s*UpdateCohesiveAfterMechanics\(lev", 1572),
    ("ResetCohesiveDiagnostics", "M", r"void ResetCohesiveDiagnostics\(", 5748),
    ("PrincipalStressRatio", "M", r"Set::Scalar PrincipalStressRatio\(", 6123),
    ("PrincipalStressRatio caller", "R", r"ratio = PrincipalStressRatio\(sigma\);", 1201),
    ("FiniteMatrix", "M", r"static bool FiniteMatrix\(", 7088),
    ("SpallSampleMode", "M", r"enum class SpallSampleMode", 8592),
    ("spallation.model default damage_law", "M", r'"spallation\.model", model_str, std::string\("damage_law"\)', 9110),
    ("enum SpallationModel", "M", r"enum class SpallationModel \{ DamageLaw, SpWeibull \};", 8714),
    ("InitializeMicrostructure reset_damage", "M", r"if \(reset_damage\)", 7525),
    ("k_eff block (kappa_eff/k_eff)", "M", r"k_eff_mf", 7549),
    ("UpdateModel", "M", r"void UpdateModel\(int a_step, Set::Scalar a_time\) override", 1300),
    ("reset_damage call (Initialize)", "M", r"/\*reset_damage=\*/true\);", 1239),
    ("reset_damage call (Regrid)", "M", r"/\*reset_damage=\*/!damage_enabled\);", 1292),
    ("D ghost fill in AdvanceMicrostructure", "M", r"Util::RealFillBoundary\(\*D_mf\[lev\],\s+geom\[lev\]\);", 1963),
    ("GB conductance h_gb0*(1 - D_face)", "M", r"h_eff\s+= h_gb0 \* \(1\.0 - D_face\)", 2142),
    ("registration D", "M", r'RegisterNewFab\(value\.D_mf,\s+nullptr, 1, 1, "D"', 570),
    ("registration kappa_eff", "M", r'RegisterNewFab\(value\.kappa_eff_mf, nullptr, 1, 0, "kappa_eff"', 571),
    ("registration k_eff", "M", r'RegisterNewFab\(value\.k_eff_mf,\s+nullptr, 1, 0, "k_eff"', 572),
    ("registration damage_F_DP", "M", r'"damage_F_DP"', 575),
    ("registration gb_delta", "M", r'"gb_delta"', 581),
    ("Sp_field registration guard == SpWeibull", "M", r"if \(value\.spallation_model == SpallationModel::SpWeibull\)$", 632),
    ("sp.a0 required under sp_weibull", "M", r'query_required\("spallation\.sp\.a0", value\.sp_a0\)', 9156),
    ("kappa factor site 2110", "M", r"kp\(i, j, k\) \* std::exp\(-a_D \* Dc\)", 2110),
    ("kappa factor site 2128", "M", r"kp\(ii, jj, kk\) \* std::exp\(-a_D \* Df\(ii, jj, kk\)\)", 2128),
    ("kappa factor site 2640", "M", r"jkc\[col\] = kp\(i, j, k\) \* std::exp\(-kappa_damage_alpha", 2640),
    ("kappa factor site 3201", "M", r"const Set::Scalar kc = kp\(i, j, k\) \* std::exp\(-kappa_damage_alpha", 3201),
    ("kappa factor site 3760", "M", r"\(a_D != 0\.0\) \? kp\(i, j, k\) \* std::exp\(-a_D", 3760),
    ("kappa factor site 4189", "M", r"const Set::Scalar kc = kp\(i, j, k\) \* std::exp\(-a_D \* Df", 4189),
    ("kappa factor site 5482 (RefreshKappaEff)", "M", r"ke\(i, j, k\) = kp\(i, j, k\) \* std::exp\(-a_D \* Dv\)", 5482),
    ("kappa factor site F:455", "F", r"std::exp\(-kappa_damage_alpha \* Df\(p\[0\], p\[1\], p\[2\]\)\)", 455),
    ("R damage-law h_spall (C_h sqrt(alpha t))", "R", r"h_spall = C_h \* std::sqrt\(alpha_eff \* t_spall_eff\);", 56),
    ("R PASS 1.5 (spall.sample window)", "R", r"PASS 1\.5: per-column max D", 251),
    ("R damage-law firing (Dv >= D_th)", "R", r"if \(Dv >= D_th\)", 1189),
    ("R bulk-averaged vapour (damage_law)", "R", r"damage_law: BULK-averaged RoP_vap path", 1264),
    ("R winner else (damage_law)", "R", r"if \(hv > hs\) \{ reg_choice = 2; h_c = hv; \}", 1441),
    ("D write R:208 (pad clearance)", "R", r"^\s*D\(i, j, k\)\s+= 0\.0;", 208),
    ("D dead bind R:715", "R", r"amrex::Array4<const Set::Scalar> const& D   = \(\*D_mf\[lev\]\)\.array\(mfi\);", 715),
    ("abort: spall_event message", "M", r"flake is the intended answer \(unit/spall_event", 864),
    ("DruckerPrager include", "M", r'#include "Numeric/DruckerPrager\.H"', 36),
    ("CohesiveZone include", "M", r'#include "Numeric/CohesiveZone\.H"', 37),
    # stage 5 (G)
    ("FeetDescent", "M", r"Set::Scalar FeetDescent\(int lev, Set::Scalar dt,", 2443),
    ("foot_clear_k armed in FeetDescent", "M", r"foot_clear_k\.assign\(|foot_clear_k\[", 2551),
    ("FeetDescent call in BuildFlameColumns", "M", r"FeetDescent\(lev, dt, k_top_col\)", 2304),
    ("refusal M:1053 path heat refuses feet (STAYS)", "M", r'ph_abort\("surface_patch\.nozzle_descent = feet"\)', 1053),
    ("refusal M:1114 regime 4 (LIFT)", "M", r'fr_abort\("regime 4 \(surface_patch\.foot_body_clearance\)", "P2-11"\)', 1115),
    ("column clearance block mech_on", "R", r"const bool mech_on = surface_patch\.enabled && surface_patch\.feet &&", 159),
    ("face-list dispatch RecedeAlongLines", "R", r"RecedeAlongLines\(lev, time, dt, k_top_global\);", 247),
    ("RecedeAlongLines", "F", r"void RecedeAlongLines\(", 340),
    ("FrInitFraction in RecedeAlongLines", "F", r"FrInitFraction\(lev\);", 356),
    ("FrAssignTs", "F", r"void FrAssignTs\(", 162),
    ("vapour rule booking ledger_E_vapor", "F", r"ledger_E_vapor \+=", 943),
    ("vapour CSV regime literal ,2,", "F", r'",2,"', 1107),
    ("ledger_err_solid credits foot_mech_E", "M", r"foot_mech_E", 4353),
    ("foot_mech registration gate", "M", r'RegisterIntegratedVariable\(&value\.foot_mech_E,', 970),
]


def section1():
    print("== 1. Sources re-resolved against the step-zero source (src_pre_1/ when present; cited -> found; drift "
          "is reported, a missing anchor FAILS)")
    cache = {k: lines(k) for k in SRC}
    for label, fk, rx, cited in ANCHORS:
        r = re.compile(rx)
        hits = [i + 1 for i, ln in enumerate(cache[fk]) if r.search(ln)]
        if not hits:
            need(False, "%-52s %s:%d -> NOT FOUND (/%s/)" % (label, fk, cited, rx))
            continue
        best = min(hits, key=lambda h: abs(h - cited))
        drift = best - cited
        extra = "" if len(hits) == 1 else " (%d hits: %s)" % (len(hits), ",".join(map(str, hits[:6])))
        need(True, "%-52s %s:%d -> %s:%d%s%s" % (label, fk, cited, fk, best,
                                                 "" if drift == 0 else " (drift %+d)" % drift, extra))
    for h in ("src/Numeric/DruckerPrager.H", "src/Numeric/CohesiveZone.H"):
        users = []
        for p in glob.glob(os.path.join(ROOT, "src", "**", "*.H"), recursive=True) + \
                glob.glob(os.path.join(ROOT, "src", "**", "*.cpp"), recursive=True) + \
                glob.glob(os.path.join(ROOT, "src", "*.cc"), recursive=False):
            if os.path.abspath(p) == os.path.join(ROOT, h):
                continue
            if os.path.basename(h) in open(p, errors="replace").read():
                users.append(os.path.relpath(p, ROOT))
        need(users == [SRC["M"]], "%s included by %s (expected MMWSpalling.H only)" % (h, users))


# ------------------------------------------------------------- deck resolution ----
def resolve_deck(path, argv=()):
    """effective key -> value: FILE = includes in order (relative to the repo root), last definition wins,
    then argv tokens key=value"""
    kv = {}

    def walk(p, depth=0):
        for ln in open(p, errors="replace").read().split("\n"):
            m = LINE.match(ln)
            if not m:
                continue
            k, v = m.group(1), m.group(2).strip()
            if k == "FILE":
                walk(os.path.join(ROOT, v), depth + 1)
                continue
            kv[k] = v
    walk(path)
    for t in argv:
        if "=" in t:
            k, v = t.split("=", 1)
            kv[k] = v
    return kv


def read_metadata(path):
    kv = {}
    for ln in open(path, errors="replace"):
        if " = " not in ln:
            continue
        k, v = ln.rstrip("\n").split(" = ", 1)
        kv[k.strip()] = v.strip()      # last occurrence wins (ParmParse LAST)
    return kv


def nums(s):
    return [float(x) for x in s.replace(",", " ").split()]


def grimsel_alpha(kv, i, T):
    g = "microstructure.phase%d.grimsel." % i
    cp_a = float(kv.get(g + "cp.a", 1370.0))
    cp_b = float(kv.get(g + "cp.b", 178000.0))
    cp_th0 = float(kv.get(g + "cp.th0", 271.0))
    lam_a = float(kv.get(g + "lambda.a", 2000.0))
    lam_th0 = float(kv.get(g + "lambda.th0", 563.0))
    cap = float(kv.get(g + "lambda.cap_theta", -1.0))
    rho = float(kv.get(g + "rho", 2650.0))
    th = T - 273.15
    cp = cp_a - cp_b / (th + cp_th0)
    thl = cap if (cap > 0.0 and th > cap) else th
    lam = lam_a / (thl + lam_th0)
    return lam / (rho * cp)


def dt_rows(kv):
    """[(level, dt, dx_min, dt_lim, ratio, phase, T_cold, law)] or None (no microstructure)"""
    if "microstructure.nphases" not in kv:
        return None
    n = int(float(kv["microstructure.nphases"]))
    t_amb = float(kv.get("losses.T_amb", 0.0))
    alpha_max, arg = -1.0, None
    for i in range(n):
        p = "microstructure.phase%d." % i
        T_ref = float(kv.get(p + "T_ref", 273.15))
        T_cold = min(T_ref, t_amb) if t_amb > 0.0 else T_ref
        if kv.get(p + "law", "constant") == "grimsel":
            a, law = grimsel_alpha(kv, i, T_cold), "grimsel"
        else:
            a, law = float(kv[p + "kappa"]) / (float(kv[p + "rho"]) * float(kv[p + "Cp"])), "constant"
        if a > alpha_max:
            alpha_max, arg = a, (i, T_cold, law)
    nc = nums(kv["amr.n_cell"])
    lo, hi = nums(kv["geometry.prob_lo"]), nums(kv["geometry.prob_hi"])
    dx0 = min((hi[d] - lo[d]) / nc[d] for d in range(3))
    maxl = int(float(kv.get("amr.max_level", 0)))
    rr = nums(kv.get("amr.ref_ratio", "2")) or [2.0]
    rr = (rr + [rr[-1]] * maxl)[:max(maxl, 1)]
    ns_raw = nums(kv.get("amr.nsubsteps", "")) if kv.get("amr.nsubsteps", "").strip() else []
    if not ns_raw:
        ns = [1] + [rr[l - 1] for l in range(1, maxl + 1)]
    elif len(ns_raw) == 1:
        ns = [1] + [ns_raw[0]] * maxl
    else:
        ns = [1] + ns_raw[:maxl]
    dt0 = float(kv["timestep"])
    out, dt, dx = [], dt0, dx0
    for l in range(maxl + 1):
        if l > 0:
            dt /= ns[l]
            dx /= rr[l - 1]
        lim = dx * dx / (6.0 * alpha_max)
        out.append((l, dt, dx, lim, dt / lim, arg[0], arg[1], arg[2]))
    return out


def live_decks():
    return [p for p in a2.deck_files() if not p.endswith(".inputs")]


TREE_OUT = ["unit/robin_face/output", "unit/robin_pinned/output", "unit/robin_jet/output", "unit/jet_enthalpy/output",
            "unit/beam_void_closure/output", "unit/scalar_flaw/output", "unit/spall_event/output",
            "unit/spall_event/output_no_event", "unit/robin_feet/output", "unit/jet_d2c/output",
            "validation/meier/sp_meier_pilot/output/dev2d", "studies/s1_surface_resolution/output_smoke"]


def section2():
    print("\n== 2. dt table PER LEVEL (L3): dt(l) <= dx(l)^2/(6 alpha_max(T_cold)); STOP on any ratio > 1")
    worst = (0.0, None)
    warn = []
    nmicro = 0
    print("  -- live decks (resolved; %d files)" % len(live_decks()))
    for p in live_decks():
        try:
            kv = resolve_deck(p)
            rows = dt_rows(kv)
        except (KeyError, ValueError, ZeroDivisionError) as e:
            print("     ?    %-70s unresolvable (%s)" % (os.path.relpath(p, T), e))
            continue
        if rows is None:
            continue
        nmicro += 1
        fl = a2.same(kv.get("surface_patch.face_list", "0"), "1")
        for (l, dt, dx, lim, r, ph, Tc, law) in rows:
            flag = ("STOP" if fl else "WARN") if r > 1.0 else ""
            print("     %-70s L%d dt %-9.4g dx %-9.4g lim %-9.4g ratio %.3f  phase%d %s T_cold %.2f %s%s"
                  % (os.path.relpath(p, T), l, dt, dx, lim, r, ph, law, Tc, "FL " if fl else "", flag))
            if fl and r > worst[0]:
                worst = (r, "%s L%d" % (os.path.relpath(p, T), l))
            if r > 1.0 and not fl:
                warn.append(os.path.relpath(p, T))
    print("  -- gate legs as they last ran (metadata, last occurrence wins)")
    metas = []
    for d in TREE_OUT:
        metas += glob.glob(os.path.join(T, d, "**", "metadata"), recursive=True)
    metas += glob.glob(os.path.join(W, "face_p2_26", "**", "metadata"), recursive=True)
    for mp in sorted(metas):
        kv = read_metadata(mp)
        try:
            rows = dt_rows(kv)
        except (KeyError, ValueError, ZeroDivisionError) as e:
            print("     ?    %-70s unresolvable (%s)" % (os.path.relpath(mp, T), e))
            continue
        if rows is None:
            continue
        nmicro += 1
        fl = a2.same(kv.get("surface_patch.face_list", "0"), "1")
        for (l, dt, dx, lim, r, ph, Tc, law) in rows:
            flag = ("STOP" if fl else "WARN") if r > 1.0 else ""
            print("     %-70s L%d dt %-9.4g dx %-9.4g lim %-9.4g ratio %.3f  phase%d %s T_cold %.2f %s%s"
                  % (os.path.relpath(os.path.dirname(mp), T), l, dt, dx, lim, r, ph, law, Tc, "FL " if fl else "", flag))
            if fl and r > worst[0]:
                worst = (r, "%s L%d" % (os.path.relpath(os.path.dirname(mp), T), l))
            if r > 1.0 and not fl:
                warn.append(os.path.relpath(os.path.dirname(mp), T))
    need(worst[0] <= 1.0, "max dt/dt_lim over the FACE-LIST decks/legs = %.4f at %s (> 1 is a STOP: the check aborts"
         " there)" % (worst[0], worst[1]))
    warn = sorted(set(warn))
    print("     WARN (not on the face list; the check warns): %d decks/legs over the limit:" % len(warn))
    for w in warn:
        print("        %s%s" % (w, "  (retiring in this packet)" if retiring(w) else ""))
    live_warn = [w for w in warn if not retiring(w)]
    need(len(live_warn) == 7, "%d live (non-retiring) decks warn: the 7 the user's decision named (%s)"
         % (len(live_warn), ", ".join(live_warn)))
    print("     %d microstructure decks/legs checked" % nmicro)
    # the packet's spot value: alpha(297.5 K) = 1.6746e-6 m^2/s => 0.2430 s at 1.5625 mm
    a = grimsel_alpha({}, 0, 297.5)
    need(abs(a - 1.6746e-6) < 5e-10 and abs(0.0015625 ** 2 / (6 * a) - 0.2430) < 5e-4,
         "spot check: Grimsel alpha(297.5 K) = %.5e m^2/s, dt_lim(1.5625 mm) = %.4f s (packet 1.6746e-6, 0.2430)"
         % (a, 0.0015625 ** 2 / (6 * a)))
    a200 = grimsel_alpha({}, 0, 200.0)
    print("     (the table floor would give alpha(200 K) = %.4e => %.4f s at 1.5625 mm -- never used)"
          % (a200, 0.0015625 ** 2 / (6 * a200)))


# ------------------------------------------------------------- 3. schedule ----
SCHED_FILE = os.path.join(T, "validation", "oglesby", "fig34_schedule.txt")


def read_schedule_file(p):
    pairs = []
    for ln in open(p):
        s = ln.split("#", 1)[0].strip()
        if not s:
            continue
        a, b = s.split()
        pairs.append((a, b))
    return pairs


def section3():
    print("\n== 3. schedule proof (L4)")
    og = os.path.join(T, "validation", "oglesby", "input")
    pre = os.path.join(HERE, "src_pre_1", "tests", "MMWSpalling", "validation", "oglesby", "input")
    if os.path.exists(pre):
        live = {}
        for ln in open(og):
            m = LINE.match(ln)
            if m:
                live[m.group(1)] = m.group(2).strip()
        print("     the live deck now reads beam.schedule.file = %s; the arrays below are the PRE-IMAGE's (src_pre_1/)"
              % live.get("beam.schedule.file", "(not set)"))
        og = pre
    kv = {}
    for ln in open(og):
        m = LINE.match(ln)
        if m:
            kv[m.group(1)] = m.group(2).strip()
    if "beam.schedule.t_end" in kv:
        t_tok, p_tok = kv["beam.schedule.t_end"].split(), kv["beam.schedule.P"].split()
        need(len(t_tok) == len(p_tok) == 159, "oglesby arrays: %d t_end, %d P (159 each)" % (len(t_tok), len(p_tok)))
        tv = [float(x) for x in t_tok]
        need(all(b > a for a, b in zip(tv, tv[1:])), "oglesby t_end strictly increasing (%g .. %g)" % (tv[0], tv[-1]))
        need(all(float(x) == int(float(x)) for x in t_tok + p_tok),
             "every token is an integer (text -> double is exact on both parsers)")
        fig = open(os.path.join(T, "studies", "p2_21_loss_table", "doe", "FIG34_FULL.txt")).read().split("\n")
        srow = [ln for ln in fig if ln.startswith("the schedule (P_use")][0]
        fig_p = srow.split(":", 1)[1].split()
        need(fig_p == p_tok[:158], "the first 158 P equal FIG34_FULL.txt's schedule (%d values), the 159th is %s at"
             " t_end %s" % (len(fig_p), p_tok[158], t_tok[158]))
        if os.path.exists(SCHED_FILE):
            pairs = read_schedule_file(SCHED_FILE)
            need(pairs == list(zip(t_tok, p_tok)), "fig34_schedule.txt: %d pairs, equal as TEXT to the deck's arrays"
                 % len(pairs))
        else:
            print("     fig34_schedule.txt not written yet (stage 1 writes it from these arrays)")
    else:
        need(os.path.exists(SCHED_FILE) and "beam.schedule.file" in kv,
             "oglesby deck already carries beam.schedule.file (stage 1 applied)")
    print("  -- live schedule decks")
    n = 0
    for p in live_decks():
        kvd = resolve_deck(p)
        if "beam.schedule.t_end" not in kvd and "beam.schedule.file" not in kvd:
            continue
        n += 1
        if "beam.schedule.t_end" in kvd:
            tv = [float(x) for x in kvd["beam.schedule.t_end"].split()]
            pv = kvd.get("beam.schedule.P", "").split()
        else:
            pr = read_schedule_file(os.path.join(ROOT, kvd["beam.schedule.file"]))
            tv, pv = [float(a) for a, _ in pr], [b for _, b in pr]
        mono = all(b > a for a, b in zip(tv, tv[1:]))
        need(mono and len(tv) == len(pv), "%-66s %3d intervals, monotone %s" % (os.path.relpath(p, T), len(tv), mono))
    own = [p for p in live_decks() if any(LINE.match(l) and LINE.match(l).group(1).startswith("beam.schedule.") for l in open(p, errors="replace"))]
    need(len(own) == 12, "%d live decks SET beam.schedule.* (packet: 12); %d resolve one (unit/face_beam/input includes oglesby)" % (len(own), n))


# --------------------------------------------------------------- 4. headers ----
def ncomp(hdr):
    L = open(hdr).read().split("\n")
    nc = int(L[1])
    return nc, L[2:2 + nc]


def section4():
    print("\n== 4. L6: robin_face form_* plotfile headers (the `none` default keeps them)")
    for d in sorted(glob.glob(os.path.join(T, "unit", "robin_face", "output", "*"))):
        h = os.path.join(d, "00000cell", "Header")
        if not os.path.exists(h):
            continue
        nc, names = ncomp(h)
        tag = os.path.basename(d)
        print("     %-24s ncomp %2d  Sp_field %s" % (tag, nc, "Sp_field" in names))
        if tag.startswith("form_"):
            need(nc == 23 and "Sp_field" not in names, "%s: 23 components, no Sp_field (spallation.model default)" % tag)
    kv = resolve_deck(os.path.join(T, "unit", "robin_face", "input_form"))
    need(kv.get("spallation.model", "damage_law") == "damage_law",
         "robin_face/input_form takes the default spallation.model (%s)" % kv.get("spallation.model", "default"))


# ------------------------------------------------------------- 5. deck edits ----
ABORT_KEYS = [r"damage\..*", r"cohesive\..*", r"spall\.alpha", r"spall\.C_h", r"spall\.t_spall",
              r"spall\.damage_threshold", r"spall\.xi_s", r"spall\.shallow_depth", r"spall\.prescribed_ratio\..*"]
TABLE_KEYS = {"damage.enabled": "0", "conductivity.kappa_damage_alpha": "0.0", "spall.sample": "top_cell"}
GRIM6 = ["kappa", "beta", "E", "mu", "rho", "Cp"]
RETIRING = ["unit/dp_yield", "unit/gb_cohesive", "unit/spall_event", "unit/regime_low_high_power",
            "unit/heterogeneous_kappa", "validation/hu/hu_end_to_end", "validation/hu/hu_spall_onset",
            "extensions/hu_spall_onset_mmwbeam", "extensions/hu_spall_onset_voronoi"]
LIVE_TESTS = None


def retiring(rel):
    return any(rel.startswith(r.rstrip("/") + "/") for r in RETIRING)


def section5():
    print("\n== 5. deck-edit list (L13): DEFAULT = no behaviour; LIVE = named")
    a_d = {"DEFAULT": [], "LIVE": []}
    abort_set, table_set, model_dl, grim = [], [], [], []
    for p in live_decks() + glob.glob(os.path.join(T, "validation", "shared", "*.inputs")):
        rel = os.path.relpath(p, T)
        own = {}
        for i, ln in enumerate(open(p, errors="replace").read().split("\n"), 1):
            m = LINE.match(ln)
            if m:
                own[m.group(1)] = (i, m.group(2).strip())
        eff = resolve_deck(p) if not p.endswith(".inputs") else {k: v for k, (i, v) in own.items()}
        for k, (i, v) in own.items():
            if re.fullmatch(r"microstructure\.phase\d+\.A_D", k):
                dmg = a2.same(eff.get("damage.enabled", "0"), "1")
                a_d["DEFAULT" if (a2.same(v, "0.0") or not dmg) else "LIVE"].append(
                    (rel, i, k, v + ("" if a2.same(v, "0.0") else " (inert: damage.enabled off)")))
            if any(re.fullmatch(x, k) for x in ABORT_KEYS) and k not in TABLE_KEYS:
                abort_set.append((rel, i, k, v))
            if k in TABLE_KEYS:
                table_set.append((rel, i, k, v, "DEFAULT" if a2.same(v, TABLE_KEYS[k]) else "LIVE"))
            if k == "spallation.model" and v == "damage_law":
                model_dl.append((rel, i))
            mm = re.fullmatch(r"microstructure\.phase(\d+)\.(kappa|beta|E|mu|rho|Cp)", k)
            if mm:
                law = eff.get("microstructure.phase%s.law" % mm.group(1), "constant")
                if law == "grimsel":
                    grim.append((rel, i, k, v))
    print("  -- microstructure.phase<i>.A_D lines (stage 2: optional, must equal 0.0)")
    for tag in ("DEFAULT", "LIVE"):
        rows = a_d[tag]
        print("     %s: %d lines in %d files" % (tag, len(rows), len(set(r[0] for r in rows))))
        for r in rows:
            if tag == "LIVE" or len(rows) < 12:
                print("        %s:%d %s = %s %s" % (r[0], r[1], r[2], r[3], "(retiring)" if retiring(r[0]) else ""))
    need(all(retiring(r[0]) for r in a_d["LIVE"]), "every LIVE A_D line (non-zero AND damage.enabled = 1) is on a retiring deck")
    nz = [r for r in a_d["DEFAULT"] if "inert" in r[3]]
    print("     of the DEFAULT lines, %d are NON-ZERO but inert (damage off) -- the packet said non-zero only on retired decks:" % len(nz))
    for r in nz:
        print("        %s:%d %s = %s %s" % (r[0], r[1], r[2], r[3], "(retiring)" if retiring(r[0]) else "<- EDIT (delete: inert)"))
    print("  -- live setters of keys that ABORT after stage 2 (damage.*, cohesive.*, spall.{alpha,C_h,...})")
    for r in abort_set:
        print("        %s:%d %s = %s %s" % (r[0], r[1], r[2], r[3], "(retiring)" if retiring(r[0]) else "<- EDIT"))
    keep = [r for r in abort_set if not retiring(r[0])]
    need(all(r[0].startswith("unit/amr_microstructure_regrid") for r in keep),
         "non-retiring setters are only amr_microstructure_regrid's (LIVE: its damage IC, deleted with its D checks)"
         " -- %d line(s)" % len(keep))
    print("  -- setters of the three table keys (retired at default: run when equal)")
    for r in table_set:
        print("        %s:%d %s = %s %s%s" % (r[0], r[1], r[2], r[3], r[4], " (retiring)" if retiring(r[0]) else ""))
    need(all(r[4] == "DEFAULT" or retiring(r[0]) or r[0].startswith("unit/amr_microstructure_regrid")
             for r in table_set), "every live table-key setter is at its default")
    print("  -- explicit spallation.model = damage_law (aborts after stage 2)")
    for r in model_dl:
        print("        %s:%d %s" % (r[0], r[1], "(retiring)" if retiring(r[0]) else "<- EDIT"))
    need(all(retiring(r[0]) for r in model_dl), "explicit damage_law only on retiring decks")
    print("  -- the six grimsel-inert scalars on a law = grimsel phase (stage 4: refused when present)")
    for r in grim:
        print("        %s:%d %s = %s" % r)
    files = sorted(set(r[0] for r in grim))
    need(files == ["unit/face_list/input", "unit/face_removal/input", "validation/shared/rock_grimsel.inputs"],
         "grimsel scalar setters = the three packet files (%s)" % files)
    print("  -- argv in live test scripts naming a retiring key")
    for t in sorted(glob.glob(os.path.join(T, "unit", "*", "test")) + glob.glob(os.path.join(T, "validation", "*", "*", "test*"))
                    + glob.glob(os.path.join(T, "extensions", "*", "test*"))):
        rel = os.path.relpath(t, T)
        txt = open(t, errors="replace").read()
        hits = sorted(set(re.findall(r"(damage\.[A-Za-z_.]+|cohesive\.[A-Za-z_.]+|spall\.(?:alpha|C_h|t_spall|"
                                     r"damage_threshold|xi_s|shallow_depth|sample|prescribed_ratio\.\w+)|"
                                     r"spallation\.model=damage_law|kappa_damage_alpha|phase\d\.A_D)", txt)))
        if hits:
            print("        %-60s %s %s" % (rel, " ".join(hits), "(retiring)" if retiring(rel) else ""))


# --------------------------------------------------------------- 6. movers ----
HASH_PAT = re.compile(r"(Level_0/Cell_D_[^/]*$)|(_removal_events\.csv$)|(_h_col_events\.csv$)|(_clusters\.csv$)|"
                      r"(_jet_profile\.csv$)|(^|/)thermo\.dat$")


def read_md5(p):
    out = {}
    for ln in open(p):
        h, rel = ln.rstrip("\n").split("  ", 1) if "  " in ln else ln.rstrip("\n").split(" ", 1)
        out[rel.strip()] = h
    return out


def section6():
    print("\n== 6. mover lists (declared before any build)")
    jd = os.path.join(T, "unit", "jet_d2c", "output")
    allf = [os.path.relpath(os.path.join(dp, f), T) for dp, dn, fs in os.walk(jd)
            if ".old" not in dp for f in fs]
    hashed = sorted(f for f in allf if HASH_PAT.search(f))
    open(os.path.join(HERE, "MOVERS_1_jet_d2c.txt"), "w").write("\n".join(hashed) + "\n")
    print("  stage 1 only-in-new: unit/jet_d2c/output under hash.sh's patterns: %d of %d files on disk"
          " (MOVERS_1_jet_d2c.txt; ⚠ the sweep re-runs the test, so the count is the re-run's)" % (len(hashed), len(allf)))
    ref = read_md5(os.path.join(W, "p2_26.md5"))
    se = sorted(r for r in ref if "unit/spall_event/" in r)
    need(len(se) == 8, "stage 2 only-in-ref: spall_event's %d paths (packet: 8) %s" % (len(se), se))
    cd = [r for r in ref if "/Level_0/Cell_D_" in r]
    other = [r for r in ref if r not in cd]
    print("  stage 3 tree: %d of %d paths are Cell_D (every one moves: three components removed); %d others"
          " (CSVs + thermo) must not move; minus spall_event's %d; plus jet_d2c's Cell_D"
          % (len(cd), len(ref), len(other), len(se)))
    need(len(cd) == 2466, "tree Cell_D count %d (packet: 2,466)" % len(cd))
    fref = read_md5(os.path.join(W, "face_p2_24.md5"))
    fcd = [r for r in fref if "/Cell_D_" in r]
    fch = [r for r in fref if r.endswith("/Cell_H")]
    fhd = [r for r in fref if r.endswith("/Header")]
    print("  stage 3 face: %d files; %d Cell_D + %d Cell_H + %d Header move (components gone); thermo/CSV/"
          "Checkpoint must not" % (len(fref), len(fcd), len(fch), len(fhd)))
    # stage 4: the t = 0 plotfile of every grimsel face leg
    t0 = sorted(r for r in fref if "/00000cell/" in r and ("/Cell_D_" in r or r.endswith("/Cell_H")))
    legs = sorted(set(r.split("/run/")[0] for r in t0))
    grim_legs = []
    for lg in legs:
        mp = glob.glob(os.path.join(W, "face_p2_26", lg, "run", "metadata"))
        kv = read_metadata(mp[0]) if mp else {}
        if kv.get("microstructure.phase0.law") == "grimsel":
            grim_legs.append(lg)
    need(grim_legs == legs, "stage 4: every face leg is law = grimsel (%d of %d)" % (len(grim_legs), len(legs)))
    open(os.path.join(HERE, "MOVERS_4_face_t0.txt"), "w").write("\n".join(t0) + "\n")
    print("  stage 4 face movers: %d paths = the t = 0 Cell_D + Cell_H of %d legs (MOVERS_4_face_t0.txt); every"
          " later plotfile, Header, thermo and CSV identical" % (len(t0), len(legs)))
    print("  stage 4 tree movers: none (L8: no tree-gate deck is law = grimsel -- checked below)")
    tree_grim = []
    for d in TREE_OUT:
        for mp in glob.glob(os.path.join(T, d, "**", "metadata"), recursive=True):
            kv = read_metadata(mp)
            if any(kv.get("microstructure.phase%d.law" % i) == "grimsel" for i in range(4)):
                tree_grim.append(os.path.relpath(mp, T))
    need(not tree_grim, "no tree-gate leg runs law = grimsel (%s)" % tree_grim[:5])


# ---------------------------------------------------------------- 7. units ----
def section7():
    print("\n== 7. units (L12) and archives")
    for k, u in (("weibull.V0", "m^3"), ("dt_lim", "s"), ("alpha (diffusivity)", "m^2/s"), ("h_applied", "m"),
                 ("foot_mech_E", "J"), ("foot_mech_vol", "m^3"), ("conductivity.h_gb0", "W/(m^2 K)"),
                 ("kappa (phase conductivity)", "W/(m K)")):
        print("     %-28s [%s]" % (k, u))
    bp = os.path.join(HERE, "bin_pre", "mmwspalling-3d-g++")
    if os.path.exists(bp):
        h = hashlib.sha256(open(bp, "rb").read()).hexdigest()
        need(h.startswith("e0551036"), "bin_pre/ is the P2-26 binary (%s)" % h[:16])
    else:
        need(False, "bin_pre/ missing")
    sp = os.path.join(HERE, "src_pre_1")
    if os.path.isdir(sp):
        for f in sorted(glob.glob(os.path.join(sp, "**", "*"), recursive=True)):
            if os.path.isfile(f):
                rel = os.path.relpath(f, sp)
                live = os.path.join(ROOT, rel)
                same = os.path.exists(live) and open(f, "rb").read() == open(live, "rb").read()
                print("     src_pre_1/%-50s %s" % (rel, "== src" if same else "DIFFERS from the tree (edited since)"))
    else:
        print("     src_pre_1/ not yet copied")


if __name__ == "__main__":
    print("P2-27 z0 -- %s" % ROOT)
    section1()
    section2()
    section3()
    section4()
    section5()
    section6()
    section7()
    print("\nz0: %s (%d failed bar(s))" % ("PASS" if not BAD else "FAIL", len(BAD)))
    for b in BAD:
        print("  FAILED: " + b)
    sys.exit(1 if BAD else 0)
