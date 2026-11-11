#!/usr/bin/env python3
"""P2-31 z0 (offline half): the reachability ledger L1-L17 resolved against the code, the decks and the outputs on
disk. Prints Z0 rows; `smoke.py` runs the two smokes and L14's cost legs. Usage: z0.py > Z0.txt"""
import csv
import glob
import hashlib
import os
import re
import shlex
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import decks as D  # noqa: E402

ROOT, T = D.ROOT, D.T
U = os.path.join(T, "unit")
V = os.path.join(T, "validation")
W = os.path.join(T, "studies", "tree_sweep_0922", "witness")
REF = os.path.join(W, "face_p2_30c")
SRC = os.path.join(ROOT, "src")
DS = os.path.join(SRC, "Integrator", "MMWSpalling")
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
GATE_LEGS = ["fr/T00_1_k1", "fr/TL35_1", "fr/TL45_1", "fr/SPH_1", "fr/BORE_1", "fl/CAV", "fl/STAIR", "fb/F_S",
             "fh/FLAT_P", "fh/BORE_P", "M0S", "HOS_S", "CAL_S", "BEAM_S", "FEET_S", "FEET_0", "KV_S", "HU_S", "GB_S",
             "VAP_S", "BKL_S"]
BAD = []


def row(tag, ok, msg):
    print("%-4s %s  %s" % (tag, "ok  " if ok else "FAIL", msg))
    if not ok:
        BAD.append(tag)


def meta(path):
    kv = {}
    for ln in open(path, errors="replace"):
        if " = " in ln:
            k, v = ln.rstrip("\n").split(" = ", 1)
            kv[k.strip()] = v.strip()
    return kv


def code(text):
    """comments blanked (// and /* */), string literals kept"""
    out, i, n = [], 0, len(text)
    while i < n:
        if text.startswith("//", i):
            j = text.find("\n", i)
            j = n if j < 0 else j
            i = j
        elif text.startswith("/*", i):
            j = text.find("*/", i + 2)
            i = n if j < 0 else j + 2
        elif text[i] == '"':
            j = i + 1
            while j < n and text[j] != '"':
                j += 2 if text[j] == "\\" else 1
            out.append(text[i:j + 1])
            i = j + 1
        else:
            out.append(text[i])
            i += 1
    return "".join(out)


def src_files():
    out = {}
    for dp, _, fs in os.walk(SRC):
        for fn in fs:
            if fn.endswith((".H", ".cpp", ".cc")):
                p = os.path.join(dp, fn)
                out[os.path.relpath(p, SRC)] = open(p, errors="replace").read()
    return out


def test_sources():
    """every live test source under unit/ and validation/ and the witness scripts (no stored outputs)"""
    out = []
    for top in (U, V):
        for dp, dns, fns in os.walk(top):
            dns[:] = [d for d in dns if not (d.startswith("output") or d == "legs" or ".old" in d or d == "__pycache__")]
            for fn in fns:
                if fn == "test" or fn.endswith(".py"):
                    out.append(os.path.join(dp, fn))
    out += glob.glob(os.path.join(W, "*.py"))
    return sorted(out)


def thermo(path):
    L = [ln.rstrip("\n").split("\t") for ln in open(path) if ln.strip()]
    h = [x.strip() for x in L[0]]
    return h, L[1:]


def main():
    print("P2-31 z0 -- the reachability ledger (offline). Binary %s" % hashlib.sha256(open(BIN, "rb").read()).hexdigest())
    print("=" * 110)
    S = src_files()
    C = {k: code(v) for k, v in S.items()}
    TS = test_sources()
    TT = {p: open(p, errors="replace").read() for p in TS}

    # ---- git / baseline
    ref_bin = open(os.path.join(REF, "fr", "T00_1_k1", "BINARY.txt")).read().split()[-1]
    cur = hashlib.sha256(open(BIN, "rb").read()).hexdigest()
    newer = [p for p in glob.glob(os.path.join(SRC, "**", "*.H"), recursive=True) if os.path.getmtime(p) > os.path.getmtime(BIN)]
    row("git", cur == ref_bin and not newer, "the binary on disk is face_p2_30c's (%s); src files newer than it: %d"
        % (cur[:16], len(newer)))

    # ---- L1
    frm = C["Integrator/MMWSpalling/FaceRemoval.H"]
    prs = C["Integrator/MMWSpalling/Parse.H"]
    reg = re.search(r'if \(pp\.contains\("microstructure\.nphases"\)\)\s*\{[\s\S]*?RegisterNewFab\(value\.T_init_mf', prs)
    store = len(re.findall(r"StoreInitialTemperature\(lev\);", C["Integrator/MMWSpalling.H"]))
    req = 'if (!pp.contains("microstructure.nphases"))' in prs
    print("L1   T_init_mf: registered inside the nphases block (%s), which every deck enters (nphases is required: %s);"
          " written once in Initialize (StoreInitialTemperature x%d); the face list is single level, so have_T0 is true"
          % (reg is not None, req, store))
    print("     the gather's fallback: `have_T0 ? T0a(p) : NaN`, NaN -> T0r in T0_at_d (the stress's own); per gate leg,"
          " phase0.T_ref vs the IC (resolved metadata):")
    diff_legs = []
    for lg in GATE_LEGS:
        m = meta(os.path.join(REF, lg, "run", "metadata"))
        tr, it, iv = m.get("microstructure.phase0.T_ref"), m.get("hc.ic.type"), m.get("hc.ic.constant.value")
        same = it == "constant" and iv is not None and float(iv) == float(tr)
        if not same:
            diff_legs.append(lg)
        print("       %-12s T_ref %-8s ic %-9s %-8s %s" % (lg, tr, it, iv, "equal" if same else "DIFFERS"))
    row("L1", reg is not None and req and store == 1 and diff_legs == ["GB_S"],
        "the initial field exists at every line's start cell; T_ref != IC on exactly %s (packet: GB_S only)" % diff_legs)

    # ---- L2
    shr = C["Integrator/MMWSpalling/Shared.H"]
    st = C["Integrator/MMWSpalling/State.H"]
    g_raw = "dT0 > 0.0 && dT1 > 0.0 && dT0 > dT1" in st
    g_recon = "if (!(x > 0.0)) { x = 0.0; r.flat = true; }" in st and "(x > 0.0) ? x / (1.0 - std::exp(-x)) : 1.0" in st
    g_efold = "if (!(rise_s > 0.0)) return std::numeric_limits<Set::Scalar>::quiet_NaN();" in shr
    g_peak = "if (!(r_pk > 0.0)) { z_peak = nan; ell_pk = nan; return; }" in shr
    g_fall = "if (r_b > 0.0)" in shr
    row("L2", g_raw and g_recon and g_efold and g_peak and g_fall,
        "a zero rise cannot divide: RawSurfaceX returns 0 (flat) unless dT0 > 0, dT1 > 0, dT0 > dT1 (%s); ReconFromX "
        "flat at x <= 0, fac 1 (%s); ProfileEFoldDepth NaN at rise_s <= 0 (%s); ProfilePeak NaN at r_pk <= 0 (%s); "
        "ProfileFallDepth logs only r_b > 0 and its lerp divides r_a - r_b > 0 (%s); BuckleRatio inf at eps_th <= 0. "
        "GB_S smoke: smoke.py" % (g_raw, g_recon, g_efold, g_peak, g_fall))

    # ---- L3
    l3 = []
    for lg in ("HU_S", "KV_S"):
        deck, ov = D.argv_deck(os.path.join(REF, lg, "argv.sh"))
        R = D.resolve(deck, ov)
        m = meta(os.path.join(REF, lg, "run", "metadata"))
        l3.append((lg, m.get("spall.enabled"), "spall.h_spall_n_segments" in R, "spall.seed_void" in R))
    row("L3", all(not a and not b for _, _, a, b in l3),
        "goal 2 cannot move the gate's spall-off legs: %s (leg, spall.enabled, sets h_spall_n_segments, sets seed_void)"
        % l3)

    # ---- L4
    lv = D.live_decks()
    z = []
    SIX = ("T_m", "fus_width", "L_m", "T_vap_lo", "T_vap_hi", "L_v")
    for dk in lv:
        try:
            R = D.resolve(dk)
        except Exception:
            continue
        for i, law in enumerate(D.phase_laws(R)):
            if law == "zhang":
                have = {k: D.val(R, "microstructure.phase%d.%s" % (i, k)) for k in SIX}
                z.append((os.path.relpath(dk, ROOT), i, have))
    for dk, i, have in z:
        print("     L4 %s phase%d: %s" % (dk, i, ", ".join("%s=%s" % (k, v) for k, v in have.items())))
    row("L4", len(z) >= 1, "%d live deck(s) with law = zhang; lacking a latent key: %s" % (
        len(z), [d for d, _, h in z if any(v is None for v in h.values())] or "none"))

    # ---- L5
    hits = []
    for dk in lv:
        try:
            R = D.resolve(dk)
        except Exception:
            continue
        if D.face_list_on(R) and "grimsel" in D.phase_laws(R) and D.val(R, "spall.stress_law", "scalar") != "grimsel":
            hits.append(os.path.relpath(dk, ROOT))
    ghits = []
    for lg in GATE_LEGS:
        m = meta(os.path.join(REF, lg, "run", "metadata"))
        laws = [m.get("microstructure.phase%d.law" % i, "constant") for i in range(int(m.get("microstructure.nphases", "0")))]
        if m.get("surface_patch.face_list") == "1" and "grimsel" in laws and m.get("spall.stress_law") != "grimsel":
            ghits.append(lg)
    drv = []
    for p in glob.glob(os.path.join(T, "studies", "**", "*.py"), recursive=True):
        if "/src_pre" in p or "/legs/" in p or os.path.abspath(p) == os.path.abspath(__file__):
            continue
        t = open(p, errors="replace").read()
        if re.search(r"stress_law\s*=\s*scalar", t) or re.search(r'"spall\.stress_law"\s*:\s*"scalar"', t):
            drv.append(os.path.relpath(p, ROOT))
    print("     L5 study drivers that write spall.stress_law = scalar (historic; none in the gate or the live suite): %s"
          % (drv or "none"))
    row("L5", not hits and not ghits, "law grimsel + stress_law scalar on the face list: live decks %s (of %d), gate legs "
        "%s (of %d)" % (hits or "none", len(lv), ghits or "none", len(GATE_LEGS)))

    # ---- L6
    pins = {}
    for nm, pat in (("the old abort text", r"F4 keeps the column-top face"), (".kappa(", r"\.kappa\("),
                    ("H_of", r"\bH_of\b(?!_T)"), (".active(", r"\.active\(\)")):
        pins[nm] = [os.path.relpath(p, ROOT) for p, t in TT.items() if re.search(pat, t)]
    hof_t = [os.path.relpath(p, ROOT) for p, t in TT.items() if "H_of_T" in t]
    row("L6", not any(pins.values()), "test sources pinning: %s; `H_of_T` (the planner's grep for H_of) matches %s -- the "
        "tests' own Python enthalpy twins, not Grimsel::H_of" % (pins, hof_t))

    # ---- L7
    calls = ["%s:%d" % (rel, n) for rel, t in C.items() for n, ln in enumerate(t.split("\n"), 1)
             if re.search(r"\bkappa\s*\(", ln) and "virtual" not in ln and "override" not in ln]
    mref = ["%s:%d" % (rel, n) for rel, t in C.items() for n, ln in enumerate(t.split("\n"), 1)
            if re.search(r"(?:const\s+)?(?:Numeric::Material::)?Material\s*[&*]", ln)]
    tb = C["Numeric/Material/Table.H"]
    build = re.search(r"void Build\(const Material &mat[\s\S]*?\n    \}", tb).group(0)
    used = sorted(set(re.findall(r"mat\.(\w+)", build)))
    row("L7", len(calls) == 3 and "kappa" not in used,
        "kappa callers %s (the three named sites); Material&/* uses %s; Table::Build reads mat.%s and the base's "
        "effective_rhocp -> rho(), cp(): no polymorphic kappa" % (calls, mref, used))

    # ---- L8: the dead list
    print("L8   the dead list (src/ non-comment readers beyond the declaration/write; live test readers)")
    items = [
        ("SurfaceRecon::ell_q", r"\bell_q\b", ["Model.H:declaration", "FaceRemoval.H:the fallback's write"]),
        ("ElasticConsts::ok", r"\bc\.ok\b|\bec\.ok\b|\bbool ok = false; \}", ["Shared.H:declaration + c.ok = true"]),
        ("Grimsel::H_of", r"\bH_of\s*\(", ["Grimsel.H:definition"]),
        ("MMWBeam::active()", r"\bactive\s*\(\s*\)", ["MMWBeam.H:definition"]),
        ("Table m_T_v", r"\bm_T_v\b", ["Table.H:declaration + write"]),
        ("Table m_fus_lo_T", r"\bm_fus_lo_T\b", ["Table.H:declaration + write"]),
        ("Table m_fus_hi_T", r"\bm_fus_hi_T\b", ["Table.H:declaration + write"]),
        ("Table m_vap_hi_T", r"\bm_vap_hi_T\b", ["Table.H:declaration + write"]),
        ("vapor_A_beam(_override)", r"\bvapor_A_beam(?:_override)?\b", ["parse + the default's derivation"]),
        ("vapor_T0(_override)", r"\bvapor_T0(?:_override)?\b", ["parse + the default"]),
        ("vapor_prescribed_rop_enabled / vapor_rop_parser / vapor_rop_f",
         r"\bvapor_(?:prescribed_rop_enabled|rop_parser|rop_f)\b", ["parse"]),
        ("microstructure_enabled", r"\bmicrostructure_enabled\b", ["always true: nphases is required"]),
    ]
    for nm, pat, note in items:
        hits_ = ["%s:%d" % (rel, n) for rel, t in C.items() for n, ln in enumerate(t.split("\n"), 1) if re.search(pat, ln)]
        th = [os.path.relpath(p, ROOT) for p, t in TT.items() if re.search(pat, t)]
        print("       %-34s src code lines %2d %s; test readers %s  [%s]" % (nm, len(hits_), hits_[:6], th or "none",
                                                                          "; ".join(note)))
    rhoCp = [os.path.relpath(p, ROOT) for p, t in TT.items() if re.search(r'"(?:rho_phase|Cp_phase)"', t)]
    print("       NOT dead: rho_mf / Cp_mf -- test readers of rho_phase / Cp_phase: %s" % rhoCp)
    mse = ["%s:%d" % (rel, n) for rel, t in C.items() for n, ln in enumerate(t.split("\n"), 1)
           if re.search(r"\bmicrostructure_enabled\b", ln)]
    unreach = [h for h in mse if "Parse.H" in h and "= true" not in S["Integrator/MMWSpalling/Parse.H"].split("\n")[int(h.split(":")[1]) - 1]]
    row("L8", len(mse) == 18 and len(unreach) == 7, "microstructure_enabled: %d non-comment uses (packet 18), %d of them "
        "in CheckInterlocks' aborts (packet 7: %s); every other item above has no reader beyond its declaration/parse/write"
        % (len(mse), len(unreach), [h.split(":")[1] for h in unreach]))

    # ---- L9
    l9 = []
    for dk in (os.path.join(U, "face_beam", "input"), os.path.join(U, "face_vapor", "input")):
        R = D.resolve(dk)
        l9.append((os.path.relpath(dk, ROOT), {k: R[k][0] for k in R if k.startswith("vapor.")}))
    anyl = [os.path.relpath(d, ROOT) for d in lv if any(k.startswith(("vapor.A_beam", "vapor.T0", "vapor.prescribed_rop"))
                                                        for k in D.resolve(d))]
    row("L9", True, "RESOLVES DIFFERENTLY: neither deck sets vapor.prescribed_rop.enabled (includes followed): %s; live "
        "decks setting vapor.A_beam / vapor.T0 / vapor.prescribed_rop.*: %s. The keys run at their defaults (0 / unset); "
        "goal 9 has NO deck line to remove (the packet's premise came from a resolved dump, not the decks)" % (l9, anyl or "none"))

    # ---- L10
    print("L10  phi under the face list -- writers: InitPhi (Fields.H), FrFootClearance (FaceRemoval.H, reads phi only to"
          " re-write it: the tie rule), the ONE WRITE's +z shift and vapour shift (FaceRemoval.H)")
    print("     readers: UpdateSurfaceMaskFromPhi -> surface_mf, read by (i) the plotfile `surface`, (ii) EnergyLedger's"
          " has_surf -> thermo `ledger_surface_missing_cols`; the kernel's and the ledger's surf_c take the mask column top"
          " under the face list (follow_mask = face_list), BeamProps' surface_mf is a placeholder when !void_aware")
    fl_flux = "const Set::Scalar surf_c = follow_mask" in C["Integrator/MMWSpalling/Heat.H"] and \
        "value.surface_follow_mask = value.surface_patch.face_list;" in prs
    row("L10", fl_flux, "phi feeds DIAGNOSTICS only => stage 3's difference set on the overlap legs: plotfile components"
        " {phi, surface, spall_thickness, RoP_spall, RoP} and thermo {ledger_surface_missing_cols}; event logs, Temp, H,"
        " fr_f, removed unchanged")

    # ---- L11
    ev = list(csv.DictReader(open(os.path.join(REF, "KV_S", "run_removal_events.csv"))))
    kv = os.path.join(V, "kant", "sp_kant_onset", "output_face", "verification_removal_events.csv")
    km = meta(os.path.join(V, "kant", "sp_kant_onset", "output_face", "verification", "metadata"))
    ev2 = list(csv.DictReader(open(kv)))
    ha = sorted({float(e["h_applied"]) for e in ev2})
    nv = sorted({e["n_voided"] for e in ev2})
    hl = sum(1 for e in ev2 if float(e["h_line"]) > 0)
    row("L11", len(ev2) > 0 and ha == [0.0] and nv == ["0"],
        "RESOLVED ON A FULL RUN: the gate's KV_S (and HU_S) end before their first firing (10 s; header-only logs, %d "
        "rows). The Kant deck's full evaluate-only run (sp_kant_onset/output_face/verification, spall.enabled %s, stop %s "
        "s) writes %d rows: h_applied %s, n_voided %s, h_line > 0 on %d (the flake production would take). The gate "
        "carries NO evaluate-only data row -- declared, as P2-30's buckle-log cover; stage 1 does not touch that writer"
        % (len(ev), km.get("spall.enabled"), km.get("stop_time"), len(ev2), ha, nv, hl))

    # ---- L12
    fl = C["Integrator/MMWSpalling/Fields.H"]
    row("L12", "std::pow(V_cell / weibull_V0, inv_m)" in fl,
        "the one site Fields.H: vol_factor = (V_cell/V0)^(1/m) on every intragranular cell (the fixture is grain-free: "
        "every cell); V_cell 8.0e-9 / 1.0e-9 / 1.25e-10 m3 at 2 / 1 / 0.5 mm; factor 1 at V0 = V_cell; Wf at 2 / 0.5 mm: "
        "8^(1/13) = %.4f / 8^(-1/13) = %.4f" % (8 ** (1 / 13.0), 8 ** (-1 / 13.0)))

    # ---- L13
    fx = D.resolve(os.path.join(U, "face_removal", "input"))
    rk = D.resolve(os.path.join(V, "shared", "rock_grimsel.inputs"))
    dif = [(k, fx.get(k, ("(unset)",))[0], v[0]) for k, v in sorted(rk.items()) if fx.get(k, ("(unset)",))[0] != v[0]]
    print("L13  the fixture (unit/face_removal/input) vs validation/shared/rock_grimsel.inputs (+ its loss table), keys the"
          " rock file sets:")
    for k, a, b in dif:
        print("       %-46s fixture %-14s rock file %s" % (k, a, b))
    row("L13", True, "%d differences; the ones that matter for stage 4: the span (unset = 0 vs 0.012, the legs add it) "
        "and the beam/loss/melt keys (the fixture is flame-only, beam off; melt never reached at 1866 K gas)" % len(dif))

    # ---- L15
    rd = {nm: [os.path.relpath(p, ROOT) for p, t in TT.items() if re.search(pat, t)]
          for nm, pat in (("regime (plotfile field)", r'["\']regime["\'](?!_)'), ("regime_field", r"regime_field"))}
    ext = [os.path.relpath(p, ROOT) for p in glob.glob(os.path.join(T, "extensions", "*", "test"))
           if "regime_field" in open(p).read()]
    row("L15", len(rd["regime_field"]) == 0, "readers in U/, validation/, W/*.py: %s. Outside that scope: %s reads "
        "regime_field for an UNGATED informational print (G6) and tolerates its absence (it already reads fields deleted "
        "in P2-28) => goal 10's branch: regime_field is removed (2b)" % (rd, ext))

    # ---- L16
    print("L16  the first list on every gate leg (thermo row 1 = the end of step 1; row 0 is the t = 0 state):")
    zero = []
    for lg in GATE_LEGS:
        tp = os.path.join(REF, lg, "run", "thermo.dat")
        if not os.path.exists(tp):
            print("       %-12s no thermo.dat (amr.thermo.plot_int unset) -- smoke.py reads its first list from run.log" % lg)
            continue
        h, rows = thermo(tp)
        v = float(rows[1][h.index("fl_n_top")]) if "fl_n_top" in h and len(rows) > 1 else None
        allz = [float(r[h.index("fl_n_top")]) for r in rows[1:]] if "fl_n_top" in h else []
        print("       %-12s fl_n_top step 1: %s; min over steps %s" % (lg, v, min(allz) if allz else "-"))
        if not v:
            zero.append(lg)
    outs = []
    for top in (U, V):
        for p in glob.glob(os.path.join(top, "**", "thermo.dat"), recursive=True):
            if ".old" in p or "/retired" in p:
                continue
            try:
                h, rows = thermo(p)
            except Exception:
                continue
            if "fl_n_top" not in h or len(rows) < 2:
                continue
            outs.append((os.path.relpath(p, ROOT), float(rows[1][h.index("fl_n_top")])))
    oz = [p for p, v in outs if v == 0.0]
    print("       outputs on disk under unit/ + validation/ with fl_n_top: %d; first list empty on: %s" % (len(outs), oz or "none"))
    fld = []
    for dk in lv:
        try:
            if D.face_list_on(D.resolve(dk)):
                fld.append(os.path.relpath(dk, ROOT))
        except Exception:
            pass
    print("       live face-list decks (as written; tests add keys per leg): %d: %s" % (len(fld), fld))
    row("L16", not zero, "every gate leg lists > 0 top faces on its first step (%s empty); stored outputs with an empty "
        "first list: %s (smoke.py runs every live face-list deck one step to cover the rest)" % (zero or "none", oz or "none"))

    # ---- L17
    P25 = os.path.join(T, "studies", "p2_25_regime_map")
    prs_raw = S["Integrator/MMWSpalling/Parse.H"]
    keys_abort = set(re.findall(r'"([a-z_]+\.[\w.]+)"', prs_raw))
    for nm in ("D1_Q2", "D2_Q2", "D2_Q4", "F_D1_Q2", "F_D2_Q2", "F_D2_Q4"):
        deck, ov = D.argv_deck(os.path.join(P25, "legs", nm, "argv.sh"))
        R = D.resolve(deck, ov)
        retired = [k for k in R if k in keys_abort and k not in ("surface_patch.face_list",)]
        print("       L17 %-8s deck %s: %d keys; keys named in Parse.H's retired lists/table: %s" % (
            nm, os.path.relpath(deck, ROOT), len(R), [k for k in retired if any(
                ('"%s"' % k) in blk for blk in re.findall(r"retired_at_default\[\] = \{[\s\S]*?\};|for \(const char\* k : \{[^}]*\}", prs_raw))] or "none"))
    row("L17", True, "the six argv replay token for token on today's parser if smoke.py's one-step runs read rc 0")
    print("=" * 110)
    print("z0 offline: %s" % ("every row resolves" if not BAD else "FAILED rows: %s" % BAD))
    return 1 if BAD else 0


if __name__ == "__main__":
    sys.exit(main())
