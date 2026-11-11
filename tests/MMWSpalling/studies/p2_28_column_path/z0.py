#!/usr/bin/env python3
"""P2-28 z0 (offline gate) -- everything step zero must establish before the first edit.

  1. the packet's Sources re-resolved against src_pre_1/ (the code BEFORE any edit): label, symbol, cited -> found; a
     missing anchor FAILS (the inventory drifted)
  2. L2: no existing face-gate leg (witness/face_p2_27/, its `metadata`, last occurrence wins) has a prescribed patch
     or removal off -- so the stage-1 branches cannot reach a gate leg
  3. the removed-component and removed-column sets per face leg (stage 2's audit declaration): every leg's cell-plotfile
     Header and its thermo.dat header (face_p2_27_headers.txt)
  4. the test classification (KEEP / RETIRE / MOVE / RECORD / EDIT) with the SILENT-INERTNESS check: every deck's
     effective keys (FILE = includes resolved, last definition wins) -> its fate after stage 2 (ABORTS, naming the
     rule, or RUNS). A retiring deck that would RUN and stays in place is a FAIL (SR 1: an onset deck would run as a
     thermal-only deck with its field gone); a KEEP deck that would ABORT is a FAIL (a survivor would break)
  5. the deck-edit list (stage 2): every live setter of a key stage 2 retires, DEFAULT or LIVE
  6. the new decks' explicit limits (stage 1): dt vs dx^2/(n alpha), n = 2d (+1 with a prescribed top)
  7. units (L10); bin_pre/ = the P2-27 binary; src_pre_1/ equal to src/
Usage: z0.py > Z0.txt   (exit 1 on a failed bar)
"""
import glob
import hashlib
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
T = os.path.join(ROOT, "tests", "MMWSpalling")
W = os.path.join(T, "studies", "tree_sweep_0922", "witness")
SRC = {"M": "src/Integrator/MMWSpalling.H", "R": "src/Integrator/MMWSpalling/Removal.H",
       "F": "src/Integrator/MMWSpalling/FaceRemoval.H", "S": "src/Integrator/MMWSpalling/Shared.H"}
LINE = re.compile(r"^\s*([A-Za-z0-9_.]+)\s*=\s*(.*?)\s*(#.*)?$")
BAD = []


def need(ok, msg):
    print("  %s  %s" % ("ok  " if ok else "FAIL", msg))
    if not ok:
        BAD.append(msg)


def lines(k):
    pre = os.path.join(HERE, "src_pre_1", SRC[k])
    return open(pre if os.path.exists(pre) else os.path.join(ROOT, SRC[k])).read().split("\n")


# ------------------------------------------------------------------ 1. sources ----
# (label, file key, regex, cited line). Regexes name a SYMBOL or a distinctive phrase, never a line number.
ANCHORS = [
    # the dispatch
    ("Advance: UpdateRemovalAfterCohesive call", "M", r"^\s*UpdateRemovalAfterCohesive\(lev, time, dt\);", 1767),
    ("Advance: UpdateSpClusters call", "M", r"^\s*UpdateSpClusters\(lev", 1763),
    ("Advance: UpdateSpAfterMechanics call", "M", r"^\s*UpdateSpAfterMechanics\(lev", 1759),
    ("Advance: BuildSurfaceTColumns call", "M", r"^\s*BuildSurfaceTColumns\(lev", 1749),
    ("R: UpdateRemovalAfterCohesive", "R", r"void UpdateRemovalAfterCohesive\(", 25),
    ("R: RemovalStepBegins call", "R", r"if \(!RemovalStepBegins\(lev\)\) return;", 27),
    ("R: ResetRemovalStepFields call", "R", r"^\s*ResetRemovalStepFields\(lev\);", 79),
    ("R: column clearance guard !face_list", "R", r"!surface_patch\.face_list", 145),
    ("R: k_top_global = ColumnTopSolid", "R", r"k_top_global = ColumnTopSolid\(lev\)", 222),
    ("R: THE DISPATCH RecedeAlongLines", "R", r"RecedeAlongLines\(lev, time, dt, k_top_global\);", 230),
    ("S: RemovalStepBegins", "S", r"bool RemovalStepBegins\(int lev\)", 340),
    ("S: RemovalStepBegins !removal_enabled return", "S", r"if \(!removal_enabled\(\)\) return false;", 342),
    ("S: ResetRemovalStepFields", "S", r"void ResetRemovalStepFields\(int lev\)", 351),
    ("S: reset regime_mf", "S", r"regime_mf\[lev\]->setVal\(0\.0\);", 353),
    ("M: ColumnTopSolid", "M", r"amrex::Vector<int> ColumnTopSolid\(int lev\)", 4792),
    ("F: ElasticConstants in the scan", "F", r"const ElasticConsts ec = ElasticConstants\(\);", 670),
    ("F: conf_off = ec.conf", "F", r"conf_off = ec\.conf", 671),
    ("F: ScanTAtZv", "F", r"return ScanTAtZv\(T_at_d,", 793),
    ("F: ScanSigmaT (criterion)", "F", r"return ScanSigmaT\(scan_g, T_z, T0_z, conf_off", 803),
    ("M: SurfaceGasTempCap", "M", r"Set::Scalar SurfaceGasTempCap\(", 8602),
    ("M: RawSurfaceX", "M", r"Set::Scalar RawSurfaceX\(", 8501),
    ("M: ReconFromX", "M", r"ReconFromX\(Set::Scalar", 8516),
    ("M: MeltFraction", "M", r"Set::Scalar MeltFraction\(Set::Scalar", 5084),
    ("M: UpdateSurfaceMaskFromPhi", "M", r"void UpdateSurfaceMaskFromPhi\(int lev\)", 6616),
    # column-only in M (stage 2)
    ("site A UpdateSpAfterMechanics", "M", r"void UpdateSpAfterMechanics\(int lev", 5844),
    ("UpdateSpClusters", "M", r"void UpdateSpClusters\(int lev", 6172),
    ("BuildSurfaceTColumns", "M", r"void BuildSurfaceTColumns\(int lev", 3163),
    ("SharedSurfaceRecon", "M", r"SharedSurfaceRecon\(", 3255),
    ("SurfaceRobinT", "M", r"Set::Scalar SurfaceRobinT\(", 3277),
    ("BuildPinEffective", "M", r"void BuildPinEffective\(int lev", 3292),
    ("PinRule", "M", r"PinRule\(", 3314),
    ("PatchDiagnostics", "M", r"void PatchDiagnostics\(int lev", 3338),
    ("PatchDiagnostics call", "M", r"^\s*PatchDiagnostics\(lev, k_top_col\);", 2435),
    ("GatherColumnSurface", "M", r"void GatherColumnSurface\(int lev", 5802),
    ("ColumnSolidRunBelowTop", "M", r"ColumnSolidRunBelowTop\(", 5154),
    ("RepairRemovalStateAfterRegrid", "M", r"void RepairRemovalStateAfterRegrid\(int lev\)", 4551),
    ("RepairRemovalStateAfterRegrid call", "M", r"RepairRemovalStateAfterRegrid\(lev\);", 1489),
    ("CheckSurfaceCovered", "M", r"void CheckSurfaceCovered\(", 4644),
    ("surface band in TagCellsForRefinement", "M", r"removal_enabled\(\) && surface_refine_depth > 0\.0", 5508),
    ("ParseBitSettings", "M", r"static void ParseBitSettings\(", 5771),
    ("ParseBitSettings call", "M", r"^\s*ParseBitSettings\(value, pp\);", 345),
    ("JetEnthalpyMarch", "M", r"void JetEnthalpyMarch\(int lev", 2778),
    ("march fl_on_m", "M", r"const bool fl_on_m", 2916),
    ("spall.surface_T parse", "M", r'"spall\.surface_T"', 5594),
    ("spall.removal_energy parse", "M", r'"spall\.removal_energy"', 5657),
    ("DetachMode enum", "M", r"enum class DetachMode", 8390),
    ("sp_emit_clusters_csv", "M", r"bool\s+sp_emit_clusters_csv", 8406),
    ("sp_emit_h_col_events_csv", "M", r"bool\s+sp_emit_h_col_events_csv", 8407),
    ("prescribed_stress parse", "M", r'"spallation\.sp\.prescribed_stress\.enabled"', 8730),
    ("prescribed_T_surface parse", "M", r'"spallation\.sp\.prescribed_T_surface\.enabled"', 8742),
    ("robin_form enum", "M", r"enum class RobinForm \{ Cell, Face, Pinned, Surface \};", 7539),
    ("robin_form default Cell", "M", r"RobinForm robin_form = RobinForm::Cell;", 7545),
    ("pinned_idle_cycles member", "M", r"Set::Scalar pinned_idle_cycles = 2\.0;", 7546),
    # registrations leaving the face legs (stage 2)
    ("reg Sp_field", "M", r'"Sp_field"', 663),
    ("reg flaw_a (STAYS)", "M", r'"flaw_a"', 669),
    ("reg Sp_cluster_id", "M", r'"Sp_cluster_id"', 674),
    ("reg h_spall_field", "M", r'"h_spall_field"', 683),
    ("reg crack_event_count", "M", r'"crack_event_count"', 686),
    ("reg T_surf_field", "M", r'"T_surf_field"', 696),
    ("reg ell_eff_field", "M", r'"ell_eff_field"', 698),
    ("reg recon_cols", "M", r'"recon_cols"', 703),
    ("reg surf_cols_robin", "M", r'"surf_cols_robin"', 768),
    ("reg surf_q_mean", "M", r'"surf_q_mean"', 769),
    ("reg surf_side_share (STAYS)", "M", r'"surf_side_share"', 770),
    ("reg patch_min_standoff (STAYS)", "M", r'"patch_min_standoff"', 809),
    ("reg patch_P_robin", "M", r'"patch_P_robin"', 810),
    # evaluate-only (stage 1)
    ("reg phi/removed/regime/RoP guard removal_enabled", "M", r"^\s*if \(value\.removal_enabled\(\)\)$", 615),
    ("follow_mask requires removal (LIFT)", "M", r'"surface\.follow_mask = 1 requires spall\.enabled "', 1120),
    ("follow_mask refuses prescribed_T (LIFT)", "M", r'"surface\.follow_mask = 1 is not supported with "', 1124),
    ("face list requires convective+h_expr (LIFT)", "M",
     r"sp\.enabled && !\(sp\.mode == SurfacePatch::Mode::ConvectiveFlame && sp\.use_expr\)", 1287),
    ("face list refuses removal off (LIFT)", "M",
     r"value\.spallation_model != SpallationModel::SpWeibull \|\| !value\.spall_enabled", 1334),
    ("kernel void_aware = removal_enabled()", "M", r"const bool void_aware = removal_enabled\(\);", 2159),
    ("fl_on refuses prescribed (LIFT)", "M", r"&& !\(surface_patch\.enabled && sp_presc\);", 2237),
    ("F: DepthScanH call", "F", r"o\.h_col = DepthScanH\(nprof", 828),
    ("F: o.fire = (o.h_col > 0)", "F", r"o\.fire = \(o\.h_col > 0\.0\) \? 1 : 0;", 832),
    ("F: buckling record BEFORE the withdrawal", "F", r"P2-23: the record, taken BEFORE the withdrawal", 881),
    ("F: withdrawal o.fire = 0", "F", r"^\s*o\.fire = 0;", 898),
    ("F: bk_log guard (thermo.plot_int stride)", "F", r"const bool bk_log = bk_on && sp_buckle_log_csv && thermo\.plot_int > 0", 703),
    ("F: win[] per-cell winner", "F", r"std::vector<int> win\(cells\.size\(\), -1\);", 959),
    ("F: ONE WRITE", "F", r"---- ONE WRITE: f, removed, H, phi", 1117),
    ("F: fr_score on every listed cell", "F", r"fsc\(c\.i, c\.j, c\.k\) = \(l >= 0\)", 1202),
    ("F: event log", "F", r"---- the event log: one row per firing", 1212),
    ("F: FrEventsHeader", "F", r"static std::string FrEventsHeader\(bool bk_on\)", 330),
    ("F: FrAssignTs", "F", r"void FrAssignTs\(int lev\)", 162),
    ("F: FrFaceFlux", "F", r"void FrFaceFlux\(const FaceEntry& e", 138),
    # prescribed patch (stage 1)
    ("kernel sp_T = surface_patch.T_at(time)", "M", r"\? surface_patch\.T_at\(time\) : 0\.0;", 2229),
    ("kernel prescribed term (factor 1, CORRECT)", "M", r"div_kgrad \+= kc \* \(sp_T - Tc\) \* inv_dz2;", 2347),
    ("kernel listed-face bill FrFaceFlux", "M", r"FrFaceFlux\(flp\[q\], fTg \? flp\[q\]\.Tg", 2382),
    ("ledger prescribed mirror (factor 1, CORRECT)", "M", r"P_robin \+= kc \* \(sp_T - Tc\) \* inv_dz2 \* V;", 4361),
    ("ledger listed-face mirror", "M", r"FrFaceFlux\(e, e\.Tg, kc, Tc, DX\[0\], a_in, a_out, a_rob, a_los\);", 4386),
    ("FaceListDiagnostics", "M", r"void FaceListDiagnostics\(int lev\)", 4075),
    ("FaceListHeat", "M", r"void FaceListHeat\(int lev, Set::Scalar z_n", 4010),
    ("FaceEntry", "M", r"struct FaceEntry$", 3440),
    ("BuildFlameColumns constant-flame early return", "M", r"if \(!surface_patch\.use_expr\)$", 2457),
    ("material kernel prescribed term (NOT touched, D)", "M", r"Dirichlet face value at distance dz/2", 1936),
    ("doc comment kappa*(T_patch - T_cell)/dz", "M", r"kappa\*\(T_patch - T_cell\)/dz", 7520),
    ("UpdateModel top-node override", "M", r"Override\s*$|Override the top-node temperature inside the patch", 1511),
    ("dt check CheckTimeStepLimit", "M", r"void CheckTimeStepLimit\(\)", 1434),
    # confining.p
    ("confining.p parse", "M", r'"confining\.p", value\.confining_p', 8641),
    ("S: c.conf = p nu/(1 - nu)", "S", r"c\.conf = confining_p \* c\.nu / \(1\.0 - c\.nu\);", 56),
]


def section1():
    print("== 1. Sources re-resolved against src_pre_1/ (cited -> found; drift reported; a missing anchor FAILS)")
    cache = {k: lines(k) for k in SRC}
    for label, fk, rx, cited in ANCHORS:
        r = re.compile(rx)
        hits = [i + 1 for i, ln in enumerate(cache[fk]) if r.search(ln)]
        if not hits:
            need(False, "%-50s %s:%d -> NOT FOUND (/%s/)" % (label, fk, cited, rx))
            continue
        best = min(hits, key=lambda h: abs(h - cited))
        drift = best - cited
        extra = "" if len(hits) == 1 else " (%d hits: %s)" % (len(hits), ",".join(map(str, hits[:6])))
        need(True, "%-50s %s:%d -> %s:%d%s%s" % (label, fk, cited, fk, best,
                                                 "" if drift == 0 else " (drift %+d)" % drift, extra))
    # the test-side anchors (files, not code)
    V = os.path.join(T, "validation")
    tchecks = [
        ("Kant Eq. 15 constants (sp_kant_onset/test)", os.path.join(V, "kant", "sp_kant_onset", "test"),
         [r"TOL_EQ15\s*=\s*0\.08", r"I_FK\s*=\s*1\.763", r"KIC0_PA\s*=\s*1\.5e6"]),
        ("pressure sweep slope (sp_kant_pressure_sweep/test)", os.path.join(V, "kant", "sp_kant_pressure_sweep", "test"),
         [r"confining\.p=", r"0\.08"]),
        ("Hu deck prescribed_T_surface (input_granite2_lefm)", os.path.join(V, "hu", "hu_spall_onset", "input_granite2_lefm"),
         [r"^spallation\.sp\.prescribed_T_surface\.enabled = 1", r"^spall\.buckle_span = 0\.0",
          r"^# check_shared: allow spall\.buckle_span", r"^surface_patch\.T_expr\s+= 638\.22\+3\.74\*t"]),
        ("Hu test reads Sp_field (test_lefm)", os.path.join(V, "hu", "hu_spall_onset", "test_lefm"), [r"Sp_field"]),
        ("reduce.py onset statistic (first dir 5, ts_fb 0)", os.path.join(V, "kant", "grimsel_calibration", "reduce.py"),
         [r"dir", r"ts_fb"]),
        ("planner's 1-D Dirichlet check", os.path.join(ROOT, "Claude_markdowns", "scripts", "p2_28_dirichlet_1d.txt"), [r"."]),
    ]
    for label, p, rxs in tchecks:
        if not os.path.exists(p):
            need(False, "%s: %s missing" % (label, os.path.relpath(p, ROOT)))
            continue
        txt = open(p, errors="replace").read()
        miss = [rx for rx in rxs if not re.search(rx, txt, re.M)]
        need(not miss, "%s: %s%s" % (label, os.path.relpath(p, ROOT), (" MISSING " + str(miss)) if miss else ""))


# ------------------------------------------------------------- deck resolution ----
def resolve_deck(path, argv=()):
    kv = {}

    def walk(p):
        for ln in open(p, errors="replace").read().split("\n"):
            m = LINE.match(ln)
            if not m:
                continue
            k, v = m.group(1), m.group(2).strip()
            if k == "FILE":
                walk(os.path.join(ROOT, v))
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
        kv[k.strip()] = v.strip()
    return kv


def on(kv, k):
    try:
        return float(kv.get(k, "0").split()[0]) != 0.0
    except ValueError:
        return False


FACE_LEGS = ["fr/T00_1_k1", "fr/TL35_1", "fr/TL45_1", "fr/SPH_1", "fr/BORE_1", "fl/CAV", "fl/STAIR", "fb/F_S",
             "fh/FLAT_P", "fh/BORE_P", "M0S", "HOS_S", "CAL_S", "BEAM_S", "FEET_S", "FEET_0"]


def section2():
    print("\n== 2. L2: every face-gate leg (face_p2_27, metadata) has removal ON and NO prescribed patch")
    for leg in FACE_LEGS:
        mp = os.path.join(W, "face_p2_27", leg, "run", "metadata")
        if not os.path.exists(mp):
            need(False, "%-14s metadata missing" % leg)
            continue
        kv = read_metadata(mp)
        en = on(kv, "surface_patch.enabled")
        mode = kv.get("surface_patch.mode", "prescribed_T")
        rem = on(kv, "spall.enabled") or on(kv, "vapor.enabled")
        presc = en and mode != "convective_flame"
        need(rem and not presc and on(kv, "surface_patch.face_list"),
             "%-14s removal %s, patch %s mode %s, face_list %s" % (leg, rem, en, mode, kv.get("surface_patch.face_list")))


# ---------------------------------------------------- 3. removed sets per leg ----
REMOVED_FIELDS = ["Sp_field", "Sp_cluster_id", "h_spall_field", "crack_event_count", "T_surf_field", "ell_eff_field"]
REMOVED_COLS_RX = re.compile(r"^(recon_\w+|patch_P_robin|surf_cols_robin|surf_q_mean|pinned_\w+)$")


def section3():
    print("\n== 3. stage 2's declared removed sets per face leg (cell Header components; thermo columns)")
    hdr = {}
    for ln in open(os.path.join(W, "face_p2_27_headers.txt")):
        rel, h = ln.rstrip("\n").split("\t", 1)
        hdr[rel] = [c.strip() for c in h.split("\t")]
    for leg in FACE_LEGS:
        hs = sorted(glob.glob(os.path.join(W, "face_p2_27", leg, "run", "*cell", "Header")))
        if not hs:
            need(False, "%-14s no cell Header" % leg)
            continue
        comps = [ln.strip() for ln in open(hs[-1])][2:2 + int(open(hs[-1]).read().split("\n")[1])]
        rf = [c for c in comps if c in REMOVED_FIELDS]
        th = hdr.get(os.path.join(leg, "run", "thermo.dat"), [])
        rc = [c for c in th if REMOVED_COLS_RX.match(c)]
        print("     %-14s ncomp %3d  fields -%d %s" % (leg, len(comps), len(rf), rf))
        print("     %-14s ncol  %3d  cols   -%d %s; ledger_err %s" % ("", len(th), len(rc), rc,
                                                                      "present" if "ledger_err" in th else "absent"))
        need(set(rf) >= {"Sp_field", "h_spall_field", "crack_event_count", "T_surf_field", "ell_eff_field"},
             "%-14s the five everywhere-fields present (Sp_cluster_id %s)" % (leg, "Sp_cluster_id" in rf))


# --------------------------------------- 4. test classification + silent inertness ----
RETIRE_UNITS = ["robin_pinned", "robin_jet", "jet_enthalpy", "robin_feet", "jet_d2c", "pad_clear_ledger", "pin_clock",
                "removal_increment", "scalar_flaw", "beam_void_closure", "seed_void", "surface_T_recon",
                "sp_v_n_regime", "sp_lefm"]
MOVE_VALID = ["kant/sp_kant_onset/input_verification", "kant/sp_kant_onset/input_verification_hfl",
              "kant/sp_kant_onset/input_model", "kant/sp_kant_onset/input_weibull",
              "kant/sp_kant_onset/input_weibull_cluster", "kant/sp_kant_pressure_sweep/input",
              "kant/sp_kant_onset_meiercfg", "p0_steady_column", "rossi"]
RECORDS = ["studies/s1_surface_resolution/input", "validation/meier/sp_meier_pilot/input_2d_dev",
           "validation/meier/sp_meier_pilot/input_2d_dev_melt", "validation/meier/sp_meier_pilot/input_energy",
           "validation/meier/sp_meier_pilot/input_drilling", "validation/meier/sp_meier_pilot/input_drilling_25mm",
           "validation/meier/sp_meier_pilot/input_feet", "validation/meier/sp_meier_pilot/input_feet_d2c",
           "extensions/hu_spall_onset_mmwbeam_lefm"]
EDIT_UNITS = ["face_removal", "face_list", "face_buckle", "face_heat", "face_beam", "gradient_blind",
              "lateral_flux_lock", "buckle_lock", "interlocks", "retired_keys"]
ABORT_PRESENCE = ["spall.clusters_csv", "spall.h_col_events_csv", "weibull.A_crit", "surface.refine_depth",
                  "surface.refine_height"]
ABORT_PREFIX = ["bit.", "spallation.sp.prescribed_stress.", "spallation.sp.prescribed_T_surface."]
AT_DEFAULT = {"spall.surface_T": "reconstructed", "spall.removal_energy": "flake", "weibull.detach_mode": "per_face"}


def fate(kv):
    """the deck's fate after stage 2: ('ABORT', reason) or ('RUNS', note)"""
    rem = on(kv, "spall.enabled") or on(kv, "vapor.enabled")
    fl = on(kv, "surface_patch.face_list")
    if rem and not fl:
        return "ABORT", "removal without the face list (one removal path)"
    if kv.get("surface_patch.robin_form", "cell") == "pinned":
        return "ABORT", "robin_form = pinned"
    for k in kv:
        if k in ABORT_PRESENCE or any(k.startswith(p) for p in ABORT_PREFIX):
            return "ABORT", "retired key %s" % k
    for k, v in AT_DEFAULT.items():
        if k in kv and kv[k] != v:
            return "ABORT", "%s = %s (retired at %s)" % (k, kv[k], v)
    if "surface.follow_mask" in kv and not on(kv, "surface.follow_mask") and rem:
        return "ABORT", "surface.follow_mask = 0 with removal"
    if rem and kv.get("surface_patch.robin_form", "cell") in ("cell", "face"):
        return "ABORT", "robin_form cell/face with removal (the face list requires surface)"
    onset = kv.get("spallation.model", "none") == "sp_weibull" and not rem and not fl
    return "RUNS", ("SILENT: sp_weibull onset deck without removal or the face list -- Sp_field gone"
                    if onset else "")


def decks_in(rel):
    p = os.path.join(T, rel)
    if os.path.isfile(p):
        return [p]
    return sorted(q for q in glob.glob(os.path.join(p, "**", "input*"), recursive=True)
                  if os.path.isfile(q) and not q.endswith(".inputs") and ".old" not in q and "/output" not in q)


# the decks the packet EDITS (their fate is judged on the edited keys): stage 1 the Hu deck (item 3), stage 2 the
# declared deck edits (item 8)
EDITED = {
    "validation/hu/hu_spall_onset/input_granite2_lefm": (
        "stage 1 item 3", ["spallation.sp.prescribed_T_surface.enabled", "spallation.sp.prescribed_T_surface.expr",
                           "spall.buckle_span"],
        {"surface_patch.face_list": "1", "surface_patch.robin_form": "surface", "surface.follow_mask": "1",
         "spall.surface_T": "reconstructed", "spall.removal_energy": "flake", "spall.removal_events_csv": "1",
         "spall.buckle_log_csv": "1"}),
    "unit/face_feet/input": ("stage 2 item 8", ["spall.h_col_events_csv"], {}),
    "validation/meier/sp_meier_pilot/input_feet_face": ("stage 2 item 8", ["bit.*"], {}),
}


def edited_kv(rel, kv):
    _, drop, add = EDITED[rel]
    out = {k: v for k, v in kv.items()
           if not any(k == d or (d.endswith("*") and k.startswith(d[:-1])) for d in drop)}
    out.update(add)
    return out


def section4():
    print("\n== 4. test classification + silent inertness (fate after stage 2, from the deck's own keys; argv-built "
          "legs inside a test script are listed in section 5)")
    seen = set()

    def show(cls, rel, must):
        for d in decks_in(rel):
            seen.add(os.path.abspath(d))
            try:
                f, why = fate(resolve_deck(d))
            except (OSError, KeyError) as e:
                f, why = "?", str(e)
            ok = {"RUNS_OK": f == "RUNS" and not why.startswith("SILENT"), "ABORT": f == "ABORT",
                  "ANY": True}[must]
            need(ok, "%-7s %-62s %-5s %s" % (cls, os.path.relpath(d, T), f, why))

    for u in RETIRE_UNITS:
        show("RETIRE", "unit/" + u, "ANY")
    show("RETIRE", "unit/robin_face/input_ledger", "ANY")
    show("RETIRE", "unit/gradient_blind/input", "ANY")
    for v in MOVE_VALID:
        show("MOVE", "validation/" + v, "ANY")
    for r in RECORDS:
        show("RECORD", r, "ABORT")
    for u in EDIT_UNITS:
        if u == "gradient_blind":
            continue
        show("EDIT", "unit/" + u, "ANY")
    # every other live deck is KEEP: it must RUN (no new abort) and must not become a silent onset deck
    others = []
    for base in ("unit", "validation", "extensions"):
        for q in glob.glob(os.path.join(T, base, "**", "input*"), recursive=True):
            if not os.path.isfile(q) or q.endswith(".inputs") or ".old" in q or "/output" in q:
                continue
            if "retired_" in q or "/legs/" in q:
                continue
            if os.path.abspath(q) in seen:
                continue
            others.append(q)
    for q in sorted(others):
        rel = os.path.relpath(q, T)
        if "/pre/" in rel:
            # an archived pre-image (records: a header line at most); it must not run silently
            try:
                f, why = fate(resolve_deck(q))
            except (OSError, KeyError) as e:
                f, why = "?", str(e)
            need(not why.startswith("SILENT"), "%-7s %-62s %-5s %s" % ("RECORD", rel, f, why))
            continue
        if rel in EDITED:
            kv = edited_kv(rel, resolve_deck(q))
            f, why = fate(kv)
            need(f == "RUNS" and not why.startswith("SILENT"),
                 "%-7s %-62s %-5s %s (edited, %s)" % ("EDITED", rel, f, why, EDITED[rel][0]))
            continue
        try:
            f, why = fate(resolve_deck(q))
        except (OSError, KeyError) as e:
            print("     ?       %-62s unresolvable (%s)" % (rel, e))
            continue
        tag = "KEEP"
        if rel.startswith("validation/kant/sp_kant_onset/") or rel.startswith("validation/kant/sp_kant_pressure"):
            tag = "NEW?"
        need(f == "RUNS" and not why.startswith("SILENT"), "%-7s %-62s %-5s %s" % (tag, rel, f, why))


# -------------------------------------------------------------- 5. deck edits ----
RETIRED_STAGE2 = ["surface_patch.pinned_idle_cycles", "spall.h_col_events_csv", "spall.clusters_csv",
                  "spall.surface_T", "spall.removal_energy", "surface.follow_mask", "weibull.detach_mode",
                  "weibull.A_crit", "surface.refine_depth", "surface.refine_height"]
TEST_TOKENS = [r"face_list=0", r"robin_form=pinned", r"robin_form=cell", r"robin_form=face", r"\bbit\.",
               r"h_col_events_csv", r"clusters_csv", r"prescribed_T_surface", r"prescribed_stress",
               r"detach_mode=connected_cluster", r"surface_T=cell", r"removal_energy=cell", r"follow_mask=0",
               r"Removal\.H", r"Sp_field", r"pinned_idle_cycles"]


def section5():
    print("\n== 5. stage 2's deck-edit list: every live setter of a retired key (DEFAULT = the value it retires at)")
    files = []
    for base in ("unit", "validation", "extensions", "studies/tree_sweep_0922/witness"):
        files += [q for q in glob.glob(os.path.join(T, base, "**", "*"), recursive=True)
                  if os.path.isfile(q) and (os.path.basename(q).startswith("input") or q.endswith(".inputs"))
                  and ".old" not in q and "/output" not in q and "retired_" not in q and "/legs/" not in q]
    n = 0
    for q in sorted(files):
        for i, ln in enumerate(open(q, errors="replace").read().split("\n")):
            m = LINE.match(ln)
            if not m:
                continue
            k, v = m.group(1), m.group(2).strip()
            if k in RETIRED_STAGE2 or k.startswith("bit.") or k.startswith("spallation.sp.prescribed_"):
                dflt = {"spall.surface_T": "reconstructed", "spall.removal_energy": "flake",
                        "surface.follow_mask": "1", "weibull.detach_mode": "per_face"}.get(k)
                cls = "DEFAULT" if (dflt is not None and v == dflt) else "LIVE"
                print("     %-7s %s:%d  %s = %s" % (cls, os.path.relpath(q, T), i + 1, k, v))
                n += 1
    print("     (%d lines)" % n)
    print("  -- argv tokens in KEEP/EDIT test scripts that stage 2 makes abort (the scripts' legs to drop or edit)")
    for u in EDIT_UNITS + ["face_feet", "face_vapor", "face_melt", "robin_face", "convective_patch",
                           "amr_microstructure_regrid", "latent_heat", "beam", "beam_tau", "beam_schedule",
                           "profile_depth"]:
        p = os.path.join(T, "unit", u, "test")
        if not os.path.exists(p):
            continue
        txt = open(p, errors="replace").read().split("\n")
        hits = {}
        for rx in TEST_TOKENS:
            r = re.compile(rx)
            ls = [i + 1 for i, ln in enumerate(txt) if r.search(ln)]
            if ls:
                hits[rx] = ls[:8]
        if hits:
            print("     %-26s %s" % (u, "; ".join("%s @ %s" % (k, ",".join(map(str, v))) for k, v in hits.items())))


# ---------------------------------------------------------------- 6. dt ----
def section6():
    print("\n== 6. stage 1's new decks: dt vs the explicit limit dx^2/(n alpha), n = 6 (+1 on a prescribed top)")
    rows = [
        # name, dx [m], alpha [m^2/s] at T_cold, dt [s], prescribed
        ("Kant face (constant law: 2.5/(2630*800))", 1.0e-3, 2.5 / (2630.0 * 800.0), 0.1, False),
        ("face_dirichlet coarse (2.5/(2640*827.2))", 3.125e-3, 2.5 / (2640.0 * 827.2), 0.5, True),
        ("face_dirichlet fine", 1.5625e-3, 2.5 / (2640.0 * 827.2), 0.125, True),
        ("Hu face (Grimsel at 280.15 K: 1.8145e-6)", 3.125e-3, 1.8145e-6, 0.5, True),
    ]
    for name, dx, a, dt, presc in rows:
        n = 7 if presc else 6
        lim = dx * dx / (n * a)
        need(dt <= lim, "%-44s dx %.6g dt %.4g lim dx^2/(%d a) = %.4g ratio %.3f" % (name, dx, dt, n, lim, dt / lim))


# --------------------------------------------------------------- 7. units, pre ----
def sha(p):
    h = hashlib.sha256()
    h.update(open(p, "rb").read())
    return h.hexdigest()


def section7():
    print("\n== 7. units (L10); bin_pre; src_pre_1")
    for u in ["T_presc [K]", "kc [W/(m K)]", "inv_dz2 [1/m^2]", "the Dirichlet term 2 kc (T_presc - Tc) inv_dz2 [W/m^3]",
              "confining.p [Pa]", "Eq. 15 Theta [K] (Delta T above T_I 293.15 K)", "slope [K/MPa]",
              "ledger_err_solid, ledger_err_raw [-] (relative to ledger_E_src)"]:
        print("     %s" % u)
    pre = os.path.join(HERE, "bin_pre", "mmwspalling-3d-g++")
    p27 = open(os.path.join(T, "studies", "p2_27_cleanup", "BINARY_5.txt")).read()
    need(os.path.exists(pre) and sha(pre) in p27, "bin_pre/ = the P2-27 stage-5 binary (%s)" % sha(pre)[:16])
    for root, _, fs in os.walk(os.path.join(HERE, "src_pre_1")):
        for f in fs:
            p = os.path.join(root, f)
            rel = os.path.relpath(p, os.path.join(HERE, "src_pre_1"))
            need(open(p, "rb").read() == open(os.path.join(ROOT, rel), "rb").read(), "src_pre_1/%s == %s" % (rel, rel))


if __name__ == "__main__":
    section1()
    section2()
    section3()
    section4()
    section5()
    section6()
    section7()
    print("\nz0: %s (%d failed)" % ("PASS" if not BAD else "FAIL", len(BAD)))
    sys.exit(1 if BAD else 0)
