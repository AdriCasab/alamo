#!/usr/bin/env python3
"""Census helper (read-only): classify an MMWSpalling input deck + argv overrides by lineage.

usage: classify.py <input> [key=val ...]
Follows `FILE = path` includes (ParmParse), last occurrence wins, argv last.
Prints one line: path|removal|lineages|closed-switches|aborts.
"""
import os
import re
import sys

ROOT = "/Users/tzetze20/amr_tools/alamo"
HERE = os.path.dirname(os.path.abspath(__file__))
SP_KEYS = set(open(os.path.join(HERE, "sp_keys.txt")).read().split())

P210 = ["spall.lateral_drive", "spall.lateral_removal", "spall.lateral_subcell", "spall.lateral_thermal",
        "spall.per_cell_removal", "surface_patch.side_face_flux", "surface_patch.side_face_h",
        "surface_patch.side_face_T", "surface_patch.side_face_h_r_min", "surface_patch.side_face_h_law",
        "surface_patch.side_face_h_datum", "surface_patch.side_face_factor_above",
        "surface_patch.side_face_area_law", "surface_patch.side_face_step_height",
        "surface_patch.side_face_staircase", "surface_patch.step_diag_radius",
        "surface_patch.jet_wall_march", "surface_patch.jet_wall_mdot", "surface_patch.jet_body",
        "surface_patch.jet_body_T", "surface_patch.jet_body_h_scale"]
P216 = ["spallation.sp.constraint_depth"]
DAMAGE_SPALL = ["spall.prescribed_ratio", "spall.alpha", "spall.t_spall", "spall.xi_s", "spall.C_h",
                "spall.damage_threshold"]


def read(path, kv, depth=0):
    if not os.path.isabs(path):
        cand = os.path.join(ROOT, path)
        path = cand if os.path.exists(cand) else path
    for raw in open(path, errors="replace"):
        line = raw.split("#", 1)[0].strip()
        if "=" not in line:
            continue
        k, v = line.split("=", 1)
        k, v = k.strip(), v.strip().strip('"')
        if k == "FILE" and depth < 5:
            read(v, kv, depth + 1)
            continue
        kv[k] = v


def on(kv, k):
    v = kv.get(k, "0").strip().lower()
    return v not in ("0", "false", "")


def classify(inp, argv):
    kv = {}
    read(inp, kv)
    for a in argv:
        if "=" in a:
            k, v = a.split("=", 1)
            kv[k.strip()] = v.strip()
    # thermal path
    if "microstructure.nphases" in kv:
        mode = kv.get("microstructure.mode", "voronoi")
        path = "MICRO-" + ("vor" if mode == "voronoi" else "expr") + "(n%s)" % kv["microstructure.nphases"]
    elif "material.type" in kv:
        path = ("MAT-IMPL" if kv.get("solver", "explicit") == "implicit" else "MAT-EXPL") + "(%s)" % kv["material.type"]
    else:
        path = "CONST-ALPHA"
    lin = []
    model = kv.get("spallation.model", "damage_law")
    spall = on(kv, "spall.enabled")
    vapor = on(kv, "vapor.enabled")
    fl = on(kv, "surface_patch.face_list")
    removal = "none"
    if spall or vapor:
        removal = "FACE" if fl else "COLUMN"
    # damage lineage
    dmg = []
    if on(kv, "damage.enabled"):
        dmg.append("damage.enabled")
    if "damage.ic.type" in kv:
        dmg.append("damage.ic")
    if on(kv, "cohesive.enabled"):
        dmg.append("cohesive")
    if spall and model == "damage_law":
        dmg.append("spall@damage_law")
    if kv.get("spall.sample") in ("column_max", "shallow_max"):
        dmg.append("sample=" + kv["spall.sample"])
    for k in DAMAGE_SPALL:
        if k in kv:
            dmg.append(k)
    if any(k.startswith("conductivity.kappa_damage_alpha") for k in kv) and kv.get("conductivity.kappa_damage_alpha", "0") not in ("0", "0.0"):
        dmg.append("kappa_damage_alpha")
    if dmg:
        lin.append("DAMAGE[" + ",".join(dmg) + "]")
    if model == "sp_weibull":
        lin.append("LEFM" + ("+weibull" if on(kv, "weibull.enabled") else ""))
    # flame / patch
    if on(kv, "surface_patch.enabled"):
        m = kv.get("surface_patch.mode", "prescribed_T")
        if m == "convective_flame":
            extra = []
            if kv.get("surface_patch.nozzle_descent") == "feet":
                extra.append("feet")
            if kv.get("surface_patch.jet_closure"):
                extra.append("jet=" + kv["surface_patch.jet_closure"])
            extra.append("robin=" + kv.get("surface_patch.robin_form", "cell"))
            lin.append("FLAME(" + ",".join(extra) + ")")
        else:
            lin.append("PATCH-prescribed_T")
    # beam
    p0 = kv.get("beam.P0", "0")
    try:
        bon = float(re.sub(r"[^0-9.eE+-].*$", "", p0) or 0) != 0.0
    except ValueError:
        bon = True
    if bon or "beam.schedule.P" in kv:
        lin.append("BEAM")
    if path.startswith("MICRO-vor") or "microstructure.grain_expr" in kv or on(kv, "weibull.enabled") and False:
        lin.append("MICRO-voronoi")
    elif path.startswith("MICRO-expr") and kv.get("microstructure.nphases", "1") not in ("1",):
        lin.append("MICRO-multiphase")
    if "microstructure.grain_expr" in kv:
        lin.append("MICRO-grain_expr")
    # stress source
    if spall or model == "sp_weibull" or on(kv, "damage.enabled") or on(kv, "cohesive.enabled"):
        ss = kv.get("spall.stress_source", "mechanics")
        lin.append("MECH" if ss == "mechanics" else "LTE")
    elif kv.get("el.type", "") not in ("", "disable"):
        lin.append("MECH(el)")
    if path.startswith("MAT") or path == "CONST-ALPHA":
        lin.append("ZHANG/global-path" if path.startswith("MAT") else "CONST-ALPHA-path")
    if vapor:
        lin.append("VAPOR")
    if on(kv, "bit.enabled"):
        lin.append("BIT")
    # closed study switches (C.7.3) still parsed
    cl = []
    if on(kv, "spall.surface_normal"):
        cl.append("surface_normal")
    if on(kv, "spall.melt_aware"):
        cl.append("melt_aware")
    if kv.get("spall.flake_coherence_length", "0") not in ("0", "0.0"):
        cl.append("flake_coherence_length")
    if kv.get("spall.flake_length") == "flux":
        cl.append("flake_length=flux")
    if kv.get("spall.flake_fraction") == "weibull":
        cl.append("flake_fraction=weibull")
    if kv.get("surface_patch.pad_shield_open", "-1") not in ("-1", "-1.0"):
        cl.append("pad_shield_open")
    if on(kv, "surface_patch.jet_above_h"):
        cl.append("jet_above_h")
    if kv.get("surface_patch.jet_T_ent_mode") == "rec_fixed":
        cl.append("jet_T_ent_mode=rec_fixed")
    if on(kv, "surface_patch.enabled") and kv.get("surface_patch.mode") == "convective_flame" and \
            kv.get("surface_patch.robin_form", "cell") in ("cell", "face"):
        cl.append("robin_form=" + kv.get("surface_patch.robin_form", "cell"))
    if kv.get("weibull.detach_mode") == "connected_cluster":
        cl.append("connected_cluster")
    if "spall.seed_void" in kv:
        cl.append("seed_void")
    try:
        if int(kv.get("amr.max_level", "0")) > 0:
            cl.append("amr.max_level=%s" % kv["amr.max_level"])
    except ValueError:
        pass
    ab = [k for k in P210 if k in kv] + [k for k in P216 if k in kv]
    unk = [k for k in kv if k.startswith("surface_patch.") and k.split(".", 1)[1] not in SP_KEYS
           and k not in P210]
    abort = ""
    if ab:
        abort = "ABORTS(retired:" + ",".join(ab) + ")"
    elif unk:
        abort = "ABORTS?(strict surface_patch unknown:" + ",".join(unk) + ")"
    return path, removal, " ".join(lin), " ".join(cl), abort


if __name__ == "__main__":
    r = classify(sys.argv[1], sys.argv[2:])
    print(" | ".join(r))
