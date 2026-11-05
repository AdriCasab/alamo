"""Render ../1_keys.md from keys.py + annot.py + readers.json + setters.json."""
import json, collections, os, re
import keys, annot
H = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(H, "readers.json")))
S = json.load(open(os.path.join(H, "setters.json")))
ABBR = {"MMWSpalling.H": "M", "Removal.H": "R", "FaceRemoval.H": "FR"}

def esc(t):
    return str(t).replace("|", "\\|").replace("\n", " ")
def uniq(vals, n=3):
    v = [x if len(x) < 40 else x[:37] + "..." for x in dict.fromkeys(vals)]
    return ",".join(v[:n]) + (",..." if len(v) > n else "")
def readers(key, kind):
    if kind == "retired":
        return "none (contains() only)"
    if kind == "grim":
        return "Grimsel.H methods (see notes)"
    h = [x for x in R[key] if x["tag"] == "read"]
    if not h:
        return "**none**"
    c = collections.OrderedDict()
    for x in h:
        c.setdefault(x["func"], []).append("%s:%d" % (ABBR[x["file"]], x["line"]))
    if key == "microstructure.nphases":
        return "%d reads of `phases`/`microstructure_enabled` in %d functions (path selector)" % (len(h), len(c))
    parts = []
    for f, ls in c.items():
        lst = ",".join(ls[:6]) + (",..." if len(ls) > 6 else "")
        parts.append("%s×%d (%s)" % (f, len(ls), lst))
    return "%d: " % len(h) + "; ".join(parts)
def prod(key):
    p = S[key]["prod"]
    out = []
    for leg in ("HOS_R", "CAL"):
        if p.get(leg) is not None:
            v = p[leg]
            if len(v) > 40: v = v[:37] + "..."
            out.append("%s=`%s`" % (leg, v))
    s = "; ".join(out) or "—"
    if p.get("CAL_all"): s += " (CAL file order: %s)" % "→".join(p["CAL_all"])
    if p.get("HOS_R_all"): s += " (HOS_R order: %s)" % "→".join(p["HOS_R_all"])
    return s
def units(key):
    d = collections.OrderedDict()
    for x in S[key]["unit"]:
        d.setdefault(x["file"].split("/")[1], []).append(x["val"])
    if not d: return "—"
    return "%d: " % len(d) + " ".join("%s(%s)" % (k, uniq(v)) for k, v in d.items())
def vals(key):
    d = collections.OrderedDict()
    for x in S[key]["val"]:
        f = x["file"].replace("validation/", "")
        parts = f.split("/")
        g = "/".join(parts[:-1]) if len(parts) > 1 else f
        d.setdefault(g, []).append((parts[-1], x["val"]))
    if not d: return "—"
    nfiles = len(set(x["file"] for x in S[key]["val"]))
    out = []
    for g, lst in d.items():
        files = list(dict.fromkeys(f for f, _ in lst))
        out.append("%s×%d(%s)" % (g, len(files), uniq([v for _, v in lst])))
    return "%d files: " % nfiles + " ".join(out)
def gate(key, default):
    lst = [x for x in S[key]["gate"] if x["nd"] is not False]
    if not lst: return "—"
    d = collections.OrderedDict()
    for x in lst:
        f = (x["file"].replace("unit/", "").replace("validation/meier/sp_meier_pilot/", "meier/")
             .replace("studies/s1_surface_resolution/", "s1/"))
        d.setdefault(f, []).append(x["val"])
    req = default.upper().startswith("REQUIRED")
    return ("(required) " if req else "") + " ".join("%s(%s)" % (k, uniq(v)) for k, v in d.items())

CLASSES = ["PROD", "FLAME-CAP", "BEAM", "ZHANG", "MICRO", "MECH", "DAMAGE", "CLOSED", "DIAG", "DEAD", "RETIRED-ABORT"]
rows = collections.defaultdict(list)
for r in keys.K:
    key = r["key"]
    if r["kind"] == "retired":
        cls, guard, pst, note = ("RETIRED-ABORT", "-", "not set (would abort)",
            ("P2-10 retired list. spall.* is top-level, so this list is the key ONLY guard (comment :89-96)" if key.startswith("spall.") else "P2-10 retired list. surface_patch.* is a strict queryclass scope, so the key would abort anyway as specified-but-not-used; the list only gives the better message") if key != "spallation.sp.constraint_depth" else "P2-16 retired list (constraint depth deleted; spall.buckle_span kept)")
    else:
        cls, guard, pst, note = annot.A[key]
    mem = ", ".join("`%s`" % m for m in r["members"]) or "—"
    if r["kind"] == "sp": mem = mem.replace("`", "`SurfacePatch::", 1) if r["members"] else "— (consumed, no member)"
    if r["kind"] == "ph": mem = "`MineralPhase::" + ", ".join(r["members"]) + "`"
    if r["kind"] == "grim": mem = "`Grimsel::%s`" % r["members"][0]
    if r["kind"] == "beam": mem = "`MMWBeam::%s`" % r["members"][0]
    parse = "; ".join(r["lines"])
    rd = readers(key, r["kind"])
    if key in ("beam.P0", "beam.schedule.t_end", "beam.schedule.P"):
        rd = ("7 (indirect, via MMWBeam::power_at, MMWBeam.H:93-98): AdvanceConstAlpha (M:1477); AdvanceMaterial (M:1518); "
              "AdvanceMaterialImplicit (M:1650); AdvanceMicrostructure (M:1826); EnergyLedger (M:4226); BeamEnergyBalance (M:4762); "
              "UpdateRemovalAfterCohesive (R:327, vapour only)")
    if r["kind"] == "grim":
        if key.endswith("grimsel.T_ref"):
            rd = "**none** (see notes)"
        else:
            rd, note = note.split(". Strict", 1)
            note = "Strict" + note
    rows[cls].append([key, parse, r["default"], mem, rd, r["guard"] + (" ‖ readers: " if r["guard"] else "readers: ") + guard if guard != "-" else (r["guard"] or "-"),
                      prod(key) + " — " + pst, units(key), vals(key), gate(key, r["default"]), note])

# grouped rows (select<> sub-keys) -- written by hand, see work/README
EXTRA = {
 "ZHANG": [
  ["material.zhang.{rho.a,rho.b,rho.cap,rho.T_cap,kappa.a,kappa.b,kappa.T_cap,cp.a_lo,cp.b_lo,cp.c_lo,cp.a_hi,cp.b_hi,cp.T_break,T_ref,T_m,fus_width,L_m,T_vap_lo,T_vap_hi,L_v} (20 keys)",
   "Numeric/Material/Zhang.H:39-64 via select_default M:180", "Zhang.H defaults", "`Zhang::*`", "via `material` (Table build M:6629, AdvanceMaterial*)",
   "material.type = zhang (default type once the block exists)", "—", "—", "zhang/zhang_oglesby uses type only (all defaults)", "—", "only reachable when material.type is given (contains() M:179)"],
  ["material.constant.{rho,cp,kappa (REQUIRED), T_ref,T_m,fus_width,L_m,T_vap_lo,T_vap_hi,L_v} (10 keys)",
   "Numeric/Material/Constant.H:33-45 via select_default M:180", "REQUIRED / Constant.H defaults", "`Constant::*`", "via `material`",
   "material.type = constant", "—", "convective_patch", "hu/{hu_breakage_index×3, hu_conduction×2, hu_thermoelastic×2, hu_thermoelastic_amr×2} (rho,cp,kappa,T_ref)", "—", ""],
 ],
 "DAMAGE": [
  ["damage.ic.{constant.value | expression.coord, .unit, .region<n>, .constant.*}",
   "IC/Constant.H:71, IC/Expression.H:150-190 via select M:496", "IC defaults", "`ic_damage`", "InitializeMicrostructure (damage IC)",
   "microstructure.nphases present & damage.ic.type given", "—", "4: amr_microstructure_regrid heterogeneous_kappa regime_low_high_power spall_event (expression)", "—", "spall_event/input, input_no_event (expression)", ""],
 ],
}
for c, lst in EXTRA.items():
    rows[c] += lst

HDR = ["Key", "Parse site", "Default", "Member(s) set", "Readers outside Parse (n: Func×count (file:line))", "Guard (parse branch ‖ reader path)",
       "(a) Production (effective value) — status", "(b) Kept unit tests (dirs)", "(c) Validation decks", "(d) Tree gate: non-default setters", "Notes / evidence"]
out = []
counts = {c: len(rows[c]) for c in CLASSES}
out.append("COUNTS " + json.dumps(counts))
for c in CLASSES:
    ng = len(EXTRA.get(c, []))
    out.append("\n### %s (%d keys%s)\n" % (c, len(rows[c]) - ng, (" + %d grouped rows" % ng) if ng else ""))
    out.append("| " + " | ".join(HDR) + " |")
    out.append("|" + "---|" * len(HDR))
    for row in sorted(rows[c], key=lambda x: x[0].lower()):
        out.append("| " + " | ".join(esc(x) for x in row) + " |")
open(os.path.join(H, "tables.md"), "w").write("\n".join(out) + "\n")
print(counts, sum(counts.values()))
