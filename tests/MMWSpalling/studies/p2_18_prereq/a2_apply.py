#!/usr/bin/env python3
"""P2-18 A2 -- (i) the "retired at default" treatment in MMWSpalling::Parse, (ii) the six DEAD
keys' parses deleted, (iii) the kept decks cleaned of lines that set a retiring key to its
default (the user, 2026-11-05: "remove from the input files inputs which are no longer in
use").

  usage: a2_apply.py --dry      list every deck line that would go, and the src edits, change nothing
         a2_apply.py --apply    do it (src + decks); writes A2_DECK_EDITS.txt beside this script

DECKS = unit/, validation/, extensions/ (input*, *.inputs) + studies/s1_surface_resolution/input
(a gate deck). NOT touched: output dirs, retired_units, */pre/*, src_pre, */legs/*, argv.sh
records, study run.py scripts (records).

The line goes when the key is DEAD (parse deleted today; any value) or when it is a RETIRING
key set to its parse default (numeric-aware). A retiring key at a LIVE value stays.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
T = os.path.join(ROOT, "tests", "MMWSpalling")
M = os.path.join(ROOT, "src", "Integrator", "MMWSpalling.H")
G = os.path.join(ROOT, "src", "Numeric", "Material", "Grimsel.H")

# DEAD keys: the parse goes today; any value on a deck was never read -> the line goes
DEAD = [r"microstructure\.phase\d+\.name", r"microstructure\.phase\d+\.alpha_attenuation",
        r"surface_patch\.face", r"surface_patch\.T_ext", r"surface_patch\.type",
        r"microstructure\.phase\d+\.grimsel\.T_ref"]
# RETIRING keys (census B-D; user's decisions 1, 2, 6, 7): the line goes only AT the parse default
AT_DEFAULT = {
    "damage.enabled": "0",                       # C1a (DAMAGE lineage retires)
    "conductivity.kappa_damage_alpha": "0.0",    # C1a
    "spall.sample": "top_cell",                  # C1a
    "solver": "explicit",                        # D2 (implicit solver deleted, decision 1)
    "spall.surface_normal": "0",                 # C2 (decision 7)
    "spall.flake_coherence_length": "0.0",       # C2
    "spall.melt_aware": "0",                     # C3 (decision 6)
    "bit.enabled": "0",                          # H  (decision 2)
    "weibull.size_volume": "cell",               # B3 (closed)
    "spall.flake_length": "twocell",             # B3 (closed)
    "spall.flake_fraction": "cf",                # B3 (closed)
}
SKIP_DIRS = ("/output", "retired_units", "/pre/", "src_pre", "/legs/", "__pycache__", "inputs_pre")


def same(a, b):
    try:
        return float(a) == float(b)
    except ValueError:
        return a == b


def deck_files():
    out = []
    for top in ("unit", "validation", "extensions"):
        for dp, _d, fs in os.walk(os.path.join(T, top)):
            if any(s in dp + "/" for s in SKIP_DIRS):
                continue
            for fn in fs:
                if fn.startswith("input") or fn.endswith(".inputs"):
                    out.append(os.path.join(dp, fn))
    out.append(os.path.join(T, "studies", "s1_surface_resolution", "input"))
    return sorted(out)


LINE = re.compile(r"^\s*([A-Za-z0-9_.]+)\s*=\s*(.*?)\s*(#.*)?$")


def plan_decks():
    """[(path, lineno, line, reason)]"""
    plan = []
    for p in deck_files():
        for i, ln in enumerate(open(p).read().split("\n"), 1):
            m = LINE.match(ln)
            if not m:
                continue
            k, v = m.group(1), m.group(2).strip()
            if any(re.fullmatch(d, k) for d in DEAD):
                plan.append((p, i, ln, "DEAD"))
            elif k in AT_DEFAULT and same(v, AT_DEFAULT[k]):
                plan.append((p, i, ln, "at default %s" % AT_DEFAULT[k]))
    return plan


# ------------------------------------------------------------------------ src ----
TOLERANCE = '''        // P2-18 A2: KEYS RETIRED AT A DEFAULT. A deleted parse is not a guard: a
        // top-level key (`spall.*`, `microstructure.*`, ...) is then silently
        // ignored, and a strict-scope key (`surface_patch.*`, `beam.*`, the
        // per-phase `grimsel.*`) aborts as "specified but not used" whatever its
        // value. Kept decks name retiring keys at their old defaults far more
        // often than at live values (census 2026-11-05 §4 A2), so each entry
        // here is read ONCE, which marks it used for the strict scopes, and
        // compared with the value the deleted code accepted: equal (or "*" = the
        // key was never read, any value is correct) -> the deck runs unchanged;
        // anything else -> abort naming the key, the value and the packet that
        // retired it. Deleting the parse of a key with a single old default ADDS
        // A ROW HERE; it does not go in the abort list above.
        struct RetiredAtDefault { const char* key; const char* value; const char* note; };
        static const RetiredAtDefault retired_at_default[] = {
            {"surface_patch.face",  "zhi",
             "P2-18: consumed, checked == zhi, read by nothing (only the +z patch exists)"},
            {"surface_patch.T_ext", "*",
             "P2-18: parsed into a member nothing read"},
            {"surface_patch.type",  "prescribed_temperature",
             "P2-18: the legacy Step 6d key; surface_patch.mode decides"},
        };
        for (const RetiredAtDefault& r : retired_at_default)
        {
            if (!pp.contains(r.key)) continue;
            std::vector<std::string> toks;
            static_cast<amrex::ParmParse&>(pp).queryarr(r.key, toks);   // marks the entry used
            std::string got;
            for (const std::string& t : toks) got += (got.empty() ? "" : " ") + t;
            if (std::string(r.value) == "*" || got == r.value) continue;
            Util::Abort(INFO, r.key, " = ", got, " was retired at its default (", r.value,
                        ") -- ", r.note, "; remove it from the input or restore the code");
        }
'''
ANCHOR = '''        for (const char* k : {"spallation.sp.constraint_depth"})
            if (pp.contains(k))
                Util::Abort(INFO, k, " was retired in P2-16 (the constraint depth was "
                            "deleted; spall.buckle_span carries the minimum flux); remove "
                            "it from the input");
'''

SRC_EDITS = [  # (file, old, new, what)
    (M, ANCHOR, ANCHOR + TOLERANCE, "the retired-at-default table + helper after the P2-16 list"),
    (M, '''                pp.query_default (pfx + ".name",     value.phases[i].name,
                                  std::string("phase") + std::to_string(i));
''', "", "delete the parse of microstructure.phase<i>.name (DEAD)"),
    (M, '''                pp.query_required(pfx + ".alpha_attenuation",
                                  value.phases[i].alpha_atten);
''', "", "delete the parse of microstructure.phase<i>.alpha_attenuation (DEAD, REQUIRED)"),
    (M, "        std::string name;\n", "", "delete the Phase member `name`"),
    (M, "        Set::Scalar alpha_atten = 0.0;\n", "", "delete the Phase member `alpha_atten`"),
    (M, '''            // Legacy `type` key — always consumed so the strict parser does
            // not choke on existing Step 6d-style inputs.
            std::string type;
            pp.query_default("type", type, std::string("prescribed_temperature"));
''', "", "delete the parse of surface_patch.type (DEAD; retired at default)"),
    (M, '''            if (v.mode == Mode::PrescribedTemperature &&
                type != "prescribed_temperature")
                Util::Abort(INFO, "surface_patch.type must be "
                                  "'prescribed_temperature', got: ", type);
''', "", "delete the type check"),
    (M, '''            std::string face;
            pp.query_default("face", face, std::string("zhi"));
            if (face != "zhi")
                Util::Abort(INFO, "surface_patch.face must be 'zhi' (only +z "
                                  "patch implemented), got: ", face);
''', "", "delete the parse of surface_patch.face (DEAD; retired at default)"),
    (M, '''            pp.query_default("T_ext",   v.T_ext,   0.0);
''', "", "delete the parse of surface_patch.T_ext (DEAD; retired at default)"),
    (M, "        Set::Scalar T_ext   = 0.0;\n", "", "delete the SurfacePatch member `T_ext`"),
    (G, '''        pp.query_default("T_ref", v.T_ref, 293.15);
''', "", "delete the parse of microstructure.phase<i>.grimsel.T_ref (DEAD: overwritten by phase<i>.T_ref)"),
]


def apply_src(dry):
    for path, old, new, what in SRC_EDITS:
        s = open(path).read()
        n = s.count(old)
        print("  %-5s %s  [%d match]  %s" % ("src", os.path.relpath(path, ROOT), n, what))
        if n != 1:
            raise SystemExit("src edit did not match exactly once: " + what)
        if not dry:
            open(path, "w").write(s.replace(old, new, 1))


def apply_decks(dry):
    plan = plan_decks()
    by_file = {}
    for p, i, ln, why in plan:
        by_file.setdefault(p, []).append((i, ln, why))
    rows = []
    for p in sorted(by_file):
        rel = os.path.relpath(p, T)
        for i, ln, why in by_file[p]:
            rows.append("%-70s %4d  %-28s %s" % (rel, i, why, ln.strip()))
        if not dry:
            L = open(p).read().split("\n")
            gone = {i for i, _, _ in by_file[p]}
            open(p, "w").write("\n".join(l for j, l in enumerate(L, 1) if j not in gone))
    txt = "\n".join(rows)
    print(txt)
    print("deck lines %s: %d in %d files" % ("that would go" if dry else "removed", len(rows), len(by_file)))
    if not dry:
        open(os.path.join(HERE, "A2_DECK_EDITS.txt"), "w").write(
            "P2-18 A2 deck cleanup -- every line removed (file, line, reason, text)\n" + txt + "\n")
    return len(rows), len(by_file)


if __name__ == "__main__":
    dry = "--apply" not in sys.argv
    print("P2-18 A2 %s" % ("DRY RUN" if dry else "APPLY"))
    print("== src")
    apply_src(dry)
    print("== decks")
    apply_decks(dry)
