#!/usr/bin/env python3
"""P2-26 -- the deck lines that name the twelve retiring keys, and (stage 2) their removal.

  usage: decks.py --dry      list every setter (DECK and ARGV), DEFAULT/LIVE, EDIT/RECORD; change nothing
         decks.py --apply    stage 2: remove the C2/C3 lines from the EDIT decks and stamp each LIVE deck's
                             header; writes DECK_EDITS.txt beside this script

The deck walk is P2-18's (`studies/p2_18_prereq/a2_apply.py` deck_files(): unit/, validation/, extensions/
input* and *.inputs, plus the gate deck studies/s1_surface_resolution/input; output dirs, retired_units,
*/pre/*, src_pre, legs, inputs_pre skipped). P2-18 already removed every line at a default there, so what
remains is LIVE by construction -- this script says so per line rather than assuming it.
RECORD decks (other studies/, retired_units/, */pre/*) are listed and never edited: a study deck is a record.
ARGV = live test scripts that pass a key on the command line (study run.py scripts are records).
`unit/sp_melt_skeleton/input` is not edited: the whole test moves to retired_units/ (Goal 5).

DEFAULT/LIVE: a line is DEFAULT when its value equals the deleted parse's default (numeric-aware), or when the
key is read only under another key that the same deck leaves off (melt.lambda_crit without spall.melt_aware = 1,
weibull.V_ref without size_volume = flux, jet_T_rec_fixed without rec_fixed, jet_T_above without jet_above_h);
otherwise LIVE.
"""
import importlib.util
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
T = os.path.join(ROOT, "tests", "MMWSpalling")
_spec = importlib.util.spec_from_file_location(
    "a2_apply", os.path.join(T, "studies", "p2_18_prereq", "a2_apply.py"))
a2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a2)

# key -> (stage/census item, the deleted parse's default, unit)
KEYS = {
    "spall.surface_normal":          ("C2", "0", "-"),
    "spall.flake_coherence_length":  ("C2", "0", "m"),
    "spall.melt_aware":              ("C3", "0", "-"),
    "melt.lambda_crit":              ("C3", "0.4", "-"),
    "spall.flake_length":            ("B3", "twocell", "-"),
    "spall.flake_fraction":          ("B3", "cf", "-"),
    "weibull.size_volume":           ("B3", "cell", "-"),
    "weibull.V_ref":                 ("B3", "-1.0", "m^3"),
    "surface_patch.pad_shield_open": ("B3", "-1.0", "-"),
    "surface_patch.jet_above_h":     ("B3", "0", "-"),
    "surface_patch.jet_T_above":     ("B3", "-1.0", "K"),
    "surface_patch.jet_T_rec_fixed": ("B3", "0.0", "K"),
}
# a value of another key that is itself the switch (jet_T_ent_mode keeps fixed / exhaust)
VALUE_SWITCH = {("surface_patch.jet_T_ent_mode", "rec_fixed"): "B3"}
# a key read only under its parent: DEFAULT unless the deck turns the parent on
PARENT = {"melt.lambda_crit": ("spall.melt_aware", "1"),
          "weibull.V_ref": ("weibull.size_volume", "flux"),
          "surface_patch.jet_T_rec_fixed": ("surface_patch.jet_T_ent_mode", "rec_fixed"),
          "surface_patch.jet_T_above": ("surface_patch.jet_above_h", "1")}
REMOVE_STAGE2 = ("spall.surface_normal", "spall.flake_coherence_length", "spall.melt_aware", "melt.lambda_crit")
NOT_EDITED = (os.path.join("unit", "sp_melt_skeleton", "input"),)
ARGV_SCRIPTS = ["unit/face_heat/test", "unit/removal_increment/test", "unit/jet_d2c/test",
                "unit/lateral_flux_lock/test", "unit/sp_melt_skeleton/test"]
RECORD_SKIP = ("/output", "/legs/", "_wit/", "src_pre", "inputs_pre", "__pycache__", ".old")


def kv(path):
    out = {}
    for ln in open(path, errors="replace").read().split("\n"):
        m = a2.LINE.match(ln)
        if m:
            out[m.group(1)] = m.group(2).strip()
    return out


def classify(k, v, deck_kv):
    if k in PARENT:
        pk, pv = PARENT[k]
        if not a2.same(deck_kv.get(pk, ""), pv):
            return "DEFAULT"
    return "DEFAULT" if a2.same(v, KEYS[k][1]) else "LIVE"


def record_decks(edit):
    """study / retired / pre decks that name a key: records, never edited"""
    out = []
    for dp, _d, fs in os.walk(T):
        if any(s in dp + "/" for s in RECORD_SKIP):
            continue
        for fn in fs:
            p = os.path.join(dp, fn)
            if (fn.startswith("input") or fn.endswith(".inputs")) and p not in edit:
                out.append(p)
    return sorted(out)


def setters():
    """[(kind, relpath, lineno, key, value, DEFAULT/LIVE, item)]"""
    edit = set(a2.deck_files())
    rows = []
    for kind, files in (("EDIT", sorted(edit)), ("RECORD", record_decks(edit))):
        for p in files:
            dk = kv(p)
            for i, ln in enumerate(open(p, errors="replace").read().split("\n"), 1):
                m = a2.LINE.match(ln)
                if not m:
                    continue
                k, v = m.group(1), m.group(2).strip()
                rel = os.path.relpath(p, T)
                if k in KEYS:
                    kd = kind if not (kind == "EDIT" and rel in NOT_EDITED) else "MOVED"
                    rows.append((kd, rel, i, k, v, classify(k, v, dk), KEYS[k][0]))
                elif (k, v) in VALUE_SWITCH:
                    rows.append((kind, rel, i, k, v, "LIVE", VALUE_SWITCH[(k, v)]))
    for rel in ARGV_SCRIPTS:
        p = os.path.join(T, rel)
        if not os.path.exists(p):
            continue
        for i, ln in enumerate(open(p).read().split("\n"), 1):
            for k in list(KEYS) + ["surface_patch.jet_T_ent_mode"]:
                # argv form `key=value` and dict form `"key": "value"` (face_heat's REFUSE legs)
                for m in re.finditer(r'%s(?:=|": ")([^"\s,\]]*)' % re.escape(k), ln):
                    v = m.group(1)
                    if k == "surface_patch.jet_T_ent_mode":
                        if v == "rec_fixed":
                            rows.append(("ARGV", rel, i, k, v, "LIVE", "B3"))
                        continue
                    lv = "DEFAULT" if a2.same(v, KEYS[k][1]) else ("LIVE" if v else "LIVE?")
                    rows.append(("ARGV", rel, i, k, v, lv, KEYS[k][0]))
            if "RECF" in ln and "rec_fixed" not in ln and re.search(r'"(KT|K15|A1[789])"', ln):
                rows.append(("ARGV", rel, i, "surface_patch.jet_T_ent_mode", "rec_fixed (RECF)", "LIVE", "B3"))
    return rows


STAMP = ("# P2-26 (2026-11-07): %s deleted from this deck -- the key is retired (census %s). "
         "The deleted line was LIVE, so this deck's next run differs from its earlier outputs BY DESIGN.\n")


def apply(dry):
    rows = [r for r in setters() if r[0] == "EDIT" and r[3] in REMOVE_STAGE2]
    by_file = {}
    for r in rows:
        by_file.setdefault(r[1], []).append(r)
    out = []
    for rel in sorted(by_file):
        p = os.path.join(T, rel)
        L = open(p).read().split("\n")
        gone = {r[2] for r in by_file[rel]}
        live = [r for r in by_file[rel] if r[5] == "LIVE" and r[6] == "C2"]
        for r in by_file[rel]:
            out.append("%-66s %4d  %-7s %-3s %s = %s" % (rel, r[2], r[5], r[6], r[3], r[4]))
        if not dry:
            body = "\n".join(l for j, l in enumerate(L, 1) if j not in gone)
            if live:
                what = ", ".join("`%s = %s`" % (r[3], r[4]) for r in live)
                body = STAMP % (what, "C2") + body
            open(p, "w").write(body)
    txt = "\n".join(out)
    print(txt)
    print("deck lines %s: %d in %d files" % ("that would go" if dry else "removed", len(out), len(by_file)))
    if not dry:
        open(os.path.join(HERE, "DECK_EDITS.txt"), "w").write(
            "P2-26 stage 2 deck edits -- every line removed (file, line, DEFAULT/LIVE, item, text); "
            "LIVE C2 decks carry a header stamp\n" + txt + "\n")


if __name__ == "__main__":
    if "--apply" in sys.argv or "--dry-apply" in sys.argv:
        apply("--apply" not in sys.argv)
    else:
        for r in setters():
            print("%-6s %-66s %4d  %-30s %-12s %-7s %s" % r)
