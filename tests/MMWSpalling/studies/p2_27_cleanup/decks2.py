#!/usr/bin/env python3
"""P2-27 stage 2 (C1a) -- the deck lines of the retired damage lineage, removed; writes DECK_EDITS_2.txt.

  usage: decks2.py --dry | --apply
The deck walk is z0's (P2-18's a2_apply deck_files + the shared rock file); decks on the retiring list move to
retired_units/ or retired_validation/ unchanged (records) and are not edited.
  * microstructure.phase<i>.A_D = <v>: DEFAULT when v == 0.0 (the new default) or when the deck leaves damage.enabled
    off (A_D was read only by the damage law's growth rate, so a non-zero value was inert); otherwise LIVE (none
    expected outside the retiring decks: z0 section 5)
  * unit/amr_microstructure_regrid/input: conductivity.kappa_damage_alpha (ln 2) and damage.ic.* -- LIVE: they made
    the deck's D field and kappa_eff non-trivial; the test's D / kappa_eff checks go with them (Goal 10)
"""
import importlib.util
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("z0", os.path.join(HERE, "z0.py"))
z0 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(z0)
T = z0.T
REGRID = os.path.join(T, "unit", "amr_microstructure_regrid", "input")


def plan():
    rows = []
    for p in z0.live_decks() + [os.path.join(T, "validation", "shared", "rock_grimsel.inputs")]:
        rel = os.path.relpath(p, T)
        if z0.retiring(rel):
            continue
        eff = z0.resolve_deck(p) if not p.endswith(".inputs") else {}
        dmg = z0.a2.same(eff.get("damage.enabled", "0"), "1")
        for i, ln in enumerate(open(p, errors="replace").read().split("\n"), 1):
            m = z0.LINE.match(ln)
            if not m:
                continue
            k, v = m.group(1), m.group(2).strip()
            if re.fullmatch(r"microstructure\.phase\d+\.A_D", k):
                tag = "DEFAULT" if (z0.a2.same(v, "0.0") or not dmg) else "LIVE"
                why = "== 0.0" if z0.a2.same(v, "0.0") else "inert: damage.enabled off"
                rows.append((p, i, ln, tag, why))
            elif p == REGRID and (k == "conductivity.kappa_damage_alpha" or k.startswith("damage.")):
                rows.append((p, i, ln, "LIVE", "the deck's damage IC / factor (its D and kappa_eff checks go too)"))
    return rows


def main():
    rows = plan()
    out = ["P2-27 stage 2 deck edits (decks2.py): %d line(s) in %d file(s)"
           % (len(rows), len(set(r[0] for r in rows)))]
    for p, i, ln, tag, why in rows:
        out.append("%-7s %s:%d  %s   [%s]" % (tag, os.path.relpath(p, T), i, ln.strip()[:90], why))
    live = [r for r in rows if r[3] == "LIVE" and r[0] != REGRID]
    out.append("LIVE outside amr_microstructure_regrid: %d (expected 0)" % len(live))
    txt = "\n".join(out)
    print(txt)
    if "--apply" in sys.argv:
        if live:
            print("REFUSED: a LIVE A_D line outside the retiring decks")
            return 1
        by = {}
        for p, i, ln, tag, why in rows:
            by.setdefault(p, set()).add(i)
        for p, drop in by.items():
            L = open(p).read().split("\n")
            open(p, "w").write("\n".join(l for n, l in enumerate(L, 1) if n not in drop))
        open(os.path.join(HERE, "DECK_EDITS_2.txt"), "w").write(txt + "\n")
        print("applied; DECK_EDITS_2.txt written")
    return 0


if __name__ == "__main__":
    sys.exit(main())
