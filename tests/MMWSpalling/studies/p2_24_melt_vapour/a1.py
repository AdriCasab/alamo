#!/usr/bin/env python3
"""P2-24 a1 -- the face gate, scored as the packet allows (L9): face_p2_23 vs face_p2_24.

  the 13 FLAME legs: every common file byte-identical, EXCEPT thermo.dat, whose common columns (the wall-clock timers
     excluded) must be identical and whose only new columns are fr_n_molten and fr_n_lid, 0 on every row
  BEAM_S (the face list without a flame on the #19 deck, 1 s): the deck now carries the melt keys and vapor.enabled = 1,
     so its thermo gains the beam's melt columns and the vapour books, and its plotfiles gain the vapour fields; it is
     scored on what both builds write: the common thermo columns identical, the new count columns 0, every plotfile's
     COMMON fields identical in double precision (yt), the CSVs byte-identical
Paths only on one side: none allowed outside BEAM_S/. Usage: a1.py  (prints the verdict line "a1 PASS|FAIL")
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
W = os.path.abspath(os.path.join(HERE, "..", "tree_sweep_0922", "witness"))
A, B = "p2_23", "p2_24"
TIMERS = {"fr_t_removal", "fr_t_step", "fl_t_path"}
NEW_FLAME = {"fr_n_molten", "fr_n_lid"}
ZERO_BEAM = {"fr_n_molten", "fr_n_lid", "fr_n_vapor", "fr_n_vapor_int", "ledger_E_vapor"}


def md5s(label):
    out = {}
    for ln in open(os.path.join(W, "face_%s.md5" % label)):
        h, rel = ln.rstrip("\n").split("  ", 1)
        out[rel] = h
    return out


def cols(path):
    rows = [ln.rstrip("\n").split("\t") for ln in open(path) if ln.strip()]
    hdr = [h.strip() for h in rows[0]]
    c = {h: [] for h in hdr}
    for r in rows[1:]:
        for h, v in zip(hdr, r):
            c[h].append(v.strip())
    return hdr, c


def fields_equal(pa, pb):
    import numpy as np
    import yt
    yt.set_log_level(50)
    da, db = yt.load(pa), yt.load(pb)
    fa = {f for t, f in da.field_list}
    fb = {f for t, f in db.field_list}
    common = sorted(fa & fb)
    a, b = da.all_data(), db.all_data()
    oa = np.lexsort((a["x"].to_value(), a["y"].to_value(), a["z"].to_value()))
    ob = np.lexsort((b["x"].to_value(), b["y"].to_value(), b["z"].to_value()))
    moved = [f for f in common if not np.array_equal(a[("boxlib", f)].to_value()[oa], b[("boxlib", f)].to_value()[ob],
                                                     equal_nan=True)]
    return common, sorted(fb - fa), sorted(fa - fb), moved


def main():
    a, b = md5s(A), md5s(B)
    only_a, only_b = sorted(set(a) - set(b)), sorted(set(b) - set(a))
    diff = sorted(r for r in set(a) & set(b) if a[r] != b[r])
    bad = []
    print("a1: face_%s %d files, face_%s %d; only in %s %d, only in %s %d; differing common paths %d" % (
        A, len(a), B, len(b), A, len(only_a), B, len(only_b), len(diff)))
    for r in only_a + only_b:
        if not r.startswith("BEAM_S/"):
            bad.append("only on one side: " + r)
    beam_pf = set()
    for r in diff:
        leg = r.split("/run/")[0]
        if r.endswith("thermo.dat"):
            ha, ca = cols(os.path.join(W, "face_" + A, r))
            hb, cb = cols(os.path.join(W, "face_" + B, r))
            moved = [h for h in ha if h not in TIMERS and h in cb and ca[h] != cb[h]]
            gone = [h for h in ha if h not in cb]
            new = [h for h in hb if h not in ca]
            if leg == "BEAM_S":
                nz = [h for h in new if h in ZERO_BEAM and any(float(v) != 0.0 for v in cb[h])]
                ok = not moved and not gone and not nz
                print("  %-14s thermo: common %d columns %s; new %d: %s; the counts and the vapour book all 0: %s" % (
                    leg, len([h for h in ha if h in cb]), "identical" if not moved else "MOVED " + " ".join(moved),
                    len(new), " ".join(new), not nz))
                for h in new:
                    if h not in ZERO_BEAM:
                        v = [float(x) for x in cb[h]]
                        fin = [x for x in v if math.isfinite(x)]
                        print("  %-14s   %-12s %s" % ("", h, "%d rows, %.6g .. %.6g%s" % (
                            len(v), min(fin), max(fin), "" if len(fin) == len(v) else " (%d NaN)" % (len(v) - len(fin)))
                            if fin else "%d rows, all NaN" % len(v)))
            else:
                nz = [h for h in new if any(float(v) != 0.0 for v in cb[h])]
                ok = not moved and not gone and set(new) <= NEW_FLAME and not nz
                print("  %-14s thermo: common columns %s; new %s (all 0: %s)%s" % (
                    leg, "identical" if not moved else "MOVED " + " ".join(moved), " ".join(new) or "none", not nz,
                    "" if set(new) <= NEW_FLAME else "  ⚠ unexpected new columns"))
            if not ok:
                bad.append("thermo " + r)
        elif leg == "BEAM_S" and "cell/" in r:
            beam_pf.add(r.split("cell/")[0] + "cell")
        else:
            bad.append("file differs: " + r)
    for pf in sorted(beam_pf):
        common, only_new, only_old, moved = fields_equal(os.path.join(W, "face_" + A, pf), os.path.join(W, "face_" + B, pf))
        print("  BEAM_S         %s: %d common fields %s; new fields %s; gone %s" % (
            pf.split("/run/")[1], len(common), "identical" if not moved else "MOVED " + " ".join(moved),
            " ".join(only_new) or "none", " ".join(only_old) or "none"))
        if moved or only_old:
            bad.append("BEAM_S plotfile " + pf)
    for r in only_a + only_b:
        if r.startswith("BEAM_S/") and not any(r.startswith(p) for p in beam_pf):
            print("  BEAM_S         only in %s: %s" % (A if r in a else B, r))
            bad.append("BEAM_S file on one side, outside a field-compared plotfile: " + r)
    for x in bad:
        print("  ⚠ " + x)
    print("a1 %s (%d failure(s))" % ("PASS" if not bad else "FAIL", len(bad)))
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
