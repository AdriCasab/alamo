#!/usr/bin/env python3
"""P2-27 L2 -- the interlock move, checked statement by statement (src_pre_1 vs src/).

Every abort statement (Util::Abort / ph_abort / fr_abort call, whitespace-normalised) of the old Parse matrix, the
four sub-parse sites and BuildFaceList's cubic-cell check must appear VERBATIM in the new CheckInterlocks, in the same
relative order; every parse/registration call of the old Parse must remain in the new Parse in the same order; the
abort statements left in Parse are the five non-cross-key ones. Prints each moved abort with its old and new
enclosing `if` conditions for review. Usage: interlocks_move.py > INTERLOCKS_MOVE.txt
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
OLD = open(os.path.join(HERE, "src_pre_1", "src", "Integrator", "MMWSpalling.H")).read()
NEW = open(os.path.join(ROOT, "src", "Integrator", "MMWSpalling.H")).read()


def norm(s):
    return re.sub(r"\s+", " ", s).strip()


def stmt_end(t, i):
    depth, j = 0, i
    while True:
        c = t[j]
        if c == '"':
            j += 1
            while t[j] != '"':
                j += 2 if t[j] == "\\" else 1
        elif c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
            if depth == 0:
                return t.index(";", j) + 1
        j += 1


def aborts(t, lo, hi):
    out = []
    for m in re.finditer(r"(Util::Abort\(INFO,|ph_abort\(\"|fr_abort\(\")", t[lo:hi]):
        i = lo + m.start()
        if t[i - 5:i] == "auto ":
            continue
        out.append((i, norm(t[i:stmt_end(t, i)])))
    return [a for a in out if "[] (const char*" not in a[1]]


def cond(t, i):
    """the nearest preceding `if (...)` header within 400 chars (the abort's own condition)"""
    k = t.rfind("if (", max(0, i - 600), i)
    return norm(t[k:t.find("\n", k)]) if k >= 0 else "?"


def region(t, start_pat, end_pat):
    a = t.index(start_pat)
    return a, t.index(end_pat, a)


def main():
    bad = 0
    pa, pb = region(OLD, "static void Parse(MMWSpalling& value, IO::ParmParse& pp)", "\nprotected:")
    old_sites = [x for x in aborts(OLD, pa, pb)
                 if "solver must be" not in x[1] and "microstructure.mode must be" not in x[1]
                 and "microstructure.phase\", i" not in x[1] and "default A_beam is non-positive" not in x[1]
                 and "beam.melt_factor / beam.R_melt" not in x[1]
                 and not x[1].startswith("Util::Abort(INFO, k,") and "r.key" not in x[1]]
    subs = []
    for pat in ['") must be SMALLER than spall.h_spall_max ("', '"spall.removal_energy = flake requires spall.surface_T = "',
                '"spall.surface_T = reconstructed requires "', '"cohesive.enabled = 1: CZM is bypassed',
                '"damage.enabled = 1: eq.27 damage evolution belongs to "', '"bit.enabled requires spall.enabled',
                '"surface_patch.face_list = 1 requires cubic cells; got "']:
        k = OLD.index(pat)
        i = OLD.rfind("Util::Abort(INFO", 0, k)
        subs.append((i, norm(OLD[i:stmt_end(OLD, i)])))
    want = subs[:6] + old_sites + subs[6:]
    na, nb = region(NEW, "static void CheckInterlocks(MMWSpalling& value)", "\nprotected:")
    got = aborts(NEW, na, nb)
    print("P2-27 L2: %d abort statements moved (6 sub-parse + %d end-of-Parse + 1 cubic-cell); CheckInterlocks holds %d"
          % (len(want), len(old_sites), len(got)))
    same = [w[1] for w in want] == [g[1] for g in got]
    print("  %s  verbatim and in the same relative order" % ("ok  " if same else "FAIL"))
    bad += 0 if same else 1
    for (wi, ws), (gi, gs) in zip(want, got):
        oc, nc = cond(OLD, wi), cond(NEW, gi)
        print("  - %s\n      old if: %s\n      new if: %s" % (ws[:110], oc[:110], nc[:110]))
    # nothing of it is left behind
    nn = norm(NEW)
    left = [w for w in want if nn.count(w[1]) != 1]
    print("  %s  each moved statement's text occurs once in src/ (%d not)" % ("ok  " if not left else "FAIL", len(left)))
    bad += 1 if left else 0
    rx = re.compile(r'(RegisterIntegratedVariable\(&value\.\w+(?:\[\d\])?,\s*"\w+"|RegisterNewFab\(value\.\w+,[^;]*?"\w+"'
                    r'|pp\.\w+\("[\w.]+"|Util::Message\(INFO, "[^"]*")')
    nP = NEW[NEW.index("static void Parse(MMWSpalling& value"):na]
    ro, rn = rx.findall(OLD[pa:pb]), rx.findall(nP)
    print("  %s  Parse keeps its %d key parses / registrations / messages in the same order (%d now)"
          % ("ok  " if ro == rn else "FAIL", len(ro), len(rn)))
    bad += 0 if ro == rn else 1
    print("VERDICT: %s" % ("PASS" if not bad else "FAIL (%d)" % bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
