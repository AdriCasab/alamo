#!/usr/bin/env python3
"""Comment review check: is the CODE (everything that is not a comment) of each file identical to the pre-image?
Strips // and /* */ comments outside string/char literals, then compares the remaining text with all whitespace
runs collapsed. Also lists the comment lines removed and added (as multisets)."""
import collections, os, re, sys
ROOT = "/Users/tzetze20/amr_tools/alamo"
PRE = os.path.join(ROOT, "tests/MMWSpalling/studies/p2_30_structure/src_pre_5")
PAIRS = [("src/Integrator/MMWSpalling.H", "Integrator/MMWSpalling.H")]
for f in sorted(os.listdir(os.path.join(ROOT, "src/Integrator/MMWSpalling"))):
    if f.endswith(".H"): PAIRS.append(("src/Integrator/MMWSpalling/" + f, "Integrator/MMWSpalling/" + f))
PAIRS.append(("src/Numeric/MMWBeam.H", "Numeric/MMWBeam.H"))
for f in sorted(os.listdir(os.path.join(ROOT, "src/Numeric/Material"))):
    if f.endswith(".H"): PAIRS.append(("src/Numeric/Material/" + f, "Numeric/Material/" + f))

def split(text):
    """-> (code text, list of comment texts)"""
    code, com, i, n = [], [], 0, len(text)
    while i < n:
        c = text[i]
        if c == '"' or c == "'":
            j = i + 1
            while j < n and text[j] != c:
                j += 2 if text[j] == "\\" else 1
            code.append(text[i:j + 1]); i = j + 1
        elif text.startswith("//", i):
            j = text.find("\n", i); j = n if j < 0 else j
            com.append(text[i:j].rstrip()); i = j
        elif text.startswith("/*", i):
            j = text.find("*/", i); j = n if j < 0 else j + 2
            com.append(text[i:j]); i = j
        else:
            code.append(c); i += 1
    return re.sub(r"\s+", " ", "".join(code)).strip(), com

bad = 0; tot_rm = tot_add = 0; removed_all = []
for live, pre in PAIRS:
    a = open(os.path.join(PRE, pre)).read(); b = open(os.path.join(ROOT, live)).read()
    ca, ma = split(a); cb, mb = split(b)
    same = ca == cb
    A, B = collections.Counter(ma), collections.Counter(mb)
    rm, add = A - B, B - A
    nrm, nadd = sum(rm.values()), sum(add.values())
    tot_rm += nrm; tot_add += nadd
    print("%-46s code %s   comments: %4d removed/changed-from, %4d added/changed-to   lines %d -> %d" % (
        live, "IDENTICAL" if same else "**DIFFERENT**", nrm, nadd, a.count("\n"), b.count("\n")))
    if not same:
        bad += 1
        k = next((i for i, (x, y) in enumerate(zip(ca, cb)) if x != y), min(len(ca), len(cb)))
        print("     first difference at code char %d:\n       pre : %r\n       live: %r" % (k, ca[max(0, k - 80):k + 80], cb[max(0, k - 80):k + 80]))
    if "--dump" in sys.argv and rm:
        removed_all.append((live, [l for l in ma if rm[l] > 0 and not (rm.subtract([l]) or False)]))
print("TOTAL comment texts out %d, in %d;  files with a code difference: %d" % (tot_rm, tot_add, bad))
if "--dump" in sys.argv:
    out = sys.argv[sys.argv.index("--dump") + 1]
    with open(out, "w") as f:
        for live, ls in removed_all:
            f.write("\n### %s\n\n```\n" % live); f.write("\n".join(ls)); f.write("\n```\n")
    print("removed comment texts ->", out)
sys.exit(1 if bad else 0)
