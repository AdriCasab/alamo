#!/usr/bin/env python3
"""P2-30: THE IDENTITY INSTRUMENT -- is the compiler's input the same program?

Preprocesses src/mmwspalling.cc with the build's OWN compile line (read from `make -n`, -c replaced by -E -P) on the
live tree and on a pre-image S/<name>/ (a copy of mmwspalling.cc and the class's headers), tokenises both, and compares.

The one thing a pure move changes is what the `INFO` macro (`__FILE__, __func__, __LINE__`, src/Util/Util.H:21) expands
to: the file string and the line integer. The mask replaces ONLY those two tokens of an expansion, and only where the
triple `"<string>" , __func__ , <integer>` occurs; the message text after it is untouched.

Reported: the tokens compared; the UNMASKED difference count and whether every differing token is an INFO file string
or line integer (L3: "the unmasked diff lists ONLY such tokens"); the masked difference count (must be 0).

Usage:  pp_identity.py <pre_name>            live tree vs S/<pre_name>
        pp_identity.py --selftest <pre_name> the instrument itself: a tree vs itself is identical; a planted one-token
                                             edit, and a planted edit inside an INFO message, are both caught; L2
Exit 0 identical, 1 not."""
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
ENV = dict(os.environ, EIGEN=os.path.join(ROOT, "ext"), CPLUS_INCLUDE_PATH="/opt/homebrew/include",
           LIBRARY_PATH="/opt/homebrew/lib:/opt/homebrew/Cellar/gcc/15.2.0_1/lib/gcc/current")
OBJ = "obj/obj-3d-g++/mmwspalling.cc.o"
# one C++ token: a string or char literal, an identifier, a number (with any exponent letters), or one other character
TOK = re.compile(r'"(?:\\.|[^"\\\n])*"|\'(?:\\.|[^\'\\\n])*\'|[A-Za-z_]\w*|\d[\w.]*|\S')


def compile_line():
    """the build's own line for the main object, from `make -n` with src/mmwspalling.cc marked new"""
    out = subprocess.run(["make", "-n", "-W", "src/mmwspalling.cc", OBJ], cwd=ROOT, env=ENV, capture_output=True,
                         text=True).stdout
    lines = [ln for ln in out.splitlines() if " -c src/mmwspalling.cc -o %s" % OBJ in ln]
    if len(lines) != 1:
        sys.exit("pp_identity: expected one compile line from make -n, found %d" % len(lines))
    return lines[0]


def preprocess(line, source, extra_inc, out):
    """-E -P of `source`, the build's flags; extra_inc goes FIRST on the include path so a pre-image shadows src/"""
    cmd = line.replace(" -c src/mmwspalling.cc -o %s" % OBJ, " -E -P %s -o %s%s" % (
        source, out, "".join(" -I" + d for d in extra_inc)))
    # -I order: the pre-image must beat -I./src/, so move ours in front of it
    for d in extra_inc:
        cmd = cmd.replace(" -I" + d, "", 1).replace(" -I./src/", " -I%s -I./src/" % d, 1)
    r = subprocess.run(["bash", "-c", cmd], cwd=ROOT, env=ENV, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit("pp_identity: the preprocessor failed:\n" + r.stderr[-3000:])
    return cmd


def tokens(path):
    return TOK.findall(open(path, errors="replace").read())


def info_positions(t):
    """indices of the file string and of the line integer of every `"..." , __func__ , N` triple"""
    pos = set()
    for i in range(len(t) - 4):
        if t[i + 2] == "__func__" and t[i + 1] == "," and t[i + 3] == "," and t[i][0] == '"' and t[i + 4].isdigit():
            pos.add(i)
            pos.add(i + 4)
    return pos


def compare(a, b):
    """(identical_after_mask, n_tokens, n_unmasked_diffs, all_diffs_are_info, first_bad)"""
    if len(a) != len(b):
        k = next((i for i, (x, y) in enumerate(zip(a, b)) if x != y), min(len(a), len(b)))
        return False, (len(a), len(b)), None, False, (k, a[max(0, k - 8):k + 8], b[max(0, k - 8):k + 8])
    pa, pb = info_positions(a), info_positions(b)
    diffs = [i for i, (x, y) in enumerate(zip(a, b)) if x != y]
    bad = [i for i in diffs if not (i in pa and i in pb)]
    first = (bad[0], a[max(0, bad[0] - 8):bad[0] + 8], b[max(0, bad[0] - 8):bad[0] + 8]) if bad else None
    return (not bad), len(a), len(diffs), (not bad), first


def run_pp(line, pre, tmp, tag):
    """the token list of the live tree (pre None) or of a pre-image directory"""
    out = os.path.join(tmp, tag + ".ii")
    if pre is None:
        cmd = preprocess(line, "src/mmwspalling.cc", [], out)
    else:
        # the including file lives in the pre-image, so a quoted include resolves there first
        cmd = preprocess(line, os.path.join(pre, "mmwspalling.cc"), [pre], out)
    return tokens(out), cmd


def report(res):
    same, n, nd, only_info, first = res
    if isinstance(n, tuple):
        print("  token counts DIFFER: %d vs %d; first difference at %d:\n    %s\n    %s" % (n[0], n[1], first[0], first[1], first[2]))
        return
    print("  %d tokens compared; unmasked differences %d (every one an INFO file string or line integer: %s); "
          "after the mask: %s" % (n, nd, "yes" if only_info else "NO", "IDENTICAL" if same else "DIFFERENT"))
    if first:
        print("  first non-INFO difference at token %d:\n    %s\n    %s" % first)


def copy_pre(pre, dst):
    shutil.copytree(pre, dst)


def selftest(name):
    pre = os.path.join(HERE, name)
    line = compile_line()
    tmp = tempfile.mkdtemp(prefix="ppid_")
    a, cmd = run_pp(line, pre, tmp, "a")
    print("the build's compile line, -c replaced by -E -P:\n  " + cmd)
    b, _ = run_pp(line, pre, tmp, "b")
    r = compare(a, b)
    print("self-identity (the pre-image against itself):")
    report(r)
    ok = r[0] and r[2] == 0
    target = os.path.join("Integrator", "MMWSpalling.H")
    tree = os.path.join(tmp, "plant")
    copy_pre(pre, tree)
    src = open(os.path.join(tree, target)).read()
    planted = 'const Set::Scalar inv_dz2 = 1.0 / (DX[2] * DX[2]);'
    if src.count(planted) != 1:
        sys.exit("planted target not unique")
    open(os.path.join(tree, target), "w").write(src.replace(planted, planted.replace(" / ", " * "), 1))
    c, _ = run_pp(line, tree, tmp, "c")
    r = compare(a, c)
    print("planted one-token edit (`/` -> `*` in the kernel's inv_dz2) must be caught:")
    report(r)
    ok = ok and not r[0]
    tree2 = os.path.join(tmp, "plant2")
    copy_pre(pre, tree2)
    msg = 'was retired in P2-10'
    src2 = open(os.path.join(tree2, target)).read()
    if src2.count(msg) < 1:
        sys.exit("planted message not found")
    open(os.path.join(tree2, target), "w").write(src2.replace(msg, 'was retired in P2-11', 1))
    d, _ = run_pp(line, tree2, tmp, "d")
    r = compare(a, d)
    print("planted edit inside an INFO abort's MESSAGE must be caught (the mask replaces only the file and line):")
    report(r)
    ok = ok and not r[0]
    n_info = len(info_positions(a)) // 2
    print("INFO expansions in the translation unit: %d" % n_info)
    shutil.rmtree(tmp, ignore_errors=True)
    print("SELFTEST", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main():
    if sys.argv[1] == "--selftest":
        return selftest(sys.argv[2])
    pre = os.path.join(HERE, sys.argv[1])
    line = compile_line()
    tmp = tempfile.mkdtemp(prefix="ppid_")
    a, _ = run_pp(line, pre, tmp, "pre")
    b, cmd = run_pp(line, None, tmp, "live")
    print("pp_identity: live tree vs %s" % sys.argv[1])
    print("  command: " + cmd[:300] + " ...")
    r = compare(a, b)
    report(r)
    shutil.rmtree(tmp, ignore_errors=True)
    print("VERDICT", "IDENTICAL" if r[0] else "DIFFERENT")
    return 0 if r[0] else 1


if __name__ == "__main__":
    sys.exit(main())
