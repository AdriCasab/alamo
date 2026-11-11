#!/usr/bin/env python3
"""P2-30 stage 1: the six source-reading tests move onto tests/MMWSpalling/unit/_src.py. Plumbing only; no clause changes.
Every replacement must match exactly once (a `want` of 2 where the text is repeated) or the script stops before writing."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
U = os.path.join(ROOT, "tests", "MMWSpalling", "unit")
IMP = 'sys.path.insert(0, os.path.join(HERE, ".."))   # tests/MMWSpalling/unit/_src.py, the one reader of src/ (P2-30)\nimport _src  # noqa: E402\n'

EDITS = {
    "retired_keys": [
        ('SRC = os.path.join(ROOT, "src")\nKEY = "spallation.sp.constraint_depth"',
         IMP + 'KEY = "spallation.sp.constraint_depth"'),
        ('''SUBST = {"MMWBeam.H": os.path.join("Numeric", "MMWBeam.H"),
         "MMWSpalling.H": os.path.join("Integrator", "MMWSpalling.H"),
         # P2-28: Removal.H is deleted; under --pre a pre-image that still has it is read back in (R11 must go RED)
         "Removal.H": os.path.join("Integrator", "MMWSpalling", "Removal.H"),
         "FaceRemoval.H": os.path.join("Integrator", "MMWSpalling", "FaceRemoval.H")}
''', '''SUBST = {"MMWBeam.H": "Numeric/MMWBeam.H",
         "MMWSpalling.H": _src.PARENT,
         # P2-28: Removal.H is deleted; under --pre a pre-image that still has it is read back in (R11 must go RED)
         "Removal.H": _src.SLICE_PREFIX + "Removal.H",
         "FaceRemoval.H": _src.SLICE_PREFIX + "FaceRemoval.H"}
'''),
        ('''    """{relpath under src/: text} for every C++ source; under --pre the three headers come from <pre>"""
    out = {}
    for dp, _d, fs in os.walk(SRC):
        for fn in fs:
            if fn.endswith((".H", ".h", ".cpp", ".cc")):
                p = os.path.join(dp, fn)
                out[os.path.relpath(p, SRC)] = open(p, errors="replace").read()
    if pre:
        for fn, rel in SUBST.items():
            if os.path.exists(os.path.join(pre, fn)):
                out[rel] = open(os.path.join(pre, fn)).read()
    return out
''', '''    """{relpath under src/: text} for every C++ source (_src's per-file map); under --pre the headers come from
    <pre>, and slices its parent no longer includes drop out"""
    out = _src.tree()
    if pre:
        for fn, rel in SUBST.items():
            if os.path.exists(os.path.join(pre, fn)):
                out[rel] = open(os.path.join(pre, fn)).read()
        out = _src.prune(out)
    return out
'''),
        ('    m = src[MMW]\n    total = sum(t.count(KEY) for t in src.values())',
         '    m = _src.header(src)              # P2-30: the header\'s text, slices expanded\n    total = sum(t.count(KEY) for t in src.values())'),
        ('''def _mmw(fn):
    return lambda s: dict(s, **{MMW: fn(s[MMW])})
''', '''def _mmw(fn):
    """a mutation is an edit of ONE file: apply `fn` to the class file (the parent, then each slice) holding its target"""
    def go(s):
        for nm in [MMW] + [_src.SLICE_PREFIX + n for n in _src.slices(s) if n not in _src.FRAGMENTS]:
            t = fn(s[nm])
            if t != s[nm]:
                return dict(s, **{nm: t})
        return s
    return go
'''),
    ],
    "lateral_flux_lock": [
        ('''MMW = os.path.join(ROOT, "src", "Integrator", "MMWSpalling.H")
REM = os.path.join(ROOT, "src", "Integrator", "MMWSpalling", "Removal.H")   # P2-28: deleted (read as empty when absent)
''', IMP + '''REM = _src.SLICE_PREFIX + "Removal.H"   # P2-28: deleted (read as empty when absent)
'''),
        ('SHR = os.path.join(ROOT, "src", "Integrator", "MMWSpalling", "Shared.H")\n',
         'SHR = _src.SLICE_PREFIX + "Shared.H"\nFRM = _src.SLICE_PREFIX + "FaceRemoval.H"\n'),
        ('''    mmw, rem = strip_comments(open(MMW).read()), \\
        (strip_comments(open(REM).read()) if os.path.exists(REM) else "") + "\\n" + strip_comments(open(SHR).read())
''', '''    mmw, rem = strip_comments(_src.header()), \\
        (strip_comments(_src.read(REM)) if _src.exists(REM) else "") + "\\n" + strip_comments(_src.read(SHR))
'''),
        ('    frm28 = strip_comments(open(os.path.join(ROOT, "src", "Integrator", "MMWSpalling", "FaceRemoval.H")).read())\n',
         '    frm28 = strip_comments(_src.read(FRM))\n'),
        ('    frm = open(os.path.join(ROOT, "src", "Integrator", "MMWSpalling", "FaceRemoval.H")).read()\n',
         '    frm = _src.read(FRM)\n'),
        ('''    src_all = ""
    for dp, _, fs in os.walk(os.path.join(ROOT, "src")):
        for f in fs:
            if f.endswith((".H", ".cpp", ".cc")):
                src_all += open(os.path.join(dp, f), errors="replace").read()
''', '''    src_all = "".join(t for r, t in _src.tree().items() if r.endswith((".H", ".cpp", ".cc")))
'''),
        ('need(not left and not os.path.exists(REM),', 'need(not left and not _src.exists(REM),'),
        ('    hc = open(os.path.join(ROOT, "src", "Integrator", "HeatConduction.H")).read()\n',
         '    hc = _src.read("Integrator/HeatConduction.H")\n'),
        ('''    for dp, _d, fs in os.walk(os.path.join(ROOT, "src")):
        for fn in fs:
            if not fn.endswith((".H", ".h", ".cpp", ".cc")):
                continue
            p = os.path.join(dp, fn)
            in_list = False
            for ln, line in enumerate(open(p, errors="replace"), 1):
                if "P2-10: THE RETIRED KEYS" in line:
                    in_list = True
                if in_list and "remove it from the input" in line:
                    in_list = False
                    continue
                if in_list:
                    continue
                if DELETED.search(line):
                    left.append("%s:%d %s" % (os.path.relpath(p, ROOT), ln, line.strip()[:70]))
''', '''    for rel, txt in sorted(_src.tree().items()):
        in_list = False
        for ln, line in enumerate(txt.split("\\n"), 1):
            if "P2-10: THE RETIRED KEYS" in line:
                in_list = True
            if in_list and "remove it from the input" in line:
                in_list = False
                continue
            if in_list:
                continue
            if DELETED.search(line):
                left.append("src/%s:%d %s" % (rel, ln, line.strip()[:70]))
'''),
        ('''    nlist = sum(open(os.path.join(ROOT, "src", "Integrator", "MMWSpalling.H")).read()
                .count(k) for k in''', '''    nlist = sum(_src.header().count(k) for k in'''),
        ('''        for dp, _d, fs in os.walk(os.path.join(ROOT, "src")):
            for fn in fs:
                if fn.endswith((".H", ".h", ".cpp", ".cc")):
                    q = os.path.join(dp, fn)
                    for ln, line in enumerate(open(q, errors="replace"), 1):
                        if pat.search(line):
                            hits.append("%s:%d %s" % (os.path.relpath(q, ROOT), ln, line.strip()[:70]))
''', '''        for rel, txt in sorted(_src.tree().items()):
            for ln, line in enumerate(txt.split("\\n"), 1):
                if pat.search(line):
                    hits.append("src/%s:%d %s" % (rel, ln, line.strip()[:70]))
'''),
    ],
    "buckle_lock": [
        ('SRC = os.path.join(ROOT, "src", "Integrator")\n', IMP),
        ('''    mmw, frm, rem = read(SRC) if False else (open(os.path.join(SRC, "MMWSpalling.H")).read(),
                                            open(os.path.join(SRC, "MMWSpalling", "FaceRemoval.H")).read(),
                                            (open(os.path.join(SRC, "MMWSpalling", "Removal.H")).read()
                                             if os.path.exists(os.path.join(SRC, "MMWSpalling", "Removal.H"))
                                             else "") + "\\n"
                                            + open(os.path.join(SRC, "MMWSpalling", "Shared.H")).read())
''', '''    REM = _src.SLICE_PREFIX + "Removal.H"
    mmw, frm, rem = (_src.header(), _src.read(_src.SLICE_PREFIX + "FaceRemoval.H"),
                     (_src.read(REM) if _src.exists(REM) else "") + "\\n" + _src.read(_src.SLICE_PREFIX + "Shared.H"))
'''),
    ],
    "profile_depth": [
        ('''SHARED = os.path.join(ROOT, "src", "Integrator", "MMWSpalling", "Shared.H")
FR = os.path.join(ROOT, "src", "Integrator", "MMWSpalling", "FaceRemoval.H")
''', IMP + '''SHARED = _src.SLICE_PREFIX + "Shared.H"
FR = _src.SLICE_PREFIX + "FaceRemoval.H"
'''),
        ('    src = open(SHARED).read()\n', '    src = _src.read(SHARED)\n'),
        ('''    fr = re.sub(r"//[^\\n]*", "", open(FR).read())
    sh = re.sub(r"//[^\\n]*", "", open(SHARED).read())
''', '''    fr = re.sub(r"//[^\\n]*", "", _src.read(FR))
    sh = re.sub(r"//[^\\n]*", "", _src.read(SHARED))
'''),
    ],
    "face_heat": [
        ('''    for fn in ("MMWSpalling.H", os.path.join("MMWSpalling", "FaceRemoval.H")):
        code = re.sub(r"//[^\\n]*", "", open(os.path.join(ROOT, "src", "Integrator", fn)).read())
''', '''    for fn, txt in (("MMWSpalling.H", _src.header()), ("MMWSpalling/FaceRemoval.H", _src.read(_src.SLICE_PREFIX + "FaceRemoval.H"))):
        code = re.sub(r"//[^\\n]*", "", txt)
'''),
        ('''    code_fr = re.sub(r"//[^\\n]*", "", open(os.path.join(ROOT, "src", "Integrator", "MMWSpalling",
                                                         "FaceRemoval.H")).read())
''', '''    code_fr = re.sub(r"//[^\\n]*", "", _src.read(_src.SLICE_PREFIX + "FaceRemoval.H"))
'''),
    ],
    "face_removal": [
        ('''    rd = lambda *q: open(os.path.join(ROOT, "src", "Integrator", *q)).read()
    code = lambda t: "\\n".join(l.split("//")[0] for l in t.split("\\n"))
    frm, mmw, shr = code(rd("MMWSpalling", "FaceRemoval.H")), code(rd("MMWSpalling.H")), \\
        code(rd("MMWSpalling", "Shared.H"))
''', '''    code = lambda t: "\\n".join(l.split("//")[0] for l in t.split("\\n"))
    frm, mmw, shr = code(_src.read(_src.SLICE_PREFIX + "FaceRemoval.H")), code(_src.header()), \\
        code(_src.read(_src.SLICE_PREFIX + "Shared.H"))
'''),
    ],
}
# the import line goes after the ROOT definition in the tests that had no src-path line to replace there
AFTER_ROOT = {"face_heat", "face_removal"}
ROOT_LINE = 'ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))\n'


def main():
    new = {}
    for name, edits in EDITS.items():
        p = os.path.join(U, name, "test")
        t = open(p).read()
        for old, rep in edits:
            n = t.count(old)
            if n != 1:
                sys.exit("%s: expected 1 match, found %d of:\n%s" % (name, n, old[:200]))
            t = t.replace(old, rep, 1)
        if name in AFTER_ROOT:
            if t.count(ROOT_LINE) != 1:
                sys.exit("%s: ROOT line" % name)
            t = t.replace(ROOT_LINE, ROOT_LINE + IMP, 1)
        new[p] = t
    for p, t in new.items():
        open(p, "w").write(t)
        print("edited", os.path.relpath(p, ROOT))


if __name__ == "__main__":
    main()
