"""P2-30: THE ONE READER OF src/ for the tests that assert on the source.

class Integrator::MMWSpalling is one header cut into consecutive slices, each #included in place inside the class body
(src/Integrator/MMWSpalling/*.H). A test that opens MMWSpalling.H reads only the short parent; every assertion on it
would then pass or fail on the wrong text. So the tests read through here:

  tree()        the PER-FILE MAP {path relative to src/: text} of every C++ source under src/
  text(src)     the class's source as a SINGLE TEXT: the parent with every in-place include expanded, in order
  header(src)   what MMWSpalling.H held before the split: the slices expanded, the older fragments Shared.H and
                FaceRemoval.H left as their #include lines (the locks read those two under their own names)
  read(rel)     one file's text; exists(rel)
  slices(src)   the fragment names the parent includes, in order
  prune(src)    a map whose parent was substituted (a --pre pre-image): slices that parent no longer includes drop out

`src` is always optional and is a per-file map, so a harness can substitute ONE file and re-expand. Nothing else in
the live tests opens a path under src/ (unit/retired_keys and z0's L5 grep hold that).
"""
import os
import re

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
SRC = os.path.join(ROOT, "src")
PARENT = "Integrator/MMWSpalling.H"
SLICE_PREFIX = "Integrator/MMWSpalling/"
# read under their own names; Removal.H was deleted in P2-28 and survives only in old --pre pre-images
FRAGMENTS = ("Shared.H", "FaceRemoval.H", "Removal.H")
EXT = (".H", ".h", ".cpp", ".cc")
_INC = re.compile(r'^[ \t]*#include "Integrator/MMWSpalling/([^"\n]+)"[ \t]*\n', re.M)


def tree():
    out = {}
    for dp, _d, fs in os.walk(SRC):
        for fn in fs:
            if fn.endswith(EXT):
                p = os.path.join(dp, fn)
                out[os.path.relpath(p, SRC)] = open(p, errors="replace").read()
    return out


def read(rel):
    return open(os.path.join(SRC, rel), errors="replace").read()


def exists(rel):
    return os.path.exists(os.path.join(SRC, rel))


def slices(src=None):
    s = src if src is not None else tree()
    return _INC.findall(s[PARENT])


def text(src=None, keep=()):
    s = src if src is not None else tree()

    def one(m):
        return m.group(0) if m.group(1) in keep else s[SLICE_PREFIX + m.group(1)]
    return _INC.sub(one, s[PARENT])


def header(src=None):
    return text(src, keep=FRAGMENTS)


def prune(src):
    used = set(slices(src))
    return {k: v for k, v in src.items() if not (k.startswith(SLICE_PREFIX) and k[len(SLICE_PREFIX):] not in used)}
