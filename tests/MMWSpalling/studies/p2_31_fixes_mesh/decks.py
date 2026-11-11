#!/usr/bin/env python3
"""P2-31: the deck reader step zero uses (L4, L5, L13, L16). An ALAMO input is `key = value` lines, `#` comments and
`FILE = <path>` includes (paths from the repo root); a later definition of a key wins, and command-line `key=value`
tokens come after the file. `resolve(deck, argv_overrides)` returns {key: (value, file, line)}.
`live_decks()` lists every input under unit/ and validation/ that a test or a study can run today (not output*/ copies,
not legs/, not retired_units/, not *.old*)."""
import os
import re
import shlex

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
T = os.path.join(ROOT, "tests", "MMWSpalling")
KV = re.compile(r"^\s*([A-Za-z_][\w.\[\]]*)\s*=\s*(.*?)\s*$")


def _read(path, out, depth=0):
    if depth > 8:
        raise RuntimeError("include depth > 8 at %s" % path)
    for n, raw in enumerate(open(path, errors="replace"), 1):
        ln = raw.split("#", 1)[0].rstrip()
        m = KV.match(ln)
        if not m:
            continue
        k, v = m.group(1), m.group(2)
        if k == "FILE":
            _read(os.path.join(ROOT, v), out, depth + 1)
            continue
        out[k] = (v, os.path.relpath(path, ROOT), n)


def resolve(deck, overrides=()):
    out = {}
    _read(deck, out)
    for t in overrides:
        if "=" in t:
            k, v = t.split("=", 1)
            out[k] = (v, "argv", 0)
    return out


def argv_deck(argv_sh):
    """(deck, overrides) from a recorded mpirun argv"""
    toks = shlex.split(open(argv_sh).read())
    i = next(i for i, t in enumerate(toks) if t.endswith("mmwspalling-3d-g++"))
    deck = toks[i + 1]
    if not os.path.isabs(deck):
        deck = os.path.join(ROOT, deck)
    return deck, [t for t in toks[i + 2:] if "=" in t]


def live_decks():
    out = []
    for top in ("unit", "validation"):
        for dp, dns, fns in os.walk(os.path.join(T, top)):
            dns[:] = [d for d in dns if not (d.startswith("output") or d == "legs" or ".old" in d
                                             or d == "__pycache__" or d.startswith("covers"))]
            for fn in fns:
                if fn.startswith("input") and not fn.endswith((".py", ".md", ".txt", ".csv")):
                    out.append(os.path.join(dp, fn))
    return sorted(out)


def val(R, k, default=None):
    return R[k][0] if k in R else default


def face_list_on(R):
    return val(R, "surface_patch.face_list", "0").strip() in ("1", "true")


def phase_laws(R):
    n = int(float(val(R, "microstructure.nphases", "0")))
    return [val(R, "microstructure.phase%d.law" % i, "constant") for i in range(n)]
