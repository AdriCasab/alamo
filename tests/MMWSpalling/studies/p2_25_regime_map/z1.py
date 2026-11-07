#!/usr/bin/env python3
"""P2-25 z1 -- the face gate's comparator (DR-1). face_p2_24's BEAM_S thermo.dat compared with ITSELF, column by column:
under the old rule (`ca[h] != cb[h]` on float lists: NaN != NaN) and under face_gate.col_moved (NaN equal to NaN), the
wall-clock timers excluded as cmp excludes them, the file parsed TWICE as cmp parses its two sides; then
`face_gate.py cmp p2_24 p2_24`.
⚠ `cmp p2_24 p2_24` alone is vacuous: identical md5 lists leave no thermo.dat to column-compare, so it reads IDENTICAL
before and after the fix. The self-comparison is the witness. Bar: 0 columns moved under the fix; cmp IDENTICAL.
Usage: z1.py > Z1.txt   (exit 1 on a failed bar)"""
import importlib.machinery
import importlib.util
import io
import os
import sys
from contextlib import redirect_stdout

HERE = os.path.dirname(os.path.abspath(__file__))
W = os.path.abspath(os.path.join(HERE, "..", "tree_sweep_0922", "witness"))
loader = importlib.machinery.SourceFileLoader("face_gate", os.path.join(W, "face_gate.py"))
spec = importlib.util.spec_from_loader("face_gate", loader)
fg = importlib.util.module_from_spec(spec)
loader.exec_module(fg)


def main():
    path = os.path.join(W, "face_p2_24", "BEAM_S", "run", "thermo.dat")
    # TWO parses, as cmp makes them (one per side): Python's list comparison tests identity first, so one parsed list
    # against itself hides the bug -- float("nan") is a fresh object on each parse, and only then is NaN != NaN
    c = fg.thermo_cols(path)
    c2 = fg.thermo_cols(path)
    cols = [h for h in c if h not in fg.TIMERS]
    old = [h for h in cols if c[h] != c2[h]]
    new = [h for h in cols if fg.col_moved(c[h], c2[h])[0]]
    nanc = [h for h in cols if any(x != x for x in c[h])]
    print("P2-25 z1: face_p2_24/BEAM_S/run/thermo.dat against itself (%d columns, timers excluded; %d rows)" % (
        len(cols), len(next(iter(c.values())))))
    print("  columns holding a NaN: %d (%s)" % (len(nanc), " ".join(nanc)))
    print("  the OLD rule (list !=): %d moved (%s)" % (len(old), " ".join(old)))
    print("  face_gate.col_moved (NaN == NaN): %d moved" % len(new))
    buf = io.StringIO()
    with redirect_stdout(buf):
        ok = fg.cmp("p2_24", "p2_24")
    verdict = [ln for ln in buf.getvalue().splitlines() if ln.startswith("VERDICT")]
    print("  cmp p2_24 p2_24: %s" % (verdict[0] if verdict else "no verdict line"))
    good = not new and ok
    print("z1 %s (bar: 0 moved under the fix; cmp IDENTICAL)" % ("PASS" if good else "FAIL"))
    return 0 if good else 1


if __name__ == "__main__":
    sys.exit(main())
