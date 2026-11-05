#!/usr/bin/env python3
"""P2-9 (a1, the unit half): the unit tests the 2639 sweep misses (DR-2), on this
binary (bin/*.log) vs P2-8's shipping binary (studies/p2_8_face_list/bin2/*.log).
Text-identical after masking what any edit moves without changing behaviour:
source line numbers (`File.H:123`) and wall-clock seconds. face_list is a key-ON
test -- its legs run the line rule now -- so its PASS lines are compared, not
its numbers. Writes UNIT_PARITY.txt."""
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
OLD = os.path.join(HERE, "..", "p2_8_face_list", "bin2")
NEW = os.path.join(HERE, "bin")
TESTS = ["side_face_flux", "side_face_area", "side_face_step_height", "lateral_drive",
         "wall_march", "face_list", "lateral_flux_lock"]


def norm(t):
    t = re.sub(r"\.H\s*:\s*\d+", ".H:N", t)
    t = re.sub(r"\.H(\x1b\[0m)?:\d+", ".H:N", t)
    t = re.sub(r"\b\d+(\.\d+)?\s*s\b", "Ns", t)
    t = re.sub(r"rc=\d+\s+\d+s", "rc", t)
    return [ln.rstrip() for ln in t.split("\n")]


def main():
    out = [__doc__.strip().split("\n")[0], "=" * 88]
    for t in TESTS:
        po, pn = os.path.join(OLD, t + ".log"), os.path.join(NEW, t + ".log")
        if not (os.path.exists(po) and os.path.exists(pn)):
            out.append("  %-22s missing (%s)" % (t, "old" if not os.path.exists(po) else "new"))
            continue
        a, b = norm(open(po).read()), norm(open(pn).read())
        if t in ("face_list", "lateral_flux_lock"):
            a = [x for x in a if x.strip().startswith(("PASS", "FAIL"))]
            b = [x for x in b if x.strip().startswith(("PASS", "FAIL"))]
        d = [i for i in range(max(len(a), len(b))) if (a[i] if i < len(a) else None) != (b[i] if i < len(b) else None)]
        verdict = "IDENTICAL" if not d else "%d line(s) differ" % len(d)
        out.append("  %-22s %s (%d vs %d lines)" % (t, verdict, len(a), len(b)))
        for i in d[:6]:
            out.append("      old: %s" % (a[i] if i < len(a) else "<none>")[:150])
            out.append("      new: %s" % (b[i] if i < len(b) else "<none>")[:150])
    s = "\n".join(out)
    print(s)
    open(os.path.join(HERE, "UNIT_PARITY.txt"), "w").write(s + "\n")


if __name__ == "__main__":
    main()
