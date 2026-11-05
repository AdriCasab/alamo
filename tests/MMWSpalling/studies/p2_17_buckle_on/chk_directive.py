#!/usr/bin/env python3
"""P2-17 z0 -- the `# check_shared: allow KEY` directive, two-sided (check (i) only; no run needed).

  D1 the edited CAL deck (grimsel_calibration/input)            -> (i) PASS, exactly 1 declared
  D2 a scratch copy WITHOUT the directive line                  -> (i) FAIL naming `spall.buckle_span (deck)`
  D3 a scratch copy whose extra directive names an UNSET key    -> (i) FAIL naming that key
(D3 keeps the real directive, so its only fault is the unset key.) check_shared's (ii) needs a run; here it is
pointed at an empty --run dir, so only the incl and (i) lines are scored. Writes Z0.txt; exit 1 on a failed row.
"""
import os
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
V = os.path.join(ROOT, "tests", "MMWSpalling", "validation")
CHECK = os.path.join(V, "shared", "check_shared.py")
DECK = os.path.join(V, "kant", "grimsel_calibration", "input")
DIRECTIVE = "# check_shared: allow spall.buckle_span\n"
UNSET = "weibull.m"


def check(deck):
    empty = tempfile.mkdtemp(prefix="z0run_")
    p = subprocess.run([sys.executable, CHECK, deck, "--run", empty], cwd=ROOT, capture_output=True, text=True)
    out = p.stdout + p.stderr
    li = [l for l in out.splitlines() if l.strip().startswith("(i)")]
    inc = [l for l in out.splitlines() if l.strip().startswith("incl")]
    return (li[0] if li else "(no (i) line)"), (inc[0] if inc else "(no incl line)"), out


def main():
    src = open(DECK).read()
    assert src.count(DIRECTIVE) == 1, "the CAL deck does not carry exactly one directive"
    d = tempfile.mkdtemp(prefix="z0deck_")
    d2, d3 = os.path.join(d, "no_directive"), os.path.join(d, "unset_directive")
    open(d2, "w").write(src.replace(DIRECTIVE, ""))
    open(d3, "w").write(src.replace(DIRECTIVE, DIRECTIVE + "# check_shared: allow %s\n" % UNSET))
    rows = []
    l1, i1, _ = check(DECK)
    rows.append(("D1 edited CAL deck: (i) PASS, 1 declared",
                 " PASS " in l1 and "(1 declared" in l1 and " PASS " in i1, l1))
    l2, _, _ = check(d2)
    rows.append(("D2 no directive: (i) FAIL naming spall.buckle_span (deck)",
                 " FAIL " in l2 and "spall.buckle_span (deck)" in l2, l2))
    l3, _, _ = check(d3)
    rows.append(("D3 directive naming an unset key (%s): (i) FAIL naming it" % UNSET,
                 " FAIL " in l3 and re.search(r"nobody sets: .*%s" % re.escape(UNSET), l3) is not None, l3))
    L = ["z0 -- the check_shared directive, two-sided", "=" * 100]
    for nm, ok, ln in rows:
        L.append("  %s  %s" % ("PASS" if ok else "FAIL", nm))
        L.append("        %s" % ln.strip())
    ok = all(r[1] for r in rows)
    L.append("z0 %s (gate)" % ("PASS" if ok else "FAIL"))
    t = "\n".join(L)
    print(t)
    open(os.path.join(HERE, "Z0.txt"), "w").write(t + "\n")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
