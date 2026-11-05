#!/usr/bin/env python3
"""P2-16 a0: the mutation control for unit/lateral_flux_lock's restated P2-2a clause.
Plants a WRONG-CELL reference in the Sp_field cell-sample branch's bracket --
(T_a(i_d, j_d, k_d) - T0_at(i, j, k)) -- in a scratch copy of MMWSpalling.H, runs the
lock (--check: source clauses only) against it, and writes CONTROL.txt. The control
holds iff the lock FAILS and the only failed check is the restated clause."""
import os, re, subprocess, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", "..", ".."))
LOCK = os.path.join(ROOT, "tests", "MMWSpalling", "unit", "lateral_flux_lock", "test")
MMW = os.path.join(ROOT, "src", "Integrator", "MMWSpalling.H")
GOOD = "(T_a(i_d, j_d, k_d) - T0_at(i_d, j_d, k_d))"
BAD = "(T_a(i_d, j_d, k_d) - T0_at(i, j, k))"
src = open(MMW).read()
assert src.count(GOOD) == 1, "the bracket to mutate is not unique"
d = tempfile.mkdtemp(prefix="lflmut_")
mut = os.path.join(d, "MMWSpalling.H")
open(mut, "w").write(src.replace(GOOD, BAD))
lock = open(LOCK).read()
line = 'MMW = os.path.join(ROOT, "src", "Integrator", "MMWSpalling.H")'
assert lock.count(line) == 1
drv = os.path.join(d, "lock_mut.py")
open(drv, "w").write(lock.replace(line, "MMW = %r" % mut).replace(
    'HERE = os.path.dirname(os.path.abspath(__file__))', "HERE = %r" % os.path.dirname(LOCK)))
r = subprocess.run([sys.executable, drv, "--check"], cwd=ROOT, capture_output=True, text=True)
out = r.stdout + r.stderr
fails = [l for l in out.splitlines() if l.startswith("  FAIL")]
ok = r.returncode != 0 and len(fails) == 1 and "BOTH Sp_field branches subtract the STORED field" in fails[0]
with open(os.path.join(HERE, "CONTROL.txt"), "w") as f:
    f.write("lateral_flux_lock with a planted WRONG-CELL reference in the Sp_field cell-sample bracket:\n")
    f.write("  %s  ->  %s\n" % (GOOD, BAD))
    f.write("the restated P2-2a clause must FAIL, and nothing else\n\n")
    f.write("\n".join(fails) + "\n")
    f.write([l for l in out.splitlines() if l.startswith(("FAIL:", "PASS:"))][-1] + "\n\n")
    f.write("control %s (rc=%d, %d failed check(s))\n" % ("HOLDS" if ok else "BROKEN", r.returncode, len(fails)))
print(open(os.path.join(HERE, "CONTROL.txt")).read())
sys.exit(0 if ok else 1)
