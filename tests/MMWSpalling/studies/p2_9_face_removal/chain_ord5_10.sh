#!/bin/bash
# P2-9 rev 2.1, ord 5-10 + a5b, in the packet's order with its gates enforced:
#   b1 (GATE_b.txt "b1  PASS") before the mesh legs; b2 (both mesh legs .done)
#   before C15n; d1 (SSTAR.txt S_STAR=, at most one added point) before P47n.
# Every leg under caffeinate (the harness also asserts it). Writes RUN_ord5_10.log.
ROOT=/Users/tzetze20/amr_tools/alamo
cd "$ROOT" || exit 1
PY=/Users/tzetze20/Desktop/code/.venv/bin/python
S=tests/MMWSpalling/studies/p2_9_face_removal
U=tests/MMWSpalling/unit/face_removal
say() { echo "$(date +%T) $*"; }

say "ord 5: K0g M0"
caffeinate -dis $PY $S/run_p2_9.py K0g M0
$PY $S/reduce_p2_9.py GATE > /dev/null
grep -q "b1  PASS" $S/GATE_b.txt || { say "STOP: b1 did not pass (GATE_b.txt)"; exit 1; }

say "ord 6: M0c M0f"
caffeinate -dis $PY $S/run_p2_9.py M0c M0f
$PY $S/reduce_p2_9.py MESH > /dev/null
[ -f $S/legs/M0c/.done ] && [ -f $S/legs/M0f/.done ] || { say "STOP: b2 -- a mesh leg did not reach 30 s"; exit 1; }

say "ord 7: C15n"
caffeinate -dis $PY $S/run_p2_9.py C15n
$PY $S/reduce_p2_9.py SSTAR > /dev/null
if ! grep -q "^S_STAR=" $S/SSTAR.txt; then
    add=$(grep -o "C20n\|C07n" $S/SSTAR.txt | head -1)
    say "d1 not bracketed: the ONE added point $add"
    caffeinate -dis $PY $S/run_p2_9.py $add
    $PY $S/reduce_p2_9.py SSTAR > /dev/null
    grep -q "^S_STAR=" $S/SSTAR.txt || { say "STOP: d1 still not bracketed"; exit 1; }
fi

say "ord 8: P47r"
caffeinate -dis $PY $S/run_p2_9.py P47r
if grep -q "P47n RUNS" $S/SSTAR.txt; then
    say "ord 9: P47n"
    caffeinate -dis $PY $S/run_p2_9.py --probe P47n
    caffeinate -dis $PY $S/run_p2_9.py P47n
fi
$PY $S/reduce_p2_9.py > $S/RESULTS.txt 2>&1
say "RESULTS written"

say "a5b: the slow tilts at 0.5 mm"
caffeinate -dis $PY $U/test --a5b > $U/a5b.log 2>&1
cp $U/a5b.log $S/A5B.txt
say "CHAIN DONE"
