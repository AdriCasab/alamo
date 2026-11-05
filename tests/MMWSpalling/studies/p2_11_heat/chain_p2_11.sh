#!/bin/bash
# P2-11 ord 2-5 after the two 30 s legs: a1 gate -> b1 (S*, at most ONE added point) -> the 120 s leg
# -> a1 on it -> c1-c3, e1. Stops at the first failed gate. Log: CHAIN.log.
set -u
ROOT=/Users/tzetze20/amr_tools/alamo; cd $ROOT
S=tests/MMWSpalling/studies/p2_11_heat
PY=/Users/tzetze20/Desktop/code/.venv/bin/python; LOG=$S/CHAIN.log
say()  { echo "$(date +%T) $*" | tee -a $LOG; }
stop() { say "STOP: $*"; exit 1; }
caffeinate -dis -w $$ &
say "P2-11 chain (binary $(shasum -a 256 bin/mmwspalling-3d-g++ | cut -c1-16))"
$PY $S/run_p2_11.py G1_S14 G1_S20 >> $S/RUN_ord2.log 2>&1 || stop "a 30 s leg failed (RUN_ord2.log)"
$PY $S/reduce_p2_11.py A1 > /dev/null; tail -1 $S/A1.txt | tee -a $LOG
grep -q "^a1 PASS" $S/A1.txt || stop "a1 ($S/A1.txt)"
$PY $S/reduce_p2_11.py SSTAR > /dev/null; cat $S/SSTAR.txt >> $LOG
if ! grep -q "^S_STAR=" $S/SSTAR.txt; then
    add=$(grep -o "G1_S[0-9]*" <(grep "adds ONE point" $S/SSTAR.txt) | head -1)
    [ -n "$add" ] || stop "b1: not bracketed after the added point ($S/SSTAR.txt)"
    say "b1: adding $add"
    $PY $S/run_p2_11.py $add >> $S/RUN_ord2.log 2>&1 || stop "$add failed"
    $PY $S/reduce_p2_11.py SSTAR > /dev/null; cat $S/SSTAR.txt >> $LOG
    grep -q "^S_STAR=" $S/SSTAR.txt || stop "b1 ($S/SSTAR.txt)"
fi
say "ord 4: G1 at S* (120 s)"
$PY $S/run_p2_11.py G1 > $S/RUN_ord4.log 2>&1 || stop "G1 failed (RUN_ord4.log)"
cat $S/RUN_ord4.log >> $LOG
$PY $S/reduce_p2_11.py A1 > /dev/null; tail -1 $S/A1.txt | tee -a $LOG
$PY $S/reduce_p2_11.py > $S/RESULTS.txt 2>&1
grep "e1" $S/RESULTS.txt | tee -a $LOG
say "CHAIN DONE"
