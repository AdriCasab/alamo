#!/bin/bash
# P2-11b ord 4-7 on the Stage B binary: X1_S14, X1_S20 -> a6 (gate) -> b1 (S*, at most ONE added point,
# gate) -> X1 (120 s) -> c1-c5 into RESULTS.txt. Stops at a failed gate. Log: CHAIN.log.
set -u
ROOT=/Users/tzetze20/amr_tools/alamo; cd $ROOT
S=tests/MMWSpalling/studies/p2_11b_path
PY=/Users/tzetze20/Desktop/code/.venv/bin/python; LOG=$S/CHAIN.log
say()  { echo "$(date +%T) $*" | tee -a $LOG; }
stop() { say "STOP: $*"; exit 1; }
caffeinate -dis -w $$ &
say "ord 4 (binary $(shasum -a 256 bin/mmwspalling-3d-g++ | cut -c1-16)): X1_S14, X1_S20"
$PY $S/run_p2_11b.py X1_S14 X1_S20 >> $S/RUN_ord4.log 2>&1 || stop "a 30 s X leg failed (RUN_ord4.log)"
$PY $S/reduce_p2_11b.py A6 > /dev/null; tail -1 $S/A6.txt | tee -a $LOG
grep -q "^a6 PASS" $S/A6.txt || stop "a6 ($S/A6.txt)"
$PY $S/reduce_p2_11b.py SSTAR > /dev/null; cat $S/SSTAR.txt >> $LOG
if ! grep -q "^S_STAR=" $S/SSTAR.txt; then
    add=$(grep "adds ONE point" $S/SSTAR.txt | grep -o "X1_S[0-9]*" | head -1)
    [ -n "$add" ] || stop "b1: not bracketed after the added point ($S/SSTAR.txt)"
    say "b1: adding $add"
    $PY $S/run_p2_11b.py $add >> $S/RUN_ord4.log 2>&1 || stop "$add failed"
    $PY $S/reduce_p2_11b.py SSTAR > /dev/null; cat $S/SSTAR.txt >> $LOG
    grep -q "^S_STAR=" $S/SSTAR.txt || stop "b1 ($S/SSTAR.txt)"
fi
say "ord 6: X1 at S* (120 s)"
$PY $S/run_p2_11b.py X1 >> $S/RUN_ord4.log 2>&1 || stop "X1 failed (RUN_ord4.log)"
$PY $S/reduce_p2_11b.py A6 > /dev/null; tail -1 $S/A6.txt | tee -a $LOG
$PY $S/reduce_p2_11b.py > $S/RESULTS.txt 2>&1
grep "e1 inputs" $S/RESULTS.txt | tee -a $LOG
say "CHAIN DONE"
