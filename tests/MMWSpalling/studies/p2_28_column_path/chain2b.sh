#!/bin/bash
# P2-28 2b -- after the stage-2 chain: the b1 audit re-run with its declaration corrected (the packet's one declared
# mover, FEET_S's h_col CSV, never existed on the face list: GONE_PATHS = none, audit_components.py), its planted
# mutation, then b5 (CAL after vs before), which the chain skipped because of that audit's first FAIL. Same binary
# (BINARY_2.txt); same log (CHAIN.log).
set -u
ROOT=/Users/tzetze20/amr_tools/alamo; cd $ROOT
S=tests/MMWSpalling/studies/p2_28_column_path
V=tests/MMWSpalling/validation
BIN=bin/mmwspalling-3d-g++
PY=/Users/tzetze20/Desktop/code/.venv/bin/python; LOG=$S/CHAIN.log
say() { echo "$(date +%T) [s2b] $*" | tee -a $LOG; }
caffeinate -dis -w $$ &
B=$(shasum -a 256 $BIN | cut -c1-64)
grep -q "$B" $S/BINARY_2.txt || { say "STOP: bin/ is not BINARY_2.txt's binary"; exit 1; }
FAIL=0
$PY $S/audit_components.py face p2_28a p2_28 > $S/K2b_audit.log 2>&1; rc=$?
say "b1 audit (declaration corrected: no mover) rc $rc: $(tail -1 $S/K2b_audit.log)"
[ $rc = 0 ] || FAIL=1
if $PY $S/audit_components.py face p2_28a p2_28 --mutate > $S/K2b_audit_mutate.log 2>&1; then
  say "FAIL: the planted mutation PASSED"; FAIL=1
else say "b1 --mutate: FAILS as it must ($(tail -1 $S/K2b_audit_mutate.log))"; fi
$PY $V/kant/grimsel_calibration/run.py CAL > $S/K2_cal_after.log 2>&1
say "b5 CAL after: rc $?"
$PY $S/cal_compare.py > $S/K2_cal_compare.log 2>&1; rc=$?
say "b5 compare rc $rc: $(grep -E 'event log|c1 level|GATE' $S/K2_cal_compare.log | tr '\n' ' ' | cut -c1-300)"
[ $rc = 0 ] || FAIL=1
[ $FAIL = 0 ] && say "2b DONE -- b1 and b5 hold" || say "2b DONE WITH A FAILURE -- see above"
