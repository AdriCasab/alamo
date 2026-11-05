#!/bin/bash
# continuation of chain_restore.sh after the r5 window fix (r1, r1b, r2 PASSED at 19:06-19:10, CHAIN.log)
set -u
ROOT=/Users/tzetze20/amr_tools/alamo; cd $ROOT
S=tests/MMWSpalling/studies/p2_12_retire_path/restore
PY=/Users/tzetze20/Desktop/code/.venv/bin/python; LOG=$S/CHAIN.log
say()  { echo "$(date +%T) $*" | tee -a $LOG; }
stop() { say "STOP: $*"; exit 1; }
caffeinate -dis -w $$ &
say "continuation (binary $(shasum -a 256 bin/mmwspalling-3d-g++ | cut -c1-16)); r5 re-scored on the existing 15 s rerun (prefix rule, declared)"
$PY $S/gates_restore.py X1CMP > $S/x1repro.log 2>&1 || stop "r5 ($S/X1REPRO.txt)"
say "r5 PASS"
$PY $S/gates_restore.py W > $S/W.log 2>&1 || stop "r3 ($S/W_p2_12r_vs_c2.txt)"
say "r3 PASS"
$PY $S/gates_restore.py UNITS > $S/units.log 2>&1 || stop "r4 ($S/UNITS_p2_12r_vs_p2_11b.txt)"
say "r4 PASS"
say "ALL RESTORE GATES PASS"
