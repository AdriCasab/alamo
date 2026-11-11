#!/bin/bash
# P2-13 ord 2: the gates a2, a3, a4 (no build; the binary is restore/BINARY.txt 0a96a095...). Stops at a failed gate.
set -u
ROOT=/Users/tzetze20/amr_tools/alamo; cd $ROOT
K=tests/MMWSpalling/validation/kant/grimsel_calibration
PY=/Users/tzetze20/Desktop/code/.venv/bin/python; LOG=$K/CHAIN.log
say()  { echo "$(date +%T) $*" | tee -a $LOG; }
stop() { say "STOP: $*"; exit 1; }
caffeinate -dis -w $$ &
B=$(shasum -a 256 bin/mmwspalling-3d-g++ | cut -c1-64)
say "ord 2 (binary ${B:0:16})"
grep -q "$B" tests/MMWSpalling/studies/p2_12_retire_path/restore/BINARY.txt || stop "binary is not restore/BINARY.txt"
$PY $K/run.py SMOKE HOS_R >> $K/RUN_ord2.log 2>&1 || stop "SMOKE or HOS_R failed (RUN_ord2.log)"
$PY $K/reduce.py A2 > /dev/null; tail -1 $K/A2.txt | tee -a $LOG
grep -q "^a2 PASS" $K/A2.txt || stop "a2 ($K/A2.txt)"
$PY $K/reduce.py A3 > /dev/null; tail -1 $K/A3.txt | tee -a $LOG
grep -q "^a3 PASS" $K/A3.txt || stop "a3 ($K/A3.txt)"
$PY $K/run.py UNITS > $K/units.log 2>&1
tail -1 $K/UNITS_p2_13_vs_p2_12r.txt | tee -a $LOG
grep -q "^a4 (every test): PASS" $K/UNITS_p2_13_vs_p2_12r.txt || stop "a4"
say "ord 3-4: SMOKE legs (10 s, witness), then the legs"
$PY $K/run.py SMOKE_CAL SMOKE_CAL_det SMOKE_CAL_h2000 SMOKE_CAL_p27 SMOKE_CAL_p48 SMOKE_CAL_2mm SMOKE_CAL_05mm >> $K/RUN_ord3.log 2>&1 || stop "a smoke leg failed"
for l in CAL CAL_det CAL_h2000 CAL_p27 CAL_p48 CAL_2mm CAL_05mm; do
  n=$(($(wc -l < $K/legs/SMOKE_$l/run_removal_events.csv) - 1)); say "  smoke $l: $n firings by 10 s"
  [ $n -gt 0 ] || stop "smoke $l: no firing (witness zero)"
done
$PY $K/run.py CAL CAL_det CAL_h2000 CAL_p27 CAL_p48 CAL_2mm CAL_05mm >> $K/RUN_ord3.log 2>&1 || stop "a leg failed"
$PY $K/reduce.py B1 | tail -1 | tee -a $LOG
$PY $K/reduce.py > $K/RESULTS.txt 2>&1
grep "^c1\|^e1" $K/RESULTS.txt | tee -a $LOG
say "CHAIN DONE (d1 runs only if c1 is OUTSIDE)"
