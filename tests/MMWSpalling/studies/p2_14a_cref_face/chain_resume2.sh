#!/bin/bash
# P2-14a resumed at a4 after a3 PASSED on the corrected comparator (no unit re-run). Log: CHAIN.log.
# read). Re-runs that unit only, re-scores a3, then the chain from a4 on. Log: CHAIN.log.
set -u
ROOT=/Users/tzetze20/amr_tools/alamo; cd $ROOT
S=tests/MMWSpalling/studies/p2_14a_cref_face
PY=/Users/tzetze20/Desktop/code/.venv/bin/python; LOG=$S/CHAIN.log
say()  { echo "$(date +%T) $*" | tee -a $LOG; }
stop() { say "STOP: $*"; exit 1; }
caffeinate -dis -w $$ &
say "RESUME: lateral_flux_lock amended (option 1); a3 re-scored with line numbers masked"
say "RESUME 2: a3 PASS on the corrected comparator (full-log check for lateral_flux_lock)"
say "a4-a8: the gate legs"
$PY $S/run.py F2_k0_pre F2_k1_pre F0_k0_pre F0_k1_pre CAL_p27_mu_pre CAL CAL_det CAL_h2000 CAL_2mm CAL_p27 CAL_p48 HOS_R CAL_p27_mu >> $S/RUN.log 2>&1
$PY $S/reduce.py GATES > /dev/null; grep "^  a\|^a4" $S/GATES_a4_a8.txt | cut -c1-220 | tee -a $LOG
grep -q "^a4-a8 PASS" $S/GATES_a4_a8.txt || stop "a4-a8 ($S/GATES_a4_a8.txt)"
say "ord 3: b legs"
$PY $S/run.py F2_k1_np1 F2_k1_g8 M2_k1 M05_k1 >> $S/RUN.log 2>&1
$PY $S/reduce.py B > /dev/null; grep "^b" $S/B.txt | cut -c1-220 | tee -a $LOG
grep -q "^b PASS" $S/B.txt || stop "b ($S/B.txt)"
say "ord 4: c legs"
$PY $S/run.py C0_h150 C2_h150 C2_h500 C0_h1000 C2_h1000 C2_h2000 HOS_d2 >> $S/RUN.log 2>&1
$PY $S/reduce.py > $S/RESULTS.txt 2>&1
grep "^c\|^   d " $S/RESULTS.txt | cut -c1-220 | tee -a $LOG
say "CHAIN DONE"
