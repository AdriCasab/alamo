#!/bin/bash
# P2-16 in the packet's run order; STOPS at the first failed gate. Log: CHAIN.log.
#   ord 1 a0 (retired_keys source half + --pre, lateral_flux_lock, its mutation control), c0/c0b (rules.py), R_30
#   ord 2 a1 tree, a2 units vs UNIT_p2_15, a3 CAL/HOS_R/F0_k1/KB150_0, a4 retired_keys run half
#   ord 3 c1 KB150_m30, KB160_m30, CAL_m30 | ord 4 RESULTS.txt (e1)
set -u
ROOT=/Users/tzetze20/amr_tools/alamo; cd $ROOT
S=tests/MMWSpalling/studies/p2_16_span_rule
U=tests/MMWSpalling/unit
W=tests/MMWSpalling/studies/tree_sweep_0922/witness
TC=tests/MMWSpalling/studies/p2_9_face_removal/tree_compare.py
PY=/Users/tzetze20/Desktop/code/.venv/bin/python; LOG=$S/CHAIN.log
say()  { echo "$(date +%T) $*" | tee -a $LOG; }
stop() { say "STOP: $*"; exit 1; }
caffeinate -dis -w $$ &
B=$(shasum -a 256 bin/mmwspalling-3d-g++ | cut -c1-64)
say "P2-16 chain (binary ${B:0:16})"
grep -q "$B" $S/BINARY.txt || stop "bin/ is not the P2-16 build"

say "ord 1: a0, c0, c0b, R_30"
$PY $U/retired_keys/test --source --pre $S/src_pre > $S/A0_retired_keys.log 2>&1
$PY $U/lateral_flux_lock/test --check > $S/A0_lfl.log 2>&1
$PY $S/lfl_control/lfl_mut.py > /dev/null 2>&1
$PY $S/reduce.py A0 > /dev/null; tail -1 $S/A0.txt | tee -a $LOG
grep -q "^a0 PASS" $S/A0.txt || stop "a0 ($S/A0.txt)"
$PY $S/rules.py > /dev/null && grep "^   R\|^c0b" $S/RULES.txt | cut -c1-200 | tee -a $LOG
$PY $S/bethel_R.py | tee -a $LOG || stop "bethel_R.py"

say "ord 2: a1 tree sweep"
bash $W/guarded_post_all_par.sh p2_16 > $S/tree.log 2>&1
[ $? = 2 ] && stop "the binary moved during the sweep"
$PY $TC p2_16 p2_12r p2_11b d2aa d2t > $S/GATE_tree.txt 2>&1; cat $S/GATE_tree.txt >> $LOG
grep -q "vs p2_12r: 0 of 2639 changed, 0 only in new, 0 only in ref" $S/GATE_tree.txt || stop "a1 vs p2_12r"
grep -q "vs p2_11b: 0 of 2639 changed, 0 only in new, 0 only in ref" $S/GATE_tree.txt || stop "a1 vs p2_11b"
grep -q "vs d2aa: 26 of 2639 changed" $S/GATE_tree.txt || stop "a1 control d2aa"
grep -q "vs d2t: 863 of 2639 changed" $S/GATE_tree.txt || stop "a1 control d2t"
say "a1 PASS; a2: every kept unit"
$PY $S/run.py UNITS > $S/units.log 2>&1
$PY $S/reduce.py A2 > /dev/null; tail -1 $S/A2.txt | tee -a $LOG
grep -q "^a2 PASS" $S/A2.txt || stop "a2 ($S/A2.txt)"
say "a3: CAL, HOS_R, F0_k1, KB150_0"
$PY $S/run.py CAL HOS_R F0_k1 KB150_0 >> $S/RUN.log 2>&1
$PY $S/reduce.py A3 > /dev/null; grep "^  \|^a3" $S/A3.txt | cut -c1-220 | tee -a $LOG
grep -q "^a3 PASS" $S/A3.txt || stop "a3 ($S/A3.txt)"
say "a4: retired_keys run half"
$PY $U/retired_keys/test --run > $S/A4_retired_keys.log 2>&1
$PY $S/reduce.py A4 > /dev/null; tail -1 $S/A4.txt | tee -a $LOG
grep -q "^a4 PASS" $S/A4.txt || stop "a4 ($S/A4.txt)"

say "ord 3: c1 legs"
$PY $S/run.py KB150_m30 KB160_m30 CAL_m30 >> $S/RUN.log 2>&1
$PY $S/reduce.py C1 > /dev/null; cut -c1-220 $S/C1.txt | tee -a $LOG
say "ord 4: RESULTS.txt"
$PY $S/reduce.py > $S/RESULTS.txt 2>&1
sed -n '/^e1 VERDICT/,$p' $S/RESULTS.txt | tee -a $LOG
say "CHAIN DONE"
