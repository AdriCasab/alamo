#!/bin/bash
# P2-14a in the packet's run order; STOPS at the first failed gate. Log: CHAIN.log.
#   ord 1 a0 (lock + --pre), a1 | ord 2 a2 tree, a3 units (face_cref included), a4-a8 legs
#   ord 3 b1-b5 | ord 4 c1-c4 legs -> RESULTS.txt
set -u
ROOT=/Users/tzetze20/amr_tools/alamo; cd $ROOT
S=tests/MMWSpalling/studies/p2_14a_cref_face
W=tests/MMWSpalling/studies/tree_sweep_0922/witness
TC=tests/MMWSpalling/studies/p2_9_face_removal/tree_compare.py
PY=/Users/tzetze20/Desktop/code/.venv/bin/python; LOG=$S/CHAIN.log
say()  { echo "$(date +%T) $*" | tee -a $LOG; }
stop() { say "STOP: $*"; exit 1; }
caffeinate -dis -w $$ &
B=$(shasum -a 256 bin/mmwspalling-3d-g++ | cut -c1-64)
say "P2-14a chain (binary ${B:0:16})"
grep -q "$B" $S/BINARY.txt || stop "bin/ is not the P2-14a build ($S/BINARY.txt)"
$PY tests/MMWSpalling/unit/constraint_ref_lock/test --pre $S/src_pre > $S/A0.txt 2>&1 || stop "a0 ($S/A0.txt)"
say "a0 PASS"
$PY $S/reduce.py A1 > /dev/null; tail -1 $S/A1.txt | tee -a $LOG
say "a2: tree sweep"
bash $W/guarded_post_all_par.sh p2_14a > $S/tree.log 2>&1
[ $? = 2 ] && stop "the binary moved during the sweep"
$PY $TC p2_14a p2_12r p2_11b d2aa d2t > $S/GATE_tree.txt 2>&1; cat $S/GATE_tree.txt >> $LOG
grep -q "vs p2_12r: 0 of 2639 changed, 0 only in new, 0 only in ref" $S/GATE_tree.txt || stop "a2 vs p2_12r"
grep -q "vs p2_11b: 0 of 2639 changed, 0 only in new, 0 only in ref" $S/GATE_tree.txt || stop "a2 vs p2_11b"
grep -q "vs d2aa: 26 of 2639 changed" $S/GATE_tree.txt || stop "a2 control d2aa"
grep -q "vs d2t: 863 of 2639 changed" $S/GATE_tree.txt || stop "a2 control d2t"
say "a2 PASS; a3: every kept unit"
$PY $S/run.py UNITS > $S/units.log 2>&1
$PY $S/reduce.py A3 > /dev/null; tail -1 $S/A3.txt | tee -a $LOG
grep -q "^a3 PASS" $S/A3.txt || stop "a3 ($S/A3.txt)"
say "a4-a8: the gate legs"
$PY $S/run.py F2_k0_pre F2_k1_pre F0_k0_pre F0_k1_pre CAL_p27_mu_pre CAL CAL_det CAL_h2000 CAL_2mm CAL_p27 CAL_p48 HOS_R CAL_p27_mu >> $S/RUN.log 2>&1
$PY $S/reduce.py GATES > /dev/null; grep "^  a\|^a4" $S/GATES_a4_a8.txt | cut -c1-200 | tee -a $LOG
grep -q "^a4-a8 PASS" $S/GATES_a4_a8.txt || stop "a4-a8 ($S/GATES_a4_a8.txt)"
say "ord 3: b legs"
$PY $S/run.py F2_k1_np1 F2_k1_g8 M2_k1 M05_k1 >> $S/RUN.log 2>&1
$PY $S/reduce.py B > /dev/null; grep "^b" $S/B.txt | cut -c1-200 | tee -a $LOG
grep -q "^b PASS" $S/B.txt || stop "b ($S/B.txt)"
say "ord 4: c legs"
$PY $S/run.py C0_h150 C2_h150 C2_h500 C0_h1000 C2_h1000 C2_h2000 HOS_d2 >> $S/RUN.log 2>&1
$PY $S/reduce.py > $S/RESULTS.txt 2>&1
grep "^c\|^   d " $S/RESULTS.txt | cut -c1-200 | tee -a $LOG
say "CHAIN DONE"
