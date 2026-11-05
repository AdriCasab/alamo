#!/bin/bash
# P2-11b ord 2 on the Stage B binary, face_heat unset: a3 (tree sweep vs p2_10_c2, controls; W vs W_c2),
# a4 (G1_C vs G1), a4b (P47n_C vs P47n). Stops at a failed gate. Log: CHAIN.log.
set -u
ROOT=/Users/tzetze20/amr_tools/alamo; cd $ROOT
S=tests/MMWSpalling/studies/p2_11b_path; W=tests/MMWSpalling/studies/tree_sweep_0922/witness
PY=/Users/tzetze20/Desktop/code/.venv/bin/python; LOG=$S/CHAIN.log
TC=tests/MMWSpalling/studies/p2_9_face_removal/tree_compare.py
say()  { echo "$(date +%T) $*" | tee -a $LOG; }
stop() { say "STOP: $*"; exit 1; }
caffeinate -dis -w $$ &
say "ord 2 (binary $(shasum -a 256 bin/mmwspalling-3d-g++ | cut -c1-16))"
bash $W/guarded_post_all_par.sh p2_11b > $S/tree.log 2>&1
[ $? = 2 ] && stop "the binary moved during the sweep"
$PY $TC p2_11b p2_10_c2 d2aa d2t | tee $S/GATE_tree.txt | tee -a $LOG
grep -q "vs p2_10_c2: 0 of 2639 changed, 0 only in new, 0 only in ref" $S/GATE_tree.txt || stop "a3 tree"
grep -q "vs d2aa: 26 of 2639 changed" $S/GATE_tree.txt || stop "a3 control d2aa"
grep -q "vs d2t: 863 of 2639 changed" $S/GATE_tree.txt || stop "a3 control d2t"
$PY $S/reduce_p2_11b.py W > $S/W.log 2>&1; tail -1 $S/W_p2_11b_vs_c2.txt | tee -a $LOG
tail -1 $S/W_p2_11b_vs_c2.txt | grep -q PASS || stop "a3 W ($S/W_p2_11b_vs_c2.txt)"
say "a3 PASS; a4: G1_C, P47n_C (120 s each)"
$PY $S/run_p2_11b.py G1_C P47n_C >> $S/RUN_ord2.log 2>&1 || stop "a control leg failed (RUN_ord2.log)"
$PY $S/reduce_p2_11b.py IDENTITY > /dev/null; grep "^a4" $S/IDENTITY.txt | tee -a $LOG
grep -q "^a4  PASS" $S/IDENTITY.txt || stop "a4 ($S/IDENTITY.txt)"
say "ORD 2 DONE"
