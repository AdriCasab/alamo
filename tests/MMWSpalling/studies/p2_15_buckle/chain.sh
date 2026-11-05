#!/bin/bash
# P2-15 in the packet's run order; STOPS at the first failed gate. Log: CHAIN.log.
#   ord 1 a0 lock (+ --pre), a1 | ord 2 a2 tree, a3 units (buckle_lock, face_buckle, face_cref re-run), a4, a5 -> SPAN.txt, a6
#   ord 3 b1-b4 (face_buckle's a3 run + the cross-binary check) | ord 4 KG*_S, KB*_S, S_h*, HOS_S -> RESULTS.txt
set -u
ROOT=/Users/tzetze20/amr_tools/alamo; cd $ROOT
S=tests/MMWSpalling/studies/p2_15_buckle
W=tests/MMWSpalling/studies/tree_sweep_0922/witness
TC=tests/MMWSpalling/studies/p2_9_face_removal/tree_compare.py
PY=/Users/tzetze20/Desktop/code/.venv/bin/python; LOG=$S/CHAIN.log
say()  { echo "$(date +%T) $*" | tee -a $LOG; }
stop() { say "STOP: $*"; exit 1; }
caffeinate -dis -w $$ &
B=$(shasum -a 256 bin/mmwspalling-3d-g++ | cut -c1-64)
say "P2-15 chain (binary ${B:0:16})"
grep -q "$B" $S/BINARY.txt || stop "bin/ is not the P2-15 build"
$PY tests/MMWSpalling/unit/buckle_lock/test --pre $S/src_pre > $S/A0.txt 2>&1 || stop "a0 ($S/A0.txt)"
say "a0 PASS"
$PY $S/reduce.py A1 > /dev/null; say "a1 written"
say "a2: tree sweep"
bash $W/guarded_post_all_par.sh p2_15 > $S/tree.log 2>&1
[ $? = 2 ] && stop "the binary moved during the sweep"
$PY $TC p2_15 p2_12r p2_11b d2aa d2t > $S/GATE_tree.txt 2>&1; cat $S/GATE_tree.txt >> $LOG
grep -q "vs p2_12r: 0 of 2639 changed, 0 only in new, 0 only in ref" $S/GATE_tree.txt || stop "a2 vs p2_12r"
grep -q "vs p2_11b: 0 of 2639 changed, 0 only in new, 0 only in ref" $S/GATE_tree.txt || stop "a2 vs p2_11b"
grep -q "vs d2aa: 26 of 2639 changed" $S/GATE_tree.txt || stop "a2 control d2aa"
grep -q "vs d2t: 863 of 2639 changed" $S/GATE_tree.txt || stop "a2 control d2t"
say "a2 PASS; a3: every kept unit"
$PY $S/run.py UNITS > $S/units.log 2>&1
$PY $S/reduce.py A3 > /dev/null; tail -1 $S/A3.txt | tee -a $LOG
grep -q "^a3 PASS" $S/A3.txt || stop "a3 ($S/A3.txt)"
say "a4: CAL, HOS_R (+ F0_k1 from a3)"
$PY $S/run.py CAL HOS_R >> $S/RUN.log 2>&1
$PY $S/reduce.py A4 > /dev/null; grep "^  \|^a4" $S/A4.txt | cut -c1-200 | tee -a $LOG
grep -q "^a4 PASS" $S/A4.txt || stop "a4 ($S/A4.txt)"
say "a5: the span-off Kant legs"
$PY $S/run.py KG130_0 KG140_0 KG150_0 KG160_0 KB150_0 KB160_0 >> $S/RUN.log 2>&1
$PY $S/derive_span.py > /dev/null; cat $S/A5.txt >> $LOG
grep -q "^a5 PASS" $S/A5.txt || stop "a5 ($S/A5.txt)"
cat $S/SPAN.txt | tee -a $LOG
$PY $S/reduce.py A6 | tee -a $LOG
say "ord 3: b"
$PY $S/reduce.py B > /dev/null; grep "^b\|FAIL" $S/B.txt | tee -a $LOG
grep -q "^b PASS" $S/B.txt || stop "b ($S/B.txt)"
say "ord 4: the span-on legs"
$PY $S/run.py KG50_S KG80_S KG110_S KG130_S KG140_S KG150_S KG160_S KB150_S KB160_S S_h150 S_h500 S_h1000 S_h2000 HOS_S >> $S/RUN.log 2>&1
$PY $S/reduce.py > $S/RESULTS.txt 2>&1
sed -n '/^c -- /,$p' $S/RESULTS.txt | cut -c1-260 | tee -a $LOG
say "CHAIN DONE"
