#!/bin/bash
# P2-19 in the packet's run order; STOPS at the first failed gate. Log: CHAIN.log.
#   ord 1 z0 (unit/profile_depth) | ord 2 a0 tree gate vs p2_18, a1 face gate cmp p2_18 p2_19, a2 T00 x3
#   ord 3 a3 CAL, then HOS and C0H150 (reports c1-c4) | ord 4 the locks, RESULTS.txt (e0)
set -u
ROOT=/Users/tzetze20/amr_tools/alamo; cd $ROOT
S=tests/MMWSpalling/studies/p2_19_flake_depth
U=tests/MMWSpalling/unit
W=tests/MMWSpalling/studies/tree_sweep_0922/witness
TC=tests/MMWSpalling/studies/p2_9_face_removal/tree_compare.py
PY=/Users/tzetze20/Desktop/code/.venv/bin/python; LOG=$S/CHAIN.log
say()  { echo "$(date +%T) $*" | tee -a $LOG; }
stop() { say "STOP: $*"; exit 1; }
caffeinate -dis -w $$ &
B=$(shasum -a 256 bin/mmwspalling-3d-g++ | cut -c1-64)
say "P2-19 chain (binary ${B:0:16})"
grep -q "$B" $S/BINARY.txt || stop "bin/ is not the P2-19 build"

say "ord 1: z0 (unit/profile_depth)"
$PY $U/profile_depth/test > $S/Z0.txt 2>&1 || { cat $S/Z0.txt >> $LOG; stop "z0 ($S/Z0.txt)"; }
grep "PASS  z0\|REPORT\|^PASS" $S/Z0.txt | cut -c1-160 | tee -a $LOG

say "ord 2: a0 tree gate"
bash $W/guarded_post_all_par.sh p2_19 > $S/tree.log 2>&1
[ $? = 2 ] && stop "the binary moved during the sweep"
$PY $TC p2_19 p2_18 d2aa d2t > $S/GATE_tree.txt 2>&1; cat $S/GATE_tree.txt >> $LOG
grep -q "vs p2_18: 0 of 2639 changed, 0 only in new, 0 only in ref" $S/GATE_tree.txt || stop "a0 vs p2_18"
grep -q "vs d2aa: 26 of 2639 changed" $S/GATE_tree.txt || stop "a0 control d2aa"
grep -q "vs d2t: 863 of 2639 changed" $S/GATE_tree.txt || stop "a0 control d2t"
say "a0 PASS; a1 face gate"
$PY $W/face_gate.py run p2_19 > $S/face.log 2>&1 || stop "face_gate run ($S/face.log)"
$PY $W/face_gate.py cmp p2_18 p2_19 > $S/GATE_face_cmp.txt 2>&1
$PY $S/reduce.py A1 > /dev/null; grep "^  \|^a1" $S/A1.txt | cut -c1-200 | tee -a $LOG
grep -q "^a1 PASS" $S/A1.txt || stop "a1 ($S/A1.txt)"
say "a2: T00 at three meshes"
$PY $S/run.py T00_2 T00_1 T00_05 >> $S/RUN.log 2>&1
$PY $S/reduce.py A2 > /dev/null; grep "^  \|^a2" $S/A2.txt | tee -a $LOG
grep -q "^a2 PASS" $S/A2.txt || stop "a2 ($S/A2.txt)"

say "ord 3: a3 CAL"
$PY $S/run.py CAL >> $S/RUN.log 2>&1
$PY $S/reduce.py A3 > /dev/null; cat $S/A3.txt | tee -a $LOG
grep -q "^a3 PASS" $S/A3.txt || stop "a3 ($S/A3.txt)"
say "HOS, C0H150 (reports)"
$PY $S/run.py HOS C0H150 >> $S/RUN.log 2>&1

say "ord 4: the locks on the shipped binary, RESULTS.txt"
for t in lateral_flux_lock buckle_lock retired_keys; do
    $PY $U/$t/test > $S/lock_$t.log 2>&1; say "  $t rc=$?"
done
$PY $U/face_removal/test --check > $S/lock_face_removal.log 2>&1
say "  face_removal --check rc=$? ($(grep -c '  FAIL' $S/lock_face_removal.log) FAIL lines; P2-18 recorded P2-9's 4 tilt-rate FAILs)"
$PY $S/reduce.py > $S/RESULTS.txt 2>&1
sed -n '/^e0 VERDICT/,$p' $S/RESULTS.txt | tee -a $LOG
say "CHAIN DONE"
