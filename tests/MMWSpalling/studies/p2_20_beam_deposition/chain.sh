#!/bin/bash
# P2-20 in the packet's run order; STOPS at the first failed gate. Log: CHAIN.log.
#   ord 1 z0 (report), z1 retired_keys (+ --pre src_pre --red R5), z2 unit/beam, z3 CAL_S/HOS_S + check_shared
#   ord 2 a0 tree gate vs p2_18, a1 face gate cmp p2_19 p2_20, a2 beam_void_closure (on the sweep's own output)
#   ord 3 Z1 (b1, b2), Z1x (b3) gates; Z1w, Z2, Z3a-c reports | ord 4 RESULTS.txt (e0)
set -u
ROOT=/Users/tzetze20/amr_tools/alamo; cd $ROOT
S=tests/MMWSpalling/studies/p2_20_beam_deposition
U=tests/MMWSpalling/unit
W=tests/MMWSpalling/studies/tree_sweep_0922/witness
TC=tests/MMWSpalling/studies/p2_9_face_removal/tree_compare.py
ZT=tests/MMWSpalling/validation/zhang/zhang_oglesby/test
CHK=tests/MMWSpalling/validation/shared/check_shared.py
BIN=bin/mmwspalling-3d-g++
PY=/Users/tzetze20/Desktop/code/.venv/bin/python; LOG=$S/CHAIN.log
say()  { echo "$(date +%T) $*" | tee -a $LOG; }
stop() { say "STOP: $*"; exit 1; }
caffeinate -dis -w $$ &
B=$(shasum -a 256 $BIN | cut -c1-64)
say "P2-20 chain (binary ${B:0:16})"
grep -q "$B" $S/BINARY.txt || stop "bin/ is not the P2-20 build"

say "ord 1: z0 (offline)"
$PY $S/z0.py > /dev/null && grep "BIT-EXACT\|NEAREST" $S/Z0.txt | cut -c1-120 | tee -a $LOG
say "z1 unit/retired_keys"
$PY $U/retired_keys/test --pre $S/src_pre --red R5 > $S/Z1_retired_keys.log 2>&1 || { tail -20 $S/Z1_retired_keys.log >> $LOG; stop "z1"; }
tail -1 $S/Z1_retired_keys.log | cut -c1-120 | tee -a $LOG
say "z2 unit/beam (its output/ archived first, SR 9)"
[ -d $U/beam/output ] && [ ! -d $S/unit_beam_output_pre ] && cp -cR $U/beam/output $S/unit_beam_output_pre
rm -rf $U/beam/output
mpirun --oversubscribe --bind-to none -np 4 $BIN $U/beam/input > $S/Z2_beam_run.log 2>&1 || stop "z2 run"
$PY $U/beam/test $U/beam/output > $S/Z2_beam.log 2>&1 || { cat $S/Z2_beam.log >> $LOG; stop "z2"; }
tail -3 $S/Z2_beam.log | tee -a $LOG
say "z3 CAL_S, HOS_S + check_shared"
$PY $S/run.py CAL_S HOS_S >> $S/RUN.log 2>&1 || stop "z3 legs"
{ $PY $CHK tests/MMWSpalling/validation/kant/grimsel_calibration/input $S/legs/CAL_S/argv.sh --run $S/legs/CAL_S/run; c1=$?
  $PY $CHK tests/MMWSpalling/validation/hoeser/input_lab $S/legs/HOS_S/argv.sh --run $S/legs/HOS_S/run; c2=$?
  grep -h "^beam.eps_r\|^beam.tan_delta" $S/legs/CAL_S/run/metadata $S/legs/HOS_S/run/metadata
  [ $c1 = 0 ] && [ $c2 = 0 ] && echo "z3 PASS" || echo "z3 FAIL"; } > $S/Z3_check_shared.txt 2>&1
cat $S/Z3_check_shared.txt | tee -a $LOG
grep -q "^z3 PASS" $S/Z3_check_shared.txt || stop "z3"

say "ord 2: a0 tree gate"
bash $W/guarded_post_all_par.sh p2_20 > $S/tree.log 2>&1
[ $? = 2 ] && stop "the binary moved during the sweep"
$PY $TC p2_20 p2_18 d2aa d2t > $S/GATE_tree.txt 2>&1; cat $S/GATE_tree.txt >> $LOG
grep -q "vs p2_18: 0 of 2639 changed, 0 only in new, 0 only in ref" $S/GATE_tree.txt || stop "a0 vs p2_18"
grep -q "vs d2aa: 26 of 2639 changed" $S/GATE_tree.txt || stop "a0 control d2aa"
grep -q "vs d2t: 863 of 2639 changed" $S/GATE_tree.txt || stop "a0 control d2t"
say "a0 PASS; a1 face gate"
$PY $W/face_gate.py run p2_20 > $S/face.log 2>&1 || stop "face_gate run ($S/face.log)"
$PY $W/face_gate.py cmp p2_19 p2_20 > $S/GATE_face_cmp.txt 2>&1
tail -3 $S/GATE_face_cmp.txt | tee -a $LOG
grep -q "VERDICT: IDENTICAL" $S/GATE_face_cmp.txt || stop "a1 ($S/GATE_face_cmp.txt)"
say "a2 beam_void_closure on the sweep's own output (this binary)"
$PY $U/beam_void_closure/test $U/beam_void_closure/output > $S/A2_beam_void_closure.log 2>&1
{ grep "^\[PASS\] T\|^\[FAIL\] T" $S/A2_beam_void_closure.log
  nbad=$(grep "^\[FAIL\] T" $S/A2_beam_void_closure.log | grep -vc "^\[FAIL\] T3 ")
  npass=$(grep -c "^\[PASS\] T[124]" $S/A2_beam_void_closure.log)
  if [ "$nbad" = 0 ] && [ "$npass" -ge 7 ] && grep -q "^\[FAIL\] T3 " $S/A2_beam_void_closure.log
  then echo "a2 PASS (T1a-d, T2a-b, T4 PASS: $npass rows; T3 FAIL by design, unchanged)"; else echo "a2 FAIL ($nbad non-T3 FAIL, $npass PASS)"; fi; } > $S/A2_beam_void_closure.txt
cut -c1-160 $S/A2_beam_void_closure.txt | tee -a $LOG
grep -q "^a2 PASS" $S/A2_beam_void_closure.txt || stop "a2"

say "ord 3: Z1 (b1, b2)"
$PY $S/run.py Z1 >> $S/RUN.log 2>&1 || stop "Z1 run"
$PY $ZT $S/legs/Z1/run > $S/TEST_Z1.log 2>&1; rz=$?
grep "  PASS\|  FAIL\|^PASS\|^FAIL" $S/TEST_Z1.log | cut -c1-200 | tee -a $LOG
[ $rz = 0 ] || stop "b1/b2 ($S/TEST_Z1.log)"
say "Z1x (b3)"
$PY $S/run.py Z1x >> $S/RUN.log 2>&1 || stop "Z1x run"
$PY $ZT $S/legs/Z1/run --explicit $S/legs/Z1x/run > $S/TEST_Z1x.log 2>&1
grep "b3" $S/TEST_Z1x.log | tee -a $LOG
grep -q "PASS  b3" $S/TEST_Z1x.log || stop "b3 ($S/TEST_Z1x.log)"
say "reports: Z1w, Z2, Z3a, Z3b, Z3c"
$PY $S/run.py Z1w Z2 Z3a Z3b Z3c >> $S/RUN.log 2>&1
for l in Z1w Z2; do $PY $ZT $S/legs/$l/run > $S/TEST_$l.log 2>&1; done

say "ord 4: RESULTS.txt"
$PY $S/reduce.py > $S/RESULTS.txt 2>&1
sed -n '/^c1\/c2 --/,$p' $S/RESULTS.txt | cut -c1-240 | tee -a $LOG
say "CHAIN DONE"
