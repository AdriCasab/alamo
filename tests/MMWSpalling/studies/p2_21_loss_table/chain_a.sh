#!/bin/bash
# P2-21 part A, the packet's run order 1-2 (the gates the Fig. 34 power reading cannot affect); STOPS at the
# first failed gate. Log: CHAIN.log.   Part B (chain_b.sh): the legs, ord 3-4.
#   ord 1  z1-z3 unit/beam_tau (+ --bin_pre bin_pre: the P2-20 binary aborts on the key), z4 unit/retired_keys
#   ord 2  a0 tree gate vs p2_18, a1 face gate cmp p2_20 p2_21, a2 beam_void_closure (the sweep's own output) +
#          unit/beam (its output archived first), a3 C1 vs P2-20's Z1 byte for byte (then C1 pruned)
set -u
ROOT=/Users/tzetze20/amr_tools/alamo; cd $ROOT
S=tests/MMWSpalling/studies/p2_21_loss_table
U=tests/MMWSpalling/unit
W=tests/MMWSpalling/studies/tree_sweep_0922/witness
TC=tests/MMWSpalling/studies/p2_9_face_removal/tree_compare.py
BIN=bin/mmwspalling-3d-g++
PY=/Users/tzetze20/Desktop/code/.venv/bin/python; LOG=$S/CHAIN.log
say()  { echo "$(date +%T) $*" | tee -a $LOG; }
stop() { say "STOP: $*"; exit 1; }
caffeinate -dis -w $$ &
B=$(shasum -a 256 $BIN | cut -c1-64)
say "P2-21 chain A (binary ${B:0:16})"
grep -q "$B" $S/BINARY.txt || stop "bin/ is not the P2-21 build"

say "ord 1: z1-z3 unit/beam_tau (U1 + bin_pre abort, U2, U3)"
$PY $U/beam_tau/test --bin_pre $S/bin_pre/mmwspalling-3d-g++ > $S/Z1_beam_tau.log 2>&1 || { tail -20 $S/Z1_beam_tau.log >> $LOG; stop "z1-z3"; }
grep "PASS  \|FAIL  " $S/Z1_beam_tau.log | cut -c1-170 | tee -a $LOG
say "z4 unit/retired_keys"
$PY $U/retired_keys/test > $S/Z4_retired_keys.log 2>&1 || { tail -30 $S/Z4_retired_keys.log >> $LOG; stop "z4"; }
grep "R6\|reassociate\|^PASS" $S/Z4_retired_keys.log | cut -c1-170 | tee -a $LOG

say "ord 2: a0 tree gate"
bash $W/guarded_post_all_par.sh p2_21 > $S/tree.log 2>&1
[ $? = 2 ] && stop "the binary moved during the sweep"
$PY $TC p2_21 p2_18 d2aa d2t > $S/GATE_tree.txt 2>&1; cat $S/GATE_tree.txt >> $LOG
grep -q "vs p2_18: 0 of 2639 changed, 0 only in new, 0 only in ref" $S/GATE_tree.txt || stop "a0 vs p2_18"
grep -q "vs d2aa: 26 of 2639 changed" $S/GATE_tree.txt || stop "a0 control d2aa"
grep -q "vs d2t: 863 of 2639 changed" $S/GATE_tree.txt || stop "a0 control d2t"
say "a0 PASS; a1 face gate"
$PY $W/face_gate.py run p2_21 > $S/face.log 2>&1 || stop "face_gate run ($S/face.log)"
$PY $W/face_gate.py cmp p2_20 p2_21 > $S/GATE_face_cmp.txt 2>&1
tail -3 $S/GATE_face_cmp.txt | tee -a $LOG
grep -q "VERDICT: IDENTICAL" $S/GATE_face_cmp.txt || stop "a1 ($S/GATE_face_cmp.txt)"
say "a2 beam_void_closure on the sweep's own output (this binary) + unit/beam"
$PY $U/beam_void_closure/test $U/beam_void_closure/output > $S/A2_beam_void_closure.log 2>&1
{ grep "^\[PASS\] T\|^\[FAIL\] T" $S/A2_beam_void_closure.log
  nbad=$(grep "^\[FAIL\] T" $S/A2_beam_void_closure.log | grep -vc "^\[FAIL\] T3 ")
  npass=$(grep -c "^\[PASS\] T[124]" $S/A2_beam_void_closure.log)
  if [ "$nbad" = 0 ] && [ "$npass" -ge 7 ] && grep -q "^\[FAIL\] T3 " $S/A2_beam_void_closure.log
  then echo "a2 beam_void_closure PASS (T1a-d, T2a-b, T4 PASS: $npass rows; T3 FAIL by design, unchanged)"; else echo "a2 beam_void_closure FAIL ($nbad non-T3 FAIL, $npass PASS)"; fi; } > $S/A2_beam_void_closure.txt
cut -c1-160 $S/A2_beam_void_closure.txt | tee -a $LOG
grep -q "^a2 beam_void_closure PASS" $S/A2_beam_void_closure.txt || stop "a2 beam_void_closure"
[ -d $U/beam/output ] && [ ! -d $S/unit_beam_output_p2_20 ] && cp -cR $U/beam/output $S/unit_beam_output_p2_20
rm -rf $U/beam/output
mpirun --oversubscribe --bind-to none -np 4 $BIN $U/beam/input > $S/A2_beam_run.log 2>&1 || stop "a2 unit/beam run"
$PY $U/beam/test $U/beam/output > $S/A2_beam.log 2>&1 || { cat $S/A2_beam.log >> $LOG; stop "a2 unit/beam"; }
tail -2 $S/A2_beam.log | tee -a $LOG
say "a3 C1 vs P2-20's Z1"
$PY $S/run.py C1 >> $S/RUN.log 2>&1 || stop "C1 run"
$PY $S/a3_compare.py > $S/A3.txt 2>&1; cat $S/A3.txt | tee -a $LOG
grep -q "^a3 PASS" $S/A3.txt || stop "a3"
# C1 pruned after a3 (guardrail): keep the first and the last cell plotfile, drop the node plotfiles
( cd $S/legs/C1/run && ls -d *cell *node | grep -v -e '^00000cell$' -e '^00075cell$' | xargs rm -rf )
echo "C1 pruned after a3 PASSED: kept 00000cell 00075cell + thermo.dat (a3 compared all 16 cell plotfiles, 112 files)" > $S/legs/C1/PRUNED.txt
say "CHAIN A DONE (ord 1-2 PASS)"
