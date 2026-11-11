#!/bin/bash
# P2-26 stage 2 (C2 + C3) gates, in the packet's run order; STOPS at a failed gate. Log: CHAIN_2.log.
#   g2  b1 tree gate p2_26 vs p2_18 (the d2aa / d2t counts informational); classify.py over every differing path
#          (MOVERS_PRE.txt, pre2_out/); b2 face gate run p2_26, cmp p2_24 p2_26; rankcheck p2_26
#   k2  b3 retired_keys (all clauses, MUTATIONS, K legs; --pre src_pre_2 --red R8); lateral_flux_lock; buckle_lock;
#          face_removal --check (unchanged from P2-19's log)
#       b4 removal_increment (SN gone), face_heat, beam_void_closure (the sweep's run; T3 FAIL by design)
#   m   m1 line counts; m2/m3 ab.py (every edited out-of-gate deck before/after) + the three edited unit tests in full
# FROM=k2 resumes after the gates (they are memoised on the binary hash anyway).
set -u
ROOT=/Users/tzetze20/amr_tools/alamo; cd $ROOT
S=tests/MMWSpalling/studies/p2_26_closed_switches
U=tests/MMWSpalling/unit
W=tests/MMWSpalling/studies/tree_sweep_0922/witness
TC=tests/MMWSpalling/studies/p2_9_face_removal/tree_compare.py
BIN=bin/mmwspalling-3d-g++
PY=/Users/tzetze20/Desktop/code/.venv/bin/python; LOG=$S/CHAIN_2.log
say()  { echo "$(date +%T) $*" | tee -a $LOG; }
stop() { say "STOP: $*"; exit 1; }
caffeinate -dis -w $$ &
B=$(shasum -a 256 $BIN | cut -c1-64)
say "P2-26 stage 2 chain (binary ${B:0:16})"
grep -q "$B" $S/BINARY_2.txt || stop "bin/ is not the stage-2 binary (BINARY_2.txt)"
FROM=${FROM:-g2}

if [ "$FROM" = "g2" ]; then
say "g2/b1: tree gate p2_26 (pre2_out/ holds the two C2 legs' old bytes)"
bash $W/guarded_post_all_par.sh p2_26 > $S/tree_2.log 2>&1
[ $? = 2 ] && stop "the binary moved during the sweep"
grep "=== .* rc=" $S/tree_2.log | tee -a $LOG
$PY $TC p2_26 p2_18 d2aa d2t --list > $S/GATE_tree_2.txt 2>&1; cat $S/GATE_tree_2.txt >> $LOG
$PY $S/classify.py p2_26 p2_18 > $S/GATE_b1.txt 2>&1; rc=$?
grep "^P2-26 b1\|archive\|class counts\|moved per leg\|^b1 " $S/GATE_b1.txt | tee -a $LOG
[ $rc = 0 ] || stop "b1 ($S/GATE_b1.txt)"
say "b1 PASS"
fi
if [ "$FROM" = "g2" ] || [ "$FROM" = "b2" ]; then
say "b2 face gate"
$PY $W/face_gate.py run p2_26 > $S/face_2.log 2>&1 || stop "face_gate run ($S/face_2.log)"
$PY $W/face_gate.py cmp p2_24 p2_26 > $S/GATE_face_2.txt 2>&1
tail -3 $S/GATE_face_2.txt | cut -c1-220 | tee -a $LOG
grep -q "VERDICT: IDENTICAL" $S/GATE_face_2.txt || stop "b2 is not IDENTICAL ($S/GATE_face_2.txt)"
$PY $W/face_gate.py rankcheck p2_26 > $S/RANKCHECK_2.txt 2>&1; rc=$?
tail -3 $S/RANKCHECK_2.txt | cut -c1-200 | tee -a $LOG
[ $rc = 0 ] || stop "b4 rankcheck ($S/RANKCHECK_2.txt)"
fi

say "k2/b3: locks"
rm -rf $S/pre_flat_2; mkdir -p $S/pre_flat_2
cp $S/src_pre_2/src/Integrator/MMWSpalling.H $S/src_pre_2/src/Integrator/MMWSpalling/Removal.H \
   $S/src_pre_2/src/Integrator/MMWSpalling/FaceRemoval.H $S/src_pre_2/src/Numeric/MMWBeam.H $S/pre_flat_2/
for t in retired_keys lateral_flux_lock buckle_lock; do
  $PY $U/$t/test > $S/K2_$t.log 2>&1 || { tail -15 $S/K2_$t.log >> $LOG; stop "b3 $t"; }
  echo "  $t: $(grep -m1 '^PASS' $S/K2_$t.log | cut -c1-140)" | tee -a $LOG
done
$PY $U/retired_keys/test --source --pre $S/pre_flat_2 --red R8 > $S/K2_retired_keys_pre.log 2>&1 \
  || { tail -15 $S/K2_retired_keys_pre.log >> $LOG; stop "b3 retired_keys --pre (R8 must be RED on src_pre_2)"; }
echo "  retired_keys --pre src_pre_2 --red R8: $(grep -m1 '^PASS' $S/K2_retired_keys_pre.log | cut -c1-60)" | tee -a $LOG
$PY $U/face_removal/test --check > $S/K2_face_removal.log 2>&1
diff -q $S/K2_face_removal.log tests/MMWSpalling/studies/p2_19_flake_depth/lock_face_removal.log > /dev/null \
  || { tail -15 $S/K2_face_removal.log >> $LOG; stop "face_removal --check differs from P2-19's recorded log"; }
echo "  face_removal --check: $(tail -1 $S/K2_face_removal.log) -- byte-identical to P2-19's log (the standing state)" | tee -a $LOG

say "b4: removal_increment, face_heat, beam_void_closure"
for t in removal_increment face_heat; do
  $PY $U/$t/test > $S/K2_$t.log 2>&1 || { tail -15 $S/K2_$t.log >> $LOG; stop "b4 $t"; }
  echo "  $t: $(tail -1 $S/K2_$t.log | cut -c1-140)" | tee -a $LOG
done
cp $W/p2_26_wit/beam_void_closure.log $S/K2_beam_void_closure.log 2>/dev/null \
  || $PY $U/beam_void_closure/test > $S/K2_beam_void_closure.log 2>&1
grep "^\[\(PASS\|FAIL\)\]\|^  \(PASS\|FAIL\)\|T[1-4]" $S/K2_beam_void_closure.log | cut -c1-160 | tee -a $LOG

say "m1: line counts"
wc -l $S/src_pre_1/src/Integrator/MMWSpalling.H $S/src_pre_1/src/Integrator/MMWSpalling/*.H \
      $S/src_pre_1/src/Numeric/MMWBeam.H | tail -1 | sed 's/^/  before: /' | tee -a $LOG
wc -l src/Integrator/MMWSpalling.H src/Integrator/MMWSpalling/*.H src/Numeric/MMWBeam.H | tail -1 \
  | sed 's/^/  after:  /' | tee -a $LOG

say "m2/m3: the edited out-of-gate decks, before/after (ab.py), then their unit tests in full"
$PY $S/ab.py > $S/AB.txt 2>&1; rc=$?
grep "^[a-z].*horizon\|=>\|^m3" $S/AB.txt | cut -c1-200 | tee -a $LOG
for t in seed_void surface_T_recon gradient_blind; do
  $PY $U/$t/test > $S/M2_$t.log 2>&1; r=$?
  echo "  $t (new binary): rc $r, $(tail -1 $S/M2_$t.log | cut -c1-140)" | tee -a $LOG
done
[ $rc = 0 ] || stop "m3: an edited deck does not run on the new binary (AB.txt)"
say "STAGE 2 DONE"
