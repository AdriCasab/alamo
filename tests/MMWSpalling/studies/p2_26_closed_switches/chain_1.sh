#!/bin/bash
# P2-26 stage 1 (B1 + B2 + B3) gates, in the packet's run order; STOPS at a failed gate. Log: CHAIN_1.log.
#   g1  a1 tree gate p2_26a vs p2_18 (0/2639; controls d2aa 26, d2t 863); a2 face gate run p2_26a, cmp p2_24 p2_26a
#   k1  a3 retired_keys (+ --pre src_pre_1 --red R7), lateral_flux_lock, buckle_lock, face_removal --check, jet_d2c
#       (its sweep log: the RECF legs are gone); face_heat (R_pad / R_above restated) and removal_increment as the
#       per-site checks of WHICH_TESTS.md for the march / FaceHr and the column path; a4 probes.py
set -u
ROOT=/Users/tzetze20/amr_tools/alamo; cd $ROOT
S=tests/MMWSpalling/studies/p2_26_closed_switches
U=tests/MMWSpalling/unit
W=tests/MMWSpalling/studies/tree_sweep_0922/witness
TC=tests/MMWSpalling/studies/p2_9_face_removal/tree_compare.py
BIN=bin/mmwspalling-3d-g++
PY=/Users/tzetze20/Desktop/code/.venv/bin/python; LOG=$S/CHAIN_1.log
say()  { echo "$(date +%T) $*" | tee -a $LOG; }
stop() { say "STOP: $*"; exit 1; }
caffeinate -dis -w $$ &
B=$(shasum -a 256 $BIN | cut -c1-64)
say "P2-26 stage 1 chain (binary ${B:0:16})"
grep -q "$B" $S/BINARY_1.txt || stop "bin/ is not the stage-1 binary (BINARY_1.txt)"

FROM=${FROM:-g1}
if [ "$FROM" = "k1b" ]; then
say "resume at k1b: a1 / a2 / the locks before face_removal already PASSED on this binary"
else
if [ "$FROM" = "g1" ]; then
say "g1/a1: tree gate p2_26a"
bash $W/guarded_post_all_par.sh p2_26a > $S/tree_1.log 2>&1
[ $? = 2 ] && stop "the binary moved during the sweep"
grep "=== .* rc=" $S/tree_1.log | tee -a $LOG
$PY $TC p2_26a p2_18 d2aa d2t > $S/GATE_tree_1.txt 2>&1; cat $S/GATE_tree_1.txt >> $LOG
grep -q "vs p2_18: 0 of 2639 changed, 0 only in new, 0 only in ref" $S/GATE_tree_1.txt || stop "a1 vs p2_18"
grep -q "vs d2aa: 26 of 2639 changed" $S/GATE_tree_1.txt || stop "a1 control d2aa"
grep -q "vs d2t: 863 of 2639 changed" $S/GATE_tree_1.txt || stop "a1 control d2t"
say "a1 PASS; a2 face gate"
$PY $W/face_gate.py run p2_26a > $S/face_1.log 2>&1 || stop "face_gate run ($S/face_1.log)"
$PY $W/face_gate.py cmp p2_24 p2_26a > $S/GATE_face_1.txt 2>&1
tail -3 $S/GATE_face_1.txt | cut -c1-220 | tee -a $LOG
grep -q "VERDICT: IDENTICAL" $S/GATE_face_1.txt || stop "a2 is not IDENTICAL ($S/GATE_face_1.txt)"
else
say "resume at k1 (FROM=$FROM): a1 / a2 already PASSED on this binary (GATE_tree_1.txt, GATE_face_1.txt)"
grep -q "vs p2_18: 0 of 2639 changed" $S/GATE_tree_1.txt && grep -q "VERDICT: IDENTICAL" $S/GATE_face_1.txt \
  || stop "resume: the gate files do not show a1 / a2 PASS"
fi

say "k1/a3: locks"
rm -rf $S/pre_flat_1; mkdir -p $S/pre_flat_1
cp $S/src_pre_1/src/Integrator/MMWSpalling.H $S/src_pre_1/src/Integrator/MMWSpalling/Removal.H \
   $S/src_pre_1/src/Integrator/MMWSpalling/FaceRemoval.H $S/src_pre_1/src/Numeric/MMWBeam.H $S/pre_flat_1/
for t in retired_keys lateral_flux_lock buckle_lock; do
  $PY $U/$t/test > $S/K1_$t.log 2>&1 || { tail -15 $S/K1_$t.log >> $LOG; stop "a3 $t"; }
  echo "  $t: $(tail -1 $S/K1_$t.log | cut -c1-140)" | tee -a $LOG
done
$PY $U/retired_keys/test --source --pre $S/pre_flat_1 --red R7 > $S/K1_retired_keys_pre.log 2>&1 \
  || { tail -15 $S/K1_retired_keys_pre.log >> $LOG; stop "a3 retired_keys --pre (R7 must be RED on src_pre_1)"; }
echo "  retired_keys --pre src_pre_1 --red R7: $(tail -1 $S/K1_retired_keys_pre.log | cut -c1-60)" | tee -a $LOG
fi
# face_removal --check re-reduces unit/face_removal/output (written 2026-10-27, P2-9's binary): it reads FAIL (4) --
# the tilted-plane rate rows TL20/TL35/TL45/TR -- in P2-18 and P2-19's lock logs too. The bar is "unchanged":
# the whole log byte-identical to studies/p2_19_flake_depth/lock_face_removal.log.
$PY $U/face_removal/test --check > $S/K1_face_removal.log 2>&1
diff -q $S/K1_face_removal.log tests/MMWSpalling/studies/p2_19_flake_depth/lock_face_removal.log > /dev/null \
  || { tail -15 $S/K1_face_removal.log >> $LOG; stop "a3 face_removal --check differs from P2-19's recorded log"; }
echo "  face_removal --check: $(tail -1 $S/K1_face_removal.log) -- byte-identical to P2-19's log (the standing state)" | tee -a $LOG
cp $W/p2_26a_wit/jet_d2c.log $S/K1_jet_d2c.log
grep -q "^PASS" $S/K1_jet_d2c.log || { tail -15 $S/K1_jet_d2c.log >> $LOG; stop "a3 jet_d2c (the sweep's run)"; }
echo "  jet_d2c (the sweep's run, p2_26a_wit/jet_d2c.log): $(tail -1 $S/K1_jet_d2c.log)" | tee -a $LOG
for t in face_heat removal_increment; do
  $PY $U/$t/test > $S/K1_$t.log 2>&1; rc=$?
  echo "  $t: rc $rc, $(tail -1 $S/K1_$t.log | cut -c1-140)" | tee -a $LOG
done
$PY $S/probes.py > /dev/null 2>&1 || stop "a4 probes (PROBES.txt)"
tail -1 $S/PROBES.txt | tee -a $LOG
say "STAGE 1 DONE"
