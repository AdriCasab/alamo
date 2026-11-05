#!/bin/bash
# P2-12 reversed: the restore gates (gates_restore.py docstring). Stops at the first failed gate. Log: CHAIN.log.
set -u
ROOT=/Users/tzetze20/amr_tools/alamo; cd $ROOT
S=tests/MMWSpalling/studies/p2_12_retire_path/restore; B=tests/MMWSpalling/studies/p2_12_retire_path/src_p2_11b
W=tests/MMWSpalling/studies/tree_sweep_0922/witness
PY=/Users/tzetze20/Desktop/code/.venv/bin/python; LOG=$S/CHAIN.log
TC=tests/MMWSpalling/studies/p2_9_face_removal/tree_compare.py
say()  { echo "$(date +%T) $*" | tee -a $LOG; }
stop() { say "STOP: $*"; exit 1; }
caffeinate -dis -w $$ &
say "P2-12 restore gates (binary $(shasum -a 256 bin/mmwspalling-3d-g++ | cut -c1-16))"
cmp -s $B/MMWSpalling.H src/Integrator/MMWSpalling.H && cmp -s $B/FaceRemoval.H src/Integrator/MMWSpalling/FaceRemoval.H \
  && cmp -s $B/lateral_flux_lock.test tests/MMWSpalling/unit/lateral_flux_lock/test || stop "r1 source differs from src_p2_11b"
[ -f tests/MMWSpalling/unit/face_heat/test ] || stop "r1 unit/face_heat not restored"
say "r1 PASS"
n=$(strings bin/mmwspalling-3d-g++ | grep -c face_heat); [ "$n" -gt 0 ] || stop "r1b binary has no face_heat ($n)"
say "r1b PASS (face_heat strings: $n)"
bash $W/guarded_post_all_par.sh p2_12r > $S/tree.log 2>&1
[ $? = 2 ] && stop "the binary moved during the sweep"
$PY $TC p2_12r p2_11b p2_10_c2 d2aa d2t | tee $S/GATE_tree.txt | tee -a $LOG
grep -q "vs p2_11b: 0 of 2639 changed, 0 only in new, 0 only in ref" $S/GATE_tree.txt || stop "r2 tree vs p2_11b"
grep -q "vs p2_10_c2: 0 of 2639 changed, 0 only in new, 0 only in ref" $S/GATE_tree.txt || stop "r2 tree vs p2_10_c2"
grep -q "vs d2aa: 26 of 2639 changed" $S/GATE_tree.txt || stop "r2 control d2aa"
grep -q "vs d2t: 863 of 2639 changed" $S/GATE_tree.txt || stop "r2 control d2t"
say "r2 PASS"
$PY $S/gates_restore.py X1REPRO > $S/x1repro.log 2>&1 || stop "r5 ($S/X1REPRO.txt)"
say "r5 PASS"
$PY $S/gates_restore.py W > $S/W.log 2>&1 || stop "r3 ($S/W_p2_12r_vs_c2.txt)"
say "r3 PASS"
$PY $S/gates_restore.py UNITS > $S/units.log 2>&1 || stop "r4 ($S/UNITS_p2_12r_vs_p2_11b.txt)"
say "r4 PASS"
say "ALL RESTORE GATES PASS"
