#!/bin/bash
# P2-12: the no-change gates of the path-mode deletion (gates.py docstring). Stops at the first failed gate.
# Log: CHAIN.log. Holds the host awake for its own lifetime.
set -u
ROOT=/Users/tzetze20/amr_tools/alamo; cd $ROOT
S=tests/MMWSpalling/studies/p2_12_retire_path; W=tests/MMWSpalling/studies/tree_sweep_0922/witness
PY=/Users/tzetze20/Desktop/code/.venv/bin/python; LOG=$S/CHAIN.log
TC=tests/MMWSpalling/studies/p2_9_face_removal/tree_compare.py
PRE=tests/MMWSpalling/studies/p2_11b_path/src_pre
say()  { echo "$(date +%T) $*" | tee -a $LOG; }
stop() { say "STOP: $*"; exit 1; }
caffeinate -dis -w $$ &
say "P2-12 gates (binary $(shasum -a 256 bin/mmwspalling-3d-g++ | cut -c1-16))"
cmp -s $PRE/MMWSpalling.H src/Integrator/MMWSpalling.H && cmp -s $PRE/FaceRemoval.H src/Integrator/MMWSpalling/FaceRemoval.H \
  || stop "d1 source differs from src_pre"
say "d1 PASS (both headers byte-identical to studies/p2_11b_path/src_pre)"
# d1b: the BINARY carries none of the retired code (a cp -p restore once left make with a stale object:
# stale_binary/ is that run's record). Strings unique to P2-11b's source.
for w in face_heat barrier_expr "heat along the surface path" fl_n_unreached; do
  n=$(strings bin/mmwspalling-3d-g++ | grep -c "$w"); [ "$n" = 0 ] || stop "d1b binary still carries '$w' ($n)"
done
say "d1b PASS (the binary carries no P2-11b string)"
bash $W/guarded_post_all_par.sh p2_12 > $S/tree.log 2>&1
[ $? = 2 ] && stop "the binary moved during the sweep"
$PY $TC p2_12 p2_10_c2 p2_11b d2aa d2t | tee $S/GATE_tree.txt | tee -a $LOG
grep -q "vs p2_10_c2: 0 of 2639 changed, 0 only in new, 0 only in ref" $S/GATE_tree.txt || stop "d2 tree vs p2_10_c2"
grep -q "vs p2_11b: 0 of 2639 changed, 0 only in new, 0 only in ref" $S/GATE_tree.txt || stop "d2 tree vs p2_11b"
grep -q "vs d2aa: 26 of 2639 changed" $S/GATE_tree.txt || stop "d2 control d2aa"
grep -q "vs d2t: 863 of 2639 changed" $S/GATE_tree.txt || stop "d2 control d2t"
say "d2 PASS"
$PY $S/gates.py RETIRED > $S/retired.log 2>&1 || stop "d5 ($S/RETIRED.txt)"
say "d5 PASS"
$PY $S/gates.py W > $S/W.log 2>&1 || stop "d3 ($S/W_p2_12_vs_c2.txt)"
say "d3 PASS"
$PY $S/gates.py UNITS > $S/units.log 2>&1 || stop "d4 ($S/UNITS_p2_12_vs_c2.txt)"
say "d4 PASS"
say "ALL GATES PASS"
