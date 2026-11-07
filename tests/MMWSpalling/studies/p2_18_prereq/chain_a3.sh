#!/bin/bash
# P2-18 A3 gates on the helper-move binary: tree gate (0 of 2639 vs p2_16), face gate (0 of 1172 vs pre), source locks.
ROOT=/Users/tzetze20/amr_tools/alamo; S=$ROOT/tests/MMWSpalling/studies/p2_18_prereq; W=$ROOT/tests/MMWSpalling/studies/tree_sweep_0922/witness
PY=/Users/tzetze20/Desktop/code/.venv/bin/python
cd $ROOT
shasum -a 256 bin/mmwspalling-3d-g++ > $S/BINARY_a3.txt
echo "== tree gate"; bash $W/guarded_post_all_par.sh p2_18_a3 > $S/tree_a3.log 2>&1; echo "guard rc $?"
echo "tree gate p2_18_a3 vs p2_16: $(diff <(awk '{print $2, $1}' $W/p2_16.md5 | sort) <(awk '{print $2, $1}' $W/p2_18_a3.md5 | sort) | grep -c '^[<>]') differing lines ($(wc -l < $W/p2_16.md5) / $(wc -l < $W/p2_18_a3.md5) files); controls: d2aa $(diff <(awk '{print $2, $1}' $W/d2aa.md5 | sort) <(awk '{print $2, $1}' $W/p2_18_a3.md5 | sort) | grep -c '^>') d2t $(diff <(awk '{print $2, $1}' $W/d2t.md5 | sort) <(awk '{print $2, $1}' $W/p2_18_a3.md5 | sort) | grep -c '^>')" | tee $S/GATE_tree_a3.txt
echo "== face gate"; $PY $W/face_gate.py run p2_18_a3 > $S/face_a3.log 2>&1; echo "face run rc $?"
$PY $W/face_gate.py cmp p2_18_pre p2_18_a3 | tee $S/GATE_face_a3.txt
echo "== locks"
for t in buckle_lock retired_keys lateral_flux_lock; do echo "--- $t"; $PY tests/MMWSpalling/unit/$t/test > $S/lock_${t}_a3.log 2>&1; echo "$t rc $?"; grep -c "FAIL" $S/lock_${t}_a3.log; done
echo "--- face_removal --check"; $PY tests/MMWSpalling/unit/face_removal/test --check > $S/lock_face_removal_a3.log 2>&1; echo "face_removal rc $?"; grep -n "a3 i\|FAIL" $S/lock_face_removal_a3.log | head -12
echo "== chain done"
