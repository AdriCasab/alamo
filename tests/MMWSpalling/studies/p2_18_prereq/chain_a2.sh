#!/bin/bash
# P2-18 A2 gates on the final binary: tree gate (0 of 2639 vs p2_16), face gate (0 of 1172 vs pre), source locks,
# the retired-at-default legs (K1-K6), every kept deck still starts (stop_time 1e-9), rank check on the final binary.
ROOT=/Users/tzetze20/amr_tools/alamo; S=$ROOT/tests/MMWSpalling/studies/p2_18_prereq; W=$ROOT/tests/MMWSpalling/studies/tree_sweep_0922/witness
T=$ROOT/tests/MMWSpalling; PY=/Users/tzetze20/Desktop/code/.venv/bin/python; BIN=$ROOT/bin/mmwspalling-3d-g++
SCR=/private/tmp/claude-501/-Users-tzetze20-amr-tools-alamo/79a0a434-469e-4979-a234-ccc76ca8c68d/scratchpad/p2_18
mkdir -p $SCR; cd $ROOT
shasum -a 256 $BIN > $S/BINARY_a2.txt
echo "== tree gate"; bash $W/guarded_post_all_par.sh p2_18 > $S/tree_a2.log 2>&1; echo "guard rc $?"
echo "tree gate p2_18 vs p2_16: $(diff <(awk '{print $2, $1}' $W/p2_16.md5 | sort) <(awk '{print $2, $1}' $W/p2_18.md5 | sort) | grep -c '^[<>]') differing lines ($(wc -l < $W/p2_16.md5) / $(wc -l < $W/p2_18.md5) files); controls: d2aa $(diff <(awk '{print $2, $1}' $W/d2aa.md5 | sort) <(awk '{print $2, $1}' $W/p2_18.md5 | sort) | grep -c '^>') d2t $(diff <(awk '{print $2, $1}' $W/d2t.md5 | sort) <(awk '{print $2, $1}' $W/p2_18.md5 | sort) | grep -c '^>')" | tee $S/GATE_tree_a2.txt
echo "== face gate"; $PY $W/face_gate.py run p2_18 > $S/face_a2.log 2>&1; echo "face run rc $?"
$PY $W/face_gate.py cmp p2_18_pre p2_18 | tee $S/GATE_face_a2.txt
echo "== locks"
for t in buckle_lock retired_keys lateral_flux_lock; do echo "--- $t"; $PY $T/unit/$t/test > $S/lock_${t}_a2.log 2>&1; echo "$t rc $?"; grep -c "FAIL" $S/lock_${t}_a2.log; done
echo "--- face_removal --check (a3 i block)"; $PY $T/unit/face_removal/test --check > $S/lock_face_removal_a2.log 2>&1; echo "rc $? (4 tilt-rate FAILs are P2-9's record, pre-existing)"; sed -n 2,12p $S/lock_face_removal_a2.log | grep -c PASS
echo "== retired-at-default legs (grimsel_calibration deck, stop_time 1e-9)"
leg() { n=$1; shift; mpirun --oversubscribe --bind-to none -np 4 $BIN $T/validation/kant/grimsel_calibration/input stop_time=1e-9 plot_file=$SCR/$n "$@" > $S/K_$n.log 2>&1; rc=$?; echo "K $n rc=$rc  $(grep -o 'surface_patch\.[a-z_A-Z]* = [^ ]* was retired at its default\|specified but not used[^$]*' $S/K_$n.log | head -1)  $(grep -c 'retired at its default' $S/K_$n.log) msg"; }
leg K1_face_default surface_patch.face=zhi
leg K2_face_other   surface_patch.face=xlo
leg K3_Text_any     surface_patch.T_ext=293.15
leg K4_type_other   surface_patch.type=foo
leg K5_atten_toplevel microstructure.phase0.alpha_attenuation=2000.0
leg K6_grimsel_Tref microstructure.phase0.grimsel.T_ref=293.15
leg K0_plain
echo "== every kept deck starts (stop_time 1e-9)"; $PY $S/runnable.py | tee $S/RUNNABLE_a2.txt | tail -3
echo "== rank check on the final binary"; $PY $W/face_gate.py rankcheck p2_18 > $S/face_a2_rankcheck.log 2>&1; tail -1 $S/face_a2_rankcheck.log
echo "== chain done"
