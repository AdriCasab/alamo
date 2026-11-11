#!/bin/bash
# P2-27 -- the gate chain, one stage per call, in the packet's run order; STOPS at a failed gate. Log: CHAIN.log.
#   usage: chain.sh <stage 1..5>          (FROM=k<n> resumes after the gates of stage n; they are hash-memoised anyway)
# Every stage: checks bin/ is BINARY_<n>.txt's binary, runs its gates and tests, then prints the macOS sleep/wake events
# inside its own window beside its wall time (pmset -g log): P2-26's "slow" stages were clamshell sleep on battery,
# not slow code, and a slow stage must never again be read as slow code.
set -u
ROOT=/Users/tzetze20/amr_tools/alamo; cd $ROOT
S=tests/MMWSpalling/studies/p2_27_cleanup
U=tests/MMWSpalling/unit
W=tests/MMWSpalling/studies/tree_sweep_0922/witness
TC=tests/MMWSpalling/studies/p2_9_face_removal/tree_compare.py
BIN=bin/mmwspalling-3d-g++
PY=/Users/tzetze20/Desktop/code/.venv/bin/python; LOG=$S/CHAIN.log
STAGE=${1:?usage: chain.sh <stage>}
say()  { echo "$(date +%T) [s$STAGE] $*" | tee -a $LOG; }
T0=$(date -u '+%Y-%m-%d %H:%M:%S'); S0=$SECONDS          # UTC: pmset -g log stamps +0000
sleeps() {
  local n
  n=$(pmset -g log | awk -v s="$T0" '($4=="Sleep"||$4=="Wake"||$4=="DarkWake") && substr($0,1,19) >= s' \
      | tee $S/SLEEP_$STAGE.txt | wc -l | tr -d ' ')
  say "wall $((SECONDS - S0)) s since $T0; sleep/wake events in this window: $n (SLEEP_$STAGE.txt)$( [ "$n" != 0 ] && echo ' -- READ THEM before reading the wall time')"
  pmset -g batt | grep -q "AC Power" || say "NOTE: on BATTERY now"
}
stop() { say "STOP: $*"; sleeps; exit 1; }
caffeinate -dis -w $$ &
B=$(shasum -a 256 $BIN | cut -c1-64)
say "P2-27 stage $STAGE chain (binary ${B:0:16}); $(pmset -g batt | head -1)"
grep -q "$B" $S/BINARY_$STAGE.txt || stop "bin/ is not the stage-$STAGE binary (BINARY_$STAGE.txt)"
FROM=${FROM:-g$STAGE}

# the tree gate (hash-memoised by binary sha; no --force) + a direct manifest diff against the previous label
tree_gate() {   # $1 label  $2 previous label
  bash $W/guarded_post_all_par.sh $1 > $S/tree_$STAGE.log 2>&1
  [ $? = 2 ] && stop "the binary moved during the sweep"
  grep "=== .* rc=\|GUARD: gate already" $S/tree_$STAGE.log | tee -a $LOG
  $PY $TC $1 $2 d2aa d2t > $S/GATE_tree_$STAGE.txt 2>&1
  cat $S/GATE_tree_$STAGE.txt | tee -a $LOG
}
face_gate() {   # $1 label  $2 reference label
  $PY $W/face_gate.py run $1 > $S/face_$STAGE.log 2>&1 || stop "face_gate run ($S/face_$STAGE.log)"
  $PY $W/face_gate.py cmp $2 $1 > $S/GATE_face_$STAGE.txt 2>&1
  tail -3 $S/GATE_face_$STAGE.txt | cut -c1-220 | tee -a $LOG
}
dt_messages() { # the dt check's Util::Message must appear in every face-gate leg's log, with no abort and no warning
  local n m w
  n=$(find $W/face_$1 -name run.log | wc -l | tr -d ' ')
  m=$(grep -l "explicit time-step limit, level 0" $(find $W/face_$1 -name run.log) | wc -l | tr -d ' ')
  w=$(grep -l "exceeds the explicit limit" $(find $W/face_$1 -name run.log) /dev/null | wc -l | tr -d ' ')
  say "dt check: message in $m of $n face-gate leg logs, over-limit lines in $w"
  [ "$m" = "$n" ] && [ "$w" = 0 ] || stop "the dt message is missing from a face leg, or a face leg is over the limit"
  m=$(grep -l "explicit time-step limit, level 0" $U/jet_d2c/output/*.log | wc -l | tr -d ' ')
  n=$(ls $U/jet_d2c/output/*.log | wc -l | tr -d ' ')
  say "dt check: message in $m of $n jet_d2c leg logs (abort legs stop before the first step)"
}
lock() {        # $1 test name  $2 log tag ; extra args after
  local t=$1 tag=$2; shift 2
  $PY $U/$t/test "$@" > $S/${tag}_$t.log 2>&1 || { tail -15 $S/${tag}_$t.log >> $LOG; stop "$tag $t"; }
  echo "  $t $*: $(tail -1 $S/${tag}_$t.log | cut -c1-150)" | tee -a $LOG
}

case $STAGE in
1)
  if [ "$FROM" = "g1" ]; then
    say "g1/a1: tree gate p2_27a vs p2_26 (jet_d2c enters the oracle: only-in-new = MOVERS_1_jet_d2c.txt)"
    tree_gate p2_27a p2_26
    $PY $S/tree_bar.py p2_27a p2_26 --only-new $S/MOVERS_1_jet_d2c.txt > $S/GATE_a1.txt 2>&1; rc=$?
    cat $S/GATE_a1.txt | tee -a $LOG
    [ $rc = 0 ] || stop "a1 ($S/GATE_a1.txt)"
    say "a1 PASS; a2 face gate p2_27a vs p2_24 (M0S now replays its record with the two C2 tokens)"
    face_gate p2_27a p2_24
    grep -q "VERDICT: IDENTICAL" $S/GATE_face_$STAGE.txt || stop "a2 is not IDENTICAL ($S/GATE_face_$STAGE.txt)"
    grep -c "spall.flake_coherence_length=0.0" $W/face_p2_27a/M0S/argv.sh | sed 's/^/  M0S argv carries spall.flake_coherence_length=0.0: /' | tee -a $LOG
    dt_messages p2_27a
  fi
  say "k1/a3: tests"
  lock interlocks K1 --pre $S/bin_pre/mmwspalling-3d-g++
  lock beam_schedule K1
  lock surface_T_recon K1
  lock retired_keys K1
  lock lateral_flux_lock K1
  lock buckle_lock K1
  $PY $U/face_removal/test --check > $S/K1_face_removal.log 2>&1
  diff -q $S/K1_face_removal.log tests/MMWSpalling/studies/p2_19_flake_depth/lock_face_removal.log > /dev/null \
    || { tail -15 $S/K1_face_removal.log >> $LOG; stop "a3 face_removal --check differs from P2-19's recorded log"; }
  echo "  face_removal --check: $(tail -1 $S/K1_face_removal.log) -- byte-identical to P2-19's log (the standing state)" | tee -a $LOG
  cp $W/p2_27a_wit/jet_d2c.log $S/K1_jet_d2c.log 2>/dev/null && \
    echo "  jet_d2c (the sweep's run): $(tail -1 $S/K1_jet_d2c.log | cut -c1-120)" | tee -a $LOG
  ;;
2)
  if [ "$FROM" = "g2" ]; then
    say "g2/b1: tree gate p2_27b vs p2_27a (spall_event's 8 paths leave: MOVERS_2_only_ref.txt)"
    tree_gate p2_27b p2_27a
    $PY $S/tree_bar.py p2_27b p2_27a --only-ref $S/MOVERS_2_only_ref.txt > $S/GATE_b1.txt 2>&1; rc=$?
    cat $S/GATE_b1.txt | tee -a $LOG
    [ $rc = 0 ] || stop "b1 ($S/GATE_b1.txt)"
    say "b1 PASS; b2 face gate p2_27b vs p2_24"
    face_gate p2_27b p2_24
    grep -q "VERDICT: IDENTICAL" $S/GATE_face_$STAGE.txt || stop "b2 is not IDENTICAL ($S/GATE_face_$STAGE.txt)"
    dt_messages p2_27b
  fi
  say "k2/b3: tests"
  rm -rf $S/pre_flat_2; mkdir -p $S/pre_flat_2
  cp $S/src_pre_2/src/Integrator/MMWSpalling.H $S/src_pre_2/src/Integrator/MMWSpalling/Removal.H \
     $S/src_pre_2/src/Integrator/MMWSpalling/FaceRemoval.H $S/pre_flat_2/
  lock retired_keys K2
  $PY $U/retired_keys/test --source --pre $S/pre_flat_2 --red R9 > $S/K2_retired_keys_pre.log 2>&1 \
    || { tail -15 $S/K2_retired_keys_pre.log >> $LOG; stop "b3 retired_keys --pre (R9 must be RED on src_pre_2)"; }
  echo "  retired_keys --pre src_pre_2 --red R9: $(tail -1 $S/K2_retired_keys_pre.log | cut -c1-60)" | tee -a $LOG
  lock lateral_flux_lock K2
  lock buckle_lock K2
  lock interlocks K2 --pre $S/bin_pre/mmwspalling-3d-g++
  O=$S/regrid_after; rm -rf $O; mkdir -p $O
  mpirun --oversubscribe --bind-to none -np 4 $BIN $U/amr_microstructure_regrid/input plot_file=$O/run > $O/run.log 2>&1
  echo "  amr_microstructure_regrid deck: rc $?" | tee -a $LOG
  $PY $U/amr_microstructure_regrid/test $O/run > $S/K2_amr_regrid_after.log 2>&1
  echo "  amr_microstructure_regrid (after; before = K2_amr_regrid_before.log rc 1): rc $?, $(grep -c FAIL: $S/K2_amr_regrid_after.log) FAIL line(s): $(grep FAIL: $S/K2_amr_regrid_after.log | cut -c1-60 | tr '\n' ';')" | tee -a $LOG
  for t in face_heat removal_increment; do
    $PY $U/$t/test > $S/K2_$t.log 2>&1 || { tail -15 $S/K2_$t.log >> $LOG; stop "b3 $t"; }
    echo "  $t: $(tail -1 $S/K2_$t.log | cut -c1-140)" | tee -a $LOG
  done
  cp $W/p2_27b_wit/scalar_flaw.log $S/K2_scalar_flaw.log 2>/dev/null && \
    echo "  scalar_flaw (the sweep's run): $(tail -1 $S/K2_scalar_flaw.log | cut -c1-120)" | tee -a $LOG
  say "m4: the Hu LEFM onset (Grimsel under Hu's granite-2 ramp), a record"
  $PY tests/MMWSpalling/validation/hu/hu_spall_onset/test_lefm > $S/M4_hu_lefm.log 2>&1; rc=$?
  tail -6 $S/M4_hu_lefm.log | tee -a $LOG
  [ $rc = 0 ] || say "m4: test_lefm rc $rc (a RECORD: reported, not a gate)"
  ;;
3)
  if [ "$FROM" = "g3" ]; then
    [ -d $S/tree_archive_p2_27b ] || stop "c1 needs tree_archive_p2_27b/ (audit_components.py archive p2_27b, BEFORE this sweep)"
    say "g3/c1: tree gate p2_27c vs p2_27b + the component audit (removed: D, kappa_eff, k_eff)"
    tree_gate p2_27c p2_27b
    $PY $S/audit_components.py tree p2_27b p2_27c > $S/AUDIT_3_tree.txt 2>&1; rc=$?
    cat $S/AUDIT_3_tree.txt | tee -a $LOG
    [ $rc = 0 ] || stop "c1 audit ($S/AUDIT_3_tree.txt)"
    say "c1 PASS; c2 face gate p2_27c + the component audit vs p2_27b, and cmp p2_24 p2_27c (expected DIFFERENT)"
    $PY $W/face_gate.py run p2_27c > $S/face_$STAGE.log 2>&1 || stop "face_gate run ($S/face_$STAGE.log)"
    $PY $W/face_gate.py cmp p2_24 p2_27c > $S/GATE_face_$STAGE.txt 2>&1
    tail -1 $S/GATE_face_$STAGE.txt | cut -c1-200 | tee -a $LOG
    $PY $S/audit_components.py face p2_27b p2_27c > $S/AUDIT_3_face.txt 2>&1; rc=$?
    cat $S/AUDIT_3_face.txt | tee -a $LOG
    [ $rc = 0 ] || stop "c2 audit ($S/AUDIT_3_face.txt)"
    $PY $W/face_gate.py rankcheck p2_27c > $S/RANKCHECK_3.txt 2>&1; rc=$?
    tail -1 $S/RANKCHECK_3.txt | cut -c1-200 | tee -a $LOG
    [ $rc = 0 ] || stop "c2 rankcheck ($S/RANKCHECK_3.txt)"
    dt_messages p2_27c
  fi
  say "k3: tests"
  rm -rf $S/pre_flat_3; mkdir -p $S/pre_flat_3
  cp $S/src_pre_3/src/Integrator/MMWSpalling.H $S/src_pre_3/src/Integrator/MMWSpalling/Removal.H \
     $S/src_pre_3/src/Integrator/MMWSpalling/FaceRemoval.H $S/pre_flat_3/
  lock retired_keys K3
  $PY $U/retired_keys/test --source --pre $S/pre_flat_3 --red R10 > $S/K3_retired_keys_pre.log 2>&1 \
    || { tail -15 $S/K3_retired_keys_pre.log >> $LOG; stop "k3 retired_keys --pre (R10 must be RED on src_pre_3)"; }
  echo "  retired_keys --pre src_pre_3 --red R10: $(tail -1 $S/K3_retired_keys_pre.log | cut -c1-60)" | tee -a $LOG
  lock lateral_flux_lock K3
  lock buckle_lock K3
  O=$S/regrid_after3; rm -rf $O; mkdir -p $O
  mpirun --oversubscribe --bind-to none -np 4 $BIN $U/amr_microstructure_regrid/input plot_file=$O/run > $O/run.log 2>&1
  $PY $U/amr_microstructure_regrid/test $O/run > $S/K3_amr_regrid.log 2>&1
  echo "  amr_microstructure_regrid: rc $?, $(grep -c FAIL: $S/K3_amr_regrid.log) FAIL line(s): $(grep FAIL: $S/K3_amr_regrid.log | cut -c1-60 | tr '\n' ';')" | tee -a $LOG
  ;;
*)
  stop "stage $STAGE is not written yet"
  ;;
esac
sleeps
say "STAGE $STAGE DONE"
