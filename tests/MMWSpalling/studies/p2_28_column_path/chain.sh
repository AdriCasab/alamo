#!/bin/bash
# P2-28 -- the gate chain, one stage per call, in the packet's run order. Log: CHAIN.log (one line per result).
#   usage: chain.sh 1        stage 1 (p2_28a): g1 (tree + face gates) then k1 (the stage-1 tests), then STOP for the user
# A failed GATE is logged and every later LONG leg (the old Kant site-A test with its 25-case sweep, the CAL leg) is
# skipped; the short tests still run and are reported. Every run is held awake (caffeinate -w on this shell).
set -u
ROOT=/Users/tzetze20/amr_tools/alamo; cd $ROOT
S=tests/MMWSpalling/studies/p2_28_column_path
U=tests/MMWSpalling/unit
V=tests/MMWSpalling/validation
W=tests/MMWSpalling/studies/tree_sweep_0922/witness
TC=tests/MMWSpalling/studies/p2_9_face_removal/tree_compare.py
BIN=bin/mmwspalling-3d-g++
PY=/Users/tzetze20/Desktop/code/.venv/bin/python; LOG=$S/CHAIN.log
STAGE=${1:?usage: chain.sh <stage>}
say()  { echo "$(date +%T) [s$STAGE] $*" | tee -a $LOG; }
T0=$(date -u '+%Y-%m-%d %H:%M:%S'); S0=$SECONDS
GATE_FAIL=0
sleeps() {
  local n
  n=$(pmset -g log | awk -v s="$T0" '($4=="Sleep"||$4=="Wake"||$4=="DarkWake") && substr($0,1,19) >= s' \
      | tee $S/SLEEP_$STAGE.txt | wc -l | tr -d ' ')
  say "wall $((SECONDS - S0)) s since $T0; sleep/wake events in this window: $n (SLEEP_$STAGE.txt)"
  pmset -g batt | grep -q "AC Power" || say "NOTE: on BATTERY now"
}
gate_fail() { GATE_FAIL=1; say "GATE FAILED: $*"; }
# one test: $1 tag (log name), $2 gate|report, then the command; logs $S/<tag>.log, one CHAIN line from its tail
t() {
  local tag=$1 kind=$2; shift 2
  "$@" > $S/$tag.log 2>&1; local rc=$?
  say "  $tag rc $rc: $(grep -E 'PASS|FAIL|RECORD' $S/$tag.log | tail -1 | cut -c1-170)"
  [ $rc != 0 ] && [ $kind = gate ] && gate_fail "$tag (rc $rc, $S/$tag.log)"
  return $rc
}
long_ok() { [ $GATE_FAIL = 0 ] || { say "  SKIPPED (a gate failed): $*"; return 1; }; }
caffeinate -dis -w $$ &
if [ $STAGE = 1 ]; then
  # r0: the first stage-1 chain ran on 8a1bfe93 (before the one-step-lag fix in RecedeAlongLines' prescribed T_s);
  # its results are kept, not overwritten (SR 9)
  if [ -f $S/GATE_a1.txt ] && [ ! -d $S/r0 ]; then
    mkdir -p $S/r0
    mv $S/GATE_a1.txt $S/GATE_a2.txt $S/tree_1.log $S/face_1.log $S/K1_*.log $S/K1b_*.log $S/r0/ 2>/dev/null
    cp $W/face_p2_28a.md5 $S/r0/face_p2_28a_r0.md5 2>/dev/null; cp $W/p2_28a.md5 $S/r0/p2_28a_r0.md5 2>/dev/null
    say "r0: the first chain's results (binary 8a1bfe93) moved to r0/"
  fi
  # the build (a no-op when bin/ is current); the stage-1 binary is whatever this leaves
  env EIGEN=$ROOT/ext CPLUS_INCLUDE_PATH=/opt/homebrew/include \
      LIBRARY_PATH=/opt/homebrew/lib:/opt/homebrew/Cellar/gcc/15.2.0_1/lib/gcc/current make -j8 > $S/build_1.log 2>&1 \
    || { say "STOP: the build failed ($S/build_1.log)"; exit 1; }
fi
B=$(shasum -a 256 $BIN | cut -c1-64)
say "P2-28 stage $STAGE chain (binary ${B:0:16}); $(pmset -g batt | head -1)"

if [ $STAGE = 1 ]; then
  # the stage-1 binary is the one on disk: every stage-1 source edit is in (src_pre_1/ holds the pre-image)
  mkdir -p $S/bin_1 && cp $BIN $S/bin_1/
  { echo "P2-28 stage 1 (p2_28a): the prescribed patch + evaluate-only on the face list; the Dirichlet factor 2"
    shasum -a 256 $BIN
    shasum -a 256 src/Integrator/MMWSpalling.H src/Integrator/MMWSpalling/*.H; } > $S/BINARY_1.txt
  # SR 9: the outputs these tests rmtree are the evidence a3 (ii) and m4 cite -- archive (move) once
  mkdir -p $S/archive
  for d in $V/kant/sp_kant_onset/output $V/kant/sp_kant_pressure_sweep/output $V/hu/hu_spall_onset/output_lefm; do
    n=$S/archive/$(echo $d | tr / _)_pre_p2_28
    if [ -d $d ] && [ ! -d $n ]; then mv $d $n && say "archived $d -> $n"; fi
  done

  # ---- smoke before campaign: every new deck one short run (rc 0, no abort); a failure stops here ----------------
  mkdir -p $S/smoke
  smoke() { local tag=$1; shift
    rm -rf $S/smoke/$tag; mkdir -p $S/smoke/$tag
    mpirun --oversubscribe --bind-to none -np 4 $BIN "$@" plot_file=$S/smoke/$tag/run > $S/smoke/$tag/run.log 2>&1
    local rc=$?
    say "  smoke $tag rc $rc $(grep -a 'ABORT' $S/smoke/$tag/run.log | head -1 | cut -c1-160)"
    [ $rc = 0 ] || { say "STOP: smoke $tag failed"; sleeps; exit 1; }
  }
  smoke GB_face $U/gradient_blind/input_face
  smoke FD_F $U/face_dirichlet/input stop_time=1.0
  smoke FD_C $U/face_dirichlet/input stop_time=1.0 surface_patch.face_list=0 surface_patch.robin_form=cell \
        surface.follow_mask=0 spall.surface_T=cell_average spall.removal_energy=cell
  smoke KP_face $V/kant/sp_kant_pressure_sweep/input_face stop_time=1.0 confining.p=48.0e6
  smoke KV_hfl $V/kant/sp_kant_onset/input_verification_hfl_face stop_time=1.0
  # (a separate site-A run at spall.surface_T = reconstructed segfaults at step 1, P2-27 binary too: not run)

  # ---- g1: the gates --------------------------------------------------------------------------------------------
  bash $W/guarded_post_all_par.sh p2_28a > $S/tree_1.log 2>&1
  [ $? = 2 ] && gate_fail "the binary moved during the tree sweep"
  $PY $TC p2_28a p2_27 d2aa d2t > $S/GATE_a1.txt 2>&1
  say "a1 tree p2_28a vs p2_27: $(grep -iE 'changed|only' $S/GATE_a1.txt | head -4 | tr '\n' ';' | cut -c1-220)"
  $PY $W/face_gate.py run p2_28a > $S/face_1.log 2>&1 || gate_fail "face_gate run (face_1.log)"
  $PY $W/face_gate.py cmp p2_27 p2_28a > $S/GATE_a2.txt 2>&1
  say "a2 face p2_27 -> p2_28a: $(tail -4 $S/GATE_a2.txt | tr '\n' ';' | cut -c1-260)"
  n=$(find $W/face_p2_28a -name run.log | wc -l | tr -d ' ')
  m=$(grep -l "explicit time-step limit, level 0" $(find $W/face_p2_28a -name run.log) | wc -l | tr -d ' ')
  say "dt check: message in $m of $n face-gate leg logs"

  # ---- k1: the stage-1 tests, in the packet's order -----------------------------------------------------------
  t K1_face_dirichlet gate $PY $U/face_dirichlet/test                                  # a7
  t K1_kant_face gate $PY $V/kant/sp_kant_onset/test_face                              # a3 (i)
  t K1_siteA_vs_face gate $PY $S/siteA_vs_face.py                                       # a3 (ii)
  t K1_pressure_face gate $PY $V/kant/sp_kant_pressure_sweep/test_face                 # a4
  t K1_dirichlet_identity gate $PY $S/dirichlet_identity.py                             # a5
  t K1_hu_lefm report $PY $V/hu/hu_spall_onset/test_lefm                                # m4a, m4b
  t K1_gradient_blind gate $PY $U/gradient_blind/test                                   # a6
  t K1_retired_keys gate $PY $U/retired_keys/test
  t K1_interlocks gate $PY $U/interlocks/test
  t K1_lateral_flux_lock gate $PY $U/lateral_flux_lock/test
  t K1_buckle_lock gate $PY $U/buckle_lock/test
  $PY $U/face_removal/test --check > $S/K1_face_removal.log 2>&1
  if diff -q $S/K1_face_removal.log tests/MMWSpalling/studies/p2_19_flake_depth/lock_face_removal.log > /dev/null; then
    say "  K1_face_removal --check: byte-identical to P2-19's log ($(tail -1 $S/K1_face_removal.log | cut -c1-60))"
  else gate_fail "face_removal --check differs from P2-19's recorded log"; fi
  t K1_pressure_siteA gate $PY $V/kant/sp_kant_pressure_sweep/test                      # a6 (the last site-A runs)
  if long_ok "the old sp_kant_onset/test (5 cases + the 25-case sweep)"; then
    t K1_kant_siteA gate $PY $V/kant/sp_kant_onset/test                                 # a6
  fi
  if long_ok "b5's BEFORE half (CAL 60 s)"; then
    $PY $V/kant/grimsel_calibration/run.py CAL > $S/K1_cal_before.log 2>&1
    rc=$?
    rm -rf $S/cal_before && cp -R $V/kant/grimsel_calibration/legs/CAL $S/cal_before
    say "  b5 BEFORE: CAL rc $rc on the stage-1 binary, archived to cal_before/ ($(ls $S/cal_before | tr '\n' ' '))"
  fi
  sleeps
  [ $GATE_FAIL = 0 ] && say "STAGE 1 DONE -- every gate held; STOP for the user (stage 2 needs the go)" \
                     || say "STAGE 1 DONE WITH A FAILED GATE -- see above"
fi

# 2: stage 2 (p2_28) -- the column path deleted (stage2_src.py, stage2_tests.py ran before this chain); the face gate
# re-pins under the component + thermo audit; every kept test; CAL after vs before; the tree gate is NOT run (a record)
if [ $STAGE = 2 ]; then
  mkdir -p $S/bin_2 && cp $BIN $S/bin_2/
  { echo "P2-28 stage 2 (p2_28): the column path deleted (one removal path, the face list)"
    shasum -a 256 $BIN
    shasum -a 256 src/Integrator/MMWSpalling.H src/Integrator/MMWSpalling/*.H; } > $S/BINARY_2.txt
  mkdir -p $S/pre2_flat
  cp $S/src_pre_2/src/Integrator/MMWSpalling.H $S/src_pre_2/src/Integrator/MMWSpalling/Removal.H \
     $S/src_pre_2/src/Integrator/MMWSpalling/FaceRemoval.H $S/pre2_flat/
  # ---- g2 / b1 / b2 ----
  $PY $W/face_gate.py run p2_28 > $S/face_2.log 2>&1 || gate_fail "face_gate run (face_2.log)"
  $PY $W/face_gate.py cmp p2_28a p2_28 > $S/GATE_b1_cmp.txt 2>&1
  say "b1 cmp p2_28a -> p2_28 (DIFFERENT by construction): $(tail -1 $S/GATE_b1_cmp.txt | cut -c1-200)"
  t K2_audit gate $PY $S/audit_components.py face p2_28a p2_28
  if $PY $S/audit_components.py face p2_28a p2_28 --mutate > $S/K2_audit_mutate.log 2>&1; then
    gate_fail "the planted mutation (Sp_field undeclared) PASSED the audit"
  else say "  K2_audit --mutate: FAILS as it must ($(tail -1 $S/K2_audit_mutate.log))"; fi
  t K2_rankcheck gate $PY $W/face_gate.py rankcheck p2_28
  n=$(find $W/face_p2_28 -name run.log | wc -l | tr -d ' ')
  m=$(grep -l "explicit time-step limit, level 0" $(find $W/face_p2_28 -name run.log) | wc -l | tr -d ' ')
  say "dt check: message in $m of $n face-gate leg logs"
  # ---- b3: the keys and the locks ----
  t K2_retired_keys gate $PY $U/retired_keys/test --pre $S/pre2_flat --red R11
  t K2_lateral_flux_lock gate $PY $U/lateral_flux_lock/test
  t K2_buckle_lock gate $PY $U/buckle_lock/test
  t K2_interlocks gate $PY $U/interlocks/test
  # ---- b4: every kept test ----
  for x in face_feet face_vapor face_melt face_beam face_heat face_list face_buckle face_dirichlet gradient_blind \
           robin_face convective_patch latent_heat beam beam_tau beam_schedule profile_depth; do
    t K2_$x gate $PY $U/$x/test
  done
  t K2_face_removal_check report $PY $U/face_removal/test --check
  t K2_amr_regrid report $PY $U/amr_microstructure_regrid/test
  t K2_kant_face report $PY $V/kant/sp_kant_onset/test_face
  t K2_pressure_face report $PY $V/kant/sp_kant_pressure_sweep/test_face
  t K2_hu_lefm report $PY $V/hu/hu_spall_onset/test_lefm
  # ---- b5: CAL after vs before ----
  if long_ok "b5's AFTER half (CAL 60 s)"; then
    $PY $V/kant/grimsel_calibration/run.py CAL > $S/K2_cal_after.log 2>&1
    t K2_cal_compare gate $PY $S/cal_compare.py
  fi
  # ---- m1 ----
  wc -l src/Integrator/MMWSpalling.H src/Integrator/MMWSpalling/*.H src/Numeric/MMWBeam.H > $S/LINES_2.txt
  say "m1 lines: $(tail -1 $S/LINES_2.txt) (12,696 + 272 at stage 1)"
  sleeps
  [ $GATE_FAIL = 0 ] && say "STAGE 2 DONE -- every gate held" || say "STAGE 2 DONE WITH A FAILED GATE -- see above"
fi

# 1b: the resume after the first stage-1 chain's one failed gate (lateral_flux_lock counted FaceListHeat's call
# sites: 3, not 2 -- the P2-28 prescribed branch; the lock restated structurally, declared): the lock on the SAME
# binary, then the two long legs the failure skipped. The gates and the other k1 tests stand from the first chain.
if [ $STAGE = 1b ]; then
  grep -q "$B" $S/BINARY_1.txt || { say "STOP: bin/ is not BINARY_1.txt's binary"; exit 1; }
  t K1b_lateral_flux_lock gate $PY $U/lateral_flux_lock/test
  if long_ok "the old sp_kant_onset/test (5 cases + the 25-case sweep)"; then
    t K1_kant_siteA gate $PY $V/kant/sp_kant_onset/test                                 # a6
  fi
  if long_ok "b5's BEFORE half (CAL 60 s)"; then
    $PY $V/kant/grimsel_calibration/run.py CAL > $S/K1_cal_before.log 2>&1
    rc=$?
    rm -rf $S/cal_before && cp -R $V/kant/grimsel_calibration/legs/CAL $S/cal_before
    say "  b5 BEFORE: CAL rc $rc on the stage-1 binary, archived to cal_before/ ($(ls $S/cal_before | tr '\n' ' '))"
  fi
  sleeps
  [ $GATE_FAIL = 0 ] && say "STAGE 1 DONE (1b) -- every gate held; STOP for the user (stage 2 needs the go)" \
                     || say "STAGE 1b DONE WITH A FAILED GATE -- see above"
fi
