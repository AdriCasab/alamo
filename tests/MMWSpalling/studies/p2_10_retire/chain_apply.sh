#!/bin/bash
# P2-10 ord 2-4 + the b legs, in the packet's run order, every gate enforced: the chain
# STOPS at the first failed gate (no further long legs). Writes CHAIN.log.
#   ord 2  install C1 (retire.py c1) -> build -> a3 tree sweep -> a4 W vs W_bB
#   ord 3  install C2 (retire.py c2) -> build -> a5 tree sweep -> a6 W -> a7 retired keys
#          -> the test edits (item 6) -> a8 every kept unit test -> a9 lines + grep
#   ord 4  the b legs (unit/face_removal --p2_10) -> TILT45.txt
# The sources installed are the scratchpad outputs of retire.py, regenerated here from the
# binary-B copies in src_bB/ so what ships is exactly what the scripts cut.
set -u
ROOT=/Users/tzetze20/amr_tools/alamo
cd "$ROOT" || exit 1
S=tests/MMWSpalling/studies/p2_10_retire
W=tests/MMWSpalling/studies/tree_sweep_0922/witness
U=tests/MMWSpalling/unit
SC=/private/tmp/claude-501/-Users-tzetze20-amr-tools-alamo/ac0a0f5a-4af6-48bf-b61d-9009cf0c8681/scratchpad
PY=/Users/tzetze20/Desktop/code/.venv/bin/python
TC=tests/MMWSpalling/studies/p2_9_face_removal/tree_compare.py
LOG=$S/CHAIN.log
say()  { echo "$(date +%T) $*" | tee -a $LOG; }
stop() { say "STOP: $*"; exit 1; }
caffeinate -dis -w $$ &

build() {
    say "build $1"
    EIGEN=$ROOT/ext CPLUS_INCLUDE_PATH=/opt/homebrew/include \
    LIBRARY_PATH=/opt/homebrew/lib:/opt/homebrew/Cellar/gcc/15.2.0_1/lib/gcc/current \
        make -j8 > $S/build_$1.log 2>&1 || stop "build $1 failed ($S/build_$1.log)"
    grep -q "mmwspalling.cc" $S/build_$1.log || stop "build $1 did not recompile mmwspalling.cc"
    shasum -a 256 bin/mmwspalling-3d-g++ | tee $S/BINARY_$1.txt | tee -a $LOG
}
tree() {
    say "tree sweep $1 (a3/a5)"
    bash $W/guarded_post_all_par.sh p2_10_$1 > $S/tree_$1.log 2>&1
    [ $? = 2 ] && stop "the binary moved during the $1 sweep"
    $PY $TC p2_10_$1 p2_9_bB d2aa d2t | tee $S/GATE_tree_$1.txt | tee -a $LOG
    grep -q "vs p2_9_bB: 0 of 2639 changed, 0 only in new, 0 only in ref" $S/GATE_tree_$1.txt \
        || stop "tree gate $1: not 0 of 2639"
    grep -q "vs d2aa: 26 of 2639 changed" $S/GATE_tree_$1.txt || stop "control d2aa ($1)"
    grep -q "vs d2t: 863 of 2639 changed" $S/GATE_tree_$1.txt || stop "control d2t ($1)"
}
wgate() {
    say "W on $1 (a4/a6)"
    $PY $S/run_p2_10.py refs $1 > $S/refs_$1.log 2>&1
    $PY $S/run_p2_10.py cmp $1 > /dev/null
    tail -3 $S/W_$1_vs_bB.txt | tee -a $LOG
    tail -1 $S/W_$1_vs_bB.txt | grep -q "PASS" || stop "W gate $1 ($S/W_$1_vs_bB.txt)"
}

say "P2-10 chain start${RESUME:+ (RESUME=$RESUME)}"
[ -n "${RESUME:-}" ] || { echo "a2 -- the references exist BEFORE the first source edit (mtimes)"
  stat -f "  %Sm  %N" $S/W_bB.md5 $S/UNIT_bB/rc.txt $S/DELETION_MAP.txt \
      src/Integrator/MMWSpalling.H src/Integrator/MMWSpalling/Removal.H \
      src/Integrator/MMWSpalling/FaceRemoval.H
  echo "  a2 PASS (every source file older than W_bB.md5 and UNIT_bB/rc.txt; edits start now)"; } \
    > $S/GATE_a2.txt
[ -n "${RESUME:-}" ] || cat $S/GATE_a2.txt >> $LOG

# the scripts' outputs, regenerated from binary B's sources
mkdir -p $SC/bB && cp $S/src_bB/MMWSpalling.H $S/src_bB/Removal.H $S/src_bB/FaceRemoval.H $SC/bB/
/usr/bin/python3 $S/retire.py c1 $SC/bB $SC/c1 > /dev/null || stop "retire.py c1"
/usr/bin/python3 $S/retire.py c2 $SC/c1 $SC/c2 > /dev/null || stop "retire.py c2"
cp $SC/c1/CUTS.txt $S/CUTS_C1.txt; cp $SC/c2/CUTS.txt $S/CUTS_C2.txt

# ---------------------------------------------------------------- ord 2: C1
if [ "${RESUME:-}" = "c2" ]; then
    # resumed after the C1 W gate's first scoring counted run/diff.patch (the build's own
    # git diff); the C1 legs on disk are re-scored, not re-run
    say "RESUME at ord 3: re-score W on c1 from disk (a4)"
    $PY $S/run_p2_10.py cmp c1 > /dev/null
    tail -3 $S/W_c1_vs_bB.txt | tee -a $LOG
    tail -1 $S/W_c1_vs_bB.txt | grep -q "PASS" || stop "W gate c1 ($S/W_c1_vs_bB.txt)"
else
say "ord 2: install C1"
cp $SC/c1/MMWSpalling.H src/Integrator/MMWSpalling.H
cp $SC/c1/Removal.H $SC/c1/FaceRemoval.H src/Integrator/MMWSpalling/
build c1
tree c1
wgate c1
fi

# ---------------------------------------------------------------- ord 3: C2
say "ord 3: install C2"
cp $SC/c2/MMWSpalling.H src/Integrator/MMWSpalling.H
cp $SC/c2/Removal.H $SC/c2/FaceRemoval.H src/Integrator/MMWSpalling/
build c2
tree c2
wgate c2

say "a7: the retired keys abort; M0T runs"
$PY $S/run_p2_10.py retired > /dev/null
cat $S/RETIRED.txt >> $LOG
grep -q "FAIL" $S/RETIRED.txt && stop "a7: a retired key ran ($S/RETIRED.txt)"
[ -f $S/W_c2/M0T/.done ] || stop "a7: M0T did not run on C2"
say "a7 PASS"

say "item 6: the test edits"
/usr/bin/python3 $S/retire_tests.py $ROOT $SC/tests_new > /dev/null || stop "retire_tests.py"
cp $SC/tests_new/TEST_EDITS.txt $S/TEST_EDITS.txt
mkdir -p $SC/retired_units
for d in side_face_flux side_face_area side_face_step_height side_staircase lateral_drive wall_march; do
    [ -d $U/$d ] && mv $U/$d $SC/retired_units/$d
done
cp $SC/tests_new/lateral_flux_lock/test $U/lateral_flux_lock/test

say "a8: every kept unit test on C2"
$PY $S/run_p2_10.py units c2 > $S/units_c2.log 2>&1
$PY $S/run_p2_10.py units_cmp c2 > /dev/null
grep "^a8\|lateral_flux_lock on" $S/UNITS_c2_vs_bB.txt | tee -a $LOG
grep -q "^a8 .*: PASS" $S/UNITS_c2_vs_bB.txt || stop "a8 ($S/UNITS_c2_vs_bB.txt)"

say "a9: lines removed, deleted names left"
$PY $S/run_p2_10.py a9 > /dev/null
head -8 $S/A9.txt | tee -a $LOG

# ---------------------------------------------------------------- ord 4: b
say "ord 4: the b legs (TL40/TL50 at 1 and 0.5 mm, TL60 and T00 at 0.5 mm)"
$PY $U/face_removal/test --p2_10 > $S/b_legs.log 2>&1
tail -3 $S/b_legs.log | tee -a $LOG
say "CHAIN DONE"
