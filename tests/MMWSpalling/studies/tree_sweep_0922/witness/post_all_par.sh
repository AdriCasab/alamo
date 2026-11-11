#!/bin/bash
# P2-28 RECORD: the tree gate retired with the column removal path (one removal path, the face list); p2_27.md5 is the column era's last list, its outputs moved with their tests to retired_units/. Kept, not run: face_gate.py is the gate.
# Parallel tree sweep. Same 11 tests, same hash.sh, same 2639-file list as
# post_all.sh + wit.sh -- only the SCHEDULE differs.
#
#   usage: guarded_post_all_par.sh <label>     (or call this directly)
#
# WHY THIS IS SAFE. Byte-identity is the oracle, and nothing in it depends on
# timing: each test keeps its own rank count and its own decomposition, and this
# line already proved the solver bit-identical across decompositions (the
# removed_mf 0-ghost fix). Concurrency changes only when a test runs, never what
# it computes. Every test writes to its own output directory, so there is no
# collision:
#   jet_d2c -> unit/jet_d2c/output          dev2d   -> validation/.../output/dev2d
#   <unit>  -> unit/<name>/output           s1_smoke-> studies/s1_.../output_smoke
#   (P2-27: spall_event, which also wrote output_no_event, retired with the damage lineage.)
#
# ⚠ It is validated by REPRODUCTION, not by argument: run it, then diff its md5
# list against the serial list. If a single hash moves, the schedule is not
# innocent and this script is wrong.
#
# LANES. macOS ships bash 3.2, which has no `wait -n`, so the pool is three
# static lanes balanced longest-processing-time-first from measured times
# (D2z-2 sweep, seconds):
#   jet_enthalpy 239, robin_feet 217, jet_d2c 189, dev2d 64, robin_jet 25,
#   robin_pinned 19, beam_void_closure 13, scalar_flaw 11, s1_smoke 9,
#   robin_face 7, spall_event 1.   total 794 s serial.
# Lane totals below are 265 / 293 / 236 s, so the critical path is ~293 s and
# the sweep drops from ~13.2 min to ~5 min. Rebalance here if a test's cost
# changes; the times are in the log lines this script prints.
#
# ⚠ Three lanes, not eleven: each test runs mpirun -np 4, so eleven at once is
# 44 ranks on a machine that has far fewer cores. Oversubscription past this
# point costs more in contention than it saves in scheduling.
set -u
S=/Users/tzetze20/amr_tools/alamo/tests/MMWSpalling/studies/tree_sweep_0922/witness
LABEL=${1:?usage: post_all_par.sh <label>}
W=$S/${LABEL}_wit; mkdir -p "$W"
cd /Users/tzetze20/amr_tools/alamo || exit 1
VENV=/Users/tzetze20/Desktop/code/.venv/bin/python
RUN="mpirun --oversubscribe --bind-to none -np 4"
T=tests/MMWSpalling

# rc is REPORTED, never fatal: beam_void_closure exits 1 by design (its T3 is
# D2z-1's documented debt) and the sweep must still produce its hash list.
run_test () {
  name=$1; shift
  local t0=$SECONDS
  echo "=== $name start $(date +%T)"
  "$@" > "$W/$name.log" 2>&1
  echo "=== $name rc=$? $((SECONDS - t0))s"
}

lane_a () {
  run_test jet_enthalpy      $VENV $T/unit/jet_enthalpy/test
  run_test robin_pinned      $VENV $T/unit/robin_pinned/test
  run_test robin_face        $VENV $T/unit/robin_face/test
}
lane_b () {
  run_test robin_feet        $VENV $T/unit/robin_feet/test
  run_test dev2d             $RUN bin/mmwspalling-3d-g++ \
                             $T/validation/meier/sp_meier_pilot/input_2d_dev
  run_test scalar_flaw       $VENV $T/unit/scalar_flaw/test
}
lane_c () {
  run_test jet_d2c           $VENV $T/unit/jet_d2c/test
  run_test robin_jet         $VENV $T/unit/robin_jet/test
  run_test beam_void_closure $VENV $T/unit/beam_void_closure/test
  run_test s1_smoke          $VENV $T/studies/s1_surface_resolution/run_sweep.py \
                             --smoke --force
}

SECONDS=0
lane_a & A=$!
lane_b & B=$!
lane_c & C=$!
wait $A; wait $B; wait $C
echo "=== all lanes joined at ${SECONDS}s"

# hash.sh LAST, after every lane has written its outputs.
$S/hash.sh "$S/$LABEL.md5"
echo "=== hashed $(wc -l < "$S/$LABEL.md5" | tr -d ' ') files in ${SECONDS}s total"
echo ALLDONE
