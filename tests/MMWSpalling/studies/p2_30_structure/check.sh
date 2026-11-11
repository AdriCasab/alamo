#!/usr/bin/env bash
# P2-30: the ONE driver. `caffeinate -dis bash tests/MMWSpalling/studies/p2_30_structure/check.sh <1|2|3|4a|4b>` runs
# that stage's Criteria rows in table order, appends to S/RESULTS.txt, and exits non-zero at the first red GATE row.
# A red verdict / report row is recorded and the chain continues. Output also goes to S/CHECK_<stage>.log.
set -u
ROOT=/Users/tzetze20/amr_tools/alamo
S=$ROOT/tests/MMWSpalling/studies/p2_30_structure
W=$ROOT/tests/MMWSpalling/studies/tree_sweep_0922/witness
T=$ROOT/tests/MMWSpalling
PY=/Users/tzetze20/Desktop/code/.venv/bin/python
cd "$ROOT" || exit 2
STAGE=$1
LOG=$S/CHECK_$STAGE.log
: > "$LOG"
say() { echo "$*" | tee -a "$LOG" | tee -a "$S/RESULTS.txt" >/dev/null; echo "$*"; }
say ""
say "== stage $STAGE started $(date '+%H:%M:%S')"

# a row: row <gate|verdict|report> <name> <command...>   the command's output lands in the log; its rc is the row's result
row() {
  local kind=$1 name=$2; shift 2
  local t0=$SECONDS
  "$@" >> "$LOG" 2>&1
  local rc=$?
  say "  $name: $([ $rc -eq 0 ] && echo PASS || echo RED) (rc $rc, $((SECONDS - t0)) s) [$kind]"
  if [ $rc -ne 0 ] && [ "$kind" = gate ]; then
    say "== RED GATE $name: stopped at $(date '+%H:%M:%S')"; exit 1
  fi
  return $rc
}
build() {
  EIGEN=$ROOT/ext CPLUS_INCLUDE_PATH=/opt/homebrew/include \
    LIBRARY_PATH=/opt/homebrew/lib:/opt/homebrew/Cellar/gcc/15.2.0_1/lib/gcc/current make -j8
}
binary_note() { shasum -a 256 "$ROOT/bin/mmwspalling-3d-g++" | tee -a "$LOG" > "$S/BINARY_$1.txt"; say "  binary $1: $(cut -c1-16 "$S/BINARY_$1.txt")"; }

case "$STAGE" in
1)
  row gate a1_identity  $PY "$S/pp_identity.py" src_pre_1
  row gate a1_edges     $PY "$S/slice_edges.py"
  row gate a2_capture   $PY "$S/locks.py" capture stage1
  row gate a2_lines     $PY "$S/locks.py" cmp pre_1 stage1
  row gate a2_bytes     $PY "$S/locks.py" bytes src_pre_1
  grep -E 'tokens compared|VERDICT' "$LOG" | tee -a "$S/RESULTS.txt"
  grep -E 'every line identical|byte-identical|DIFFERENT' "$LOG" | tee -a "$S/RESULTS.txt"
  ;;
*)
  $PY "$S/stages.py" "$STAGE"
  rc=$?
  if [ $rc -ne 0 ]; then say "== stage $STAGE RED: the driver exited $rc at $(date '+%H:%M:%S')"; exit 1; fi
  ;;
esac
say "== stage $STAGE GREEN $(date '+%H:%M:%S')"
