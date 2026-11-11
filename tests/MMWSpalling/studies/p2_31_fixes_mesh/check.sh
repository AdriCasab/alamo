#!/usr/bin/env bash
# P2-31: the ONE driver. `caffeinate -dis bash tests/MMWSpalling/studies/p2_31_fixes_mesh/check.sh <1|2a|2b|3|4>` runs
# that stage's Criteria rows in table order (stages.py), appends to S/RESULTS.txt, honours the driver's exit status and
# exits non-zero at the first red GATE row. A red verdict / report row is recorded and the chain continues.
set -u
ROOT=/Users/tzetze20/amr_tools/alamo
S=$ROOT/tests/MMWSpalling/studies/p2_31_fixes_mesh
PY=/Users/tzetze20/Desktop/code/.venv/bin/python
cd "$ROOT" || exit 2
STAGE=$1
: > "$S/CHECK_$STAGE.log"
say() { echo "$*" | tee -a "$S/RESULTS.txt"; }
say ""
say "== stage $STAGE started $(date '+%H:%M:%S') on $(shasum -a 256 bin/mmwspalling-3d-g++ | cut -c1-16)"
$PY "$S/stages.py" "$STAGE"
rc=$?
if [ $rc -ne 0 ]; then say "== stage $STAGE RED: the driver exited $rc at $(date '+%H:%M:%S')"; exit 1; fi
say "== stage $STAGE GREEN $(date '+%H:%M:%S')"
