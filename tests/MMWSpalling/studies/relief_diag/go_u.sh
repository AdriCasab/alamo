#!/bin/bash
S=/Users/tzetze20/amr_tools/alamo/tests/MMWSpalling/studies/relief_diag
PY=/Users/tzetze20/Desktop/code/.venv/bin/python
cd /Users/tzetze20/amr_tools/alamo || exit 1
caffeinate -dis $PY "$S/run_uniform.py" > "$S/RUN_U.out" 2>&1
$PY "$S/reduce_uniform.py" > "$S/REDUCE_U.out" 2>&1
echo done
