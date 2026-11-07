#!/bin/bash
# one-command driver: clean the comparison legs and run them (Lc, then M), then reduce
S=/Users/tzetze20/amr_tools/alamo/tests/MMWSpalling/studies/relief_diag
PY=/Users/tzetze20/Desktop/code/.venv/bin/python
cd /Users/tzetze20/amr_tools/alamo || exit 1
rm -rf "$S/Lc" "$S/M"
caffeinate -dis $PY "$S/run.py" Lc M > "$S/RUN.out" 2>&1
$PY "$S/reduce.py" > "$S/REDUCE.out" 2>&1
echo "go.sh done rc=$?"
