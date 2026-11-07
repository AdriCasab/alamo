#!/bin/bash
# the 1 mm twins of Lc and M (mesh sensitivity of the ratio)
S=/Users/tzetze20/amr_tools/alamo/tests/MMWSpalling/studies/relief_diag
PY=/Users/tzetze20/Desktop/code/.venv/bin/python
cd /Users/tzetze20/amr_tools/alamo || exit 1
rm -rf "$S/Lc1" "$S/M1"
caffeinate -dis $PY "$S/run.py" Lc1 M1 > "$S/RUN1.out" 2>&1
echo "go1.sh done rc=$?"
