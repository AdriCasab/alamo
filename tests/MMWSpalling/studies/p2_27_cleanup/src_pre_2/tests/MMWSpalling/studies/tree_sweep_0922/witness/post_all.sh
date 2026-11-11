#!/bin/bash
# Copied from the D2c scratchpad on 2026-09-18 (D2e Goal 0); paths point here.
S=/Users/tzetze20/amr_tools/alamo/tests/MMWSpalling/studies/tree_sweep_0922/witness
cd /Users/tzetze20/amr_tools/alamo
SECONDS=0
/Users/tzetze20/Desktop/code/.venv/bin/python tests/MMWSpalling/unit/jet_d2c/test > tests/MMWSpalling/unit/jet_d2c/output/test_pass.log 2>&1
echo "=== jet_d2c rc=$? ${SECONDS}s"
$S/wit.sh ${1:?label, e.g. perf or post_d2e}
