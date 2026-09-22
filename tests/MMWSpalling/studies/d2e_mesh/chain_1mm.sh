#!/bin/bash
# D2e: start Pb_1mm when Pa_1mm is done and Pc_1mm when B0_1mm is done (at
# most two 1 mm runs at a time), only if the 1 mm step gate PASSED.
cd "$(dirname "$0")"
PY=/Users/tzetze20/Desktop/code/.venv/bin/python
until grep -q "GATE" gate_watch.log 2>/dev/null; do sleep 60; done
if ! grep -q "GATE PASSED" gate_watch.log; then echo "gate failed: chain not started"; exit 1; fi
( until [ -f output/Pa_1mm.done ]; do sleep 60; done
  echo "Pb_1mm start $(date +%T)"; $PY run.py --cases Pb_1mm --jobs 1 > campaign_1mm_b.log 2>&1; echo "Pb_1mm end $(date +%T)" ) &
( until [ -f output/B0_1mm.done ]; do sleep 60; done
  echo "Pc_1mm start $(date +%T)"; $PY run.py --cases Pc_1mm --jobs 1 > campaign_1mm_c.log 2>&1; echo "Pc_1mm end $(date +%T)" ) &
wait
echo CHAIN_DONE
