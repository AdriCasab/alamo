#!/bin/bash
# Stops four legs as soon as the pools start them, and records it. They are the 0.098 mm legs at absorption depths 4 and
# 8 mm: 41 and 82 cells per absorption depth, where the 0.195 mm legs already have 20 and 41. They answer nothing the
# coarser legs do not, and each would cost hours (the line scan runs every step while a firing is being withdrawn).
# A stopped leg keeps a non-zero rc; reduce.py prints it as "not run".
S=/Users/tzetze20/amr_tools/alamo/tests/MMWSpalling/studies/beam_column_1110
while pgrep -f "beam_column_1110/run.py" > /dev/null; do
  for L in D4_Q2_N512_A D4_Q2_N512_H D8_Q2_N512_A D8_Q2_N512_H; do
    if pkill -f "mmwspalling-3d-g.. .*legs/$L/run"; then echo "$(date '+%F %T') stopped $L (cost; see skip_watch.sh)" >> $S/SKIPPED.txt; fi
  done
  sleep 15
done
