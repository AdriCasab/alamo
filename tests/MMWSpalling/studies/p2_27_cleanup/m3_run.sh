#!/bin/bash
# P2-27 m3 (the stage-5 chain's last step, run on its own after k5's interlocks re-run): a 20 s prefix of
# validation/meier/sp_meier_pilot/input_feet_face on the stage-5 binary -- recorded, not scored.
set -u
cd /Users/tzetze20/amr_tools/alamo
S=tests/MMWSpalling/studies/p2_27_cleanup; LOG=$S/CHAIN.log; BIN=bin/mmwspalling-3d-g++
grep -q "$(shasum -a 256 $BIN | cut -c1-64)" $S/BINARY_5.txt || { echo "not the stage-5 binary" | tee -a $LOG; exit 1; }
caffeinate -dis -w $$ &
echo "$(date +%T) [s5] k5 interlocks re-run with the removed scalars on the pre binary's argv: $(tail -1 $S/K5_interlocks.log | cut -c1-80)" >> $LOG
echo "$(date +%T) [s5] m3: a 20 s prefix of input_feet_face (recorded, not scored)" >> $LOG
O=$S/m3_feet_face; rm -rf $O; mkdir -p $O; t0=$SECONDS
mpirun --oversubscribe --bind-to none -np 4 $BIN tests/MMWSpalling/validation/meier/sp_meier_pilot/input_feet_face \
  stop_time=20.0 plot_file=$O/run > $S/M3_feet_face.log 2>&1
rc=$?
n4=$(awk -F, 'NR>1 && $4==4' $O/run_removal_events.csv 2>/dev/null | wc -l | tr -d ' ')
n1=$(awk -F, 'NR>1 && $4==1' $O/run_removal_events.csv 2>/dev/null | wc -l | tr -d ' ')
echo "  m3 input_feet_face 20 s: rc $rc, wall $((SECONDS - t0)) s; regime-4 rows $n4, spall rows $n1; $(grep -a 'ABORT' $S/M3_feet_face.log | head -1 | cut -c1-160)" | tee -a $LOG
echo "$(date +%T) [s5] STAGE 5 DONE" >> $LOG
