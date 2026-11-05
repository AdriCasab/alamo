#!/bin/bash
# P2-17 RESUME from a4 after the a3 comparator restatement (chain.sh otherwise); STOPS at the first failed gate.
#   ord 1 z0 (gate), z1, z2 | ord 2 CAL1 -> a1, face_heat -> a2, O_S164 -> a3, O_S140 -> a4, N_S070 -> a5 (gates)
#   ord 3 the remaining 16 ladder legs | ord 4 RESULTS.txt (c1-c4, e0, e1)
# No src/ change and no build: bin/ must still be P2-16's binary.
set -u
ROOT=/Users/tzetze20/amr_tools/alamo; cd $ROOT
S=tests/MMWSpalling/studies/p2_17_buckle_on
U=tests/MMWSpalling/unit
PY=/Users/tzetze20/Desktop/code/.venv/bin/python; LOG=$S/CHAIN.log
say()  { echo "$(date +%T) $*" | tee -a $LOG; }
stop() { say "STOP: $*"; exit 1; }
caffeinate -dis -w $$ &
B=$(shasum -a 256 bin/mmwspalling-3d-g++ | cut -c1-64)
say "P2-17 chain (binary ${B:0:16})"
[ "$B" = "38c8623885c5b8b61d10597453aa380882cf22ea6c4c75274fc6cb7bdf1075ba" ] || stop "bin/ is not P2-16's binary"
echo "$B  bin/mmwspalling-3d-g++" > $S/BINARY.txt
say "RESUME at a4 (a3 restated by step: A3.first.txt is the first scoring; O_S164 not re-run)"
tail -1 $S/A3.txt | tee -a $LOG
grep -q "^a3 PASS" $S/A3.txt || stop "a3 ($S/A3.txt)"
say "a4 O_S140"
$PY $S/run.py O_S140 >> $S/RUN.log 2>&1
$PY $S/reduce.py A4 > /dev/null; cat $S/A4.txt | tee -a $LOG
grep -q "^a4 PASS" $S/A4.txt || stop "a4 ($S/A4.txt)"
say "a5 N_S070"
$PY $S/run.py N_S070 >> $S/RUN.log 2>&1
$PY $S/reduce.py A5 > /dev/null; tail -2 $S/A5.txt | tee -a $LOG
grep -q "^a5 PASS" $S/A5.txt || stop "a5 ($S/A5.txt)"

say "ord 3: the remaining ladder legs"
$PY $S/run.py O_S200 O_S120 O_S100 O_S087 O_S070 N_S200 N_S164 N_S140 N_S120 N_S100 N_S087 \
    A_S140 A_S120 A_S100 A_S087 A_S070 >> $S/RUN.log 2>&1
say "ord 4: RESULTS.txt"
$PY $S/reduce.py > $S/RESULTS.txt 2>&1
sed -n '/^e1 VERDICT/,$p' $S/RESULTS.txt | tee -a $LOG
say "CHAIN DONE"
