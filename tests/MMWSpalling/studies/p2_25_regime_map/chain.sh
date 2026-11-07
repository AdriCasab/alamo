#!/bin/bash
# P2-25 in the packet's run order; STOPS at the first failed gate. Log: CHAIN.log.
#   ord 1  z0 (the ladders, the stop rule, the binary), z1 (the comparator), z2 (the vapour mass on P2-24's legs)
#   ord 2  a0 tree gate vs p2_18 (controls d2aa 26, d2t 863); a1 face gate run p2_25, cmp p2_24 p2_25 (IDENTICAL or
#          STOP -- no re-pin); a2 unit/{retired_keys,latent_heat,beam_tau,face_beam,face_vapor}, face_melt --bin_pre
#   ord 3  b0 P400R (gate); the 42 map legs by flux row, shallow to deep; reduce; the K legs (14); reduce; the W (4) and
#          F (5) legs; check_shared on D1_Q2 and D71_Q2; reduce
#   ord 4  RESULTS.txt; prune every leg (first and last plotfile, MANIFEST.md5 first)
set -u
ROOT=/Users/tzetze20/amr_tools/alamo; cd $ROOT
S=tests/MMWSpalling/studies/p2_25_regime_map
U=tests/MMWSpalling/unit
W=tests/MMWSpalling/studies/tree_sweep_0922/witness
TC=tests/MMWSpalling/studies/p2_9_face_removal/tree_compare.py
OT=tests/MMWSpalling/validation/oglesby/test
CS=tests/MMWSpalling/validation/shared/check_shared.py
P24=tests/MMWSpalling/studies/p2_24_melt_vapour
BIN=bin/mmwspalling-3d-g++
PY=/Users/tzetze20/Desktop/code/.venv/bin/python; LOG=$S/CHAIN.log
say()  { echo "$(date +%T) $*" | tee -a $LOG; }
stop() { say "STOP: $*"; exit 1; }
caffeinate -dis -w $$ &
B=$(shasum -a 256 $BIN | cut -c1-64)
say "P2-25 chain (binary ${B:0:16})"
grep -q "$B" $S/BINARY.txt || stop "bin/ is not the P2-25 binary (BINARY.txt)"

say "ord 1: z0, z1, z2"
$PY $S/z0.py > $S/Z0.txt 2>&1 || { tail -5 $S/Z0.txt >> $LOG; stop "z0"; }
tail -1 $S/Z0.txt | tee -a $LOG
$PY $S/run.py --check | tee -a $LOG | grep -q EQUALS || stop "run.py's literal table differs from STOPS.json"
$PY $S/z1.py > $S/Z1.txt 2>&1 || { cat $S/Z1.txt >> $LOG; stop "z1"; }
grep "OLD rule\|col_moved\|^z1" $S/Z1.txt | cut -c1-200 | tee -a $LOG
$PY $OT $P24/legs/M0/run --onset --alt $P24/legs/M100/run --flat $P24/legs/M0F/run > $S/Z2.txt 2>&1
grep "vapour mass" $S/Z2.txt | tee -a $LOG
grep -q "vapour mass 3.28 g" $S/Z2.txt && grep -q "vapour mass 0 g" $S/Z2.txt && grep -q "vapour mass 14.4 g" $S/Z2.txt || stop "z2"

say "ord 2: a0 tree gate"
bash $W/guarded_post_all_par.sh p2_25 > $S/tree.log 2>&1
[ $? = 2 ] && stop "the binary moved during the sweep"
$PY $TC p2_25 p2_18 d2aa d2t > $S/GATE_tree.txt 2>&1; cat $S/GATE_tree.txt >> $LOG
grep -q "vs p2_18: 0 of 2639 changed, 0 only in new, 0 only in ref" $S/GATE_tree.txt || stop "a0 vs p2_18"
grep -q "vs d2aa: 26 of 2639 changed" $S/GATE_tree.txt || stop "a0 control d2aa"
grep -q "vs d2t: 863 of 2639 changed" $S/GATE_tree.txt || stop "a0 control d2t"
say "a0 PASS; a1 face gate"
$PY $W/face_gate.py run p2_25 > $S/face.log 2>&1 || stop "face_gate run ($S/face.log)"
$PY $W/face_gate.py cmp p2_24 p2_25 > $S/GATE_face_cmp.txt 2>&1
tail -3 $S/GATE_face_cmp.txt | cut -c1-220 | tee -a $LOG
grep -q "VERDICT: IDENTICAL" $S/GATE_face_cmp.txt || stop "a1 is not IDENTICAL -- the finding; no re-pin ($S/GATE_face_cmp.txt)"
say "a2 unit tests"
for t in retired_keys latent_heat beam_tau face_beam face_vapor; do
  $PY $U/$t/test > $S/A2_$t.log 2>&1 || { tail -15 $S/A2_$t.log >> $LOG; stop "a2 $t"; }
  echo "  $t: $(tail -1 $S/A2_$t.log | cut -c1-120)" | tee -a $LOG
done
$PY $U/face_melt/test --bin_pre $P24/bin_p2_24/mmwspalling-3d-g++ > $S/A2_face_melt.log 2>&1 || { tail -15 $S/A2_face_melt.log >> $LOG; stop "a2 face_melt"; }
echo "  face_melt: $(tail -1 $S/A2_face_melt.log)" | tee -a $LOG

say "ord 3: b0 P400R"
$PY $S/run.py P400R >> $S/RUN.log 2>&1 || stop "P400R run"
$PY $S/reduce.py > $S/RESULTS.txt 2>&1
grep "b0\|L10" $S/RESULTS.txt | cut -c1-220 | tee -a $LOG
grep -q "PASS b0" $S/RESULTS.txt || stop "b0: the rule does not remove on this binary"
say "  the 42 map legs (by flux row, shallow to deep)"
for q in 0p3 0p5 1 2 4 8; do
  for d in 0p5 1 2 4 8 16 71; do
    $PY $S/run.py D${d}_Q${q} >> $S/RUN.log 2>&1 || say "  D${d}_Q${q} rc != 0 (legs/D${d}_Q${q}/run.log)"
  done
  say "  flux row q $q done"
done
$PY $S/reduce.py > $S/RESULTS.txt 2>&1
sed -n '/letters only/,/C cells/p' $S/RESULTS.txt | tee -a $LOG
say "  K legs (14)"
for k in K05 K2; do for d in 0p5 1 2 4 8 16 71; do
  $PY $S/run.py ${k}_D${d}_Q2 >> $S/RUN.log 2>&1 || say "  ${k}_D${d}_Q2 rc != 0"
done; done
$PY $S/reduce.py > $S/RESULTS.txt 2>&1
say "  W (4) and F (5) legs"
for l in W6_D1_Q2 W6_D4_Q2 W19_D1_Q2 W19_D4_Q2 F_D1_Q2 F_D2_Q2 F_D4_Q2 F_D2_Q1 F_D2_Q4; do
  $PY $S/run.py $l >> $S/RUN.log 2>&1 || say "  $l rc != 0"
done
say "check_shared on D1_Q2 (--allow beam.tan_delta, with its argv.sh) and D71_Q2"
{ $PY $CS $S/input $S/legs/D1_Q2/argv.sh --allow beam.tan_delta --run $S/legs/D1_Q2/run
  $PY $CS $S/input $S/legs/D71_Q2/argv.sh --run $S/legs/D71_Q2/run; } > $S/CHECK_SHARED.txt 2>&1
grep "check_shared" $S/CHECK_SHARED.txt | cut -c1-160 | tee -a $LOG

say "ord 4: RESULTS.txt, prune"
$PY $S/reduce.py > $S/RESULTS.txt 2>&1
sed -n '/^e0 VERDICT/,$p' $S/RESULTS.txt | cut -c1-240 | tee -a $LOG
$PY $S/prune.py --all
du -sh $S/legs | tee -a $LOG
say "CHAIN DONE"
