#!/bin/bash
# P2-22 -- the in-rock beam spread + test #19 at its true spot. STOPS at the first failed GATE; the reports run on.
#   ord 1  unit/beam_tau (U1-U3; --bin_pre = P2-20's binary, which must abort on the key), unit/retired_keys
#   ord 2  tree gate vs p2_18 (0/2639 expected: every tree deck has eps_r 1), face gate cmp p2_21 p2_22 (identical:
#          beam-off), beam_void_closure + unit/beam (eps_r 1: unchanged)
#   ord 3  C1 vs P2-20's Z1 (REPORT: expected to MOVE, the deck has eps_r 5.02; the magnitude is recorded),
#          A (table, true spot), X (A explicit dt 0.03), F (flat control); each scored by the Zhang test
#   ord 4  RESULTS.txt
set -u
ROOT=/Users/tzetze20/amr_tools/alamo; cd $ROOT
S=tests/MMWSpalling/studies/p2_22_inrock_spread
U=tests/MMWSpalling/unit
W=tests/MMWSpalling/studies/tree_sweep_0922/witness
TC=tests/MMWSpalling/studies/p2_9_face_removal/tree_compare.py
ZT=tests/MMWSpalling/validation/zhang/zhang_oglesby/test
BIN=bin/mmwspalling-3d-g++
PY=/Users/tzetze20/Desktop/code/.venv/bin/python; LOG=$S/CHAIN.log
say()  { echo "$(date +%T) $*" | tee -a $LOG; }
stop() { say "STOP: $*"; exit 1; }
caffeinate -dis -w $$ &
B=$(shasum -a 256 $BIN | cut -c1-64)
say "P2-22 chain (binary ${B:0:16})"
grep -q "$B" $S/BINARY.txt || stop "bin/ is not the P2-22 build"

prune() {   # keep the first plotfile at or after the onset time (if any) and the last
  $PY - "$S/legs/$1" <<'PYEOF'
import os, sys, math, shutil
d = sys.argv[1]; run = os.path.join(d, "run")
L = open(os.path.join(run, "thermo.dat")).read().splitlines()
h = [c.strip() for c in L[0].split("\t")]
rows = [[float(x) for x in r.split("\t")] for r in L[1:] if r.strip()]
t_on = next((r[0] for r in rows if r[0] > 0 and "bc_z_melt" in h and math.isfinite(r[h.index("bc_z_melt")])), None)
pfs = sorted(p for p in os.listdir(run) if p.endswith("cell"))
def ptime(pf):
    H = open(os.path.join(run, pf, "Header")).read().split("\n"); return float(H[2 + int(H[1]) + 1])
times = {pf: ptime(pf) for pf in pfs}
keep = {pfs[-1]}
if t_on is not None:
    after = [pf for pf in pfs if times[pf] >= t_on - 1e-9]
    if after: keep.add(after[0])
gone = [pf for pf in pfs if pf not in keep]
for pf in gone: shutil.rmtree(os.path.join(run, pf))
for pf in [p for p in os.listdir(run) if p.endswith("node")]: shutil.rmtree(os.path.join(run, pf))
open(os.path.join(d, "PRUNED.txt"), "w").write(
    "pruned after its test (guardrail): kept %s (onset row at %s s; the last plotfile); removed %d cell plotfiles\n"
    % (" ".join("%s (t %g s)" % (pf, times[pf]) for pf in sorted(keep)), t_on, len(gone)))
PYEOF
}

say "ord 1: unit/beam_tau (U1-U3)"
$PY $U/beam_tau/test --bin_pre tests/MMWSpalling/studies/p2_21_loss_table/bin_pre/mmwspalling-3d-g++ > $S/Z1_beam_tau.log 2>&1 || { tail -20 $S/Z1_beam_tau.log >> $LOG; stop "z1-z3"; }
grep "PASS  \|FAIL  " $S/Z1_beam_tau.log | cut -c1-170 | tee -a $LOG
say "unit/retired_keys"
$PY $U/retired_keys/test > $S/Z4_retired_keys.log 2>&1 || { tail -30 $S/Z4_retired_keys.log >> $LOG; stop "z4"; }
grep "^PASS\|^FAIL" $S/Z4_retired_keys.log | cut -c1-170 | tee -a $LOG

say "ord 2: tree gate"
bash $W/guarded_post_all_par.sh p2_22 > $S/tree.log 2>&1
[ $? = 2 ] && stop "the binary moved during the sweep"
$PY $TC p2_22 p2_18 d2aa d2t > $S/GATE_tree.txt 2>&1; cat $S/GATE_tree.txt >> $LOG
grep -q "vs p2_18: 0 of 2639 changed, 0 only in new, 0 only in ref" $S/GATE_tree.txt || stop "a0 vs p2_18"
grep -q "vs d2aa: 26 of 2639 changed" $S/GATE_tree.txt || stop "a0 control d2aa"
grep -q "vs d2t: 863 of 2639 changed" $S/GATE_tree.txt || stop "a0 control d2t"
say "a0 PASS; face gate"
$PY $W/face_gate.py run p2_22 > $S/face.log 2>&1 || stop "face_gate run ($S/face.log)"
$PY $W/face_gate.py cmp p2_21 p2_22 > $S/GATE_face_cmp.txt 2>&1
tail -3 $S/GATE_face_cmp.txt | tee -a $LOG
grep -q "VERDICT: IDENTICAL" $S/GATE_face_cmp.txt || stop "a1 ($S/GATE_face_cmp.txt)"
say "beam_void_closure (the sweep's own output) + unit/beam"
$PY $U/beam_void_closure/test $U/beam_void_closure/output > $S/A2_beam_void_closure.log 2>&1
{ grep "^\[PASS\] T\|^\[FAIL\] T" $S/A2_beam_void_closure.log
  nbad=$(grep "^\[FAIL\] T" $S/A2_beam_void_closure.log | grep -vc "^\[FAIL\] T3 ")
  npass=$(grep -c "^\[PASS\] T[124]" $S/A2_beam_void_closure.log)
  if [ "$nbad" = 0 ] && [ "$npass" -ge 7 ] && grep -q "^\[FAIL\] T3 " $S/A2_beam_void_closure.log
  then echo "a2 beam_void_closure PASS (T1a-d, T2a-b, T4 PASS: $npass rows; T3 FAIL by design, unchanged)"; else echo "a2 beam_void_closure FAIL ($nbad non-T3 FAIL, $npass PASS)"; fi; } > $S/A2_beam_void_closure.txt
cut -c1-160 $S/A2_beam_void_closure.txt | tee -a $LOG
grep -q "^a2 beam_void_closure PASS" $S/A2_beam_void_closure.txt || stop "a2 beam_void_closure"
rm -rf $U/beam/output
mpirun --oversubscribe --bind-to none -np 4 $BIN $U/beam/input > $S/A2_beam_run.log 2>&1 || stop "a2 unit/beam run"
$PY $U/beam/test $U/beam/output > $S/A2_beam.log 2>&1 || { cat $S/A2_beam.log >> $LOG; stop "a2 unit/beam"; }
tail -2 $S/A2_beam.log | tee -a $LOG

say "ord 3: C1 vs P2-20's Z1 (report: expected to move)"
$PY $S/run.py C1 >> $S/RUN.log 2>&1 || stop "C1 run"
$PY $S/c1_compare.py > $S/C1_vs_Z1.txt 2>&1; cat $S/C1_vs_Z1.txt | tee -a $LOG
( cd $S/legs/C1/run && ls -d *cell *node 2>/dev/null | grep -v -e '^00000cell$' -e '^00075cell$' | xargs rm -rf )
say "A (table, true spot)"
$PY $S/run.py A >> $S/RUN.log 2>&1 || stop "A run"
$PY $ZT $S/legs/A/run > $S/TEST_A.log 2>&1; rz=$?
grep "  PASS\|  FAIL\|b1r\|^PASS\|^FAIL" $S/TEST_A.log | cut -c1-230 | tee -a $LOG
[ $rz = 0 ] || say "  A: the Zhang test exit $rz (see TEST_A.log)"
say "X (A explicit dt 0.03)"
$PY $S/run.py X >> $S/RUN.log 2>&1 || stop "X run"
$PY $ZT $S/legs/A/run --explicit $S/legs/X/run > $S/TEST_X.log 2>&1
grep "b3" $S/TEST_X.log | cut -c1-230 | tee -a $LOG
prune A; prune X
say "F (flat control, true spot)"
$PY $S/run.py F >> $S/RUN.log 2>&1 || stop "F run"
$PY $ZT $S/legs/F/run > $S/TEST_F.log 2>&1
grep "  FAIL" $S/TEST_F.log | cut -c1-200 | tee -a $LOG
prune F

say "ord 4: RESULTS.txt"
$PY $S/reduce.py > $S/RESULTS.txt 2>&1
sed -n '/^legs/,$p' $S/RESULTS.txt | cut -c1-240 | tee -a $LOG
du -sh $S/legs | tee -a $LOG
say "CHAIN DONE"
