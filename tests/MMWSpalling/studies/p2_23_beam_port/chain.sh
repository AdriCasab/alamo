#!/bin/bash
# P2-23 in the packet's run order; STOPS at the first failed gate. Log: CHAIN.log.
#   ord 1  z1, z1r, z2 unit/face_beam (U1, U2; U3 on bin_pre = the P2-22 build), z3 unit/retired_keys
#   ord 2  a0 tree gate vs p2_18 (controls d2aa 26, d2t 863); a1 face gate run p2_23, cmp p2_22 p2_23 (the 13
#          legs identical, BEAM_S only in new); a2 beam_void_closure (the sweep's own output), unit/beam, unit/beam_tau
#   ord 3  B0 (b1 gates the rest), P, B0F, B1 (the fine leg, last); the full test; each leg pruned after it to its
#          first-fire and last plotfiles with an md5 manifest of every plotfile file first
#   ord 4  RESULTS.txt
set -u
ROOT=/Users/tzetze20/amr_tools/alamo; cd $ROOT
S=tests/MMWSpalling/studies/p2_23_beam_port
U=tests/MMWSpalling/unit
W=tests/MMWSpalling/studies/tree_sweep_0922/witness
TC=tests/MMWSpalling/studies/p2_9_face_removal/tree_compare.py
OT=tests/MMWSpalling/validation/oglesby/test
BIN=bin/mmwspalling-3d-g++
PY=/Users/tzetze20/Desktop/code/.venv/bin/python; LOG=$S/CHAIN.log
say()  { echo "$(date +%T) $*" | tee -a $LOG; }
stop() { say "STOP: $*"; exit 1; }
caffeinate -dis -w $$ &
B=$(shasum -a 256 $BIN | cut -c1-64)
say "P2-23 chain (binary ${B:0:16})"
grep -q "$B" $S/BINARY.txt || stop "bin/ is not the P2-23 build"

say "ord 1: z1, z1r, z2 unit/face_beam"
$PY $U/face_beam/test --bin_pre $S/bin_pre/mmwspalling-3d-g++ > $S/Z1_face_beam.log 2>&1 || { tail -15 $S/Z1_face_beam.log >> $LOG; stop "z1/z2"; }
grep "PASS  \|FAIL  \|z1r" $S/Z1_face_beam.log | cut -c1-200 | tee -a $LOG
say "z3 unit/retired_keys"
$PY $U/retired_keys/test > $S/Z3_retired_keys.log 2>&1 || { tail -25 $S/Z3_retired_keys.log >> $LOG; stop "z3"; }
tail -1 $S/Z3_retired_keys.log | cut -c1-160 | tee -a $LOG

say "ord 2: a0 tree gate"
bash $W/guarded_post_all_par.sh p2_23 > $S/tree.log 2>&1
[ $? = 2 ] && stop "the binary moved during the sweep"
$PY $TC p2_23 p2_18 d2aa d2t > $S/GATE_tree.txt 2>&1; cat $S/GATE_tree.txt >> $LOG
grep -q "vs p2_18: 0 of 2639 changed, 0 only in new, 0 only in ref" $S/GATE_tree.txt || stop "a0 vs p2_18"
grep -q "vs d2aa: 26 of 2639 changed" $S/GATE_tree.txt || stop "a0 control d2aa"
grep -q "vs d2t: 863 of 2639 changed" $S/GATE_tree.txt || stop "a0 control d2t"
say "a0 PASS; a1 face gate"
$PY $W/face_gate.py run p2_23 > $S/face.log 2>&1 || stop "face_gate run ($S/face.log)"
$PY $W/face_gate.py cmp p2_22 p2_23 > $S/GATE_face_cmp.txt 2>&1
$PY - $W > $S/A1.txt 2>&1 <<'EOF'
import sys, os
W = sys.argv[1]
rd = lambda l: dict(reversed(ln.rstrip("\n").split("  ", 1)) for ln in open(os.path.join(W, "face_%s.md5" % l)))
a, b = rd("p2_22"), rd("p2_23")
only_a = sorted(set(a) - set(b)); only_b = sorted(set(b) - set(a))
cmpt = open(os.path.join(W, "..", "..", "p2_23_beam_port", "GATE_face_cmp.txt")).read()
real = int(cmpt.split("VERDICT:")[1].split("(")[1].split()[0])
ok = real == 0 and not only_a and only_b and all(p.startswith("BEAM_S/") for p in only_b)
print("a1: %d substantive difference(s) on the common paths (thermo timers excluded); only in p2_22: %d; only in p2_23: "
      "%d, all under BEAM_S/: %s (%s)" % (real, len(only_a), len(only_b), all(p.startswith("BEAM_S/") for p in only_b),
                                          ", ".join(only_b[:6])))
print("a1 %s" % ("PASS" if ok else "FAIL"))
EOF
cat $S/A1.txt | cut -c1-240 | tee -a $LOG
grep -q "^a1 PASS" $S/A1.txt || stop "a1 ($S/A1.txt, $S/GATE_face_cmp.txt)"
say "a2 beam_void_closure (the sweep's own output), unit/beam, unit/beam_tau"
$PY $U/beam_void_closure/test $U/beam_void_closure/output > $S/A2_beam_void_closure.log 2>&1
{ grep "^\[PASS\] T\|^\[FAIL\] T" $S/A2_beam_void_closure.log
  nbad=$(grep "^\[FAIL\] T" $S/A2_beam_void_closure.log | grep -vc "^\[FAIL\] T3 ")
  npass=$(grep -c "^\[PASS\] T[124]" $S/A2_beam_void_closure.log)
  if [ "$nbad" = 0 ] && [ "$npass" -ge 7 ] && grep -q "^\[FAIL\] T3 " $S/A2_beam_void_closure.log
  then echo "a2 beam_void_closure PASS ($npass rows; T3 FAIL by design)"; else echo "a2 beam_void_closure FAIL ($nbad non-T3 FAIL)"; fi; } > $S/A2_beam_void_closure.txt
tail -1 $S/A2_beam_void_closure.txt | tee -a $LOG
grep -q "^a2 beam_void_closure PASS" $S/A2_beam_void_closure.txt || stop "a2 beam_void_closure"
[ -d $U/beam/output ] && [ ! -d $S/unit_beam_output_p2_22 ] && cp -cR $U/beam/output $S/unit_beam_output_p2_22
rm -rf $U/beam/output
mpirun --oversubscribe --bind-to none -np 4 $BIN $U/beam/input > $S/A2_beam_run.log 2>&1 || stop "a2 unit/beam run"
$PY $U/beam/test $U/beam/output > $S/A2_beam.log 2>&1 || { cat $S/A2_beam.log >> $LOG; stop "a2 unit/beam"; }
grep "integral" $S/A2_beam.log | tee -a $LOG
$PY $U/beam_tau/test > $S/A2_beam_tau.log 2>&1 || { tail -20 $S/A2_beam_tau.log >> $LOG; stop "a2 unit/beam_tau"; }
tail -1 $S/A2_beam_tau.log | tee -a $LOG

prune() {   # keep the first-fire and the last cell plotfile; md5 manifest of every plotfile file first
  $PY - "$S/legs/$1" <<'EOF'
import os, sys, shutil, hashlib
d = sys.argv[1]; run = os.path.join(d, "run")
L = open(os.path.join(run, "thermo.dat")).read().splitlines()
h = [c.strip() for c in L[0].split("\t")]
rows = [[float(x) for x in r.split("\t")] for r in L[1:] if r.strip()]
t_fire = next((r[0] for r in rows if r[0] > 0 and r[h.index("bk_n_eval")] > 0), None)
man = []
for p in sorted(os.listdir(run)):
    if p.endswith(("cell", "node")) and os.path.isdir(os.path.join(run, p)):
        for dp, _, fs in os.walk(os.path.join(run, p)):
            for fn in sorted(fs):
                q = os.path.join(dp, fn)
                man.append("%s  %s" % (hashlib.md5(open(q, "rb").read()).hexdigest(), os.path.relpath(q, run)))
open(os.path.join(d, "MANIFEST.md5"), "w").write("\n".join(man) + "\n")
pfs = sorted(p for p in os.listdir(run) if p.endswith("cell"))
ptime = lambda pf: float(open(os.path.join(run, pf, "Header")).read().split("\n")[2 + int(open(os.path.join(run, pf, "Header")).read().split("\n")[1]) + 1])
keep = {pfs[-1]}
if t_fire is not None:
    after = [pf for pf in pfs if ptime(pf) >= t_fire - 1e-9]
    if after: keep.add(after[0])
gone = [pf for pf in pfs if pf not in keep]
for pf in gone: shutil.rmtree(os.path.join(run, pf))
for pf in [p for p in os.listdir(run) if p.endswith("node")]: shutil.rmtree(os.path.join(run, pf))
open(os.path.join(d, "PRUNED.txt"), "w").write("pruned after its test: kept %s (first fire %s s; the last); removed %d "
    "cell plotfiles and the node plotfiles; MANIFEST.md5 lists every plotfile file before the prune\n"
    % (" ".join(sorted(keep)), t_fire, len(gone)))
EOF
}

say "ord 3: B0 (b1 gates P, B0F, B1)"
$PY $S/run.py B0 >> $S/RUN.log 2>&1 || stop "B0 run"
$PY $OT $S/legs/B0/run > $S/TEST_B0.log 2>&1; rb=$?
grep "b1\|b2 \|b2r" $S/TEST_B0.log | cut -c1-230 | tee -a $LOG
[ $rb = 0 ] || stop "b1 ($S/TEST_B0.log): the scan never reached a beam line"
for l in P B0F B1; do
  say "  $l"
  $PY $S/run.py $l >> $S/RUN.log 2>&1 || say "  $l rc != 0 (legs/$l/run.log)"
done
$PY $OT $S/legs/B0/run --control $S/legs/P/run --flat $S/legs/B0F/run --fine $S/legs/B1/run > $S/TEST_all.log 2>&1
sed -n '/^b2 /,$p' $S/TEST_all.log | cut -c1-240 | tee -a $LOG
for l in B0 P B0F B1; do prune $l; done

say "ord 4: RESULTS.txt"
$PY $S/reduce.py > $S/RESULTS.txt 2>&1
sed -n '/^e0 VERDICT/,$p' $S/RESULTS.txt | cut -c1-240 | tee -a $LOG
du -sh $S/legs | tee -a $LOG
say "CHAIN DONE"
