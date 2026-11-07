#!/bin/bash
# P2-24 in the packet's run order (U, gates, M0, M100, M0F, M1); STOPS at the first failed gate. Log: CHAIN.log.
#   ord 1  z1 unit/face_melt, z2 unit/face_vapor (U3 on bin_pre = the P2-23 build), z3 unit/retired_keys +
#          unit/latent_heat (z0 is offline: Z0.txt, written before the chain)
#   ord 2  a0 tree gate vs p2_18 (controls d2aa 26, d2t 863); a1 face gate run p2_24, cmp p2_23 p2_24 (the record) and
#          a1.py (the packet's L9 scoring); a2 beam_void_closure (the sweep's own output), unit/beam, unit/beam_tau,
#          unit/face_beam
#   ord 3  M0 (b1, b2 gate the rest), check_shared on the deck (M0's metadata), M100, M0F, M1; the full test; each leg
#          pruned after it to its onset / 600 s / 1,580 s / last plotfiles with an md5 manifest of every plotfile first
#   ord 4  RESULTS.txt
# FROM=a1 bash chain.sh  resumes at the face gate (ord 1 and a0 re-checked from their logs, not re-run)
set -u
ROOT=/Users/tzetze20/amr_tools/alamo; cd $ROOT
S=tests/MMWSpalling/studies/p2_24_melt_vapour
U=tests/MMWSpalling/unit
W=tests/MMWSpalling/studies/tree_sweep_0922/witness
TC=tests/MMWSpalling/studies/p2_9_face_removal/tree_compare.py
OT=tests/MMWSpalling/validation/oglesby/test
CS=tests/MMWSpalling/validation/shared/check_shared.py
BIN=bin/mmwspalling-3d-g++
PY=/Users/tzetze20/Desktop/code/.venv/bin/python; LOG=$S/CHAIN.log
say()  { echo "$(date +%T) $*" | tee -a $LOG; }
stop() { say "STOP: $*"; exit 1; }
caffeinate -dis -w $$ &
B=$(shasum -a 256 $BIN | cut -c1-64)
say "P2-24 chain (binary ${B:0:16})"
grep -q "$B" $S/BINARY.txt || stop "bin/ is not the P2-24 build"

if [ "${FROM:-}" = a1 ]; then
  # a resume after the a1 stop of run 1 (run1_a1stop/: the BEAM_S fixture lost its buckle log to the #19 deck's
  # buckle_log_csv = 0 and turned it back on). ord 1 and a0 ran on THIS binary in run 1; their logs are re-checked, not
  # re-run (one full tree gate per packet).
  say "resume at a1: ord 1 and a0 from run 1 (run1_a1stop/CHAIN.log), the same binary"
  for f in Z1_face_melt.log Z2_face_vapor.log; do grep -q "^PASS$" $S/$f || stop "run 1's $f"; done
  grep -q "^PASS" $S/Z3_retired_keys.log && grep -q "PASS" $S/Z3_latent_heat.log || stop "run 1's z3"
  grep -q "vs p2_18: 0 of 2639 changed, 0 only in new, 0 only in ref" $S/GATE_tree.txt || stop "run 1's a0"
  cat $S/GATE_tree.txt >> $LOG
else
say "ord 1: z1 unit/face_melt"
$PY $U/face_melt/test --bin_pre $S/bin_pre/mmwspalling-3d-g++ > $S/Z1_face_melt.log 2>&1 || { tail -15 $S/Z1_face_melt.log >> $LOG; stop "z1"; }
grep "PASS  \|FAIL  " $S/Z1_face_melt.log | cut -c1-200 | tee -a $LOG
say "z2 unit/face_vapor"
$PY $U/face_vapor/test --bin_pre $S/bin_pre/mmwspalling-3d-g++ > $S/Z2_face_vapor.log 2>&1 || { tail -20 $S/Z2_face_vapor.log >> $LOG; stop "z2"; }
grep "PASS  \|FAIL  \|U3 H_vap" $S/Z2_face_vapor.log | cut -c1-200 | tee -a $LOG
say "z3 unit/retired_keys, unit/latent_heat"
$PY $U/retired_keys/test > $S/Z3_retired_keys.log 2>&1 || { tail -25 $S/Z3_retired_keys.log >> $LOG; stop "z3 retired_keys"; }
tail -1 $S/Z3_retired_keys.log | cut -c1-160 | tee -a $LOG
$PY $U/latent_heat/test > $S/Z3_latent_heat.log 2>&1 || { tail -25 $S/Z3_latent_heat.log >> $LOG; stop "z3 latent_heat"; }
tail -1 $S/Z3_latent_heat.log | cut -c1-160 | tee -a $LOG

say "ord 2: a0 tree gate"
bash $W/guarded_post_all_par.sh p2_24 > $S/tree.log 2>&1
[ $? = 2 ] && stop "the binary moved during the sweep"
$PY $TC p2_24 p2_18 d2aa d2t > $S/GATE_tree.txt 2>&1; cat $S/GATE_tree.txt >> $LOG
grep -q "vs p2_18: 0 of 2639 changed, 0 only in new, 0 only in ref" $S/GATE_tree.txt || stop "a0 vs p2_18"
grep -q "vs d2aa: 26 of 2639 changed" $S/GATE_tree.txt || stop "a0 control d2aa"
grep -q "vs d2t: 863 of 2639 changed" $S/GATE_tree.txt || stop "a0 control d2t"
say "a0 PASS"
fi
say "a1 face gate"
$PY $W/face_gate.py run p2_24 > $S/face.log 2>&1 || stop "face_gate run ($S/face.log)"
$PY $W/face_gate.py cmp p2_23 p2_24 > $S/GATE_face_cmp.txt 2>&1
$PY $S/a1.py > $S/A1.txt 2>&1
cut -c1-240 $S/A1.txt | tee -a $LOG
grep -q "^a1 PASS" $S/A1.txt || stop "a1 ($S/A1.txt, $S/GATE_face_cmp.txt)"
say "a2 beam_void_closure (the sweep's own output), unit/beam, unit/beam_tau, unit/face_beam"
$PY $U/beam_void_closure/test $U/beam_void_closure/output > $S/A2_beam_void_closure.log 2>&1
{ grep "^\[PASS\] T\|^\[FAIL\] T" $S/A2_beam_void_closure.log
  nbad=$(grep "^\[FAIL\] T" $S/A2_beam_void_closure.log | grep -vc "^\[FAIL\] T3 ")
  npass=$(grep -c "^\[PASS\] T[124]" $S/A2_beam_void_closure.log)
  if [ "$nbad" = 0 ] && [ "$npass" -ge 7 ] && grep -q "^\[FAIL\] T3 " $S/A2_beam_void_closure.log
  then echo "a2 beam_void_closure PASS ($npass rows; T3 FAIL by design)"; else echo "a2 beam_void_closure FAIL ($nbad non-T3 FAIL)"; fi; } > $S/A2_beam_void_closure.txt
tail -1 $S/A2_beam_void_closure.txt | tee -a $LOG
grep -q "^a2 beam_void_closure PASS" $S/A2_beam_void_closure.txt || stop "a2 beam_void_closure"
rm -rf $U/beam/output
mpirun --oversubscribe --bind-to none -np 4 $BIN $U/beam/input > $S/A2_beam_run.log 2>&1 || stop "a2 unit/beam run"
$PY $U/beam/test $U/beam/output > $S/A2_beam.log 2>&1 || { cat $S/A2_beam.log >> $LOG; stop "a2 unit/beam"; }
grep "integral" $S/A2_beam.log | tee -a $LOG
$PY $U/beam_tau/test > $S/A2_beam_tau.log 2>&1 || { tail -20 $S/A2_beam_tau.log >> $LOG; stop "a2 unit/beam_tau"; }
tail -1 $S/A2_beam_tau.log | tee -a $LOG
$PY $U/face_beam/test > $S/A2_face_beam.log 2>&1 || { tail -20 $S/A2_face_beam.log >> $LOG; stop "a2 unit/face_beam"; }
grep "PASS  \|FAIL  \|z1r" $S/A2_face_beam.log | cut -c1-200 | tee -a $LOG

prune() {   # keep the onset, 600 s, 1,580 s and last cell plotfiles; md5 manifest of every plotfile file first
  $PY - "$S/legs/$1" <<'EOF'
import os, sys, shutil, hashlib, math
d = sys.argv[1]; run = os.path.join(d, "run")
L = open(os.path.join(run, "thermo.dat")).read().splitlines()
h = [c.strip() for c in L[0].split("\t")]
rows = [[float(x) for x in r.split("\t")] for r in L[1:] if r.strip()]
t_on = next((r[0] for r in rows if math.isfinite(r[h.index("bc_z_melt")])), None)
man = []
for p in sorted(os.listdir(run)):
    if p.endswith(("cell", "node")) and os.path.isdir(os.path.join(run, p)):
        for dp, _, fs in os.walk(os.path.join(run, p)):
            for fn in sorted(fs):
                q = os.path.join(dp, fn)
                man.append("%s  %s" % (hashlib.md5(open(q, "rb").read()).hexdigest(), os.path.relpath(q, run)))
open(os.path.join(d, "MANIFEST.md5"), "w").write("\n".join(man) + "\n")
pfs = sorted(p for p in os.listdir(run) if p.endswith("cell"))
def ptime(pf):
    H = open(os.path.join(run, pf, "Header")).read().split("\n")
    return float(H[2 + int(H[1]) + 1])
keep = {pfs[-1]}
for t in [x for x in (t_on, 600.0, 1580.0) if x is not None]:
    after = [pf for pf in pfs if ptime(pf) >= t - 1e-9]
    if after: keep.add(after[0])
gone = [pf for pf in pfs if pf not in keep]
for pf in gone: shutil.rmtree(os.path.join(run, pf))
for pf in [p for p in os.listdir(run) if p.endswith("node")]: shutil.rmtree(os.path.join(run, pf))
open(os.path.join(d, "PRUNED.txt"), "w").write("pruned after its test: kept %s (the first at/after onset %s s, 600 s, "
    "1580 s; the last); removed %d cell plotfiles and the node plotfiles; MANIFEST.md5 lists every plotfile file before "
    "the prune\n" % (" ".join(sorted(keep)), t_on, len(gone)))
EOF
}

say "ord 3: M0 (b1, b2 gate M100, M0F, M1)"
$PY $S/run.py M0 >> $S/RUN.log 2>&1 || stop "M0 run"
$PY $OT $S/legs/M0/run --onset > $S/TEST_M0.log 2>&1; rb=$?
grep "b1\|b2\|b1r" $S/TEST_M0.log | cut -c1-230 | tee -a $LOG
[ $rb = 0 ] || stop "b1/b2 ($S/TEST_M0.log)"
$PY $CS tests/MMWSpalling/validation/oglesby/input --run $S/legs/M0/run > $S/CHECK_SHARED.txt 2>&1
tail -6 $S/CHECK_SHARED.txt | cut -c1-200 | tee -a $LOG
for l in M100 M0F M1; do
  say "  $l"
  $PY $S/run.py $l >> $S/RUN.log 2>&1 || say "  $l rc != 0 (legs/$l/run.log)"
done
$PY $OT $S/legs/M0/run --onset --alt $S/legs/M100/run --flat $S/legs/M0F/run --fine $S/legs/M1/run > $S/TEST_all.log 2>&1
sed -n '/^m1 /,$p' $S/TEST_all.log | cut -c1-240 | tee -a $LOG
for l in M0 M100 M0F M1; do prune $l; done

say "ord 4: RESULTS.txt"
$PY $S/reduce.py > $S/RESULTS.txt 2>&1
sed -n '/^e0 VERDICT/,$p' $S/RESULTS.txt | cut -c1-240 | tee -a $LOG
du -sh $S/legs | tee -a $LOG
say "CHAIN DONE"
