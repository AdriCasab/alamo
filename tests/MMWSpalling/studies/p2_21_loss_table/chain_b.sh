#!/bin/bash
# P2-21 part B, the packet's run order 3-4 (after chain_a.sh passed ord 1-2); STOPS at the first failed gate.
#   ord 3  T1b (b1, b1r, b2), T1x (b3: T1b's onset read first; within 10 %), then the reports T1a, T1c, T0, M2,
#          M3.5, T1d, T2; each leg's test log, then the leg pruned to its onset plotfile and its last
#   ord 4  RESULTS.txt (reduce.py: c5-c10, the scan c8, e0)
# Every T-leg runs at Fig. 34's MEASURED incident power (run.py ramp_schedule; the user's decision 2026-11-06).
set -u
ROOT=/Users/tzetze20/amr_tools/alamo; cd $ROOT
S=tests/MMWSpalling/studies/p2_21_loss_table
ZT=tests/MMWSpalling/validation/zhang/zhang_oglesby/test
BIN=bin/mmwspalling-3d-g++
PY=/Users/tzetze20/Desktop/code/.venv/bin/python; LOG=$S/CHAIN.log
say()  { echo "$(date +%T) $*" | tee -a $LOG; }
stop() { say "STOP: $*"; exit 1; }
caffeinate -dis -w $$ &
B=$(shasum -a 256 $BIN | cut -c1-64)
say "P2-21 chain B (binary ${B:0:16})"
grep -q "$B" $S/BINARY.txt || stop "bin/ is not the P2-21 build"
grep -q "CHAIN A DONE" $LOG || stop "chain A has not passed"

prune() {   # keep the first plotfile at or after the onset time (if any) and the last; record it
  $PY - "$S/legs/$1" <<'EOF'
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
EOF
}

say "ord 3: T1b (b1, b1r, b2)"
$PY $S/run.py T1b >> $S/RUN.log 2>&1 || stop "T1b run"
$PY $ZT $S/legs/T1b/run > $S/TEST_T1b.log 2>&1; rz=$?
grep "  PASS\|  FAIL\|b1r\|^PASS\|^FAIL" $S/TEST_T1b.log | cut -c1-230 | tee -a $LOG
[ $rz = 0 ] || stop "b1/b2 ($S/TEST_T1b.log)"
say "T1x (b3): stop_time from T1b's onset"
$PY $S/run.py T1x >> $S/RUN.log 2>&1 || stop "T1x run"
$PY $ZT $S/legs/T1b/run --explicit $S/legs/T1x/run > $S/TEST_T1x.log 2>&1
grep "b3" $S/TEST_T1x.log | tee -a $LOG
grep -q "PASS  b3" $S/TEST_T1x.log || stop "b3 ($S/TEST_T1x.log)"
prune T1b; prune T1x
say "b3 PASS; the reports: T1a T1c T0 M2 M3.5 T1d T2"
for l in T1a T1c T0 M2 M3.5 T1d T2; do
  $PY $S/run.py $l >> $S/RUN.log 2>&1 || say "  $l rc != 0 (see legs/$l/run.log)"
  $PY $ZT $S/legs/$l/run > $S/TEST_$l.log 2>&1
  grep "  FAIL" $S/TEST_$l.log | cut -c1-200 | tee -a $LOG
  [ $l = T2 ] || prune $l   # T2 keeps its 4 plotfiles (c9 reads 50 / 100 / 150 s)
  say "  $l done"
done
$PY $ZT $S/legs/T1c/run --scan m1=$S/legs/T0/run m2=$S/legs/M2/run m3.5=$S/legs/M3.5/run m8.88=$S/legs/T1c/run \
  > $S/TEST_scan.log 2>&1

say "ord 4: RESULTS.txt"
$PY $S/reduce.py > $S/RESULTS.txt 2>&1
sed -n '/^c5-c7/,$p' $S/RESULTS.txt | cut -c1-240 | tee -a $LOG
du -sh $S/legs | tee -a $LOG
say "CHAIN B DONE"
