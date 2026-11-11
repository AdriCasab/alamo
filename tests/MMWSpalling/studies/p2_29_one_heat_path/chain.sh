#!/bin/bash
# P2-29 chain -- one call per stage, serial, held awake (caffeinate -w on this shell); log CHAIN.log.
#   bash chain.sh 1h  stage 1's Hu section after the user's trim (7 new decks; see the block below)
#   bash chain.sh 1   stage 1 (p2_29a, additive): z0 -> build -> a1 face gate -> a2 latent_heat -> a7/m2 fixtures ->
#                     the Zhang legs ZO, ZF, ZQ, ZQ1, ZL -> a6, a3, a4, a5, m3, m8 -> Hu old (bin_pre) + new (bin_1),
#                     9 decks each -> a8, m4-m6 -> STOP for the user.
# A failed GATE skips every later simulation leg (the reductions of runs already made still report).
set -u
STAGE=${1:?usage: chain.sh <stage>}
ROOT=/Users/tzetze20/amr_tools/alamo; cd $ROOT
S=tests/MMWSpalling/studies/p2_29_one_heat_path
W=tests/MMWSpalling/studies/tree_sweep_0922/witness
H=tests/MMWSpalling/validation/hu
Z=tests/MMWSpalling/validation/zhang/zhang_oglesby
BIN=bin/mmwspalling-3d-g++
BPRE=$S/bin_pre/mmwspalling-3d-g++
PY=/Users/tzetze20/Desktop/code/.venv/bin/python
LOG=$S/CHAIN.log
MPI="mpirun --oversubscribe --bind-to none -np 4"
say()  { echo "$(date +%T) [s$STAGE] $*" | tee -a $LOG; }
T0=$(date "+%Y-%m-%d %H:%M:%S"); S0=$SECONDS
FAILED=0
finish() {
  n=$(pmset -g log | awk -v s="$T0" '($4=="Sleep"||$4=="Wake"||$4=="DarkWake") && substr($0,1,19) >= s' \
      | tee $S/SLEEP_$STAGE.txt | wc -l | tr -d ' ')
  say "wall $((SECONDS - S0)) s since $T0; sleep/wake events in this window: $n (SLEEP_$STAGE.txt)"
  pmset -g batt | grep -q "AC Power" || say "NOTE: on BATTERY now"
}
gate() {   # name rc log
  if [ "$2" = 0 ]; then say "  $1 rc 0: $(grep -E 'RESULT|VERDICT|^PASS|^FAIL' $3 | tail -1 | cut -c1-200)"
  else say "  GATE FAILED: $1 (rc $2, $3): $(grep -E 'RESULT|VERDICT|FAIL' $3 | tail -1 | cut -c1-200)"; FAILED=1; fi
}
report() { say "  $1 (report): $(grep -E 'RESULT' $2 | tail -1 | cut -c1-200)"; }
# one simulation leg: name binary deck plot_file [argv ...]; argv + binary recorded (SR 8)
leg() {
  local name=$1 bin=$2 deck=$3 pf=$4; shift 4
  if [ $FAILED = 1 ]; then say "  SKIPPED (a gate failed): leg $name"; return 1; fi
  mkdir -p $(dirname $pf)
  echo "$MPI $bin $deck $* plot_file=$pf" > $pf.argv.sh     # <run>.argv.sh / .BINARY.txt / .log beside the run
  shasum -a 256 $bin > $pf.BINARY.txt
  local t=$SECONDS
  $MPI $bin $deck "$@" plot_file=$pf > $pf.log 2>&1; local rc=$?
  say "  leg $name rc $rc ($((SECONDS - t)) s)"
  [ $rc = 0 ] || { FAILED=1; say "  GATE FAILED: leg $name did not run ($pf.log)"; }
  return $rc
}
caffeinate -dis -w $$ &

if [ $STAGE = 2 ]; then
  # P2-29 stage 2 (the user's go, 2026-11-10): the separate heat kernels deleted (stage2_src.py), the phase decks
  # renamed and the old ones / the AMR tests moved (stage2_tests.py) -- both applied before this chain
  env EIGEN=$ROOT/ext CPLUS_INCLUDE_PATH=/opt/homebrew/include \
      LIBRARY_PATH=/opt/homebrew/lib:/opt/homebrew/Cellar/gcc/15.2.0_1/lib/gcc/current make -j8 > $S/build_2b.log 2>&1 \
    || { say "STOP: the build failed ($S/build_2b.log)"; finish; exit 1; }
  mkdir -p $S/bin_2 && cp $BIN $S/bin_2/
  { echo "P2-29 stage 2 (p2_29): the separate heat kernels deleted; retired keys; foot_rule max/mean deleted"
    shasum -a 256 $BIN; shasum -a 256 src/Integrator/MMWSpalling.H src/Numeric/Material/Table.H; } > $S/BINARY_2.txt
  say "P2-29 stage 2 chain (binary $(shasum -a 256 $BIN | cut -c1-16), bin_2/); $(pmset -g batt | head -1)"
  # b1: the face gate
  $PY $W/face_gate.py run p2_29 > $S/face_2.log 2>&1
  $PY $W/face_gate.py cmp p2_29a p2_29 > $S/GATE_b1.txt 2>&1
  grep -q "^VERDICT: IDENTICAL" $S/GATE_b1.txt; gate b1 $? $S/GATE_b1.txt
  # b2: ZQ and hu_conduction re-run on bin_2 (the stage-1 argv; the renamed decks)
  leg ZQ $BIN $Z/input $S/legs2/ZQ/run energy_ledger.enabled=1
  leg hu_conduction/granite2 $BIN $H/hu_conduction/input_granite2 $S/legs2/hu_conduction/granite2
  leg hu_conduction/sandstone2 $BIN $H/hu_conduction/input_sandstone2 $S/legs2/hu_conduction/sandstone2
  $PY $S/b2.py > $S/K2_b2.log 2>&1; gate b2 $? $S/K2_b2.log
  # b3: retired_keys (R12 RED on the stage-2 pre-image) + the locks + interlocks
  $PY tests/MMWSpalling/unit/retired_keys/test --pre $S/pre2_flat --red R12 > $S/K2_retired_keys.log 2>&1
  gate "b3 retired_keys" $? $S/K2_retired_keys.log
  for t in lateral_flux_lock buckle_lock interlocks; do
    $PY tests/MMWSpalling/unit/$t/test > $S/K2_$t.log 2>&1; gate "b3 $t" $? $S/K2_$t.log
  done
  # b4: source assertions (R12 carries them; restated here as a grep)
  n=$(grep -cE "AdvanceConstAlpha|AdvanceMaterial|SolverKind|material_table|therm_alpha" src/Integrator/MMWSpalling.H)
  h=$(grep -rnE "[^0] \* kc \* \(sp_T - Tc\) \* inv_dz2|[^*] kc \* \(sp_T - Tc\) \* inv_dz2" src/ | grep -vc "2.0 \* kc" )
  say "  b4 source: deleted names $n hit(s); a half-conductance term without its 2.0: $h hit(s)"
  [ "$n" = 0 ] && [ "$h" = 0 ] || { FAILED=1; say "  GATE FAILED: b4"; }
  # per-edit checks on the renamed decks (WHICH_TESTS: the phase law, the one heat path, the feet rule)
  for t in latent_heat beam_tau face_feet; do
    $PY tests/MMWSpalling/unit/$t/test > $S/K2_$t.log 2>&1; gate "check $t" $? $S/K2_$t.log
  done
  finish
  if [ $FAILED = 0 ]; then say "STAGE 2 DONE -- b1-b4 hold"
  else say "STAGE 2 DONE WITH A FAILED GATE -- see above"; fi
  exit 0
fi

if [ $STAGE = 1b ]; then
  # the user (2026-11-10): the breakage decks SINGLE LEVEL on their 40^3 base grid -- AMR + mechanics fails on the
  # production path (MLMG at the first regrid, pre-existing on bin_pre); a8 is dropped (no level 1)
  B=$(shasum -a 256 $BIN | cut -c1-64)
  grep -q "$B" $S/BINARY_1.txt || { say "STOP: bin/ is not BINARY_1.txt's binary"; exit 1; }
  say "P2-29 stage 1, the breakage decks single level (binary ${B:0:16}); $(pmset -g batt | head -1)"
  for c in granite2_37 sandstone2_37 sandstone2_89; do
    leg "hu_new/hu_breakage_index/$c" $BIN $H/hu_breakage_index/input_${c}_phase $S/legs/hu_new/hu_breakage_index/$c
  done
  if [ -d $S/legs/hu_new/hu_breakage_index/sandstone2_89 ]; then
    $PY $S/hu_compare.py report > $S/K1_m4_m6.log 2>&1; report "m4-m6" $S/K1_m4_m6.log
  fi
  finish
  if [ $FAILED = 0 ]; then say "STAGE 1 (breakage) DONE -- every leg ran; STOP for the user (stage 2 needs the go)"
  else say "STAGE 1 (breakage) DONE WITH A FAILED LEG -- see above"; fi
  exit 0
fi

if [ $STAGE = 1h ]; then
  # the user's trim (2026-11-10, mid-chain): the Hu section again, the NEW loop only (7 decks; hu_thermoelastic_amr
  # retired, not ported; no old side: m4/m5 read the archived 09-25 outputs, m6 Hu's data)
  B=$(shasum -a 256 $BIN | cut -c1-64)
  grep -q "$B" $S/BINARY_1.txt || { say "STOP: bin/ is not BINARY_1.txt's binary"; exit 1; }
  say "P2-29 stage 1, the Hu section after the trim (binary ${B:0:16}); $(pmset -g batt | head -1)"
  for d in hu_conduction/input_granite2 hu_conduction/input_sandstone2 hu_thermoelastic/input_granite2 \
           hu_thermoelastic/input_sandstone2 hu_breakage_index/input_granite2_37 hu_breakage_index/input_sandstone2_37 \
           hu_breakage_index/input_sandstone2_89; do
    t=$(dirname $d); c=$(basename $d); c=${c#input_}
    leg "hu_new/$t/$c" $BIN $H/${d}_phase $S/legs/hu_new/$t/$c
  done
  if [ -d $S/legs/hu_new/hu_breakage_index/sandstone2_89 ]; then
    $PY $S/hu_compare.py a8 > $S/K1_a8.log 2>&1; gate a8 $? $S/K1_a8.log
    $PY $S/hu_compare.py report > $S/K1_m4_m6.log 2>&1; report "m4-m6" $S/K1_m4_m6.log
  fi
  finish
  if [ $FAILED = 0 ]; then say "STAGE 1 (Hu section) DONE -- every gate held; STOP for the user (stage 2 needs the go)"
  else say "STAGE 1 (Hu section) DONE WITH A FAILED GATE -- see above"; fi
  exit 0
fi

if [ $STAGE = 1 ]; then
  say "P2-29 stage 1 chain; $(pmset -g batt | head -1)"
  $PY $S/z0.py > $S/Z0.txt 2>&1; gate z0 $? $S/Z0.txt
  [ $FAILED = 1 ] && { say "STOP: z0 failed"; finish; exit 1; }
  env EIGEN=$ROOT/ext CPLUS_INCLUDE_PATH=/opt/homebrew/include \
      LIBRARY_PATH=/opt/homebrew/lib:/opt/homebrew/Cellar/gcc/15.2.0_1/lib/gcc/current make -j8 > $S/build_1b.log 2>&1 \
    || { say "STOP: the build failed ($S/build_1b.log)"; finish; exit 1; }
  mkdir -p $S/bin_1 && cp $BIN $S/bin_1/
  { echo "P2-29 stage 1 (p2_29a): microstructure.phaseN.law = zhang; Table.H's sharp band (the user, 2026-11-10)"
    shasum -a 256 $BIN; shasum -a 256 src/Integrator/MMWSpalling.H src/Numeric/Material/Table.H; } > $S/BINARY_1.txt
  say "binary $(shasum -a 256 $BIN | cut -c1-16) (bin_1/, BINARY_1.txt); bin_pre $(shasum -a 256 $BPRE | cut -c1-16)"
  # SR 9: the outputs this chain overwrites
  [ -d $S/archive/unit_beam_tau_output ] || { mv tests/MMWSpalling/unit/beam_tau/output $S/archive/unit_beam_tau_output \
      && say "archived (moved) unit/beam_tau/output"; }
  for t in hu_conduction hu_thermoelastic; do
    [ -d $S/archive/${t}_output ] || { cp -rp $H/$t/output $S/archive/${t}_output && say "archived (copied) $H/$t/output"; }
  done
  # a1: the face gate (every grimsel/constant leg identical; nothing in it is zhang or a sharp band)
  $PY $W/face_gate.py run p2_29a > $S/face_1.log 2>&1
  $PY $W/face_gate.py cmp p2_28 p2_29a > $S/GATE_a1.txt 2>&1
  grep -q "^VERDICT: IDENTICAL" $S/GATE_a1.txt; gate a1 $? $S/GATE_a1.txt
  # a2: latent_heat + its zhang leg
  $PY tests/MMWSpalling/unit/latent_heat/test > $S/K1_latent_heat.log 2>&1; gate a2 $? $S/K1_latent_heat.log
  # the source locks, on the shipped binary (cheap; not a packet row)
  for t in retired_keys lateral_flux_lock buckle_lock interlocks; do
    $PY tests/MMWSpalling/unit/$t/test > $S/K1_$t.log 2>&1; gate "lock $t" $? $S/K1_$t.log
  done
  # a7 + m2: the five fixtures, old (bin_pre) vs new (bin_1)
  if [ $FAILED = 0 ]; then
    $PY $S/fixtures.py run > $S/K1_fix_run.log 2>&1; say "  fixtures run rc $? ($S/K1_fix_run.log)"
    $PY $S/fixtures.py score > $S/K1_fixtures.log 2>&1; gate a7 $? $S/K1_fixtures.log
    say "  m2: $(grep '^m2' $S/K1_fixtures.log | cut -c1-220)"
  else say "  SKIPPED (a gate failed): a7 / m2"; fi
  # the Zhang legs
  leg ZO $BPRE $Z/input $S/legs/ZO/run solver=explicit timestep=0.05 amr.plot_int=200
  leg ZF $BIN $Z/input_phase $S/legs/ZF/run "geometry.prob_lo=0.0 0.0 0.0" "geometry.prob_hi=0.10 0.10 0.025" \
      "amr.n_cell=128 128 32" energy_ledger.enabled=1
  leg ZI $BPRE $Z/input $S/legs/ZI/run   # the implicit record: the old deck as written (m3)
  leg ZQ $BIN $Z/input_phase $S/legs/ZQ/run energy_ledger.enabled=1
  leg ZQ1 $BIN $Z/input_phase $S/legs/ZQ1/run energy_ledger.enabled=1 max_step=2999
  leg ZL $BIN $S/ZL.input $S/legs/ZL/run energy_ledger.enabled=1
  if [ -d $S/legs/ZQ1/run ] && [ -d $S/legs/ZL/run ]; then
    for r in a6 a3 a4; do $PY $S/zcompare.py $r > $S/K1_$r.log 2>&1; gate $r $? $S/K1_$r.log; done
    $PY $Z/test $S/legs/ZQ/run > $S/K1_a5_zhang.log 2>&1; gate "a5 zhang_oglesby/test on ZQ" $? $S/K1_a5_zhang.log
    for r in m3 m8; do $PY $S/zcompare.py $r > $S/K1_$r.log 2>&1; report $r $S/K1_$r.log; done
  fi
  # the Hu decks: old on bin_pre, new on bin_1
  for side in old new; do
    for d in hu_conduction/input_granite2 hu_conduction/input_sandstone2 hu_thermoelastic/input_granite2 \
             hu_thermoelastic/input_sandstone2 hu_thermoelastic_amr/input_granite2 hu_thermoelastic_amr/input_sandstone2 \
             hu_breakage_index/input_granite2_37 hu_breakage_index/input_sandstone2_37 hu_breakage_index/input_sandstone2_89; do
      t=$(dirname $d); c=$(basename $d); c=${c#input_}
      if [ $side = old ]; then leg "hu_old/$t/$c" $BPRE $H/$d $S/legs/hu_old/$t/$c
      else leg "hu_new/$t/$c" $BIN $H/${d}_phase $S/legs/hu_new/$t/$c; fi
    done
  done
  if [ -d $S/legs/hu_new/hu_breakage_index/sandstone2_89 ]; then
    $PY $S/hu_compare.py a8 > $S/K1_a8.log 2>&1; gate a8 $? $S/K1_a8.log
    $PY $S/hu_compare.py report > $S/K1_m4_m6.log 2>&1; report "m4-m6" $S/K1_m4_m6.log
  fi
  finish
  if [ $FAILED = 0 ]; then say "STAGE 1 DONE -- every gate held; STOP for the user (stage 2 needs the go)"
  else say "STAGE 1 DONE WITH A FAILED GATE -- see above"; fi
fi
