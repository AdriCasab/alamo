#!/bin/bash
# Copied from the D2c scratchpad on 2026-09-18 (D2e Goal 0); paths point here.
S=/Users/tzetze20/amr_tools/alamo/tests/MMWSpalling/studies/tree_sweep_0922/witness
W=$S/$1_wit; mkdir -p $W
cd /Users/tzetze20/amr_tools/alamo
VENV=/Users/tzetze20/Desktop/code/.venv/bin/python
RUN="mpirun --oversubscribe --bind-to none -np 4"
T=tests/MMWSpalling
run_test () { name=$1; shift; echo "=== $name start $(date +%T)"; SECONDS=0; "$@" > $W/$name.log 2>&1; echo "=== $name rc=$? ${SECONDS}s"; }
run_test dev2d $RUN bin/mmwspalling-3d-g++ $T/validation/meier/sp_meier_pilot/input_2d_dev
for t in robin_face robin_pinned robin_jet jet_enthalpy robin_feet beam_void_closure scalar_flaw spall_event; do
  run_test $t $VENV $T/unit/$t/test
done
run_test s1_smoke $VENV $T/studies/s1_surface_resolution/run_sweep.py --smoke --force
$S/hash.sh $S/$1.md5
echo ALLDONE
