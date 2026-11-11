#!/bin/bash
# P2-28 RECORD: the tree gate retired with the column removal path (one removal path, the face list); p2_27.md5 is the column era's last list, its outputs moved with their tests to retired_units/. Kept, not run: face_gate.py is the gate.
# Copied from the D2c scratchpad on 2026-09-18 (D2e Goal 0); paths point here.
# targeted D2b witnesses: Level_0 Cell_D + event/cluster CSVs + thermo
# P2-27: + unit/jet_d2c/output (the sweep always ran the test; its rec_fixed legs retired at P2-26, so it enters
#        the oracle -- declared only-in-new at p2_27a)
# P2-27 stage 2: - unit/spall_event/output{,_no_event} (the damage-law removal's test, retired with the lineage;
#        its 8 paths are declared only-in-ref at p2_27b)
cd /Users/tzetze20/amr_tools/alamo/tests/MMWSpalling
D="unit/robin_face/output unit/robin_pinned/output unit/robin_jet/output unit/jet_enthalpy/output unit/beam_void_closure/output unit/scalar_flaw/output unit/robin_feet/output validation/meier/sp_meier_pilot/output/dev2d studies/s1_surface_resolution/output_smoke unit/jet_d2c/output"
for d in $D; do
  find $d -path '*.old*' -prune -o \( -path '*Level_0/Cell_D_*' -o -name '*_removal_events.csv' -o -name '*_h_col_events.csv' -o -name '*_clusters.csv' -o -name '*_jet_profile.csv' -o -name thermo.dat \) -print0 | xargs -0 md5 -r
done | sort -u -k2 > $1
