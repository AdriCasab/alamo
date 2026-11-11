#!/bin/bash
# Copied from the D2c scratchpad on 2026-09-18 (D2e Goal 0); paths point here.
# targeted D2b witnesses: Level_0 Cell_D + event/cluster CSVs + thermo
cd /Users/tzetze20/amr_tools/alamo/tests/MMWSpalling
D="unit/robin_face/output unit/robin_pinned/output unit/robin_jet/output unit/jet_enthalpy/output unit/beam_void_closure/output unit/scalar_flaw/output unit/spall_event/output unit/spall_event/output_no_event unit/robin_feet/output validation/meier/sp_meier_pilot/output/dev2d studies/s1_surface_resolution/output_smoke"
for d in $D; do
  find $d -path '*.old*' -prune -o \( -path '*Level_0/Cell_D_*' -o -name '*_removal_events.csv' -o -name '*_h_col_events.csv' -o -name '*_clusters.csv' -o -name '*_jet_profile.csv' -o -name thermo.dat \) -print0 | xargs -0 md5 -r
done | sort -u -k2 > $1
