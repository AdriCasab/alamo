#!/bin/bash
# Census helper (2026-11-05): for each tree-gate input, print the keys that pick the removal path.
cd /Users/tzetze20/amr_tools/alamo/tests/MMWSpalling
for f in validation/meier/sp_meier_pilot/input_2d_dev unit/robin_face/input_form unit/robin_face/input_ledger unit/robin_pinned/input_pinned unit/robin_jet/input_jet unit/jet_enthalpy/input_quarter unit/jet_enthalpy/input_strip unit/robin_feet/input_column unit/robin_feet/input_feet unit/beam_void_closure/input unit/scalar_flaw/input unit/spall_event/input unit/spall_event/input_no_event unit/removal_increment/input unit/removal_increment/input_presc studies/s1_surface_resolution/input; do
  echo "## $f"
  grep -E "^\s*(surface_patch\.(enabled|mode|face_list|feet|foot_body_clearance|robin_form|nozzle_descent)|spallation\.model|spall\.(enabled|stress_source|surface_T|removal_energy|melt_aware|surface_normal|flake_coherence_length|buckle_span)|vapor\.enabled|bit\.enabled|beam\.(P0|power|P_expr|power_expr)|energy_ledger\.enabled|el\.type|amr\.max_level|weibull\.enabled)\b" "$f" | sed 's/#.*//' | tr -s ' ' | tr '\n' ';'
  echo
done
