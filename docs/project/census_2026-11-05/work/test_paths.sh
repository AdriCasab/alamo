#!/bin/bash
# Census helper (2026-11-05): per test dir, which column-path-only features its inputs/test switch on.
cd /Users/tzetze20/amr_tools/alamo/tests/MMWSpalling
for d in unit/* validation/*/* validation/p0_steady_column studies/s1_surface_resolution; do
  [ -d "$d" ] || continue
  files=$(ls $d/input* $d/test* $d/*.py 2>/dev/null)
  [ -z "$files" ] && continue
  c=$(cat $files 2>/dev/null)
  tags=""
  echo "$c" | grep -qE "vapor\.enabled\s*=\s*1|vapor\.enabled=1" && tags="$tags vapor"
  echo "$c" | grep -qE "melt_aware\s*=\s*1|melt_aware=1" && tags="$tags melt_aware"
  echo "$c" | grep -qE "bit\.enabled\s*=\s*1|bit\.enabled=1" && tags="$tags bit"
  echo "$c" | grep -qE "surface_normal\s*=\s*1|surface_normal=1" && tags="$tags surface_normal"
  echo "$c" | grep -qE "foot_body_clearance\s*=\s*1|foot_body_clearance=1" && tags="$tags feet_clearance"
  echo "$c" | grep -qE "robin_form\s*=\s*pinned|robin_form=pinned" && tags="$tags pinned"
  echo "$c" | grep -qE "mode\s*=\s*prescribed_T|mode=prescribed_T" && tags="$tags prescribed_T"
  echo "$c" | grep -qE "spallation\.model\s*=\s*sp_weibull|spallation\.model=sp_weibull" || { echo "$c" | grep -qE "spall\.enabled\s*=\s*1" && tags="$tags damage_law"; }
  echo "$c" | grep -qE "flake_coherence_length\s*=\s*0\.[0-9]*[1-9]" && tags="$tags coherence"
  echo "$c" | grep -qE "stress_source\s*=\s*mechanics|el\.type\s*=\s*static" && tags="$tags stress_mf/static"
  echo "$c" | grep -qE "beam\.P0\s*=\s*[1-9]|beam\.P0=[1-9]" && tags="$tags BEAM"
  echo "$c" | grep -qE "face_list\s*=\s*1|face_list=1|face_list\": \"1|face_list\": str" && tags="$tags FACE_LIST"
  echo "$c" | grep -qE "spall\.enabled\s*=\s*1|spall\.enabled=1|vapor\.enabled\s*=\s*1" || tags="$tags (no-removal?)"
  echo "$d:$tags"
done
