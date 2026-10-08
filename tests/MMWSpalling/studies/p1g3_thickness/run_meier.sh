#!/bin/bash
# P1g-3 deliverable 10 / (c1): the two Meier legs.
#
#   usage: run_meier.sh MEIER_CTRL | MEIER_BOTH
#
# ⚠⚠ THE ARGV IS NOT RETYPED -- IT IS READ FROM LEG B'S OWN argv.sh. The packet
# says "Leg B's argv verbatim", and that argv is ~2 kB of nested h_expr with
# balanced parentheses. Retyping or regenerating it would be a silent chance to
# differ from the leg every offline number in this step was measured on. So the
# script strips only the launcher prefix and the trailing plot_file from the
# archived argv and keeps everything between them byte for byte.
#
# ⚠⚠ `spall.lateral_subcell=1` GOES ON BOTH LEGS. Leg B's argv has
# `lateral_drive=1` and `lateral_removal=1` but NOT `lateral_subcell`, and
# without it regime 6 voids a whole 2 mm cell, so a derived wall thickness never
# reaches the geometry and (c1) would compare two whole-cell runs. That is also
# why MEIER_CTRL must be a NEW RUN and not the archived Leg B.
#
# ⚠ The ONLY difference between the two legs is the two new key lines, so any
# bore difference is attributable to the derived thickness and nothing else.
set -u
ROOT=/Users/tzetze20/amr_tools/alamo
S=$ROOT/tests/MMWSpalling/studies/p1g3_thickness
LEGB=$ROOT/tests/MMWSpalling/studies/p1f_lateral/output_legB_fix/argv.sh
NAME="${1:?usage: run_meier.sh MEIER_CTRL|MEIER_BOTH}"

[ -f "$LEGB" ] || { echo "FATAL: no archived argv at $LEGB"; exit 1; }

case "$NAME" in
  MEIER_CTRL) NEWKEYS="" ;;
  MEIER_BOTH) NEWKEYS="spall.flake_length=flux spall.flake_fraction=weibull" ;;
  # ⚠⚠ MEIER_LAT carries the SAME KEYS as MEIER_BOTH -- the difference is the
  # BINARY. MEIER_BOTH ran before the LATERAL site was wired, so its wall
  # thickness stayed pinned at the 1.500 mm cap on 68 % of firings and it is a
  # VERTICAL-ONLY leg. MEIER_LAT runs with both sites derived. Separate names so
  # the two are never confused, and each leg records its own binary sha.
  MEIER_LAT)  NEWKEYS="spall.flake_length=flux spall.flake_fraction=weibull" ;;
  # ⚠ CTRL2 re-runs the control on the SAME binary as MEIER_LAT. The toolchain is
  # not byte-reproducible, so scoring a leg from one build against a control from
  # another is the provenance hole the P1g-0b review caught. They run in parallel,
  # so this costs no wall time.
  MEIER_CTRL2) NEWKEYS="" ;;
  *) echo "FATAL: unknown leg $NAME"; exit 1 ;;
esac

mkdir -p "$S/$NAME"

# Everything after the input path and before plot_file=, verbatim.
# ⚠ python is used only to split the line -- it does not rewrite any token.
MID=$(python3 - "$LEGB" <<'PY'
import re, sys
s = open(sys.argv[1]).read().strip()
i = s.index("input_drilling") + len("input_drilling")
j = s.rindex("plot_file=")
sys.stdout.write(s[i:j].strip())
PY
)
[ -n "$MID" ] || { echo "FATAL: could not extract the argv body"; exit 1; }

INP=$ROOT/tests/MMWSpalling/validation/meier/sp_meier_pilot/input_drilling

# ⚠ caffeinate: these runs are NOT checkpointed, so a host sleep loses the whole
# leg and it restarts from t = 0.
CMD="caffeinate -dis mpirun --oversubscribe --bind-to none -np 4 \
  $ROOT/bin/mmwspalling-3d-g++ $INP $MID \
  spall.lateral_subcell=1 $NEWKEYS \
  plot_file=$S/$NAME/FH"

# record exactly what ran, beside the output, before running it
printf '%s\n' "$CMD" > "$S/$NAME/argv.sh"
shasum -a 256 "$ROOT/bin/mmwspalling-3d-g++" >> "$S/$NAME/argv.sh"

echo "=== $NAME start $(date +%T)"
eval "$CMD" > "$S/$NAME.log" 2>&1
RC=$?
echo "=== $NAME rc=$RC  $(date +%T)"
exit $RC
