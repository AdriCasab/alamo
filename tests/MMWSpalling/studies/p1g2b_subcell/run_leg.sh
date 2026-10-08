#!/bin/bash
# P1g-2b deliverable 3: one Meier leg, Leg B's argv VERBATIM plus the ledger's keys.
#
#   usage: run_leg.sh CTRL3 | CTRL4
#
# ⚠⚠ THE ARGV IS NOT RETYPED -- IT IS READ FROM LEG B'S OWN argv.sh, between
# `input_drilling` and `plot_file=`. That argv is ~2 kB of nested h_expr with
# balanced parentheses and a hand-set standoff decay; retyping or regenerating it
# would be a silent chance to differ from the leg every offline number in (a1)-(a4)
# was measured on. Proven approach, carried from p1g3_thickness/run_meier.sh.
#
# ⚠⚠ NO SOURCE IS EDITED BY THIS STEP, so the binary must hash the SAME before and
# after. Both hashes are recorded per leg and (b1) asserts them. If they differ,
# something rebuilt mid-campaign and the leg is void -- the D2q-2b failure mode.
#
# ⚠ `spall.lateral_subcell=0` is written EXPLICITLY on CTRL3 rather than omitted, so
# the diff against CTRL2 is one visible key. It is byte-identical to omission:
# MMWSpalling.H:5834 is `query_default(..., false)` and the member default at :9573
# is `false`, with NO `requires` abort (removed after P1g-1's negative control
# returned rc 6). ⚠ But the thermo HEADER still differs -- :875-884 registers
# lat_bank, lat_bank_hf, lat_accum_sum AND lat_void_colmax only under the key -- so
# score_subcell.py must read thermo.dat BY COLUMN NAME.
#
# ⚠ caffeinate: these runs are NOT checkpointed; a host sleep loses the whole leg.
set -u
ROOT=/Users/tzetze20/amr_tools/alamo
S=$ROOT/tests/MMWSpalling/studies/p1g2b_subcell
LEGB=$ROOT/tests/MMWSpalling/studies/p1f_lateral/output_legB_fix/argv.sh
BIN=$ROOT/bin/mmwspalling-3d-g++
NAME="${1:?usage: run_leg.sh CTRL3|CTRL4}"

[ -f "$LEGB" ] || { echo "FATAL: no archived argv at $LEGB"; exit 1; }
[ -x "$BIN" ]  || { echo "FATAL: no binary at $BIN"; exit 1; }

case "$NAME" in
  # banking OFF on the SAME binary as CTRL2 => the attribution, and the drift arm.
  CTRL3) NEWKEYS="spall.lateral_subcell=0" ;;
  # banking + THE THERMAL HALF. P1g-2's shift at Meier scale for the FIRST time.
  # ⚠ lateral_thermal needs the banking conjunction (bank_on && sp_lateral_thermal),
  # so subcell=1 is carried too -- thermal alone would be inert.
  CTRL4) NEWKEYS="spall.lateral_subcell=1 spall.lateral_thermal=1" ;;
  *) echo "FATAL: unknown leg $NAME"; exit 1 ;;
esac

mkdir -p "$S/$NAME"

# Everything after the input path and before plot_file=, verbatim.
# ⚠ python only SPLITS the line -- it rewrites no token.
MID=$(python3 - "$LEGB" <<'PY'
import sys
s = open(sys.argv[1]).read().strip()
i = s.index("input_drilling") + len("input_drilling")
j = s.rindex("plot_file=")
sys.stdout.write(s[i:j].strip())
PY
)
[ -n "$MID" ] || { echo "FATAL: could not extract the argv body"; exit 1; }

INP=$ROOT/tests/MMWSpalling/validation/meier/sp_meier_pilot/input_drilling
CMD="caffeinate -dis mpirun --oversubscribe --bind-to none -np 4 \
  $BIN $INP $MID $NEWKEYS plot_file=$S/$NAME/FH"

shasum -a 256 "$BIN" > "$S/$NAME.sha_before"
printf '%s\n' "$CMD" > "$S/$NAME.argv.sh"
cat "$S/$NAME.sha_before" >> "$S/$NAME.argv.sh"

echo "=== $NAME start $(date +%T)   keys: $NEWKEYS"
eval "$CMD" > "$S/$NAME.log" 2>&1
RC=$?
shasum -a 256 "$BIN" > "$S/$NAME.sha_after"
echo "=== $NAME rc=$RC  $(date +%T)"
if ! cmp -s "$S/$NAME.sha_before" "$S/$NAME.sha_after"; then
  echo "GUARD FAILURE: the binary CHANGED during $NAME -- the leg is VOID."
  echo "   before: $(cat "$S/$NAME.sha_before")"
  echo "   after : $(cat "$S/$NAME.sha_after")"
  exit 2
fi
echo "=== $NAME binary unchanged across the leg"
exit $RC
